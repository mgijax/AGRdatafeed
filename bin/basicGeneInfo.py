#!/usr/bin/env python2.7 
#
# basicGeneInfo.py
#
# Script to dump basic gene information from MGI in the AGR standard JSON format.
# The format is described here:
#       https://github.com/alliance-genome/agr_schemas
#
# Usage:
# To dump all genes and pseudogenes:
#       % python basicGeneInfo.py > FILE
# To dump specific genes/pseudogenes:
#       % python basicGeneInfo.py MGI:96449 MGI:96677 MGI:2685845
# 
# This script uses MouseMine webservices API.
# Implementation approach. A fair amount of detail is needed for each gene.
# In the db, this info requires traversing to subobjects (and beyond).
# One option is to specify 'jsonobjects' as the return type, and the server
# will put things together for you in nice object form. My experience is that
# this breaks the server for very large queries like these.
#
# This leads to the approach taken here, which is to break up the query 
# into simple(r) pieces each of which returns a stream of parts sorted by object id.
# The streams are exectuted in parallel, then merged based on object id.
# 
# Example: for each gene we need its synonyms (among many other parts). These are
# stored in a separate table, and a given gene can have zero or more synonyms.
# In the one-giant-query approach, we'd tack on another join. 
# Because a given gene can have 0 synonyms, the join would have to be outer.
# In the parallel streams approach, we add another query that returns each synonym with
# its gene's ID, sorted on that ID. All the other streams are also composed of pieces that
# returns streams of parts sorted by gene ID.
# The bottom line is that the server can execute the collection of simple queries much faster and
# with less memory then the one-giant-query approach.
#
# Author: Joel Richardson
#
import sys
import argparse
import json
import heapq
import itertools
import re
import os
from AGRlib import stripNulls, buildMetaObject, sql, makeOneOfConstraint
from AGRqlib import qMcvTerms, qGenes, qGeneHasPhenotype, qGeneHasImpc, qGeneSynonyms, qGeneHasExpression, qGeneHasExpressionImage, qGeneLocations, qGeneProteinIds, qGeneXrefs, qGeneSecondaryIds

#-----------------------------------
MCV2SO_AUX = {
    #"complex/cluster/region"
    6238175 : "SO:0000110",
    #"cytogenetic marker"
    6238176 : "SO:0000110",
    #"BAC/YAC end"
    6238177 : "SO:0000110",
    #"other genome feature"
    6238178 : "SO:0000110",
    #"DNA segment"
    6238179 : "SO:0000110",
    #"unclassified gene"
    6238184 : "SO:0000704",
    #"other feature type"
    6238185 : "SO:0000110",
    #"unclassified non-coding RNA gene"
    6238186 : "SO:0000704",
    #"unclassified cytogenetic marker"
    7222413 : "SO:0000110",
    #"unclassified other genome feature"
    7648969 : "SO:0000110",
    #"mutation defined region"
    11928467 : "SO:0000110",
}

#----------------------------------
XREF_DBS = {
    "Entrez Gene": "NCBI_Gene",
    "Ensembl Gene Model": "ENSEMBL"
}

#-----------------------------------

taxon         = os.environ["TAXONID"]
GLOBALTAXONID = os.environ["GLOBALTAXONID"]
MYGENEURL     = os.environ["MYGENEURL"]
SAMPLEIDS     = os.environ["SAMPLEIDS"].split()
MGD_OLD_PREFIX= os.environ["MGD_OLD_PREFIX"]

# For AGR, old style MGD ids must be given a distinct global prefix (we're using 'MGD_old:') and have 
# a stanza describing that kind of ID in the resourceDescriptors.yaml file. This routine adds the 
# prefix, if appropriate,
def formatSecondary(identifier):
  return (MGD_OLD_PREFIX if identifier.startswith("MGD-") else "") + identifier

# In the MGI fewi, mouse genes link to a MyGenes wiki page which is a human readable description.
# The MyGenes page is for the HUMAN ortholog of the mouse gene.
# In the database, human genes have a cross reference to MyGenes, where the "id" is the part needed
# to fill out the complete URL (the base part if constant). 
# The link is only displayed if the mouse/human gene have a 1:1 orthology relationship.
# Here we use the same logic to construct a link (or not) for the given mouse gene (obj).
#
def formatMyGeneLink(obj):
    mgl = obj.get("myGeneLink", [None])[0]
    if not mgl:
        return None
    symbol = mgl['homologues.homologue.crossReferences.identifier']
    # FIXME: Currently, the agr schema for global id's doesn't allow ":" in the suffix part.
    # A few of these wiki links do, so we'll filter them out for now.
    # TODO: Ask DQMs to change the globalId pattern. Then remove this filter.
    return None if ":" in symbol else {"id": "WIKIP:"+symbol}

# Selects the xrefs to be exported for the object and formats them according to the spec.
#       - restricts which xrefs are exported
#       - translates provider name
#       - packs provider name and id into a object
#       - ensures uniqueness 
# Returns a list of cross reference objects.
#
def formatXrefs(obj):
    xrefs = set()
    for x in obj.get("xrefs",[]):
      dp = XREF_DBS.get(x["ldbName"], None)
      if dp:
        xrefs.add((dp, x["accid"]))
    for x in obj.get('proteinIds', []):
        p = x.get('proteinId','')
        if p:
            xrefs.add(("UniProtKB", p))
    pid = obj.get('pantherId', [None])[0]
    if pid:
      xrefs.add(('PANTHER', pid['homologues.crossReferences.identifier']))
    xrefs = list(xrefs)
    xrefs.sort()
    # new xref format for 1.0.0.0. Includes 2 parts: the id, and a list of 
    # page-tags (see resourceDescriptors.yaml)
    xrs = [{"id": x[0]+":"+x[1]} for x in xrefs]
    # add xrefs to MGI pages for this gene
    pgs = ["gene","gene/references"]
    if obj.get('expressed', None):
        pgs.append("gene/expression")
    if obj.get('expressedImages', None):
        pgs.append("gene/expression_images")
    if obj.get('hasPheno', False):
        pgs.append('gene/phenotype')
    if obj.get('hasImpc', False):
        pgs.append('gene/phenotypes_impc')
    xrs.append({"id": obj["markerId"], "pages":pgs })
    # add xref to MyGene page (if applicable)
    mgl = formatMyGeneLink(obj)
    if mgl: xrs.append(mgl)
    #
    return xrs


# Format the genome location for the obj. The agr standard is to allow multiple
# locations per object, but we will only have one. Even so, we have to return a list.
# 
def formatGenomeLocation(obj):

    loc = obj["location"][0]
    if loc["genomicchromosome"]:
        return [{
            "assembly"          : loc['assembly'],
            "chromosome"        : loc["genomicchromosome"],
            "startPosition"     : int(loc['startcoordinate']),
            "endPosition"       : int(loc['endcoordinate']),
            "strand"            : loc['strand']
        }]
    else:
        if loc["chromosome"] and loc["chromosome"] != 'UN':
          return [{
            "assembly"          : '',
            "chromosome"        : loc["chromosome"]
        }]

def formatDescription (obj) :
    d = obj["description"]
    if not d:
      return None
    return "PHENOTYPE: %s [provided by MGI curators]" % d

def getJsonObj(obj):
      #try:
          synonyms = []
          if 'synonyms' in obj:
              synonyms = [ r['synonym'] for r in obj['synonyms'] ]
          #
          secondaryIds = []
          if 'secondaryIds' in obj:
              secondaryIds = [ formatSecondary(r['accid']) for r in obj['secondaryIds'] ]
          #
          basicGeneticEntity = stripNulls({
            "primaryId"         : obj["markerId"],
            "taxonId"           : GLOBALTAXONID,
            "secondaryIds"      : secondaryIds,
            "synonyms"          : [ s for s in synonyms if s != obj["symbol"] and s != obj["name"] ],
            "crossReferences"   : formatXrefs(obj),
            "genomeLocations"   : formatGenomeLocation(obj),
          })
          return stripNulls({
            "basicGeneticEntity": basicGeneticEntity,
            "symbol"            : obj["symbol"],
            "name"              : obj["name"],
            "geneSynopsis"      : formatDescription(obj),
            "soTermId"          : obj["soTermId"],
          })
      #except:
          #sys.stderr.write('ERROR in getJsonObj. obj=' + str(obj) + '\n')
          #sys.exit(1)
#
def parseCmdLine():
    parser = argparse.ArgumentParser(description='Dumps basic gene information to a JSON file.')

    parser.add_argument(
      '-s','--sample',
      action='store_true',
      default=False,
      help='Generate sample output',
      required=False)

    parser.add_argument(
      'identifiers', 
      metavar='ids',
      nargs='*',
      help='Specific MGI ids to dump.')

    args = parser.parse_args()
    if args.sample:
      args.identifiers.extend(SAMPLEIDS)
    return args
  
#
def main(args):
    ##
    qs = [
        ('gene',            qGenes),
        ('synonyms',        qGeneSynonyms),
        ('secondaryIds',    qGeneSecondaryIds),
        ('expressed',       qGeneHasExpression),
        ('expressedImages', qGeneHasExpressionImage),
        ('location',        qGeneLocations),
        ('proteinIds',      qGeneProteinIds),
        ('xrefs',           qGeneXrefs),
        #('pantherId',       mousePantherIds),
    ]

    # set of markers with phenotype annotations
    hasPheno = set()
    for r in sql(qGeneHasPhenotype):
        hasPheno.add(r['_marker_key'])

    # set of markers with alleles in the IMPC collection
    hasImpc = set()
    for r in sql(qGeneHasImpc):
        hasImpc.add(r['_marker_key'])

    # Mapping from MCV term key to SO id.
    # Initialize from the hard coded mappings
    # then load what's in the db (which is incomplete).
    mcv2so = MCV2SO_AUX.copy()
    so_re = re.compile(r'SO:[0-9]+')
    for r in sql(qMcvTerms):
        m = so_re.search(r['note'])
        if m:
            mcv2so[r['_term_key']] = m.group(0)

    id2gene = {}
    for label, q in qs:
        if label == 'gene':
            for r in sql(q):
                r = dict(r)
                r['soTermId'] = mcv2so[r['_mcv_term_key']]
                id2gene[r['_marker_key']] = r
        else:
            for r in sql(q):
                r = dict(r)
                obj = id2gene.get(r['_marker_key'], None)
                if obj:
                    obj.setdefault(label,[]).append(r)
    print('{\n  "metaData": %s,\n  "data": [' % json.dumps(buildMetaObject(), indent=2))
    first=True
    for i in id2gene:
        obj = id2gene[i]
        if not first: print(',', end='')
        obj["hasPheno"] = obj["_marker_key"] in hasPheno
        obj["hasImpc"] = obj["_marker_key"] in hasImpc
        print(json.dumps(getJsonObj(obj), indent=2))
        first = False
    print(']\n}')

mousePantherIds = '''
    <query
      name="BGI_mousePantherIds"
      model="genomic"
      view="
        Gene.primaryIdentifier
        Gene.homologues.crossReferences.identifier
        "
      sortOrder="Gene.primaryIdentifier asc"
      >
      %(extraConstraint)s
      <constraint path="Gene.organism.taxonId" op="=" value="10090"/>
      <constraint path="Gene.homologues.dataSets.name" op="=" value="Panther data set"/>
      <constraint path="Gene.homologues.type" op="!=" value="paralogue"/>
      <constraint path="Gene.dataSets.name" op="=" value="Mouse Gene Catalog from MGI"/>
      <constraint path="Gene.homologues.homologue.organism.taxonId" op="=" value="9606"/>
    </query>
    '''

main(parseCmdLine())
