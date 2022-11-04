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
import urllib.request, urllib.parse, urllib.error
import json
import heapq
import itertools
import re
import os
from AGRlib import stripNulls, buildMetaObject, doQuery, makeOneOfConstraint, dataProviders

#-----------------------------------

MOUSEMINE     = os.environ["MOUSEMINEURL"]
taxon         = os.environ["TAXONID"]
GLOBALTAXONID = os.environ["GLOBALTAXONID"]
MYGENEURL     = os.environ["MYGENEURL"]
SAMPLEIDS     = os.environ["SAMPLEIDS"].split()
MGD_OLD_PREFIX= os.environ["MGD_OLD_PREFIX"]

# In MouseMine, synonyms and secondary ids are lumped together as "synonyms". 
# This function distinguishes a synonym value as being either a secondary id or not.
#
def isSecondaryId(identifier):
        return identifier.startswith("MGI:") or identifier.startswith("MGD-")

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
      dp = dataProviders.get(x["crossReferences.source.name"], None)
      if dp:
        xrefs.add((dp, x["crossReferences.identifier"]))
    for x in obj.get('proteinIds', []):
        p = x.get('proteins.uniprotAccession','')
        if p:
            xrefs.add(("UniProtKB", p))
    pid = obj.get('pantherId', [None])[0]
    if pid:
      xrefs.add(('PANTHER', pid['homologues.crossReferences.identifier']))
    xrefs = list(xrefs)
    xrefs.sort()
    # new xref format for 1.0.0.0. Includes 2 parts: the id, and a list of page-tags (see resourceDescriptors.yaml)
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
    xrs.append({"id": obj["primaryIdentifier"], "pages":pgs })
    # add xref to MyGene page (if applicable)
    mgl = formatMyGeneLink(obj)
    if mgl: xrs.append(mgl)
    #
    return xrs


# Convert strand value as stored in MouseMine (+1/-1/0) to the AGR standard (+/-/.)
#
def convertStrand(s):
    return "+" if s in ["+1","1"] else "-" if s == "-1" else "."

# Format the genome location for the obj. The agr standard is to allow multiple
# locations per object, but we will only have one. Even so, we have to return a list.
# 
def formatGenomeLocation(chrom, loc):
    if loc:
        return [{
            "assembly"          : loc['chromosomeLocation.assembly'],
            "chromosome"        : chrom,
            "startPosition"     : int(loc['chromosomeLocation.start']),
            "endPosition"       : int(loc['chromosomeLocation.end']),
            "strand"            : convertStrand(loc['chromosomeLocation.strand'])
        }]
    #elif obj.chromosome and obj.chromosome.primaryIdentifier != "UN":
    else:
        if chrom and chrom != 'UN':
          return [{
            "assembly"          : '',
            "chromosome"        : chrom
        }]

def formatDescription (obj) :
    d = obj["description"]
    if not d:
      return None
    i = d.find("PHENOTYPE")
    if i == -1:
      return None
    else:
      return d[i:]

def getJsonObj(obj):
      #try:
          synonyms = [ r['synonyms.value'] for r in obj.get('synonyms',[]) ]
          basicGeneticEntity = stripNulls({
            "primaryId"         : obj["primaryIdentifier"],
            "taxonId"           : GLOBALTAXONID,
            "secondaryIds"      : [ formatSecondary(s) for s in synonyms if isSecondaryId(s) ],
            "synonyms"          : [ s for s in synonyms if not isSecondaryId(s) and s != obj["symbol"] and s != obj["name"] ],
            "crossReferences"   : formatXrefs(obj),
            "genomeLocations"   : formatGenomeLocation(obj.get('chromosome.primaryIdentifier', None), obj.get('location', [None])[0]),
          })
          return stripNulls({
            "basicGeneticEntity": basicGeneticEntity,
            "symbol"            : obj["symbol"],
            "name"              : obj["name"],
            "geneSynopsis"      : formatDescription(obj),
            "soTermId"          : obj["sequenceOntologyTerm.identifier"],
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
    qmods = {
      'extraConstraint' : makeOneOfConstraint('Gene.primaryIdentifier', args.identifiers),
    }
    ##
    qs = [
        ('gene', mouseGenes),
        ('synonyms', mouseSynonyms),
        ('expressed', mouseExpressedGenes),
        ('expressedImages', mouseExpressedGenesWithImages),
        ('location', mouseLocations),
        ('proteinIds', mouseProteinIds),
        ('xrefs', mouseXrefs),
        ('pantherId', mousePantherIds),
        ('myGeneLink', mouseMyGeneLinks),
    ]

    hasPheno = set()
    for r in doQuery(mouseHasPheno, MOUSEMINE):
        hasPheno.add(r['primaryIdentifier'])

    hasImpc = set()
    for r in doQuery(mouseHasImpc, MOUSEMINE):
        hasImpc.add(r['primaryIdentifier'])

    id2gene = {}
    for label, q in qs:
        if label == 'gene':
            for r in doQuery (q % qmods, MOUSEMINE):
                r['mgiid'] = r['primaryIdentifier']
                id2gene[r['primaryIdentifier']] = r
        else:
            for r in doQuery (q % qmods, MOUSEMINE):
                obj = id2gene.get(r['primaryIdentifier'], None)
                if obj:
                    obj.setdefault(label,[]).append(r)
    print('{\n  "metaData": %s,\n  "data": [' % json.dumps(buildMetaObject(), indent=2))
    first=True
    for i in id2gene:
        obj = id2gene[i]
        if not first: print(',', end='')
        obj["hasPheno"] = obj["primaryIdentifier"] in hasPheno
        obj["hasImpc"] = obj["primaryIdentifier"] in hasImpc
        print(json.dumps(getJsonObj(obj), indent=2))
        first = False
    print(']\n}')

mouseGenes = '''
    <query
      name="BGI_mouseGenes"
      model="genomic"
      view="
        Gene.primaryIdentifier
        Gene.symbol
        Gene.name
        Gene.description
        Gene.sequenceOntologyTerm.identifier
        Gene.chromosome.primaryIdentifier
        "
      sortOrder="Gene.primaryIdentifier asc"
      >
      %(extraConstraint)s
      <constraint path="Gene.organism.taxonId" op="=" value="10090"/>
      <constraint path="Gene.dataSets.name" op="=" value="Mouse Gene Catalog from MGI"/>
      <constraint path="Gene.sequenceOntologyTerm.identifier" op="!=" value="SO:0000902"/>
    </query>
    '''

mouseSynonyms = '''
    <query
      name="BGI_mouseSynonyms"
      model="genomic"
      view="
        Gene.primaryIdentifier
        Gene.synonyms.value
        "
      sortOrder="Gene.primaryIdentifier asc Gene.synonyms.value asc"
      >
      %(extraConstraint)s
      <constraint path="Gene.organism.taxonId" op="=" value="10090"/>
      <constraint path="Gene.dataSets.name" op="=" value="Mouse Gene Catalog from MGI"/>
    </query>
    '''

mouseLocations = '''
    <query
      name="BGI_mouseLocations"
      model="genomic"
      view="
        Gene.primaryIdentifier
        Gene.chromosomeLocation.locatedOn.primaryIdentifier
        Gene.chromosomeLocation.start
        Gene.chromosomeLocation.end Gene.chromosomeLocation.strand
        Gene.chromosomeLocation.assembly
        "
      sortOrder="Gene.primaryIdentifier asc"
      >
      %(extraConstraint)s
      <constraint path="Gene.organism.taxonId" op="=" value="10090"/>
      <constraint path="Gene.dataSets.name" op="=" value="Mouse Gene Catalog from MGI"/>
    </query>
    '''

mouseProteinIds = '''
    <query
      name="BGI_mouseProteinIds"
      model="genomic"
      view="
        Gene.primaryIdentifier
        Gene.proteins.uniprotAccession
        "
      sortOrder="Gene.primaryIdentifier asc Gene.proteins.uniprotAccession asc"
      >
      %(extraConstraint)s
      <constraint path="Gene.organism.taxonId" op="=" value="10090"/>
      <constraint path="Gene.dataSets.name" op="=" value="Mouse Gene Catalog from MGI"/>
    </query>
    '''

mouseExpressedGenes = '''
    <query
      name="BGI_mouseExpressedGenes"
      model="genomic"
      view="
        Gene.primaryIdentifier
        "
      sortOrder="Gene.primaryIdentifier asc"
      >
      %(extraConstraint)s
      <constraint path="Gene.organism.taxonId" op="=" value="10090"/>
      <constraint path="Gene.dataSets.name" op="=" value="Mouse Gene Catalog from MGI"/>
      <constraint path="Gene.expression" op="IS NOT NULL"/>
    </query>
    '''

mouseExpressedGenesWithImages = '''
    <query
      name="BGI_mouseExpressedGenesWithImages"
      model="genomic"
      view="
        Gene.primaryIdentifier
        "
      sortOrder="Gene.primaryIdentifier asc"
      >
      %(extraConstraint)s
      <constraint path="Gene.organism.taxonId" op="=" value="10090"/>
      <constraint path="Gene.dataSets.name" op="=" value="Mouse Gene Catalog from MGI"/>
      <constraint path="Gene.expression.image" op="IS NOT NULL"/>
    </query>
    '''

mouseXrefs = '''
    <query
      name="BGI_mouseXrefs"
      model="genomic"
      view="
        Gene.primaryIdentifier
        Gene.crossReferences.source.name
        Gene.crossReferences.identifier
        "
      sortOrder="Gene.primaryIdentifier asc"
      >
      %(extraConstraint)s
      <constraint path="Gene.organism.taxonId" op="=" value="10090"/>
      <constraint path="Gene.dataSets.name" op="=" value="Mouse Gene Catalog from MGI"/>
    </query>
    '''

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

mouseMyGeneLinks = '''
    <query
      name="BGI_mouseMyGeneLinks"
      model="genomic"
      view="
        Gene.primaryIdentifier
        Gene.homologues.homologue.primaryIdentifier
        Gene.homologues.homologue.symbol
        Gene.homologues.homologue.crossReferences.identifier
        "
      sortOrder="Gene.primaryIdentifier asc"
      >
      %(extraConstraint)s
      <constraint path="Gene.organism.taxonId" op="=" value="10090"/>
      <constraint path="Gene.homologues.dataSets.name" op="=" value="Mouse/Human Orthologies from MGI"/>
      <constraint path="Gene.homologues.homologue.crossReferences.source.name" code="E" op="=" value="MyGene"/>
      <constraint path="Gene.dataSets.name" op="=" value="Mouse Gene Catalog from MGI"/>
    </query>
    '''

mouseHasPheno = '''
    <query 
        model="genomic"
        view="Gene.primaryIdentifier"
        >
        <constraint path="Gene.ontologyAnnotations.ontologyTerm" type="MPTerm"/>
        <constraint path="Gene.ontologyAnnotations.ontologyTerm.identifier" op="IS NOT NULL"/>
        </query>
    '''

mouseHasImpc = '''
    <query
        model="genomic"
        view="Gene.primaryIdentifier"
        >
        <constraint path="Gene.alleles.projectCollection" op="=" value="IMPC"/>
        </query>
        '''
main(parseCmdLine())
