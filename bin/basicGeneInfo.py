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
# 
# Author: Joel Richardson
#
import sys
from subprocess import Popen
import argparse
import json
import heapq
import itertools
import re
import os
from AGRlib import stripNulls, buildMetaObject, sql
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
# Goes to PantherDB to get the download file, then parses the file to generate a mapping from MGI ids to Panther IDs.
# Returns the map.
PANTHERURL="ftp://ftp.pantherdb.org/ortholog/current_release/RefGenomeOrthologs.tar.gz"
def getPantherIds () :
    def parseMouseId (s) :
        idPart = s.split("|")[1]
        return "MGI:" + idPart.split("=")[-1]

    def parseLine(line):
        parts = line.split()
        pthrId = parts[-1]
        if parts[0].startswith('MOUSE'):
            return (parseMouseId(parts[0]), pthrId)
        elif parts[1].startswith('MOUSE'):
            return (parseMouseId(parts[1]), pthrId)

    cmd = 'curl -o "RefGenomeOrthologs.tar.gz" -z "RefGenomeOrthologs.tar.gz" "%s"' % PANTHERURL
    sp = Popen(cmd, shell=True)
    rc = sp.wait()
    # tar outputs file names to stdout. Redirect to /dev/null so these don't end up
    # in the json file.
    cmd = 'tar -xvf RefGenomeOrthologs.tar.gz > /dev/null'
    sp = Popen(cmd, shell=True)
    rc = sp.wait()
    mgi2panther = {}
    with open('RefGenomeOrthologs','r') as fd:
        for line in fd:
            res = parseLine(line[:-1])
            if res:
                mgi2panther[res[0]] = res[1]
    return mgi2panther

#-----------------------------------

GLOBALTAXONID = os.environ["GLOBALTAXONID"]
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
    pid = obj['pantherId']
    if pid:
      xrefs.add(('PANTHER', pid))
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

#
def main():
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
    ]

    mgi2panther = getPantherIds()

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
                r['soTermId'] = mcv2so[r['_mcv_term_key']]
                r['pantherId'] = mgi2panther.get(r['markerId'], None)
                id2gene[r['_marker_key']] = r
        else:
            for r in sql(q):
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

main()
