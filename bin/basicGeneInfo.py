#!/usr/bin/env python2.7 
#
# basicGeneInfo.py
#
# Script to dump basic gene information from MGI in the AGR standard JSON format.
# The format is described here:
#	https://github.com/alliance-genome/agr_schemas
#
# Usage:
# To dump all genes and pseudogenes:
#	% python basicGeneInfo.py > FILE
# To dump specific genes/pseudogenes:
#	% python basicGeneInfo.py MGI:96449 MGI:96677 MGI:2685845
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
import urllib
import json
import heapq
import itertools
import re
from AGRlib import getConfig, stripNulls, buildMetaObject, doQuery, makeOneOfConstraint

#-----------------------------------
# load config settings
cp = getConfig()

MOUSEMINE     = cp.get("DEFAULT","MOUSEMINEURL")
taxon         = cp.get("DEFAULT","TAXONID")
GLOBALTAXONID = cp.get("DEFAULT","GLOBALTAXONID")
GENELITURL    = cp.get("DEFAULT","GENELITURL")
MYGENEURL     = cp.get("DEFAULT","MYGENEURL")
SAMPLEIDS     = cp.get("DEFAULT","SAMPLEIDS").split()
MGD_OLD_PREFIX= cp.get("DEFAULT","MGD_OLD_PREFIX")

# Mapping from data provider name as stored in MGI to name as needed by AGR
# Cross references exported to the file are limited to those where the provider's name
# has an entry in this map.
dataProviders	= {}
for n in cp.options("dataProviders"):
    dataProviders[n] = cp.get("dataProviders", n)

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
#	- restricts which xrefs are exported
#	- translates provider name
#	- packs provider name and id into a object
#	- ensures uniqueness 
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
	    "assembly"		: loc['chromosomeLocation.assembly'],
	    "chromosome"	: chrom,
	    "startPosition"	: int(loc['chromosomeLocation.start']),
	    "endPosition"	: int(loc['chromosomeLocation.end']),
	    "strand"		: convertStrand(loc['chromosomeLocation.strand'])
	}]
    #elif obj.chromosome and obj.chromosome.primaryIdentifier != "UN":
    else:
        if chrom and chrom != 'UN':
	  return [{
	    "assembly"		: '',
	    "chromosome"	: chrom
	}]

def getJsonObj(obj):
      try:
	  return stripNulls({
	    "primaryId"		: obj["primaryIdentifier"],
	    "symbol"		: obj["symbol"],
	    "name"		: obj["name"],
	    "geneSynopsis"	: obj["description"],
	    "soTermId"		: obj["sequenceOntologyTerm.identifier"],
	    "taxonId"		: GLOBALTAXONID,
	    "synonyms"		: [ s for s in obj["synonyms"] if not isSecondaryId(s) and s != obj["symbol"] and s != obj["name"] ],
	    "secondaryIds"	: [ formatSecondary(s) for s in obj["synonyms"] if isSecondaryId(s) ],
	    "crossReferences"	: formatXrefs(obj),
	    "genomeLocations"	: formatGenomeLocation(obj.get('chromosome.primaryIdentifier', None), obj.get('location', [None])[0])
	  })
      except:
          sys.stderr.write('ERROR in getJsonObj. obj=' + str(obj) + '\n')
	  sys.exit(1)
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
    qs2 = []
    for label, q in qs:
	qiter = doQuery (q % qmods, MOUSEMINE)
	qiter = itertools.groupby(qiter, lambda x: x['primaryIdentifier'])
	qiter = itertools.imap(lambda x, y=label: (x[0], y, list(x[1])), qiter)
	qs2.append(qiter)

    print '{\n  "metaData": %s,\n  "data": [' % json.dumps(buildMetaObject(MOUSEMINE), indent=2)
    first=True
    for x in itertools.groupby(heapq.merge(*qs2), lambda x: x[0]):
      obj = { 'mgiid' : x[0] }
      for y in list(x[1]):
	if y[1] == 'gene':
	    # copy in all the basic gene attrs (symbol, name, etc)
	    obj.update(y[2][0])
	elif y[1] == 'synonyms':
	    # make a simple list of synonyms
	    obj['synonyms'] = map(lambda x: x['synonyms.value'], y[2])
	else:
	    obj[y[1]] = y[2]

      if not obj.get("primaryIdentifier"):
        continue
      if not first: print ",",
      print json.dumps(getJsonObj(obj), indent=2)
      first=False
    print ']\n}'

mouseGenes = '''
    <query
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

main(parseCmdLine())
