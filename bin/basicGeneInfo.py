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
#
# Author: Joel Richardson
#

# standard libs
import sys
import json
import itertools
import time
import types
import argparse

# nonstandard dependencies
from AGRlib import getConfig, stripNulls, buildMetaObject
from intermine.webservice import Service

#-----------------------------------
# load config settings
cp = getConfig()

MOUSEMINEURL  = cp.get("DEFAULT","MOUSEMINEURL")
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

#-----------------------------------
# MouseMine connection

mousemine = Service(MOUSEMINEURL)

#-----------------------------------
def buildExpressedGeneQuery(service):
    query = service.new_query('Gene')
    query.add_view("primaryIdentifier")
    query.add_constraint("expression.assayId", "IS NOT NULL", code="A")
    return query

#-----------------------------------
# Constructs and returns the core of the query, suitable for any SequenceFeature subclass.
#
def buildSequenceFeatureQuery(service, subclassName, ids):
    query = service.new_query(subclassName)
    #
    query.add_view(
	"primaryIdentifier", "symbol", "name", "description",
	"sequenceOntologyTerm.identifier", "synonyms.value",
	"crossReferences.source.name", "crossReferences.identifier",
	"chromosomeLocation.locatedOn.primaryIdentifier",
	"chromosomeLocation.start", "chromosomeLocation.end",
	"chromosomeLocation.strand",
	"chromosomeLocation.assembly",
	"chromosome.primaryIdentifier"

    )
    #
    query.add_sort_order("primaryIdentifier", "ASC")
    #
    query.add_constraint("organism.taxonId", "=", taxon, code = "A")
    query.add_constraint("dataSets.name", "=", "Mouse Gene Catalog from MGI", code = "B")
    if len(ids):
	query.add_constraint("primaryIdentifier", "ONE OF", ids, code = "C")
    #
    query.outerjoin("synonyms")
    query.outerjoin("crossReferences")
    query.outerjoin("chromosomeLocation")
    query.outerjoin("chromosome")
    #
    return query

# Builds/returns the query for class Gene
#
def buildGeneQuery(service, ids):
    # start w/ the seq feature query, then add Gene-specific parts
    query = buildSequenceFeatureQuery(service, 'Gene', ids)
    query.add_view(
	"proteins.uniprotAccession",
	"homologues.homologue.primaryIdentifier", "homologues.homologue.symbol",
	"homologues.homologue.crossReferences.identifier"
    )
    query.add_constraint("homologues.dataSets.name", "=", "Mouse/Human Orthologies from MGI", code = "D")
    query.add_constraint("homologues.homologue.crossReferences.source.name", "=", "MyGene", code = "E")
    query.add_constraint("sequenceOntologyTerm.identifier", "!=", "SO:0000902", code = "F") # no transgenes
    query.outerjoin("proteins")
    query.outerjoin("homologues")
    return query

# Builds/returns the query for class Pseudogene
#
def buildPseudogeneQuery(service, ids):
    return buildSequenceFeatureQuery(service, "Pseudogene", ids)

#
# Queries MouseMine for Panther ids, which are attached to Homologue objects.
# Returns a dict mapping MGI ids to Panther ids.
def cachePantherIds(service, ids):
  d = {}
  query = mousemine.new_query("Gene")
  query.add_view("primaryIdentifier", "homologues.crossReferences.identifier")
  query.add_constraint("organism.taxonId", "=", "10090", code = "A")
  query.add_constraint("homologues.dataSets.name", "=", "Panther data set", code = "B")
  query.add_constraint("homologues.type", "!=", "paralogue", code="C")
  if len(ids):
    query.add_constraint("primaryIdentifier", "ONE OF", ids, code = "D")

  for row in query.rows():
    d[row["primaryIdentifier"]] = row["homologues.crossReferences.identifier"]

  return d

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
    if not hasattr(obj, "homologues") or len(obj.homologues) != 1:
        return None
    xrs = obj.homologues[0].homologue.crossReferences
    if len(xrs) == 0:
        return None
    symbol = xrs[0].identifier
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
    for x in obj.crossReferences:
      dp = dataProviders.get(x.source.name, None)
      if dp:
        xrefs.add((dp, x.identifier))
    if hasattr(obj,"proteins"):
	for x in obj.proteins:
	    if x.uniprotAccession:
	        xrefs.add(("UniProtKB", x.uniprotAccession))
    if hasattr(obj,"pantherId") and obj.pantherId:
      xrefs.add(("PANTHER", obj.pantherId))
    xrefs = list(xrefs)
    xrefs.sort()
    # new xref format for 1.0.0.0. Includes 2 parts: the id, and a list of page-tags (see resourceDescriptors.yaml)
    xrs = [{"id": x[0]+":"+x[1]} for x in xrefs]
    # add xrefs to MGI pages for this gene
    pgs = ["gene","gene/references"]
    if obj.primaryIdentifier in expressed:
        pgs.append("gene/expression")
    xrs.append({"id": obj.primaryIdentifier, "pages":pgs })
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
def formatGenomeLocation(obj):
    locations = []
    if obj.chromosomeLocation:
	loc = obj.chromosomeLocation
        locations.append({
	    "assembly"		: loc.assembly,
	    "chromosome"	: loc.locatedOn.primaryIdentifier,
	    "startPosition"	: loc.start,
	    "endPosition"	: loc.end,
	    "strand"		: convertStrand(loc.strand)
	})
    elif obj.chromosome and obj.chromosome.primaryIdentifier != "UN":
        locations.append({
	    "assembly"		: '',
	    "chromosome"	: obj.chromosome.primaryIdentifier,
	})
    return locations

# Here is the magic by which an object returned by the query is converted to an object
# conforming to the spec.
#
def getJsonObj(obj):
  return stripNulls({
    "primaryId"		: obj.primaryIdentifier,
    "symbol"		: obj.symbol,
    "name"		: obj.name,
    "geneSynopsis"	: obj.description,
    "soTermId"		: obj.sequenceOntologyTerm.identifier,
    "taxonId"		: GLOBALTAXONID,
    "synonyms"		: [ s.value for s in obj.synonyms if not isSecondaryId(s.value) and s.value != obj.symbol and s.value != obj.name ],
    "secondaryIds"	: [ formatSecondary(s.value) for s in obj.synonyms if isSecondaryId(s.value) ],
    "crossReferences"	: formatXrefs(obj),
    "genomeLocations"	: formatGenomeLocation(obj)
  })

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
  
# Main prog. Build the query, run it, and output 
def main():
  args = parseCmdLine()
  ids = args.identifiers
  #
  global expressed
  expressed = set()
  for x in buildExpressedGeneQuery(mousemine).rows():
      expressed.add(x["primaryIdentifier"])
  #
  mgi2panther = cachePantherIds(mousemine, ids)
  def addPantherId(obj):
    obj.pantherId = mgi2panther.get(obj.primaryIdentifier, None)
    return obj
  #
  query = itertools.chain(buildGeneQuery(mousemine, ids), buildPseudogeneQuery(mousemine, ids))
  jobj = {
    "metaData" : buildMetaObject(mousemine),
    "data" : [ getJsonObj(addPantherId(x)) for x in query ]
  }
  print json.dumps(jobj, sort_keys=True, indent=2, separators=(',', ': ')),


#
main()
