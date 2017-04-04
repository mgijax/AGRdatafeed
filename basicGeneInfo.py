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
from ConfigParser import ConfigParser

from AGRlib import stripNulls, buildMetaObject

# nonstandard dependencies

# See: http://henry.precheur.org/projects/rfc3339 
from rfc3339 import rfc3339
from intermine.webservice import Service


#-----------------------------------
# Load config
cp = ConfigParser()
cp.optionxform = str # make keys case sensitive
cp.read("config.cfg")

MOUSEMINEURL	= cp.get("DEFAULT","MOUSEMINEURL")
TAXONID		= cp.get("DEFAULT","TAXONID")
GENELITURL	= cp.get("DEFAULT","GENELITURL")
MYGENEURL	= cp.get("DEFAULT","MYGENEURL")
SAMPLEIDS	= cp.get("DEFAULT","SAMPLEIDS").split()
PREFIX		= "MGI:"

# Mapping from data provider name as stored in MGI to name as needed by AGR
# Cross references exported to the file are limited to those where the provider's name
# has an entry in this map.
dataProviders	= {}
for n in cp.options("dataProviders"):
    dataProviders[n] = cp.get("dataProviders", n)

#-----------------------------------
# RFC 3339 timestamps

# Returns the current date-time in RFC-3339 format
# Example: "2017-01-26T15:00:42-05:00"
#
def getTimeStamp():
    return rfc3339(time.time())

#-----------------------------------
# MouseMine connection

mousemine = Service(MOUSEMINEURL)

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
    query.add_constraint("organism.taxonId", "=", TAXONID, code = "A")
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

# Build/returns the query for clas Pseudogene
#
def buildPseudogeneQuery(service, ids):
    return buildSequenceFeatureQuery(service, "Pseudogene", ids)

# In MouseMine, synonyms and secondary ids are lumped together as "synonyms". 
# This function distinguishes a synonym value as being either a secondary id or not.
#
def isSecondaryId(identifier):
    return identifier.startswith("MGI:") or identifier.startswith("MGD-")

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
    xrefs = list(xrefs)
    xrefs.sort()
    return [{"dataProvider":x[0],"id":x[1]} for x in xrefs]

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
    return MYGENEURL + obj.homologues[0].homologue.crossReferences[0].identifier


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
    "primaryId"		: PREFIX + obj.primaryIdentifier,
    "symbol"		: obj.symbol,
    "name"		: obj.name,
    "geneSynopsis"	: obj.description,
    "geneSynopsisUrl"	: formatMyGeneLink(obj),
    "geneLiteratureUrl"	: GENELITURL % obj.primaryIdentifier,
    "soTermId"		: obj.sequenceOntologyTerm.identifier,
    "taxonId"		: TAXONID,
    "synonyms"		: [ s.value for s in obj.synonyms if not isSecondaryId(s.value) ],
    "secondaryIds"	: [ PREFIX + s.value for s in obj.synonyms if isSecondaryId(s.value) ],
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
  query = itertools.chain(buildGeneQuery(mousemine, ids), buildPseudogeneQuery(mousemine, ids))
  jobj = {
    "metaData" : buildMetaObject(mousemine),
    "data" : [ getJsonObj(x) for x in query ]
  }
  print json.dumps(jobj, sort_keys=True, indent=2, separators=(',', ': ')),


#
main()
