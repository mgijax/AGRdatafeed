#!/usr/bin/env python 
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
# ISSUES:
# - GeneSegments not included
# - Pseudogenes have no homology or protein associations
#

import sys
import json
import itertools
import time
import types


# The following library in not Python std lib.
# See: http://henry.precheur.org/projects/rfc3339 
#
from rfc3339 import rfc3339

# Returns the current date-time in RFC-3339 format
# Example: "2017-01-26T15:00:42-05:00"
#
def getTimeStamp():
    return rfc3339(time.time())

#
from intermine.webservice import Service
#mousemine = Service("http://www.mousemine.org/mousemine/service")
mousemine = Service("http://beta.mousemine.org/mousemine/service")

# for each xref provider name that we want to export, add an entry mapping our 
# name to the AGR standard.
#
dataProviders = {
  "Entrez Gene" : "NCBIGene",
  "Ensembl Gene Model" : "Ensembl",
}

# Base URL for MyGene wiki pages
myGeneBaseUrl = "https://en.wikipedia.org/wiki/"

# Constructs and returns the metaData (header) for the dump file.
#
def buildMetaObject(service):
    # get current date
    currentDate = getTimeStamp()

    # Get the MGI release date, which is embedded in the description field of the MGI DataSource obj
    # For example, "Mouse Genome Informatics [MGI 6.07 2017-01-24]"
    release = None
    query = service.new_query("DataSource")
    query.add_view("name", "description")
    query.add_constraint("name", "=", "MGI", code = "B")
    for r in query:
      i = r.description.find('[')
      release = r.description[i:].strip()[1:-1].strip()

    return {
    "dataProvider" : "MGI",
    "dateProduced" : currentDate,
    "release" : release
    }

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
	"chromosomeLocation.strand"

    )
    #
    query.add_sort_order("primaryIdentifier", "ASC")
    #
    query.add_constraint("organism.taxonId", "=", "10090", code = "A")
    query.add_constraint("primaryIdentifier", "CONTAINS", "MGI:", code = "B")
    if len(ids):
	query.add_constraint("primaryIdentifier", "ONE OF", ids, code = "C")
    #
    query.outerjoin("synonyms")
    query.outerjoin("crossReferences")
    query.outerjoin("chromosomeLocation")
    #
    return query

# Builds/returns the query for class Gene
#
def buildGeneQuery(service, ids):
    # start w/ the seq feature query, then add Gene-specific parts
    query = buildSequenceFeatureQuery(service, 'Gene', ids)
    query.add_view(
	"proteins.primaryAccession",
	"homologues.homologue.primaryIdentifier", "homologues.homologue.symbol",
	"homologues.homologue.crossReferences.identifier"
    )
    query.add_constraint("homologues.dataSets.name", "=", "Mouse/Human Orthologies from MGI", code = "D")
    query.add_constraint("homologues.homologue.crossReferences.source.name", "=", "MyGene", code = "E")
    query.add_constraint("sequenceOntologyTerm.identifier", "!=", "SO:0000902", code = "F")
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
# Returns a list of cross reference objects.
#
def formatXrefs(obj):
    xrefs = []
    for x in obj.crossReferences:
      dp = dataProviders.get(x.source.name, None)
      if dp:
        xrefs.append({"dataProvider":dp, "id":x.identifier})
    if hasattr(obj,"proteins"):
	for x in obj.proteins:
	    xrefs.append({"dataProvider":"UniProtKB", "id":x.primaryAccession})
    return xrefs

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
    return myGeneBaseUrl + obj.homologues[0].homologue.crossReferences[0].identifier


def formatGenomeLocation(loc):
    if loc:
	s = loc.strand
	s = "+" if s in ["+1","1"] else "-" if s == "-1" else "."
        return {
	    "assembly"		: "GRCm38",
	    "chromosome"	: loc.locatedOn.primaryIdentifier,
	    "startPosition"	: loc.start,
	    "endPosition"	: loc.end,
	    "strand"		: s
	}
    return None

# The AGR spec is to leave out attributes that would otherwise have a null value or an empty list value.
# The implementation constructs complete objects, then calls this function to strip out any nulls/empty lists.
#
def stripNulls(obj):
  for k,v in obj.items():
    if v is None or type(v) is types.ListType and len(v) == 0:
      del obj[k]
  return obj

# Here is the magic by which an object returned by the query is converted to an object
# conforming to the spec.
#
def getJsonObj(obj):
  return stripNulls({
    "primaryId"		: obj.primaryIdentifier,
    "symbol"		: obj.symbol,
    "name"		: obj.name,
    "geneSynopsis"	: obj.description,
    "geneSynopsisUrl"	: formatMyGeneLink(obj),
    "geneLiteratureUrl"	: "http://www.informatics.jax.org/reference/marker/%s?typeFilter=Literature" % obj.primaryIdentifier,
    "soTermId"		: obj.sequenceOntologyTerm.identifier,
    "taxonId"		: "10090",
    "synonyms"		: [ s.value for s in obj.synonyms if not isSecondaryId(s.value) ],
    "secondaryIds"	: [ s.value for s in obj.synonyms if isSecondaryId(s.value) ],
    "crossReferences"	: formatXrefs(obj),
    "genomeLocations"	: [ formatGenomeLocation(obj.chromosomeLocation) ]
  })

# Main prog. Build the query, run it, and output 
def main(ids):
  #
  query = itertools.chain(buildGeneQuery(mousemine, ids), buildPseudogeneQuery(mousemine, ids))
  jobj = {
    "metaData" : buildMetaObject(mousemine),
    "data" : [ getJsonObj(x) for x in query ]
  }
  print json.dumps(jobj, sort_keys=True, indent=2, separators=(',', ': ')),

#
main(sys.argv[1:])



