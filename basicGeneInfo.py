!/usr/bin/env python

import sys
import json
import time
from intermine.webservice import Service

# for each xref provider name that we want to export, add an entry mapping our
# name to the AGR standard.
#
dataProviders = {
  "NCBI Gene Model" : "NCBI Gene",
  "Ensembl Gene Model" : "Ensembl",
  "VEGA Gene Model" : "Vega"
}

# Base URL for MyGene wiki pages
myGeneBaseUrl = "https://en.wikipedia.org/wiki/"

#
def buildMetaObject(service):
    # get current date
    currentDate = time.asctime(time.localtime(time.time()))

    # Get the MGI release date, which is embedded in the description field of the MGI DataSource obj
    # For example, "Mouse Genome Informatics [MGI 6.07 2017-01-24]"
    release = None
    query = service.new_query("DataSource")
    query.add_view("name", "description")
    query.add_constraint("name", "=", "MGI", code = "B")
    for r in query:
      release = r.description
      break
      i = r.description.find('[')
      release = r.description[i:].strip()[1:-1].strip()

    return {
    "dataProvider" : "MGI",
    "dateProduced" : currentDate,
    "release" : release
    }
#
def buildQuery(service, ids):
    query = service.new_query("Gene")
    #
    query.add_view(
	"primaryIdentifier", "symbol", "name", "description",
	"sequenceOntologyTerm.identifier", "synonyms.value",
	"crossReferences.source.name", "crossReferences.identifier",
	"proteins.primaryAccession",
	"homologues.homologue.primaryIdentifier", "homologues.homologue.symbol",
	"homologues.homologue.crossReferences.identifier"
    )
    #
    query.add_sort_order("Gene.primaryIdentifier", "ASC")
    #
    query.add_constraint("organism.taxonId", "=", "10090", code = "A")
    query.add_constraint("primaryIdentifier", "CONTAINS", "MGI:", code = "B")
    query.add_constraint("homologues.dataSets.name", "=", "Mouse/Human Orthologies from MGI", code = "C")
    query.add_constraint("homologues.homologue.crossReferences.source.name", "=", "MyGene", code = "D")
    if len(ids):
	query.add_constraint("primaryIdentifier", "ONE OF", ids, code = "E")
    #
    query.outerjoin("synonyms")
    query.outerjoin("crossReferences")
    query.outerjoin("proteins")
    query.outerjoin("homologues")
    #
    return query

# In MouseMine, synonyms and secondary ids are lumped together as "synonyms".
# This function distinguishes a synonym value as being either a secondary id or not.
def isSecondary(identifier):
    return identifier.startswith("MGI:") or identifier.startswith("MGD-")

# Selects the xrefs to be exported for the object and formats them according to the spec.
#	- translates provider name
#	- packs provider name and id into a object
# Returns a list of cross reference objects.
def formatXrefs(obj):
    xrefs = []
    for x in obj.crossReferences:
      dp = dataProviders.get(x.source.name, None)
      if dp:
        xrefs.append({"dataProvider":dp, "id":x.identifier})
    for x in obj.proteins:
      xrefs.append({"dataProvider":"UniProtKB", "id":x.primaryAccession})
    return xrefs

# In the MGI fewi, mouse genes link to a MyGenes wiki page which is a human readable description.
# The MyGenes page is for the HUMAN ortholog of the mouse gene.
# In the database, human genes have a cross reference to MyGenes, where the "id" is the part needed
# to fill out the complete URL (the base part if constant).
# The link is only displayed if the mouse/human gene have a 1:1 orthology relationship.
# Here we use the same logic to construct a link (or not) for the given mouse gene (obj).
def formatMyGeneLink(obj):
    if len(obj.homologues) != 1:
        return None
    return myGeneBaseUrl + obj.homologues[0].homologue.crossReferences[0].identifier

def getJsonObj(obj):
  jobj = {
    "primaryId"	: obj.primaryIdentifier,
    "symbol"	: obj.symbol,
    "name"		: obj.name,
    "geneSynopsis"	: obj.description,
    "geneSynopsisUrl": formatMyGeneLink(obj),
    "soTermId"	: obj.sequenceOntologyTerm.identifier,
    "taxonId"	: "10090",
    "synonyms"	: [ s.value for s in obj.synonyms if not isSecondary(s.value) ],
    "secondaryIds"	: [ s.value for s in obj.synonyms if isSecondary(s.value) ],
    "crossReferences": formatXrefs(obj)
  }
  return jobj

# Main prog. Build the query, run it, and output
def main(ids):
    service = Service("http://www.mousemine.org/mousemine/service")

    print "{"

    #
    meta = buildMetaObject(service)
    print '"meta":', json.dumps(meta, indent = 2, separators=(',', ': ')), ","

    #
    print '"data": [',
    query = buildQuery(service, ids)
    for i,obj in enumerate(query):
      jobj = getJsonObj(obj)
      if i > 0:
          print ",",
      print json.dumps(jobj, indent=2, separators=(',', ': ')),

    #
    print "]"
    print "}"

main(sys.argv[1:])
