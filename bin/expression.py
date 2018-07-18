#!/usr/bin/env python2.7 
#
# expression.py
#
# Script to dump expression info in AGR json format.
# The format is described here:
#	https://github.com/alliance-genome/agr_schemas
#
# See: https://github.com/alliance-genome/agr_schemas/releases/tag/1.0.0.4
# Of particular note (from the doc):
#   4. an expression annotation has one gene, one where, one pub, one assay and one when 
#
# The data to be dumped include positve expression results in wild type specimens.
# NOTE: In situ reporter knockin hets are considered wild type.
#
# Assay types mapped to MMO terms.
#
# TODO: anatomy terms mapped to AGR grouping labels.
# See: https://docs.google.com/spreadsheets/d/1_UcKTq7y-wsQ83_kJlP6X5mQdgKykaPWlDdkZzCuygI/edit#gid=313336440
#
# Theiler stages mapped to:
# See https://docs.google.com/spreadsheets/d/1_UcKTq7y-wsQ83_kJlP6X5mQdgKykaPWlDdkZzCuygI/edit#gid=463437957
#
# Within a given assay, only report one result per anatomy/stage combo. (Ie, if the assay contains the
# same time/space result in multiple specimens, only output one row.)
#
# Author: Joel Richardson
#

# standard libs
import sys
import re
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

#-----------------------------------
# MouseMine connection

MOUSEMINEURL  = cp.get("DEFAULT","MOUSEMINEURL")
mousemine = Service(MOUSEMINEURL)

#-----------------------------------
def log(msg):
    sys.stderr.write(msg + '\n')

#-----------------------------------
# Returns a mapping from EMAPA id to EMAPA structure name
def loadEMAPA (service):
    id2emapa = {}
    query = service.new_query("EMAPATerm")
    query.add_view("identifier", "name")
    for t in query:
        id2emapa[t.identifier] = t.name
    return id2emapa

#-----------------------------------
# Returns unique-ified expression assay results. 
#
def getExpressionData(service,ids):
    query = service.new_query("GXDExpression")
    #
    query.add_view(
        "assayId", "assayType", "feature.primaryIdentifier",
	"stage", "structure.identifier", "publication.mgiJnum",
	"publication.pubMedId"
    )
    query.add_constraint("detected", "=", "true", code = "B")
    query.add_constraint("genotype.hasMutantAllele", "=", "false", code = "C")
    query.add_constraint("assayType", "=", "In situ reporter (knock in)", code = "D")
    query.add_constraint("genotype.zygosity", "=", "ht", code = "E")

    lexp = "B and (C or (D and E))"
    if len(ids):
	lexp = "A and " + lexp
	query.add_constraint("feature.primaryIdentifier", "ONE OF", ids, code = "A")
    query.set_logic(lexp)

    # Grrr. Because of a bug in the python IM client lib, sort specifications are ignored.
    # Have to read in the whole thing and sort in memory. Yuck.
    data = list(query.rows())
    data.sort(key=lambda r: (r["assayId"],r["structure.identifier"],r["stage"]))
    #
    prev = None
    unique = []
    for r in data:
	if not prev or r["assayId"] != prev["assayId"] or r["stage"] != prev["stage"] or r["structure.identifier"] != prev["structure.identifier"]:
	    unique.append(r)
	#
	prev = r
    #
    #log('getExpressionData: %d results => %d unique results' % (len(data), len(unique)))
    return unique

# Here is the magic by which an object returned by the query is converted to an object
# conforming to the spec.
#
def getJsonObj(obj, structureName):
  mkid = lambda i,p: None if i is None else p+i
  return stripNulls({
      'geneId': obj['feature.primaryIdentifier'],
      'evidence' : {
          'modPublicationId': obj['publication.mgiJnum'],
	  'pubMedId': mkid(obj['publication.pubMedId'], 'PMID:')
      },
      'whenExpressedStage': obj['stage'],
      'assay': obj['assayType'],
      'dateAssigned' : '??',
      'wildtypeExpressionTermIdentifiers' : {
          'anatomicalStructureTermId' : obj['structure.identifier'],
	  'whereExpressedStatement' : structureName
      },
      'crossReference' : {
          'id' : obj['assayId'],
	  'pages' : [ 'gene/expression/annotation' ]
      }
  })

#
def parseCmdLine():
    parser = argparse.ArgumentParser(description='Dumps basic allele information to a JSON file.')

    parser.add_argument(
      'identifiers', 
      metavar='ids',
      nargs='*',
      help='Specific MGI ids to dump.')

    args = parser.parse_args()
    return args
  
# Main prog. Build the query, run it, and output 
def main():
  args = parseCmdLine()
  ids = args.identifiers
  #
  id2emapa = loadEMAPA(mousemine)
  exprData = getExpressionData(mousemine, ids)
  print '{ "metaData" : %s, ' % json.dumps(buildMetaObject(mousemine))
  print '  "data" : ['
  for r in exprData:
      print json.dumps(getJsonObj(r, id2emapa[r['structure.identifier']]), sort_keys=True, indent=2, separators=(',', ': '))
  print ']}'

#
main()
