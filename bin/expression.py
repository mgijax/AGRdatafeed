#!/usr/bin/env python2.7 
#
# Expression.py
#
# Script to dump expression info in AGR json format.
# The format is described here:
#	https://github.com/alliance-genome/agr_schemas
#
# The data to be dumped include positve expression results in wild type specimens.
# (NOTE: In situ reporter knockin hets are considered wild type.)
#
# Assay types mapped to MMO terms.
#
# TODO: anatomy terms mapped to AGR grouping labels.
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

MOUSEMINEURL  = cp.get("DEFAULT","MOUSEMINEURL")
taxon         = cp.get("DEFAULT","TAXONID")
SAMPLEIDS     = cp.get("DEFAULT","SAMPLEALLELEIDS").split()

#-----------------------------------
# MouseMine connection

mousemine = Service(MOUSEMINEURL)

#-----------------------------------
# Constructs and returns the core of the query, suitable for any SequenceFeature subclass.
#
def buildExpressionQuery(service,ids):
    query = service.new_query("GXDExpression")
    #
    query.add_view(
        "assayId", "assayType", "feature.primaryIdentifier", "feature.symbol",
	"stage", "structure.identifier", "structure.name", "publication.mgiJnum",
	"publication.pubMedId", "genotype.hasMutantAllele", "genotype.zygosity",
	"genotype.symbol"
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

    #
    query.add_sort_order("GXDExpression.assayId", "ASC")

    return query

# Here is the magic by which an object returned by the query is converted to an object
# conforming to the spec.
#
def getJsonObj(obj):
  return stripNulls({
      'geneId': obj.feature.primaryIdentifier,
      'evidence' : {
          'modPublicationId': obj.publication.mgiJnum,
	  'pubMedId': obj.publication.pubMedId
      },
      'whenExpressedStage': obj.stage,
      'assay': obj.assayType,
      'dateAssigned' : '??',
      'wildtypeExpressionTermIdentifiers' : {
          'anatomicalStructureTermId' : obj.structure.identifier,
	  'whereExpressedStatement' : obj.structure.name
      },
      'crossReference' : {
          'id' : obj.assayId,
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
  query = buildExpressionQuery(mousemine, ids)
  #jobj = {
    #"metaData" : buildMetaObject(mousemine),
    #"data" : [ getJsonObj(x) for x in query ]
  #}

  print '{ "metaData" : %s, ' % json.dumps(buildMetaObject(mousemine))
  print '  "data" : ['
  for obj in query:
      print json.dumps(getJsonObj(obj), sort_keys=True, indent=2, separators=(',', ': '))
  print ']}'

#
main()
