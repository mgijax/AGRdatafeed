#!/usr/bin/env python2.7 
#
# basicAlleleInfo.py
#
# Script to dump basic allele information from MGI in the AGR standard JSON format.
# The format is described here:
#	https://github.com/alliance-genome/agr_schemas
#
# Example allele IDs: MGI:3603004 MGI:2680557
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
from ConfigParser import ConfigParser

# nonstandard dependencies
from AGRlib import stripNulls, buildMetaObject
from intermine.webservice import Service

#-----------------------------------
# Load config
cp = ConfigParser()
cp.optionxform = str # make keys case sensitive
cp.read("config.cfg")

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
# Constructs and returns the core of the query, suitable for any SequenceFeature subclass.
#
def buildAlleleQuery(service,ids):
    query = service.new_query("Allele")
    query.add_constraint("feature", "Gene")
    query.add_view("primaryIdentifier", "symbol", "name", "feature.primaryIdentifier", "synonyms.value")
    query.add_constraint("organism.taxonId", "=", "10090", code = "A")
    query.add_constraint("alleleType", "NONE OF", ["QTL", "Transgenic"], code = "B")
    query.add_constraint("alleleType", "IS NULL", code = "E")
    query.add_constraint("isWildType", "=", "false", code = "C")
    query.add_constraint("ontologyAnnotations.ontologyTerm.id", "IS NOT NULL", code = "D")
    #
    query.outerjoin("synonyms")
    #
    query.add_sort_order("symbol", "ASC")
    #
    lexp = "A and (B or E) and C and D"
    if len(ids):
	query.add_constraint("primaryIdentifier", "ONE OF", ids, code = "F")
        lexp += " and F"
    query.set_logic(lexp)
    #
    return query

#
sup_re = re.compile(r'([<>])')
def insertSups (s) :
    def replace (s) :
        if s == "<":
            return "<sup>"
        elif s == ">":
            return "</sup>"
        else:
            return s
    return ''.join(map(replace, sup_re.split(s)))

#
def formatXrefs(obj):
    return [{"id":obj.primaryIdentifier, "pages":["allele"]}]

# Here is the magic by which an object returned by the query is converted to an object
# conforming to the spec.
#
def getJsonObj(obj):
  syns = list(set(insertSups(s.value) for s in obj.synonyms ))
  syns.sort()
  return stripNulls({
    "primaryId"		: obj.primaryIdentifier,
    "symbol"		: insertSups(obj.symbol),
    "taxonId"           : "NCBITaxon:10090",
    "gene"	        : obj.feature.primaryIdentifier,
    "synonyms"          : syns,
    "secondaryIds"      : [],
    "crossReferences"   : formatXrefs(obj)
  })

#
def parseCmdLine():
    parser = argparse.ArgumentParser(description='Dumps basic allele information to a JSON file.')

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
  query = buildAlleleQuery(mousemine, ids)
  jobj = {
    "metaData" : buildMetaObject(mousemine),
    "data" : [ getJsonObj(x) for x in query ]
  }
  print json.dumps(jobj, sort_keys=True, indent=2, separators=(',', ': ')),


#
main()
