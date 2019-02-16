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

# nonstandard dependencies
from AGRlib import getConfig, stripNulls, buildMetaObject, makeOneOfConstraint, doQuery

#-----------------------------------
# load config settings
cp = getConfig()

MOUSEMINE     = cp.get("DEFAULT","MOUSEMINEURL")
taxon         = cp.get("DEFAULT","TAXONID")
GLOBALTAXONID = cp.get("DEFAULT","GLOBALTAXONID")
GENELITURL    = cp.get("DEFAULT","GENELITURL")
MYGENEURL     = cp.get("DEFAULT","MYGENEURL")
SAMPLEIDS     = cp.get("DEFAULT","SAMPLEALLELEIDS").split()
MGD_OLD_PREFIX= cp.get("DEFAULT","MGD_OLD_PREFIX")

# Mapping from data provider name as stored in MGI to name as needed by AGR
# Cross references exported to the file are limited to those where the provider's name
# has an entry in this map.
dataProviders	= {}
for n in cp.options("dataProviders"):
    dataProviders[n] = cp.get("dataProviders", n)

#-----------------------------------
# Constructs and returns the core of the query, suitable for any SequenceFeature subclass.
#
def buildAlleleQuery(url,ids):
    qopts = {
      'xtraConstraint': makeOneOfConstraint('Allele.feature.primaryIdentifier', ids)
    }
    aid2syns = {}
    qsynonyms = '''<query
      model="genomic"
      view="
      Allele.primaryIdentifier
      Allele.synonyms.value
      "
      constraintLogic="A and (B or E) and C and D"
      sortOrder="Allele.primaryIdentifier asc"
      >
      <constraint code="A" path="Allele.organism.taxonId" op="=" value="10090" />
      <constraint code="B" path="Allele.alleleType" op="NONE OF"><value>QTL</value><value>Transgenic</value></constraint>
      <constraint code="E" path="Allele.alleleType" op="IS NULL" />
      <constraint code="C" path="Allele.isWildType" op="=" value="false" />
      <constraint code="D" path="Allele.ontologyAnnotations.ontologyTerm.id" op="IS NOT NULL" />
      %(xtraConstraint)s
      </query>
    ''' % qopts
    for r in doQuery(qsynonyms, url):
      aid2syns.setdefault(r['primaryIdentifier'], set()).add(r['synonyms.value'])

    qalleles = '''<query
      model="genomic"
      view="
      Allele.primaryIdentifier
      Allele.symbol
      Allele.name
      Allele.feature.primaryIdentifier
      "
      constraintLogic="A and (B or E) and C and D"
      sortOrder="Allele.primaryIdentifier asc"
      >
      <constraint code="A" path="Allele.organism.taxonId" op="=" value="10090" />
      <constraint code="B" path="Allele.alleleType" op="NONE OF"><value>QTL</value><value>Transgenic</value></constraint>
      <constraint code="E" path="Allele.alleleType" op="IS NULL" />
      <constraint code="C" path="Allele.isWildType" op="=" value="false" />
      <constraint code="D" path="Allele.ontologyAnnotations.ontologyTerm.id" op="IS NOT NULL" />
      %(xtraConstraint)s
      </query>
    ''' % qopts
    
    for r in doQuery(qalleles, url):
      r['synonyms'] = list(aid2syns.get(r['primaryIdentifier'],[]))
      yield r

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
    return [{"id":obj["primaryIdentifier"], "pages":["allele"]}]

# Here is the magic by which an object returned by the query is converted to an object
# conforming to the spec.
#
def getJsonObj(obj):
  ###
  # MouseMine is loading things as allele synonyms that should not be included. 
  # Until this is fixed, weed out the unwanted ones here. Basically, the allele's
  # symbol and name are both being included and need to be screen out.
  # Note that in mousemine, allele.name includes the gene's name followed by a semicolon followed
  # by the allele's name.
  nn = obj["name"].rsplit(';')[-1].strip()
  syns = set()
  for s in obj["synonyms"]:
      if s != obj["symbol"] and s != nn:
          syns.add(s)
  syns = map(insertSups, list(syns))
  syns.sort()
  ###
  return stripNulls({
    "primaryId"		: obj["primaryIdentifier"],
    "symbol"		: insertSups(obj["symbol"]),
    "symbolText"	: obj["symbol"],
    "taxonId"           : "NCBITaxon:10090",
    "gene"	        : obj["feature.primaryIdentifier"],
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
  query = buildAlleleQuery(MOUSEMINE, ids)
  print '{\n  "metaData": %s,\n  "data": [' % json.dumps(buildMetaObject(MOUSEMINE), indent=2)
  first = True
  for a in query:
    if not first: print ",",
    print json.dumps(getJsonObj(a), indent=2)
    first=False
  print "]}"

#
main()
