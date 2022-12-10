#!/usr/bin/env python2.7 
#
# basicAlleleInfo.py
#
# Script to dump basic allele information from MGI in the AGR standard JSON format.
# The format is described here:
#       https://github.com/alliance-genome/agr_schemas
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
import os
import argparse

#
from AGRlib import stripNulls, buildMetaObject, sql
from AGRqlib import qAlleles, qAlleleSynonyms, qExpressors

#-----------------------------------
GLOBALTAXONID = os.environ["GLOBALTAXONID"]

#-----------------------------------
#
def getAlleles():
    # Query for alleles that have expressed component 
    expressors = set()
    for r in sql(qExpressors):
        expressors.add(r['_allele_key'])

    # Query allele synonyms, build index of id -> synonyms
    aid2syns = {}
    for r in sql(qAlleleSynonyms):
      aid2syns.setdefault(r['_allele_key'], set()).add(r['synonym'])

    # Main allele query.
    for r in sql(qAlleles):
      ak = r['_allele_key']
      aid = r['alleleId']
      r['synonyms'] = list(aid2syns.get(ak, []))
      # If the allele has a driver or has expressed components, then the allele has a
      # "construct". At the Alliance, constructs must have an ID, but at MGI they don't (they're not objects).
      # So we create a fake ID for it. These are not displayed and are not used to create links.
      # Constructs are dumped separately. Here we just need the ID
      if r['drivenBy'] or ak in expressors:
          r['construct'] = aid + "_con"
      #
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
    return [{"id":obj["alleleId"], "pages":["allele", "allele/references"]}]

# Here is the magic by which an object returned by the query is converted to an object
# conforming to the spec.
#
def getJsonObj(obj):
  syns = list(map(insertSups, obj["synonyms"]))
  syns.sort()
  ###
  isTgAllele = obj["alleleType"] == "Transgenic"
  isTgMarker = obj["markerType"] == "transgene"
  ###
  gene = None if isTgAllele and isTgMarker else obj["markerId"]

  aors = [] # alleleObjectRelations
  if gene:
      aors.append({"objectRelation": {"associationType":"allele_of","gene":gene }})
  if 'construct' in obj:
      aors.append({"objectRelation": {"associationType":"contains","construct":obj['construct'] }})
  #
  return stripNulls({
    "primaryId"         : obj["alleleId"],
    "symbol"            : insertSups(obj["symbol"]),
    "symbolText"        : obj["symbol"],
    "taxonId"           : GLOBALTAXONID,
    "synonyms"          : syns,
    "secondaryIds"      : [],
    "crossReferences"   : formatXrefs(obj),
    "alleleObjectRelations": aors,
    "description" : obj["molecularNote"]
  })

# Main prog. Build the query, run it, and output 
def main():
  #
  query = getAlleles()
  print('{\n  "metaData": %s,\n  "data": [' % json.dumps(buildMetaObject(), indent=2))
  first = True
  for a in query:
    if not first: print(",", end=' ')
    print(json.dumps(getJsonObj(a), indent=2))
    first=False
  print("]}")

#
main()
