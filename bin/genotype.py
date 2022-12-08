#
# genotype.py
#
# Want all genotypes that have disease or pheno annotations.
#

import sys
import os
import re
import json
import itertools
import time
import types
import argparse

from AGRlib import stripNulls, buildMetaObject, makeOneOfConstraint, sql, makePubRef
from AGRqlib import qSubmittedAlleleIds, qSubmittedGenotypes, qGenotypeAllelePair

GLOBALTAXONID = os.environ["GLOBALTAXONID"]

# IDs of genotypes to omit (the "Not applicable" and the "Not specified" genotypes)
SKIP = ["MGI:2166309", "MGI:2166310" ]

#
def parseCmdLine():
    parser = argparse.ArgumentParser(description='Dumps genotype information to a JSON file.')

    parser.add_argument(
      'identifiers', 
      metavar='ids',
      nargs='*',
      help='Specific MGI ids to dump.')

    args = parser.parse_args()
    return args
  
bracketed_re = re.compile(r'<([^>]+)>')
def htmlify (s) :
  parts = bracketed_re.split(s)
  parts2 = []
  for i,part in enumerate(parts):
    if i % 2:
      part = '<sup>%s</sup>' % part
    parts2.append(part)
  return ''.join(parts2)

# Returns a JSON object for the given genotype
def getJsonObj (g, includedAlleles) :
    # Turn AllelePairs  into genotypeComponents.
    # Have to guard against genotypes having duplicate AllelePair records (it happens).
    pts = set([(c["alleleId"], mgi2geno.get(c["pairState"], None)) for c in g["components"]])
    #
    comps = [{
      "alleleID" : p[0],
      "zygosity" : p[1]
    } for p in pts ]
    #
    for c in comps:
      if c["zygosity"] is None:
        #return None
        pass
    #
    name = ("%s [background:] %s" % (g["alleles"], g["strain"])).replace("\n", " ")
    return {
      "primaryID" : g["genotypeId"],
      "subtype" : "genotype",
      "name" :  htmlify(name),
      #"nameText" : g["name"].strip(),
      "taxonId" : GLOBALTAXONID,
      "crossReference" : {
        "id" : g["genotypeId"],
        "pages" : [ "genotype" ]
      },
      "affectedGenomicModelComponents" :  comps
    }

# Main prog. Build the query, run it, and output 
def main():
  args = parseCmdLine()
  ids = args.identifiers
  xtra = makeOneOfConstraint('Genotype.primaryIdentifier', ids)

  # Process Genotype-AllelePairs. Build index from genotype id to list of component (allele+state)
  id2components = {}
  toDelete = set(SKIP)
  for r in sql(qGenotypeAllelePair):
    gid = r["genotypeId"]
    id2components.setdefault(gid, []).append(r)

  # Build set of MGI ids of alleles being sent to the alliance
  includedAlleles = set()
  for r in sql(qSubmittedAlleleIds):
    includedAlleles.add(r['mgiid'])

  # Process genotypes. For each one, find / attach its components if any and output.
  # Screen for genotypes to be deleted.
  #
  print('{\n  "metaData": %s,\n  "data": [' % json.dumps(buildMetaObject(), indent=2))
  first = True
  for g in sql(qSubmittedGenotypes):
    gid = g["genotypeId"]
    g["components"] = id2components.get(gid,[])
    gobj = getJsonObj(g, includedAlleles)
    if gobj:
        if not first: print(",", end=' ')
        print(json.dumps(gobj, indent=2))
        first=False
  print("]}")

## ----------------------------------------------------

#
# valid GENO term ids for Alliance submissions
# See: http://www.ontobee.org/ontology/GENO
#
genoZygosityTerms = [
("GENO:0000602", "homoplasmic"),
("GENO:0000603", "heteroplasmic"),
("GENO:0000604", "hemizygous X-linked"),
("GENO:0000605", "hemizygous Y-linked"),
("GENO:0000606", "hemizygous insertion-linked"),
("GENO:0000135", "heterozygous"),
("GENO:0000136", "homozygous"),
("GENO:0000137", "unspecified zygosity"),
]

mgi2geno = {
 "Homozygous" : "GENO:0000136", 
 "Heterozygous" : "GENO:0000135",
 "Hemizygous Insertion" : "GENO:0000606",
 "Indeterminate" : "GENO:0000137",
 "Hemizygous X-linked" : "GENO:0000604",
 "Hemizygous Y-linked" : "GENO:0000605",
 "Hemizygous Deletion" : "GENO:0000134", # no term for hemizygous deletion. Use plain hemizygous.
 "Homoplasmic" : "GENO:0000602",
 "Heteroplasmic" : "GENO:0000603"
}

main()
