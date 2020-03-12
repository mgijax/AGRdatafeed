#
# genotype.py
#
# Only want genotypes where any/all alleles are also in the alleleInfo file.
#

import sys
import re
import json
import itertools
import time
import types
import argparse

from AGRlib import getConfig, stripNulls, buildMetaObject, makeOneOfConstraint, doQuery, makePubRef

cp = getConfig()

MOUSEMINE     = cp.get("DEFAULT","MOUSEMINEURL")
GLOBALTAXONID = cp.get("DEFAULT","GLOBALTAXONID")

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
  try:
    # Turn AllelePairs  into genotypeComponents.
    # Have to guard against genotypes having duplicate AllelePair records (it happens).
    pts = set([(c["allelePairs.allele1.primaryIdentifier"],
            mgi2geno[c["allelePairs.pairState"]]) for c in g["components"]])
    #
    comps = [{
      "alleleID" : p[0],
      "zygosity" : p[1]
    } for p in pts ]
    #
    # Test that all component alleles are in the included set. 
    # Skip if that is not the case.
    for c in comps:
      if c["alleleID"] not in includedAlleles:
        return None
    #
    return {
      "primaryID" : g["primaryIdentifier"],
      "subtype" : "genotype",
      "name" :  htmlify(g["name"].strip()),
      #"nameText" : g["name"].strip(),
      "taxonId" : GLOBALTAXONID,
      "crossReference" : {
        "id" : g["primaryIdentifier"],
        "pages" : [ "genotype" ]
      },
      "affectedGenomicModelComponents" :  comps
    }
  except:
    return None

# Main prog. Build the query, run it, and output 
def main():
  args = parseCmdLine()
  ids = args.identifiers
  xtra = makeOneOfConstraint('Genotype.alleles.feature.primaryIdentifier', ids)
  xtra2 = makeOneOfConstraint('Allele.feature.primaryIdentifier', ids)

  # Process Genotype-AllelePairs. Build index from genotype id to list of component (allele+state)
  id2components = {}
  toDelete = set(SKIP)
  for r in doQuery(q_genotypeAlleles % xtra, MOUSEMINE):
    gid = r["primaryIdentifier"]
    try:
      id2components.setdefault(gid, []).append(r)
    except:
      toDelete.add(gid)

  # Build set of MGI ids of alleles being sent to the alliance
  includedAlleles = set()
  for r in doQuery(q_alleles % xtra2, MOUSEMINE):
    includedAlleles.add(r['primaryIdentifier'])

  # Process genotypes. For each one, find / attach its components if any and output.
  # Screen for genotypes to be deleted.
  #
  print '{\n  "metaData": %s,\n  "data": [' % json.dumps(buildMetaObject(MOUSEMINE), indent=2)
  first = True
  for g in doQuery(q_genotypes % xtra, MOUSEMINE):
    gid = g["primaryIdentifier"]
    if gid in toDelete:
      continue
    g["components"] = id2components.get(gid,[])
    gobj = getJsonObj(g, includedAlleles)
    if gobj:
	if not first: print ",",
	print json.dumps(gobj, indent=2)
	first=False
  print "]}"

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
 # "Hemizygous Deletion" : ??? - no geno term for this one (as of 18-03-2019)
 "Homoplasmic" : "GENO:0000602",
 "Heteroplasmic" : "GENO:0000603"
}

q_alleles = '''<query
  model="genomic"
  view="Allele.primaryIdentifier"
  constraintLogic="A and (B or E) and C and D"
  sortOrder="Allele.primaryIdentifier asc"
  >
  <constraint code="A" path="Allele.organism.taxonId" op="=" value="10090" />
  <constraint code="B" path="Allele.alleleType" op="NONE OF"><value>QTL</value></constraint>
  <constraint code="E" path="Allele.alleleType" op="IS NULL" />
  <constraint code="C" path="Allele.isWildType" op="=" value="false" />
  <constraint code="D" path="Allele.ontologyAnnotations" op="IS NOT NULL" />
  %s
</query>
'''
q_genotypes = '''<query
  model="genomic"
  view="
    Genotype.primaryIdentifier
    Genotype.name
    Genotype.background.primaryIdentifier
    "
  sortOrder="Genotype.primaryIdentifier asc"
  >
  %s
</query>
'''

q_genotypeAlleles = '''<query
  model="genomic"
  view="
    Genotype.primaryIdentifier
    Genotype.allelePairs.allele1.primaryIdentifier
    Genotype.allelePairs.allele1.symbol
    Genotype.allelePairs.pairState
    Genotype.allelePairs.allele2.primaryIdentifier
    Genotype.allelePairs.allele2.symbol
    Genotype.allelePairs.allele2.isWildType
    "
  sortOrder="Genotype.primaryIdentifier asc Genotype.allelePairs.allele1.primaryIdentifier asc"
  >
    <join path="Genotype.allelePairs.allele2" style="OUTER"/>
    %s
  </query>
'''

main()
