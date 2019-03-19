#
# variant.py
#
#
import sys
import json
from AGRlib import getConfig, stripNulls, buildMetaObject, makeOneOfConstraint, doQuery

cp = getConfig()
MOUSEMINE     = cp.get("DEFAULT","MOUSEMINEURL")
GLOBALTAXONID = cp.get("DEFAULT","GLOBALTAXONID")

chr2accid = {
  "GRCm38" : {
    "chr1" : "NC_000067.6",
    "chr2" : "NC_000068.7",
    "chr3" : "NC_000069.6",
    "chr4" : "NC_000070.6",
    "chr5" : "NC_000071.6",
    "chr6" : "NC_000072.6",
    "chr7" : "NC_000073.6",
    "chr8" : "NC_000074.6",
    "chr9" : "NC_000075.6",
    "chr10" : "NC_000076.6",
    "chr11" : "NC_000077.6",
    "chr12" : "NC_000078.6",
    "chr13" : "NC_000079.6",
    "chr14" : "NC_000080.6",
    "chr15" : "NC_000081.6",
    "chr16" : "NC_000082.6",
    "chr17" : "NC_000083.6",
    "chr18" : "NC_000084.6",
    "chr19" : "NC_000085.6",
    "chrX" : "NC_000086.7",
    "chrY" : "NC_000087.7",
  }
}

#
def getJsonObj(r) :
  vtype = r["varTypeId"]
  rr = {
    "alleleId": r["alleleId"],
    "assembly" : r["build"],
    "chromosome" : r["chr"],
    "start" : int(r["start"]),
    "end" : int(r["end"]),
    "genomicReferenceSequence" : r["refAllele"],
    "genomicVariantSequence" : r["varAllele"],
    "type" : vtype,
    "consequence" : r["varEffectId"],
    "sequenceOfReferenceAccessionNumber": "RefSeq:" + chr2accid[r["build"]][r["chr"]],
    "references" : [{ "publicationId" : r["pubmedId"] }]
  }
  #
  if vtype ==   "SO:0000159":
    #deletion
    if len(rr["genomicVariantSequence"]) != 1:
      log("Skipping deletion because var seq len != 1: " + str(rr))
      return None
    if rr["genomicReferenceSequence"][0] != rr["genomicVariantSequence"]:
      log("Skipping deletion because cannot compute padding base: " + str(rr))
      return None
    rr["paddedBase"] = rr["genomicVariantSequence"]
    rr["genomicVariantSequence"] = "N/A"
    rr["genomicReferenceSequence"] = rr["genomicReferenceSequence"][1:]
  elif vtype == "SO:0002007":
    # MNV
    pass
  elif vtype == "SO:1000008":
    # point mutation
    if len(rr["genomicReferenceSequence"]) != 1 or len(rr["genomicVariantSequence"]) != 1:
      log("Skipping point mutation because of lengths. " + str(rr))
      return None
  else:
    # error
    raise RuntimeError("Unknown or unhandled vtype: " + vtype)
  return rr

#
def log (s) :
  sys.stderr.write(s)
  sys.stderr.write('\n')

#
def main () :
  colNames = []
  i = 0
  n = 0
  print '{\n  "metaData": %s,\n  "data": [' % json.dumps(buildMetaObject(MOUSEMINE), indent=2)
  for line in sys.stdin:
    row = line[:-1].split('\t')
    if i == 0:
      colNames = row[:]
    else:
      r = dict([ (colNames[j], row[j]) for j in range(len(colNames)) ])
      rr = getJsonObj(r)
      if rr:
	  if n: print ",",
	  print json.dumps(rr, indent=2)
	  n += 1
    i += 1
  print "]}"
#
main ()
