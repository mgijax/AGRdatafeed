#
# variant.py
#
#
import sys
import re
import subprocess
import json
from AGRlib import getConfig, stripNulls, buildMetaObject, makeOneOfConstraint, doQuery, sql

cp = getConfig()
MOUSEMINE     = cp.get("DEFAULT","MOUSEMINEURL")
GLOBALTAXONID = cp.get("DEFAULT","GLOBALTAXONID")

chr2accid = {
  "GRCm38" : {
    "1" : "NC_000067.6",
    "2" : "NC_000068.7",
    "3" : "NC_000069.6",
    "4" : "NC_000070.6",
    "5" : "NC_000071.6",
    "6" : "NC_000072.6",
    "7" : "NC_000073.6",
    "8" : "NC_000074.6",
    "9" : "NC_000075.6",
    "10" : "NC_000076.6",
    "11" : "NC_000077.6",
    "12" : "NC_000078.6",
    "13" : "NC_000079.6",
    "14" : "NC_000080.6",
    "15" : "NC_000081.6",
    "16" : "NC_000082.6",
    "17" : "NC_000083.6",
    "18" : "NC_000084.6",
    "19" : "NC_000085.6",
    "X" : "NC_000086.7",
    "Y" : "NC_000087.7",
    "MT": "NC_005089.1"
  }
}

#
S_RE = re.compile('[^a-zA-Z-]')
def cleanseSequence (s) :
   return S_RE.sub('', s)
#
def getJsonObj(r) :
  vtype = r["type"]
  rr = stripNulls({
    "alleleId": r["allele_id"],
    "assembly" : r["build"],
    "chromosome" : r["chromosome"],
    "start" : int(r["startcoordinate"]),
    "end" : int(r["endcoordinate"]),
    "genomicReferenceSequence" : cleanseSequence(r["referencesequence"]),
    "genomicVariantSequence" : cleanseSequence(r["variantsequence"]),
    "type" : vtype,
    "consequence" : r["effect"],
    "sequenceOfReferenceAccessionNumber": "RefSeq:" + chr2accid[r["build"]][r["chromosome"]],
    "references" : [{"publicationId" : x} for x in r['refs']],
    "note" : r["note"],
    "crossReferences" : [{ "id" : r["allele_id"], "pages" : ["allele"] }]
  })
  #
  grs = rr["genomicReferenceSequence"]
  gvs = rr["genomicVariantSequence"]
  if vtype == "SO:0000159":
    # deletion
    # In MGI, the coding convention for deletions is that the genomic sequence contains the 
    # deleted part plus a padding base, while the variant sequence contains just the padding base.
    # In AGR, the genomic sequence contains just the deleted part, the variant sequence is 'N/A',
    # and the padding base is stored separately.
    if len(gvs) != 1:
      log("\nSkipping deletion because var seq len != 1: " + str(rr))
      return None
    if grs[0] != gvs:
      log("\nSkipping deletion because cannot compute padding base: " + str(rr))
      return None
    rr["start"] += 1
    rr["paddedBase"] = gvs
    rr["genomicVariantSequence"] = "N/A"
    rr["genomicReferenceSequence"] = rr["genomicReferenceSequence"][1:]
  elif vtype == "SO:0000667":
    # insertion
    # In MGI, the coding convention for insertions is that the variant sequence contains the 
    # inserted part plus a padding base, while the genomic sequence contains just the padding base.
    # In AGR, the variant sequence contains just the inserted part, the genomic sequence is 'N/A',
    # and the padding base is stored separately.
    if len(grs) != 1:
      log("\nSkipping insertion because genomic seq len != 1: " + str(rr))
      return None
    if grs != gvs[0]:
      log("\nSkipping insertion because cannot compute padding base: " + str(rr))
      return None
    rr["paddedBase"] = grs
    rr["genomicReferenceSequence"] = "N/A"
    rr["genomicVariantSequence"] = gvs[1:]
  elif vtype == "SO:0002007":
    # MNV
    if len(grs) != len(gvs):
      log("\nSkipping MNV because lengths are != : " + str(rr))
      return None
  elif vtype == "SO:1000032":
     # delins
     pass
  elif vtype == "SO:1000008":
    # point mutation
    if len(grs) != 1 or len(gvs) != 1:
      log("\nSkipping point mutation because of lengths. " + str(rr))
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
#
def main () :
    #
    vk2effects = {}
    for x in sql(Q_EFFECTS):
      vk2effects[x['_variant_key']] = x['accid']

    vk2types = {}
    for x in sql(Q_TYPES):
      vk2types[x['_variant_key']] = x['accid']

    vk2refs = {}
    for x in sql(Q_REFS):
      rid = ('PMID:' + x['pubmedid']) if x['pubmedid'] else x['mgiid']
      vk2refs.setdefault(x['_variant_key'], []).append(rid)

    vk2notes = {}
    for x in sql(Q_VARIANT_NOTES):
        vk2notes[x['_variant_key']] = { 'note' : x['note'] }

    n = 0
    print('{\n  "metaData": %s,\n  "data": [' % json.dumps(buildMetaObject(MOUSEMINE), indent=2))
    for x in sql(Q_VARIANTS):
      x['build'] = "GRCm38" # FIXME: should get this from the DB
      x['type'] = vk2types.get(x['_variant_key'], None)
      x['effect'] = vk2effects.get(x['_variant_key'], None)
      x['refs'] = vk2refs.get(x['_variant_key'], [])
      x['note'] = vk2notes.get(x['_variant_key'], None)
      try:
          j = getJsonObj(x)
          if not j: continue
      except:
          log("\nSkipping variant because of error: key=" +x['_variant_key'] + " " + str(x))
          log("Error=" + str(sys.exc_info()[1]))
          continue
      if n > 0: sys.stdout.write(",")
      try:
          sys.stdout.write(json.dumps(j, indent=2))
          n += 1
      except:
          log("\nSkipping variant because of encoding error: key=" +x['_variant_key'] + " " + str(j))
          log("Error=" + str(sys.exc_info()[1]))
    print("]}")

#
Q_VARIANTS = '''
  select
      v._variant_key,
      m.chromosome,
      aa.accid as allele_id,
      vs.startcoordinate,
      vs.endcoordinate,
      vs.referencesequence,
      vs.variantsequence
  from
      mrk_marker m,
      all_allele a,
      acc_accession aa,
      all_variant v,
      all_variant_sequence vs
  where m._marker_key = a._marker_key
      and a._allele_key = v._allele_key
      and a._allele_key = aa._object_key
      and aa._mgitype_key = 11
      and aa._logicaldb_key = 1
      and aa.preferred = 1
      and v._variant_key = vs._variant_key
      and v.isreviewed = 1
      and vs._sequence_type_key = 316347
  order by v._variant_key
  '''

#
Q_VARIANT_NOTES = '''
    select nc.note, n._object_key as _variant_key
    from mgi_notechunk nc, mgi_note n
    where nc._note_key = n._note_key
    and n._notetype_key = 1051
    '''
#
Q_TYPES = '''
  select
      v._variant_key,
      aa.accid
  from
      all_variant v,
      voc_annot va,
      voc_term vt,
      acc_accession aa
  where v._variant_key = va._object_key
      and va._annottype_key = 1026
      and va._term_key = vt._term_key
      and vt._term_key = aa._object_key
      and aa._mgitype_key = 13
      and aa.preferred = 1
  order by v._variant_key
  '''

#
Q_EFFECTS = '''
  select
      v._variant_key,
      aa.accid
  from
      all_variant v,
      voc_annot va,
      voc_term vt,
      acc_accession aa
  where v._variant_key = va._object_key
      and va._annottype_key = 1027
      and va._term_key = vt._term_key
      and vt._term_key = aa._object_key
      and aa._mgitype_key = 13
      and aa.preferred = 1
  order by v._variant_key
  '''

Q_REFS = '''
  select
      v._variant_key,
      aa.accid as mgiid,
      aa2.accid as pubmedid
  from
      all_variant v
      join mgi_reference_assoc ra
        on v._variant_key = ra._object_key
        and ra._mgitype_key = 45
      join acc_accession aa
        on ra._refs_key = aa._object_key
        and aa._mgitype_key = 1
        and aa._logicaldb_key = 1
        and aa.prefixpart = 'MGI:'
      left outer join acc_accession aa2
        on aa2._object_key = ra._refs_key
        and aa2._mgitype_key = 1
        and aa2._logicaldb_key = 29
  '''

#
main ()
#
