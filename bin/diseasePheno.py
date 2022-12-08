#!/usr/bin/env python2.7

# Create JSON file describing genotype-disease annotations for AGR data ingest.
# Usage:
#       python diseasePheno.py -d > MGI_0.6.2_diseaseAnnotations.json
# Writes JSON to stdout.
# See JSON specs at:
#    https://github.com/alliance-genome/agr_schemas/blob/development/disease/diseaseModelAnnotation.json
#    https://github.com/alliance-genome/agr_schemas/blob/release-1.0.0.3/phenotype/phenotypeModelAnnotation.json
# For debugging, accepts optional MGI IDs for genotypes/alleles/genes on command line
#   to restrict output to DAF records for those genotypes.
# Example genotype ID:  MGI:5526095 MGI:2175208 MGI:5588576
# Example gene IDs: MGI:97490 MGI:99607
# Example allele IDs: MGI:3603004 MGI:2680557
#       python diseasePheno.py -d MGI:5526095 MGI:2175208 MGI:5588576 MGI:99607 MGI:97490 MGI:97490 MGI:99607 > disease.sample.json
#       python diseasePheno.py -p MGI:5526095 MGI:2175208 MGI:5588576 MGI:99607 MGI:97490 MGI:97490 MGI:99607 > pheno.sample.json
#
# Original author: Jim Kadin
# Revisions: Joel Richardson
#

import sys
import os
import argparse
import json
from itertools import chain
from AGRlib import stripNulls, buildMgiDataProviderObject, buildMetaObject, getTimeStamp, sql, makePubRef
from AGRqlib import tAnnots, tAnnotEvidence, tAnnotBaseAnnots, tGenotypeLabels

mpGeneCfg = {
    "okind" : "MP",
    "skind" : "gene",
    "_annottype_key" : 1015, # MP-Gene
    "_baseannottype_key" : 1002, # MP-Genotype
    "subj_tblname"   : "MRK_Marker",
    "subj_keycol"    : "_marker_key",
    "subj_labelcol"  : "symbol",
    "_mgitype_key"   : 2,
    "voc_ldbkey"    : 34,
    }
mpAlleleCfg = {
    "okind" : "MP",
    "skind" : "allele",
    "_annottype_key" : 1028, # MP-Allele
    "_baseannottype_key" : 1002, # MP-Genotype
    "subj_tblname"   : "ALL_Allele",
    "subj_keycol"    : "_allele_key",
    "subj_labelcol"  : "symbol",
    "_mgitype_key"   : 11,
    "voc_ldbkey"    : 34,
    }
mpGenoCfg = {
    "okind" : "MP",
    "skind" : "genotype",
    "_annottype_key" : 1002, # MP-Genotype
    "_baseannottype_key" : None,
    "subj_tblname"   : "GXD_Genotype",
    "subj_keycol"    : "_genotype_key",
    "subj_labelcol"  : "_genotype_key", # special handling
    "_mgitype_key"   : 12, # genotype
    "voc_ldbkey"    : 34, # MP
    }
doGeneCfg = {
    "okind" : "DO",
    "skind" : "gene",
    "_annottype_key" : 1023, # DO-Gene
    "_baseannottype_key" : 1020, # DO-Genotype
    "subj_tblname"   : "MRK_Marker",
    "subj_keycol"    : "_marker_key",
    "subj_labelcol"  : "symbol",
    "_mgitype_key"   : 2,  # marker
    "voc_ldbkey"    : 191, # DO
    }
doAlleleCfg = {
    "okind" : "DO",
    "skind" : "allele",
    "_annottype_key" : 1029, # DO-Allele (derived)
    "_baseannottype_key" : 1020, # DO-Genotype
    "subj_tblname"   : "ALL_Allele",
    "subj_keycol"    : "_allele_key",
    "subj_labelcol"  : "symbol",
    "_mgitype_key"   : 11, # allele
    "voc_ldbkey"    : 191, # DO
    }
doAlleleDirectCfg = {
    "okind" : "DO",
    "skind" : "allele",
    "_annottype_key" : 1021, # DO-Allele (direct)
    "_baseannottype_key" : None,
    "subj_tblname"   : "ALL_Allele",
    "subj_keycol"    : "_allele_key",
    "subj_labelcol"  : "symbol",
    "_mgitype_key"   : 11, # allele
    "voc_ldbkey"    : 191, # DO
    }
doGenoCfg = {
    "okind" : "DO",
    "skind" : "genotype",
    "_annottype_key" : 1020, # DO-Genotype
    "_baseannottype_key" : None,
    "subj_tblname"   : "GXD_Genotype",
    "subj_keycol"    : "_genotype_key",
    "subj_labelcol"  : "_genotype_key", # special handling
    "_mgitype_key"   : 12, # genotype
    "voc_ldbkey"    : 191, # DO
    }

PHENO_CFGS = [mpGeneCfg, mpAlleleCfg, mpGenoCfg]
DISEASE_CFGS = [doGeneCfg, doAlleleCfg, doAlleleDirectCfg, doGenoCfg]

#
code2eco = {
  "TAS" : "ECO:0000033"
}

def getAnnotations (cfg) :
    LIMIT = "" # " limit 10"
    
    gk2label = {}
    if (cfg["skind"] == "genotype") :
        for r in sql(tGenotypeLabels % cfg):
            label = "%s [background:] %s" % (r["alleles"],r["strain"])
            gk2label[r["_genotype_key"]] = label

    ak2evs = {}
    for r in sql(tAnnotEvidence % cfg + LIMIT):
        ak2evs.setdefault(r["_annot_key"],[]).append(r)
   
    ak2bas = {}
    if (cfg["_baseannottype_key"]):
        for r in sql(tAnnotBaseAnnots % cfg + LIMIT):
            ak2bas.setdefault(r["_annot_key"],[]).append(r)

    for r in sql(tAnnots % cfg + LIMIT):
        r = dict(r)
        r["evidence"] = ak2evs.get(r["_annot_key"],[])
        r["baseAnnots"] = ak2bas.get(r["_annot_key"],[])
        rr = applyConversions(r, cfg["okind"], cfg["skind"])
        if rr:
            if (cfg["skind"] == "genotype") :
                rr["subjectLabel"] = gk2label[r["subjectKey"]]
            for e in rr["evidence"]:
                rr2 = formatDafJsonRecord(rr, e, "disease" if cfg["okind"] == "DO" else "phenotype", cfg["skind"])
                yield rr2
      
########
# Applies various transformations to a 'raw' annotation returned from the db to prepare it
# for export to AGR. Example: AGR dates must be in rfc3339 format. See comments for specific 
# transforms.
# Args:
#   a - one annotation
#   okind - kind of annotation ontology term, "DO" or "MP"
#   skind - kind of annotation subject, "gene", "allele", or "genotype"
# Returns:
#   The annotation, with conversions applied
#
geno2genes = {}   # genotype->rolled up genes index
def applyConversions(a, okind, skind):
    global geno2genes

    # screen out "normal" phenotype annotations
    # but allow "not" disease annotations. Not sure why the alliance 
    # handles one and not the other.
    if okind == "MP" and a["qualifier"] == "normal":
        return None
    elif okind == "DO":
        a["qualifier"] = "not" if a["qualifier"] == "NOT" else None
    ##
    # 1. Attach base annotations to their proper evidence records
    # 2. Set evidence annot date to earliest base annotation's date
    # 3. Collect base annotataion genotype ids
    for e in a["evidence"]:
        #
        if okind == "DO":
            e["codes"] = [code2eco[e["code"]]]
        else:
            e["codes"] = [e["code"]]
        # Attach base annotations to the inverted evidence records they belong to
        e["baseAnnots"] = []
        e["baseSubjectIds"] = []
        dt = None
        for ba in a.get("baseAnnots", []):
            if ba["_annotevidence_key"] == e["_annotevidence_key"]:
               e["baseAnnots"].append(ba)
               e["baseSubjectIds"].append(ba["genotypeId"])
               if dt is None or ba["_basemodification_date"] < dt:
                   dt = ba["_basemodification_date"]
        if dt:
            e["annotationDate"] = dt
    # Now merge evidence recs having the same publication
    prev = None
    elist = []
    for e in a["evidence"]:
        if prev and prev["_refs_key"] == e["_refs_key"]:
            prev["baseAnnots"] += e["baseAnnots"]
            prev["baseSubjectIds"] += e["baseSubjectIds"]
            prev["codes"] += e["codes"]
        else:
            if (prev) : elist.append(prev)
            prev = e
    if prev:
        elist.append(prev)
    a["evidence"] = elist

    a["baseAnnots"] = []
    
    # object relation for genes and alleles
    if skind == "gene":
        a["objectName"] = a["subjectLabel"]
        a["objectRelation"] = {
            "objectType" : "gene",
            "associationType" : "is_implicated_in"
        }
        for ba in a["baseAnnots"]:
            # record for later use: which genotypes roll up to this gene
            geno2genes.setdefault(ba["genotypeId"], set()).add(a["subjectId"])
    elif skind == "allele":
        a["objectName"] = a["subjectLabel"]
        a["objectRelation"] = {
            "objectType" : "allele",
            "associationType" : "is_implicated_in"
        }
    else: # skind == "genotype"
        a["objectName"] = "???"
        a["objectRelation"] = {
            "objectType" : "genotype",
            "associationType" : "is_model_of",
            # The following line is the reason genes-disease annotation must be processed first...
            "inferredGeneAssociation": list(geno2genes.get(a["subjectId"],[]))
        }
    return a

#
def buildDataProviderObject(a, kind, skind):
    if kind != "disease":
        raise RuntimeError("Cannot build data provider object for this kind: " + kind)
    ident = a["subjectId"]
    if skind == "gene":
        page = "gene"
    elif skind == "allele":
        page = "allele"
    elif skind == "genotype":
        page = "genotype"
    else:
        raise RuntimeError("Cannot build data provider object for this skind: " + skind)
    return {
        "crossReference" : {
            "id" : ident,
            "pages" : [page]
        },
        "type" : "curated"
    }
#
def log (s):
    sys.stderr.write(s + '\n')

########
# Returns a JSON object for one annotation formatted according to the AGR disease or pheno annotation spec.
# Args:
#    annot
#    kind (string) Either "disease" or "phenotype"
#
def formatDafJsonRecord (annot, evidence, kind, skind):
    #
    adate = getTimeStamp(evidence['annotationDate'])
    try:
      if kind == "disease":
        return stripNulls({
            'objectId':                     annot["subjectId"],
            'objectName':                   annot["objectName"],
            'objectRelation':               annot["objectRelation"],
            'negation':                     annot["qualifier"],
            'DOid':                         annot["termId"],
            'evidence':                     {
                "publication" : makePubRef(evidence["refPmid"], evidence["refMgiId"]),
                "evidenceCodes" :  list(set(evidence["codes"])),
                },
            'primaryGeneticEntityIDs':      evidence["baseSubjectIds"],
            'dateAssigned':                 adate,
            'dataProvider':                 [ buildDataProviderObject(annot, kind, skind) ],
        })
      else:
        # guard against data issue: MP record w/o a name. 
        if not annot["term"]:
            return None
        #
        return stripNulls({
            'objectId':                     annot["subjectId"],
            'phenotypeTermIdentifiers':     [{ "termId" : annot["termId"], 'termOrder' : 1 }],
            'phenotypeStatement':           annot["term"],
            'evidence':                     makePubRef(evidence["refPmid"], evidence["refMgiId"]),
            'primaryGeneticEntityIDs':      evidence["baseSubjectIds"],
            'dateAssigned':                 adate,
            })
    except:
        log('ERROR in annotation record: ' + str(annot))
        raise

#####
def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--diseases",
        dest="doDiseases",
        action="store_true",
        default=False,
        help="Output disease annotations.")

    parser.add_argument(
        "-p", "--phenotypes",
        dest="doPhenotypes",
        action="store_true",
        default=False,
        help="Output phenotype annotations.")

    return parser.parse_args()

#####
def main():
    args = getArgs()
    mdo = buildMetaObject()
    print('{"metaData": %s,' % json.dumps(mdo)) 
    print(' "data"    : [')
    if args.doDiseases:
        #
        # IMPORTANT! Gene annotations must be retrieved *before* Genotype annots. 
        geneAnnots = getAnnotations(doGeneCfg)
        alleleAnnots = getAnnotations(doAlleleCfg)
        alleleDirectAnnots = getAnnotations(doAlleleDirectCfg)
        genotypeAnnots = getAnnotations(doGenoCfg)
        #
        for i,ga in enumerate(chain(geneAnnots, alleleAnnots, alleleDirectAnnots, genotypeAnnots)):
            print("," if i>0 else "", json.dumps(ga, indent = 2))
    elif args.doPhenotypes:
        #
        geneAnnots = getAnnotations(mpGeneCfg)
        alleleAnnots = getAnnotations(mpAlleleCfg)
        genotypeAnnots = getAnnotations(mpGenoCfg)
        #
        for i,ga in enumerate(chain(geneAnnots, alleleAnnots, genotypeAnnots)):
            print("," if i>0 else "", json.dumps(ga, indent=2))

    print("]}")

main()  # optional geno and/or gene IDs on the cmd line to restrict query
