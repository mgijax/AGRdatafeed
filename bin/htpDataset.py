#
# htxMetaData.py
#
# Dumps curated metadata for high throughput expression studies 
# (experiments + samples).
#

import sys
import os
import re
import json
import itertools
import time
import types
import argparse
from emapa_lib import mkWhereExpressedObj

# nonstandard dependencies
from AGRlib import stripNulls, buildMetaObject, sql, makePubRef, getTimeStamp

TIMESTAMP = getTimeStamp()

AGE2BIN = [
    [42.01, 99999, "UBERON:0000113"], # post-juvenile adult
    [21.01, 42.0, "post embryonic, pre-adult"],
    [0.0, 21.0, "UBERON:0000068"],  # embryo
]

#-----------------------------------

def parseCmdLine():
    parser = argparse.ArgumentParser(description='Dumps metadata for high-throughput experiments and samples to a JSON file.')

    parser.add_argument(
      '-d', '--dump',
      choices=['experiments','samples'],
      help='Specify what to dump.',
      required=True)

    args = parser.parse_args()
    return args
  
#-----------------------------------
def getReferences():
    #
    pmid2mgi = {}
    for r in sql(htPmid2Mgi):
        pmid2mgi[r['pmid']] = r['mgiid']
    #
    eid2refs = {}
    for r in sql(htPmids):
        pmid = r['pmid']
        mgiid = pmid2mgi.get(pmid,None)
        rr = { "pmid": pmid, "mgiid": mgiid }
        eid2refs.setdefault(r['experimentId'], []).append(rr)
    return eid2refs
    
#-----------------------------------
def getVariables():
    eid2vars = {}
    for r in sql(htVariables):
        if r["variable"]:
            eid2vars.setdefault(r['experimentId'], []).append(r)
    return eid2vars
    
#-----------------------------------
def getSamples():
    eid2samples = {}
    for r in sql(htSamples):
        rkey = (r["name"],r["age"],r["sex"],r["structureId"])
        eid2samples.setdefault(r['experimentId'], {})[rkey] = r
    return eid2samples

#-----------------------------------
def getExperiments () :
    for r in sql(htExperiments):
        yield r

#-----------------------------------
def getCatTags (obj) :
    vs = [ obj["studyType"] ] 
    vs += [ v["variable"] for v in obj["variables"]]
    vs = list(map(lambda x: TAG_MAP.get(x, x), vs))
    for x in vs:
        if not x in TAGS:
            sys.stderr.write('Tag not recognized: %s\n' % x)
    return vs

#-----------------------------------
def getSex (obj) :
    # Male, Female, Pooled, Not Specified
    #  =>
    # male, female, pooled, unknown
    if obj["sex"] is None or obj["sex"] == "Not Specified":
        return "unknown"
    return obj["sex"].lower()

#-----------------------------------
def getAssayType (exptType) :
    if exptType == "transcription profiling by array":
        return "MMO:0000648"
    elif exptType == "RNA-Seq":
        return "MMO:0000659"
    else:
        raise RuntimeError("Unknown experiment type: " + str(exptType))

#-----------------------------------
def getHTdata (kind):
    eid2samples = getSamples()
    eid2refs = getReferences()
    eid2vars = getVariables()
    for e in getExperiments():
        eid = e['experimentId']
        e['samples'] = (eid2samples.get(eid,{}).values())
        e['variables'] = eid2vars.get(eid,[])
        e['references'] = eid2refs.get(eid,[])
        e['curationDate'] = getTimeStamp(e['curationDate'])
        #
        if kind == "experiments":
            yield getExptJsonObj(e)
        else:
            for s in e['samples']:
                s['experimentId'] = e['experimentId']
                s['curationDate'] = e['curationDate']
                yield getSampleJsonObj(s)
            
#-----------------------------------
def findUberonTerm (ageMin, ageMax) :
    ageMin = float(ageMin)
    ageMax = float(ageMax)
    for i in range(len(AGE2BIN)) :
        bMin, bMax, uberonTerm = AGE2BIN[i]
        if ageMin <= bMax and ageMax >= bMin:
            return uberonTerm
    raise RuntimeError("Could not find Uberon term for range (%1.2f, %1.2f)" % (ageMin, ageMax))

#-----------------------------------
def getStageObj (obj) :
    try:
        u = findUberonTerm(obj["ageMin"], obj["ageMax"])
    except:
        raise RuntimeError("Could not find Uberon age term for:" + str(obj))

    return {
        "stageName" : "TS" + str(obj["stage"]),
        "stageUberonSlimTerm": {"uberonTerm":u}
    }
#-----------------------------------
def getSampleJsonObj (obj) :
    return stripNulls({
        "sampleTitle" : obj["name"],
        "datasetIds" : [ "ArrayExpress:" + obj["experimentId"] ],
        "assayType" : getAssayType(obj["experimentType"]),
        "sampleType": "OBI:0000880", # "RNA extract"
        "sampleAge" : { "age" : obj["age"], "stage" : getStageObj(obj) },
        "sampleLocations" : [ mkWhereExpressedObj(obj["structureId"], int(obj["stage"])) ],
        "sex" : getSex(obj),
        "taxonId" : "NCBITaxon:10090",
        "assemblyVersions" : [ "GRCm39" ],
        "dateAssigned" : obj['curationDate']
    })

#-----------------------------------
GEO_re = re.compile('^E-GEOD-(\d+)$')
def getExptJsonObj(obj):
    # If it's a GEO experiment, also include the GEO id. We should get this from MouseMine,
    # but MM doesn't have this. So construct it from the AE version of the ID. 
    # Example: E-GEOD-33885 (ArrayExpress) -> GSE33885 (GEO)
    pid = "ArrayExpress:" + obj["experimentId"]
    xrefs = [{ "id" : pid, "pages": ["htp/dataset"] }]
    geoXref = None
    geoidmatch = GEO_re.match(obj["experimentId"])
    if geoidmatch:
        geoid = 'GEO:GSE'+geoidmatch.group(1)
        xrefs.append({ "id" : geoid, "pages": ["htp/dataset"] })

    return stripNulls({
        "datasetId" : {
            "primaryId" : pid,
            "preferredCrossReference" : {
                 "id" : "MGI:" + obj["experimentId"],
                 "pages": ["htp/dataset"] ,
            },
            "crossReferences" : xrefs,
        },
        "title" : obj["name"],
        "summary" : obj["description"],
        "categoryTags" : getCatTags(obj),
        "publications" : [
            makePubRef(p["pmid"], p["mgiid"]) for p in obj["references"]],
        "dateAssigned" : obj['curationDate'],
    })
#-----------------------------------
def main() :
    args = parseCmdLine()

    first = True
    mdo = json.dumps(buildMetaObject(), indent=2)
    print('{\n  "metaData": %s,\n  "data": [' % mdo)
    for obj in getHTdata(args.dump):
        if not first: print(",", end=' ')
        print(json.dumps(obj, indent=2))
        first=False
    print("]}")


#-----------------------------------
TAGS=set('''
age
aging
cell aging
amino acid metabolism
amino acid utilization
anatomical structure
baseline
bioinformatics and computational biology
biotic stimulus
carbon utilization
cell cycle
cell morphogenesis
cell type
cell wall organization
cellular ion homeostasis
chemical stimulus
chromatin organization
chromosome segregation
time of day
cofactor metabolism
colony morphology
cytoskeleton and molecular motors
dauer larval development
developmental stage
developmental time course
diauxic shift
disease
DNA binding
DNA damage stimulus
DNA replication, recombination and repair
dosage compensation
embryo development
environmental-sensing
evolution
fermentation
filamentous growth
flocculation
gene silencing by miRNA
gene
genetic interaction
genome variation
genotype
histone modification
innate immune response
lipid metabolic process
mating
metabolism
metabolism and metabolic regulation
metal or metalloid ion stress
metal or metalloid ion utilization
mitotic cell cycle
mRNA processing
mouse species
natural variant
nitrogen utilization
nuclear structure and function
nucleotide metabolism
nutrient utilization
subcellular component
oxidative stress
oxygen level alteration
phosphorus utilization
physical interaction
physical stimulus
ploidy
pregnancy
prions
programmed cell death
protein dephosphorylation
protein glycosylation
protein modification
protein phosphorylation
protein structure and folding
protein trafficking, localization and degradation
proteolysis
QTL
respiration
response to osmotic stress
response to radiation
response to starvation
response to temperature stimulus
response to unfolded protein
RNA catabolism
RNA interference
RNA processing and metabolism
RNA structure
RNAi gene function study
sex
signaling
single cell variation
space flight
species
sporulation
stationary phase
stationary phase entry
stationary phase maintenance
strain study
stress
sulfur utilization
synthetic biology
transcription
transcriptional regulation
transcriptome
translational regulation
transposons
ubiquitin or ULP modification
vulval development
WT vs. mutant
unclassified
'''.strip().split('\n'))
#-----------------------------------
TAG_MAP = {
    "WT vs. Mutant" : "WT vs. mutant",
    "Baseline" : "baseline",
    "mouse strain" : "strain study"
}
#-----------------------------------
htExperiments = '''
    SELECT 
      a.accid as "experimentId",
      name, 
      et.term as "experimentType", 
      e.description, 
      st.term as "studyType", 
      src.term as source, 
      e.last_curated_date as "curationDate"
    FROM 
       GXD_HTExperiment e,
       VOC_Term et,
       VOC_Term st,
       VOC_Term src,
       ACC_Accession a
    WHERE e._curationstate_key = 20475421 /* Done */
    AND e._experimenttype_key =  et._term_key
    AND e._studytype_key = st._term_key
    AND e._source_key = src._term_key
    AND e._experiment_key = a._object_key
    AND a._mgitype_key = 42 /* ht experiment type */
    AND a.preferred = 1
    '''
#-----------------------------------
htSamples = '''
    SELECT 
      ea.accid as "experimentId", 
      et.term as "experimentType",
      s.name, 
      s.age, 
      s.agemin, 
      s.agemax, 
      st.term as sex, 
      s._stage_key as stage, 
      em.term as "structureName", 
      ema.accid as "structureId",
      ga.accid as "genotypeId"
    FROM 
      GXD_HTSample s, 
      GXD_HTExperiment e,
      ACC_Accession ea,
      VOC_Term st,
      VOC_Term em,
      ACC_Accession ema,
      ACC_Accession ga,
      VOC_Term et
    WHERE e._experiment_key = s._experiment_key
    AND e._experimenttype_key = et._term_key
    AND s._relevance_key = 20475450 /* Yes */
    AND s._experiment_key = ea._object_key
    AND ea._mgitype_key = 42 /* ht experiment type */
    AND ea.preferred = 1
    AND s._sex_key = st._term_key
    AND s._emapa_key = em._term_key
    AND s._emapa_key = ema._object_key
    AND ema._mgitype_key = 13
    AND ema._logicaldb_key = 169 /* EMAPA */
    AND ema.preferred = 1
    AND s._genotype_key = ga._object_key
    AND ga._mgitype_key = 12 
    AND ga._logicaldb_key = 1
    AND ga.preferred = 1
    '''

#-----------------------------------
htPmids = '''
    SELECT 
      p.value as pmid,
      a.accid as "experimentId"
    FROM
      GXD_HTExperiment e,
      MGI_Property p,
      ACC_Accession a
    WHERE p._propertytype_key = 1002
    AND p._object_key = e._experiment_key
    AND p._propertyterm_key = 20475430 /* PMID */
    AND e._experiment_key = a._object_key
    AND a._mgitype_key = 42 /* ht experiment type */
    AND a.preferred = 1
    '''

htPmid2Mgi = '''
    SELECT 
      a.accid as pmid, a2.accid as mgiid
    FROM
      ACC_Accession a,
      ACC_Accession a2
    WHERE a._logicaldb_key = 29 /* Pubmed */
    AND a._mgitype_key = 1
    AND a.preferred = 1
    AND a._object_key = a2._object_key
    AND a2._mgitype_key = 1
    AND a2._logicaldb_key = 1
    AND a2.prefixpart = 'MGI:'
    AND a2.preferred = 1
    '''

#-----------------------------------
htVariables = '''
    SELECT
      a.accid as "experimentId", t.term as variable
    FROM
      GXD_HtExperimentVariable v,
      ACC_Accession a,
      VOC_Term t
    WHERE v._experiment_key = a._object_key 
    AND a._mgitype_key = 42 /* ht experiment type */
    AND a.preferred = 1
    AND v._term_key = t._term_key
    AND t.term != 'Not Applicable'
    '''
#-----------------------------------
main()
