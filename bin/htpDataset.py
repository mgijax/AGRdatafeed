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
from AGRlib import stripNulls, buildMetaObject, doQuery, makePubRef, getTimeStamp

TIMESTAMP = getTimeStamp()

AGE2BIN = [
    [42.01, 99999, "UBERON:0000113"], # post-juvenile adult
    [21.01, 42.0, "post embryonic, pre-adult"],
    [0.0, 21.0, "UBERON:0000068"],  # embryo
]

#-----------------------------------

MOUSEMINE     = os.environ["MOUSEMINEURL"]

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
def getReferences(url):
    eid2refs = {}
    for r in doQuery(htReferences, url):
        eid2refs.setdefault(r['experimentId'], []).append(r)
    return eid2refs
    
#-----------------------------------
def getVariables(url):
    eid2vars = {}
    for r in doQuery(htVariables, url):
        if r["variables.name"]:
            eid2vars.setdefault(r['experimentId'], []).append(r)
    return eid2vars
    
#-----------------------------------
def getSamples(url):
    eid2samples = {}
    for r in doQuery(htSamples, url):
        rkey = (r["samples.name"],r["samples.age"],r["samples.sex"],r["samples.structure.identifier"])
        eid2samples.setdefault(r['experimentId'], {})[rkey] = r
    return eid2samples

#-----------------------------------
def getExperiments (url) :
    for r in doQuery(htExperiments, url):
        yield r

#-----------------------------------
def getCatTags (obj) :
    vs = [ obj["studyType"] ] 
    vs += [ v["variables.name"] for v in obj["variables"]]
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
    if obj["samples.sex"] is None or obj["samples.sex"] == "Not Specified":
        return "unknown"
    return obj["samples.sex"].lower()

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
    eid2samples = getSamples(MOUSEMINE)
    eid2refs = getReferences(MOUSEMINE)
    eid2vars = getVariables(MOUSEMINE)
    for e in getExperiments(MOUSEMINE):
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
        u = findUberonTerm(obj["samples.ageMin"], obj["samples.ageMax"])
    except:
        raise RuntimeError("Could not find Uberon age term for:" + str(obj))

    return {
        "stageName" : "TS" + obj["samples.stage"],
        "stageUberonSlimTerm": {"uberonTerm":u}
    }
#-----------------------------------
def getSampleJsonObj (obj) :
    return stripNulls({
        "sampleTitle" : obj["samples.name"],
        "datasetIds" : [ "ArrayExpress:" + obj["experimentId"] ],
        "assayType" : getAssayType(obj["experimentType"]),
        "sampleType": "OBI:0000880", # "RNA extract"
        "sampleAge" : { "age" : obj["samples.age"], "stage" : getStageObj(obj) },
        "sampleLocations" : [ mkWhereExpressedObj(obj["samples.structure.identifier"], int(obj["samples.stage"])) ],
        "sex" : getSex(obj),
        "taxonId" : "NCBITaxon:" + obj["samples.organism.taxonId"],
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
            makePubRef(p["publications.pubMedId"], p["publications.mgiId"]) for p in obj["references"]],
        "dateAssigned" : obj['curationDate'],
    })
#-----------------------------------
def main() :
    args = parseCmdLine()

    first = True
    mdo = json.dumps(buildMetaObject(MOUSEMINE), indent=2)
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
    <query
        model="genomic"
        view="
            HTExperiment.experimentId
            HTExperiment.seriesId
            HTExperiment.name
            HTExperiment.experimentType
            HTExperiment.studyType
            HTExperiment.source
            HTExperiment.description
            HTExperiment.curationDate
            "
        >
    </query>
    '''
#-----------------------------------
htSamples = '''
    <query
        model="genomic"
        view="
            HTExperiment.experimentId
            HTExperiment.experimentType
            HTExperiment.samples.name
            HTExperiment.samples.sex
            HTExperiment.samples.age
            HTExperiment.samples.ageMin
            HTExperiment.samples.ageMax
            HTExperiment.samples.stage
            HTExperiment.samples.structure.identifier
            HTExperiment.samples.structure.name
            HTExperiment.samples.genotype.primaryIdentifier
            HTExperiment.samples.organism.taxonId
            "
        >
        <constraint path="HTExperiment.samples.organism.taxonId" op="=" value="10090"/>
        <constraint path="HTExperiment.samples.ageMin" op="!=" value="-1"/>
    </query>
    '''
#-----------------------------------
htReferences = '''
    <query
        model="genomic"
        view="
            HTExperiment.experimentId
            HTExperiment.publications.mgiId
            HTExperiment.publications.pubMedId
            "
        >
    </query>
    '''

#-----------------------------------
htVariables = '''
    <query
        model="genomic"
        view="HTExperiment.experimentId
        HTExperiment.variables.name"
        >
    </query>
'''
#-----------------------------------
main()
