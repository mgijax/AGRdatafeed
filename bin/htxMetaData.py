#
# htxMetaData.py
#
# Dumps curated metadata for high throughput expression studies 
# (experiments + samples).
#

import sys
import re
import json
import itertools
import time
import types
import argparse
from emapa_lib import mkWhereExpressedObj

# nonstandard dependencies
from AGRlib import stripNulls, buildMetaObject, doQuery, makePubRef, getConfig, getTimeStamp

TIMESTAMP = getTimeStamp()

#-----------------------------------
# load config settings

cp = getConfig()

MOUSEMINE     = cp.get("DEFAULT","MOUSEMINEURL")
GLOBALTAXONID = cp.get("DEFAULT","GLOBALTAXONID")

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
    
def getVariables(url):
    eid2vars = {}
    for r in doQuery(htVariables, url):
        if r["variables.name"]:
            eid2vars.setdefault(r['experimentId'], []).append(r)
    return eid2vars
    
def getSamples(url):
    eid2samples = {}
    for r in doQuery(htSamples, url):
        eid2samples.setdefault(r['experimentId'], []).append(r)
    return eid2samples

def getExperiments (url) :
    for r in doQuery(htExperiments, url):
        yield r

#-----------------------------------
def getCatTags (obj) :
    vs = [ obj["studyType"] ] 
    vs += [ v["variables.name"] for v in obj["variables"]]
    return vs

def getSex (obj) :
    # Male, Female, Pooled, Not Specified
    #  =>
    # male, female, pooled, unknown
    if obj["samples.sex"] == "Not Specified":
        return "unknown"
    return obj["samples.sex"].lower()

def getAssayType (exptType) :
    if exptType == "transcription profiling by array":
        return "MMO:0000648"
    elif exptType == "RNA-Seq":
        return "MMO:0000659"
    else:
        raise RuntimeError("Unknown experiment type: " + str(exptType))

def getSampleJsonObj (obj) :
    return stripNulls({
        "sampleTitle" : obj["samples.name"],
        "datasetId" : [ obj["experimentId"] ],
        "assayType" : getAssayType(obj["experimentType"]),
        "sampleAge" : { "age" : obj["samples.age"] },
        "sampleLocation" : [ mkWhereExpressedObj(obj["samples.structure.identifier"], int(obj["samples.stage"])) ],
        "sex" : getSex(obj),
        "taxonId" : GLOBALTAXONID,
        "assemblyVersion" : [ "GRCm38.p6" ],
        "dateAssigned" : TIMESTAMP, #FIXME
    })

def getExptJsonObj(obj):
    return stripNulls({
        "datasetId" : { "primaryId" : obj["experimentId"] },
        "title" : obj["name"],
        "summary" : obj["description"],
        "categoryTags" : getCatTags(obj),
        "publication" : [
            makePubRef(p["publications.pubMedId"], p["publications.mgiId"]) for p in obj["references"]],
        "dateAssigned" : TIMESTAMP, #FIXME
    })
#-----------------------------------
def getHTdata (kind):
    eid2samples = getSamples(MOUSEMINE)
    eid2refs = getReferences(MOUSEMINE)
    eid2vars = getVariables(MOUSEMINE)
    for e in getExperiments(MOUSEMINE):
        e['samples'] = eid2samples.get(e['experimentId'],[])
        e['variables'] = eid2vars.get(e['experimentId'],[])
        e['references'] = eid2refs.get(e['experimentId'],[])
        if kind == "experiments":
            yield getExptJsonObj(e)
        else:
            for s in e['samples']:
                yield getSampleJsonObj(s)
            
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
            "
        >
    </query>
    '''
htSamples = '''
    <query
        model="genomic"
        view="
            HTExperiment.experimentId
            HTExperiment.experimentType
            HTExperiment.samples.name
            HTExperiment.samples.sex
            HTExperiment.samples.age
            HTExperiment.samples.stage
            HTExperiment.samples.structure.identifier
            HTExperiment.samples.structure.name
            HTExperiment.samples.genotype.primaryIdentifier
            "
        >
    </query>
    '''
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
