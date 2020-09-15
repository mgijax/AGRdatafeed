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
        eid2samples.setdefault(r['experimentId'], []).append(r)
    return eid2samples

#-----------------------------------
def getExperiments (url) :
    for r in doQuery(htExperiments, url):
        yield r

#-----------------------------------
def getCatTags (obj) :
    vs = [ obj["studyType"] ] 
    vs += [ v["variables.name"] for v in obj["variables"]]
    return vs

#-----------------------------------
def getSex (obj) :
    # Male, Female, Pooled, Not Specified
    #  =>
    # male, female, pooled, unknown
    if obj["samples.sex"] == "Not Specified":
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
        e['samples'] = eid2samples.get(eid,[])
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
def getSampleJsonObj (obj) :
    return stripNulls({
        "sampleTitle" : obj["samples.name"],
        "datasetIds" : [ "ArrayExpress:" + obj["experimentId"] ],
        "assayType" : getAssayType(obj["experimentType"]),
        "sampleType": "OBI:0000880", # "RNA extract"
        "sampleAge" : { "age" : obj["samples.age"] },
        "sampleLocations" : [ mkWhereExpressedObj(obj["samples.structure.identifier"], int(obj["samples.stage"])) ],
        "sex" : getSex(obj),
        "taxonId" : "NCBITaxon:" + obj["samples.organism.taxonId"],
        "assemblyVersions" : [ "GRCm38.p6" ],
        "dateAssigned" : obj['curationDate']
    })

#-----------------------------------
def getExptJsonObj(obj):
    pid = "ArrayExpress:" + obj["experimentId"]
    return stripNulls({
        "datasetId" : {
            "primaryId" : pid,
            "crossReferences" : [
                { "id" : "MGI:" + obj["experimentId"], "pages": ["htp/dataset"] },
                { "id" : pid, "pages": ["htp/dataset"] },
                ]
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
            HTExperiment.samples.stage
            HTExperiment.samples.structure.identifier
            HTExperiment.samples.structure.name
            HTExperiment.samples.genotype.primaryIdentifier
            HTExperiment.samples.organism.taxonId
            "
        >
        <constraint path="HTExperiment.samples.organism.taxonId" op="=" value="10090"/>
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
