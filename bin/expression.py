#!/usr/bin/env python2.7 
#
# expression.py
#
# Script to dump expression info in AGR json format.
# The format is described here:
#       https://github.com/alliance-genome/agr_schemas
#
# See: https://github.com/alliance-genome/agr_schemas/releases/tag/1.0.0.4
# Of particular note (from the doc):
#   4. an expression annotation has one gene, one where, one pub, one assay and one when 
#
# The data to be dumped include positve expression results in wild type specimens.
# NOTE: In situ reporter knockin hets are considered wild type.
#
# Assay types mapped to MMO terms.
#
# See: https://docs.google.com/spreadsheets/d/1_UcKTq7y-wsQ83_kJlP6X5mQdgKykaPWlDdkZzCuygI/edit#gid=313336440
#
# FIXME: Need annotation dates!
#
# Theiler stages mapped to:
# See https://docs.google.com/spreadsheets/d/1_UcKTq7y-wsQ83_kJlP6X5mQdgKykaPWlDdkZzCuygI/edit#gid=463437957
#
# Within a given assay, only report one result per anatomy/stage combo. (Ie, if the assay contains the
# same time/space result in multiple specimens, only output one row.)
#
# Author: Joel Richardson
#

# standard libs
import sys
import os
import re
import json
import itertools
import time
import types
import argparse

# nonstandard dependencies
from AGRlib import stripNulls, buildMetaObject, doQuery, makeOneOfConstraint, makePubRef

#-----------------------------------
# Mapping from our assay type to MMO ids
# courtesy of Connie Smith @MGI.
assayType2mmo = dict([
    ('Immunohistochemistry','MMO:0000498'),
    ('In situ reporter (knock in)','MMO:0000672'),
    ('Northern blot','MMO:0000647'),
    ('Nuclease S1','MMO:0000653'),
    ('RT-PCR','MMO:0000655'),
    ('RNase protection','MMO:0000652'),
    ('RNA in situ','MMO:0000658'),
    ('Western blot','MMO:0000669'),
])
 
#-----------------------------------
# mappings from high level EMAPA terms to UBERON terms
# mappings from: https://docs.google.com/spreadsheets/d/1_UcKTq7y-wsQ83_kJlP6X5mQdgKykaPWlDdkZzCuygI/edit#gid=313336440
uberonEmapaTbl = [row.split('\t') for row in ('''
UBERON:0001009\tEMAPA:16104\tcardiovascular system
UBERON:0005409\tEMAPA:16246\talimentary system
UBERON:0000949\tEMAPA:35306\tendocrine system
UBERON:0001008\tEMAPA:17366\turinary system
UBERON:0002330\tEMAPA:35329\texocrine system
UBERON:0002193\tEMAPA:18765\themolymphoid system
UBERON:0002416\tEMAPA:17524\tintegumental system
UBERON:0002423\tEMAPA:16840\tliver and biliary system
UBERON:0002204\tEMAPA:32714\tmusculoskeletal system
UBERON:0001016\tEMAPA:16469\tnervous system
UBERON:0000990\tEMAPA:17381\treproductive system
UBERON:0001004\tEMAPA:16727\trespiratory system
UBERON:0001032\tEMAPA:16192\tsensory organ system
UBERON:0005726\tEMAPA:36004\tolfactory system
UBERON:0005726\tEMAPA:36885\tgustatory system
UBERON:0007037\tN/A
UBERON:0002105\tEMAPA:37985\tvestibulo-auditory system
UBERON:0002104\tEMAPA:36003\tvisual system
UBERON:0000924\tEMAPA:35985\tectoderm
UBERON:0000925\tEMAPA:35986\tendoderm
UBERON:0000926\tEMAPA:35987\tmesoderm
UBERON:0003104\tEMAPA:16097\tmesenchyme
UBERON:0001013\tEMAPA:35112\tadipose tissue
UBERON:0000026\tEMAPA:37283\tappendage
UBERON:0016887\tEMAPA:16042\textraembryonic component
UBERON:6005023\tN/A
UBERON:0002539\tEMAPA:16117\tbranchial arch
'''.strip().split('\n'))]

# index the table
emapa2uberon = {}
for r in uberonEmapaTbl:
    if r[1] != 'N/A':
        emapa2uberon[r[1]] = r[0]

# the set of high level EMAPA term IDs
highlevelemapa = set(emapa2uberon.keys())

#-----------------------------------
# mappings from Theiler stages to UBERON stage term IDs
ts2uberon = dict( \
    [('TS%02d'%(i+1),'UBERON:0000068') for i in range(26)] + \
    [('TS27','post embryonic, pre-adult')] + \
    [('TS28', 'UBERON:0000113')])

#-----------------------------------
# MouseMine connection

MOUSEMINE  = os.environ["MOUSEMINEURL"]

#-----------------------------------
def log(msg):
    sys.stderr.write(msg + '\n')

#-----------------------------------
# Returns a mapping from EMAPA id to EMAPA term obj
# Each term has attributes: name, startsAt, endsAt/
def loadEMAPA (url):
    log('Loading EMAPA...')
    id2emapa = {}
    q = '''
    <query
      model="genomic"
      view="
      EMAPATerm.identifier
      EMAPATerm.name
      EMAPATerm.startsAt
      EMAPATerm.endsAt
      "
      >
      </query>
    '''
    for t in doQuery(q, url):
        t["startsAt"] = int(t["startsAt"])
        t["endsAt"] = int(t["endsAt"])
        id2emapa[t["identifier"]] = t
    log('Loaded %d EMAPA terms.'%len(id2emapa))
    return id2emapa

#-----------------------------------
# Loads/returns a mapping from EMAPA id to the IDs of its immediate parents
def loadEMAPAParents(url):
    log('Loading EMAPA parents...')

    q = '''<query
    name=""
    model="genomic"
    view="OntologyRelation.childTerm.identifier
    OntologyRelation.parentTerm.identifier"
    longDescription=""
    sortOrder="OntologyRelation.childTerm.identifier asc"
    >
        <constraint path="OntologyRelation.childTerm" type="EMAPATerm"/>
        <constraint path="OntologyRelation.parentTerm" type="EMAPATerm"/>
        <constraint path="OntologyRelation.direct" op="=" value="true"/>
    </query>'''
    id2pids = {}
    for i,r in enumerate(doQuery(q, url)):
        id2pids.setdefault(r["childTerm.identifier"], []).append(r["parentTerm.identifier"])
    log('Loaded %d parent/child relations.' % i)
    return id2pids


#-----------------------------------
# Returns the IDs of all ancestors of the given term at the given stage.
# Includes the term's own ID (ie, reflexive transitive closure)
#
def ancestorsAt (termId, stage, id2emapa, id2pids) :
    ancestors = set([termId])
    def _(t):
        for pid in id2pids.get(t["identifier"],[]):
            p = id2emapa[pid]
            if p["startsAt"] <= stage and p["endsAt"] >= stage:
                ancestors.add(pid)
                _(p)
        
    _(id2emapa[termId])
    return ancestors

#-----------------------------------
# Returns unique-ified expression assay results. 
#
def getExpressionData(url, ids):
  log('Getting expression data...')
  q = '''<query
    model="genomic"
    view="
        GXDExpression.assayId
        GXDExpression.assayType
        GXDExpression.feature.primaryIdentifier
        GXDExpression.stage
        GXDExpression.structure.identifier
        GXDExpression.publication.mgiId
        GXDExpression.publication.pubMedId"
    sortOrder="GXDExpression.assayId asc GXDExpression.structure.identifier asc GXDExpression.stage asc"
    constraintLogic="A and (B or (C and D)) and E"
    >
      <constraint path="GXDExpression.detected" code="A" op="=" value="true"/>
      <constraint path="GXDExpression.genotype.hasMutantAllele" code="B" op="=" value="false"/>
      <constraint path="GXDExpression.assayType" code="C" op="=" value="In situ reporter (knock in)"/>
      <constraint path="GXDExpression.genotype.zygosity" code="D" op="=" value="ht"/>
      %s
    </query>
  ''' % makeOneOfConstraint('GXDExpression.feature.primaryIdentifier', ids)
  prev = None
  qcount = 0
  ycount = 0
  for r in doQuery(q, MOUSEMINE):
      qcount += 1
      if not prev \
      or r["assayId"] != prev["assayId"] \
      or r["stage"] != prev["stage"] \
      or r["structure.identifier"] != prev["structure.identifier"]:
          ycount += 1
          yield r
      #
      prev = r
      #
  log('getExpressionData: %d results => %d unique results' % (qcount, ycount))

#
def getWhereExpressedObj (eid) :
    pass

# Here is the magic by which an object returned by the query is converted to an object
# conforming to the spec.
#
def getJsonObj(obj, structureName, uids):
  mkid = lambda i,p: None if i is None else p+i
  try:
      return stripNulls({
          'geneId': obj['feature.primaryIdentifier'],
          'evidence' : makePubRef(obj['publication.pubMedId'], obj['publication.mgiId']),
          'assay': assayType2mmo[obj['assayType']],
          'dateAssigned' : '2018-07-18T13:27:43-04:00', # FIXME
          'whereExpressed': {
              'anatomicalStructureTermId' : obj['structure.identifier'],
              'anatomicalStructureUberonSlimTermIds': [{'uberonTerm':u} for u in uids],
              'whereExpressedStatement' : structureName
          },  
          'whenExpressed': {
              'stageName': obj['stage'],
              'stageUberonSlimTerm': {'uberonTerm':ts2uberon[obj['stage']]}
          },  
          'crossReference' : {
              'id' : obj['assayId'],
              'pages' : [ 'gene/expression/annotation/detail' ]
          }
      })
  except:
    log('ERROR in expression object: ' + str(obj))
    raise

#
def parseCmdLine():
    parser = argparse.ArgumentParser(description='Dumps basic allele information to a JSON file.')

    parser.add_argument(
      'identifiers', 
      metavar='ids',
      nargs='*',
      help='Specific MGI ids to dump.')

    args = parser.parse_args()
    return args
  
# Main prog. Build the query, run it, and output 
def main():
  noMapping = set()
  args = parseCmdLine()
  ids = args.identifiers
  #
  id2emapa = loadEMAPA(MOUSEMINE)
  id2pids = loadEMAPAParents(MOUSEMINE)

  '''
  print 'Telencephalon ancestors at stage 15:'
  ancs = ancestorsAt('EMAPA:16910', 15, id2emapa, id2pids)
  for aid in ancs:
      a = id2emapa[aid]
      print a.identifier, a.name, a.startsAt, a.endsAt
      
  print 'Telencephalon ancestors at stage 18:'
  ancs = ancestorsAt('EMAPA:16910', 18, id2emapa, id2pids)
  for aid in ancs:
      a = id2emapa[aid]
      print a.identifier, a.name, a.startsAt, a.endsAt
      
  '''

  exprData = getExpressionData(MOUSEMINE, ids)
  print('{ "metaData" : %s, ' % json.dumps(buildMetaObject(MOUSEMINE)))
  print('  "data" : [')
  for i,r in enumerate(exprData):
      if i: print(",", end=' ')
      eid = r['structure.identifier']
      structureName = id2emapa[eid]["name"]
      s = int(r['stage'][2:])
      ancs = ancestorsAt(eid, s, id2emapa, id2pids)
      # the high level EMAPA ids this annot rolls up to (intersect my ancestors with the HL EMAPA set)
      hla = highlevelemapa & ancs
      # the uberon IDs these map to
      uids = set([emapa2uberon.get(a,'Other') for a in hla])
      '''
      if eid == 'EMAPA:35177':
          log('\n' + eid + ' ' + structureName + ' TS ' + str(s))
          log('ancs=' + str(ancs))
          log('hla=' + str(hla))
          log('uids=' + str(uids))
      '''
      if len(uids) == 0:
          noMapping.add((eid, structureName))
          uids = ['Other']
      # get the JSON object for this annotation
      jobj = getJsonObj(r, structureName, uids)
      #
      print(json.dumps(jobj, sort_keys=True, indent=2, separators=(',', ': ')))
  print(']}')
  noMapping = list(noMapping)
  noMapping.sort()
  for p in noMapping:
      log('No UBERON mapping for: %s %s' % p)

#
if __name__ == "__main__":
    main()
