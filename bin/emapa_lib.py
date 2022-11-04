# standard libs
import sys
import os
import json
import argparse

# nonstandard dependencies
from AGRlib import stripNulls, buildMetaObject, sql, makeOneOfConstraint, makePubRef

#-----------------------------------
def log(msg):
    sys.stderr.write(msg + '\n')

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
# Returns a mapping from EMAPA id to EMAPA term obj
# Each term has attributes: term, startstage, endstage/
def loadEMAPA ():
    log('Loading EMAPA...')
    id2emapa = {}
    q = '''
    SELECT aa.accid, vt.term, vte.startstage, vte.endstage
    FROM VOC_Term_Emapa vte, VOC_Term vt, ACC_Accession aa
    WHERE vt._term_key = vte._term_key
    AND vt._term_key = aa._object_key
    AND aa._mgitype_key = 13
    AND aa._logicaldb_key = 169
    AND aa.preferred = 1
    '''
    for t in sql(q):
        t["startstage"] = int(t["startstage"])
        t["endstage"] = int(t["endstage"])
        id2emapa[t["accid"]] = t
    log('Loaded %d EMAPA terms.'%len(id2emapa))
    return id2emapa

#-----------------------------------
# Loads/returns a mapping from EMAPA id to the IDs of its immediate parents
def loadEMAPAParents():
    log('Loading EMAPA parents...')

    q = '''
        SELECT ca.accid as childid, pa.accid as parentid
        FROM 
          DAG_Edge e, 
          DAG_Node cn, 
          VOC_Term ct, 
          ACC_Accession ca,
          DAG_Node pn, 
          VOC_Term pt,
          ACC_Accession pa
        WHERE e._child_key = cn._node_key
        AND e._parent_key = pn._node_key
        AND cn._object_key = ct._term_key
        AND pn._object_key = pt._term_key
        AND pt._vocab_key = 90
        AND pt._term_key = pa._object_key
        AND pa._mgitype_key = 13
        AND pa._logicaldb_key = 169
        AND pa.preferred = 1
        AND ct._term_key = ca._object_key
        AND ca._mgitype_key = 13
        AND ca._logicaldb_key = 169
        AND ca.preferred = 1
        '''
    id2pids = {}
    for i,r in enumerate(sql(q)):
        id2pids.setdefault(r["childid"], []).append(r["parentid"])
    log('Loaded %d parent/child relations.' % i)
    return id2pids


#-----------------------------------
# Returns the IDs of all ancestors of the given term at the given stage.
# Includes the term's own ID (ie, reflexive transitive closure)
#
def ancestorsAt (termId, stage) :
    global id2emapa, id2pids
    ancestors = set([termId])
    def _(t):
        for pid in id2pids.get(t["accid"],[]):
            p = id2emapa[pid]
            if p["startstage"] <= stage and p["endstage"] >= stage:
                ancestors.add(pid)
                _(p)
        
    _(id2emapa[termId])
    return ancestors

#
def mkWhereExpressedObj (eid, stage) :
    global id2emapa, id2pids
    structureName = id2emapa[eid]["term"]
    ancs = ancestorsAt(eid, stage)
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
        uids = ['Other']
    return {
      'anatomicalStructureTermId' : eid,
      'anatomicalStructureUberonSlimTermIds': [{'uberonTerm':u} for u in uids],
      'whereExpressedStatement' : structureName
    }

#
id2emapa = loadEMAPA()
id2pids = loadEMAPAParents()

