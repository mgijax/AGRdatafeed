#
# contruct.py
#
# Script to generate construct objects for alleles in MGI.
#
# Alleles that have expressed components and/or driver genes are considered to
# have constructs. One construct is created per allele, which  may contain multiple components.
# (e.g., multiple expressed components, or an expressed component and a driver, or ...)
#
# The alliance ingest schema wants an identifier for each construct. MGI does not have
# actual construct objects, much less identifiers for them. So this script creates ersatz
# IDs for the constructs, based on the owning allele's id. If an allele's id is MGI:123456,
# its contruct has the id MGI:123456_con
#

import sys
import os
import json
from AGRlib import stripNulls, buildMetaObject, sql
from AGRqlib import qSubmittedAlleleIds, tConstructRelationships, qConstructNonMouseComponents

EXPRESSES_cat_key = 1004
DRIVER_cat_key = 1006

def log (s) :
    sys.stderr.write(s + '\n')

def loadSubmittedAlleles () :
    aids = set()
    for a in sql(qSubmittedAlleleIds):
        aids.add(a["mgiid"])
    return aids
    
mk2nmdId = {} # marker key -> non-mouse ID
def loadNonMouseGeneIds () :
    for r in sql(qConstructNonMouseComponents):
        mk2nmdId[r['_marker_key']] = r['accid']

def loadRelationship (key) :
    rels = []
    # read the relationships 
    for r in sql(tConstructRelationships % key):
        rk = r['_relationship_key']
        rels.append(r)
    return rels

def rel2constrComp (r) :
    symbol = r["genesymbol"]
    gid = r["gene"]
    if gid is None:
        gid = mk2nmdId.get(r["_marker_key"],None)
        if gid :
            if gid.startswith('ZDB'):
                gid = 'ZFIN:' + gid
            elif gid.startswith('XB-'):
                gid = 'Xenbase:' + gid

    if r["_category_key"] == EXPRESSES_cat_key:
        reln = "expresses"
    elif r["_category_key"] == DRIVER_cat_key:
        reln = "is_regulated_by"
    else:
        raise RuntimeError("Internal error: unknown _category_key: " + str(r))

    return {
        "componentSymbol" : symbol,
        "componentID" : gid,
        "componentRelation" : reln
    }
    
def main () :
    loadNonMouseGeneIds()
    submittedIds = loadSubmittedAlleles()
    aid2rels = {}
    for r in loadRelationship(EXPRESSES_cat_key):
        aid2rels.setdefault(r['allele'],[]).append(r)
    for r in loadRelationship(DRIVER_cat_key):
        aid2rels.setdefault(r['allele'],[]).append(r)
    aids = list(aid2rels.keys())
    #
    print('{\n  "metaData": %s,\n  "data": [' % json.dumps(buildMetaObject(), indent=2))
    first = True
    for aid in aids:
        if not aid in submittedIds:
            continue
        arels = aid2rels[aid]
        ccomps = list(filter(lambda x:x, map(rel2constrComp, arels)))
        obj = {
          "primaryId" : aid + '_con',
          "name" : arels[0]["allelesymbol"] + ' construct',
          "crossReferences" : [],
          "synonyms" : [],
          "constructComponents": ccomps
        }
        if not first: print(",", end=' ')
        print(json.dumps(stripNulls(obj), indent=2))
        first=False
    print("]}")

#
main()

