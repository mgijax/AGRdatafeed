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
from AGRlib import stripNulls, buildMetaObject, makeOneOfConstraint, sql
from AGRqlib import qSubmittedAlleleIds, tConstructRelationships, tConstructProperties, qConstructNonMouseDrivers

EXPRESSES_cat_key = 1004
DRIVER_cat_key = 1006

def log (s) :
    sys.stderr.write(s + '\n')

def loadSubmittedAlleles () :
    aids = set()
    for a in sql(qSubmittedAlleleIds):
        aids.add(a["mgiid"])
    return aids
    
ak2nmdId = {} # allele_key -> non-mouse-driver ID
def loadNonMouseDrivers () :
    for r in sql(qConstructNonMouseDrivers):
        ak2nmdId[r['_allele_key']] = r['accid']

def loadRelationship (key) :
    rels = []
    # read the properties and create index from rel key -> prop obj
    rk2props = {}
    for rp in sql(tConstructProperties % key):
        rk2props.setdefault(rp['_relationship_key'], {})[rp['property']] = rp['value']
    # read the relationships and add properties
    for r in sql(tConstructRelationships % key):
        rk = r['_relationship_key']
        if rk in rk2props:
            r['properties'] = rk2props[rk]
        rels.append(r)
    return rels

def rel2constrComp (r) :
    if r["_category_key"] == EXPRESSES_cat_key:
        reln = "expresses"
        symbol = r["genesymbol"]
        try:
            if r["relationship"] == "expresses_an_orthologous_gene":
                symbol = r["properties"]["Non-mouse_Gene_Symbol"]
                organism = r["properties"]["Non-mouse_Organism"]
                if organism == "human":
                    gid = r["properties"].get("Non-mouse_HGNC_Gene_ID",None)
                elif organism == "rat":
                    gid = r["properties"].get("Non-mouse_RGD_Gene_ID",None)
                elif organism == "zebrafish":
                    gid = "ZFIN:" + r["properties"].get("Non-mouse_ZFIN_Gene_ID",None)
                else:
                    gid = None
            else:
                gid = r["gene"]
        except:
            log("Error processing: " + str(r))
            return None
    elif r["_category_key"] == DRIVER_cat_key:
        reln = "is_regulated_by"
        symbol = r["genesymbol"]
        if r["_organism_key"] == "1":
            gid = r["gene"]
        else:
            gid = ak2nmdId.get(r["_allele_key"], None)
            if gid and gid.startswith('ZDB'):
                gid = 'ZFIN:' + gid
    else:
        raise RuntimeError("Internal error: unknown _category_key: " + str(r))

    return {
        "componentSymbol" : symbol,
        "componentID" : gid,
        "componentRelation" : reln
    }
    
def main () :
    loadNonMouseDrivers()
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

