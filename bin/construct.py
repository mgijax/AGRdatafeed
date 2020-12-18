#
# contruct.py
#
# Script to generate construct objects for alleles in MGI.
# Alleles that have expressed components or that have driver genes are considered to
# have constructs.
#

import sys
import json
from AGRlib import getConfig, stripNulls, buildMetaObject, makeOneOfConstraint, sql

cp = getConfig()
MOUSEMINE = cp.get("DEFAULT","MOUSEMINEURL")

EXPRESSES_cat_key = 1004
DRIVER_cat_key = 1006

def loadRelationship (key) :
    rels = []
    # read the properties and create index from rel key -> prop obj
    rk2props = {}
    for rp in sql(qProperties % key):
        rk2props.setdefault(rp['_relationship_key'], {})[rp['property']] = rp['value']
    # read the relationships and add properties
    for r in sql(qRelationships % key):
        rk = r['_relationship_key']
        if rk in rk2props:
            r['properties'] = rk2props[rk]
        rels.append(r)
    return rels

def rel2constrComp (r) :
    gid = r["gene"]
    symbol = r["genesymbol"]
    if "properties" in r:
        symbol = r["properties"]["Non-mouse_Gene_Symbol"]
        gid = r["properties"].get("Non-mouse_NCBI_Gene_ID", None)
        if gid:
            gid = "NCBI_Gene:"+gid
    if r["relationship"].startswith("express"):
        reln = "expresses"
    elif r["relationship"] == "has_driver":
        reln = "is_regulated_by"
    return {
        "componentSymbol" : symbol,
        "componentID" : gid,
        "componentRelation" : reln
    }
    
def main () :
    aid2rels = {}
    for r in loadRelationship(EXPRESSES_cat_key):
        aid2rels.setdefault(r['allele'],[]).append(r)
    for r in loadRelationship(DRIVER_cat_key):
        aid2rels.setdefault(r['allele'],[]).append(r)
    aids = list(aid2rels.keys())
    print('{\n  "metaData": %s,\n  "data": [' % json.dumps(buildMetaObject(MOUSEMINE), indent=2))
    first = True
    for aid in aids:
        arels = aid2rels[aid]
        ccomps = list(map(rel2constrComp, arels))
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


# query MGI_Relationship for "expressed_component" (1004) and "has_driver" (1006)
qRelationships = '''
    SELECT 
        r._relationship_key,
        r._category_key,
        aa.accid as allele,
        al.symbol as alleleSymbol,
        am.accid as gene,
        mm.symbol as genesymbol,
        rr.term as relationship,
        q.term as qualifier,
        e.abbreviation as evidencecode,
        r._refs_key
    FROM
        MGI_Relationship r,
        VOC_Term q,
        VOC_Term e,
        VOC_Term rr,
        ACC_Accession aa,
        ACC_Accession am,
        ALL_Allele al,
        MRK_Marker mm
    WHERE r._category_key = %d
    AND r._relationshipterm_key = rr._term_key
    AND r._qualifier_key = q._term_key
    AND r._evidence_key = e._term_key
    AND r._object_key_1 = aa._object_key
    AND aa._mgitype_key = 11
    AND aa._logicaldb_key = 1
    AND aa.preferred = 1
    AND r._object_key_2 = am._object_key
    AND am._mgitype_key = 2
    AND am._logicaldb_key = 1
    AND am.preferred = 1
    AND r._object_key_1 = al._allele_key
    AND r._object_key_2 = mm._marker_key
    ORDER BY r._relationship_key
    '''

# query for relationship properties. Will get attached as list to relationship.
qProperties = '''
    SELECT r._relationship_key, t.term as property, p.value
    FROM MGI_Relationship r
      JOIN MGI_Relationship_Property p
        ON r._relationship_key = p._relationship_key
      JOIN VOC_Term t
        ON p._propertyname_key = t._term_key
    WHERE r._category_key = %d
    ORDER BY r._relationship_key, p.sequenceNum
    '''

# go!
main()

