#
# contruct.py
#
# Script to generate construct objects for alleles in MGI.
# Alleles that have expressed components or that have driver genes are considered to
# have constructs.
#

import sys
import json
from AGRlib import getConfig, stripNulls, buildMetaObject, makeOneOfConstraint, sql, doQuery

cp = getConfig()
MOUSEMINE = cp.get("DEFAULT","MOUSEMINEURL")

EXPRESSES_cat_key = 1004
DRIVER_cat_key = 1006

def loadSubmittedAlleles () :
    aids = set()
    for a in doQuery(qAlleles, MOUSEMINE):
        aids.add(a["primaryIdentifier"])
    return aids
    
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

def loadNcbi2HgncMapping () :
    global ncbi2hgnc
    ncbi2hgnc = {}
    for r in sql(qNcbi2Hgnc):
        ncbi2hgnc["NCBI_Gene:" + r['NCBI']] = r['HGNC']
    return ncbi2hgnc

def rel2constrComp (r) :
    global ncbi2hgnc
    gid = r["gene"]
    symbol = r["genesymbol"]
    if "properties" in r:
        symbol = r["properties"]["Non-mouse_Gene_Symbol"]
        gid = r["properties"].get("Non-mouse_NCBI_Gene_ID", None)
        if gid:
            gid = "NCBI_Gene:"+gid
            gid = ncbi2hgnc.get(gid,gid)
    if r["relationship"].startswith("express"):
        reln = "expresses"
    elif r["relationship"] == "has_driver":
        reln = "is_regulated_by"
    if gid == "":
        gid = None
    return {
        "componentSymbol" : symbol,
        "componentID" : gid,
        "componentRelation" : reln
    }
    
def main () :
    submittedIds = loadSubmittedAlleles()
    aid2rels = {}
    for r in loadRelationship(EXPRESSES_cat_key):
        aid2rels.setdefault(r['allele'],[]).append(r)
    for r in loadRelationship(DRIVER_cat_key):
        aid2rels.setdefault(r['allele'],[]).append(r)
    aids = list(aid2rels.keys())
    #
    ncbi2hgnc = loadNcbi2HgncMapping()
    print('{\n  "metaData": %s,\n  "data": [' % json.dumps(buildMetaObject(MOUSEMINE), indent=2))
    first = True
    for aid in aids:
        if not aid in submittedIds:
            continue
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
        MGI_Relationship r
        LEFT JOIN ACC_Accession am 
            /* Have to outer join to get the MGI id since not all of these will be mouse genes. */
            ON r._object_key_2 = am._object_key
            AND am._mgitype_key = 2
            AND am._logicaldb_key = 1
            AND am.preferred = 1
        JOIN VOC_Term q
            ON r._qualifier_key = q._term_key
        JOIN VOC_Term e
            ON r._evidence_key = e._term_key
        JOIN VOC_Term rr
            ON r._relationshipterm_key = rr._term_key
        JOIN ACC_Accession aa
            ON r._object_key_1 = aa._object_key
            AND aa._mgitype_key = 11
            AND aa._logicaldb_key = 1
            AND aa.preferred = 1
        JOIN ALL_Allele al
            ON r._object_key_1 = al._allele_key
        JOIN MRK_Marker mm
            ON r._object_key_2 = mm._marker_key
    WHERE r._category_key = %d
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

# query for alleles being submitted to the Alliance. This is a MOUSEMINE query!
qAlleles = '''<query
      model="genomic"
      view="Allele.primaryIdentifier"
      constraintLogic="A and B and C and (D or E)"
      >
      <constraint code="A" path="Allele.organism.taxonId" op="=" value="10090" />
      <constraint code="B" path="Allele.alleleType" op="NONE OF">
        <value>QTL</value>
      </constraint>
      <constraint code="C" path="Allele.isWildType" op="=" value="false" />
      <constraint code="D" path="Allele.glTransmission" op="!=" value="Cell Line"/>
      <constraint code="E" path="Allele.glTransmission" op="IS NULL" />
      </query>
    '''
# query to return mapping from NCBI gene ids to HGNC ids for human genes.
# (The VAST majority of non-mouse constructs contain human genes.)
qNcbi2Hgnc = '''
    select a.accid as "NCBI", a2.accid as "HGNC"
    from acc_accession a, mrk_marker m, acc_accession a2
    where a._logicaldb_key = 55
    and a._object_key = m._marker_key
    and a._mgitype_key = 2
    and m._organism_key != 1
    and a2._object_key = m._marker_key
    and a2._logicaldb_key = 64
    and a2._mgitype_key = 2
    '''

# go!
main()

