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
import json
from AGRlib import getConfig, stripNulls, buildMetaObject, makeOneOfConstraint, sql, doQuery

cp = getConfig()
MOUSEMINE = cp.get("DEFAULT","MOUSEMINEURL")

EXPRESSES_cat_key = "1004"
DRIVER_cat_key = "1006"

def log (s) :
    sys.stderr.write(s + '\n')

def loadSubmittedAlleles () :
    aids = set()
    for a in doQuery(qAlleles, MOUSEMINE):
        aids.add(a["primaryIdentifier"])
    return aids
    
ak2nmdId = {} # allele_key -> non-mouse-driver ID
def loadNonMouseDrivers () :
    for r in sql(qNonMouseDrivers):
        ak2nmdId[r['_allele_key']] = r['accid']

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
    else:
        raise RuntimeError("Internal error: unknown _category_key: " + str(r))

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
    print('{\n  "metaData": %s,\n  "data": [' % json.dumps(buildMetaObject(MOUSEMINE), indent=2))
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


# query MGI_Relationship for "expressed_component" (1004) and "has_driver" (1006)
qRelationships = '''
    SELECT 
        r._relationship_key,
        r._category_key,
        al._allele_key,
        aa.accid as allele,
        al.symbol as alleleSymbol,
        mm._organism_key,
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
    WHERE r._category_key = %s
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
    WHERE r._category_key = %s
    ORDER BY r._relationship_key, p.sequenceNum
    '''

# query for non-mouse drivers used in recombinase alleles
qNonMouseDrivers = '''
    SELECT a.accid, r._object_key_1 as _allele_key
    FROM
        MGI_Relationship r,
        MRK_Marker m,
        ACC_Accession a
    WHERE
        r._category_key = %s
    AND r._object_key_2 = m._marker_key
    AND m._organism_key != 1
    AND a._object_key = m._marker_key
    AND a._mgitype_key = 2
    AND a._logicaldb_key in (64,47,172) /* HGNC, RGD, ZFIN */
    ''' % DRIVER_cat_key 

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

# go!
main()

