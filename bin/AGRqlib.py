# Collect the SQL in one place so that queries can be shared

# query for MGI ids of alleles being submitted to the Alliance. 
qSubmittedAlleleIds = '''
    SELECT aa.accid as mgiid
    FROM ALL_Allele a, ACC_Accession aa
    WHERE a._allele_key = aa._object_key
    AND aa._mgitype_key = 11
    AND aa._logicaldb_key = 1
    AND aa.preferred = 1
    AND a._transmission_key != 3982953 /* Cell line */
    AND a._allele_type_key != 847130 /* QTL */
    AND a._allele_status_key in (847114, 3983021) /* approved, autoload */
    AND a.iswildtype = 0
    '''

qAlleles = '''
    SELECT
      aa.accid as "alleleId",
      a.symbol,
      a.name,
      t.term as "alleleType",
      n.note,
      ma.accid as "markerId",
      d.symbol as "drivenBy",
      mc.directterms as "markerType"
    FROM
      ALL_Allele a
        LEFT JOIN MGI_Note n
          ON a._allele_key = n._object_key
          AND n._notetype_key = 1021 /* molecular */
        LEFT JOIN MGI_Relationship r
          ON r._category_key = 1006 /* allele-driver */
          AND a._allele_key = r._object_key_1
          LEFT JOIN MRK_Marker d
            ON r._object_key_2 = d._marker_key,
      ACC_Accession aa,
      VOC_Term t,
      ACC_Accession ma,
      MRK_MCV_Cache mc
    WHERE a._transmission_key != 3982953 /* Cell line */
    AND a._allele_type_key != 847130 /* QTL */
    AND a._allele_status_key in (847114, 3983021) /* approved, autoload */
    AND a.iswildtype = 0
    AND a._allele_key = aa._object_key
    AND aa._mgitype_key = 11
    AND aa._logicaldb_key = 1
    AND aa.preferred = 1
    AND a._allele_type_key = t._term_key
    AND a._marker_key = ma._object_key
    AND ma._mgitype_key = 2
    AND ma._logicaldb_key = 1
    AND ma.preferred = 1
    AND a._marker_key = mc._marker_key
    AND mc.qualifier = 'D'
   '''

# Return genotypes submitted to the alliance
qSubmittedGenotypes = '''
    SELECT
      a.accid as "genotypeId",
      n.note as alleles,
      s.strain,
      sa.accid as "backgroundId"
    FROM 
      GXD_Genotype g,
      ACC_Accession a,
      MGI_Note n,
      ACC_Accession sa,
      PRB_Strain s
    WHERE exists (
      /* the genotyeps has an annotation */
      SELECT _annot_key
      FROM VOC_Annot
      WHERE _annottype_key in (1002,1020) /* MP-geno, DO-geno */
      AND _object_key = g._genotype_key
      )
    AND g._genotype_key = a._object_key
    AND a._mgitype_key = 12
    AND a._logicaldb_key = 1
    AND a.preferred = 1
    AND g._genotype_key = n._object_key
    AND n._notetype_key = 1016 /* Combo type 1 */
    AND g._strain_key = sa._object_key
    AND sa._mgitype_key = 10 /* strain */
    AND sa._logicaldb_key = 1
    AND sa.preferred = 1
    AND g._strain_key = s._strain_key
    '''

# Return the first allele and zygosity from each allelepair record.
qGenotypeAllelePair = '''
    SELECT
      ga.accid as "genotypeId",
      aa.accid as "alleleId",
      ps.term as "pairState" 
    FROM 
      GXD_AllelePair ap,
      ACC_Accession ga,
      ALL_Allele a1,
      ACC_Accession aa,
      VOC_Term ps
    WHERE ap._genotype_key = ga._object_key
    AND ga._mgitype_key = 12
    AND ga._logicaldb_key = 1
    AND ga.preferred = 1
    AND ap._allele_key_1 = a1._allele_key
    AND a1._allele_key = aa._object_key
    AND aa._mgitype_key = 11
    AND aa._logicaldb_key = 1
    AND aa.preferred = 1
    AND ap._pairstate_key = ps._term_key
    '''

# Returns all the terms in the EMAPA.
qEmapaTerms = '''
    SELECT aa.accid, vt.term, vte.startstage, vte.endstage
    FROM VOC_Term_Emapa vte, VOC_Term vt, ACC_Accession aa
    WHERE vt._term_key = vte._term_key
    AND vt._term_key = aa._object_key
    AND aa._mgitype_key = 13
    AND aa._logicaldb_key = 169
    AND aa.preferred = 1
    '''

# Returns ids of EMAPA terms and their immediate parents.
qEmapaTermsAndParents = '''
    SELECT 
      ca.accid as childid,
      pa.accid as parentid
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

# old mousemine expression query for comparison
'''<query
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
    </query>
  '''

qMutantGenotypes = '''
    SELECT g._genotype_key
    FROM GXD_Genotype g
    WHERE EXISTS(
      SELECT 1
      FROM gxd_allelegenotype ag, all_allele a
      WHERE g._genotype_key = ag._genotype_key
      AND ag._allele_key = a._allele_key 
      AND a.iswildtype = 0
    )
    AND g._genotype_key > 0
    '''

qNonMutantGenotypes = '''
    SELECT g._genotype_key
    FROM GXD_Genotype g
    WHERE NOT EXISTS(
      SELECT 1
      FROM gxd_allelegenotype ag, all_allele a
      WHERE g._genotype_key = ag._genotype_key
      AND ag._allele_key = a._allele_key 
      AND a.iswildtype = 0
    )
    /* AND g._genotype_key > 0 */
    '''

qHomozygousGenotypes = '''
    SELECT _genotype_key
    FROM gxd_allelepair
    WHERE _Compound_key = 847167 and _PairState_key = 847138
    '''

qHeterozygousGenotypes = '''
    SELECT _genotype_key
    FROM gxd_allelepair
    WHERE _Compound_key = 847167 and _PairState_key = 847137
    '''

qGxdExpression = '''
    SELECT
      a.accid as "assayId",
      t.assaytype as "assayType",
      fa.accid as "geneId",
      ex._stage_key as stage,
      sa.accid as "structureId",
      ra.accid as "refMgiId",
      pa.accid as "refPubmedId"
    FROM
      GXD_Expression ex,
      ACC_Accession a,
      GXD_AssayType t,
      ACC_Accession fa,
      ACC_Accession sa,
      ACC_Accession ra 
        LEFT JOIN ACC_Accession pa 
        ON ra._object_key = pa._object_key
        AND pa._mgitype_key = 1
        AND pa._logicaldb_key = 29 /* pubmed */
    WHERE ex.isforgxd = 1
    /* detected */
    AND ex.expressed = 1
    AND (
      /* non-mutant genotypes */
      ex._genotype_key in (%s)
      OR (
        ex._assaytype_key = 9 /* in situ reporter (knockin) */
        AND
        ex._genotype_key IN (%s))) /* heterozygote */
    /* assayId */
    AND ex._assay_key = a._object_key
    AND a._mgitype_key = 8 /* GXD assay */
    AND a._logicaldb_key = 1
    AND a.preferred = 1
    /* assayType */
    AND ex._assaytype_key = t._assaytype_key
    AND ex._marker_key = fa._object_key
    AND fa._mgitype_key = 2
    AND fa._logicaldb_key = 1
    AND fa.preferred = 1
    /* structureId */
    AND ex._emapa_term_key = sa._object_key
    AND sa._mgitype_key = 13
    AND sa._logicaldb_key = 169 /* emapa */
    AND sa.preferred = 1
    /* refMgiId */
    AND ex._refs_key = ra._object_key
    AND ra._mgitype_key = 1
    AND ra._logicaldb_key = 1
    AND ra.prefixpart = 'MGI:'
    AND ra.preferred = 1
    /**/
    ORDER BY a.accid, sa.accid, ex._stage_key
    ''' % (qNonMutantGenotypes,qHeterozygousGenotypes)
