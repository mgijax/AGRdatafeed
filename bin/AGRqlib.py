# Collect all the SQL in one place so that queries can be shared
# Naming convention:
# - fully formed queries start with "q"
# - templates (things with '%(name)s' substitutions) start with "t"


#--------------------------------------------------------------------------
# GENES
#--------------------------------------------------------------------------

#
qMcvTerms = '''
    SELECT _term_key, term, note
    FROM VOC_Term
    WHERE _vocab_key = 79 /* MCV */
    '''

# genes submitted to Alliance
qGenes = '''
    SELECT 
        m._marker_key, 
        a.accid as "markerId",
        m.symbol,
        m.name,
        n.note as description,
        m.chromosome,
        va._term_key as _mcv_term_key
    FROM 
        MRK_Marker m
        LEFT JOIN MRK_Notes n
            ON n._marker_key = m._marker_key
            AND m._organism_key = 1,
        ACC_Accession a,
        VOC_Annot va
    WHERE m._marker_status_key = 1  /* official */
    AND m._marker_type_key in (1,7) /* genes and pseudogenes */
    AND m._marker_key = a._object_key
    AND a._mgitype_key = 2
    AND a._logicaldb_key = 1
    AND a.preferred = 1
    AND m._marker_key = va._object_key
    AND va._annottype_key = 1011 /* Marker-MCV */
    '''

#
qGeneSynonyms = '''
    SELECT distinct _object_key as _marker_key, synonym
    FROM MGI_Synonym
    WHERE _synonymtype_key = 1004 /* exact marker synonyms */
    '''

#
qGeneSecondaryIds = '''
    SELECT a._object_key as _marker_key, a.accid
    FROM ACC_Accession a
    WHERE a._mgitype_key = 2
    AND a._logicaldb_key = 1
    AND a.preferred = 0
    '''

#
qGeneLocations = '''
    SELECT
        lc._marker_key,
        lc.chromosome,
        lc.genomicchromosome,
        lc.startcoordinate,
        lc.endcoordinate,
        lc.strand,
        lc.version as assembly
    FROM MRK_Location_Cache lc
    WHERE lc._organism_key = 1
    '''

#
qGeneXrefs = '''
    SELECT
      a._object_key as _marker_key,
      a.accid,
      ldb._logicaldb_key,
      CASE WHEN ldb._logicaldb_key = 9 THEN 'GenBank' ELSE ldb.name END AS "ldbName"
    FROM 
      MRK_Marker m,
      ACC_Accession a,
      ACC_LogicalDB ldb
    WHERE m._marker_status_key = 1 /* official */
    AND m._marker_type_key in (1,7) /* genes, pseudogenes */
    AND m._organism_key = 1
    AND m._marker_key = a._object_key
    AND a._mgitype_key = 2
    AND a._logicaldb_key not in (1,171,146,19,212) /* exclude MGI logical dbs */
    AND a._logicaldb_key = ldb._logicaldb_key
    '''

#
qGeneProteinIds = '''
    SELECT _object_key as _marker_key, accid as "proteinId"
    FROM acc_accession 
    WHERE _mgitype_key = 2
    AND _logicaldb_key in (13,41)
    '''

# genes with phenotype annots
qGeneHasPhenotype = '''
    SELECT distinct _object_key as _marker_key
    FROM VOC_Annot
    WHERE _annottype_key = 1015 /* MP-Gene */
    '''

# genes for alleles in the IMPC collection
qGeneHasImpc = '''
    SELECT distinct _marker_key
    FROM ALL_Allele
    WHERE _collection_key = 24755824 /* IMPC */
    '''

# genes that have expression data
qGeneHasExpression = '''
    SELECT distinct _marker_key
    FROM GXD_Expression
    '''

# genes that have expression data
qGeneHasExpressionImage = '''
    SELECT distinct _marker_key
    FROM GXD_Expression
    WHERE hasimage = 1
    '''

#--------------------------------------------------------------------------
# ALLELES
#--------------------------------------------------------------------------

# MGI ids of alleles being submitted to the Alliance. 
qSubmittedAlleleIds = '''
    SELECT aa.accid as mgiid
    FROM ALL_Allele a, ACC_Accession aa
    WHERE a._allele_key = aa._object_key
    AND aa._mgitype_key = 11
    AND aa._logicaldb_key = 1
    AND aa.preferred = 1
    AND a._transmission_key != 3982953 /* not Cell line */
    AND a._allele_type_key != 847130   /* not QTL */
    AND a._allele_status_key in (847114, 3983021) /* status= approved or autoload */
    AND a.iswildtype = 0               /* not wildtype */
    '''

# basic allele info
qAlleles = '''
    SELECT
      a._allele_key,
      aa.accid as "alleleId",
      a.symbol,
      a.name,
      t.term as "alleleType",
      n.note as "molecularNote",
      ma.accid as "markerId",
      mc.directterms as "markerType",
      d.symbol as "drivenBy"
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
    WHERE a._transmission_key != 3982953 /* not Cell line */
    AND a._allele_type_key != 847130 /* not QTL */
    AND a._allele_status_key in (847114, 3983021) /* approved, autoload */
    AND a.iswildtype = 0        /* not wild type */
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
#
qExpressors = '''
    SELECT distinct _object_key_1 as _allele_key
    FROM MGI_Relationship
    WHERE _category_key = 1004 /* expresses component */
    '''

#
qAlleleSynonyms = '''
    SELECT _object_key as _allele_key, synonym
    FROM MGI_Synonym
    WHERE _synonymtype_key = 1016 /* allele synonyms */
    '''

#--------------------------------------------------------------------------
# GENOTYPES
#--------------------------------------------------------------------------

# Return genotypes submitted to the alliance
qSubmittedGenotypes = '''
    SELECT
      g._genotype_key,
      a.accid as "genotypeId",
      n.note as alleles,
      s.strain,
      sa.accid as "backgroundId"
    FROM 
      GXD_Genotype g
      LEFT JOIN MGI_Note n
        ON g._genotype_key = n._object_key
        AND n._notetype_key = 1016 /* Combo type 1 */
      LEFT JOIN PRB_Strain s
        ON g._strain_key = s._strain_key,
      ACC_Accession a,
      ACC_Accession sa
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
    AND g._strain_key = sa._object_key
    AND sa._mgitype_key = 10 /* strain */
    AND sa._logicaldb_key = 1
    AND sa.preferred = 1
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

#--------------------------------------------------------------------------
# EXPRESSION
#--------------------------------------------------------------------------

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

#
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

#
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

#
qHomozygousGenotypes = '''
    SELECT _genotype_key
    FROM gxd_allelepair
    WHERE _Compound_key = 847167 and _PairState_key = 847138
    '''

#
qHeterozygousGenotypes = '''
    SELECT _genotype_key
    FROM gxd_allelepair
    WHERE _Compound_key = 847167 and _PairState_key = 847137
    '''

#
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

#--------------------------------------------------------------------------
# ANNOTATIONS
#--------------------------------------------------------------------------

#
# Params:
#  _annottype_key : 1015=MP/Marker, 1023=DO/Marker, 1029=DO/Allele, 1028=MP/Allele, 1002=MP/Genotype, 1020=DO/Genotype
#  subj_keycol : _marker_key , _allele_key , _genotype_key
#  subj_labelcol : symbol . Genotypes require special handling to generate their labels. For this query, use '_genotype_key'
#  subj_tblname : MRK_Marker , ALL_Allele , GXD_Genotype
#  _mgitype_key : 2=marker, 11=allele , 12=genotype
#  voc_ldbkey : 34=MP, 191=DO
tAnnots = '''
   SELECT
     va._annot_key,
     CASE WHEN qt.term is null THEN '' ELSE qt.term END as qualifier,
     aa.accid as "subjectId",
     subj.%(subj_keycol)s as "subjectKey",
     subj.%(subj_labelcol)s as "subjectLabel",
     vt.term,
     vta.accid as "termId"
   FROM
     VOC_Annot va,
     VOC_Term qt,
     %(subj_tblname)s subj, /* MRK_Marker or ALL_Allele */
     ACC_Accession aa,
     VOC_Term vt,
     ACC_Accession vta
   WHERE va._annottype_key = %(_annottype_key)d
   AND va._qualifier_key = qt._term_key
   AND va._object_key = subj.%(subj_keycol)s
   AND va._object_key = aa._object_key
   AND aa._mgitype_key = %(_mgitype_key)d
   AND aa._logicaldb_key = 1
   AND aa.preferred = 1
   AND va._term_key = vt._term_key
   AND vt._term_key = vta._object_key
   AND vta._mgitype_key = 13
   AND vta._logicaldb_key = %(voc_ldbkey)d
   AND vta.preferred = 1
   ORDER BY va._annot_key
   '''
# Params:
#  _annottype_key: one of 1015, 1023, 1028, 1029
tAnnotEvidence = '''
    SELECT
      va._annot_key,
      ve._annotevidence_key,
      et.abbreviation as code,
      br._refs_key,
      ra.accid as "refMgiId",
      raj.accid as "refJnum",
      rap.accid as "refPmid",
      ve.modification_date as "annotationDate"
    FROM
      VOC_Annot va,
      VOC_Evidence ve,
      VOC_Term et,
      BIB_Refs br,
      ACC_Accession ra,
      ACC_Accession raj
         LEFT JOIN ACC_Accession rap
         ON raj._object_key = rap._object_key
         AND rap._mgitype_key = 1
         AND rap._logicaldb_key = 29 /* pubmed */
         AND rap.preferred = 1
    WHERE va._annottype_key = %(_annottype_key)d
    AND va._annot_key = ve._annot_key
    AND ve._evidenceterm_key = et._term_key
    AND ve._refs_key = br._refs_key
    /* ref's MGI id */
    AND ve._refs_key = ra._object_key
    AND ra._mgitype_key = 1
    AND ra._logicaldb_key = 1
    AND ra.prefixpart = 'MGI:'
    AND ra.preferred = 1
    /* ref's J# id */
    AND ve._refs_key = raj._object_key
    AND raj._mgitype_key = 1
    AND raj._logicaldb_key = 1
    AND raj.prefixpart = 'J:'
    AND raj.preferred = 1
    ORDER BY va._annot_key, br._refs_key
    '''

# Params:
#  _annottype_key : one of 1015, 1023, 1028, 1029
tAnnotBaseAnnots = '''
    select distinct
      va._annot_key,
      ve._annotevidence_key,
      ve._refs_key,
      ba._annot_key as _baseannot_key,
      be._annotevidence_key as _baseannotevidence_key,
      be._refs_key as _baserefs_key,
      be.modification_date as _basemodification_date,
      aa.accid as "genotypeId"
    from 
      voc_annot va,
      voc_evidence ve,
      voc_evidence_property vep,
      voc_term vept,
      voc_annot ba,
      voc_evidence be,
      acc_accession aa
    where va._annot_key = ve._annot_key
    and ve._annotevidence_key = vep._annotevidence_key 
    and vep._propertyterm_key = vept._term_key
    and vept.term = '_SourceAnnot_key'
    and va._annottype_key = %(_annottype_key)d
    and cast(vep.value as integer) = ba._annot_key
    and ba._annottype_key = %(_baseannottype_key)d
    and ba._annot_key = be._annot_key
    and ba._object_key = aa._object_key
    and aa._mgitype_key = 12
    and aa._logicaldb_key = 1
    and aa.preferred = 1
    and ve._refs_key = be._refs_key
    order by va._annot_key
    '''

# Params:
#   _annottype_key: one of 1002, 1020
tGenotypeLabels = '''
    SELECT
      g._genotype_key,
      n.note as alleles,
      s.strain
    FROM 
      GXD_Genotype g
          LEFT JOIN MGI_Note n
          ON g._genotype_key = n._object_key
          AND n._notetype_key = 1016 /* Combo type 1 */
      
          LEFT JOIN PRB_Strain s
          ON g._strain_key = s._strain_key
    WHERE exists (
      /* the genotyeps has an annotation */
      SELECT _annot_key
      FROM VOC_Annot
      WHERE _annottype_key = %(_annottype_key)d
      AND _object_key = g._genotype_key
      )
    '''

#--------------------------------------------------------------------------
# High-Throughput Expression Metadata
#--------------------------------------------------------------------------

#-----------------------------------
qHTExperiments = '''
    SELECT 
      a.accid as "experimentId",
      name, 
      et.term as "experimentType", 
      e.description, 
      st.term as "studyType", 
      src.term as source, 
      e.last_curated_date as "curationDate"
    FROM 
       GXD_HTExperiment e,
       VOC_Term et,
       VOC_Term st,
       VOC_Term src,
       ACC_Accession a
    WHERE e._curationstate_key = 20475421 /* Done */
    AND e._experimenttype_key =  et._term_key
    AND e._studytype_key = st._term_key
    AND e._source_key = src._term_key
    AND e._experiment_key = a._object_key
    AND a._mgitype_key = 42 /* ht experiment type */
    AND a.preferred = 1
    '''
#-----------------------------------
qHTSamples = '''
    SELECT 
      ea.accid as "experimentId", 
      et.term as "experimentType",
      s.name, 
      s.age, 
      s.agemin, 
      s.agemax, 
      st.term as sex, 
      s._stage_key as stage, 
      em.term as "structureName", 
      ema.accid as "structureId",
      ga.accid as "genotypeId"
    FROM 
      GXD_HTSample s, 
      GXD_HTExperiment e,
      ACC_Accession ea,
      VOC_Term st,
      VOC_Term em,
      ACC_Accession ema,
      ACC_Accession ga,
      VOC_Term et
    WHERE e._experiment_key = s._experiment_key
    AND e._experimenttype_key = et._term_key
    AND s._relevance_key = 20475450 /* Yes */
    AND s._experiment_key = ea._object_key
    AND ea._mgitype_key = 42 /* ht experiment type */
    AND ea.preferred = 1
    AND s._sex_key = st._term_key
    AND s._emapa_key = em._term_key
    AND s._emapa_key = ema._object_key
    AND ema._mgitype_key = 13
    AND ema._logicaldb_key = 169 /* EMAPA */
    AND ema.preferred = 1
    AND s._genotype_key = ga._object_key
    AND ga._mgitype_key = 12 
    AND ga._logicaldb_key = 1
    AND ga.preferred = 1
    '''

#-----------------------------------
qHTPmids = '''
    SELECT 
      p.value as pmid,
      a.accid as "experimentId"
    FROM
      GXD_HTExperiment e,
      MGI_Property p,
      ACC_Accession a
    WHERE p._propertytype_key = 1002
    AND p._object_key = e._experiment_key
    AND p._propertyterm_key = 20475430 /* PMID */
    AND e._experiment_key = a._object_key
    AND a._mgitype_key = 42 /* ht experiment type */
    AND a.preferred = 1
    '''

qHTPmid2Mgi = '''
    SELECT 
      a.accid as pmid, a2.accid as mgiid
    FROM
      ACC_Accession a,
      ACC_Accession a2
    WHERE a._logicaldb_key = 29 /* Pubmed */
    AND a._mgitype_key = 1
    AND a.preferred = 1
    AND a._object_key = a2._object_key
    AND a2._mgitype_key = 1
    AND a2._logicaldb_key = 1
    AND a2.prefixpart = 'MGI:'
    AND a2.preferred = 1
    '''

#-----------------------------------
qHTVariables = '''
    SELECT
      a.accid as "experimentId", t.term as variable
    FROM
      GXD_HtExperimentVariable v,
      ACC_Accession a,
      VOC_Term t
    WHERE v._experiment_key = a._object_key 
    AND a._mgitype_key = 42 /* ht experiment type */
    AND a.preferred = 1
    AND v._term_key = t._term_key
    AND t.term != 'Not Applicable'
    '''
