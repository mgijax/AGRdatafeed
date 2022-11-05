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
