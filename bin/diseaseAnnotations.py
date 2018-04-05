#!/usr/bin/env python2.7

# Create JSON file describing genotype-disease annotations for AGR data ingest.
# Usage:
#       python diseaseAnnotations.py > MGI_0.6.2_diseaseAnnotations.json
# Pulls data from MouseMine.
# Writes JSON to stdout.
# See JSON spec at:
#https://github.com/alliance-genome/agr_schemas/blob/development/disease/diseaseModelAnnotation.json
# (this is development branch, might want to see master branch too)
# For debugging, accepts optional MGI IDs for genotypes on command line
#   to restrict output to DAF records for those genotypes.
# Example genotype ID:  MGI:5526095 MGI:2175208 MGI:5588576
# Example gene IDs: MGI:97490 MGI:99607
# Example allele IDs: MGI:3603004 MGI:2680557
#       python diseaseAnnotations.py MGI:5526095 MGI:2175208 MGI:5588576 MGI:99607 MGI:97490 MGI:97490 MGI:99607 > sample.json
#
# Original author: Jim Kadin
# Revisions: Joel Richardson
#

import sys
import json
import itertools
from intermine.webservice import Service
from AGRlib import getConfig, stripNulls, buildMetaObject, getTimeStamp

# load config settings
cp = getConfig()

MOUSEMINEURL    = cp.get("DEFAULT","MOUSEMINEURL")
TAXONID         = cp.get("DEFAULT","TAXONID")
MOUSETAXONID    = cp.get("DEFAULT","GLOBALTAXONID")
DO_GENES        = cp.getboolean("dafFile","DO_GENES")
DO_ALLELES      = cp.getboolean("dafFile","DO_ALLELES")
DO_GENOS        = cp.getboolean("dafFile","DO_GENOS")

########
# Returns an iterator over the disease annotation for sequence features or genotypes, depending
# on the kind argument. Optionally restrict to a specific set of ids.
# Args:
#   service - a connection to MouseMine
#   kind - either "SequenceFeature" or "Genotype"
#   ids - either None, or a list of MGI ids. Annotations will be returned only for objects in this list
# Returns:
#   An iterator that will yield annotations. Annotations are objects - use dot notation to access parts,
#   e.g., a.subject.symbol
#
def annotations(service, kind, ids = None):
      query = service.new_query("OntologyAnnotation")
      #
      query.add_constraint("subject", kind)
      query.add_constraint("ontologyTerm", "DOTerm")
      #
      query.add_view(
          "subject.primaryIdentifier",
          "subject.symbol",
          "subject.name",
          "ontologyTerm.identifier",
          "ontologyTerm.name",
          "evidence.annotationDate",
          "evidence.code.code",
          "evidence.publications.pubMedId",
          "evidence.publications.mgiJnum",
          "evidence.publications.mgiId"
      )
      if kind == "Allele":
          query.add_view("subject.feature.primaryIdentifier")

      query.add_constraint("subject.organism.taxonId", "=", "10090")
      #
      if kind == "SequenceFeature": # FIXME: Alleles should also have baseAnnotations. 
          query.add_constraint("evidence.baseAnnotations.subject", "Genotype")
          query.add_view(
            "evidence.baseAnnotations.subject.primaryIdentifier",
            "evidence.baseAnnotations.evidence.annotationDate"
          )
          query.outerjoin("evidence.baseAnnotations")
      #
      if kind == "xAllele": # disabled
          query.add_constraint("subject.feature.ontologyAnnotations.ontologyTerm", "DOTerm")
          query.add_constraint("ontologyTerm", "IS", "OntologyAnnotation.subject.feature.ontologyAnnotations.ontologyTerm")
      #
      if ids and len(ids):
          query.add_constraint("subject.primaryIdentifier", "ONE OF", ids)
      #
      query.add_sort_order("subject.primaryIdentifier", "ASC")
      query.add_sort_order("ontologyTerm.identifier", "ASC")
      #
      for a in query:
          if not applyConversions(a, kind):
              continue
          for e in a.invevidence:
            a.agrevidence = e
            yield formatDafJsonRecord(a)

########
# Applies various transformations to a 'raw' annotation returned from the db to prepare it
# for export to AGR. Example: AGR dates must be in rfc3339 format. See comments for specific 
# transforms.
# Args:
#   a - one annotation
#   kind - either "SequenceFeature" or "Genotype"
# Returns:
#   The annotation, with conversions applied
#
geno2genes = {}   # genotype->rolled up genes index
gene2disease = {} # gene id -> set of disease ids
def applyConversions(a, kind):
    global geno2genes
    #
    # "not" is the only recognized qualifier for 1.0
    a.qualifier = "not" if a.qualifier == "NOT" else None
    #
    # MGI (and MM) store evidence as one evidence code with >= 1 ref.
    # AGR inverts this: one reference with >= 1 evidence code.
    ref2codes = {}
    for e in a.evidence:
        for p in e.publications:
            # FIXME: temporary tweak for MouseMine annotations. Remove once allele-annotations are being loaded from MGI
            if e.code.code == "DOA": e.code.code = "TAS"
            ##
            ref2codes.setdefault((p.mgiId,p.pubMedId), set()).add(e.code.code)
    a.invevidence = []
    for (k,es) in ref2codes.items():
        (mgiId, pubMedId) = k
        p = { "modPublicationId" : mgiId }
        if pubMedId: p["pubMedId"] = "PMID:" + pubMedId
        a.invevidence.append({ "publication" : p, "evidenceCodes" : list(es) })
    #
    # object relation for genes and alleles
    if kind == "SequenceFeature":
        # Record for later use
        gene2disease.setdefault(a.subject.primaryIdentifier,set()).add(a.ontologyTerm.identifier)
        #
        a.objectName = a.subject.symbol
        a.objectRelation = {
            "objectType" : "gene",
            "associationType" : "is_implicated_in"
        }
        setAnnotationDate(a, kind)
        for e in a.evidence:
            for ba in e.baseAnnotations:
                # record for later use: which genotypes roll up to this gene
                geno2genes.setdefault(ba.subject.primaryIdentifier, set()).add(a.subject.primaryIdentifier)
    elif kind == "Allele":
        # FIXME FIXME FIXME
        # Rolled up allele-disease annotations have never been implemented in MGI. The rules implemented for mousmine are
        # ancient and incomplete. The following is a TEMPORARY measure to filter the annotations to a more acceptible set.
        # We will only include a rolled up allele-disease annotation if the allele's gene has a rolled-up annotation to 
        # the same disease. (Sue Bello suggested this rule.) 
        # The ultimate solution is to compute the rolled up annotations in MGI (as with gene-disease annotations) and simply
        # load them into MouseMine. There is a TR for this. When that TR is complete, the following can go away...
        # 
        if not a.ontologyTerm.identifier in gene2disease.get(a.subject.feature.primaryIdentifier, []):
            return None
        a.objectName = a.subject.symbol
        a.objectRelation = {
            "objectType" : "allele",
            "associationType" : "is_implicated_in"
        }
        setAnnotationDate(a, kind)
    else: # kind == "Genotype"
        a.objectName = a.subject.name
        a.objectRelation = {
            "objectType" : "genotype",
            "associationType" : "is_model_of",
            # The following line is the reason genes-disease annotation must be processed first...
            "inferredGeneAssociation": list(geno2genes.get(a.subject.primaryIdentifier,[]))
        }
        setAnnotationDate(a, kind)
    return a

# Sets the annotation date for a gene or allele -to-disease annotation.
# This is the min annotation date from associated base (genotype) evidence recs
def setAnnotationDate(a, kind):
    d = None
    for e in a.evidence:
        if kind == "Genotype" or kind == "Allele":  # FIXME: allele should have base annots
            if d is None or e.annotationDate < d:
                d = e.annotationDate
        else:
            for ba in e.baseAnnotations:
                for be in ba.evidence:
                  if d is None or be.annotationDate < d:
                      d = be.annotationDate
    a.annotationDate = getTimeStamp(d)

########
# Returns a JSON object for one annotation formatted according to the AGR disease annotation spec.
#
def formatDafJsonRecord (annot):
    return stripNulls({
        'taxonId':                      MOUSETAXONID,
        'objectId':                     annot.subject.primaryIdentifier,
        'objectName':                   annot.objectName,
        'objectRelation':               annot.objectRelation,
        #'experimentalConditions':      [],
        'qualifier':                    annot.qualifier,
        'DOid':                         annot.ontologyTerm.identifier,
        #'with':                        [],
        #'modifier':                    None,
        'evidence':                     annot.agrevidence,
        #'geneticSex':                  '',
        'dateAssigned':                 annot.annotationDate,
        'dataProvider':                 'MGI',
    })

#####
def main(ids):
    service = Service(MOUSEMINEURL)

    mdo = buildMetaObject(service)
    print '{"metaData": %s,' % json.dumps(mdo) 
    print ' "data"    : ['
    #
    # IMPORTANT! Gene annotations must be retrieved *before* Genotype annots. 
    geneAnnots = annotations(service, "SequenceFeature", ids)
    alleleAnnots = annotations(service, "Allele", ids)
    genotypeAnnots = annotations(service, "Genotype", ids)
    #
    for i,ga in enumerate(itertools.chain(geneAnnots, alleleAnnots, genotypeAnnots)):
        print "," if i>0 else "", json.dumps(ga)

    print "]}"

#####
main(sys.argv[1:]) 	# optional geno and/or gene IDs on the cmd line to restrict query
