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
#       python diseaseAnnotations.py MGI:5526095 MGI:2175208 MGI:5588576 MGI:99607 MGI:97490 > sample.json
#
# Original author: Jim Kadin
# Revisions: Joel Richardson
#

import sys
import json
from ConfigParser import ConfigParser
from intermine.webservice import Service
from AGRlib import stripNulls, buildMetaObject, getTimeStamp

##### Load config
cp = ConfigParser()
cp.optionxform = str
cp.read("config.cfg")

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
          "qualifier",
          "evidence.annotationDate",
          "evidence.code.code",
          "evidence.publications.pubMedId",
          "evidence.publications.mgiJnum",
          "evidence.publications.mgiId"
      )
      query.add_constraint("subject.organism.taxonId", "=", "10090")
      #
      if kind == "Gene":
          query.add_constraint("evidence.baseAnnotations.subject", "Genotype")
          query.add_view(
            "evidence.baseAnnotations.subject.symbol",
            "evidence.baseAnnotations.subject.background.name",
            "evidence.baseAnnotations.subject.primaryIdentifier",
            "evidence.baseAnnotations.evidence.annotationDate"
          )
          query.outerjoin("evidence.baseAnnotations")
      #
      if ids and len(ids):
          query.add_constraint("subject.primaryIdentifier", "ONE OF", ids)
      #
      query.add_sort_order("subject.primaryIdentifier", "ASC")
      query.add_sort_order("ontologyTerm.identifier", "ASC")
      #
      for a in query:
          applyConversions(a, kind)
          for e in a.invevidence:
            a.agrevidence = e
            yield formatDafJsonRecord(a)

########
# Applies various transformations to the 'raw' annotations returned by from the db to prepare it
# for export to AGR. Example: AGR dates must be in rfc3339 format. See comments for specific 
# transforms.
# Args:
#   a - one annotation
#   kind - either "SequenceFeature" or "Genotype"
# Returns:
#   The annotation, with conversions applied
#
geno2genes = {}   # genotype->rolled up genes index
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
            ref2codes.setdefault((p.mgiId,p.pubMedId), set()).add(e.code.code)
    a.invevidence = []
    for (k,es) in ref2codes.items():
        (mgiId, pubMedId) = k
        p = { "modPublicationId" : mgiId }
        if pubMedId: p["pubMedId"] = "PMID:" + pubMedId
        a.invevidence.append({ "publication" : p, "evidenceCodes" : list(es) })
    #
    # object relation for genes and alleles
    if kind == "SequenceFeature" or kind == "Allele":
        a.objectName = a.subject.symbol
        a.objectRelation = {
            "objectType" : "gene" if kind == "SequenceFeature" else "allele",
            "associationType" : "is_implicated_in"
        }
        #
        # annotationDate for gene is min annotation date from associated base (genotype) evidence recs
        d = None
        for e in a.evidence:
            for ba in e.baseAnnotations:
                geno2genes.setdefault(ba.subject.primaryIdentifier, set()).add(a.subject.primaryIdentifier)
                for be in ba.evidence:
                  if d is None or be.annotationDate < d:
                      d = be.annotationDate
        a.annotationDate = getTimeStamp(d)
    else: # kind == "Genotype"
        a.objectName = a.subject.name
        a.objectRelation = {
            "objectType" : "genotype",
            "associationType" : "is_model_of",
            "inferredGeneAssociation": list(geno2genes.get(a.subject.primaryIdentifier,[]))
        }
        #
        # annotationDate for genotype is min annotation date from associated evidence recs
        d = None
        for e in a.evidence:
            if d is None or e.annotationDate < d:
                d = e.annotationDate
        a.annotationDate = getTimeStamp(d)
    return a

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
    # IMPORTANT! Gene annotations must be retrieved *before* Genotype annots. 
    geneAnnots = list(annotations(service, "SequenceFeature", ids))
    alleleAnnots = []
    if DO_ALLELES: alleleAnnots = list(annotations(service, "Allele", ids))
    genoAnnots = list(annotations(service, "Genotype", ids))
    annots = (geneAnnots if DO_GENES else []) + alleleAnnots + genoAnnots
    jobj = {
      "metaData" : buildMetaObject(service),
      "data"     : annots
      }
    print json.dumps(jobj, sort_keys=True, indent=2, separators=(',', ': ')),

#####
main(sys.argv[1:]) 	# optional geno and/or gene IDs on the cmd line to restrict query
