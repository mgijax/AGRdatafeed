#!/usr/bin/env python2.7

# Create JSON file describing genotype-disease annotations for AGR data ingest.
# Usage:
#       python diseasePheno.py -d > MGI_0.6.2_diseaseAnnotations.json
# Pulls data from MouseMine.
# Writes JSON to stdout.
# See JSON specs at:
#    https://github.com/alliance-genome/agr_schemas/blob/development/disease/diseaseModelAnnotation.json
#    https://github.com/alliance-genome/agr_schemas/blob/release-1.0.0.3/phenotype/phenotypeModelAnnotation.json
# For debugging, accepts optional MGI IDs for genotypes/alleles/genes on command line
#   to restrict output to DAF records for those genotypes.
# Example genotype ID:  MGI:5526095 MGI:2175208 MGI:5588576
# Example gene IDs: MGI:97490 MGI:99607
# Example allele IDs: MGI:3603004 MGI:2680557
#       python diseasePheno.py -d MGI:5526095 MGI:2175208 MGI:5588576 MGI:99607 MGI:97490 MGI:97490 MGI:99607 > disease.sample.json
#       python diseasePheno.py -p MGI:5526095 MGI:2175208 MGI:5588576 MGI:99607 MGI:97490 MGI:97490 MGI:99607 > pheno.sample.json
#
# Original author: Jim Kadin
# Revisions: Joel Richardson
#

import sys
import argparse
import json
from itertools import imap, chain, groupby
from AGRlib import getConfig, stripNulls, buildMgiDataProviderObject, buildMetaObject, getTimeStamp, makeOneOfConstraint, doQuery
import heapq

# load config settings
cp = getConfig()

MOUSEMINE    = cp.get("DEFAULT","MOUSEMINEURL")
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
#   okind - kind of annotation ontology term, "DOTerm" or "MPTerm"
#   skind - kind of annotation subject, "SequenceFeature", "Allele", or "Genotype"
#   ids - either None, or a list of MGI ids. Annotations will be returned only for objects in this list
# Returns:
#   An iterator that will yield annotations. Annotations are objects - use dot notation to access parts,
#   e.g., a.subject.symbol
#
def annotations(url, okind, skind, ids = None):
    qopts = {
      'alleleFeatView': "OntologyAnnotation.subject.feature.primaryIdentifier" if skind == "Allele" else "",
      'xtraConstraint': makeOneOfConstraint("OntologyAnnotation.subject.primaryIdentifier",  ids),
      'xtraConstraint2': makeOneOfConstraint("OntologyAnnotationEvidence.annotation.subject.primaryIdentifier",  ids),
      'okind': okind,
      'skind': skind,
    }

    qAnnots = '''<query
    model="genomic"
    view="
	OntologyAnnotation.id
	OntologyAnnotation.subject.primaryIdentifier
	OntologyAnnotation.subject.symbol
	OntologyAnnotation.subject.name
	OntologyAnnotation.ontologyTerm.identifier
	OntologyAnnotation.ontologyTerm.name
	OntologyAnnotation.qualifier
	%(alleleFeatView)s
	"
    sortOrder="OntologyAnnotation.id asc"
    >
    <constraint path="OntologyAnnotation.ontologyTerm" type="%(okind)s"/>
    <constraint path="OntologyAnnotation.subject" type="%(skind)s"/>
    <constraint path="OntologyAnnotation.subject.organism.taxonId" op="=" value="10090"/>
    %(xtraConstraint)s
    </query>
    ''' % qopts

    qEvidence = '''<query
    model="genomic"
    view="
	OntologyAnnotation.id
	OntologyAnnotation.evidence.id
	OntologyAnnotation.evidence.annotationDate
	OntologyAnnotation.evidence.code.code
	OntologyAnnotation.evidence.publications.pubMedId
	OntologyAnnotation.evidence.publications.mgiJnum
	OntologyAnnotation.evidence.publications.mgiId
	%(alleleFeatView)s
	"
    sortOrder="OntologyAnnotation.id asc OntologyAnnotation.evidence.id asc"
    >
    <constraint path="OntologyAnnotation.ontologyTerm" type="%(okind)s"/>
    <constraint path="OntologyAnnotation.subject" type="%(skind)s"/>
    <constraint path="OntologyAnnotation.subject.organism.taxonId" op="=" value="10090"/>
    %(xtraConstraint)s
    </query>
    ''' % qopts

    qBaseAnnots = '''<query
	model="genomic"
	view="
	    OntologyAnnotationEvidence.id
	    OntologyAnnotationEvidence.annotation.id
	    OntologyAnnotationEvidence.baseAnnotations.subject.primaryIdentifier
	    OntologyAnnotationEvidence.baseAnnotations.evidence.annotationDate
	    "
	sortOrder="OntologyAnnotationEvidence.annotation.id asc OntologyAnnotationEvidence.id asc"
	>
	<constraint path="OntologyAnnotationEvidence.annotation.ontologyTerm" type="%(okind)s"/>
	<constraint path="OntologyAnnotationEvidence.annotation.subject" type="%(skind)s"/>
	%(xtraConstraint2)s
	</query>
    ''' % qopts

    qs = [
      imap(lambda x: (x[0], 'annotation', list(x[1])), groupby(doQuery(qAnnots, url), lambda e: e['id'])),
      imap(lambda x: (x[0], 'evidence', list(x[1])), groupby(doQuery(qEvidence, url), lambda e: e['id'])),
      imap(lambda x: (x[0], 'baseAnnots', list(x[1])), groupby(doQuery(qBaseAnnots, url), lambda e: e['annotation.id'])),
    ]
    for x in groupby(heapq.merge(*qs), lambda x: x[0]):
	r = {}
	for y in x[1]:
	    if y[1] == 'annotation':
		r.update(y[2][0])
	    elif y[1] == 'evidence':
		r['evidence'] = y[2]
	    elif y[1] == 'baseAnnots':
	        r['baseAnnots'] = y[2]
	rr = applyConversions(r, skind)
	if rr:
	  for e in rr["invevidence"]:
	      rr["agrevidence"] = e
	      yield formatDafJsonRecord(rr, "disease" if okind == "DOTerm" else "phenotype")

########
# Applies various transformations to a 'raw' annotation returned from the db to prepare it
# for export to AGR. Example: AGR dates must be in rfc3339 format. See comments for specific 
# transforms.
# Args:
#   a - one annotation
#   kind - "SequenceFeature", "Allele", or "Genotype"
# Returns:
#   The annotation, with conversions applied
#
geno2genes = {}   # genotype->rolled up genes index
gene2term = {} # gene id -> set of disease or phenotype ids
def applyConversions(a, kind):
    global geno2genes
    #
    # "not" is the only recognized qualifier for 1.0
    if not 'qualifier' in a:
      print '\n\nERROR:', a
    a["qualifier"] = "not" if a["qualifier"] == "NOT" else None
    #
    # MGI (and MM) store evidence as one evidence code with >= 1 ref.
    # AGR inverts this: one reference with >= 1 evidence code.
    ref2codes = {}
    for e in a["evidence"]:
	# FIXME: temporary tweak for MouseMine annotations. Remove once allele-annotations are being loaded from MGI
	if e["evidence.code.code"] == "DOA": e["evidence.code.code"] = "TAS"
	##
	pmid = e["evidence.publications.pubMedId"]
	mgiid = e["evidence.publications.mgiId"]
	ref2codes.setdefault((mgiid,pmid), set()).add(e["evidence.code.code"])
    a["invevidence"] = []
    for (k,es) in ref2codes.items():
        (mgiId, pubMedId) = k
        p = { "modPublicationId" : mgiId }
        if pubMedId: p["pubMedId"] = "PMID:" + pubMedId
        a["invevidence"].append({ "publication" : p, "evidenceCodes" : list(es) })
    #
    # object relation for genes and alleles
    if kind == "SequenceFeature":
        # Record for later use
        gene2term.setdefault(a["subject.primaryIdentifier"], set()).add(a["ontologyTerm.identifier"])
        #
        a["objectName"] = a["subject.symbol"]
        a["objectRelation"] = {
            "objectType" : "gene",
            "associationType" : "is_implicated_in"
        }
        setAnnotationDate(a, kind)
        for ba in a["baseAnnots"]:
	    # record for later use: which genotypes roll up to this gene
	    geno2genes.setdefault(ba["baseAnnotations.subject.primaryIdentifier"], set()).add(a["subject.primaryIdentifier"])
    elif kind == "Allele":
        # FIXME FIXME FIXME
        # Rolled up allele-disease/pheno annotations have never been implemented in MGI. The rules implemented for mousmine are
        # ancient and incomplete. The following is a TEMPORARY measure to filter the annotations to a more acceptible set.
	#
	# To wit:
	#
        # We will only include a rolled up allele-disease/pheno annotation if the allele's gene also has a rolled-up annotation to 
        # the same disease/pheno. (Sue Bello suggested this rule.) 
	#
        # The ultimate solution is to compute the rolled up annotations in MGI (as with gene-disease/pheno annotations) and simply
        # load them into MouseMine. There is a TR for this. When that TR is complete, the following can go away...
        # 
        if not a["ontologyTerm.identifier"] in gene2term.get(a["subject.feature.primaryIdentifier"], []):
            return None
        a["objectName"] = a["subject.symbol"]
        a["objectRelation"] = {
            "objectType" : "allele",
            "associationType" : "is_implicated_in"
        }
        setAnnotationDate(a, kind)
    else: # kind == "Genotype"
        a["objectName"] = a["subject.name"]
        a["objectRelation"] = {
            "objectType" : "genotype",
            "associationType" : "is_model_of",
            # The following line is the reason genes-disease annotation must be processed first...
            "inferredGeneAssociation": list(geno2genes.get(a["subject.primaryIdentifier"],[]))
        }
        setAnnotationDate(a, kind)
    return a

# Sets the annotation date for a gene or allele -to-disease annotation.
# This is the min annotation date from associated base (genotype) evidence recs
def setAnnotationDate(a, kind):
    d = None
    for e in a["evidence"]:
	if d is None or e["evidence.annotationDate"] < d:
	    d = e["evidence.annotationDate"]
    for ba in a.get("baseAnnots", []):
	bd = ba["baseAnnotations.evidence.annotationDate"]
	if d is None or bd < d:
	    d = bd
    a["annotationDate"] = getTimeStamp(d)

########
# Returns a JSON object for one annotation formatted according to the AGR disease or pheno annotation spec.
# Args:
#    annot
#    kind (string) Either "disease" or "phenotype"
#
def formatDafJsonRecord (annot, kind):
    if kind == "disease":
        return stripNulls({
            #'taxonId':                      MOUSETAXONID,
            'objectId':                     annot["subject.primaryIdentifier"],
            'objectName':                   annot["objectName"],
            'objectRelation':               annot["objectRelation"],
            #'experimentalConditions':      [],
            'qualifier':                    annot["qualifier"],
            'DOid':                         annot["ontologyTerm.identifier"],
            #'with':                        [],
            #'modifier':                    None,
            'evidence':                     annot["agrevidence"],
            #'geneticSex':                  '',
            'dateAssigned':                 annot["annotationDate"],
            'dataProvider':                 [ buildMgiDataProviderObject() ],
        })
    else:
        return stripNulls({
            'objectId':                     annot["subject.primaryIdentifier"],
            'phenotypeTermIdentifiers':     [{ "termId" : annot["ontologyTerm.identifier"], 'termOrder' : 1 }],
            'phenotypeStatement':           annot["ontologyTerm.name"],
	    'evidence': {
		'modPublicationId':         annot["agrevidence"]['publication'].get('modPublicationId', None),
		'pubMedId':                 annot["agrevidence"]['publication'].get('pubMedId', None),
	    },
            'dateAssigned':                 annot["annotationDate"],
            })

#####
def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "ids",
        nargs="*",
        help="Specific ids to generate annotation output for."
    )
    parser.add_argument(
        "-d", "--diseases",
        dest="doDiseases",
        action="store_true",
        default=False,
        help="Output disease annotations.")

    parser.add_argument(
        "-p", "--phenotypes",
        dest="doPhenotypes",
        action="store_true",
        default=False,
        help="Output phenotype annotations.")

    return parser.parse_args()

#####
def main():
    args = getArgs()
    service = MOUSEMINE
    mdo = buildMetaObject(service)
    print '{"metaData": %s,' % json.dumps(mdo) 
    print ' "data"    : ['
    if args.doDiseases:
        #
        # IMPORTANT! Gene annotations must be retrieved *before* Genotype annots. 
        geneAnnots = annotations(service, "DOTerm", "SequenceFeature", args.ids)
        alleleAnnots = annotations(service, "DOTerm", "Allele", args.ids)
        genotypeAnnots = annotations(service, "DOTerm", "Genotype", args.ids)
        #
        for i,ga in enumerate(chain(geneAnnots, alleleAnnots, genotypeAnnots)):
            print "," if i>0 else "", json.dumps(ga)
    if args.doPhenotypes:
        geneAnnots = annotations(service, "MPTerm", "SequenceFeature", args.ids)
        alleleAnnots = annotations(service, "MPTerm", "Allele", args.ids)
        #
        for i,ga in enumerate(chain(geneAnnots, alleleAnnots)):
            print "," if i>0 else "", json.dumps(ga)

    print "]}"

#####
main() 	# optional geno and/or gene IDs on the cmd line to restrict query
