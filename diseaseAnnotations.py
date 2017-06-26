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
#       python diseaseAnnotations.py MGI:5526095 MGI:2175208 MGI:5588576 > sample.json
#
# Author: jim kadin

import sys
import json
from ConfigParser import ConfigParser
from intermine.webservice import Service
from AGRlib import AGRjsonFormatter, buildMetaObject

##### Load config
cp = ConfigParser()
cp.optionxform = str
cp.read("config.cfg")

MOUSEMINEURL    = cp.get("DEFAULT","MOUSEMINEURL")
TAXONID         = cp.get("DEFAULT","TAXONID")
GLOBALTAXONID   = cp.get("DEFAULT","GLOBALTAXONID")
DO_GENES        = cp.getboolean("dafFile","DO_GENES")
DO_GENOS        = cp.getboolean("dafFile","DO_GENOS")
RFC3339TIME	= "T10:00:00-05:00"	# add to dates to make RFC3339 date/time

########
class GenoToGeneFinder (object):
    # Builds a mapping from genotypes to rolled-up gene(s).
    def __init__(self, service):
	query = self.buildQuery(service)
	self.buildMapping(query)
    ######

    def genoToRolledUpGenes(self, genoID):
	# return list of (geneID, symbol) pairs that roll up to the genoID
	return self.genoToGene.get(genoID, {}).keys()
    ######

    def buildMapping(self, query):

	self.genoToGene = {}	# geneoToGene[genoID] == {(geneID:,symbol):1}
				# There can be multiple rolled up genes:
				#  1 example: gene + transgene if the TG
				#   expresses the gene.
	for row in query.rows():
	    genoID = row["ontologyAnnotations.evidence.baseAnnotations.subject.primaryIdentifier"]
	    geneInfo = (row["primaryIdentifier"], row["symbol"])
	    self.genoToGene.setdefault(genoID,{})[geneInfo] = 1

    ######
    def buildQuery(self,service):
	query = service.new_query("SequenceFeature")

	# Type constraints should come early 
	query.add_constraint("ontologyAnnotations.ontologyTerm", "DiseaseTerm")
	query.add_constraint("ontologyAnnotations.evidence.baseAnnotations.subject", "Genotype")

	# The view specifies the output columns
	query.add_view(
	"primaryIdentifier", "symbol",
	"ontologyAnnotations.ontologyTerm.identifier",
	"ontologyAnnotations.ontologyTerm.name", "ontologyAnnotations.qualifier",
	"ontologyAnnotations.evidence.baseAnnotations.subject.symbol",
	"ontologyAnnotations.evidence.baseAnnotations.subject.background.name",
	"ontologyAnnotations.evidence.baseAnnotations.subject.primaryIdentifier"
	)

	query.add_constraint("organism.taxonId", "=", "10090", code = "C")
	return query
    ######
####### end class GenoToGeneFinder

class GenoAnnotationQuery (object): 
    # Provide an iterator, annotations(), for all genotype disease annotations
    #   pulled from MouseMine.
    # Each returned annotation record is a dictionary with fields defined
    #    below.

    def __init__(self, service, ids):
	# ids is list of Genotype MGI ids to restrict the query to.
	#   if empty list, get them all
	self.service = service
	self.ids = ids

    def buildGenoAnnotQuery(self):

	# Get a new query on the class (table) you will be querying:
	query = self.service.new_query("Genotype")

	#Type constraints come before all mentions of the paths they constrain
	query.add_constraint("ontologyAnnotations.ontologyTerm", "DOTerm")
	if len(self.ids):
	    query.add_constraint("primaryIdentifier", "ONE OF", self.ids)

	# The view specifies the output columns
	query.add_view(
	    "primaryIdentifier", "name",
	    "ontologyAnnotations.ontologyTerm.identifier",
	    "ontologyAnnotations.ontologyTerm.name",
	    "ontologyAnnotations.qualifier",
	    "ontologyAnnotations.evidence.code.code",
	    "ontologyAnnotations.evidence.publications.pubMedId",
	    "ontologyAnnotations.evidence.publications.mgiJnum",
	    "ontologyAnnotations.evidence.publications.mgiId",
	    "ontologyAnnotations.evidence.annotationDate"
	)
	return query

    def annotations(self):
	genoAnnotQuery = self.buildGenoAnnotQuery()

	for genoAnnotRow in genoAnnotQuery.rows():
	    genoAnnot = genoAnnotRow.to_d()

	    # 'PMID:' prefix
	    pubmedID = genoAnnot['Genotype.ontologyAnnotations.evidence.publications.pubMedId']
	    if pubmedID != '' and pubmedID != None:
		pubmedID = 'PMID:' + str(pubmedID)

	    result = {
	    'genoID'  : genoAnnot['Genotype.primaryIdentifier'],
	    'genoName': genoAnnot['Genotype.name'],
	    'doID'  : genoAnnot['Genotype.ontologyAnnotations.ontologyTerm.identifier'],

	    'doTerm':
	      genoAnnot['Genotype.ontologyAnnotations.ontologyTerm.name'],

	    'qualifier':
	      genoAnnot['Genotype.ontologyAnnotations.qualifier'],

	    'evidenceCode':
	      genoAnnot['Genotype.ontologyAnnotations.evidence.code.code'],

	    'pubmedID': pubmedID,

	    'mgiJnum':
	      genoAnnot['Genotype.ontologyAnnotations.evidence.publications.mgiJnum'],
	    'mgiRefID':
	      genoAnnot['Genotype.ontologyAnnotations.evidence.publications.mgiId'],
	    'annotDate':
	      genoAnnot['Genotype.ontologyAnnotations.evidence.annotationDate'],
	    }
	    yield result

####### end class GenoAnnotationQuery

class DiseaseAnnotationFormatter(object):
    # Knows the (json) structure of disease annotation.
    # Creates/returns the json object for a disease annotation.

    def __init__(self, config, service):
	self.geneFinder = GenoToGeneFinder(service)
	self.AGRjf      = AGRjsonFormatter(config)

    def getGeneAnnotJsonObjs(self, ga):
	# return list of gene annot json objects representing
	#   the genoAnnot (ga)

	doID = ga['doID']

	rolledGenes = self.geneFinder.genoToRolledUpGenes( ga['genoID'] )

	geneJsons = []
	for geneID, symbol in rolledGenes:
	    # Fields that are commented out are ones MGI doesn't use.
	    # No reason to set them. If they were are null or [], stripNulls()
	    #   would remove them anyway.
	    geneJsons.append( self.AGRjf.stripNulls( \
	    {
	    'taxonId'			: GLOBALTAXONID,
	    'objectId'			: self.AGRjf.addIDprefix(geneID),
	    'objectName'		: symbol,
	    'objectRelation'		: self.getGeneObjectRelationObj(ga),
	    #'experimentalConditions'	: [],
	    'qualifier'			: self.getQualifier(ga),
	    'DOid'			: doID,
	    #'with'			: [],
	    #'modifier'			: None,
	    'evidence'			: self.getEvidenceObj(ga),
	    #'geneticSex'		: '',
	    'dateAssigned'		: ga['annotDate'] + RFC3339TIME,
	    'dataProvider'		: 'MGI',
	    } )
	    )
	return geneJsons
    ######

    def getGeneObjectRelationObj(self, ga):
	# see https://github.com/alliance-genome/agr_schemas/blob/development/disease/diseaseObjectRelation.json
	
	return \
	  { "associationType"	: "causes_condition",
	    "objectType"	: "gene",
	    "inferredFromID"	: ga["genoID"],		# MGI defined json attr
	    "inferredFromName"	: ga["genoName"]	# MGI defined json attr
	  }

    def getGenoAnnotJsonObjs(self, ga):
	# return list of genotype annot json objects representing
	#   the genoAnnot (ga)

	doID = ga['doID']

	# Fields that are commented out are ones MGI doesn't use.
	# No reason to set them. If they were set to null or [], stripNulls()
	#   would remove them anyway.
	return [ \
	    self.AGRjf.stripNulls( \
	    {
	    'taxonId'		: GLOBALTAXONID,
	    'objectId'		: ga['genoID'],
	    'objectName'	: ga['genoName'],
	    'objectRelation'	: self.getGenoObjectRelationObj(ga),
	    #'experimentalConditions': [],
	    'qualifier'		: self.getQualifier(ga),
	    'DOid'		: doID,
	    #'with'		: [],
	    #'modifier'		: None,
	    'evidence'		: self.getEvidenceObj(ga),
	    #'geneticSex'	: '',
	    'dateAssigned'	: ga['annotDate'] + RFC3339TIME,
	    'dataProvider'	: 'MGI',
	    } )
	    ]
	######

    def getGenoObjectRelationObj(self, ga):
	# see https://github.com/alliance-genome/agr_schemas/blob/development/disease/diseaseObjectRelation.json
	
	rolledGenes = [ self.AGRjf.addIDprefix(g[0]) for g in
		      self.geneFinder.genoToRolledUpGenes(ga['genoID']) ]

	return \
	  { "associationType" : "is_model_of",
	    "objectType"      : "genotype",
	    "inferredGeneAssociation" : rolledGenes
	  }

    def getQualifier(self,ga):
	if ga['qualifier'] is None: return None
	else: return ga['qualifier'].lower()

    def getEvidenceObj(self, ga):
	# see https://github.com/alliance-genome/agr_schemas/blob/development/disease/evidence.json
	# for us, because of the way we pull annotations from MM,
	#    we only have one pub for each annotation record.

	pmID = ga['pubmedID']
	if pmID == '': pmID = None

	return \
	[ { "evidenceCode" : ga['evidenceCode'],
	    "publications" : [ { "modPublicationId" : ga['mgiRefID'],
	                         "pubMedId"         : pmID
			       }
			     ]
	  }
	]
    
######## end Class DiseaseAnnotationFormatter

def main(ids):
    service = Service(MOUSEMINEURL)
    query   = GenoAnnotationQuery(service, ids)
    df      = DiseaseAnnotationFormatter(cp, service)

    annotJsonObjs = []
    for annot in query.annotations():
	if DO_GENOS: annotJsonObjs += df.getGenoAnnotJsonObjs(annot)
				# Can be [] if OMIMID doesn't map to DOID.
				# Once MM starts storing DO annotations, 
				#  this check won't be necessary as the query
				#  will only return annots to DO.
	if DO_GENES: annotJsonObjs += df.getGeneAnnotJsonObjs(annot)
				# Can be [] if OMIMID doesn't map to DOID.
				# OR geno annot doesn't roll up to genes

    jobj = {
      "metaData" : buildMetaObject(service),
      "data"     : annotJsonObjs
      }
    print json.dumps(jobj, sort_keys=True, indent=2, separators=(',', ': ')),

#####

main(sys.argv[1:]) 	# optional geno IDs on the cmd line to restrict query
