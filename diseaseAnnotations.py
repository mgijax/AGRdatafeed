#!/usr/bin/env python2.7

# Create JSON file describing genotype-disease annotations for AGR data ingest.
# Pulls data from MouseMine.
# Writes JSON to stdout.
# See JSON spec at:
#https://github.com/alliance-genome/agr_schemas/blob/development/disease/diseaseModelAnnotation.json
# (this is development branch, might want to see master branch too)
# For debugging, accepts optional MGI IDs for genotypes on command line
#   to restrict output to DAF records for those genotypes.
# Example genotype ID:  MGI:5526095 MGI:2175208 MGI:5588576
#
# Author: jim kadin

import sys
import json
from ConfigParser import ConfigParser
from intermine.webservice import Service
from AGRlib import stripNulls, buildMetaObject

##### Load config
cp = ConfigParser()
cp.optionxform = str
cp.read("config.cfg")

MOUSEMINEURL    = cp.get("DEFAULT","MOUSEMINEURL")
TAXONID         = cp.get("DEFAULT","TAXONID")

########
class OmimToDOfinder (object):
    # build and provide a mapping from OMIM IDs to DO IDs

    def __init__(self, service):
	query = self.buildQuery(service)
	self.buildMapping(query)
    ######
    def omimToDO(self, OmimID):
	return self.OmimToDO.get(OmimID, set())

    ######
    def buildQuery(self,service):

	# Get a new query on the class (table) you will be querying:
	query = service.new_query("OMIMTerm")

	# Type constraints should come early
	query.add_constraint("crossReferences", "DOTerm")

	# The view specifies the output columns
	query.add_view("crossReferences.identifier",
			"crossReferences.name",
			"identifier",		# from MM, these have "OMIM:"
			"name")
	return query
    ######

    def buildMapping(self, query):

	self.OmimToDO = {}	# OmimToDO[OMIM ID] == set(DO IDs)
				# could be 0 or >1 in odd cases
	for row in query.rows():
	    OmimID = str(row["identifier"])
	    DOID   = row["crossReferences.identifier"]
	    self.OmimToDO.setdefault(OmimID,set()).add(DOID)
####### end class OmimToDOfinder

class GenoToGeneFinder (object):
    # find rolled up gene(s) for a genotype to disease annotation

    def __init__(self, service):
	query = self.buildQuery(service)
	self.buildMapping(query)
    ######

    def genoToRolledUpGenes(self, genoID):
	return self.genoToGene.get(genoID, set())

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

    def buildMapping(self, query):

	self.genoToGene = {}	# geneoToGene[genoID] == set(GeneIDs)
				# There can be multiple rolled up genes:
				#  1 example: gene + transgene if the TG
				#   expresses the gene.
	for row in query.rows():
	    genoID = row["ontologyAnnotations.evidence.baseAnnotations.subject.primaryIdentifier"]
	    geneID = row["primaryIdentifier"]
	    self.genoToGene.setdefault(genoID,set()).add(geneID)
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
	query.add_constraint("ontologyAnnotations.ontologyTerm", "OMIMTerm")
	if len(self.ids):
	    query.add_constraint("primaryIdentifier", "ONE OF", self.ids)

	# The view specifies the output columns
	query.add_view(
	    "primaryIdentifier", "name",
	    "ontologyAnnotations.ontologyTerm.omimId",
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

	    # 'OMIM:' prefix
	    omimID = 'OMIM:' + \
	     str(genoAnnot['Genotype.ontologyAnnotations.ontologyTerm.omimId'])

	    result = {
	    'genoID'  : genoAnnot['Genotype.primaryIdentifier'],
	    'genoName': genoAnnot['Genotype.name'],
	    'omimID'  : omimID,

	    'omimTerm':
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

    def __init__(self, service):
	self.geneFinder = GenoToGeneFinder(service)
	self.doFinder   = OmimToDOfinder(service)

    def getAnnotJsonObj(self, ga):
	# return a json object representing a genoAnnot (ga)

	# get list of matching DO IDs based on annotation's OMIM ID
	DOIDs = [ d for d in self.doFinder.omimToDO(ga['omimID']) ]

	if len(DOIDs) != 1:
	    return None		# skip if no DO ids or more than one

	doID = DOIDs.pop()

	# Fields that are commented out are ones MGI doesn't use.
	# No reason to set them. If they were set to null or [], stripNulls()
	#   would remove them anyway.
	return stripNulls( \
	{
	'taxonID'			: TAXONID,
	'objectId'			: ga['genoID'],
	'objectName'			: ga['genoName'],
	'objectRelation'		: self.getObjectRelationObj(ga),
	#'experimentalConditions'	: [],
	'qualifier'			: ga['qualifier'],
	'DOid'				: doID,
	#'with'				: [],
	#'modifier'			: None,
	'evidence'			: self.getEvidenceObj(ga),
	#'geneticSex'			: '',
	'dateAssigned'			: ga['annotDate'],
	'dataProvider'			: 'MGI',
	} )
    ######

    def getEvidenceObj(self, ga):
	# see https://github.com/alliance-genome/agr_schemas/blob/development/disease/evidence.json
	# for us, because of the way we pull annotations from MM,
	#    we only have one pub for each annotation record.

	pmID = ga['pubmedID']
	if pmID == '': pmID = None

	return \
	[ { "evidenceCode" : ga['evidenceCode'],
	    "publications" : [ { "publicationModId" : ga['mgiRefID'],
	                         "pubMedId"         : pmID
			       }
			     ]
	  }
	]

    def getObjectRelationObj(self, ga):
	# see https://github.com/alliance-genome/agr_schemas/blob/development/disease/diseaseObjectRelation.json
	
	return \
	{ "objectRelation" :
	  { "associationType" : "is_model_of",
	    "objectType"      : "genotype",
	    "inferredGeneAssociation" :
	      [ g for g in self.geneFinder.genoToRolledUpGenes( ga['genoID'] ) ]
	  }
	}
    
######## end Class DiseaseAnnotationFormatter

def main(ids):
    service = Service(MOUSEMINEURL)
    query   = GenoAnnotationQuery(service, ids)
    df      = DiseaseAnnotationFormatter(service)

    annotJobjs = []
    for annot in query.annotations():
	annotJson = df.getAnnotJsonObj(annot)
	if annotJson != None:	# Can be None if the OMIMID doesn't map to DOID.
				# Once MM starts storing DO annotations, 
				#  this check won't be necessary as the query
				#  will only return annots to DO.
	    annotJobjs.append(annotJson)

    jobj = {
      "metaData" : buildMetaObject(service),
      "data"     : annotJobjs
      }
    print json.dumps(jobj, sort_keys=True, indent=2, separators=(',', ': ')),

#####

main(sys.argv[1:]) 	# optional geno IDs on the cmd line to restrict query
