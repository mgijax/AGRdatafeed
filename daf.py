#!/usr/bin/env python

# Create DAF (disease annot file) file for AGR data ingest.
# Pulls data from MouseMine.
# Writes DAF to stdout.
# For debugging, accepts optional MGI IDs for genotypes on command line
#   to restrict output to DAF records for those genotypes.
#
# Author: jim kadin

import sys
#import json

from intermine.webservice import Service

#########
# Constants for genotype Annotation record field names that are used in the
#   multiple places
GENO_ID  = 'Genotype.primaryIdentifier'

# Structure of a DAF file.
# Some of these fields are constant for us, some come from data pulled
#  from MouseMine.
# The "field" for a column is the genotype Annotation record field name
#  (most are determined by the fields returned from the MM query)
dafColumns = [		# ordered by the columns in the DAF file
    {'colName': 'Taxon',		       'constant': 'taxon:10090'},
    {'colName': 'DB Object Type',	       'constant': 'genotype'},
    {'colName': 'DB',		 	       'constant': 'MGI'},
    {'colName': 'DB Object ID',  	          'field': GENO_ID},
    {'colName': 'DB Object Symbol',	          'field': 'Genotype.name'},
    {'colName': 'Inferred gene association',      'field': 'RolledUpGenes'},
		#'field': is computed, not from the query output
    {'colName': 'Gene Product Form ID',	       'constant': ''},
    {'colName': 'Experimental conditions',     'constant': ''},
    {'colName': 'Association type',            'constant': 'is_model_of'},
    {'colName': 'Qualifier',
		'field': 'Genotype.ontologyAnnotations.qualifier'},
    {'colName': 'DO ID',                          'field': 'DOIDs'},
		#'field': is computed, not from the query output
    {'colName': 'With',                        'constant': ''},
    {'colName': 'Modifier - association type', 'constant': ''},
    {'colName': 'Modifier - Qualifier',	       'constant': ''},
    {'colName': 'Modifier - genetic',	       'constant': ''},
    {'colName': 'Modifier - experimental conditions',
                                               'constant': ''},
    {'colName': 'Evidence Code',
		'field': 'Genotype.ontologyAnnotations.evidence.code.code'},
    {'colName': 'genetic sex',                 'constant': ''},
    {'colName': 'DB:Reference', 		  'field': 'REF_ID' },
		#'field': is computed, not from the query output
    {'colName': 'Date',                        'constant': '2016/12/25'},
    {'colName': 'Assigned By',                 'constant': 'MGI'},
]
#####

def formatDafRow( annotRow):
    # return an annotation row formatted as a DAF row

    dafRow = []
    for dafColDesc in dafColumns:
	value = dafColDesc['constant'] if dafColDesc.has_key('constant') \
			else str(annotRow[ dafColDesc['field'] ]) 
	if value == 'None': value = ''
	dafRow.append(value)

    return '\t'.join(dafRow)
#####
class OmimToDOfinder (object):
    # find DO term for a given OMIM term ID

    def __init__(self, service):
	query = self.buildQuery(service)
	self.buildMapping(query)
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
				# in theory, should only be 1 DO ID, right?
	for row in query.rows():
	    OmimID = row["identifier"]
	    DOID = row["crossReferences.identifier"]
	    self.OmimToDO.setdefault(OmimID,set()).add(DOID)
    ######

    def omimToDO(self, OmimID):
	return self.OmimToDO.get(OmimID, set())

####### end class OmimToDOfinder

class GenoToGeneFinder (object):
    # find rolled up gene(s) for a genotype to disease annotation

    ### HMMM, seems like we are only using genotype ID to map to rolled up
    # genes. Is it possible to have a genotype, say with 2 disease annot's
    #        and both of those do not roll up to both those genes?
    # I.e, is the roll up based only on the genotype, or can it be based on
    #  the annotation itself?

    def __init__(self, service):
	query = self.buildQuery(service)
	self.buildMapping(query)
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
				# in theory, should only be 1 causitive gene
				#  per genotype annot, but allow a set for now.
	for row in query.rows():
	    genoID = row["ontologyAnnotations.evidence.baseAnnotations.subject.primaryIdentifier"]
	    geneID = row["primaryIdentifier"]
	    self.genoToGene.setdefault(genoID,set()).add(geneID)
    ######

    def genoToRolledUpGenes(self, genoID):
	return self.genoToGene.get(genoID, set())

####### end class GenoToGeneFinder

#####
def buildGenoAnnotQuery(service, ids):
    # ids is list of Genotype MGI ids to restrict the query to for testing.
    #   if empty list, get them all

    # Get a new query on the class (table) you will be querying:
    query = service.new_query("Genotype")

    #Type constraints come early before all mentions of the paths they constrain
    query.add_constraint("ontologyAnnotations.ontologyTerm", "OMIMTerm")
    if len(ids):
	query.add_constraint("primaryIdentifier", "ONE OF", ids)

    # The view specifies the output columns
    query.add_view(
	"primaryIdentifier", "name", "ontologyAnnotations.ontologyTerm.omimId",
	"ontologyAnnotations.ontologyTerm.name", "ontologyAnnotations.qualifier",
	"ontologyAnnotations.evidence.code.code",
	"ontologyAnnotations.evidence.publications.pubMedId",
	"ontologyAnnotations.evidence.publications.mgiJnum",
	"ontologyAnnotations.evidence.publications.mgiId"
    )
    return query
#####

def main(ids):
    service = Service("http://www.mousemine.org/mousemine/service")
    geneFinder = GenoToGeneFinder(service)
    doFinder = OmimToDOfinder(service)

    # print column headings
    print '\t'.join( [col['colName'] for col in dafColumns] )

    genoAnnotQuery = buildGenoAnnotQuery(service, ids)

    for genoAnnotRow in genoAnnotQuery.rows():
	#print '--------'
	genoAnnot = genoAnnotRow.to_d()

	### Clean up the genoAnnot
	# shorter field names
	OMIM_ID  = "Genotype.ontologyAnnotations.ontologyTerm.omimId"
	PM_ID    = 'Genotype.ontologyAnnotations.evidence.publications.pubMedId'
	#MGI_Jnum = 'Genotype.ontologyAnnotations.evidence.publications.mgiJnum'
	MGI_Ref_ID = 'Genotype.ontologyAnnotations.evidence.publications.mgiId'

	# rolled up genes
	genoID = genoAnnot[GENO_ID]
	genoAnnot['RolledUpGenes'] = \
			    '|'.join(geneFinder.genoToRolledUpGenes(genoID))

	# OMIM --> DO
	omimID = "OMIM:" + str(genoAnnot[OMIM_ID])
	DOIDs = doFinder.omimToDO(omimID)
	if len(DOIDs) == 0:		# skip if no DO ID
	    #print "skipping %s, OMIM: %s" % (genoID, omimID)
	    continue
	genoAnnot['DOIDs'] = '|'.join(DOIDs)

	# J# if no PMID
	pmID = str(genoAnnot[PM_ID])
	refID = genoAnnot[MGI_Ref_ID]
	genoAnnot['REF_ID'] = "PMID:" + pmID if pmID!='' and pmID!='None' \
				else  refID

	print formatDafRow(genoAnnot)
#####

main(sys.argv[1:]) 	# optional geno IDs on the cmd line to restrict query
