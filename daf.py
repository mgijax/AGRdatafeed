#!/usr/bin/env python2.7

# Create DAF (disease annot file) file for AGR data ingest.
# Pulls data from MouseMine.
# Writes DAF to stdout.
# For debugging, accepts optional MGI IDs for genotypes on command line
#   to restrict output to DAF records for those genotypes.
# Example genotype ID:  MGI:5526095 MGI:2175208
#
# Author: jim kadin

import sys
from ConfigParser import ConfigParser
from intermine.webservice import Service

##### Load config
cp = ConfigParser()
cp.optionxform = str
cp.read("config.cfg")

MOUSEMINEURL    = cp.get("DEFAULT","MOUSEMINEURL")
TAXONID         = cp.get("DEFAULT","TAXONID")
DAFHEADER	= cp.get("dafFile","DAFHEADER") 

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
				# in theory, should only be 1 DO ID, right?
	for row in query.rows():
	    OmimID = str(row["identifier"])
	    DOID   = row["crossReferences.identifier"]
	    self.OmimToDO.setdefault(OmimID,set()).add(DOID)
####### end class OmimToDOfinder

class GenoToGeneFinder (object):
    # find rolled up gene(s) for a genotype to disease annotation
    # build and provide a mapping from OMIM IDs to DO IDs

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

#####
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

class DafFormatter(object):
    # knows the structure of a DAF file, formats a line in a DAF file, etc.

    #########
    # Structure of a DAF file.
    # Some of these fields (columns) are constant for us, some come from 
    #  from MouseMine, some are computed.
    # The "field" for a column is the genotype Annotation record field name
    #	that holds the contents for the column.
    #
    dafColumns = [		# ordered by the columns in the DAF file
    {'colName': 'Taxon',		       'constant': 'taxon:'+ TAXONID},
    {'colName': 'DB Object Type',	       'constant': 'genotype'},
    {'colName': 'DB',		 	       'constant': 'MGI'},
    {'colName': 'DB Object ID',  	          'field': 'genoID'},
    {'colName': 'DB Object Symbol',	          'field': 'genoName'},
    {'colName': 'Inferred gene association',      'field': 'ROLLEDUPGENES'},
    {'colName': 'Gene Product Form ID',	       'constant': ''},
    {'colName': 'Experimental conditions',     'constant': ''},
    {'colName': 'Association type',            'constant': 'is_model_of'},
    {'colName': 'Qualifier',			  'field': 'qualifier'},
    {'colName': 'DO ID',                          'field': 'DOIDs'},
    {'colName': 'With',                        'constant': ''},
    {'colName': 'Modifier - association type', 'constant': ''},
    {'colName': 'Modifier - Qualifier',	       'constant': ''},
    {'colName': 'Modifier - genetic',	       'constant': ''},
    {'colName': 'Modifier - experimental conditions',
                                               'constant': ''},
    {'colName': 'Evidence Code',		  'field': 'evidenceCode'},
    {'colName': 'genetic sex',                 'constant': ''},
    {'colName': 'DB:Reference', 		  'field': 'REFID' },
    {'colName': 'Date',                           'field': 'annotDate'},
    {'colName': 'Assigned By',                 'constant': 'MGI'},
    ]
    #########

    def __init__(self, service):
	self.geneFinder = GenoToGeneFinder(service)
	self.doFinder   = OmimToDOfinder(service)

    def getDafHeader(self):
	return DAFHEADER + '\n' + \
	    '\t'.join([col['colName'] for col in DafFormatter.dafColumns]) +'\n'

    def omim2DO(self, omimID):
	return '|'.join( self.doFinder.omimToDO(omimID) )

    def getRolledUpGenes(self, genoID):
	return '|'.join( self.geneFinder.genoToRolledUpGenes(genoID) )
    
    def getReferenceID(self, pmID, mgiRefID):
	return pmID if pmID!='' and pmID!=None else  mgiRefID

    def formatDafRow(self, genoAnnot):
	# return as string representing the DAF row for genoAnnot

	DOIDs = self.omim2DO( genoAnnot['omimID'] )
	if len(DOIDs) == 0:		# skip if no DO ID
	    return ''

	# Compute a few DAF columns
	genoAnnot['DOIDs'] = DOIDs

	genoAnnot['ROLLEDUPGENES'] = self.getRolledUpGenes(genoAnnot['genoID'])

	genoAnnot['REFID'] = self.getReferenceID(genoAnnot['pubmedID'],
							genoAnnot['mgiRefID'])

	dafRow = []
	for dafColDesc in DafFormatter.dafColumns:
	    value = dafColDesc['constant'] if dafColDesc.has_key('constant') \
			    else str(genoAnnot[ dafColDesc['field'] ])
	    if value == 'None': value = ''
	    dafRow.append(value)

	return '\t'.join(dafRow) + '\n'

######## end Class DafFormatter

def main(ids):
    service    = Service(MOUSEMINEURL)
    annotQuery = GenoAnnotationQuery(service, ids)
    df         = DafFormatter(service)

    sys.stdout.write( df.getDafHeader() )

    for genoAnnot in annotQuery.annotations():
	sys.stdout.write( df.formatDafRow( genoAnnot ) )
#####

main(sys.argv[1:]) 	# optional geno IDs on the cmd line to restrict query

