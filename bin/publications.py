#
# publications.py
#
# Dumps all publications for loading into the ABC.
#

import sys
import os
import json
import argparse
from AGRlib import stripNulls, buildMetaObject, makeOneOfConstraint, sql, getTimeStamp

MOUSEMINE     = os.environ["MOUSEMINEURL"]

AGR_REF_CATS = [
    "Research Article",
    "Review Article",
    "Thesis",
    "Book",
    "Other",
    "Preprint",
    "Conference Publication",
    "Personal Communication",
    "Direct Data Submission",
    "Internal Process Reference",
    "Unknown",
    "Retraction"
    ]
MGI2AGR_REF_TYPES = {
    "MGI Direct Data Submission" : "Direct Data Submission",
    "Not Specified" : "Unknown",
    "MGI Curation Record" : "Internal Process Reference",
    "Unreviewed Article" : "Preprint",
    "MGI Data Load" : "Internal Process Reference",
    "Book" : "Book",
    "Newsletter" : "Other",
    "External Resource" : "Other",
    "Annual Report/Bulletin" : "Other",
    "Peer Reviewed Article" : "Research Article",
    "Dissertation/Thesis" : "Thesis",
    "Personal Communication" : "Personal Communication",
    "Conference Proceedings/Abstracts" : "Conference Publication",
    "JAX Notes" : "Other",
    }

REFS_Q = '''
    SELECT 
        b.title,
        b.authors,
        b.journal,
        b.vol,
        b.issue,
        b.pgs,
        b.date,
        b.modification_date,
        bc.citation,
        bc.mgiid,
        bc.pubmedid,
        b.abstract,
        bc.relevanceterm,
        bc.referencetype,
        bb.book_au,
        bb.book_title,
        bb.place,
        bb.publisher,
        bb.series_ed
    FROM
        BIB_Refs b
        LEFT OUTER JOIN BIB_Books bb ON b._refs_key = bb._refs_key
        JOIN BIB_Citation_cache bc ON b._refs_key = bc._refs_key

    /* optional WHERE clause */
    %s 

    '''

def getExchangeObj (r) :
    pmid = 'PMID:%s' % r['pubmedid']
# Args:
#    r - A reference record from the database
#    mode - One of all, pubmed, nonpubmed
# Returns:
#    object or None
#
def getObj (r, which) :
    pmid = 'PMID:%s' % r['pubmedid']
    primaryId = pmid if r['pubmedid'] else r['mgiid']
    #
    if which == "pubmed" and not r['pubmedid']:
        return None

    if which == "nonpubmed" and r['pubmedid']:
        return None

    #
    if which == "all" or which == "nonpubmed":
        #
        return {
            'primaryId'        : primaryId,
            'title'            : r['title'] if r['title'] else '',
            'authors'          : getAuthors(r, primaryId),
            'datePublished'    : r['date'],
            'dateLastModified' : getTimeStamp(r['modification_date']),
            'volume'           : r['vol'] if r['vol'] else '',
            'pages'            : r['pgs'] if r['pgs'] else '',
            'abstract'         : r['abstract'] if r['abstract'] else '',
            'citation'         : r['citation'] if r['citation'] else '',
            'issueName'        : r['issue'] if r['issue'] else '',
            'allianceCategory' : getAllianceCategory(r),
            'resourceAbbreviation' : r['journal'] or r['book_title'] or '',
            'MODReferenceTypes' : [{ 'referenceType' : r['referencetype'], 'source' : 'MGI' }],
            'tags'              : getTags(r, primaryId),
            'crossReferences'   : [{'id': r['mgiid'],'pages':['reference']}],
        }
    else:
        #
        return {
            'modId'            : r['mgiid'],
            'pubMedId'         : pmid,
            'allianceCategory' : getAllianceCategory(r),
            'dateLastModified' : getTimeStamp(r['modification_date']),
            'MODReferenceTypes' : [{ 'referenceType' : r['referencetype'], 'source' : 'MGI' }],
            'tags'              : getTags(r, primaryId),
        }

def getAllianceCategory (r) :
    return MGI2AGR_REF_TYPES[r['referencetype']]

def parseAuthor (s) :
    return {
        'name' : s
    }

def getAuthors (r, pid) :
    if r['authors'] is None or r['authors'] == "":
        return []
    authors = list(map(parseAuthor, r['authors'].split(';')))
    for a in authors:
        aid = (a['name'] + pid).replace(" ", "")
        a['referenceId'] = aid
    return authors

def getTags (r, pid) :
    tags = []
    if r['relevanceterm'] == "discard":
        tags.append('notRelevant')
    else:
        tags.append('inCorpus')

    return list(map(lambda t: {"referenceId":pid, "tagName":t, "tagSource":"MGI"}, tags))

def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        dest="which",
        choices=["all", "pubmed", "nonpubmed"],
        default="all",
        help="Which publications to output. all=all pubs in long form; pubmed=pubmed pubs in exchange form; nonpubmed=Non-pubmed pubs in long form."
        )
    parser.add_argument(
        "-s",
        dest="doSample",
        action="store_true",
        default=False,
        help="Just do a small sample for dev purpposes."
        )
    parser.add_argument(
        "-i",
        dest="mgiId",
        default=None,
        help="Just run for this MGI publication id (for debugging)."
        )
    return parser.parse_args()

def main () :
    opts = getArgs()
    whereClause = ''
    if opts.mgiId:
        whereClause = "WHERE bc.mgiid = '%s'" % opts.mgiId
    #
    print('{\n  "metaData": %s,\n  "data": [' % json.dumps(buildMetaObject(MOUSEMINE), indent=2))
    #
    n = 0
    for r in sql(REFS_Q % whereClause):
        obj = getObj(r, opts.which)
        if obj:
            if n > 0: sys.stdout.write(",")
            print(json.dumps(obj, indent=2))
            n += 1
            if opts.doSample and n > 100:
                break
    #
    print(']}')

main()
