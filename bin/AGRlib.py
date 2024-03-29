# helpers for AGR data feeds

import types
import time
import json
import re
import datetime
import os
import urllib.request, urllib.parse, urllib.error
import itertools
import sys
import subprocess
import db

#----------------------------------
# See: http://henry.precheur.org/projects/rfc3339 
from rfc3339 import rfc3339
date_re = re.compile(r'(\d\d\d\d)-(\d\d)-(\d\d)')

#----------------------------------
# The AGR spec is to leave out attributes that would otherwise have a null
#   value or an empty list value.
#  This function removes these from a Json obj.
def stripNulls(obj):
    for k,v in list(obj.items()):
        if v is None or type(v) is list and len(v) == 0:
            del obj[k]
        elif type(v) is list:
            for i in v:
                if type(i) is dict: stripNulls(i)
        elif type(v) is dict:
            stripNulls(v)
    return obj

#----------------------------------
#
def buildMgiDataProviderObject () :
    return {
        "crossReference" : {
            "id" : "MGI",
            "pages" : ["homepage"]
        },
        "type" : "curated"
    }
#----------------------------------
# Constructs and returns the metaData (header) for the dump file.
#
def buildMetaObject():
    # get current date
    currentDate = getTimeStamp()

    # Get the MGI release date, which is embedded in the description field
    #   of the MGI DataSource obj
    # For example, "Mouse Genome Informatics [MGI 6.07 2017-01-24]"
    release = None
    r = list(sql('\nselect * from mgi_dbinfo\n'))[0]
    release = r['public_version'] + " " + r['lastdump_date'].split()[0]
    return {
    "dataProvider" : buildMgiDataProviderObject(),
    "dateProduced" : currentDate,
    "release" : release
    }

#----------------------------------
# Returns a PublicationRef which contains a publication id (required) and an optional xref
# to the mod for that publication.
# The publication id should be the pubmed id , if available. Otherwise, use the MOD pub id.
# See: https://github.com/alliance-genome/agr_schemas/blob/release-1.0.0.8/publicationRef.json
def makePubRef (pubmedId, mgiId) :
  if pubmedId:
    pid = pubmedId
    if not pid.startswith("PMID:"): pid = "PMID:" + pid
  else:
    pid = mgiId
  #
  if pid:
    if mgiId:
        return {
          "publicationId" : pid,
          "crossReference" : {"id":mgiId ,"pages":["reference"]}
        }
    else:
        return { "publicationId" : pid }
  else:
    return None

#-----------------------------------
# RFC 3339 timestamps
#
# With no arguments, returns the current date-time in RFC-3339 format
# With a string argument of the form 'yyyy-mm-dd', converts to rfc3339 format and returns.
# Note that simply concatenating a string such as 'T00:00:00-05:00' to the date is not
# sufficient because of DST differences (i.e., for part of the year, the offset if -04:00
# rather that "-05:00"
# Examples:
#   getTimeStamp() --> "2017-01-26T15:00:42-05:00"
#   getTimeStamp("2007-01-15") --> "2007-01-15T00:00:00-05:00"
#   getTimeStamp("2014-05-01") --> "2014-05-01T00:00:00-04:00"
#
def getTimeStamp(s = None):
    if s:
        m = date_re.match(s)
        d = datetime.datetime(int(m.group(1)), int(m.group(2)), int(m.group(3)))
        return rfc3339(d)
    else:
        return rfc3339(time.time())

#---------------------------------
def sql (query) :
    server = db.get_sqlServer()
    database = db.get_sqlDatabase()
    sys.stderr.write(f"\nSQL query ({server}.{database}): {query}\n")
    for r in db.sql(query):
        yield dict(r)

