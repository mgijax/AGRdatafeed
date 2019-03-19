# helpers for AGR data feeds

import types
import time
import json
import re
import datetime
import os
import urllib
import itertools
import sys


#
def getConfig():
    from ConfigParser import ConfigParser
    DIR=os.path.dirname(__file__)
    cfname= os.path.abspath(os.path.join(DIR,"..","config.cfg"))
    cp = ConfigParser()
    cp.optionxform = str # make keys case sensitive
    cp.read(cfname)
    return cp


# See: http://henry.precheur.org/projects/rfc3339 
from rfc3339 import rfc3339
date_re = re.compile(r'(\d\d\d\d)-(\d\d)-(\d\d)')

# The AGR spec is to leave out attributes that would otherwise have a null
#   value or an empty list value.
#  This function removes these from a Json obj.
def stripNulls(obj):
    for k,v in obj.items():
        if v is None or type(v) is types.ListType and len(v) == 0:
            del obj[k]
        elif type(v) is types.ListType:
            for i in v:
                if type(i) is types.DictType: stripNulls(i)
        elif type(v) is types.DictType:
            stripNulls(v)
    return obj

#
def buildMgiDataProviderObject () :
    return {
        "crossReference" : {
	    "id" : "MGI",
	    "pages" : ["homepage"]
	},
	"type" : "curated"
    }
#
def doInterMineQuery(q, url): 
    # sys.stderr.write('Intermine query=' + q + '\n')
    fmt = 'tab'
    url = '%s/query/results?format=%s&query=%s' % (url,fmt,urllib.quote_plus(q))
    fd = urllib.urlopen(url)
    for line in fd: 
        toks = line[:-1].split('\t')
        yield toks
    fd.close()

# Returns the list of view paths from the given query.
VIEW_re=re.compile(r'view="([^"]*)"')
def getView (q, stripRoot=True):
    view = None
    m = VIEW_re.search(q)
    if m:
      view = m.group(1).strip().split()
      if stripRoot:
        view = map(lambda v: '.'.join(v.split('.')[1:]), view)
    return view

# 
def doQuery(q, url):
  view = getView(q)
  def makeObject (row) :
    return dict(zip(view, map(lambda f: f if f != '""' else None, row)))
  return itertools.imap(makeObject, doInterMineQuery(q, url))

# Constructs and returns the metaData (header) for the dump file.
#
def buildMetaObject(mouseMineUrl):
    # get current date
    currentDate = getTimeStamp()

    # Get the MGI release date, which is embedded in the description field
    #   of the MGI DataSource obj
    # For example, "Mouse Genome Informatics [MGI 6.07 2017-01-24]"
    release = None
    q = '''<query 
      model="genomic"
      view="DataSource.name DataSource.description"
      >
      <constraint path="DataSource.name" op="=" value="MGI"/>
      </query>
      '''
    for r in doInterMineQuery(q,mouseMineUrl):
      i = r[1].find('[')
      release = r[1][i:].strip()[1:-1].strip()

    return {
    "dataProvider" : buildMgiDataProviderObject(),
    "dateProduced" : currentDate,
    "release" : release
    }
#
def makeOneOfConstraint(path, vals):
  cnst = ''
  if vals:
    cvals = ''.join(map(lambda i: '<value>%s</value>' % i, vals))
    cnst = '<constraint path="%s" op="ONE OF">%s</constraint>' % (path,cvals)
  return cnst

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
    return {
      "publicationId" : pid,
      "crossReference" : {"id":mgiId ,"pages":["reference"]}
    }
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

#-----------------------------------

if __name__ == "__main__":
    obj =  \
    {   "keep1" : 23,
	"lose1" : None,
	"lose2" : [],
	"keep2" : [23, 64, 'hike', [],
		    { 'keep3': 78, 'lose3': [], 'lose4':None }],
	"prefix MGI:" : "MGI:12345",
	"no prefix"   : "MGI:foo",
	"prefix MGI: 2" : [ { "foo": "MGI:23456"} ]
    }
    print "Before"
    print json.dumps(obj, sort_keys=True, indent=2, separators=(',', ': ') )
    print "After cleansing"
    stripNulls(obj)
    print json.dumps(obj, sort_keys=True, indent=2, separators=(',', ': ') )
