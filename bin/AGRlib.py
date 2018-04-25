# helpers for AGR data feeds

import types
import time
import json
import re
import datetime
import os


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
# Constructs and returns the metaData (header) for the dump file.
#
def buildMetaObject(service):
    # get current date
    currentDate = getTimeStamp()

    # Get the MGI release date, which is embedded in the description field
    #   of the MGI DataSource obj
    # For example, "Mouse Genome Informatics [MGI 6.07 2017-01-24]"
    release = None
    query = service.new_query("DataSource")
    query.add_view("name", "description")
    query.add_constraint("name", "=", "MGI", code = "B")
    for r in query:
      i = r.description.find('[')
      release = r.description[i:].strip()[1:-1].strip()

    return {
    "dataProvider" : [buildMgiDataProviderObject()],
    "dateProduced" : currentDate,
    "release" : release
    }
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
