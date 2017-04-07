# helpers for AGR data feeds

import types
import time
import json
import re
from ConfigParser import ConfigParser

# See: http://henry.precheur.org/projects/rfc3339 
from rfc3339 import rfc3339

IDPREFIX = "MGI:"

#-----------------------------------
class AGRjsonFormatter(object):
    # some helpers that know certain details of the AGR json schemas

    def __init__(self, config):
	self.prefixIDs = config.getboolean("DEFAULT", "PREFIX_IDS")
				# MGI IDs, current and old
	self.IDregex = re.compile('^(MGI:(\d)+)|(MGD-[A-Z]+-(\d)+)$')

    # return a clean JSON obj
    def cleanse(self,obj):
	self.stripNulls(obj)
	self.fixIDs(obj)
	return obj

    # The AGR spec is to leave out attributes that would otherwise have a null
    #   value or an empty list value.
    #  This function removes these from a Json obj.
    def stripNulls(self, obj):
	for k,v in obj.items():
	    if v is None or type(v) is types.ListType and len(v) == 0:
		del obj[k]
	    elif type(v) is types.ListType:
		for i in v:
		    if type(i) is types.DictType: self.stripNulls(i)
	    elif type(v) is types.DictType:
		self.stripNulls(v)
	return obj

    # Add "MGI:" prefix to all MGI IDs in json obj (recursively).
    # Not clear we want to do this. AGR may only want to add the prefix for
    #  primary gene IDs. Not all. Still in flux.
    def fixIDs(self, obj):

	if not self.prefixIDs: return

	if type(obj) is types.ListType:
	    for i,v in enumerate(obj):
		if type(v) is types.StringType or type(v) is types.UnicodeType:
		    obj[i] = self.fixID(v)
		else: self.fixIDs(v)

	if type(obj) is types.DictType:
	    for i, v in obj.items():
		if type(v) is types.StringType or type(v) is types.UnicodeType:
		    obj[i] = self.fixID(v)
		else: self.fixIDs(v)
	return obj
	
    # add prefix only if it matches the ID template
    def fixID(self, s):
	if self.prefixIDs and re.match(self.IDregex,s):
	    s = IDPREFIX + s
	return s

    # assuming s is an MGI ID that we want to add the prefix, do so
    def addIDprefix(self,s):
	if self.prefixIDs: return IDPREFIX + s
	else:  return s

# ----- end class AGRjsonFormatter

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
    "dataProvider" : "MGI",
    "dateProduced" : currentDate,
    "release" : release
    }
#-----------------------------------
# RFC 3339 timestamps

# Returns the current date-time in RFC-3339 format
# Example: "2017-01-26T15:00:42-05:00"
#
def getTimeStamp():
    return rfc3339(time.time())

#-----------------------------------

if __name__ == "__main__":
    Ajf = AGRjsonFormatter( ConfigParser({ "PREFIX_IDS": "true"}) )
    print Ajf.fixID("MGI:12345")
    print Ajf.fixID("MGD-MRK-12345")
    print Ajf.fixID("MGD-MRK-MRK-12345")
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
    Ajf.stripNulls(obj)
    print json.dumps(obj, sort_keys=True, indent=2, separators=(',', ': ') )
