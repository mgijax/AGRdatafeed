# helpers for AGR data feeds

import types
import time
import json
# See: http://henry.precheur.org/projects/rfc3339 
from rfc3339 import rfc3339

#-----------------------------------
# The AGR spec is to leave out attributes that would otherwise have a null
#   value or an empty list value. This function removes these from a Json obj.
#
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

#-----------------------------------
# RFC 3339 timestamps

# Returns the current date-time in RFC-3339 format
# Example: "2017-01-26T15:00:42-05:00"
#
def getTimeStamp():
    return rfc3339(time.time())

#-----------------------------------
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

if __name__ == "__main__":
    obj =  \
    {   "keep1" : 23,
	"lose1" : None,
	"lose2" : [],
	"keep2" : [23, 64, 'hike', [],
		    { 'keep3': 78, 'lose3': [], 'lose4':None }]
    }
    print "Before"
    print json.dumps(obj, sort_keys=True, indent=2, separators=(',', ': ') )
    print "After stripNulls"
    print json.dumps(stripNulls(obj), sort_keys=True, indent=2, separators=(',', ': ') )


