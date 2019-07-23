#
# lsAgr.py
#
# Lists the contents of the Alliance downloads area, with optional filtering.
# Deals with (ie hides) the paginated nature of the API, so the user just sees everything.
#
# Example usage:
#       python lsAgr.py -x 10090 -s 1.0.0.8
#
import sys
import xml.etree.ElementTree as ET
try:
  # Python 3
  from urllib.request import urlopen
  from urllib.parse import quote
except:
  # Python 2
  from urllib import urlopen
  from urllib import quote
import re
import json
import argparse

AGRURL='http://download.alliancegenome.org/'

def getOptions ():
  parser = argparse.ArgumentParser(description="List the files that have been uploaded to the Alliance.")
  parser.add_argument('-x', '--taxonid',
    dest='taxonid',
    action="append",
    metavar="TAXON or PROVIDER",
    help="Select for this taxon id (eg 10090) or provider (eg MGI). Some schema versions use taxon ids, and some use provider names. To get all data submission from a given MOD, you have to specify both. Repeatable, to specify multiple.")
  parser.add_argument('-d', '--datatype',
    dest='datatype',
    action="append",
    metavar="TYPE",
    help="Select for this data type. Repeatable, to specify multiple types.")
  parser.add_argument('-s', '--schema',
    dest='schema',
    action="append",
    metavar="VERSION",
    help="Select for this schema version. Repeatable.")
  parser.add_argument('-e', '--etag',
    dest='etag',
    action="append",
    metavar="UID",
    help="Select for this ETag (unique identifier). Repeatable.")
  parser.add_argument('-m', '--modified',
    dest='modified',
    metavar="DATE",
    help="Select for this modification date. Specify dates using yyyy-mm-dd format (eg '2019-02-03'). You may specify a single date, or a date range using 'mindate..maxdate'. You may leave off maxdate to specify anything on or after mindate (eg '2017-01-01..'), or you may leave off mindate to specify anything on or before maxdate (eg, '..2018-12-31).")
  parser.add_argument('-f', '--format',
    dest='format',
    metavar="FORMAT",
    choices=['json','tab'],
    help="Output format.")
  args = parser.parse_args()
  if args.modified:
    m = args.modified
    if '..' in m:
        m = m.split('..')
        if len(m) != 2:
          parser.error('Bad date specification: ' + args.modified)
        args.modified = m
    else:
        args.modified = [m,m]
  return args

# regex to match the leading stuff in braces
RE=re.compile(r'[{].*[}]')
def processOneContentSection (c):
    #
    # quick and dirty conversion to dictionary 
    d = dict([(RE.sub('', x.tag), x.text) for x in c])
    # skip anything that's not a file entry
    if not 'ETag' in d:
      return
    # for some reason, the ETag values include the quote characters
    d['ETag'] = d['ETag'].replace('"','')
    # split the file path into parts
    kparts = d['Key'].split("/")
    if len(kparts) != 4:
      return
    (kSchema, kType, kProvider, kFile) = kparts
    # just grab the yyyy-mm-dd part of the modification date
    mdate = d['LastModified'][:10]
    #
    # Apply filters
    #
    if opts.etag and d['ETag'] not in opts.etag:
      return
    #
    if opts.taxonid and kProvider not in opts.taxonid:
      return
    #
    if opts.schema and kSchema not in opts.schema:
      return
    #
    if opts.datatype and kType not in opts.datatype:
      return
    #
    if opts.modified:
      if opts.modified[0] and mdate < opts.modified[0]:
        return
      if opts.modified[1] and mdate > opts.modified[1]:
        return
    #
    if opts.format == 'tab':
      row = [
        d['LastModified'],
        d['Key'],
        d['ETag'],
        d['Size'],
      ]
      print('\t'.join(row))
    else:
      print(json.dumps(d, indent=2))

def main():
    global opts
    opts = getOptions()
    again = True
    ctoken = None
    while again:
        again = False
        url = AGRURL
        if ctoken:
            url += '&continuation-token=' + quote(ctoken)
        fd = urlopen(url)
        s = fd.read()
        fd.close()
        root = ET.fromstring(s)
        for c in root:
          ctag = RE.sub('', c.tag)
          if ctag == 'Contents':
              processOneContentSection(c)  
          elif ctag == 'NextContinuationToken':
              ctoken = c.text
              again = True

#
main()
