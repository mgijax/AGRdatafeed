import sys
import xml.etree.ElementTree as ET
import urllib
import re
import json
import argparse

AGRURL='http://download.alliancegenome.org/'

def getOptions ():
  parser = argparse.ArgumentParser(description="List the files that have been uploaded to the Alliance.")
  parser.add_argument('-x', '--taxonid',
    metavar="TAXONID",
    help="Select for this taxon id.")
  parser.add_argument('-d', '--datatype',
    metavar="TYPE",
    choices= ["BGI","ALLELE","EXPRESSION","PHENOTYPE","DAF","GFF"],
    help="Select for this data type. One of: %(choices)s")
  parser.add_argument('-s', '--schema',
    metavar="VERSION",
    help="Select for this schema version.")
  parser.add_argument('-e', '--etag',
    metavar="UID",
    help="Select for this ETag (unique identifier).")
  parser.add_argument('-m', '--modified',
    metavar="DATE",
    help="Select for this modification date. Specify dates using yyyy-mm-dd format (eg '2019-02-03'). You may specify a single date, or a date range using 'mindate..maxdate'. You may leave off maxdate to specify anything on or after mindate (eg '2017-01-01..'), or you may leave off mindate (eg, '..2018-12-31) to spcify anything on or before maxdate.")
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

def main():
    opts = getOptions()
    RE=re.compile(r'[{].*[}]')
    fd = urllib .urlopen(AGRURL)
    s = fd.read()
    fd.close()
    root = ET.fromstring(s)
    for c in root:
      # quick and dirty conversion to dictionary 
      d = dict([(RE.sub('', x.tag), x.text) for x in c])
      # skip anything that's not a file entry
      if not 'ETag' in d:
	continue
      # for some reason, the ETag values include the quote characters
      d['ETag'] = d['ETag'].replace('"','')
      # split the file path into parts
      kparts = d['Key'].split("/")
      # just grab the yyyy-mm-dd part of the modification date
      mdate = d['LastModified'][:10]
      #
      # Apply filters
      #
      if opts.etag and opts.etag != d['ETag']:
	continue
      #
      if opts.taxonid and opts.taxonid not in kparts:
	continue
      #
      if opts.schema and opts.schema not in kparts:
	continue
      #
      if opts.datatype and opts.datatype not in kparts:
	continue
      #
      if opts.modified:
	if opts.modified[0] and mdate < opts.modified[0]:
	  continue
	if opts.modified[1] and mdate > opts.modified[1]:
	  continue
      #
      print json.dumps(d, indent=2)

#
main()
