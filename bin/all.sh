#!/bin/bash
PYTHON=python
${PYTHON} basicGeneInfo.py > bgi.json
python agr_validate.py -s gene/geneMetaData.json -d ../AGRdatafeed/bgi.json

${PYTHON} alleleInfo.py > alleles.json
python agr_validate.py  -s allele/alleleMetaData.json -d ../AGRdatafeed/alleles.json

${PYTHON} diseaseAnnotations.py > disease.json
