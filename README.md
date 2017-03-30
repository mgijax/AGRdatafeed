# AGRdatafeed

Prepare MGI data to be ingested by AGR.

To Use:
cp config.cfg.default config.cfg

alter config.cfg as needed - most likey set the MouseMine URL

run setup script - installs needed 3rd party libraries
(only need to run it once)

run rundatafeed
to generate all (both so far) data feed output files

OR to get individual data feed output files:
python diseaseAnnotations.py
python basicGeneInfo.py
