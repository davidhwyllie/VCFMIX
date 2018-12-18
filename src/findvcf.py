#!/usr/bin/env python
 
# finds vcf files to process for mixtures;
# generates a series of .toprocess files, which allow multiple threads to be assigned to each .toprocess file.

# necessary libraries
import os
import unittest
import uuid
import inspect
import datetime 
import hashlib
import random
import csv
import sys
import shutil
import glob
from vcfScan import lineageScan

# use a list of samples which are provided in a flat file.
# they are TB samples from either B'ham or Brighton.


# read these from the flat file, and divide the guids into approx equal sized files.
to_process=[]
outputdir = os.path.join('..','output')
globpattern="/mnt/microbio/ndm-hicf/ogre/pipeline_output/*/MAPPING/2e6b7bc7-f52c-4649-8538-c984ab3894bb_R00000039/STD/basecalls/*v3.vcf.gz"
v = lineageScan()
inputfiles = glob.glob(globpattern)
print("Found {0} input files".format(len(inputfiles)))
for inputfile in inputfiles:
    guid=os.path.basename(inputfile)[0:36]
    
    # test whether the file has already been parsed
    targetfile = os.path.join(outputdir,'{0}.txt'.format(guid))
    if os.path.exists(targetfile):
        print('exists {0}'.format(guid))
    else:

        res = v.parse(vcffile = inputfile, guid= guid)
        print(guid, res)
        v.region_stats.to_csv(targetfile)
        
