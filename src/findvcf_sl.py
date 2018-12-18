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

# read the reference
genbank_file_name = os.path.join("..", "testdata", "NC_0103971.1.gb")
rs = regionScan_from_genbank(genbank_file_name)

# record the regions studied	
output_excel = os.path.join('..','testdata', 'regions.xlsx')
rs.regions.to_excel(output_excel)

# define the output directory
outputdir = os.path.join('..','output')

# identify the files to study
to_process=[]
#globpattern="/mnt/microbio/ndm-hicf/ogre/pipeline_output/*/MAPPING/2a836c1b-c9d1-4c16-ab2d-2c9d26c28f3e_R00000277/STD/basecalls/*v3.vcf.gz"
#inputfiles = glob.glob(globpattern)
inputfiles = ["/mnt/microbio/ndm-hicf/ogre/pipeline_output/d9b0b4dc-1f92-4ba6-a850-9649647ff792/MAPPING/2a836c1b-c9d1-4c16-ab2d-2c9d26c28f3e_R00000277/STD/basecalls/d9b0b4dc-1f92-4ba6-a850-9649647ff792_v3.vcf.gz"]

print("Found {0} input files".format(len(inputfiles)))
for inputfile in inputfiles:
    guid=os.path.basename(inputfile)[0:36]
    
    # test whether the file has already been parsed
    targetfile = os.path.join(outputdir,'{0}.txt'.format(guid))
    if os.path.exists(targetfile):
        print('exists {0}'.format(guid))
    else:
        res = rs.parse(vcffile = inputfile, guid= guid)
        print(guid, res)
        rs.region_stats.to_csv(targetfile)
        
