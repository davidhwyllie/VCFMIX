# Find sites at which multiple high-quality sequences map within VCF files
Allows identification of, and quantification of mixtures of, high quality bases within VCF files.  It is designed for use with reference mapped bacterial sequences.

## Background
Extraction, summarising and depiction of mixed base frequencies at lineage defining positions in *Mycobacterium tuberculosis* has application to laboratory quality control and the detection of mixed sequences.
Techniques for detecting inter-lineage mixtures have been [published](https://www.ncbi.nlm.nih.gov/pubmed/30209183).  This module computes the F2 and F47 statistics described in the paper.

In cases where the phylogeny is unknown, or where variation exists at only a few bases in closely related samples, approaches exist which compare multiple samples to identify whether the positions to which mixtures of bases are mapped are likely to be of technical, or biological origin.  One such is the [mixPore algorithm](to add link to paper) implemented as part of the [findNeighbour3 server](https://github.com/davidhwyllie/findNeighbour3).

This software has predominantly been used to consume data from a [bioinformatic pipeline used for TB processing in England](https://github.com/oxfordmmm/CompassCompact).

## Implementation

VCFMIX provides python classes to extract high quality bases in defined positions, or in all positions, within VCF files.

The core class is *vcfScan*.    
The *add_roi* methods allows definition of a region of interest within which variation should be summarised.  
The *parse* method generates a summary for all defined rois within a vcf file.

More specialised classes inherit from :vcfScan:
* *lineageScan* defines variation in a set of lineage defining positions, e.g. those of [Coll F *et al* 2014](https://www.ncbi.nlm.nih.gov/pubmed/25176035).
* *FastaMixtureMarker* takes the output from vcfScan and uses it to modify a fasta files containing consensus basecalls, e.g. as produced by [bioinformatic pipeline used for TB processing in England](https://github.com/oxfordmmm/CompassCompact), replacing the called consensus with the relevant [IUPAC code for the mixture discerned](https://www.bioinformatics.org/sms/iupac.html).
If you are looking for the *regionScan_from_genbank* used to perform the calculations described in our publication on [Adaptive masking](https://www.ncbi.nlm.nih.gov/pubmed/29875188), a strategy for assessing whether mis-mapping of non-target bacterial DNA is occurring, this has been moved to the [adaptivemasking module](https://github.com/davidhwyllie/adaptivemasking).

### Obtaining test data
Note that test data in the example below is not included when the github project is cloned, because it is larger than github's limit.
you need to download a test data file manually.

```
cd testdata
wget --no-check-certificate -O 52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz "https://ora.ox.ac.uk/objects/uuid:5e4ec1f8-e212-47db-8910-161a303a0757/download_file?file_format=x-tar&safe_filename=52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz&type_of_work=Dataset"
```

### Example using lineageScan to compute F2 and F47 statistics
( see also [scientific rationale](https://www.ncbi.nlm.nih.gov/pubmed/30209183) )

```
import os
from vcfmix import lineageScan

# create a lineagescan object;
v = lineageScan()


# scan an input file
inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
res = v.parse(vcffile = inputfile, guid='528')

# print details of the regions scanned
print(v.regions_stats)

# export details of the regions scanned
outputfile = os.path.join('..','unitTest_tmp','528.txt')
v.region_stats.to_csv(outputfile)

# compute F2 and F47 statistics (see publication)
summary1 = v.f_statistics()

```

### Example using lineageScan and FastaMixtureMasker to generate a consensus fasta file
( see also [scientific rationale]( http://biorxiv.org/cgi/content/short/681502v1); the output is an input for the [findNeighbour3 server](https://github.com/davidhwyllie/findNeighbour3)).


The below example is in /src/fmm_example.py.

```
""" Example of use of FastaMixtureMarker """
import os
import pathlib
from vcfScan import FastaMixtureMarker, vcfScan
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


print('set up vcfScan object.. (only needs to be done once)')
# use the baseCounts4 tag to identify high quality bases
# don't report minor variants present at less that 5% frequency (done here simply to speed up computations)
v = vcfScan(expectedErrorRate = 0.001, infotag = 'BaseCounts4', report_minimum_maf = 0.05, compute_pvalue = False)     

# we define one region for each base of the genome
for i in range(4411532):     
    v.add_roi(str(1+i),set([1+i]))
    
print('parsing vcf')
inputfile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
if not os.path.exists(inputfile):
    self.fail("Input file does not exist.  Please see README.  You may need to install test data.")
v.parse(vcffile = inputfile)
print("Parse complete; writing output")

# make sure a target directory exists
targetdir = os.path.join("..", "unitTest_tmp")		# a writeable directory
pathlib.Path(targetdir).mkdir(parents=True, exist_ok=True)

# write mixed bases to a csv file   
mixfile = os.path.join(targetdir,'output_table.txt')
v.bases.to_csv(mixfile, index=None)

# the fasta file used contains high confidence single base base calls.  Example is from PHE TB pipeline https://github.com/oxfordmmm/CompassCompact
fastafile=os.path.join("..",'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.fasta')
fmm = FastaMixtureMarker(expectedErrorRate = 0.001, mlp_cutoff=6.65, clustering_cutoff = 10, min_maf=0)			
df, seq = fmm.mark_mixed(fastafile, mixfile)

iupac = ['A','C','G','T','r','R','w','W','y','Y','m','M','s','S','k','K']
resDict={}
for item in iupac:
    resDict[item] =  seq.count(item)
    print("There were {1} mixed bases of type {0}".format(item, resDict[item]))

# write fasta
record = SeqRecord(Seq(seq,
                       generic_dna),
                   id="test_id", name="fmm_output_test",
                   description="consensus sequence with mixed bases recorded as iupac codes")
fasta_output_filename = os.path.join(targetdir, 'test_output.fasta')
with open(fasta_output_filename,'w') as f:
    SeqIO.write(record, f, 'fasta')

```