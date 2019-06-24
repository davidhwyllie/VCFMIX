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
