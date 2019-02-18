## FishLifeExonCapture

# Tutorial for assembling exon capture data across the diversity of fishes

Software required:

Trimmomatic

aTRAM (requires Python 3)

samtools 1.8 or above

BLAST+ (command line)

bwa

Velvet

Biopython

Exonerate

Macse

CD-Hit


# Step 1: Organize your fastq files into separate directories for each species

You should establish a main working directory for all of your samples. Starting with a list of demultiplexed fastq files from the sequencer, we need to make separate directories to process our samples through the pipeline. These files are commonly gzip compressed (so you might see .fastq.gz or .fq.gz). 

Make directories for each of the fastq files, and move them to those directories:

```
for f in *_R1.fastq.gz;
do
mkdir ${f%_*.*.*};
mv ${f%_*.*.*}*gz ${f%_*.*.*}/;
done
```

# Step 2: Run Trimmomatic to quality trim the sequences

Run the trimmomatic-loop-PE.sh script in the main working directory. This calls the trimmomatic.jar file in the location where it is stored on Makaira. 

***If you are running this on another system, you may want to make a copy of this file and change the path.***

```
#module load trimmomatic
../FishLifeExonCapture/trimmomatic-loop-PE.sh
```

Note: This script will look for a file called 'adapters.fa' in the main working directory, to trim out adapter contamination. Different sets of adapters are used for different library preparations, so you should check which are appropriate. Trimmomatic has all of these sequences packaged with it, so you can just move this to the 'adapters.fa' file. The adapters are proprietary Illumina sequences, so I have not included them here.

There is a version for single-end files, trimmomatic-loop-SE.sh

When the script is finished, each .fastq.gz file will have an associated .trimmed.fastq.gz file, and two 'rem' files. These are the leftover reads that no longer have mate-pairs.


# Step 3: Map raw reads back to representative bait sequences

aTRAM gets better assemblies than other software, but for its first iteration, it is dependent on a reference sequence. We have many possible reference sequence, and if the reference sequence doesn't recruit reads during the first iteration, nothing assembles.To get around having to choose the right reference, we're generating a starting contig with the reads that mapped to the reference sequences in all_Master.fasta, which contains all the sequences the baits were designed on. This might be redundant in some cases, but tends to produce longer, more accurate assemblies, without having to choose a reference sequence for every taxon.

This script maps the raw reads with bwa against the reference sequences that all exon baits were designed on, as well as coding mitochondrial genes. These sequences are included in the all_Master.fasta file. It then makes individual fastq files for each locus, with PCR duplicates removed.

```
module load samtools/1.8
../FishLifeExonCapture/map-exons.sh
```

Another version of the script is available to work with older versions (<2) of samtools.

For the older set of Otophysi markers (see Arcila et al 2017), use:

```
#module load samtools/1.8
../FishLifeExonCapture/map-exons-Otophysi.sh
```

# Step 4: Build initial assemblies in Velvet

The previous step should have generated a .fq file for each locus (providing that some reads mapped to the reference sequences ). Now we can generate an initial assembly for each locus with Velvet. This script runs the assemblies for .fq files that contain reads and removes empty files. It then pulls out the longest assembled contig to feed to aTRAM.

```
#module load velvet
../FishLifeExonCapture/initialVelvet.sh
```



# Step 5: Run aTRAM

This script uses the default parameters for aTRAM, using Velvet again as the assembler. It uses ten iterations, which seems to be sufficient for most loci. If you're trying to assemble something longer, like a full mitogenome, you'll want to run aTRAM directly and change this.

Colonial One (C1) has a special python environment that was set up for running aTRAM, and these modules are loaded as part of the runaTRAM_c1.sh script. This version also assumes that you are working in the /lustre/ file system, which is not compatible with the sqlite files that aTRAM uses. To get around this, this version copies files to the /scratch/ space on the compute node, then back to the working directory. So if you are using C1, run:

```
../FishLifeExonCapture/runaTRAM_c1.sh
```

If you're not working in a /lustre/ system, and don't need the special python environment for C1, you can run:

```
../FishLifeExonCapture/runaTRAM.sh
```


# Step 6: Find reading frames and filter exons

To efficiently filter the reading frames for the exons, we need to collapse highly similar contigs produced by aTRAM, and then get the reading frame from exonerate. This will produce a exonerate_filtered.fa file for each exon that passes the filters. 

In this current version, the reading frames are percomorph-specific. More reading frame sets will be added, and this script will be modified to reflect that. 

```
module load cd-hit

../FishLifeExonCapture/ExonFiltering.sh
```

There is a second version to deal with the Otophysi set of markers available. The names of these loci start with 'G' instead of 'E'. It works the same way, it just calls a different set of markers and reading frames.

```
module load cd-hit
#module load exonerate

../FishLifeExonCapture/filterOtophysiExons.sh
```

# Step 7: Reading-Frame Aware Alignment

First we need to gather all the separate exon files into a single file that we can use for alignment. The following script will make a new directory called Alignments/, and .unaligned.fasta files for each exon. It's just a bash script, so it doesn't require extra software.

```
../FishLifeExonCapture/preAlignment.sh
```

Now we need to align all of the files in the new Alignments/ directory. Move into the Alignments directory. You will now run scripts out of this directory.

We will use MACSE2 to align our exons, but there are also instructions for using TranslatorX. Just be aware that TranslatorX cannot tolerate any insertions or deletions that cause frame shifts, but MACSE can. 

To run MACSE2:

```
cd Alignments/
../../FishLifeExonCapture/run_macse.sh
```

To run a stand-alone version of TranslatorX and Mafft (which I have named tx.pl, but you may have named something different):
```
cd Alignments/
for f in *.fasta;
do
perl tx.pl -i $f -o $f.tx -p F;
done
```

# Step 8: Alignment Filtering (Optional)

There are several alignment filtering scripts included in this repository. All require biopython.

***AlignmentCleaner.py***

This is a multi-purpose script that cleans:

-Single-taxon insertions (often, though certainly not always, these are caused by assembly errors)

-Gappy edges above a user-specified gap threshold

-Short sequences that fall below a user-specified threshold

If you want to clean out edges composed of more than 60% gaps, single-taxon insertions, and sequences that cover less than 50% of the alignment for all MACSE2 alignments, you could run:

```
# Nuclear Exons aligned with MACSE2
for f in E*NT_aligned.fasta;
do
python ../../FishLifeExonCapture/AlignmentCleaner.py -f $f -o $f.cleaned.fasta -c 0.5 -t 0.6;
done
```


Other Utilities:

If I want to drop the  taxa from an alignment for any reason, I use the dropTaxa.py script. It just takes a text file of names to remove from a fasta file (one name per line), the fasta file to prune, and the name of an output file. 

```
python ../../FishLifeExonCapture/dropTaxa.py -f E0001.NT_aligned.fasta -o E0001.NT_aligned.dropped.fasta -t badTaxaList.txt;
```

You can run the AlignmentCleaner.py script on the output of dropTaxa.py, to remove any gaps left in the alignment.


If you want to flag divergent sequences, you can use the AlignmentChecker.py script. This script flags taxa that fall above a certain user-specified distance threshold, and prints taxon names that fall above this threshold. Your mileage may vary depending on how divergent the taxa are in your alignment.


Other notes:

I always look at the exons that are numbered higher than E1730 by eye. These exons have known paralogs, but were included in the dataset to connect with older PCR-based datasets. They are usually fine, but you can never be to careful.

My python scripts always use fasta-formatted files because they are easier for me to work with. If you need other formats for phylogenetic analysis, the AlignmentConverter.py script is quite flexible. Just note that to specify a relaxed phylip format, you'll need to type phylip-relaxed (otherwise, it's strict phylip).


# Step 9: Clean Up (Optional)

All of this data processing leaves a lot of intermediate files for each sample. I tend to keep them for debugging purposes, but at some point they no longer are necessary. You can run a script to remove these intermediate files, but still keep the output of trimmomatic, aTRAM, and the exon filtering.

```
../FishLifeExonCapture/CleanUp.sh
```


