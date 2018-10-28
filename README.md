## FishLifeExonCapture

# Tutorial for assembling exon capture data across the diversity of fishes

Software required:

Trimmomatic

aTRAM (requires Python 3)

samtools (older versions should work, but newer versions >2 preferred)

BLAST+ (command line)

bwa

Velvet

Biopython

Exonerate

Mafft

TranslatorX (standalone, or preferred translation-aware aligner)

CD-Hit

FastUniq


# Step 1: Organize your fastq files into separate directories for each species

You should establish a main working directory for all of your samples. Starting with a list of demultiplexed fastq files from the sequencer, we need to make separate directories to process our samples through the pipeline. These files are commonly gzip compressed (so you might see .fastq.gz or .fq.gz). 

Make directories for each of the fastq files, and move them to those directories:

```
for f in *_R1.fastq.gz;
do
mkdir ${f%.*};
mv $f ${f%.*}/;
done
```

# Step 2: Run Trimmomatic to quality trim the sequences

Run the trimmomatic-loop-PE.sh script in the main working directory. This calls the trimmomatic.jar file in the location where it is stored on Colonial One. 

***If you are running this on another system, you may want to make a copy of this file and change the path.***

```
module load trimmomatic
../FishLifeExonCapture/trimmomatic-loop-PE.sh
```

Note: This script will look for a file called 'adapters.fa' in the main working directory, to trim out adapter contamination. Different sets of adapters are used for different library preparations, so you should check which are appropriate. Trimmomatic has all of these sequences packaged with it, so you can just move this to the 'adapters.fa' file. The adapters are proprietary Illumina sequences, so I have not included them here.

There is a version for single-end files, trimmomatic-loop-SE.sh

When the script is finished, each .fastq.gz file will have an associated .trimmed.fastq.gz file, and two 'rem' files. These are the leftover reads that no longer have mate-pairs.


# Step 3: Map raw reads back to representative bait sequences

aTRAM gets better assemblies than other software, but for its first iteration, it is dependent on a reference sequence. We have many possible reference sequence, and if the reference sequence doesn't recruit reads during the first iteration, nothing assembles.To get around having to choose the right reference, we're generating a starting contig with the reads that mapped to the reference sequences in all_Master.fasta, which contains all the sequences the baits were designed on. This might be redundant in some cases, but tends to produce longer, more accurate assemblies, without having to choose a reference sequence for every taxon.

This script maps the raw reads with bwa against the reference sequences that all exon baits were designed on, as well as coding mitochondrial genes. These sequences are included in the all_Master.fasta file. It then makes individual fastq files for each locus.

```
module load samtools
module load bwa
../FishLifeExonCapture/map-exons.sh
```

Another version of the script is available to work with older versions (<2) of samtools.

For the older set of Otophysi markers (see Arcila et al 2017), use:

```
module load samtools
module load bwa
../FishLifeExonCapture/map-exons-Otophysi.sh
```

# Step 4: Build initial assemblies in Velvet

The previous step should have generated a .fq file for each locus (providing that some reads mapped to the reference sequences ). Now we can generate an initial assembly for each locus with Velvet. This script runs the assemblies for .fq files that contain reads and removes empty files. It then pulls out the longest assembled contig to feed to aTRAM.

```
module load velvet
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
module load exonerate

../FishLifeExonCapture/ExonFiltering.sh
```

There is a second version to deal with the Otophysi set of markers available. The names of these loci start with 'G' instead of 'E'. It works the same way, it just calls a different set of markers and reading frames.

```
module load cd-hit
module load exonerate

../FishLifeExonCapture/filterOtophysiExons.sh
```

# Step 7: First Alignment

First we need to gather all the separate exon files into a single file that we can use for alignment. The following script will make a new directory called Alignments/, and .unaligned.fasta files for each exon. It's just a bash script, so it doesn't require extra software.

```
../FishLifeExonCapture/preAlignment.sh
```

Now we need to align all of the files in the new Alignments/ directory. I didn't make a special script for this, since people have many different alignment preferences. I typically use a standalone version of translatorX and mafft for this task. My script might look something like this:

```
cd Alignments/
for f in *.fasta;
do
perl tx.pl -i $f -o $f.tx -p F;
done
```

# Step 8: Alignment Filtering

Even with care, things we don't want in our alignments get through. There are three python scripts built to filter out some of the noisy sequences and scaffolding Ns that sometimes get into the alignments. Both require Biopython. Use '-h' to see the arguments and default settings of each script.

The first is called AlignmentChecker.py. This script builds a distance matrix of your alignment and outputs taxon names with average distances that are above a user-specified threshold. I am working on alternative versions of this script to deal with proteins, and possibly speed up the run time. The threshold you specify may vary depending on the level of diversity in your matrix. In an alignment of mostly one family, I would set this somewhat lower than I would for an alignment that included many orders. 

If I wanted this script to write a file with the names of the taxa that it flags, I might do this:

```
for f in *tx_nt_ali.fasta;
do
python ../../FishLifeExonCapture/AlignmentChecker.py -f $f -d 0.5 > $f.badTaxa.txt;
done
```

If it's flagging a lot of taxa, or a lot of different alignments, you'll need to investigate further. It might be that you have a mislabeled taxon that is quite divergent from the rest of your sequences.

If I want to prune the flagged taxa from the alignment, I use the dropTaxa.py script. It just takes a text file of names to remove from a fasta file (one name per line), the fasta file to prune, and the name of an output file.

```
for f in *tx_nt_ali.fasta;
do
python ../../FishLifeExonCapture/dropTaxa.py -f $f -o $f.dropped.fasta -t $f.badTaxa.txt;
done
```

Finally, removing these taxa often leaves gaps in the alignment, or large gaps are created from scaffolding Ns, or just single-taxon insertions. The gapCleaner.py script removes these columns from the alignment. It will also remove sequences that cover less than 50% of the alignment. If you want to change this threshold, you can change the fraction with the '-c' flag.


```
for f in *dropped.fasta;
do
python ../../FishLifeExonCapture/gapCleaner.py -f $f -o $f.cleaned.fasta -c 0.5;
done
```

I always look at the exons that are numbered higher than E1730 by eye. These exons have known paralogs, but were included in the dataset to connect with older PCR-based datasets. They are usually fine, but you can never be to careful.

My python scripts always use fasta-formatted files because they are easier for me to work with. If you need other formats for phylogenetic analysis, the AlignmentConverter.py script is quite flexible. Just note that to specify a relaxed phylip format, you'll need to type phylip-relaxed (otherwise, it's strict phylip).



