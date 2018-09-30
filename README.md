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


# Step 1: Organize your fastq files into separate directories for each species

You should establish a main working directory for all of your samples. Starting with a list of demultiplexed fastq files from the sequencer, we need to make separate directories to process our samples through the pipeline. These files are commonly gzip compressed (so you might see .fastq.gz or .fq.gz). If you need to expand the files you can run the following:

```
gunzip *gz
```

Now make directories for each of the fastq files, and move them to those directories:

```
for f in *fastq;
do
mkdir ${f%.*};
mv $f ${f%.*}/;
done
```

# Step 2: Run Trimmomatic to quality trim the sequences

Run the trimmomatic-loop.sh script in the main working directory. This calls the trimmomatic.jar file in the location where it is stored on Colonial One. If you are running this on another system, you may want to make a copy of this file and change the path.

```
module load trimmomatic
../FishLifeExonCapture/trimmomatic-loop.sh
```

Note: This script will look for a file called 'adapters.fa' in the main working directory, to trim out adapter contamination. Different sets of adapters are used for different library preparations, so you should check which are appropriate. Trimmomatic has all of these sequences packaged with it, so you can just move this to the 'adapters.fa' file. The adapters are proprietary Illumina sequences, so I have not included them here.

When the script is finished, each .fastq file will have an associated .trimmed.fastq file, and the originial untrimmed file is compressed.


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
../FishLifeExonCapture/map-exons-Otophysi.sh
```

# Step 4: Build initial assemblies in Velvet

The previous step should have generated a .fq file for each locus (providing that some reads mapped to the reference sequences ). Now we can generate an initial assembly for each locus with Velvet. This script runs the assemblies for .fq files that contain reads and removes empty files. It then pulls out the longest assembled contig to feed to aTRAM.

```
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

aTRAM tends to write all of the contigs it found in each of its iterations, regardless of whether they are identical or not. To reduce the number of contigs in the aTRAM output file, I run CD-HIT, a software that clusters sequences based on some identity threshold, and writes the longest sequences to a file. It is also quite fast. 

I set the cluster identity to 0.98 below, but you can change it where it says -c.

If I'm in my working directory with all of my species directories, I can run:

```
module load cd-hit

for directory in *;
do
if [  -d $directory  ];
then
cd $directory;
for f in *.filtered_contigs.fasta;
do
cd-hit-est -i $f -o ${f%.*}.cdhit -c 0.98;
done;
cd ../;
fi;
done
```

