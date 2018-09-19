## FishLifeExonCapture

# Tutorial for assembling exon capture data across the diversity of fishes

Software required:

aTRAM (requires Python 3)

samtools (most versions should work, but newer versions >2 preferred)

BLAST+ (command line)

bwa

Velvet

Biopython

Exonerate

Mafft

TranslatorX

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

# Step 3: Pre-process trimmed reads for aTRAM

# A note on lustre file systems: Some university servers use lustre scratch file systems. Colonial One at GWU uses this type of file system, and they generally prefer that you run jobs in the /lustre/groups/ workspace. However, lustre systems cannot use SQL-based databases. This is an ongoing issue for data analysis with aTRAM, that I am actively working on. Right now, I usually run jobs from this point in the non-lustre groups space, but the system administrators don't like it.

aTRAM uses SQLite databases and BLAST databases of the raw reads to run, and it can set these up with a script that comes with aTRAM, aTRAM-preprocessor.py. 

You can run this for all of your files like so:

```
module load aTRAM
../FishLifeExonCapture/preprocess4atram.sh
```

# Step 4: Map raw reads back to representative bait sequences

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

# Step 5: Build initial assemblies in Velvet

The previous step should have generated a .fq file for each locus (providing that some reads mapped to the reference sequences ). Now we can generate an initial assembly for each locus with Velvet:

```
../FishLifeExonCapture/initialVelvet.sh
```

# Step 6: Extract longest assembled contig

Velvet often assembles incomplete contigs at this stage, but the longest one should contain some part of the exon we are trying to assemble. To pull out that contig, run:

```
for f in *initial;
do
python ../FishLifeExonCapture/getLongest.py -f $f/contigs.fa -o $f.combined.fa;
done
```

# Step 7: Run aTRAM

This script uses the default parameters for aTRAM, using Velvet again as the assembler. It uses ten iterations, which seems to be sufficient for most loci. If you're trying to assemble something longer, like a full mitogenome, you'll want to run aTRAM directly and change this.

Colonial one has a special python environment that was set up for aTRAM. To run it under this environment, do:

```
module load velvet
module load sqlite
module load aTRAM/2.0
source $aTRAM

../FishLifeExonCapture/runaTRAM.sh
```

# Step 8: Find reading frames and filter exons

Run the thing.
