# FishLifeExonCapture

# Tutorial for assembling exon capture data across the diversity of fishes

Software required:

aTRAM

samtools (most versions should work)

BLAST+ (command line)

bwa

Biopython

Exonerate

Mafft

TranslatorX

# Step 1: Organize your fastq files into separate directories for each species

Starting with a list of demultiplexed fastq files from the sequencer, we need to make separate directories to process our samples through the pipeline. These files are commonly gzip compressed (so you might see .fastq.gz). If you need to expand the files you can run the following:

```
for f in *gz;
do
gunzip $f;
done
```

