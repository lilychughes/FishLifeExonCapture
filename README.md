## FishLifeExonCapture

# Tutorial for assembling exon capture data across the diversity of fishes

Software required:

Trimmomatic

aTRAM (requires Python 3)

samtools 1.8 or above

BLAST+ (command line)

bwa

Velvet

Trinity

Biopython

Exonerate

Macse

CD-Hit


# Step 0: Organize your fastq files into separate directories for each species

You should establish a main working directory for all of your samples. I'll call this 'Project_Directory' here, but go ahead and name it something meaningful to you. The 'Project_Directory' and the 'FishLifeExonCapture' folder should be next to each other.

Move your raw fastq files into the Project_Directory, then cd into the Project_Directory.

Starting with a list of paired-end demultiplexed fastq files from the sequencer, we need to make separate directories to process our samples through the pipeline. These files are commonly gzip compressed (so you might see .fastq.gz or .fq.gz). 

Rename these fastq files as you want your samples to be named in your alignments, before the "\_R1.fastq.gz" endings. __Don't use special characters or spaces__, but letters, numbers, and underscores are ok. I typically name my files like this:

Family_Genus_species_S1234_R1.fastq.gz
Family_Genus_species_S1234_R2.fastq.gz

Make directories for each of the fastq files, and move them to those directories:

```
for f in *_R1.fastq.gz;
do
mkdir ${f%_*.*.*};
mv ${f%_*.*.*}*gz ${f%_*.*.*}/;
done
```

# Step 1: Run Trimmomatic to quality trim the sequences

Run the trimmomatic-loop-PE.sh script in the main working directory. You will need to edit the path to the Trimmomatic jar file for your particular system.

***If you are running this on another system, you may want to make a copy of this file and change the path to the trimmomatic jar file*** 
```
../FishLifeExonCapture/trimmomatic-loop-PE.sh
```

Note: This script will look for a file called 'adapters.fa' in the Project_Directory, to trim out adapter contamination. Different sets of adapters are used for different library preparations, so you should check which are appropriate. Trimmomatic has all of these sequences packaged with it, so you can just move this to the 'adapters.fa' file. The adapters are proprietary Illumina sequences, so I have not included them here.

When the script is finished, each .fastq.gz file will have an associated .trimmed.fastq.gz file, and two 'rem' files. These are the leftover reads that no longer have mate-pairs.


# Step 2: Map raw reads back to representative bait sequences

aTRAM gets better assemblies than other software, but for its first iteration, it is dependent on a reference sequence. We have many possible reference sequence, and if the reference sequence doesn't recruit reads during the first iteration, nothing assembles.To get around having to choose the right reference, we're generating a starting contig with the reads that mapped to the reference sequences in all_Master.fasta, which contains all the sequences the baits were designed on. This might be redundant in some cases, but tends to produce longer, more accurate assemblies, without having to choose a reference sequence for every taxon.

This script maps the raw reads with bwa against the reference sequences that all exon baits were designed on, as well as coding mitochondrial genes. These sequences are included in the all_Master.fasta file. It then makes individual fastq files for each locus, with PCR duplicates removed.

***You will need bwa,samtools 1.7 or higher in your path. You will also need biopython and python2.7.***
```
../FishLifeExonCapture/map-exons.sh
```

For the older set of Otophysi markers (see Arcila et al 2017), use:

***You will need bwa,samtools 1.7 or higher in your path. You will also need biopython and python2.7.***
```
../FishLifeExonCapture/map-exons-Otophysi.sh
```

# Step 3: Build initial assemblies in Velvet

The previous step should have generated a .fq file for each locus (providing that some reads mapped to the reference sequences ). Now we can generate an initial assembly for each locus with Velvet. This script runs the assemblies for .fq files that contain reads and removes empty files. It then pulls out the longest assembled contig to feed to aTRAM.


***You will need velvet in your path. You will also need biopython and python2.7.***
```
../FishLifeExonCapture/initialVelvet.sh
```



# Step 4: Run aTRAM

This script uses the default parameters for aTRAM, using Trinity as the assembler. It uses five iterations, which seems to be sufficient for most loci. If you're trying to assemble something longer, like a full mitogenome, you may want to try MitoBIM or similar.

***Both scripts below use 6 CPUs with Trinity. You should be able to request this directly from your scheduling software***


***Requires blast+, aTRAM 2.0 and its dependencies, trinity, sqlite in your path.***

```
../FishLifeExonCapture/runaTRAM.sh
```




Colonial One (the GW cluster) prefers that we run our jobs out of the /lustre/ partition, which is incompatible with the sqlite files that aTRAM uses. So for now, there's a separate script for running aTRAM in this environment that moves the files to the /scratch/ directory of the compute nodes.

***Requires blast+, aTRAM 2.0 and its dependencies, trinity, sqlite in your path.***

```
../FishLifeExonCapture/runaTRAM_c1.sh
```




# Step 5: Find reading frames and filter exons

To efficiently filter the reading frames for the exons, we need to collapse identical contigs produced by aTRAM with CD-HIT, and then get the reading frame from exonerate. This will produce a exonerate_filtered.final_contigs.fa file for each exon that passes the filters. For a contig to pass to this final_contigs.fa file, it must have the correct reading frame as determined by exonerate. If only one contig assembled with the reading frame, the exon sequence passes directly to the final_contigs.fa file. If more than one contig assembled with the reading frame, the reading frames are compared with CD-HIT, and if they are 99% similar, the longer one will be passed onto the final_contigs.fa file. If the sequences are more divergent, they do not pass this filter.


There are several versions of reference reading frames to use with Exonerate. More can be added upon request.


***You will need cd-hit and exonerate in your path. You will also need biopython and python2.7.***


For percomorph fishes:
```
../FishLifeExonCapture/ExonFilteringPercomorph.sh
```

For elopomorph fishes:
```
../FishLifeExonCapture/ExonFilteringOsteoglossomorph.sh
```

For osteoglossomorph fishes:
```
../FishLifeExonCapture/ExonFilteringOsteoglossomorph.sh
```


There is a second version to deal with the Otophysi set of markers available (see Arcila et al. 2017). The names of these loci start with 'G' instead of 'E'. It works the same way, it just calls a different set of markers and reading frames.



***You will need cd-hit and exonerate in your path. You will also need biopython and python2.7.***
```
../FishLifeExonCapture/ExonFilteringOtophysi.sh
```

Note: Additional filtering references will be added for more divergent groups shortly! Stay tuned!


# Step 5b: Filter Exons and Flanking Introns (Optional)

This is a new feature. You must run this AFTER running Step 5. 

Instead of only including the portion of the assembled sequence that matches to the reference reading, it includes the entire contig. If more than one contig assembled with the reading frame, the contigs are compared with CD-HIT, and if they are 98% similar, the longer one will be passed to a file ending in .filtered_flanks.fa.

```
../FishLifeExonCapture/FlankFlitering.sh
```


# Step 6: Reading-Frame Aware Alignment

First we need to gather all the separate exon files into a single file that we can use for alignment. The following script will make a new directory called Alignments/, and .unaligned.fasta files for each exon. It's just a bash script, so it doesn't require extra software.

```
../FishLifeExonCapture/preAlignment.sh
```

An alternative version exists for the older set of Otophysan markers:

```
../FishLifeExonCapture/preAlignment_Otophysi.sh
```

Now we need to align all of the files in the new Alignments/ directory. Move into the Alignments directory. You will now run scripts out of this directory.

We will use MACSE2 to align our exons, but there are also instructions for using TranslatorX. Just be aware that TranslatorX cannot tolerate any insertions or deletions that cause frame shifts, but MACSE can. 

To run MACSE2:

***You may need to replace the full path to the MACSE jar file in this script.***

```
cd Alignments/
../../FishLifeExonCapture/run_macse.sh
```

As with the last step, there is an alternative script for the otophysan markers:

```
cd Alignments/
../../FishLifeExonCapture/run_macse_Otophysi.sh
```


Alternatively, if you would prefer to run TranslatorX, some options are below:
To run a stand-alone version of TranslatorX and Mafft (which I have named tx.pl, but you may have named something different):
```
cd Alignments/
for f in *.fasta;
do
perl tx.pl -i $f -o $f.tx -p F;
done
```

# Step 6b: Alignment for contigs with flanking regions

At this time, I don't use a reading-frame-aware aligner to deal with these sequences. The script below will gather the contigs with flanking regions that passed all filters into files that can be aligned for phylogenetic analysis.

```
../FishLifeExonCapture/preAlignment_Flanks.sh
```

This will create a folder called Alignments_Flanks


# Step 7: Alignment Filtering (Optional)

There are several alignment filtering scripts included in this repository. 

***All require biopython and are written for python2.7.***

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

My python scripts always use fasta-formatted files because they are easier for me to work with. If you need other formats for phylogenetic analysis, the AlignmentConverter.py script is quite flexible. Just note that to specify a relaxed phylip format, you'll need to type phylip-relaxed (otherwise, it's strict phylip).


# Step 8: Clean Up (Optional)

All of this data processing leaves a lot of intermediate files for each sample. I tend to keep them for debugging purposes, but at some point they no longer are necessary. You can run a script to remove these intermediate files, but still keep the output of trimmomatic, aTRAM, and the exon filtering.

```
../FishLifeExonCapture/CleanUp.sh
```


