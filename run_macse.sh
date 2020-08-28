#!/bin/sh

#change the path to the macse_v2.03.jar file as necessary below:
macsepath="/path/to/macse.jar"

#Run macse2 for nuclear alignments
for sequences in E*.unaligned.fasta;
do
    if [  ! -e ${sequences%.*.*}.macse.txt  ];
	then
	echo macse alignment started > ${sequences%.*.*}.macse.txt;
	java -jar $macsepath -prog trimNonHomologousFragments -seq $sequences -min_homology_to_keep_seq 0.4 -out_NT ${sequences%.*}.NT_trimNonHomologousFragments.fasta -out_AA  ${sequences%.*}.AA_trimNonHomologousFragments.fasta;
	java -jar $macsepath -prog alignSequences -seq ${sequences%.*}.NT_trimNonHomologousFragments.fasta -out_NT ${sequences%.*.*}.NT_aligned.fasta -out_AA ${sequences%.*.*}.AA_aligned.fasta;
# replace '!' insertion character with 'N', since it seems to interfere with other software
	sed -i 's/!/N/g' ${sequences%.*.*}.NT_aligned.fasta;
	fi;
done


#Run macse2 for mitochondrial alignments

for sequences in ND*unaligned.fasta C*unaligned.fasta AT*unaligned.fasta;
do
    if [  ! -e ${sequences%.*.*}.macse.txt  ];
	then
	echo macse alignment started > ${sequences%.*.*}.macse.txt;
	java -jar $macsepath -prog trimNonHomologousFragments -seq $sequences -min_homology_to_keep_seq 0.4 -out_NT ${sequences%.*}.NT_trimNonHomologousFragments.fasta -out_AA  ${sequences%.*}.AA_trimNonHomologousFragments.fasta -gc_def 2;
	java -jar $macsepath -prog alignSequences -seq ${sequences%.*}.NT_trimNonHomologousFragments.fasta -out_NT ${sequences%.*.*}.NT_aligned.fasta -out_AA ${sequences%.*.*}.AA_aligned.fasta -gc_def 2;
# replace '!' insertion character with 'N', since it seems to interfere with other software
	sed -i 's/!/N/g' ${sequences%.*.*}.NT_aligned.fasta;
	fi;
done