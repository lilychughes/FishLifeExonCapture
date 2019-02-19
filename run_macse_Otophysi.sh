#!/bin/sh

#change the path to the macse_v2.03.jar file as necessary!

#Run macse2 for nuclear alignments
for sequences in G*.unaligned.fasta;
do
    if [  ! -e ${sequences%.*.*}.macse.txt  ];
	then
	echo macse alignment started > ${sequences%.*.*}.macse.txt;
	java -jar ../../macse_v2.03.jar -prog trimNonHomologousFragments -seq $sequences -min_homology_to_keep_seq 0.4 -out_NT ${sequences%.*}.NT_trimNonHomologousFragments.fasta -out_AA  ${sequences%.*}.AA_trimNonHomologousFragments.fasta;
	java -jar ../../macse_v2.03.jar -prog alignSequences -seq ${sequences%.*}.NT_trimNonHomologousFragments.fasta -out_NT ${sequences%.*.*}.NT_aligned.fasta -out_AA ${sequences%.*.*}.AA_aligned.fasta;
	fi;
done

# replace '!' insertion character with 'N', since it seems to interfere with other software

sed -i 's/!/N/g' *NT_aligned*