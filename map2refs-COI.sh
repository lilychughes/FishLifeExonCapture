#!/bin/bash


# requires bwa and samtools

for f in *.fastq;
do
bwa mem COI_index $f | samtools view -F 4 -bS -o $f.mapped.bam -;
samtools sort $f.mapped.bam sorted.bam;
samtools index $f.mapped.sorted.bam;
samtools view -b $f "COI_Syngnathus_fuscus" "COI_Synchiropus_atrilabiatus" "COI_Fistularia_petimba" "COI_Parupeneus_porphyreus" "COI_Pegasus_laternarius" "COI_Dactyloptena_orientalis" "COI_Thunnus_tonggol" "COI_Bramidae_sp" "COI_Cubiceps_pauciradiatus" "COI_Dysalotus_alcocki" "COI_Caranx_hippos" "COI_Toxotes_jaculatrix" "COI_Pleuronectidae_sp" "COI_Sphyraena_flavicauda" "COI_Acanthurus_reversus" "COI_Cottus_rhotheus" "COI_Lutjanus_vitta" "COI_Rhinecanthus_rectangulus" "COI_Melanotaenia_papuae" "COI_Acanthemblemaria_spinosa" "COI_Chanda_nama" "COI_Nandus_nandus" - \
| samtools bam2fq - > $f.COI.fq
done