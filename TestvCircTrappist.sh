#!/bin/bash

#launching the test of vCircTrappist. It should work if all dependencies were correctly installed and your computer setup is sufficient

bash ~/vCircTrappist/vCircTrappist_full.sh -F ./TestvCircTrappist/alv.fasta -G ./TestvCircTrappist/alv.gff -Q ./TestvCircTrappist/aln_virus.fastq -S R

rm ./aln_bsj_out.txt
rm aln_bsj_sites.txt
rm aln_circ.sam
rm aln_kmer_compare.csv
rm aln_new_count.csv
rm aln_virus.bam
rm aln_virus.sam
rm aln_virus_coverage.csv
rm aln_virus_sorted.bam
rm circ_antisense_nonU2.bam
rm circ_antisense_nonU2.sam
rm circ_antisense_nonU2_a.sam
rm circ_antisense_nonU2_sorted.bam
rm circ_antisense_nonU2_coverage.csv
rm circ_antisense_U2.bam
rm circ_antisense_U2.sam
rm circ_antisense_U2_a.sam
rm circ_antisense_U2_sorted.bam
rm circ_antisense_U2_coverage.csv
rm circ_sense_nonU2.bam
rm circ_sense_nonU2.sam
rm circ_sense_nonU2_a.sam
rm circ_sense_nonU2_sorted.bam
rm circ_sense_nonU2_coverage.csv
rm circ_sense_U2.bam
rm circ_sense_U2.sam
rm circ_sense_U2_a.sam
rm circ_sense_U2_sorted.bam
rm circ_sense_U2_coverage.csv
rm aln_verifa.txt
rm aln_split.sam
rm aln_circ_list.sam
