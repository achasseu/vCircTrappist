#!/bin/bash

#launching the test of vCircTrappist. It should work if all dependencies were correctly installed and your computer setup is sufficient
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
echo "Script directory: $SCRIPT_DIR"

bash $SCRIPT_DIR/vCircTrappist_full.sh -F $SCRIPT_DIR/TestvCircTrappist/alv.fasta -G $SCRIPT_DIR/TestvCircTrappist/alv.gff -Q $SCRIPT_DIR/TestvCircTrappist/data.fastq -S R -O $SCRIPT_DIR/TestvCircTrappist

rm $SCRIPT_DIR/TestvCircTrappist/aln_bsj_out.txt
rm $SCRIPT_DIR/TestvCircTrappist/aln_bsj_sites.txt
rm $SCRIPT_DIR/TestvCircTrappist/aln_circ.sam
rm $SCRIPT_DIR/TestvCircTrappist/aln_kmer_compare.csv
rm $SCRIPT_DIR/TestvCircTrappist/aln_new_count.csv
rm $SCRIPT_DIR/TestvCircTrappist/aln_virus.bam
rm $SCRIPT_DIR/TestvCircTrappist/aln_virus.sam
rm $SCRIPT_DIR/TestvCircTrappist/aln_virus_coverage.csv
rm $SCRIPT_DIR/TestvCircTrappist/aln_virus_sorted.bam
rm $SCRIPT_DIR/TestvCircTrappist/circ_antisense_nonU2.bam
rm $SCRIPT_DIR/TestvCircTrappist/circ_antisense_nonU2.sam
rm $SCRIPT_DIR/TestvCircTrappist/circ_antisense_nonU2_a.sam
rm $SCRIPT_DIR/TestvCircTrappist/circ_antisense_nonU2_sorted.bam
rm $SCRIPT_DIR/TestvCircTrappist/circ_antisense_nonU2_coverage.csv
rm $SCRIPT_DIR/TestvCircTrappist/circ_antisense_U2.bam
rm $SCRIPT_DIR/TestvCircTrappist/circ_antisense_U2.sam
rm $SCRIPT_DIR/TestvCircTrappist/circ_antisense_U2_a.sam
rm $SCRIPT_DIR/TestvCircTrappist/circ_antisense_U2_sorted.bam
rm $SCRIPT_DIR/TestvCircTrappist/circ_antisense_U2_coverage.csv
rm $SCRIPT_DIR/TestvCircTrappist/circ_sense_nonU2.bam
rm $SCRIPT_DIR/TestvCircTrappist/circ_sense_nonU2.sam
rm $SCRIPT_DIR/TestvCircTrappist/circ_sense_nonU2_a.sam
rm $SCRIPT_DIR/TestvCircTrappist/circ_sense_nonU2_sorted.bam
rm $SCRIPT_DIR/TestvCircTrappist/circ_sense_nonU2_coverage.csv
rm $SCRIPT_DIR/TestvCircTrappist/circ_sense_U2.bam
rm $SCRIPT_DIR/TestvCircTrappist/circ_sense_U2.sam
rm $SCRIPT_DIR/TestvCircTrappist/circ_sense_U2_a.sam
rm $SCRIPT_DIR/TestvCircTrappist/circ_sense_U2_sorted.bam
rm $SCRIPT_DIR/TestvCircTrappist/circ_sense_U2_coverage.csv
rm $SCRIPT_DIR/TestvCircTrappist/aln_verifa.txt
rm $SCRIPT_DIR/TestvCircTrappist/aln_split.sam
rm $SCRIPT_DIR/TestvCircTrappist/aln_circ_list.sam
