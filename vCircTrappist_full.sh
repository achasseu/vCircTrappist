#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
echo "Script directory: $SCRIPT_DIR"

#getting options for the diverse scripts
while getopts ":F:G:Q:S:O:*" option ;
do
	case $option in
	F)
		echo received -F with $OPTARG "as a fasta reference file"
		FASTA=$OPTARG
		;;
		:)
			echo "option $OPTARG needs an argument"
		;;
	G)
		echo received -G with $OPTARG "as a gff reference file"
		GFF=$OPTARG
		;;
		:)
			echo "option $OPTARG needs an argument"
		;;
	Q)
		echo received -Q with $OPTARG "as a fastq sequences file"
		FASTQ=$OPTARG
		;;
		:)
			echo "option $OPTARG needs an argument"
		;;
	S)
		echo received -S with $OPTARG "to determine if the library is stranded or not"
		strand=$OPTARG
		;;
		:)
			echo "option $OPTARG needs an argument"
		;;
	O)
		echo received -O with $OPTARG "to determine your output file"
		output=$OPTARG
		;;
		:)
			echo "option $OPTARG needs an argument"
		;;
	*)
		echo "invalid option $OPTARG"
		;;
esac
done

mkdir -p $output

#align with bwa
bwa index $FASTA
bwa mem -a -T15 $FASTA $FASTQ > $output/aln.sam

#deleting non-aligned reads
echo "Deleting non-aligned reads"
samtools view -F4 $output/aln.sam > $output/aln_virusb.sam
python3 $SCRIPT_DIR/Correction_SAM.py -I $output/aln.sam -S $output/aln_virusb.sam -O $output/aln_virus.sam
rm $output/aln.sam
rm $output/aln_virusb.sam

samtools view -bS $output/aln_virus.sam > $output/aln_virus.bam
samtools sort $output/aln_virus.bam > $output/aln_virus_sorted.bam
samtools depth -a $output/aln_virus_sorted.bam > $output/aln_virus_coverage.csv

#looking for splitted reads
echo "Looking for splitted reads"
python3 $SCRIPT_DIR/splitfilter.py -a $output/aln_virus.sam -o $output

#looking for circRNA signatures (back-splicing)
echo "Looking for circRNA signatures (back-splicing)"
python3 $SCRIPT_DIR/circHunter3.py -f $FASTA -o $output

#comparing BS junctions and BS sites (counting)
echo "Comparing BS junctions and BS sites (counting)"
python3 $SCRIPT_DIR/bsj_id.py -o $output

#are the reads the results of polymerase jump ?
echo "Is the splicing happening in a repeated region ?"
python3 $SCRIPT_DIR/repet.py -f $FASTA -o $output

#outputting the reads
echo "Outputting the reads"
python3 $SCRIPT_DIR/circ_listing.py -o $output

#sorting the reads according to their features + outputting their features
echo "Sorting the reads according to their features + outputting their features"
python3 $SCRIPT_DIR/circ_caracterisator.py -f $FASTA -g $GFF -s $strand -o $output
python3 $SCRIPT_DIR/Correction_SAM.py -I $output/aln_virus.sam -S $output/circ_sense_U2_a.sam -O $output/circ_sense_U2.sam
python3 $SCRIPT_DIR/Correction_SAM.py -I $output/aln_virus.sam -S $output/circ_antisense_U2_a.sam -O $output/circ_antisense_U2.sam
python3 $SCRIPT_DIR/Correction_SAM.py -I $output/aln_virus.sam -S $output/circ_sense_nonU2_a.sam -O $output/circ_sense_nonU2.sam
python3 $SCRIPT_DIR/Correction_SAM.py -I $output/aln_virus.sam -S $output/circ_antisense_nonU2_a.sam -O $output/circ_antisense_nonU2.sam

#sorting the reads with SAMtools
echo "Sorting the reads with SAMtools"
samtools view -bS $output/circ_sense_U2.sam > $output/circ_sense_U2.bam
samtools sort $output/circ_sense_U2.bam > $output/circ_sense_U2_sorted.bam
samtools depth -a $output/circ_sense_U2_sorted.bam > $output/circ_sense_U2_coverage.csv

samtools view -bS $output/circ_antisense_U2.sam > $output/circ_antisense_U2.bam
samtools sort $output/circ_antisense_U2.bam > $output/circ_antisense_U2_sorted.bam
samtools depth -a $output/circ_antisense_U2_sorted.bam > $output/circ_antisense_U2_coverage.csv

samtools view -bS $output/circ_sense_nonU2.sam > $output/circ_sense_nonU2.bam
samtools sort $output/circ_sense_nonU2.bam > $output/circ_sense_nonU2_sorted.bam
samtools depth -a $output/circ_sense_nonU2_sorted.bam > $output/circ_sense_nonU2_coverage.csv

samtools view -bS $output/circ_antisense_nonU2.sam > $output/circ_antisense_nonU2.bam
samtools sort $output/circ_antisense_nonU2.bam > $output/circ_antisense_nonU2_sorted.bam
samtools depth -a $output/circ_antisense_nonU2_sorted.bam > $output/circ_antisense_nonU2_coverage.csv

#launching the graphical visualisator
python3 $SCRIPT_DIR/covvisualisator.py -o $output
python3 $SCRIPT_DIR/circ_visualisator.py -o $output

echo "Job done."
