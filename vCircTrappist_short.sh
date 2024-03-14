#!/bin/bash


#getting options for the diverse scripts
while getopts ":F:G:Q:S:*" option ;
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
	*)
		echo "invalid option $OPTARG"
		;;
esac
done



samtools view -bS ./aln_virus.sam > ./aln_virus.bam
samtools sort ./aln_virus.bam > ./aln_virus_sorted.bam
samtools depth -a ./aln_virus_sorted.bam > ./aln_virus_coverage.csv

#looking for splitted reads
echo "Looking for splitted reads"
python3 ~/vCircTrappist/splitfilter.py -a ./aln_virus.sam

#looking for circRNA signatures (back-splicing)
echo "Looking for circRNA signatures (back-splicing)"
python3 ~/vCircTrappist/circHunter3.py -f $FASTA

#comparing BS junctions and BS sites (counting)
echo "Comparing BS junctions and BS sites (counting)"
python3 ~/vCircTrappist/bsj_id.py

#are the reads the results of polymerase jump ?
echo "Are the reads the results of polymerase jump ?"
python3 ~/vCircTrappist/repet.py -f $FASTA

#outputting the reads
echo "Outputting the reads"
python3 ~/vCircTrappist/circ_listing.py

#sorting the reads according to their features + outputting their features
echo "Sorting the reads according to their features + outputting their features"
python3 ~/vCircTrappist/circ_caracterisator.py -f $FASTA -g $GFF -s $strand
python3 ~/vCircTrappist/Correction_SAM.py -I ./aln_virus.sam -S ./circ_sense_U2_a.sam -O ./circ_sense_U2.sam
python3 ~/vCircTrappist/Correction_SAM.py -I ./aln_virus.sam -S ./circ_antisense_U2_a.sam -O ./circ_antisense_U2.sam
python3 ~/vCircTrappist/Correction_SAM.py -I ./aln_virus.sam -S ./circ_sense_nonU2_a.sam -O ./circ_sense_nonU2.sam
python3 ~/vCircTrappist/Correction_SAM.py -I ./aln_virus.sam -S ./circ_antisense_nonU2_a.sam -O ./circ_antisense_nonU2.sam

#sorting the reads with SAMtools
echo "Sorting the reads with SAMtools"
samtools view -bS ./circ_sense_U2.sam > ./circ_sense_U2.bam
samtools sort ./circ_sense_U2.bam > ./circ_sense_U2_sorted.bam
samtools depth -a ./circ_sense_U2_sorted.bam > ./circ_sense_U2_coverage.csv

samtools view -bS ./circ_antisense_U2.sam > ./circ_antisense_U2.bam
samtools sort ./circ_antisense_U2.bam > ./circ_antisense_U2_sorted.bam
samtools depth -a ./circ_antisense_U2_sorted.bam > ./circ_antisense_U2_coverage.csv

samtools view -bS ./circ_sense_nonU2.sam > ./circ_sense_nonU2.bam
samtools sort ./circ_sense_nonU2.bam > ./circ_sense_nonU2_sorted.bam
samtools depth -a ./circ_sense_nonU2_sorted.bam > ./circ_sense_nonU2_coverage.csv

samtools view -bS ./circ_antisense_nonU2.sam > ./circ_antisense_nonU2.bam
samtools sort ./circ_antisense_nonU2.bam > ./circ_antisense_nonU2_sorted.bam
samtools depth -a ./circ_antisense_nonU2_sorted.bam > ./circ_antisense_nonU2_coverage.csv

#launching the graphical visualisator
python3 ~/vCircTrappist/covvisualisator.py
python3 ~/vCircTrappist/circ_visualisator.py
echo "Job done."
