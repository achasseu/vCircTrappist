# vCircTrappist

## Introduction	    
This program allows to perform circRNA alignments on viral genomes (mono- or multi- segments). The whole script might take from 15 minutes to 3 days to run completely depending on the amount of data to process. In most cases, the results should be available within an hour.

The python-based scripts were used in Chasseur _et al_, 2022 and introduced in Chasseur _et al_, 2024 (in preparation). If you use this program, you should cite the latter.

## Requirements
The set-up we used was the following :
- 64GB RAM
- 512GB SSD
- AMD® Ryzen 7 3700x
  
to run vCircTrappist with a typical sequencing dataset (20Gb).

## Dependencies
#### Python3 (>3.9.7)
  ```
  sudo apt install python3
  sudo apt install python3-pip
  ```
#### Numpy (>1.26.4)
  `pip install numpy`
#### Pandas (>2.2.1)
  `pip install pandas` 
#### Matplotlib (>v3.8.3)
  `pip install matplotlib`
#### Burrows-Wheeler Aligner (BWA) (>0.7.17)
  `sudo apt install bwa`
#### Samtools (>1.10)
  `sudo apt install samtools`
#### Cigar Python library (>0.1.3)
  `pip install cigar`

## Manual
To run vCircTrappist, just download and extract the "vCircTrappist-main" directory in any folder, at your convenience.

### Test
To test it, you should copy/paste the command `bash [path_to_vCircTrappist_folder]/TestvCircTrappist.sh` in your terminal. Both the expected and obtained results are located in the Test folder of vCircTrappist.

### Usage

To use vCircTrappist ("v" for virus ; "circ" for circular RNA ; "Trappist" for the fact that we trap circRNAs but also because of the beers), you should provide (1) a fastq file (single-end or only one mate of your paired-end reads), (2) a fasta reference file of your viral genome, (3) a GFF file of your viral genome (no need for the full annotations but the program won't run without it).


There are three commands that allow you to use vCircTrappist :


(1) `bash [Path_to_the_vCircTrappist_folder]/vCircTrappist_full.sh -F [path_to_your_fasta_genome_file] -G [path_to_your_GFF_genome_file] -Q [path_to_your_fastq_sequencing_file] -S [F or R according to the strandness of your library] -O [Path_to_your_output_folder]`

This command allows you to run the full version of vCircTrappist. This version allows you to run the full script in only one command. It goes from the alignment to the publication-ready graphical visualisation in a straight-forward manner.
We advise you to clean (trimmomatic, dedupe.sh) your data before running the script to obtain more precise results.


(2) `bash [Path_to_the_vCircTrappist_folder]/vCircTrappist_short.sh -F [path_to_your_fasta_file] -G [path_to_your_GFF_file] -Q [path_to_your_fastq_file] -S [F or R according to the strandness of your library] -O [Path_to_your_output_folder]`

This command basically provides the same output but allows you to run the script without carrying out the entire alignment procedure which takes time.


(3) `bash [Path_to_the_vCircTrappist_folder]/visu_only.sh -B [Start_Coordinates_on_the_X_axis] -E [End_Coordinates_on_the_X_axis] -O [Path_to_your_output_folder]`

This command allows you to run only the graphical visualisation program but on a short segment of the genome. It can be useful in the case of pretty long genomes, such as herpesviridae genomes.

### Output

This command generates different data files that can be used for deep circRNA characterization :
- "sites_sorting.csv" is a table of all circRNA features extracted from the alignment.
- "aln_circ_list.sam" is an alignment file of all the circRNAs mapped on the viral genome. You can visualize these using third-party programs such as IGV or Geneious.
- "aln_bsj_sites.csv" is the list of the backjunctions and back-splicing sites.
- Four alignments files "sense_U2", "antisense_U2", "sense_nU2" and "antisense_nU2" allow you to visualize the circRNAs displaying the features you expect to see. For example, you might want to display the circRNAs that only exhibit a canonical splicing signature.
- "aln_split.sam" is an alignment file of all the reads that display a "splicing" signature (linear and circular).
- Two PNG images are generated : "Circ_Coverage_genome.png" and "Circ_list_genome.png". Both of these are publication-ready figures that either depict the full coverage of circRNA reads or back-splice junctions on the genome.

### Example
If you run the test, you should obtain a perfect example of publication-ready figures. In this test, we chose a mapping on the Avian Leukosis Virus Genome, using data from Qiu _et al_, 2018 (DOI: 10.1371/journal.pone.0204931).
![Example of results obtained using vCircTrappist](./TestvCircTrappist/Expected_Coverage_Results.png)
![Example of results obtained using vCircTrappist](./TestvCircTrappist/Expected_Sashimi_Results.png)

### License
This program is the property of the University of Namur and is available under a GPLv3 license.
