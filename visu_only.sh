#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
echo "Script directory: $SCRIPT_DIR"

#getting options for the diverse scripts
while getopts ":B:E:I:O:*" option ;
do
	case $option in
	B)
		echo received -B with $OPTARG "as a limit for the graph file"
		beg=$OPTARG
		;;
		:)
			echo "option $OPTARG needs an argument"
		;;
	E)
		echo received -E with $OPTARG "as a limit for the graph file"
		en=$OPTARG
		;;
		:)
			echo "option $OPTARG needs an argument"
		;;			
	I)
		echo received -I with $OPTARG "to determine your input file"
		input=$OPTARG
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

#launching the graphical visualisator
python3 $SCRIPT_DIR/covvisualisator.py -b $beg -e $en -i $input -o $output
python3 $SCRIPT_DIR/circ_visualisator.py -b $beg -e $en -i $input -o $output
