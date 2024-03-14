#!/bin/bash


#getting options for the diverse scripts
while getopts ":B:E:*" option ;
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
	*)
		echo "invalid option $OPTARG"
		;;
esac
done

#launching the graphical visualisator
python3 ~/circJager/visualisator.py -b $beg -e $en
python3 ~/vCircTrappist/circ_visualisator.py -b $beg -e $en
