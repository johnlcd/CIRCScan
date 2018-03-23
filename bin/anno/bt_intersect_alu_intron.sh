#!/bin/bash

echo ======================================

echo " "
echo USAGE:
echo " "
echo " $0 -a < a.bed > -b < b.bed> -o < out_put > < -s / -S > "
echo " "

echo --------------------------------------

echo " "
echo EXAMPLEs:
echo " "
echo " $0 -a aaa.bed -b bbb.bed -o abc -s "
echo " "

echo ======================================

echo " "
echo Command is: $0 $*
echo " "


if [ $# == "7" ] && [ $1 == "-a" ] && [ $3 == "-b" ] && [ $5 == "-o" ]; then

	
	if [ $7 == "-s" ]; then
	
		bedtools intersect -s -F 1.00 -a $2 -b $4 -wa -wb > $6

	elif [ $7 == "-S" ]; then

		bedtools intersect -S -F 1.00 -a $2 -b $4 -wa -wb > $6

	fi

#		bedtools intersect -s -F 1.00 -a intron.bed -b alu_hg19.bed -wa -wb > fw_alu_intron_overlap

#		bedtools intersect -S -F 1.00 -a intron.bed -b alu_hg19.bed -wa -wb > rev_alu_intron_overlap

fi


## END
