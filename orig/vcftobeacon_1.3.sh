#!/bin/bash
#title           :vcftobeacon.sh
#description     :The following script takes the vcf files from the folder and generates a file with the columns for the Beacon.
#author		 :Dietmar
#date            :16/08/2019
#version         :1.3
#usage		 :bash vcftobeacon.sh
#notes           :for python pandas. Install bcftools and add plugin +fill-tags to the path.
#bash_version    :4.4.20(1)-release




ls -I list *.vcf>list
sed -i 's/.vcf//g' list

echo "The following script takes the vcf files from the folder and generates the three files needed for Beacon."
echo "---------------"
echo "Starting script"

aa=( $(wc <list))
bb=1

for i in $(cat list)
do
	echo "Analyzing" $bb "/" $aa; bb=$((bb+1))
	echo "Normalizing file" $i
	echo "Normalizing file" $i>$i.summary.log
	date>>$i.summary.log
	bcftools norm -m -both $i.vcf -o $i.norm.vcf
	echo "Getting Info from file" $i>>$i.summary.log
	bcftools +fill-tags $i.norm.vcf | bcftools query --allow-undef-tags -H -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%END\t%SVLEN\t%AC\t%AN\t%NS\t%AF\n'> $i.ann1
	echo "Sample" $i " first annotation performed">>$i.summary.log
	
	echo "Generating file XX.variants.csv"
	paste <(bcftools view $i.vcf |\
	awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} 	!/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
	\
	<(bcftools query -f '[\t%SAMPLE=%GT]\n' $i.vcf |awk 'BEGIN {print "nHet"} {print gsub(/0\|1|1\|0|0\/1|1\/0/, "")}') \
	\
	<(bcftools query -f '[\t%SAMPLE=%GT]\n' $i.vcf | awk 'BEGIN {print "nHomAlt"} {print gsub(/1\|1|1\/1/, "")}') \
	 | sed 's/,\t/\t/g' | sed 's/,$//g'>$i.ann2
	echo "Sample" $i " second annotation performed"
	echo "Sample" $i " second annotation performed">>$i.summary.log

	#now we will put in a file the column for the sum of nHet and nHomAlt by running a python script.
	python3 ./sum_var.py -i ./$i.ann2
	echo "Calculating number of Het + Homo_Alt"
	mv Sum_variants $i.Sum_variants

	#getting info for TYPE/SVTYPE:
	bcftools query --allow-undef-tags -H -f'%TYPE\t%SVTYPE\n' $i.norm.vcf > $i.ann3
	echo "Calculating SVTYPE/TYPE"
	python3 ./type_svtype.py -i ./$i.ann3
	mv type $i.type
	
	#As we already have all needed parameters, now we put them all into one file
	paste $i.ann1 $i.Sum_variants $i.type>joined
	sed -i 's/# //' joined
	echo "Generating file" $i".variants.csv"
	awk -v OFS='\t' '{print $1, $2, $3, $4, $5,$6,$13,$7,$8,$9,$10,$11,$12}' joined>$i.variants.csv

	echo "Starting generating file XX.variants.matching.sample.csv"
	paste <(bcftools view $i.vcf |\
	awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
	!/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
	\
	<(bcftools view $i.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HetSamples"} \
	!/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|1|1\|0|0\/1|1\/0/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
	\
	<(bcftools view $i.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesAlt"} \
	!/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/1\|1|1\/1/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
	\
	| sed 's/,\t/\t/g' | sed 's/,$//g'>$i.matching.sample.int
	awk '{print $6, $7}' $i.matching.sample.int> sub.$i.matching.sample.int
	sed -i 's/\t/,/g' sub.$i.matching.sample.int
	sed -i 's/ /,/g' sub.$i.matching.sample.int
	sed -i '0,/HetSamples,HomSamplesAlt/s//sampleId/' sub.$i.matching.sample.int
	awk '{print $1,$2,$3,$4,$5}' $i.matching.sample.int> $i.matching.sample.int2
	paste $i.matching.sample.int2 $i.type sub.$i.matching.sample.int>$i.variants.matching.sample.csv
	
	echo "Starting generating file XX.samples.csv"
	bcftools query -l $i.vcf >$i.samples.csv
	echo $i "DONE"
	date>>$i.summary.log
done


