#!/usr/bin/env sh

###############################################################################
#
# Author contact:
# zhaofei920810@gmail.com
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
###############################################################################

input=$(realpath $1)
output1=$(realpath $2)
output2=$(realpath $3)
dir=$(dirname $1)
id=$(basename $1 | sed 's/_genofiltered_gatk.vcf//')
pool1=$4
pool2=$5
genome=$7
depth=$6
output3=$(realpath $8)

gatk --java-options "-Xmx10G -Djava.io.tmpdir=./" \
	SelectVariants \
	-R $genome \
	-V $input \
	--exclude-filtered \
	-select-type SNP \
	-xl-select-type MNP \
	--max-nocall-number 0 \
	-sn $pool1 \
	-sn $pool2 \
	-select 'DP >= '$depth'' \
	-select 'vc.getGenotype("'$pool1'").getDP() >= '$depth' && vc.getGenotype("'$pool2'").getDP() >= '$depth'' \
	-O $dir/temp.${id}.realuse.vcf

gatk --java-options "-Xmx10G -Djava.io.tmpdir=./" \
	SelectVariants \
	-R $genome \
	-V $input \
	--exclude-filtered \
	-select-type INDEL \
	--max-nocall-number 0 \
	-sn $pool1 \
	-sn $pool2 \
	-select 'DP >= '$depth'' \
	-select 'vc.getGenotype("'$pool1'").getDP() >= '$depth' && vc.getGenotype("'$pool2'").getDP() >= '$depth'' \
	-O $dir/temp.${id}.realuse.indel.vcf

gatk --java-options "-Xmx10G -Djava.io.tmpdir=./" \
	VariantsToTable \
	-V $dir/temp.${id}.realuse.vcf \
	-R $genome \
	-F CHROM -F POS -F REF -F ALT -GF AD -GF DP \
	-O $dir/temp.${id}.realuse.gatk.table

gatk --java-options "-Xmx10G -Djava.io.tmpdir=./" \
	VariantsToTable \
	-V $dir/temp.${id}.realuse.indel.vcf \
	-R $genome \
	-F CHROM -F POS -F REF -F ALT -GF AD -GF DP \
	-O $dir/temp.${id}.realuse.indel.gatk.table

if [ -s $dir/temp.${id}.realuse.gatk.table ]; then
	cat $dir/temp.${id}.realuse.gatk.table | awk '$4!~/,/' | tr ',' '\t' | awk -v bulk1=$pool1 -v bulk2=$pool2 'BEGIN{OFS="\t";print "#CHROM\tPOSm1\tPOS\tREF\tALT\t"bulk1"_ref\t"bulk1"_alt\t"bulk1"_depth\t"bulk1"_ratio\t"bulk2"_ref\t"bulk2"_alt\t"bulk2"_depth\t"bulk2"_ratio\tsnpindex"}NR>1{print $1,$2-1,$2,$3,$4,$5,$6,$7,$6/$7,$8,$9,$10,$9/$10,$6/$7-$9/$10}' >$output1
	cat $output1 | awk -v bulk1=$pool1 -v bulk2=$pool2 'BEGIN{OFS="\t"; print "CHROM\tPOS\tAD_REF."bulk1"\tAD_ALT."bulk1"\tAD_REF."bulk2"\tAD_ALT."bulk2}NR>1{print $1,$3,$6,$7,$10,$11}' >$output2
	cat $dir/temp.${id}.realuse.indel.gatk.table | awk '$4!~/,/' | tr ',' '\t' | awk -v bulk1=$pool1 -v bulk2=$pool2 'BEGIN{OFS="\t";print "#CHROM\tPOSm1\tPOS\tREF\tALT\t"bulk1"_ref\t"bulk1"_alt\t"bulk1"_depth\t"bulk1"_ratio\t"bulk2"_ref\t"bulk2"_alt\t"bulk2"_depth\t"bulk2"_ratio\tsnpindex"}NR>1{print $1,$2-1,$2,$3,$4,$5,$6,$7,$6/$7,$8,$9,$10,$9/$10,$6/$7-$9/$10}' >$output3
	if [ -d $dir/temp_output ]; then
		rm -rf $dir/temp_output
	fi
	mkdir $dir/temp_output && mv $dir/temp.${id}* $dir/temp_output
else
	echo "No SNP in file"
	exit
fi
