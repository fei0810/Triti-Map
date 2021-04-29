#!/usr/bin/env sh
set -e
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

input1=$(realpath $1)
dir=$(dirname $1)
output1=$(realpath $2)
output2=$(realpath $3)
id=$(basename $1 | sed 's/_snpindex.input.txt//')
depth=$4
minindex=$5
genome=$6
#output3=$7

if [ -f "${genome}.fai" ]; then
	echo genome index exists
else
	samtools faidx $genome
fi

awk 'BEGIN{OFS="\t"}NR>1{print $1,$2,$3,"snp"NR-1,$14}' $input1 >$dir/temp.${id}.bed

bedops --chop 1000000 --stagger 10000 <(awk '{print $1"\t0\t"$2}' ${genome}.fai) >$dir/temp.${id}.1m10k.bed

bedmap --delim '\t' --prec 3 --skip-unmapped --echo --count --mean --echo-map-id $dir/temp.${id}.1m10k.bed $dir/temp.${id}.bed >$dir/temp.${id}.txt

awk -v depth=$depth -v minindex=$minindex '$4>depth && ($5>=minindex || $5<=-minindex)' $dir/temp.${id}.txt | bedops -m - | awk 'BEGIN{OFS="\t"}{print $1,$2,$3+1000000}' | bedops -m - | awk 'BEGIN{OFS="\t"}{print $1,$2,$3-1000000}' >${output1}

bedmap --skip-unmapped --delim '\t' --header --echo $input1 $output1 | cat <(head -n1 $input1) - >${output2}

if [ -d $dir/temp_output ]; then
	rm -rf $dir/temp_output
fi
mkdir $dir/temp_output && mv $dir/temp.${id}* $dir/temp_output
