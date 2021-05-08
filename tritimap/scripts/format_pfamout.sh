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
email=$2
scanpath=$(realpath $3)
output1=$4
output2=$5
dir=$(dirname $input)
name=$(basename $input | sed 's/fasta//')

#use conda install hmmer
esl-translate $input >${dir}/temp.pfam.${name}pep

#use conda install seqkit
seqkit fx2tab -l -n ${dir}/temp.pfam.${name}pep | tr ' ' '\t' | cut -f 1,2,7 | sort -k2,2 -k3,3nr | awk '{print $0"\t"$2}' | uniq -f 3 | cut -f1-3 | cut -f1 >${dir}/temp.pfam.${name}long.id

echo ${dir}/temp.pfam.${name}pep
echo ${dir}/temp.pfam.${name}long.id

if [ $(awk 'END{print NR}' ${dir}/temp.pfam.${name}long.id) -gt 10000 ]; then
	echo "Your function annotation input sequences is more than 10000, please filter sequences first"
	exit 1
fi

seqkit grep -f ${dir}/temp.pfam.${name}long.id ${dir}/temp.pfam.${name}pep >${dir}/temp.pfam.${name}long.pep

echo ${dir}/temp.pfam.${name}long.pep

seqkit split -s 1 ${dir}/temp.pfam.${name}long.pep

for j in $(ls ${dir}/temp.pfam.${name}long.pep.split/*.pep); do
	if [ ! -f "${j}.out.txt" ]; then
		echo $j
		python $scanpath --email $email --database pfam --outformat out,sequence --outfile ${j} --sequence $j
		while [ $? -ne 0 ]; do
			echo "rerun"
			python $scanpath --email $email --database pfam --outformat out,sequence --outfile ${j} --sequence $j
		done
		echo "pfam ${j} done"
		shuf -i 2-5 -n1 | xargs sleep
	fi
done

out2tab() {
	file=$1
	num=$(awk '$1~/^$/{print NR}' $file | head -n 1)
	seqid=$(grep '>' ${file/out.txt/}sequence.txt | cut -d ' ' -f2 | sed 's/source=//')
	paste <(head -n $num $file | egrep -v '^$|Scores|---|Description' | sed 's/[ ][ ]* / /g;s/^ //;s/ /\t/9' | awk -F'\t' -v id=$seqid '{print id"\t"$2"\t"$1}' | cut -f1,2) <(head -n $num $file | egrep -v '^$|Scores|---|Description' | sed 's/[ ][ ]* / /g;s/^ //;s/ /\t/9' | awk -F'\t' -v id=$seqid '{print id"\t"$2"\t"$1}' | cut -f3 | tr ' ' '\t' | cut -f4,9)
}

for j in $(ls ${dir}/temp.pfam.${name}long.pep.split/*.out.txt); do
	out2tab $j >${j/out.txt/}info.txt
done

cat <(awk 'BEGIN{print "seqid\tDescription\tE-value\tModel"}') ${dir}/temp.pfam.${name}long.pep.split/*info.txt >${output1}

seqkit grep -f <(cut -f1 ${output1}) $input >${output2}

if [ ! -d $dir/temp_output ]; then
	mkdir $dir/temp_output
fi
rm -rf ${dir}/temp_output/temp.pfam.${name}*
mv ${dir}/temp.pfam.${name}* ${dir}/temp_output
