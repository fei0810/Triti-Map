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

#use conda install hmmer
esl-translate $input >${input/fasta/}pep

#use conda install seqkit
seqkit fx2tab -l -n ${input/fasta/}pep | tr ' ' '\t' | cut -f 1,2,7 | sort -k2,2 -k3,3nr | awk '{print $0"\t"$2}' | uniq -f 3 | cut -f1-3 | cut -f1 >${input/fasta/}long.id

echo ${input/fasta/}pep
echo ${input/fasta/}long.id

seqkit grep -f ${input/fasta/}long.id ${input/fasta/}pep >${input/fasta/}long.pep

echo ${input/fasta/}long.pep

seqkit split -s 1 ${input/fasta/}long.pep

for j in $(ls ${input/fasta/}long.pep.split/*.pep); do
	echo $j
	python $scanpath --email $email --database pfam --outformat out,sequence --outfile ${j} --sequence $j
	echo "pfam ${j} done"
	shuf -i 2-5 -n1 | xargs sleep
done

out2tab() {
	file=$1
	num=$(awk '$1~/^$/{print NR}' $file | head -n 1)
	seqid=$(grep '>' ${file/out.txt/}sequence.txt | cut -d ' ' -f2 | sed 's/source=//')
	paste <(head -n $num $file | egrep -v '^$|Scores|---|Description' | sed 's/[ ][ ]* / /g;s/^ //;s/ /\t/9' | awk -F'\t' -v id=$seqid '{print id"\t"$2"\t"$1}' | cut -f1,2) <(head -n $num $file | egrep -v '^$|Scores|---|Description' | sed 's/[ ][ ]* / /g;s/^ //;s/ /\t/9' | awk -F'\t' -v id=$seqid '{print id"\t"$2"\t"$1}' | cut -f3 | tr ' ' '\t' | cut -f4,9)
}

for j in $(ls ${input/fasta/}long.pep.split/*.out.txt); do
	out2tab $j >${j/out.txt/}info.txt
done

cat <(awk 'BEGIN{print "seqid\tDescription\tE-value\tModel"}') ${input/fasta/}long.pep.split/*info.txt >${output1}

seqkit grep -f <(cut -f1 ${output1}) $input >${output2}
