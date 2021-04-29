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

input=$(realpath $1)
email=$2
blastpath=$(realpath $3)
transR=$(realpath $4)
database=$5
output=$6
dir=$(dirname $input)
name=$(basename $input | sed 's/fasta//')

if [ -d ${dir}/temp.blast.${name}split ]; then
	rm -rf ${dir}/temp.blast.${name}split
fi
seqkit split -s 1 --out-dir ${dir}/temp.blast.${name}split ${input}

#use EBI API get json resoults
for j in $(ls ${dir}/temp.blast.${name}split/*.fasta); do
	echo $j
	python $blastpath --email $email --stype dna --program blastn --database $database --outformat out,json --outfile ${j} $j
	echo "blast ${j} done"
	shuf -i 2-5 -n1 | xargs sleep
done

#json 2 data table
for j in $(ls ${dir}/temp.blast.${name}split/*.json.json); do
	echo $j
	Rscript $transR $j ${j/json.json/}table.txt
	head -n 6 ${j/json.json/}table.txt >${j/json.json/}table.top5hit.txt
	echo "json2table $j done"
done

cat <(awk 'BEGIN{print "seqid\thit_db\thit_id\thit_desc\thit_url\thsp_bit_score\thsp_align_len\thsp_identity\thsp_query_from\thsp_query_to\thsp_hit_from\thsp_hit_to"}') <(cat ${dir}/temp.blast.${name}split/*table.top5hit.txt | grep -v 'seqid') >${output}

if [ ! -d $dir/temp_output ]; then
	mkdir $dir/temp_output
fi

rm -rf ${dir}/temp_output/temp.blast.${name}split
mv ${dir}/temp.blast.${name}split ${dir}/temp_output
