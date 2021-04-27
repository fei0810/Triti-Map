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

R1=$(realpath $1)
R2=$(realpath $2)
id=$(basename $1 | sed 's/_fastp_R1.fq.gz//')
type=$3
thread=$4
memory=$5
dir=$6
output=$7

if [ -e ${dir}/07_assembleout/${id} ]; then
	rm -rf ${dir}/07_assembleout/${id}
fi

mkdir -p ${dir}/07_assembleout/${id}

if [ $type == "dna" ]; then
	abyss-pe -C ${dir}/07_assembleout/${id} j=$thread k=90 name=$id in=''$R1' '$R2''
	cp ${dir}/07_assembleout/${id}/${id}-scaffolds.fa ${dir}/07_assembleout/temp_${id}_denovo_scaffolds.fasta
	seqkit replace -p .+ -r ""$id"_{nr}" ${dir}/07_assembleout/temp_${id}_denovo_scaffolds.fasta >$output && rm ${dir}/07_assembleout/temp_${id}_denovo_scaffolds.fasta
elif [ $type == "rna" ]; then
	spades.py --rna -1 $R1 -2 $R2 -m $memory -t $thread -o ${dir}/07_assembleout/${id}
	cp ${dir}/07_assembleout/${id}/transcripts.fasta ${dir}/07_assembleout/temp_${id}_denovo_scaffolds.fasta
	seqkit replace -p 'NODE' -r "$id" ${dir}/07_assembleout/temp_${id}_denovo_scaffolds.fasta >$output && rm ${dir}/07_assembleout/temp_${id}_denovo_scaffolds.fasta
else
	exit
fi
