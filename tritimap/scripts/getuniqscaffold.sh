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

input1=$(realpath $1)
id1=$(basename $1 | sed 's/_denovo_scaffolds.fasta//')
input2=$(realpath $2)
id2=$(basename $2 | sed 's/_denovo_scaffolds.fasta//')
dir=$(dirname $1)
region=$(realpath $3)
length=$4
genome=$5
thread=$6
type=$7
mem=$8
output1=$9
output2=${10}
output1_summary=${13}
output2_summary=${14}
unmap1=${11}
unmap2=${12}
#done add to config file
blast_ident=${15}
blast_pre=${15}

# # get longer than 500bp scaffold
bioawk -v len=$length -c fastx '{if(length($seq)>=len){print ">"$name;print $seq}}' $input1 >${dir}/temp.${id1}.scaffolds.500.fasta
bioawk -v len=$length -c fastx '{if(length($seq)>=len){print ">"$name;print $seq}}' $input2 >${dir}/temp.${id2}.scaffolds.500.fasta

# map to ref
if [ $type == "dna" ]; then
	bwa-mem2 mem -v 1 -t $thread $genome ${dir}/temp.${id1}.scaffolds.500.fasta | samtools view -S -b - | samtools sort -@ $thread -o ${dir}/temp.${id1}.scaffolds.500.bam -
	bwa-mem2 mem -v 1 -t $thread $genome ${dir}/temp.${id2}.scaffolds.500.fasta | samtools view -S -b - | samtools sort -@ $thread -o ${dir}/temp.${id2}.scaffolds.500.bam -
elif [ $type == "rna" ]; then
	minimap2 -ax splice --split-prefix minimap_${id1} -I $mem $genome ${dir}/temp.${id1}.scaffolds.500.fasta -t 30 | samtools view -S -b - | samtools sort -o ${dir}/temp.${id1}.scaffolds.500.bam -
	minimap2 -ax splice --split-prefix minimap_${id2} -I $mem $genome ${dir}/temp.${id2}.scaffolds.500.fasta -t 30 | samtools view -S -b - | samtools sort -o ${dir}/temp.${id2}.scaffolds.500.bam -
else
	exit
fi

# get not completely matched scaffold
samtools view ${dir}/temp.${id1}.scaffolds.500.bam -F 2048 | egrep 'SA:Z|XA:Z' | awk -v name=${id1} -F "\t" '$6~/S/ {print ">"$1":"$3":"$4":"$6":"length($10)":"name"\n"$10}' >${dir}/temp.${id1}_step1_candidate.fasta
samtools view ${dir}/temp.${id2}.scaffolds.500.bam -F 2048 | egrep 'SA:Z|XA:Z' | awk -v name=${id2} -F "\t" '$6~/S/ {print ">"$1":"$3":"$4":"$6":"length($10)":"name"\n"$10}' >${dir}/temp.${id2}_step1_candidate.fasta

# get uniq to each bulk
makeblastdb -in ${dir}/temp.${id1}_step1_candidate.fasta -dbtype nucl
makeblastdb -in ${dir}/temp.${id2}_step1_candidate.fasta -dbtype nucl

blastn -db ${dir}/temp.${id1}_step1_candidate.fasta -query ${dir}/temp.${id2}_step1_candidate.fasta -out ${dir}/temp.${id2}2${id1}_step1_candidate.outfmt6 -outfmt 6 -num_threads $thread
blastn -db ${dir}/temp.${id2}_step1_candidate.fasta -query ${dir}/temp.${id1}_step1_candidate.fasta -out ${dir}/temp.${id1}2${id2}_step1_candidate.outfmt6 -outfmt 6 -num_threads $thread

# strictly filter use || ; loosely filter use &&
paste ${dir}/temp.${id1}2${id2}_step1_candidate.outfmt6 <(cat ${dir}/temp.${id1}2${id2}_step1_candidate.outfmt6 | cut -f1 | tr ':' '\t' | awk '{print $(NF-1)}') | awk '$3=$3/100 {print $0"\t"$4/$NF}' | awk -v ident=$blast_ident -v pre=$blast_pre '$3>ident || $NF > pre' | cut -f1 | sort -u >${dir}/temp.${id1}_needfilter.id1

paste ${dir}/temp.${id2}2${id1}_step1_candidate.outfmt6 <(cat ${dir}/temp.${id2}2${id1}_step1_candidate.outfmt6 | cut -f2 | tr ':' '\t' | awk '{print $(NF-1)}') | awk '$3=$3/100 {print $0"\t"$4/$NF}' | awk -v ident=$blast_ident -v pre=$blast_pre '$3>ident || $NF > pre' | cut -f2 | sort -u >${dir}/temp.${id1}_needfilter.id2

cat ${dir}/temp.${id1}_needfilter.id1 ${dir}/temp.${id1}_needfilter.id2 | sort -u >${dir}/temp.${id1}_needfilter.id

paste ${dir}/temp.${id1}2${id2}_step1_candidate.outfmt6 <(cat ${dir}/temp.${id1}2${id2}_step1_candidate.outfmt6 | cut -f2 | tr ':' '\t' | awk '{print $(NF-1)}') | awk '$3=$3/100 {print $0"\t"$4/$NF}' | awk -v ident=$blast_ident -v pre=$blast_pre '$3>ident || $NF > pre' | cut -f2 | sort -u >${dir}/temp.${id2}_needfilter.id1

paste ${dir}/temp.${id2}2${id1}_step1_candidate.outfmt6 <(cat ${dir}/temp.${id2}2${id1}_step1_candidate.outfmt6 | cut -f1 | tr ':' '\t' | awk '{print $(NF-1)}') | awk '$3=$3/100 {print $0"\t"$4/$NF}' | awk -v ident=$blast_ident -v pre=$blast_pre '$3>ident || $NF > pre' | cut -f1 | sort -u >${dir}/temp.${id2}_needfilter.id2

cat ${dir}/temp.${id2}_needfilter.id1 ${dir}/temp.${id2}_needfilter.id2 | sort -u >${dir}/temp.${id2}_needfilter.id

seqkit grep -n -v -f ${dir}/temp.${id1}_needfilter.id ${dir}/temp.${id1}_step1_candidate.fasta >${dir}/temp.${id1}_step2_candidate.fasta
seqkit grep -n -v -f ${dir}/temp.${id2}_needfilter.id ${dir}/temp.${id2}_step1_candidate.fasta >${dir}/temp.${id2}_step2_candidate.fasta

cp ${dir}/temp.${id1}_step2_candidate.fasta ${dir}/${id1}_all_denovo.fasta
cp ${dir}/temp.${id2}_step2_candidate.fasta ${dir}/${id2}_all_denovo.fasta

# get candidate region fasta
cat $region | grep -v 'start' | cut -f1,3,4 | awk '{print $1":"$2"-"$3}' >${dir}/temp.${id1}2${id2}_samtools.region.txt

samtools faidx -r ${dir}/temp.${id1}2${id2}_samtools.region.txt $genome >${dir}/temp.${id1}2${id2}_samtools.region.fasta

makeblastdb -in ${dir}/temp.${id1}2${id2}_samtools.region.fasta -dbtype nucl

# map to database
blastn -db ${dir}/temp.${id1}2${id2}_samtools.region.fasta -query ${dir}/temp.${id1}_step2_candidate.fasta -out ${dir}/temp.${id1}_step3_candidate.outfmt6 -outfmt 6 -num_threads $thread
cut -f1 ${dir}/temp.${id1}_step3_candidate.outfmt6 | sort -u >${dir}/temp.${id1}_step3_candidate.id
seqkit grep -n -f ${dir}/temp.${id1}_step3_candidate.id ${dir}/temp.${id1}_step2_candidate.fasta >${dir}/temp.${id1}_step3_candidate.fasta

cp ${dir}/temp.${id1}_step3_candidate.fasta $output1

cut -f1,2,3,4,9,10 ${dir}/temp.${id1}_step3_candidate.outfmt6 | sed 's/:/\t/6;s/-/\t/' | awk 'BEGIN{OFS="\t"}{print $2,$3+$7,$3+$8,$6,$5,$1}' | uniq -f 5 | awk '{print $NF"\t"$0}' | cut -f1-6 >${dir}/temp.${id1}_candidate_denovo2ref.info.txt

cat <(awk 'BEGIN{print "seqid\tchrom\tstart\tend\thit_length\thit_score"}') ${dir}/temp.${id1}_candidate_denovo2ref.info.txt >$output1_summary

blastn -db ${dir}/temp.${id1}2${id2}_samtools.region.fasta -query ${dir}/temp.${id2}_step2_candidate.fasta -out ${dir}/temp.${id2}_step3_candidate.outfmt6 -outfmt 6 -num_threads $thread
cut -f1 ${dir}/temp.${id2}_step3_candidate.outfmt6 | sort -u >${dir}/temp.${id2}_step3_candidate.id
seqkit grep -n -f ${dir}/temp.${id2}_step3_candidate.id ${dir}/temp.${id2}_step2_candidate.fasta >${dir}/temp.${id2}_step3_candidate.fasta

cp ${dir}/temp.${id2}_step3_candidate.fasta $output2

cut -f1,2,3,4,9,10 ${dir}/temp.${id2}_step3_candidate.outfmt6 | sed 's/:/\t/6;s/-/\t/' | awk 'BEGIN{OFS="\t"}{print $2,$3+$7,$3+$8,$6,$5,$1}' | uniq -f 5 | awk '{print $NF"\t"$0}' | cut -f1-6 >${dir}/temp.${id2}_candidate_denovo2ref.info.txt

cat <(awk 'BEGIN{print "seqid\tchrom\tstart\tend\thit_length\thit_score"}') ${dir}/temp.${id2}_candidate_denovo2ref.info.txt >$output2_summary

#### add feizhao treat unmap scaffolds ####

samtools view -f 4 ${dir}/temp.${id1}.scaffolds.500.bam | awk -v name=${id1} -F "\t" '{print ">"$1":"length($10)":"name"\n"$10}' >${dir}/temp.${id1}_unmap.fasta
samtools view -f 4 ${dir}/temp.${id2}.scaffolds.500.bam | awk -v name=${id2} -F "\t" '{print ">"$1":"length($10)":"name"\n"$10}' >${dir}/temp.${id2}_unmap.fasta

makeblastdb -in ${dir}/temp.${id1}_unmap.fasta -dbtype nucl
makeblastdb -in ${dir}/temp.${id2}_unmap.fasta -dbtype nucl

blastn -db ${dir}/temp.${id1}_unmap.fasta -query ${dir}/temp.${id2}_unmap.fasta -out ${dir}/temp.${id2}2${id1}_unmap.outfmt6 -outfmt 6 -num_threads $thread
blastn -db ${dir}/temp.${id2}_unmap.fasta -query ${dir}/temp.${id1}_unmap.fasta -out ${dir}/temp.${id1}2${id2}_unmap.outfmt6 -outfmt 6 -num_threads $thread

paste ${dir}/temp.${id1}2${id2}_unmap.outfmt6 <(cat ${dir}/temp.${id1}2${id2}_unmap.outfmt6 | cut -f1 | tr ':' '\t' | awk '{print $(NF-1)}') | awk '$3=$3/100 {print $0"\t"$4/$NF}' | awk -v ident=$blast_ident -v pre=$blast_pre '$3>ident || $NF > pre' | cut -f1 | sort -u >${dir}/temp.${id1}_unmap_needfilter.id1

paste ${dir}/temp.${id2}2${id1}_unmap.outfmt6 <(cat ${dir}/temp.${id2}2${id1}_unmap.outfmt6 | cut -f2 | tr ':' '\t' | awk '{print $(NF-1)}') | awk '$3=$3/100 {print $0"\t"$4/$NF}' | awk -v ident=$blast_ident -v pre=$blast_pre '$3>ident || $NF > pre' | cut -f2 | sort -u >${dir}/temp.${id1}_unmap_needfilter.id2

cat ${dir}/temp.${id1}_unmap_needfilter.id1 ${dir}/temp.${id1}_unmap_needfilter.id2 | sort -u >${dir}/temp.${id1}_unmap_needfilter.id

paste ${dir}/temp.${id1}2${id2}_unmap.outfmt6 <(cat ${dir}/temp.${id1}2${id2}_unmap.outfmt6 | cut -f2 | tr ':' '\t' | awk '{print $(NF-1)}') | awk '$3=$3/100 {print $0"\t"$4/$NF}' | awk -v ident=$blast_ident -v pre=$blast_pre '$3>ident || $NF > pre' | cut -f2 | sort -u >${dir}/temp.${id2}_unmap_needfilter.id1

paste ${dir}/temp.${id2}2${id1}_unmap.outfmt6 <(cat ${dir}/temp.${id2}2${id1}_unmap.outfmt6 | cut -f1 | tr ':' '\t' | awk '{print $(NF-1)}') | awk '$3=$3/100 {print $0"\t"$4/$NF}' | awk -v ident=$blast_ident -v pre=$blast_pre '$3>ident || $NF > pre' | cut -f1 | sort -u >${dir}/temp.${id2}_unmap_needfilter.id2

cat ${dir}/temp.${id2}_unmap_needfilter.id1 ${dir}/temp.${id2}_unmap_needfilter.id2 | sort -u >${dir}/temp.${id2}_unmap_needfilter.id

seqkit grep -n -v -f ${dir}/temp.${id1}_unmap_needfilter.id ${dir}/temp.${id1}_unmap.fasta >${dir}/temp.${id1}_step2_unmap.fasta
seqkit grep -n -v -f ${dir}/temp.${id2}_unmap_needfilter.id ${dir}/temp.${id2}_unmap.fasta >${dir}/temp.${id2}_step2_unmap.fasta

cp ${dir}/temp.${id1}_step2_unmap.fasta $unmap1
cp ${dir}/temp.${id2}_step2_unmap.fasta $unmap2

if [ -d $dir/temp_output ]; then
	rm -rf $dir/temp_output
fi
mkdir ${dir}/temp_output && mv ${dir}/temp.${id1}* ${dir}/temp.${id2}* ${dir}/temp_output
