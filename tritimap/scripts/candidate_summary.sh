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

blast_info=$1
pfam_info=$2
mapping_info=$3
output=$4
id=$(basename $1 | sed 's/_candidate_denovo_blast_anno.txt//')
dir=$(dirname $1)

paste <(awk -F"\t" 'BEGIN{getline;a=$1;printf ("%s\t%s",$1,$4)}{if(a==$1){printf "//"$4}else{printf "\n%s\t%s",$1,$4;a=$1}}END{printf "\n"}' $pfam_info) <(awk -F"\t" 'BEGIN{getline;a=$1;printf ("%s\t%s",$1,$2)}{if(a==$1){printf "//"$2}else{printf "\n%s\t%s",$1,$2;a=$1}}END{printf "\n"}' $pfam_info) | cut -f1,2,4 | sort -k1,1 | grep -v 'seqid' >$dir/temp_${id}_pfam.summary.txt

paste <(awk -F"\t" 'BEGIN{getline;a=$1;printf ("%s\t%s",$1,$3)}{if(a==$1){printf "//"$3}else{printf "\n%s\t%s",$1,$3;a=$1}}END{printf "\n"}' $blast_info) <(awk -F"\t" 'BEGIN{getline;a=$1;printf ("%s\t%s",$1,$4)}{if(a==$1){printf "//"$4}else{printf "\n%s\t%s",$1,$4;a=$1}}END{printf "\n"}' $blast_info) | cut -f1,2,4 | sort -k1,1 | grep -v 'seqid' >$dir/temp_${id}_blast.summary.txt

join -1 1 -2 1 $dir/temp_${id}_pfam.summary.txt $dir/temp_${id}_blast.summary.txt -t$'\t' >$dir/temp_${id}_pfam_blast.summary.txt

join -1 1 -2 1 $dir/temp_${id}_pfam_blast.summary.txt <(sort -k1,1 $mapping_info | grep -v 'seqid') -t$'\t' >$dir/temp_${id}_candidate_denovo_summary_info.txt

cat <(awk 'BEGIN{print "SeqID\tPfamModel\tPfamDescription\tBlastHitID\tBlastHitDescription\tRefChrom\tRefStart\tRefEnd\tRefHitLength\tRefHitScore"}') $dir/temp_${id}_candidate_denovo_summary_info.txt >$output && rm -f $dir/temp_${id}_pfam.summary.txt $dir/temp_${id}_blast.summary.txt $dir/temp_${id}_pfam_blast.summary.txt $dir/temp_${id}_candidate_denovo_summary_info.txt
