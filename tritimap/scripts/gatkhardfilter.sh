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
output=$(realpath $2)
genome=$3
id=$(basename $1 | sed 's/_raw_gatk.vcf.gz//')
dir=$(dirname $1)
echo $input
echo $id
echo $dir

gatk --java-options "-Xmx10G -Djava.io.tmpdir=./" \
	SelectVariants \
	-R $genome \
	-V $input \
	-select-type SNP \
	--select-type-to-exclude MNP \
	-O $dir/temp_${id}_raw_gatk.snp.vcf

gatk --java-options "-Xmx10G -Djava.io.tmpdir=./" \
	SelectVariants \
	-R $genome \
	-V $input \
	-select-type INDEL \
	-O $dir/temp_${id}_raw_gatk.indel.vcf

gatk --java-options "-Xmx10G -Djava.io.tmpdir=./" \
	VariantFiltration \
	-R $genome \
	-V $dir/temp_${id}_raw_gatk.snp.vcf \
	-window 35 -cluster 3 \
	--filter-expression "(vc.hasAttribute('QD') && QD<2.0 ) || FS > 60.0 || (vc.hasAttribute('MQ') && MQ < 40.0) || ( vc.hasAttribute('ReadPosRankSum' ) && ReadPosRankSum < -8.0 ) || ( vc.hasAttribute('MQRankSum') && MQRankSum < -12.5 ) || ( vc.hasAttribute('SOR') && SOR > 10.0 ) || QUAL < 30.0" \
	--filter-name "filter" \
	-O $dir/temp_${id}_filter_gatk.snp.vcf

gatk --java-options "-Xmx10G -Djava.io.tmpdir=./" \
	VariantFiltration \
	-R $genome \
	-V $dir/temp_${id}_raw_gatk.indel.vcf \
	-window 35 -cluster 3 \
	--filter-expression "(vc.hasAttribute('QD') && QD<2.0) || FS > 200.0 || (vc.hasAttribute ('ReadPosRankSum') && ReadPosRankSum < -20.0) || (vc.hasAttribute('SOR') && SOR > 10.0)|| QUAL < 20.0" \
	--filter-name "filter" \
	-O $dir/temp_${id}_filter_gatk.indel.vcf

gatk --java-options "-Xmx10G -Djava.io.tmpdir=./" \
	MergeVcfs \
	-I $dir/temp_${id}_filter_gatk.snp.vcf \
	-I $dir/temp_${id}_filter_gatk.indel.vcf \
	-O $dir/temp_${id}_filter_gatk.vcf

gatk --java-options "-Xmx10G -Djava.io.tmpdir=./" \
	SelectVariants \
	--exclude-filtered -R $genome \
	-V $dir/temp_${id}_filter_gatk.vcf \
	-O $output

if [ -d $dir/temp_output ]; then
	rm -rf $dir/temp_output
fi
mkdir $dir/temp_output && mv $dir/temp_${id}* $dir/temp_output
