#!/bin/bash  
set -x

#GATK 4 picard 2.20 samtools 1.16.1 bedtools v2.31.0
genome=
prefix=293T_V6

script_dir=./
#RNA_bam file name
threads=8
knownSNV=/public3/home/scg8972/genome/dbsnp150.vcf
boost_workingdir=/public3/home/scg8972/0_czhang/DVR_Analysis/2.5_boosted_scan/293T_V6/
DNA_BQSR_dir=/public3/home/scg8972/0_czhang/DVR_Analysis/4_BQSR_DNA/293T_COM/
RNA_bam_dir=/public3/home/scg8972/0_czhang/DVR_Analysis/2_BAM_RNA-seq/293T/

SampleNamelist=(C6-DMSO-1 \
    C6-DMSO-2 \
    C6-DMSO-3 \
    C6-DOX-1 \
    C6-DOX-2 \
    C6-DOX-3)

#ensure_path_ends_with_slash
ensure_path_ends_with_slash() {  
    local path="$1"  
    # 检查路径是否以/结尾，如果不是，则添加  
    if [[ ! "$path" =~ /$ ]]; then  
        echo "$path/"  
    else  
        echo "$path"  
    fi  
}
DNA_BQSR_dir=$(ensure_path_ends_with_slash "$DNA_BQSR_dir")  
boost_workingdir=$(ensure_path_ends_with_slash "$boost_workingdir")  
RNA_bam_dir=$(ensure_path_ends_with_slash "$RNA_bam_dir")
DNA_bam=/public3/home/scg8972/0_czhang/DVR_Analysis/1_BAM_WGS/293T/C6-WGSA_plasmid_sorted.bam
DNA_name=C6_WGS
RNA_bam_sample1=(${boost_workingdir}C6-DMSO-1_split.bam \
    ${boost_workingdir}C6-DMSO-2_split.bam \
    ${boost_workingdir}C6-DMSO-3_split.bam)
RNA_bam_sample2=(${boost_workingdir}C6-DOX-1_split.bam \
    ${boost_workingdir}C6-DOX-2_split.bam \
    ${boost_workingdir}C6-DOX-3_split.bam)

cd $boost_workingdir
for SamplePrefix in "${SampleNamelist[@]}";
do  
    python ${script_dir}bam_calibration_RNA_boost.py \
        --bam ${RNA_bam_dir}${SamplePrefix}.final.bam \
        --output ${SamplePrefix} \
        --genome $genome
    RNA_report_list+=("-I ${boost_workingdir}_recalibration_report.grp ")
done
#========DNA_BQSR=========
cd $DNA_BQSR_dir
python ${script_dir}bam_calibration_DNA.py \
    --bam $DNA_bam \
    --output $DNA_name \
    --genome $genome \
    --known $knownSNV
#=======================

RNA_bam=(${RNA_bam_sample1[@]} ${RNA_bam_sample2[@]})

for element in "${RNA_bam[@]}";
do  
    RNA_bam_list+=("-I ${element} ")
done
cd $boost_workingdir
gatk Mutect2 -R $genome ${RNA_bam_list[@]} \
    -I ${DNA_BQSR_dir}${DNA_name}_recalibration.bam \
    -normal $DNA_name
    -O $prefix.mutect.vcf

gatk FilterMutectCalls \
    -V $prefix.mutect.vcf \
    -R $genome \
    -O ${prefix}_boost.filter.vcf \
    --max-events-in-region 4

bcftools view -f PASS ${prefix}_boost.filter.vcf > ${prefix}_boost.known.vcf
gatk IndexFeatureFile -I ${prefix}_boost.known.vcf