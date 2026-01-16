#!/bin/bash  
set -x
conda activate CADRES
#GATK 4 picard 2.20 samtools 1.16.1 bedtools v2.31.0
genome=/home/sunjun/hdd/sunjun_file/genome/plsmid/Homo_sapiens.GRCh38.dna.primary_assembly_plsmid.fa
prefix=test_293

script_dir=/home/sunjun/cadres/CADRES/
#RNA_bam file name
threads=8
knownSNV=/home/sunjun/hdd/sunjun_file/genome/dbsnp150.vcf
boost_workingdir=/home/sunjun/cadres/3_boost_scan/
DNA_BQSR_dir=/home/sunjun/cadres/1_dna_bam/
RNA_bam_dir=/home/sunjun/cadres/2_rna_bam/

SampleNamelist=(C6-DOX-1A_chr22 \
    C6-DOX-2A_chr22 \
    C6-DOX-3A_chr22 \
    C6+DOX-1A_chr22 \
    C6+DOX-2A_chr22 \
    C6+DOX-3A_chr22)

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
DNA_bam=/home/sunjun/cadres/1_dna_bam/test_sorted_chr22.bam
DNA_name=C6_WGS
RNA_bam_sample1=(${RNA_bam_dir}C6-DOX-1A_chr22.bam \
    ${RNA_bam_dir}C6-DOX-2A_chr22.bam \
    ${RNA_bam_dir}C6-DOX-3A_chr22.bam)
RNA_bam_sample2=(${RNA_bam_dir}C6+DOX-1A_chr22.bam \
    ${RNA_bam_dir}C6+DOX-2A_chr22.bam \
    ${RNA_bam_dir}C6+DOX-3A_chr22.bam)

cd $boost_workingdir
for SamplePrefix in "${SampleNamelist[@]}";
do  
    python ${script_dir}bam_calibration_RNA_boost.py \
        --bam ${RNA_bam_dir}${SamplePrefix}.bam \
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
    -normal $DNA_name \
    -O $prefix.mutect.vcf

gatk FilterMutectCalls \
    -V $prefix.mutect.vcf \
    -R $genome \
    -O ${prefix}_boost.filter.vcf \
    --max-events-in-region 4

bcftools view -f PASS ${prefix}_boost.filter.vcf > ${prefix}_boost.known.vcf
gatk IndexFeatureFile -I ${prefix}_boost.known.vcf