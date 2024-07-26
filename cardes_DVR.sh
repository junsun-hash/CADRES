#!/bin/bash  
#SBATCH -p amd_512  
#SBATCH -N 1  
#SBATCH -n 1  
#SBATCH -c 64  
#SBATCH -o /public3/home/scg8972/slurm/slurm-%j.out  
set -x
source activate DVR
#GATK 4 picard 2.20 samtools 1.16.1 bedtools v2.31.0
genome=/public3/home/scg8972/genome/plsmid/Homo_sapiens.GRCh38.dna.primary_assembly_plsmid.fa
prefix=293T_V6

#所需脚本所在位置的目录，以"/"结尾
#注意！所有dir都要以"/"结尾！！
script_dir=./
#RNA_bam file name
threads=8
knownSNV=/public3/home/scg8972/genome/dbsnp150_reorder.vcf
genomead=/public3/home/scg8972/genome/germline_resource/af-only-gnomad_chr_3.hg38.vcf
geneanno=/public3/home/scg8972/genome/DVR_ref/refGene.txt
knownediting=/public3/home/scg8972/genome/DVR_ref/REDIportal_AtoI.txt
#========RNA_BQSR=========
#不同步骤的工作目录，boost_workingdir必须和上个脚本目录相同
boost_workingdir=/public3/home/scg8972/0_czhang/DVR_Analysis/2.5_boosted_scan/293T_V6/
RNA_BQSR_dir=/public3/home/scg8972/0_czhang/DVR_Analysis/3_BQSR_RNA/293T_V6/
DNA_BQSR_dir=/public3/home/scg8972/0_czhang/DVR_Analysis/4_BQSR_DNA/293T_COM/
cd $RNA_BQSR_dir
#gatk IndexFeatureFile -I ${boost_workingdir}${prefix}_2.known.vcf
SampleNamelist=(C6-DMSO-1 \
    C6-DMSO-2 \
    C6-DMSO-3 \
    C6-DOX-1 \
    C6-DOX-2 \
    C6-DOX-3)

for SamplePrefix in "${SampleNamelist[@]}";
do  
    gatk BaseRecalibrator -I ${boost_workingdir}${SamplePrefix}_split.bam \
        -R $genome \
        -O ${SamplePrefix}_recalibration_report.grp \
        --known-sites ${boost_workingdir}${prefix}_2.known.vcf
    RNA_report_list+=("-I ${SamplePrefix}_recalibration_report.grp ")
done

RNA_report_list_str=${RNA_report_list[@]}
gatk GatherBQSRReports $RNA_report_list_str -O ${prefix}_recalibration_report.grp

for SamplePrefix in "${SampleNamelist[@]}";
do  
    gatk ApplyBQSR -R $genome \
        -I ${boost_workingdir}${SamplePrefix}_split.bam \
        --bqsr-recal-file ${prefix}_recalibration_report.grp \
        -O ${SamplePrefix}_recalibration.bam    
done
gatk AnalyzeCovariates \
    -bqsr ${prefix}_recalibration_report.grp \
    -plots AnalyzeCovariates.pdf
ehco "========CalculateContamination========="
for SamplePrefix in "${SampleNamelist[@]}";
do  
    gatk --java-options '-Xmx200G' GetPileupSummaries \
        -I ${SamplePrefix}_recalibration.bam \
        -V $genomead \
        -L $genomead \
        -O ${SamplePrefix}_pileups.table
    gatk --java-options '-Xmx200G' CalculateContamination \
        -I ${SamplePrefix}_pileups.table \
        -O ${SamplePrefix}_contamination.table
    contamination_list+=("--contamination-table ${RNA_BQSR_dir}${SamplePrefix}_contamination.table ")    
done
contamination_list_str=${contamination_list[@]}
echo $contamination_list_str

echo "========RV&DVR_SETTING============"
RV_workingdir=/public3/home/scg8972/0_czhang/DVR_Analysis/5_mutect2/293T_V6/
DNA_bam=${DNA_BQSR_dir}C6_WGS_recalibration.bam
RNA_bam_sample1=(${RNA_BQSR_dir}C6-DMSO-1_recalibration.bam \
    ${RNA_BQSR_dir}C6-DMSO-2_recalibration.bam \
    ${RNA_BQSR_dir}C6-DMSO-3_recalibration.bam)
RNA_bam_sample2=(${RNA_BQSR_dir}C6-DOX-1_recalibration.bam \
    ${RNA_BQSR_dir}C6-DOX-2_recalibration.bam \
    ${RNA_BQSR_dir}C6-DOX-3_recalibration.bam)
normal=C6_WGS

#=======RV======
cd $RV_workingdir
RNA_bam=(${RNA_bam_sample1[@]} ${RNA_bam_sample2[@]})
RNA_bam_sample1=$(echo "${RNA_bam_sample1[@]}" | tr ' ' ',')  
RNA_bam_sample2=$(echo "${RNA_bam_sample2[@]}" | tr ' ' ',')  

for element in "${RNA_bam[@]}";
do  
    RNA_bam_list+=("-I ${element} ")
done
RNA_bam_list_str=${RNA_bam_list[@]}

gatk Mutect2 -R $genome $RNA_bam_list_str \
    -I $DNA_bam \
    -normal $normal \
    --germline-resource $genomead \
    -O $prefix.temp.vcf

gatk FilterMutectCalls \
    -V $prefix.temp.vcf \
    -R $genome \
    $contamination_list_str \
    --min-median-base-quality 12 \
    --max-events-in-region 4 \
    -O $prefix.temp2.vcf

bcftools view -f PASS $prefix.temp2.vcf > $prefix.vcf

gatk SelectVariants \
    -R $genome \
    -V $prefix.vcf  \
    --select-type-to-include SNP \
    -O $prefix.SNV.vcf

bash ${script_dir}revised_convertVCF.sh $prefix.SNV.vcf ${prefix}_SNPIR.SNV.vcf

#export PERL5LIB=/public3/home/scg8972/0_czhang/Soft/SNPiR
perl ${script_dir}filter_homopolymer_nucleotides.pl \
    -infile ${prefix}_SNPIR.SNV.vcf \
    -outfile $prefix.homo.vcf \
    -refgenome $genome
#制作用于pblat的all.bam参考文件

samtools merge -f ${prefix}all.bam ${RNA_bam[@]}
samtools sort -@ $threads -o ${prefix}all_sorted.bam ${prefix}all.bam
samtools index ${prefix}all_sorted.bam

#export PERL5LIB=/public3/home/scg8972/0_czhang/Soft/SNPiR

perl ${script_dir}pblat_candidates_ln.pl \
    -infile $prefix.homo.vcf \
    -outfile $prefix.pblat.vcf \
    -bamfile ${prefix}all_sorted.bam \
    -refgenome $genome \
    -minbasequal 5 \
    -threads 32
rm ${prefix}all_sorted.bam
rm ${prefix}all.bam
awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $prefix.pblat.vcf > $prefix.knownedit.vcf

bedtools intersect -a $prefix.SNV.vcf -b $prefix.knownedit.vcf -wa -header > $prefix.final.vcf
#==========DVR===========

for element in "${RNA_bam[@]}";
do  
    samtools view -@ $threads -f 64 -b $element > ${element}_seperate_firstEnd.bam
    samtools view -@ $threads -f 128 -b $element > ${element}_seperate_secondEnd.bam
    RNA_bam_sep_list+=("${element}_seperate_firstEnd.bam ")
    RNA_bam_sep_list+=("${element}_seperate_secondEnd.bam ")
done
RNA_bam_sep_list_str=${RNA_bam_sep_list[@]}
samtools mpileup -B -d 100000 -f $genome -l $prefix.final.vcf -q 30 -Q 17 -a -o ${prefix}.pileup ${RNA_bam_sep_list_str}
rm $RNA_bam_sep_list_str
python ${script_dir}vcf_to_mats_input_For_Mutect2.py \
    $prefix.final.vcf \
    ${prefix}.inc.txt \
    $RNA_bam_sample1 \
    $RNA_bam_sample2 \
    20 5 T ${prefix}.pileup T
python ${script_dir}MATS_LRT.py ${prefix}.inc.txt ${prefix}_rMATS-DVR_results 4 0.0001
python ${script_dir}FDR.py ${prefix}_rMATS-DVR_results_rMATS_Result_P.txt ${prefix}_rMATS-DVR_results_rMATS_Result_FDR.txt

#python '+directory+'snv_annotation.py --input '+output+'_rMATS-DVR_results/rMATS_Result_FDR.txt --output '+output+'_rMATS-DVR_results/rMATS-DVR_Result.txt --summary '+output+'_rMATS-DVR_results/rMATS-DVR_Result_summary.txt --label1 '+lab1+' --label2 '+lab2+' --snp '+known+' --repeat '+repeatmask+' --editing '+knownediting+' --gene +geneanno

python ${script_dir}snv_annotation.py --input ${prefix}_rMATS-DVR_results_rMATS_Result_FDR.txt --output ${prefix}_rMATS-DVR_Result.txt --summary ${prefix}rMATS-DVR_Result_summary.txt --label1 DMSO --label2 DOX --snp $knownSNV --repeat NA --editing $knownediting --gene $geneanno
