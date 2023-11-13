#!/bin/sh

#SBATCH --job-name=SRA001_BAM
#SBATCH --output=job%j.out
#SBATCH --error job%j.err
#SBATCH -p gcai-lab
###Number of Cores Max 20
#SBATCH -n 4
###Number of Nodes
#SBATCH -N 1

module load samtools/gcc/1.12
module load java/1.8.0_162
module load bwa/3.16

kim=/home/fqin/SRA001
fastq_dir=$kim/fastq/single
bam_dir=$kim/fastq/single_bam
cd $fastq_dir

total_files=`find -name '*.fastq' | wc -l`
arr=( $(ls *.fastq) )
echo "mapping started" >> map.log
echo "---------------" >> map.log

for ((i=0; i<$total_files; i+=1))
{
sample_name=`echo ${arr[$i]} | awk -F "_" '{print $1}'`
echo "[mapping running for] $sample_name"
echo "bwa mem -t 4 $kim/hg19.fa ${arr[$i]}  > ${sample_name:0:10}.sam" >> map.log
bwa mem -t 4 $kim/hg19.fa ${arr[$i]} > $bam_dir/${sample_name:0:10}.sam

# Convert .sam to .bam
samtools view -bS $bam_dir/${sample_name:0:10}.sam > $bam_dir/${sample_name:0:10}.bam

# Sort
java -Xmx30G -jar $kim/picard.jar SortSam \
    INPUT=$bam_dir/${sample_name:0:10}.bam OUTPUT=$bam_dir/${sample_name:0:10}.sorted.bam \
    SORT_ORDER=coordinate

# Add read group
java -Xmx40G -jar $kim/picard.jar AddOrReplaceReadGroups \
    I=$bam_dir/${sample_name:0:10}.sorted.bam O=$bam_dir/${sample_name:0:10}.sorted.rg.bam RGID=${sample_name:0:10} \
    RGLB=Chung_Et_Al RGPL=ILLUMINA RGPU=machine RGSM=${sample_name:0:10}
samtools index $bam_dir/${sample_name:0:10}.sorted.rg.bam

# Dedup
java -Xmx40G -jar $kim/picard.jar MarkDuplicates \
    REMOVE_DUPLICATES=true \
    I=$bam_dir/${sample_name:0:10}.sorted.rg.bam O=$bam_dir/${sample_name:0:10}.sorted.rg.dedup.bam \
    METRICS_FILE=$bam_dir/${sample_name:0:10}.sorted.rg.dedup.metrics.txt \
    PROGRAM_RECORD_ID= MarkDuplicates PROGRAM_GROUP_VERSION=null \
    PROGRAM_GROUP_NAME=MarkDuplicates
java -jar $kim/picard.jar BuildBamIndex \
    I=$bam_dir/${sample_name:0:10}.sorted.rg.dedup.bam 
}


