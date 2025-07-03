#!/bin/bash

## RNApipeline.sh

# Gong lab bulk RNA-Seq pipeline (strand specific & PE fastqs) 
# Mon Oct 24 21:47:59 2022 ZY

## This pipeline assume clean fastq inputs and only perform two key steps of RNA-Seq analysis, i.e., alignment and reads counting;
## input 1. fastq files (PE absolute path); 2. destinated folder (should be exist); 3. ref index path (absolute path); 4. sample ID; 5. threads (default 4; optional); 6. gtf file
## 
## For reference genome, you can used pre-indexed refs for hg38 and mm10:
## hg38 pre-built index: /public/home/gonglianggroup/Gonglab/software/hg38/STAR_index
## mm10 pre-built index: /public/home/gonglianggroup/Gonglab/software/mm10/STAR_Index

## Usage:
## RNApipeline.sh -1 fastq1 -2 fastq2 -o outdir -r hg38 -s A1 -t 8 -g hg38.gtf -a ATCAC,GGCAT (-a illumina/stranded_illumina)
## make sure STAR and htseq-count in your env path


set -e

trap 'error_handler $LINENO' ERR
error_handler() {
    local line_num=$1
    echo "Error: pipeline failed at line $line_num"
    exit 1
}

helpFunc()
{
    echo ""
    echo "Usage: ./RNApipeline.sh -1 <input abs path fastq1> -2 <input abs path fastq2> -o <output directory> -r <ref index path> -s <sample ID string> -t <thread default 4> -g <gtf annotation> -a <adapter, separate by comma> -m <miscallenous arguments pass to QC config, separate by comma> -S <specie: mm or hs> -h <showing this message>"
    echo "Heads up: pipeline script, QC summary script, and QC config file should at the same folder"
    echo -e "
   ________     \\        //
          /  /   \\      //
         /  /     \\    //
        /  /       \\  //
       /  /         \\//
      /  /           ||
     /  /            ||
    /  /             ||
   /  /              ||
   ==========        ||
"
 }



thread=8
adapter="illumina"
while getopts "1:2:o:r:s:t:g:a:m:S:h" opt; do
  case $opt in
    1 ) fastq1="$OPTARG"
    ;;
    2 ) fastq2="$OPTARG"
    ;;
    o ) output="$OPTARG"
    ;;
    r ) refindex="$OPTARG"
    ;;
    s ) sample="$OPTARG"
    ;;
    t ) thread="$OPTARG"
    ;;
    g ) gtf="$OPTARG"
    ;;
    a ) adapter="$OPTARG"
    ;;
    m ) misca="$OPTARG"
    ;;
    S ) specie="$OPTARG"
    ;;
    h ) helpFunc ; exit 0
    ;;
    \? )
     echo "Invalid Option" 
     exit 1
    ;;
  esac
done


if [[ $(($# / 2)) -lt 6 ]]; then
    echo ""
    echo "only $(($# / 2)) arguments listed"
    helpFunc
    exit 2
fi

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo ${SCRIPT_DIR}


if [[ ! -d ${output} ]]; then
  echo ""
  echo "output folder: ${output} doesn't exist
  Creating new output folder ..."
  mkdir -p ${output}
fi

module load apps/STAR/2.7.11
module load apps/python/3.10.12

echo ${fastq1}
echo ${fastq2}
echo ${output}
echo ${refindex}
echo ${sample}
echo ${thread}
echo ${gtf}
echo ${misca}
echo ${specie}

timega=`which trim_galore || true`

echo "start trimming fastq ..."
threads_trim=6
mkdir ${output}/clean_reads -p; cd ${output}/clean_reads
if [[ ${adapter} != "none" && ${adapter} == *","* ]]; then
    adapter1=$( echo ${adapter} | awk -F"," '{print $1}' | tr -d " " )
    adapter2=$( echo ${adapter} | awk -F"," '{print $2}' | tr -d " " )
    echo -e "Using adapter: ${adapter1} and adapter2: ${adapter2} ..."
    RUN_TRIM="${timega} -q 20 -j ${threads_trim} --stringency 3 --gzip --paired --length 20 -a ${adapter1} -a2 ${adapter2} --output_dir ${output}/clean_reads --paired ${fastq1} ${fastq2}"

elif [[ ${adapter} =~ ^[A|T|C|G]{3,} && ${adapter} != *","* ]]; then
    echo -e "Using adapter1: ${adapter} ..."
    RUN_TRIM="${timega} -q 20 -j ${threads_trim} --stringency 3 --gzip --paired --length 20 -a ${adapter} --output_dir ${output}/clean_reads --paired ${fastq1} ${fastq2}"

elif [[ ${adapter} != "none" && ${adapter} == "stranded_illumina" ]]; then
    echo -e "Using stranded_illumina default adapters ..."
    RUN_TRIM="${timega} -q 20 -j ${threads_trim} --stringency 3 --gzip --length 20 --paired --stranded_illumina --output_dir ${output}/clean_reads --paired ${fastq1} ${fastq2}"

elif [[ ${adapter} != "none" && ${adapter} == "BGI" ]]; then
    echo -e "Using BGI default adapters ..."
    RUN_TRIM="${timega} -q 20 -j ${threads_trim} --stringency 3 --gzip --paired --length 20 -a AAGTCGGAGGCCAAGCGG -a2 AAGTCGGATCGTAGCCATG --output_dir ${output}/clean_reads --paired ${fastq1} ${fastq2}"
    
elif [[ ${adapter} == "illumina" ]]; then
    echo -e "Using illumina default adapters ..."
    RUN_TRIM="${timega} -q 20 -j ${threads_trim} --stringency 3 --gzip --length 20 --paired --illumina --output_dir ${output}/clean_reads --paired ${fastq1} ${fastq2}"

elif [[ ${adapter} == "none" ]]; then
    echo -e "Using trim_galore searched adapters ..."
    RUN_TRIM="${timega} -q 20 -j ${threads_trim} --stringency 3 --gzip --length 20 --paired --output_dir ${output}/clean_reads --paired ${fastq1} ${fastq2}"

else
echo -e "Error: Please provide adapter sequences or sequencing platform \n"
exit 1
fi

${RUN_TRIM} && echo -e "Trim_galore is done for ${sample}! \n Now changing the input fastq for alignment"

export rawfastq1=`basename ${fastq1}`
export rawfastq2=`basename ${fastq2}`
export fastq1=$( echo `basename ${fastq1}` | sed 's/.fq.gz/_val_1.fq.gz/g' | sed 's/.fastq.gz/_val_1.fq.gz/g' )
export fastq2=$( echo `basename ${fastq2}` | sed 's/.fq.gz/_val_2.fq.gz/g' | sed 's/.fastq.gz/_val_2.fq.gz/g' )
export fastq1=${output}/clean_reads/${fastq1}
export fastq2=${output}/clean_reads/${fastq2}


echo "Start running STAR......"
mkdir ${output}/BAM -p
cd ${output}/BAM
STAR --runThreadN ${thread} --genomeDir ${refindex} \
--readFilesIn ${fastq1} ${fastq2} --readFilesCommand zcat --limitBAMsortRAM 44178937863 --outFileNamePrefix ${sample} \
--outSAMtype BAM SortedByCoordinate --outSAMattrRGline ID:Gong PU:MGI SM:${sample} --outFilterType BySJout --alignMatesGapMax 1000000 \
--outFilterIntronMotifs RemoveNoncanonical --seedSearchStartLmax 50 --alignIntronMin 21 --alignIntronMax 1000000 --alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 --chimSegmentMin 15 --chimSegmentReadGapMax 0 --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic &&

samtools index -@ ${thread} ${sample}Aligned.sortedByCoord.out.bam
echo "STAR alignment completed....."
echo "Start counting....."

Star_counts=${output}/BAM/${sample}ReadsPerGene.out.tab
strandness=(`awk '($1 !~ /^N_/) {sum1+=$2; sum2+=$3; sum3+=$4;} END {print sum2/sum1, sum3/sum1}' ${Star_counts} | xargs printf '%.2f\n'`)
strand_ratio=$( echo ${strandness[0]} / ${strandness[1]} | bc -l )
if(( $( echo "${strand_ratio} > 0.7" | bc -l ) )); then 
    echo "This is a unstranded library, running htseq-count with -s no ..."
    htseqcmd="htseq-count -f bam ${sample}Aligned.sortedByCoord.out.bam ${gtf} -s no -m union"
elif(( $( echo "${strand_ratio} < 0.2" | bc -l ) )); then
    echo "This is a stranded library and primary strand is reverse, running htseq-count with -s reverse ..."
    htseqcmd="htseq-count -f bam ${sample}Aligned.sortedByCoord.out.bam ${gtf} -s reverse -m union"
else
    echo "Not sure about the strandness of this library, better check the results manually, considering as stranded though ..."
    htseqcmd="htseq-count -f bam ${sample}Aligned.sortedByCoord.out.bam ${gtf} -s reverse -m union"
fi


${htseqcmd} > ${sample}_htseq.counts && echo "Counting finished....."



## start doing QC summary
echo "Start RNA-Seq QC summary ... "

cd ${output}
if [ "$specie" == "hs" ]; then
echo "INFO: specie is human ..."
defaul_config=${SCRIPT_DIR}/RNA-Seq_default.config
elif [ "$specie" == "mm" ]; then
echo "INFO: specie is mouse ..."
defaul_config=${SCRIPT_DIR}/RNA-Seq_default.mouse.config
else
echo -e "Error: specie not recognized! only support mm or hs \n"
exit 1
fi

cp ${defaul_config} ${output}
defaul_config=${output}/RNA-Seq_default.config

if [[ ! -n $misca ]]; then
    IFS="," read -r -a new_misca <<< "$misca"
    for item in ${new_misca[@]}; do
        echo ${item} >> ${defaul_config}
    done
fi

echo -e "ID=${sample}" >> ${defaul_config}
echo -e "trim_galore_report1=${output}/clean_reads/${rawfastq1}_trimming_report.txt" >> ${defaul_config}
echo -e "trim_galore_report2=${output}/clean_reads/${rawfastq2}_trimming_report.txt" >> ${defaul_config}
echo -e "BAM=${output}/BAM/${sample}Aligned.sortedByCoord.out.bam" >> ${defaul_config}
echo -e "star_report=${output}/BAM/${sample}Log.final.out" >> ${defaul_config}
echo -e "GTF=${gtf}" >> ${defaul_config}
${SCRIPT_DIR}/RNA-Seq_QC_summary.sh -c ${defaul_config} -o ${output} && echo -e "RNA-Seq QC summary is done"

echo -e "RNA-Seq pipeline completed successfully!"








