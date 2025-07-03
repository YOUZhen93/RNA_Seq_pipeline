#!/bin/bash

#### RNA QC pipeline
# Author Zhen Y
# require RSeQC

function usage
{
    echo -e "usage: bash RNA_Seq_QC_report_v1.sh
    --conf ${conf} --out_dir ${out_dir}
    --help <showing this message>

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

while [ "$1" != "" ]; do
    case $1 in
        -c | --conf )           shift
                                conf=$1
                                ;;
        -o | --out_dir )        shift
                                out_dir=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done


if [[ ! -d ${out_dir} ]]; then
  echo ""
  echo "output folder: ${out_dir} doesn't exist"
  echo "
  Creating new output folder:
  ${out_dir}

  "
  mkdir -p ${out_dir}
fi


source ${conf}
cd ${out_dir}
out_file=${ID}.RNAseq_QC.tsv
out_file_md=${ID}.RNAseq_QC.markdown
dates=`date`




echo -e "# RNA-Seq QC summary table" > ${out_file_md}
echo -e '
<br>
<br>
<table border="1" style="border:1px solid white; border-collapse: collapse; font-size:15px; width:60%; height=40%" class="table">
  <thead>
    <tr bgcolor="#f0a1a8" style="text-align: left;">
      <th style="padding: 10px">Title</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 5px">Author</td>
      <td>Zhen Y</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Pipeline</td>
      <td>RNA-Seq_pipeline</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Date</td>
      <td>'`date`'</td>
    </tr>
  </tbody>
</table>
<br>
<br>' >> ${out_file_md}


echo -e "RNA-Seq QC summary table" > ${out_file}
echo -e "Author\tZhen Y" >> ${out_file}
echo -e "Pipeline\tdefault RNA-Seq_pipeline.sh" >> ${out_file}
echo -e "Date\t"$( echo `date` ) >> ${out_file}

echo -e "## Basic summary" >> ${out_file_md}
echo -e '
<br>
<br>
<table border="1" style="border:1px solid white; border-collapse: collapse; font-size:15px; width:60%; height=40%" class="table">
  <thead>
    <tr bgcolor="#f0a1a8" style="text-align: left;">
      <th style="padding: 10px">Items</th>
      <th>Stats</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 5px">Sample_ID</td>
      <td>'${ID}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Library_type</td>
      <td>'${Library_type}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Tissue</td>
      <td>'${Tissue}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Condition</td>
      <td>'${Condition}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Replicates</td>
      <td>'${Replicates}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Platform</td>
      <td>'${Platform}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Reference_Genome</td>
      <td>'${Reference_Genome}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Alignment_software</td>
      <td>'${Alignment_software}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">strandness</td>
      <td>'${Strandness}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">annotation</td>
      <td>'${Annotation}'</td>
    </tr>
  </tbody>
</table>
<br>
<br>' >> ${out_file_md}


# echo -e "## Basic summary" >> ${out_file_md}
# echo -e "> > **Sample_ID:**&emsp;"${ID}"<br>" >> ${out_file_md}
# echo -e "> > **Library_type:**&emsp;"${Library_type}"<br>" >> ${out_file_md}
# echo -e "> > **Tissue:**&emsp;"${Tissue}"<br>" >> ${out_file_md}
# echo -e "> > **Condition:**&emsp;"${Condition}"<br>" >> ${out_file_md}
# echo -e "> > **Replicates:**&emsp;"${Replicates}"<br>" >> ${out_file_md}
# echo -e "> > **IP factor:**&emsp;"${IP}"<br>" >> ${out_file_md}
# echo -e "> > **Platform:**&emsp;"${Platform}"<br>" >> ${out_file_md}
# echo -e "> > **Reference_Genome:**&emsp;"${Reference_Genome}"<br>" >> ${out_file_md}
# echo -e "> > **Alignment_software:**&emsp;"${Alignment_software}"<br>" >> ${out_file_md}
# echo -e "> > **Peaks_caller:**&emsp;"${Peaks_caller} >> ${out_file_md}

echo -e "Basic summary" >> ${out_file}
echo -e "Sample_ID\t"${ID} >> ${out_file}
echo -e "Library_type\t"${Library_type} >> ${out_file}
echo -e "Tissue\t"${Tissue} >> ${out_file}
echo -e "Condition\t"${Condition} >> ${out_file}
echo -e "Replicates\t"${Replicates} >> ${out_file}
echo -e "Platform\t"${Platform} >> ${out_file}
echo -e "Reference_Genome\t"${Reference_Genome} >> ${out_file}
echo -e "Alignment_software\t"${Alignment_software} >> ${out_file}
echo -e "strandness\t"${Strandness} >> ${out_file}
echo -e "annotation\t"${Annotation} >> ${out_file}


echo -e "## Reads QC summary" >> ${out_file_md}
echo -e "Reads QC summary" >> ${out_file}

R1totalbases=$( grep "Total written (filtered):" ${trim_galore_report1} | awk '{print $4}' | sed 's/[(|)]//g' )
R2totalbases=$( grep "Total written (filtered):" ${trim_galore_report2} | awk '{print $4}' | sed 's/[(|)]//g' )
TotalbasespostQC=$( expr $(echo ${R1totalbases} | sed 's/,//g') + $(echo ${R2totalbases} | sed 's/,//g') )

echo -e '
<br>
<br>
<table border="1" style="border:1px solid white; border-collapse: collapse; font-size:15px; width:60%; height=40%" class="table">
  <thead>
    <tr bgcolor="#f0a1a8" style="text-align: left;">
      <th style="padding: 10px">Items</th>
      <th>Stats</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 5px">Total reads processed</td>
      <td>'$( grep "Total reads processed:" ${trim_galore_report1} | awk '{print $NF}' )'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Reads with adapters from R1</td>
      <td>'$( grep "Reads with adapters:" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Reads with adapters from R2</td>
      <td>'$( grep "Reads with adapters:" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g')'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Passed Reads from R1</td>
      <td>'$( grep "Reads written (passing filters):" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Passed Reads from R2</td>
      <td>'$( grep "Reads written (passing filters):" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' )'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Total bases processed from R1</td>
      <td>'$( grep "Total basepairs processed:" ${trim_galore_report1} | awk -F" " '{print $4}' )'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Total bases processed from R2</td>
      <td>'$( grep "Total basepairs processed:" ${trim_galore_report2} | awk -F" " '{print $4}' )'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Bases trimmed from R1</td>
      <td>'$( grep "Quality-trimmed:" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Bases trimmed from R2</td>
      <td>'$( grep "Quality-trimmed:" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' )'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Passed bases from R1</td>
      <td>'$( grep "Total written (filtered):" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Passed bases from R2</td>
      <td>'$( grep "Total written (filtered):" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' )'</td>
    </tr>
  </tbody>
</table>
<br>
<br>' >> ${out_file_md}

echo -e "Total reads processed\t"$( grep "Total reads processed:" ${trim_galore_report1} | awk '{print $NF}' ) >> ${out_file}
# echo -e "> > **Total reads processed:**&emsp;"$( grep "Total reads processed:" ${trim_galore_report1} | awk '{print $NF}' )"<br>" >> ${out_file_md}

echo -e "Reads with adapters from R1:\t"$( grep "Reads with adapters:" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file}
# echo -e "> > **Reads with adapters from R1:**&emsp;"$( grep "Reads with adapters:" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )"<br>" >> ${out_file_md}

echo -e "Reads with adapters from R2:\t"$( grep "Reads with adapters:" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g') >> ${out_file}
# echo -e "> > **Reads with adapters from R2:**&emsp;"$( grep "Reads with adapters:" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g')"<br>" >> ${out_file_md}

echo -e "Passed Reads from R1:\t"$( grep "Reads written (passing filters):" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file}
# echo -e "> > **Passed Reads from R1:**&emsp;"$( grep "Reads written (passing filters):" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )"<br>" >> ${out_file_md}

echo -e "Passed Reads from R2:\t"$( grep "Reads written (passing filters):" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file}
# echo -e "> > **Passed Reads from R2:**&emsp;"$( grep "Reads written (passing filters):" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' )"<br>" >> ${out_file_md}

echo -e "Total bases processed from R1:\t"$( grep "Total basepairs processed:" ${trim_galore_report1} | awk -F" " '{print $4}' ) >> ${out_file}
# echo -e "> > **Total bases processed from R1:**&emsp;"$( grep "Total basepairs processed:" ${trim_galore_report1} | awk -F" " '{print $4}' )"<br>" >> ${out_file_md}

echo -e "Total bases processed from R2:\t"$( grep "Total basepairs processed:" ${trim_galore_report2} | awk -F" " '{print $4}' ) >> ${out_file}
# echo -e "> > **Total bases processed from R2:**&emsp;"$( grep "Total basepairs processed:" ${trim_galore_report2} | awk -F" " '{print $4}' )"<br>" >> ${out_file_md}

echo -e "Bases trimmed from R1:\t"$( grep "Quality-trimmed:" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file}
# echo -e "> > **Bases trimmed from R1:**&emsp;"$( grep "Quality-trimmed:" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )"<br>" >> ${out_file_md}

echo -e "Bases trimmed from R2:\t"$( grep "Quality-trimmed:" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file}
# echo -e "> > **Bases trimmed from R2:**&emsp;"$( grep "Quality-trimmed:" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' )"<br>" >> ${out_file_md}

echo -e "Passed bases from R1:\t"$( grep "Total written (filtered):" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file}
# echo -e "> > **Passed bases from R1:**&emsp;"$( grep "Total written (filtered):" ${trim_galore_report1} | awk '{print $NF}' | sed 's/[(|)]//g' )"<br>" >> ${out_file_md}

echo -e "Passed bases from R2:\t"$( grep "Total written (filtered):" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file}
# echo -e "> > **Passed bases from R2:**&emsp;"$( grep "Total written (filtered):" ${trim_galore_report2} | awk '{print $NF}' | sed 's/[(|)]//g' ) >> ${out_file_md}


echo -e "## Alignment summary" >> ${out_file_md}
echo -e "Alignment summary" >> ${out_file}

# create index
bampath=`dirname ${BAM}`
cd ${bampath}
samtools index -@ ${cores} ${BAM}
cd ${out_dir}


# star stats output summary

UMRN=$( grep "Uniquely mapped reads number" ${star_report} | awk -F"|" '{print $2}' | tr -d " \t" | xargs printf "%'.3d" )
UMRR=$( grep "Uniquely mapped reads %" ${star_report} | awk -F"|" '{print $2}' | tr -d " \t" )
NST=$( grep "Number of splices: Total" ${star_report} | awk -F"|" '{print $2}' | tr -d " \t" | xargs printf "%'.3d" )
NSA=$( grep "Number of splices: Annotated (sjdb)" ${star_report} | awk -F"|" '{print $2}' | tr -d " \t" | xargs printf "%'.3d" )
NSGT=$( grep "Number of splices: GT/AG" ${star_report} | awk -F"|" '{print $2}' | tr -d " \t" | xargs printf "%'.3d" )
GTAR=$( echo ${NSGT} / ${NST} | tr -d "," | bc -l | xargs printf "%.3f" )
RMML=$( grep "% of reads mapped to multiple loci" ${star_report} | awk -F"|" '{print $2}' | tr -d " \t" )
CHIMERIC=$( grep "% of chimeric reads" ${star_report} | awk -F"|" '{print $2}' | tr -d " \t" )


echo -e '
<br>
<br>
<table border="1" style="border:1px solid white; border-collapse: collapse; font-size:15px; width:60%; height=40%" class="table">
  <thead>
    <tr bgcolor="#f0a1a8" style="text-align: left;">
      <th style="padding: 10px">Items</th>
      <th>Stats</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 5px">Uniquely mapped reads</td>
      <td>'${UMRN}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Uniquely mapped reads %</td>
      <td>'${UMRR}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Total number of splices</td>
      <td>'${NST}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Annotated number of splices</td>
      <td>'${NSA}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Number of GT/AG splices</td>
      <td>'${NSGT}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">GT/AG splice ratio (over annotated)</td>
      <td>'${GTAR}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">Reads with multiple hits %</td>
      <td>'${RMML}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Chimeric reads %</td>
      <td>'${CHIMERIC}'</td>
    </tr>
  </tbody>
</table>
<br>
<br>' >> ${out_file_md}

echo -e "Uniquely mapped reads\t"${UMRN} >> ${out_file}
# echo -e "> > **Overall alignment rate:**&emsp;"$( grep overal ${bowtie2_report} | awk '{print $1}' )"<br>" >> ${out_file_md}
echo -e "Uniquely mapped reads %\t"${UMRR} >> ${out_file}
# echo -e "> > **Duplicate reads:**&emsp;"${dupreads}"<br>" >> ${out_file_md}
echo -e "Total number of splices\t"${NST} >> ${out_file}
# echo -e "> > **Duplicate reads ratio:**&emsp;"${dupratio}"<br>" >> ${out_file_md}
echo -e "Annotated number of splices\t"${NSA} >> ${out_file}
# echo -e "> > **Paired mapped reads:**&emsp;"${pairedreads}"<br>" >> ${out_file_md}
echo -e "Number of GT/AG splices\t"${NSGT} >> ${out_file}
# echo -e "> > **Singleton:**&emsp;"${singletons}"<br>" >> ${out_file_md}
echo -e "GT/AG splice ratio (over annotated)\t"${GTAR} >> ${out_file}
# echo -e "> > **Insert size:**&emsp;"${insertsize}"<br>" >> ${out_file_md}
echo -e "Reads with multiple hits\t"${RMML} >> ${out_file}
# echo -e "> > **Average MAPQ:**&emsp;"${meanmapq}"<br>" >> ${out_file_md}
echo -e "Chimeric reads %\t"${CHIMERIC} >> ${out_file}
# echo -e "> > **Average coverage:**&emsp;"${meancov}"<br>" >> ${out_file_md}




echo -e "## RNA-Seq QC metrics summary" >> ${out_file_md}
echo -e "RNA-Seq QC metrics summary" >> ${out_file}


# Running RSeQC
echo "Performing RSeQC and calculating TPM ..."

bam_stat.py -i ${BAM} -q 20 > bam_stat
geneBody_coverage.py -i ${BAM} -r ${RNA_HK} -f pdf -o gene_body
split_bam.py -i ${BAM}  -r ${RNAQC_rRNA} -o rRNA 1>rRNA_report && rm -rf rRNA.junk.bam rRNA.ex.bam
junction_annotation.py -i ${BAM} -o junc_anno -r ${RNAQC_BED}
junction_saturation.py -i ${BAM} -r ${RNAQC_BED} -o junc_satura
read_distribution.py  -i ${BAM} -r ${RNAQC_BED} > reads_distribution

#cd ${bampath}
#TPMCalculator -g ${GTF} -c 150 -p -e -b ${BAM}

cd ${out_dir}
RNAR=$( grep "rRNA.in.bam" rRNA_report | awk -F":" '{print $2}' | tr -d " \t" | xargs printf "%'.3d" )
TOTR=$( grep "Total records:" rRNA_report | awk -F":" '{print $2}' | tr -d " \t" | xargs printf "%'.3d" )
RNAra=$( echo ${RNAR} / ${TOTR} | tr -d "," | bc -l | xargs printf "%.3f" )

RDTOT=$( grep "Total Assigned Tags" reads_distribution | sed 's/[^0-9]*//g' | tr -d "\t " )
CDS_Exons=$( echo "$( grep "CDS_Exons" reads_distribution | awk '{print $3}' ) / ${RDTOT}" | bc -l | xargs printf "%.3f" )
U5_Exons=$( echo "$( grep "5'UTR_Exons" reads_distribution | awk '{print $3}' ) / ${RDTOT}" | bc -l | xargs printf "%.3f" )
U3_Exons=$( echo "$( grep "3'UTR_Exons" reads_distribution | awk '{print $3}' ) / ${RDTOT}" | bc -l | xargs printf "%.3f" )
Introns=$( echo "$( grep "Introns" reads_distribution | awk '{print $3}' ) / ${RDTOT}" | bc -l | xargs printf "%.3f" )
TSS_up_1kb=$( echo "$( grep "TSS_up_1kb" reads_distribution | awk '{print $3}' ) / ${RDTOT}" | bc -l | xargs printf "%.3f" )
TSS_up_5kb=$( echo "$( grep "TSS_up_5kb" reads_distribution | awk '{print $3}' ) / ${RDTOT}" | bc -l | xargs printf "%.3f" )
TSS_up_10kb=$( echo "$( grep "TSS_up_10kb" reads_distribution | awk '{print $3}' ) / ${RDTOT}" | bc -l | xargs printf "%.3f" )
TES_down_1kb=$( echo "$( grep "TES_down_1kb" reads_distribution | awk '{print $3}' ) / ${RDTOT}" | bc -l | xargs printf "%.3f" )
TES_down_5kb=$( echo "$( grep "TES_down_5kb" reads_distribution | awk '{print $3}' ) / ${RDTOT}" | bc -l | xargs printf "%.3f" )
TES_down_10kb=$( echo "$( grep "TES_down_10kb" reads_distribution | awk '{print $3}' ) / ${RDTOT}" | bc -l | xargs printf "%.3f" )


echo -e '
<br>
<br>
<table border="1" style="border:1px solid white; border-collapse: collapse; font-size:15px; width:60%; height=40%" class="table">
  <thead>
    <tr bgcolor="#f0a1a8" style="text-align: left;">
      <th style="padding: 10px">Items</th>
      <th>Stats</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 5px">rRNA %</td>
      <td>'${RNAra}'</td>
    </tr>
    <tr>
      <td style="padding: 5px">CDS_Exons %</td>
      <td>'${CDS_Exons}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">5UTR_Exons %</td>
      <td>'${U5_Exons}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">3UTR_Exons %</td>
      <td>'${U3_Exons}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">Introns %</td>
      <td>'${Introns}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">TSS_up_1kb %</td>
      <td>'${TSS_up_1kb}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">TSS_up_5kb %</td>
      <td>'${TSS_up_5kb}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">TSS_up_10kb %</td>
      <td>'${TSS_up_10kb}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">TES_down_1kb %</td>
      <td>'${TES_down_1kb}'</td>
    </tr>
    <tr>
      <td style="padding: 5px;">TES_down_5kb %</td>
      <td>'${TES_down_5kb}'</td>
    </tr>
    <tr style="background-color: rgba(0, 0, 0, 0.2);">
      <td style="padding: 5px;">TES_down_10kb %</td>
      <td>'${TES_down_10kb}'</td>
    </tr>
  </tbody>
</table>
<br>
<br>' >> ${out_file_md}

echo -e "rRNA %\t"${RNAra} >> ${out_file}
# echo -e "> > **Fragment length:**&emsp;"${FragLen}"<br>" >> ${out_file_md}
echo -e "CDS_Exons %\t"${CDS_Exons} >> ${out_file}
# echo -e "> > **FragCC:**&emsp;"${FragCC}"<br>" >> ${out_file_md}
echo -e "5UTR_Exons %\t"${U5_Exons} >> ${out_file}
# echo -e "> > **RelCC:**&emsp;"${RelCC}"<br>" >> ${out_file_md}
echo -e "3UTR_Exons %\t"${U3_Exons} >> ${out_file}
# echo -e "> > **RiP:**&emsp;"${RiP} >> ${out_file_md}
echo -e "Introns %\t"${Introns} >> ${out_file}
echo -e "TSS_up_1kb %\t"${TSS_up_1kb} >> ${out_file}
echo -e "TSS_up_5kb %\t"${TSS_up_5kb} >> ${out_file}
echo -e "TSS_up_10kb %\t"${TSS_up_10kb} >> ${out_file}
echo -e "TES_down_1kb %\t"${TES_down_1kb} >> ${out_file}
echo -e "TES_down_5kb %\t"${TES_down_5kb} >> ${out_file}
echo -e "TES_down_10kb %\t"${TES_down_10kb} >> ${out_file}




echo -e '
### housekeeping genes body coverage
<br>
<embed src="gene_body.geneBodyCoverage.curves.pdf" alt="drawing" width="600" height="700"/>
<br>
<br>' >> ${out_file_md}
echo -e '
### splicing events
<br>
<embed src="junc_anno.splice_events.pdf" alt="drawing" width="600" height="700"/>
<br>
<br>' >> ${out_file_md}
echo -e '
### junction saturation
<br>
<embed src="junc_satura.junctionSaturation_plot.pdf" alt="drawing" width="600" height="700"/>
<br>' >> ${out_file_md}

python3 -m markdown ${out_file_md} > ${out_file_md}.html
