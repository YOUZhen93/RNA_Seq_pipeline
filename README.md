# bulk RNA-Seq analysis pipeline (work on linux)
## Author: Zhen Y

### main script: RNApipeline.sh
example: ChIP-Seq_pipeline.sh -1 1.fastq.gz -2 2.fastq.gz -M input|peak -o outdir -a illumina -r ref_genome_STAR_index -s sample_id -t 10 -P primary_assembly_bed -g gtf_annotation -S hs|mm -m tissue="kidney",sample_id="A1",group="normal"

Notes: -m: miscallenous sample information pass to config file, separated by comma.
       -a: adapters has 4 options: stranded_illumina or illumina or none or BGI; users can provide their own adapter sequences, if two ends have different adapter seqs, use comma to separate them, e.g., ATGC,GTAC
       
       
### config file: RNA-Seq_default.config/RNA-Seq_default.mouse.config (for mouse)
Notes used for QC script ChIP_Seq_QC_report_v1.sh, modified as needed; put under the same folder with main script
this also contain QC required files (mm10 as example):
RNAQC_BED=mm10_GENCODE_vm25.bed  # gencode gene bed file
RNAQC_rRNA=mm10_rRNA.bed  # rRNA gene bed file
RNA_HK=mm10.HouseKeepingGenes.bed  # house keeping gene bed file

### QC script: RNA-Seq_QC_summary.sh
Notes: this would generated .tsv QC report along with .markdown and .html QC report files; keep generated figures with the .html file; require qualimap tool and chipseq_qc.R script;
This requires RSeQC package


