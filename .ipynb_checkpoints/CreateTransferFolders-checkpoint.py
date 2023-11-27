# 2023.09.11
# This will replicate the folder structure of each sample and transfer only necessary files


# Import libraries

import os
from pathlib import Path
import glob
import shutil
import re


rnaseq_dir = "/data/heustonef/AML/fastqs/Liu-Bresciani_fastqs_071823/"
fastq_prefix = "CTRL|FPD"


# Get a list of samples to run
sample_list = [f for f in os.listdir(".") if f.startswith(("CTRL", "FPD"))]
for sampleID in sample_list:
    os.chdir(sampleID)    
    [os.rename(f, "{}{}{}".format(sampleID, "_" , f)) for f in os.listdir(".") if f.endswith(("html", "json"))]
    Path("origFastqs").mkdir(parents=True, exist_ok=True)
    [os.rename(f, os.path.join("origFastqs", f)) for f in os.listdir(".") if f.endswith("gz")]
    [os.rename(os.path.join("fastpOUT", f), os.path.join(".", f)) for f in os.listdir("fastpOUT") if f.endswith("gz")]
    if len(os.listdir("fastpOUT")) == 0:
                os.rmdir("fastpOUT")    
    os.chdir(rnaseq_dir)
        

#Pass 1
for sampleID in sample_list: 
    os.chdir(sampleID)
    fastq_files = glob.glob("*gz")
    cd_line = ''.join(("&& cd ", sampleID, " \\"))
    readfilesin_line = ''.join(("--readFilesIn ", ' '.join(fastq_files), " \\"))
    outfilenameprefix_line = ''.join(("--outFileNamePrefix ", sampleID, "_pass1 \\"))
    cmdList = [
    "GENOMEDIR=/fdb/STAR_current/UCSC/hg38/genes-100/ \\"  ,
    "&& SJDBGTFFILE=/fdb/STAR_current/UCSC/hg38/genes.gtf \\",
    cd_line,
    "&& mkdir –p star_pass1 \\",
    "&& STAR \\",
    "--runThreadN $SLURM_CPUS_PER_TASK \\",
    "--genomeDir $GENOMEDIR \\",
    readfilesin_line,
    "--readFilesCommand zcat \\",
    outfilenameprefix_line,
    "--outSAMtype BAM SortedByCoordinate \\",
    "--sjdbGTFfile $SJDBGTFFILE \\",
    "--sjdbOverhang 100 \\",
    "--outFilterMultimapScoreRange 1 \\",
    "--outFilterMultimapNmax 20 \\",
    "--outFilterMismatchNmax 10 \\",
    "--alignIntronMax 500000 \\",
    "--alignMatesGapMax 1000000 \\",
    "--sjdbScore 2 \\",
    "--alignSJDBoverhangMin 1 \\",
    "--outFilterMatchNminOverLread 0.33 \\",
    "--outFilterScoreMinOverLread 0.33 \n\n"]
    swarm_cmd = '\n'.join(cmdList)
    os.chdir(rnaseq_dir)
    with open("STARcmd_pass1.swarm", 'a') as f:
        f.write(swarm_cmd)
    
    
#pass 2
for sampleID in sample_list: 
    os.chdir(sampleID)
    fastq_files = glob.glob("*gz")
    tab_files = os.path.join(rnaseq_dir, "STAR_pass1_SJouttabFiles")
    tab_files = glob.glob('/'.join((tab_files, "*.tab")))
    cd_line = ''.join(("&& cd ", sampleID, " \\"))
    readfilesin_line = ''.join(("--readFilesIn ", ' '.join(fastq_files), " \\"))
    outfilenameprefix_line = ''.join(("--outFileNamePrefix ", sampleID, "_pass2 \\"))
    sjdbFileChrStartEnd_line = ''.join(("--sjdbFileChrStartEnd ", ' '.join(tab_files), " \\"))
    cmdList = [
    "GENOMEDIR=/fdb/STAR_current/UCSC/hg38/genes-100/ \\"  ,
    "&& SJDBGTFFILE=/fdb/STAR_current/UCSC/hg38/genes.gtf \\",
    cd_line,
    "&& mkdir –p star_pass2 \\",
    "&& STAR \\",
    "--runThreadN $SLURM_CPUS_PER_TASK \\",
    "--genomeDir $GENOMEDIR \\",
    readfilesin_line,
    "--readFilesCommand zcat \\",
    outfilenameprefix_line,
    "--outSAMtype BAM SortedByCoordinate \\",
    "--sjdbGTFfile $SJDBGTFFILE \\",
    "--sjdbOverhang 100 \\",
    "--outFilterMultimapScoreRange 1 \\",
    sjdbFileChrStartEnd_line, 
    "--outFilterMultimapNmax 20 \\",
    "--outFilterMismatchNmax 10 \\",
    "--alignIntronMax 500000 \\",
    "--alignMatesGapMax 1000000 \\",
    "--sjdbScore 2 \\",
    "--alignSJDBoverhangMin 1 \n\n"]
    swarm_cmd = '\n'.join(cmdList)
    os.chdir(rnaseq_dir)
    with open("STARcmd_pass2.swarm", 'a') as f:
        f.write(swarm_cmd)



for sampleID in sample_list: 
    os.chdir(sampleID)
    fastq_files = glob.glob("*gz")
    tab_files = os.path.join(rnaseq_dir, "STAR_pass1_SJouttabFiles")
    tab_files = glob.glob('/'.join((tab_files, "*.tab")))
    cd_line = ''.join(("&& cd ", sampleID, " \\"))
    readfilesin_line = ''.join(("--readFilesIn ", ' '.join(fastq_files), " \\"))
    outfilenameprefix_line = ''.join(("--outFileNamePrefix ", sampleID, "_pass2 \\"))
    sjdbFileChrStartEnd_line = ''.join(("--sjdbFileChrStartEnd ", ' '.join(tab_files), " \\"))
    cmdList = [
    "GENOMEDIR=/fdb/STAR_current/UCSC/hg38/genes-100/ \\"  ,
    "&& SJDBGTFFILE=/fdb/STAR_current/UCSC/hg38/genes.gtf \\",
    cd_line,
    "&& mkdir –p star_pass2 \\",
    "&& STAR \\",
    "--runThreadN $SLURM_CPUS_PER_TASK \\",
    "--genomeDir $GENOMEDIR \\",
    readfilesin_line,
    "--readFilesCommand zcat \\",
    outfilenameprefix_line,
    "--outSAMtype BAM SortedByCoordinate \\",
    "--sjdbGTFfile $SJDBGTFFILE \\",
    "--sjdbOverhang 100 \\",
    "--outFilterMultimapScoreRange 1 \\",
    sjdbFileChrStartEnd_line, 
    "--outFilterMultimapNmax 20 \\",
    "--outFilterMismatchNmax 10 \\",
    "--alignIntronMax 500000 \\",
    "--alignMatesGapMax 1000000 \\",
    "--sjdbScore 2 \\",
    "--alignSJDBoverhangMin 1 \n\n"]
      --outSAMattrRGline ID:GRPundef \
      --chimMultimapScoreRange 3 \
      --chimScoreJunctionNonGTAG -4 \
      --chimMultimapNmax 20 \
      --chimNonchimScoreDropMin 10 \
      --peOverlapNbasesMin 12 \
      --peOverlapMMp 0.1 \
      --alignInsertionFlush Right \
      --alignSplicedMateMapLminOverLmate 0 \
      --alignSplicedMateMapLmin 30
    swarm_cmd = '\n'.join(cmdList)
    os.chdir(rnaseq_dir)
    with open("STARcmd_pass2.swarm", 'a') as f:
        f.write(swarm_cmd)



     STAR --genomeDir ${star_index_dir} \                                                                                     
          --readFilesIn ${left_fq_filename} ${right_fq_filename} \                                                                      
          --outReadsUnmapped None \
          --twopassMode Basic \
          --readFilesCommand "gunzip -c" \
          --outSAMstrandField intronMotif \  # include for potential use with StringTie for assembly
          --outSAMunmapped Within 
          --chimSegmentMin 12 \  # ** essential to invoke chimeric read detection & reporting **
          --chimJunctionOverhangMin 8 \
          --chimOutJunctionFormat 1 \   # **essential** includes required metadata in Chimeric.junction.out file.
          --alignSJDBoverhangMin 10 \
          --alignMatesGapMax 100000 \   # avoid readthru fusions within 100k
          --alignIntronMax 100000 \
          --alignSJstitchMismatchNmax 5 -1 5 5 \   # settings improved certain chimera detections
          --outSAMattrRGline ID:GRPundef \
          --chimMultimapScoreRange 3 \
          --chimScoreJunctionNonGTAG -4 \
          --chimMultimapNmax 20 \
          --chimNonchimScoreDropMin 10 \
          --peOverlapNbasesMin 12 \
          --peOverlapMMp 0.1 \
          --alignInsertionFlush Right \
          --alignSplicedMateMapLminOverLmate 0 \
          --alignSplicedMateMapLmin 30




