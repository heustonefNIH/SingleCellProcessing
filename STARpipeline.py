# 2023.09.11
# This will replicate the folder structure of each sample and transfer only necessary files


# Import libraries

import os
from pathlib import Path
import glob
# import shutil
# import re


rnaseq_dir = "/data/heustonef/AML/fastqs/Liu-Bresciani_fastqs_071823/"
fastq_prefix = "CTRL|FPD"


#"rsem-prepare-reference --gtf $gr39/gencode.v39.primary_assembly.annotation.gtf --num-threads 4 $gr39/GRCh38.primary_assembly.genome.fa rsem-ref/gencode_release39"


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
    outfilenameprefix_line = ''.join(("--outFileNamePrefix ", sampleID, "_pass1- \\"))
    mvSJfile_line = ''.join(("&& mv ", sampleID, "_pass1-SJ.out.tab ../STARpass1_SJouttabFiles"))
    cmdList = [
	"GENOMEDIR=/fdb/STAR_current/GENCODE/Gencode_human/release_39/genes-100 \\"  ,
    cd_line,
    "&& mkdir -p star_pass1 \\",
    "&& STAR \\",
    "--runThreadN $SLURM_CPUS_PER_TASK \\",
    "--genomeDir $GENOMEDIR \\",
    readfilesin_line,
    "--readFilesCommand zcat \\",
    outfilenameprefix_line,
    "--outSAMtype BAM SortedByCoordinate \\",
    "--sjdbOverhang 100 \\",
    "--outFilterMultimapScoreRange 1 \\",
    "--outFilterMultimapNmax 20 \\",
    "--outFilterMismatchNmax 10 \\",
    "--alignIntronMax 500000 \\",
    "--alignMatesGapMax 1000000 \\",
    "--sjdbScore 2 \\",
    "--alignSJDBoverhangMin 1 \\",
    "--outFilterMatchNminOverLread 0.33 \\",
    "--outFilterScoreMinOverLread 0.33 \\",
	mvSJfile_line, 
	"\n\n"]
    swarm_cmd = '\n'.join(cmdList)
    os.chdir(rnaseq_dir)
    with open("STAR_pass1.swarm", 'a') as f:
        f.write(swarm_cmd)
    
    
# #Pass 2 - STAR for sambamba
# for sampleID in sample_list: 
#     os.chdir(sampleID)
#     fastq_files = glob.glob("*gz")
#     tab_files = os.path.join(rnaseq_dir, "STAR_pass1_SJouttabFiles")
#     tab_files = glob.glob('/'.join((tab_files, "*.tab")))
#     cd_line = ''.join(("&& cd ", sampleID, " \\"))
#     readfilesin_line = ''.join(("--readFilesIn ", ' '.join(fastq_files), " \\"))
#     outfilenameprefix_line = ''.join(("--outFileNamePrefix ", sampleID, "_pass2 \\"))
#     sjdbFileChrStartEnd_line = ''.join(("--sjdbFileChrStartEnd ", ' '.join(tab_files), " \\"))
#     cmdList = [
# 	"GENOMEDIR=fdb/GENCODE/Gencode_human/release_39 \\"  ,
#     cd_line,
#     "&& mkdir –p star_pass2 \\",
#     "&& STAR \\",
#     "--runThreadN $SLURM_CPUS_PER_TASK \\",
#     "--genomeDir $GENOMEDIR \\",
#     readfilesin_line,
#     "--readFilesCommand zcat \\",
#     outfilenameprefix_line,
#     "--outSAMtype BAM SortedByCoordinate \\",
#     "--sjdbGTFfile $SJDBGTFFILE \\",
#     "--sjdbOverhang 100 \\",
#     "--outFilterMultimapScoreRange 1 \\",
#     sjdbFileChrStartEnd_line, 
#     "--outFilterMultimapNmax 20 \\",
#     "--outFilterMismatchNmax 10 \\",
#     "--alignIntronMax 500000 \\",
#     "--alignMatesGapMax 1000000 \\",
#     "--sjdbScore 2 \\",
#     "--alignSJDBoverhangMin 1 \n\n"]
#     swarm_cmd = '\n'.join(cmdList)
#     os.chdir(rnaseq_dir)
#     with open("STARcmd_pass2.swarm", 'a') as f:
#         f.write(swarm_cmd)


#Pass 2 - STAR for STAR-FUSION-RSEM branch
#swarm -f STAR_for_fusion_rsem_pass2.swarm --time=24:00:00 -g 100 -t 16 --merge-output --module STAR --sbatch '--mail-type=BEGIN,END,FAIL'
for sampleID in sample_list: 
	os.chdir(sampleID)
	fastq_files = glob.glob("*gz")
	tab_files = os.path.join(rnaseq_dir, "STARpass1_SJouttabFiles")
	tab_files = glob.glob('/'.join((tab_files, "*.tab")))
	cd_line = ''.join(("&& cd ", sampleID, " \\"))
	readfilesin_line = ''.join(("--readFilesIn ", ' '.join(fastq_files), " \\"))
	outfilenameprefix_line = ''.join(("--outFileNamePrefix ", sampleID, "_starfusion \\"))
	sjdbFileChrStartEnd_line = ''.join(("--sjdbFileChrStartEnd ", ' '.join(tab_files), " \\"))
	cmdList = [
	"GENOMEDIR=/fdb/STAR_current/GENCODE/Gencode_human/release_39/genes-100 \\"  ,
	cd_line,
	"&& STAR \\",
	"--runThreadN $SLURM_CPUS_PER_TASK \\",
	"--genomeDir $GENOMEDIR \\",
	readfilesin_line,
	"--readFilesCommand zcat \\",
	outfilenameprefix_line,
	"--outSAMtype BAM SortedByCoordinate \\",
	"--sjdbOverhang 100 \\",
	"--outFilterMultimapScoreRange 1 \\",
	sjdbFileChrStartEnd_line, 
	"--outFilterMultimapNmax 20 \\",
	"--outFilterMismatchNmax 10 \\",
	"--alignIntronMax 500000 \\",
	"--alignMatesGapMax 1000000 \\",
	"--sjdbScore 2 \\",
	"--alignSJDBoverhangMin 1 \\",
	"--outSAMattrRGline ID:GRPundef \\",
	"--outSAMunmapped Within \\",
	"--chimMultimapScoreRange 3 \\",
	"--chimScoreJunctionNonGTAG -4 \\",
	"--chimMultimapNmax 20 \\",
	"--chimNonchimScoreDropMin 10 \\",
	"--peOverlapNbasesMin 12 \\",
	"--peOverlapMMp 0.1 \\",
	"--alignInsertionFlush Right \\",
	"--alignSplicedMateMapLminOverLmate 0 \\",
	"--alignSplicedMateMapLmin 30 \\",
	"--quantMode TranscriptomeSAM \\",
	"--outSAMstrandField intronMotif \\", # include for potential use with StringTie for assembly
	"--chimSegmentMin 12 \\",  # ** essential to invoke chimeric read detection & reporting **
	"--chimJunctionOverhangMin 8 \\",
	"--chimOutJunctionFormat 1 \\",   # **essential** includes required metadata in Chimeric.junction.out file.
	"--alignSJstitchMismatchNmax 5 -1 5 5 \n\n"]   # settings improved certain chimera detections
	star_pass2_for_fusion_rsem = '\n'.join(cmdList)
	os.chdir(rnaseq_dir)
	with open("STAR_for_fusion_rsem_pass2.swarm", 'a') as f:
		f.write(star_pass2_for_fusion_rsem)


#STAR for raw counts
for sampleID in sample_list: 
	os.chdir(sampleID)
	fastq_files = glob.glob("*gz")
	tab_files = os.path.join(rnaseq_dir, "STARpass1_SJouttabFiles")
	tab_files = glob.glob('/'.join((tab_files, "*.tab")))
	cd_line = ''.join(("&& cd ", sampleID, " \\"))
	readfilesin_line = ''.join(("--readFilesIn ", ' '.join(fastq_files), " \\"))
	outfilenameprefix_line = ''.join(("--outFileNamePrefix ", sampleID, "_starfusion \\"))
	sjdbFileChrStartEnd_line = ''.join(("--sjdbFileChrStartEnd ", ' '.join(tab_files), " \\"))
	cmdList = [
	"GENOMEDIR=/fdb/STAR_current/GENCODE/Gencode_human/release_39/genes-100 \\"  ,
	cd_line,
	"&& STAR \\",
	"--runThreadN $SLURM_CPUS_PER_TASK \\",
	"--genomeDir $GENOMEDIR \\",
	readfilesin_line,
	"--readFilesCommand zcat \\",
	outfilenameprefix_line,
	"--outSAMtype BAM SortedByCoordinate \\",
	"--sjdbOverhang 100 \\",
	"--outFilterMultimapScoreRange 1 \\",
	sjdbFileChrStartEnd_line, 
	"--outFilterMultimapNmax 20 \\",
	"--outFilterMismatchNmax 10 \\",
	"--alignIntronMax 500000 \\",
	"--alignMatesGapMax 1000000 \\",
	"--sjdbScore 2 \\",
	"--alignSJDBoverhangMin 1 \\",
	"--outSAMattrRGline ID:GRPundef \\",
	"--outSAMunmapped Within \\",
	"--chimMultimapScoreRange 3 \\",
	"--chimScoreJunctionNonGTAG -4 \\",
	"--chimMultimapNmax 20 \\",
	"--chimNonchimScoreDropMin 10 \\",
	"--peOverlapNbasesMin 12 \\",
	"--peOverlapMMp 0.1 \\",
	"--alignInsertionFlush Right \\",
	"--alignSplicedMateMapLminOverLmate 0 \\",
	"--alignSplicedMateMapLmin 30 \\",
	"--quantMode GeneCounts \n\n"]
	star_for_raw_counts = '\n'.join(cmdList)
	os.chdir(rnaseq_dir)
	with open("star_for_raw_counts.swarm", 'a') as f:
		f.write(star_for_raw_counts)

  
#STAR-FUSION
#swarm -f Fusion.swarm --time=24:00:00 -g 50 -t 16 --merge-output --module STAR,STAR-Fusion --sbatch '--mail-type=BEGIN,END,FAIL'
for sampleID in sample_list: 
	os.chdir(sampleID)
	cd_line = ''.join(("&& cd ", sampleID, " \\"))
	chimericFile = ''.join((sampleID, "_starfusionChimeric.out.junction"))
	J_line = ''.join((("-J ", chimericFile, " \\")))
	cmdList = [
	"FUSIONDIR=/fdb/CTAT/__genome_libs_StarFv1.10/GRCh38_gencode_v37_CTAT_lib_Mar012021/ctat_genome_lib_build_dir \\",
	cd_line,
	"&& STAR-Fusion \\",
	"--genome_lib_dir $FUSIONDIR \\",
	J_line,
	"--FusionInspector validate \\",
	"--denovo_reconstruct \\",
	"--examine_coding_effect \\",
	"--output_dir STARFUSION \n\n"]
	fusion_cmd = '\n'.join(cmdList)
	os.chdir(rnaseq_dir)
	with open("Fusion.swarm", 'a') as f:
		f.write(fusion_cmd)

# 	chimericFile = ''.join((sampleID, "_starfusionChimeric.out.junction"))
# 	cd_line = ''.join(("&& cd ", sampleID, " \\"))
# 	J_line = ''.join((("-J ", chimericFile, " \\")))

# 	cmdList = [
# 	"&& FUSIONDIR=/fdb/CTAT/__genome_libs_StarFv1.10/GRCh38_gencode_v37_CTAT_lib_Mar012021/ctat_genome_lib_build_dir \\"  ,
# 	cd_line,
# 	"&& STAR-Fusion \\",
# 	"--genome_lib_dir $FUSIONDIR \\",
# 	J_line,
# 	"--FusionInspector validate \\",
# 	"--denovo_reconstruct \\",
# 	"--examine_coding_effect \\"]
# 	fusion_cmd = '\n'.join(cmdList)


# gr39="fdb/GENCODE/Gencode_human/release_39"
# rsem-calculate-expression -p $SLURM_CPUS_PER_TASK --alignments --paired-end --output-genome-bam --estimate-rspd --calc-ci --seed 12345 v39pass2Aligned.toTranscriptome.out.bam --forward-prob 0 rsem-ref/gencode_release39 test



#rsem-calculate-expression using pre-loaded indices and STAR-aligned bams
#swarm -f RSEM.swarm --time=24:00:00 -g 50 -t 16 --merge-output --module rsem --sbatch '--mail-type=BEGIN,END,FAIL'

for sampleID in sample_list:
	os.chdir(sampleID)
	cd_line = ''.join(("cd ", sampleID, " \\"))
	alignedToTranscriptome_line = ''.join((sampleID, "_starfusionAligned.toTranscriptome.out.bam \\"))
	rsemref_line = ''.join((rnaseq_dir, "rsem-prepare-ref_GencodeV39/gencode_release39 \\"))
	# sample_info = ''.join((sampleID, ' >& ', sampleID, "rsemLOG.txt"))
	cmdList = [
		cd_line, 
		"&& rsem-calculate-expression \\", 
		"-p $SLURM_CPUS_PER_TASK \\",
		"--alignments \\" ,
		"--paired-end \\",
		"--output-genome-bam \\",
		"--estimate-rspd \\",
		"--calc-ci \\",
		"--seed 12345 \\",
		alignedToTranscriptome_line,
		"--forward-prob 0 \\",
		rsemref_line,
		sampleID,
		"\n\n"]
	rsem_cmd = '\n'.join(cmdList)
	os.chdir(rnaseq_dir)
	with open("RSEM.swarm", 'a') as f:
		f.write(rsem_cmd)




sample_list = [f for f in os.listdir(".") if f.startswith(("CTRL", "FPD"))]
# for sampleID in sample_list: 
# 	os.chdir(sampleID)
# 	fastq_files = glob.glob("*gz")
# 	tab_files = os.path.join(rnaseq_dir, "STARpass1_SJouttabFiles")
# 	tab_files = glob.glob('/'.join((tab_files, "*.tab")))
# 	cd_line = ''.join(("&& cd ", sampleID, " \\"))
# 	readfilesin_line = ''.join(("--readFilesIn ", ' '.join(fastq_files), " \\"))
# 	outfilenameprefix_line = ''.join(("--outFileNamePrefix ", sampleID, "_starfusion \\"))
# 	sjdbFileChrStartEnd_line = ''.join(("--sjdbFileChrStartEnd ", ' '.join(tab_files), " \\"))
# 	cmdList = [
# 	"GENOMEDIR=/fdb/STAR_current/GENCODE/Gencode_human/release_39/genes-100 \\"  ,
# 	cd_line,
# 	"&& mkdir -p STAR-FUSION-RSEM_branch \\",
# 	"&& STAR \\",
# 	"--runThreadN $SLURM_CPUS_PER_TASK \\",
# 	"--genomeDir $GENOMEDIR \\",
# 	readfilesin_line,
# 	"--readFilesCommand zcat \\",
# 	outfilenameprefix_line,
# 	"--outSAMtype BAM SortedByCoordinate \\",
# 	"--sjdbOverhang 100 \\",
# 	"--outFilterMultimapScoreRange 1 \\",
# 	sjdbFileChrStartEnd_line, 
# 	"--outFilterMultimapNmax 20 \\",
# 	"--outFilterMismatchNmax 10 \\",
# 	"--alignIntronMax 500000 \\",
# 	"--alignMatesGapMax 1000000 \\",
# 	"--sjdbScore 2 \\",
# 	"--alignSJDBoverhangMin 1 \\",
# 	"--outSAMattrRGline ID:GRPundef \\",
# 	"--outSAMunmapped Within \\",
# 	"--chimMultimapScoreRange 3 \\",
# 	"--chimScoreJunctionNonGTAG -4 \\",
# 	"--chimMultimapNmax 20 \\",
# 	"--chimNonchimScoreDropMin 10 \\",
# 	"--peOverlapNbasesMin 12 \\",
# 	"--peOverlapMMp 0.1 \\",
# 	"--alignInsertionFlush Right \\",
# 	"--alignSplicedMateMapLminOverLmate 0 \\",
# 	"--alignSplicedMateMapLmin 30 \\",
# 	"--quantMode TranscriptomeSAM \\",
# 	"--outSAMstrandField intronMotif \\", # include for potential use with StringTie for assembly
# 	"--chimSegmentMin 12 \\",  # ** essential to invoke chimeric read detection & reporting **
# 	"--chimJunctionOverhangMin 8 \\",
# 	"--chimOutJunctionFormat 1 \\",   # **essential** includes required metadata in Chimeric.junction.out file.
# 	"--alignSJstitchMismatchNmax 5 -1 5 5 \\"]   # settings improved certain chimera detections
# 	star_pass2_for_fusion_rsem = '\n'.join(cmdList)

# 	chimericFile = ''.join((sampleID, "_starfusionChimeric.out.junction"))
# 	J_line = ''.join((("-J ", chimericFile, " \\")))
# 	cmdList = [
# 	"&& FUSIONDIR=/fdb/CTAT/__genome_libs_StarFv1.10/GRCh38_gencode_v37_CTAT_lib_Mar012021/ctat_genome_lib_build_dir \\"  ,
# 	"&& module load STAR-Fusion \\",
# 	"&& STAR-Fusion \\",
# 	"--genome_lib_dir $FUSIONDIR \\",
# 	J_line,
# 	"--FusionInspector validate \\",
# 	"--denovo_reconstruct \\",
# 	"--examine_coding_effect \\",
# 	"--output_dir STARFUSION \\"]
# 	fusion_cmd = '\n'.join(cmdList)
# 	alignedToTranscriptome_line = ''.join((sampleID, "_starfusionAligned.toTranscriptome.out.bam \\"))
# 	sample_info = ''.join((sampleID, ' >& ', sampleID, "rsemLOG.txt"))
# 	rsemref_line = ''.join((rnaseq_dir, "/rsem-prepare-ref_GencodeV39"))
# 	cmdList = [
# 	"\n&& module purge \\",
# 	"&& module load rsem \\",
# 	"&& rsem-calculate-expression \\", 
# 	"-p $SLURM_CPUS_PER_TASK \\",
# 	"--alignments \\" ,
# 	"--paired-end \\",
# 	"--output-genome-bam \\",
# 	"--estimate-rspd \\",
# 	"--calc-ci \\",
# 	"--seed 12345 \\",
# 	alignedToTranscriptome_line,
# 	"--forward-prob 0 \\",
# 	rsemref_line,
# 	sample_info, 
# 	"\n\n"]
# 	rsem_cmd = '\n'.join(cmdList)
# 	swarm_cmd = ''.join((star_pass2_for_fusion_rsem, fusion_cmd, rsem_cmd))
# 	os.chdir(rnaseq_dir)
# 	with open("STAR-FUSION-RSEM.swarm", 'a') as f:
# 		f.write(swarm_cmd)



