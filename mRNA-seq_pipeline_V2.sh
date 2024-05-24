#!/bin/bash
# Processing reads

## Code taken from CB207V course, provided by Anniina Vihervaara, adapted by JT
## These scripts should run in the ag-vm1 in the base environment
## Running this file processes raw read files in fastq.gz format to bedgraph
## density profile
## For each step the file format is given for in- and output files

## You need the CHM13v2chromSizes.txt file and this script in the directory
## with all your replicate_condition_R1.fq and replicate_condition_R2.fq files.
## You also need the chm13v2_hiSat2 files and chm13v2_chrs.gtf.

# The following steps are done sequentially:
# 1 - Trimming reads
# 2 - Decompressing fastq.gz to fq for mapping
# 3 - Aligning reads to genome * this step takes a while ~ 15-30 min/file
# 4 - Sorting the bam files * ~ 5-10 min/file
# 5 - Indexing sorted files ~ 30 s/file
# 6 - Genome coverage bedgraph * > 1h/file
# 7 - bedgraph to bigWig
# 8 - Normalisation (FPM)
# 9 - Feature count
# 10 - Combining reads
# 11 - bedgraph to bigWig

# ------------------------------------------------------------------------------
# Input the replicates, drugs and timepoints used in the experiment here #
# This allows looping through all samples with only a few lines of code.
# You also don't need to change this other than here.

# We loop through all the replicates, drugs and timepoints
Replicates="R1 R2 R3"
Drugs="Ima Bos"
Timepoints="0h 1h 6h 24h Rec" ### These are the 5 timepoints I have

# ------------------------------------------------------------------------------
# 0 - Just a welcome message and setup for the script

printf "\n#### Welcome to the mRNA-seq data processing pipeline ####\n"

sample_start () {
	printf "$1 \t running\t"
}

sample_done () {
	printf "$1 \t done\n"
}

give_start_time () {
	printf "Step $1, started `date "+%a %H:%M:%S"`\n"
}

give_end_time () {
	printf "Step $1, ended `date "+%a %H:%M:%S"` \n \n"
}

# The locale needs to be set to US for the decimal separator to be correct
LC_NUMERIC="en_US.UTF-8"

# ------------------------------------------------------------------------------
# 1 - Trimming reads: Sometimes Illumina reads need to be trimmed
# I will however not be trimming my reads as HiSat2 performs soft-trimming
# .fq.gz

give_start_time "1 - trimming reads"

for replicate in ${Replicates}
do
	for drug in ${Drugs}
	do
		for timepoint in ${Timepoints}
		do
			sample_start ${replicate}${drug}${timepoint}

			fastp -f 10 -t 60 -F 10 -T 60 \
			-i fastq_files/${replicate}${drug}${timepoint}_R1.fq.gz \
			-I fastq_files/${replicate}${drug}${timepoint}_R2.fq.gz \
			-o ${replicate}${drug}${timepoint}_R1_trimmed.fq.gz \
			-O ${replicate}${drug}${timepoint}_R2_trimmed.fq.gz
			# -i inputs the file containing Read 1
			# -I Inputs the file containing Read 2
			# -f nucleotides to clip from the front of Read1
			# -F nucleotides to clip from the Front of Read2
			# -t nucleotides to clip from the tail of Read1
			# -T nucleotides to clip from the Tail of Read2
			# -o outputs a file for trimmed Read 1
			# -O Outputs a file for trimmed Read 2

			sample_done ${replicate}_${timepoint}
	done
done

give_end_time 1

# ------------------------------------------------------------------------------
# 2 - Decompressing fastq files
# Need to decompress files before aligning reads with hiSat2
# Moved this to a separate file: decompress_move.sh
# .fq.gz --> .fq

give_start_time "2 - Decompressing fastq files"

gunzip -rk fastq_files

give_end_time 2

give_start_time "2 - Decompressing fastq files"

reads="R1 R2"

for replicate in ${Replicates}
do
	for drug in ${Drugs}
	do
		for timepoint in ${Timepoints}
		do
			sample_start ${replicate}${drug}${timepoint}

			for read in ${reads}
			do
				gunzip -k fastq_files/${replicate}${drug}${timepoint}_${read}.fastq.gz
			done

			sample_done ${replicate}${drug}${timepoint}
		done
	done
done

give_end_time 2

# ------------------------------------------------------------------------------
# 3 - Aligning reads to genome
# Align reads with hiSat2 loop:
# .fq --> .bam

give_start_time "3 - Aligning reads to genome"

for replicate in ${Replicates}
do
	for drug in ${Drugs}
	do
		for timepoint in ${Timepoints}
		do
			sample_start ${replicate}${drug}${timepoint}

			hisat2 -k 10 --rna-strandness RF --no-mixed --no-discordant -x \
			chm13v2_hiSat2 -1 fastq_files/${replicate}${drug}${timepoint}_R1.fastq -2 \
			fastq_files/${replicate}${drug}${timepoint}_R2.fastq | samtools view -S -b '-' \
			> bams/${replicate}${drug}${timepoint}.bam

			sample_done ${replicate}${drug}${timepoint}
		done
	done
done

give_end_time 3

# # ------------------------------------------------------------------------------
# # 4 - Sorting the bam files
# # Sorts the bam files according to chromosomal positions, necessary for IGV
# # .bam

give_start_time "4 - Sorting the bam files"

for replicate in ${Replicates}
do
	for drug in ${Drugs}
	do
		for timepoint in ${Timepoints}
		do
			sample_start ${replicate}${drug}${timepoint}

			samtools sort bams/${replicate}${drug}${timepoint}.bam -o \
			bams/${replicate}${drug}${timepoint}_sorted.bam

			sample_done ${replicate}${drug}${timepoint}
		done
	done
done

give_end_time 4

# # ------------------------------------------------------------------------------
# # 5 - Indexing sorted bam files
# # Indexes the sorted reads and outputs bai index files, necessary for IGV
# # .bai

give_start_time "5 - Indexing sorted bam files"

for replicate in ${Replicates}
do
	for drug in ${Drugs}
	do
		for timepoint in ${Timepoints}
		do
			sample_start ${replicate}${drug}${timepoint}

			samtools index bams/${replicate}${drug}${timepoint}_sorted.bam

			sample_done ${replicate}${drug}${timepoint}
		done
	done
done

give_end_time 5

# # ------------------------------------------------------------------------------
# ## 6 - Genome coverage bedgraph
# ## generates density profiles with bedtools,
# ## .bam --> .bedgraph
# ## ~ 2 h/file?

give_start_time "6 - Genome coverage bedgraph"

for replicate in ${Replicates}
do
	for drug in ${Drugs}
	do
		for timepoint in ${Timepoints}
		do
			sample_start ${replicate}${drug}${timepoint}

			bedtools genomecov -split -bg -ibam bams/${replicate}${drug}${timepoint}_sorted.bam \
			> bedgraphs/${replicate}${drug}${timepoint}.bedgraph

			sample_done ${replicate}${drug}${timepoint}
		done
	done
done

give_end_time 6

# ------------------------------------------------------------------------------
## 7 - Signal counts and sequencing depth normalisation
## Generates a list of signal counts in each sample.
## .tsv

give_start_time "7 - Signal counts"

for drug in ${Drugs}
do

	#### We start by initiating a text file
	touch ${drug}_0h_1h_6h_24h_rec_mRNAseq_nfs.tsv
	#### Add a header to the file
	printf "sample\ttotal_mapping\n" >> ${drug}_0h_1h_6h_24h_rec_mRNAseq_nfs.tsv

	for timepoint in ${Timepoints}
	do
    for replicate in ${Replicates}
  	do

			sample_start ${replicate}${drug}${timepoint}

			# Count the number of reads in the bam file
			uC=$(samtools view -c bams/${replicate}${drug}${timepoint}_sorted.bam)

			# Divide the reads by 2 (since we have paired-end reads)
			c=`echo $uC / 2 | bc -l`
			echo ${c} #### prints the count on screen

			# Normalises the counts in the bed files to sequencing depth
			echo $c | awk '{c="'$c'"; printf "%s\t%s\t%s\t%s\n", $1, $2, $3, \
			($4*1000000)/c}' bedgraphs/${replicate}${drug}${timepoint}.bedgraph \
			> bedgraphs/${replicate}${drug}${timepoint}_FPM.bedgraph
			#### adjust the fourth column (signal)

			printf ${replicate}${drug}${timepoint}"\t"${c}"\n" >> ${drug}_0h_1h_6h_24h_rec_mRNAseq_nfs.tsv
			#### prints the sample name and its count

			sample_done ${replicate}${drug}${timepoint}
		done
	done
done

give_end_time 7

# ------------------------------------------------------------------------------
## 7.5- Sorting FPM normalised files
## Sorts bedfiles with sort.
## _FPM.bed --> _FPM_sorted.bed

give_start_time "7.5 - Sorting normalised bedgraphs"

for replicate in ${Replicates}
do
  for drug in ${Drugs}
  do
    for timepoint in ${Timepoints}
    do
      sample_start ${replicate}${drug}${timepoint}
      sort -k1,1 -k2,2n bedgraphs/${replicate}${drug}${timepoint}_FPM.bedgraph \
			> bedgraphs/${replicate}${drug}${timepoint}_FPM_sorted.bedgraph
      sample_done ${replicate}${drug}${timepoint}
    done
  done
done

give_end_time 7.5

# # ------------------------------------------------------------------------------
## 8 - bedgraoph to bigWig
## Converts human-readable bedgraphs to human unreadable bigWigs to save space
## .bedgraph --> .bigWig
# Here I changed the script so that it creates bigwigs from non-normalised .bed

give_start_time "8 - bedgraoph to bigWig"

for replicate in ${Replicates}
do
	for drug in ${Drugs}
	do
		for timepoint in ${Timepoints}
		do
			sample_start ${replicate}${drug}${timepoint}
			bedGraphToBigWig bedgraphs/${replicate}${drug}${timepoint}.bedgraph \
			CHM13v2chromSizes.txt bigWigs/${replicate}${drug}${timepoint}.bigWig
			sample_done ${replicate}${drug}${timepoint}
		done
	done
done
give_end_time 8

## 8.5 - bedgraoph to bigWig
## Converts human-readable bedgraphs to human unreadable bigWigs to save space
## .bedgraph --> .bigWig
# Here we convert the FPM-normalised and sorted bams to bigWigs

give_start_time "8.5 - bedGraph to BigWig 2"
for replicate in ${Replicates}
do
	for drug in ${Drugs}
	do
		for timepoint in ${Timepoints}
		do
			sample_start ${replicate}${drug}${timepoint}
			bedGraphToBigWig bedgraphs/${replicate}${drug}${timepoint}_FPM_sorted.bedgraph \
			CHM13v2chromSizes.txt bigWigs/${replicate}${drug}${timepoint}_FPM.bigWig
			sample_done ${replicate}${drug}${timepoint}
		done
	done
done
give_end_time 8.5
# # ------------------------------------------------------------------------------
## 9 - Feature counts
## Counts how many reads map to each mRNA using the exon coordinates in the
## chm13v2_chrs.gtf file and returns a text file
## .txt

give_start_time "9 - Feature counts"

for drug in ${Drugs}
do
	for timepoint in ${Timepoints}
	do
		sample_start ${drug}${timepoint}

		featureCounts -p -B -M -g gene_id -a chm13v2_chrs.gtf -o \
		${drug}${timepoint}_BRs_mRNA_featureCounts.txt bams/R1${drug}${timepoint}_sorted.bam \
		bams/R2${drug}${timepoint}_sorted.bam bams/R3${drug}${timepoint}_sorted.bam

		sample_done ${condition}
	done
done
give_end_time 9

# # ------------------------------------------------------------------------------
## 10 – Combine replicates
## Combines the replicates into one file and gives the average normalised
## expression
## .bedgraph --> .ubg --> .bedgraph

give_start_time "10 - Combining replicates"

for drug in ${Drugs}
do
	for timepoint in ${Timepoints}
	do
		sample_start ${drug}${timepoint}
		bedtools unionbedg -i bedgraphs/R*${drug}${timepoint}_FPM_sorted.bedgraph > \
		bedgraphs/${drug}${timepoint}_mRNAseq_BRs.ubg
		awk '{print $1,$2,$3,$4}' bedgraphs/${drug}${timepoint}_mRNAseq_BRs.ubg > \
		bedgraphs/${drug}${timepoint}_mRNAseq_BRs.bedgraph
		sample_done ${condition}
	done
done
give_end_time 10

# # ------------------------------------------------------------------------------
## 10.5 – Sorting replicate merged bedgraphs
## bedgraph --> bedgraph

give_start_time "7.5 - Sorting normalised bedgraphs"

for replicate in ${Replicates}
do
  for drug in ${Drugs}
  do
    for timepoint in ${Timepoints}
    do
      sample_start ${drug}${timepoint}
      sort -k1,1 -k2,2n bedgraphs/${drug}${timepoint}_mRNAseq_BRs.bedgraph \
			> bedgraphs/${drug}${timepoint}_mRNAseq_BRs_sorted.bedgraph
      sample_done ${drug}${timepoint}
    done
  done
done

give_end_time 7.5

# # ------------------------------------------------------------------------------
## 11 – BedGraph to BigWig again
## Converts human-readable bedgraphs to human unreadable bigWigs to save space
## .bedgraph --> .bigWig

give_start_time "11.5 - bedGraph to BigWig"
for drug in ${Drugs}
do
	for timepoint in ${Timepoints}
	do
		sample_start ${replicate}${drug}${timepoint}
		bedGraphToBigWig bedgraphs/${drug}${timepoint}_mRNAseq_BRs_sorted.bedgraph \
		CHM13v2chromSizes.txt bigWigs/${drug}${timepoint}_mRNAseq_BRs.bigWig
		sample_done ${replicate}${drug}${timepoint}
	done
done
give_end_time 11.5
