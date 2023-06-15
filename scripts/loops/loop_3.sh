#!/bin/bash

(

set -ue pipefail

for FNAME in $(awk -F"\t" '$2 == "3" {print $0}' ../file_info.tsv |\
				cut -f 1 |\
				sort |\
				uniq)

do

mkdir fastqc/${FNAME}

echo "FastQC on ${FNAME}"
fastqc ../raw/${FNAME}_R1.fastq \
		../raw/${FNAME}_R2.fastq \
		-outdir fastqc/${FNAME}/ \
		-d temp/ \
		-q

echo "Vsearch trimming on forward reads"
vsearch --fastq_filter ../raw/${FNAME}_R1.fastq \
		--fastq_trunclen_keep 200 \
		--fastqout ../raw/trimmed_R1/${FNAME}_R1_trim.fastq

echo "Vsearch trimming on reverse reads"
vsearch --fastq_filter ../raw/${FNAME}_R2.fastq \
		--fastq_trunclen_keep 190 \
		--fastqout ../raw/trimmed_R2/${FNAME}_R2_trim.fastq

mkdir merging/${FNAME}
mkdir merging/${FNAME}/Not_merged

echo "Vsearch mergepairs ${FNAME}"
vsearch --fastq_mergepairs ../raw/trimmed_R1/${FNAME}_R1_trim.fastq \
		--reverse ../raw/trimmed_R2/${FNAME}_R2_trim.fastq \
		--fastq_maxmergelen 313 --fastq_minmergelen 311 \
		--fastq_maxns 0 \
		--fastq_maxdiffs 2 \
		--fastqout_notmerged_fwd merging/${FNAME}/Not_merged/${FNAME}_R1_notmerged.fastq \
		--fastqout_notmerged_rev merging/${FNAME}/Not_merged/${FNAME}_R2_notmerged.fastq \
		--fastqout merging/${FNAME}/${FNAME}_merged.fastq \
		--log merging/${FNAME}/${FNAME}_merge.log

# set variable with replicate name
y=$(awk -F"\t" '$2 == "3" {print $0}' ../file_info.tsv |\
					awk -F"\t" '$1 == "'$FNAME'" {print $0}' |\
					cut -f 4)

echo "Vsearch relabel with ${FNAME}.#"
vsearch --fastq_filter merging/${FNAME}/${FNAME}_merged.fastq \
		--relabel ${y}. \
		--fastqout fastq_prefilter/${y}.fastq

done

) 2>&1 | tee log/loopoutput_3.log
