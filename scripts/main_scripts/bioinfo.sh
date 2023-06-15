#!/bin/bash

# The first part will change the file names leaving only the sample IDs. After creating a directory with the project name,
# cd in that directory and create new directories for every marker analyzed.

# set the directories structure
#mkdir -p ./{sites_selection,\
			scripts,\
			statistics_analyses,\
			raw/{trimmed_R1,\
				trimmed_R2},\
			bioinfo_analyses/{align,\
							denoising,\
							chimera_removal/chimeras,\
							clustering,\
							fastqc,\
							temp,\
							log,\
							cutadapt,\
							merging,\
							fasta/discarded,\
							fastq_prefilter,\
							fastq_prefilter/dataset_prefilter/quality_log,\
							filterlen/discarded,dereplication}}

# Download Pearman 2020 data from the sra database. The id of the project reported in the pearman paper (PRJNA557002)
# is searched in the sra database and the runInfo file is downloaded from each link. The csv file of the summary contains
# the different ids of the files and the corresponding name in the paper. However they result less then those reported in
# the supplementary files of the paper. 45 against the 49 reported in the supplementary.
# Those file are here reported: CHE.A1.Sessile, GA.A1.Sessile, GA.A2.Sessile, PAS.A1.Sessile.

# fasterq-dump download Pearman et al. 2020
#for i in $(cut -f1 scripts/data_download/Pearman_etal_2020/accn); do fasterq-dump --split-files $i -O raw/; done

# Carvalho et al 2019 will use the same method

# fasterq-dump download Carvalho et al. 2019
#for i in $(cut -f1 scripts/data_download/Carvalho_etal_2019/accn); do fasterq-dump --split-files $i -O raw/; done

# Nichols et al 2021 have to be downloaded manually from the Dryad repository

# The ARMS-MBON samples will be downloaded from the ENA database, using wget and a list of ftp links associated
# to the reads of each sample

#wget -i scripts/data_download/armsmbon/ftp -P raw/

# unzip, if needed

# Then rename each read as R1 and R2

#rename 's/_1/_R1/' raw/*.fastq
#rename 's/_2/_R2/' raw/*.fastq

# A series of loop scripts will start.
# It will use fastqc for quality control, vsearch for merging paired ends and trim the end of sequences and cutadapt to
# remove primers. If activating a conda environment whithin a script, remember to add -i to bash
# (bash -i ../scripts/loop_X.sh) to run the script interactively and allowing the activation of the conda
# environment for cutadapt.

# In this case I'll use different loops for the three set of samples: East Antarctic Wilkes land, Carvalho et al. (2019)
# (with additionally some samples from the ARMS MBON program) and the rest of the dataset.
# The loop for the East Antarctic Wilkes land samples simply include a step that reverse complement the reads before the
# other steps.
# The other two loops are designed for samples that present or not the primers, and thus have a primer removal step or
# lack it.

# Infos on the specific samples treated in the different loops are presented in the file_info.tsv file 

cd bioinfo_analyses

for j in $(cat ../file_info.tsv | cut -f 2 | tail +2 | sort | uniq)
do

bash -i ../scripts/loops/loop_"$j".sh

done

# merge all cutadapt reports into a single one in the cutadapt main directory
cat cutadapt/*/*_report.txt > cutadapt/cutadapt_reports.txt

# concatenate all samples with renamed merged reads
cat fastq_prefilter/*.fastq > fastq_prefilter/dataset_prefilter/dataset_prefilter.fastq

# remove sequences with > 1 maximum expected error and convert to fasta
vsearch --fastq_filter fastq_prefilter/dataset_prefilter/dataset_prefilter.fastq \
			--fastq_maxee 1 \
			--fastaout fasta/dataset.fasta

# Next i will use vsearch for filtering the dataset according to the lenght of the sequences, dereplicate it and
# relabel the headers and remove the singletons. The length filter is set to 313 - 311 to account for small differences
# in the length of the leray region.
vsearch --fastx_filter fasta/dataset.fasta \
			--fastq_maxlen 313 \
			--fastq_minlen 311 \
			--fastaout filterlen/dataset_filterlen.fasta \
			--log filterlen/log_filter.txt \
			--fastaout_discarded filterlen/discarded/dataset_filtlen_discarded.fasta

# the next step is the dereplication to produce the first batch of unique sequences and relabel the headers with a new ID
vsearch --derep_fulllength filterlen/dataset_filterlen.fasta \
			--sizeout \
			--relabel Uniq \
			--log dereplication/log_derep.txt \
			--output dereplication/unique_seqs.fa

# Now the uchime algorithm will be run to remove the chimeras in the dereplicated dataset.
# uchime3_denovo of vsearch will be run
vsearch --uchime3_denovo dereplication/unique_seqs.fa \
			--nonchimeras chimera_removal/unique_seqs_chim.fa \
			--chimeras chimera_removal/chimeras/chimeras.fa \
			--uchimealns chimera_removal/uchime_aln \
			--uchimeout chimera_removal/chimera_tab \
			--log chimera_removal/uchime3_log.txt

# Now the reliability of the dereplicated, non singleton sequences will be checked by aligning those seqs to COI regions of
# mitochondrial genomes from the NCBI. The accession numbers of the metazoan COI genomes are listed in the "accn" file in the
# "scripts" folder of the project. The "retrieve.sh" will run interactively and download the gene names information from the
# NCBI and ask the user to choose the different coi nomenclatures found. All the sequences for the corresponding
# accession numbers will be downloaded and trimmed for the corresponding coordinates of the COI region using "samtools"
# (http://www.htslib.org/)
# The sequences will also be checked for genetic code informations, as those will be used later by the macse program for a
# correct alignment

bash ../scripts/alignment/retrieve.sh

# simplify the accn_gencode file
cat align/retrieve/accn_gencode | grep ':' | tr ':' '\t' | sed 's/\.[0-9]//' > align/retrieve/accn_gencode_macse

# align the genomes together
java -jar ~/bin/macse/macse_v2.06.jar \
		-prog alignSequences \
		-seq align/retrieve/clean_genome_coord.fasta \
		-gc_file align/retrieve/accn_gencode_macse \
		-out_NT align/genomes_COI_NT.fasta \
		-out_AA align/genomes_COI_AA.fasta

# create a fasta file with the Leray primers and put it in the scripts folder. The next line of code will align the primers
# to the previously formed alignment
java -jar ~/bin/macse/macse_v2.06.jar \
		-prog enrichAlignment \
		-align align/genomes_COI_NT.fasta \
		-seq align/genomes_COI_NT.fasta \
		-seq_lr ../scripts/alignment/primers.fasta \
		-maxTotalINS_inSeq 0

# Check the result using MEGAX and adjust the alignment of the primers (they might be shifted of one or two bases or "cut")
# In this case (with the genomes used here), the forward primer shows a gap in the third last position (in respect to the 3'
# end) which I eliminated considering that it doesn't coindice with the bases showed in the genomes alignment
# (c homopolymers).
# The region inside the primers was extracted, comprising a couple of codons on the extremities (thus reaching a total length
# of 318 bp: 312 (leray only codons) + 1 (remaining leray) + 2 (finish codon with the remaining leray) + 3
# (last additional codon))
# Now all the region extracted is translatable as amminoacid sequence without uncomplete codons. This will help macse with
# the COI alignment. Moreover, will allow us to calculate the entropy levels for each sequence base position thanks to the
# known codon position of each single position of the sequence (i.e. they all start with a complete codon)
# Save the file as align/genomes_refalign.fas

# Align the dereplicated sequences to the reference genome alignment. The retrieve/accn_gencode_macse file will be used to
# translate each genomic sequence according to the right genetic code, whereas all the dereplicated sequences will be
# translated with the fifth code for simplicity. The command will be run by NOT allowing any insertions in the sequences
# any frameshift and any stop codon. Those seqs that have those characteristics will be discarded.

# set name of file to align for simplicity
var2=$(echo chimera_removal/unique_seqs_chim.fa)

# get number of total reads, those will be used to split the fasta file into multiple fasta file to align in parallel
num=$(grep -c '^>' $var2)

# get number of cores and remove one (change based on availability)
cr=$(echo "$(nproc) -1" | bc)

# get number of sequences to split the fasta and align for every processor available
num2=$(echo $num/$cr | bc)

# split the fasta file in $num2 number of sequences each
seqkit split $var2 -s $num2 -w 0 -o chimera_removal/

# the splitted files' names are stored in a variable
list_split=$(ls $var2.split/)

# run the macse alignment program in parallel based on the chosen number of processors available
echo $list_split | tr ' ' '\n' |\
parallel -j $cr 'java -jar ~/bin/macse/macse_v2.06.jar \
						-prog enrichAlignment \
						-align align/genomes_refalign.fas \
						-seq align/genomes_refalign.fas \
						-seq_lr '$var2'.split/{} \
						-gc_file align/retrieve/accn_gencode_macse \
						-out_NT align/alignment_NT_{}.fasta \
						-out_AA align/alignment_AA_{}.fasta \
						-gc_def 5 \
						-maxFS_inSeq 0 \
						-maxINS_inSeq 0 \
						-maxSTOP_inSeq 0 \
						-fixed_alignment_ON'

# concatenate sequences that passed were successfully aligned
cat align/alignment_NT* > align/alignment_NT.fasta

# now remove from the alignment the reference genome sequences
cat align/alignment_NT.fasta |\
	grep -A1 '>Uniq' --no-group-separator |\
	sed -e '/^>/s/$/;/g' -e 's/!/-/g' > align/align_noref.fasta

# remove gaps from the alignment, only to create the total ESV table with the sequences that passed the macse requirements
cat align/align_noref.fasta |\
	sed -e 's/-//g' > align/aligned_unique_nogap.fa

# search exact to reconstruct the otu table for the dereplicated, denoised Exact Sequence Variants
vsearch --search_exact filterlen/dataset_filterlen.fasta \
		--db align/aligned_unique_nogap.fa \
		--otutabout align/ESV_counts.tab \
		--log usearchglobal_report.txt

# remove the "#OTU ID" data item in the first (top leftmost) cell that characterize the "qiime-formatted" otu tables
# produced by vsearch as it conflicts with the following procedures
sed -i 's/#OTU ID//' align/ESV_counts.tab

# The DnoisE program from Antich et al. (2022) use a denoising procedure which is optimized for COI metabarcoding sequences.
# The entropy values for the first, second and third codon positions will be calculated.

# In order to provide an exact evaluation of the entropy values and identification of "daughter" (sequencing errors)
# sequences, the analysis must be conducted exclusively with samples that have been sequenced in the same run and
# preferably that have been processed using the same PCR settings etc.

# For this reason, a for loop will allow to conduct the analysis separately for each group of samples specified in the
# "groups" text file that must be prepared beforehand and identifies the groups of samples that were processed alltogether.

# the analyses are conducted on a separate script

bash ../scripts/main_scripts/denoising.sh

# Now the clustering can be performed on the denoised dataset. This will be done using swarm (Swarm 3.0.0) Mahé et al. (2021)

# remove size info from the denoised fasta file
sed -e 's/;size=/_/g' -e 's/;$//' denoising/denoised_total_derep_abund.fasta > clustering/denoised.fasta

# The clustering will be performed using 13 as "resolution", as suggested by many authors due to the high variability of COI
# see Antich et al. (2021) https://doi.org/10.1186/s12859-021-04115-6

cr=$(echo "$(nproc) -1" | bc)
# set the number of cores

swarm \
	-t $cr \
	-d 13 clustering/denoised.fasta \
	-s clustering/stats \
	-o clustering/clust_swarm \
	-w clustering/clusts.fa

# Remove the _"abundance" values from the cluster names
sed -e 's/_[0-9]*//g' clustering/clust_swarm > clustering/clust_swarm_mod

# Aggregate abundance info for each cluster using the recount_swarm script as above
Rscript ../scripts/recount_swarm.R clustering/clust_swarm_mod \
									denoising/ASV_counts_total.tab

# Copy the abundance table to the directory for the statistical analyses
cp clustering/clust_swarm_mod.counts.tsv statistics_analyses/

