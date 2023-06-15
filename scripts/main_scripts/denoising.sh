#!/bin/bash

# The DnoisE program take into account both the entropy values of the codon positions in the COI fragment included and
# the abundance of each dereplicated sequence. Here, the pipeline is conducted separately for each sample until the primers
# removal, and then on the entire dataset for additional quality control, dereplication and aligning. This means that the
# abundance annotation given by vsearch to the dereplicated sequences refer to the abundance in the entire dataset,
# not exclusively on each sequencing run, the correct value that DnoisE needs to identify "mother" and "daughter"
# sequences, and thus must be corrected for each group of samples defined by the "file_info" file (fifth column).

# The following for loop does that. For every group of samples, the dereplicated sequences occurring exclusively in
# those samples are subtracted from the entire dereplicated dataset, the "dataset_filterlen.fasta" file, which contains
# each (non-dereplicated) sequence with the sample annotation in the header, is subsetted to include only the reads of
# those samples and a new otu table is produced (using the same command of vsearch "search_exact") with the correct
# abundances. The entropy values are calculated on the subset of dereplicated seqs and the DnoisE program is run.
# As the dereplicated sequences are also aligned, due to incompatibilities of some vsearch commands with gaps,
# additional control steps must be conducted, including the seqkit program.

for i in $(cut -f 5 ../file_info.tsv | tail +2 | sort | uniq)

do
	
	cr=$(echo "$(nproc) -1" | bc)
	# reset the number of cores needed
	
	colnum=$(cat align/ESV_counts.tab |\
				head -n1 |\
				tr '\t' '\n' |\
				grep -n -f <(awk -F"\t" -v i="$i" '$5 == i {print $0}' ../file_info.tsv |\
								cut -f 4) - |\
				cut -d ":" -f 1 |\
				tr '\n' ',' |\
				sed 's/,$//')
	# the numbers of the columns of each group of samples is stored in a variable, this will automatically tell "cut" which 
	# columns must use to extract the samples from the total otu table
	
	mkdir denoising/$i
	# make a directory for each group of samples
	
	cut -f"$colnum" align/ESV_counts.tab |\
	tail +2 |\
	sed 's/\t/+/g' |\
	cut -d '|' -f 2- | bc |\
	paste - <(cut -f1 align/ESV_counts.tab | tail +2) |\
	grep -vw "0" |\
	cut -f2 |\
	sed 's/$/;/' |\
	grep -A1 --no-group-separator -f - align/align_noref.fasta |\
	sed 's/;.*$//' > denoising/$i/align_noref_"$i".fasta
	
	# the total otu table is filtered to include only the columns corresponding to the wanted samples and, for each row
	# (dereplicated and aligned seqs) the total amount of reads is calculated (thus, the amount of reads for those
	# particular samples). The name of the seq is pasted alongside and grep is used to extract those with at least one read
	# and those names are used by grep to subset the total dereplicated aligned seqs from all the other, here absent seqs.
	# The abundance annotation is removed from the subset as it will be later included with the correct one.
	# The sed 's/$/;/' is included as it identifies the end of the seq name, which might "confuse" grep with similar names
	# (e.g. "Uniq1" and "Uniq10", would be grepped altogether if searched with "Uniq1", but not if searched with "Uniq1;")

	seqkit seq -w0 filterlen/dataset_filterlen.fasta |\
	grep -A1 --no-group-separator -f <(awk -F"\t" -v i="$i" '$5 == i {print $0}' ../file_info.tsv |\
											cut -f4 |\
											sed -e 's/$/./' -e 's/^/>/') - > denoising/$i/dataset_filterlen_group"$i".fasta
	
	# here the sample names for a specific "group" are used to subset the total fasta file based on header names, which,
	# at the "filterlen" point include headers with the sample name separated by a number using the "." character.
	# The resulting dataset includes the reads of only those specific samples before the dereplication step
	
	vsearch --search_exact denoising/$i/dataset_filterlen_group"$i".fasta \
				--db <(seqkit seq -g denoising/$i/align_noref_"$i".fasta) \
				--sizeout \
				--otutabout denoising/$i/ASV_counts_group"$i".tab \
				--dbmatched denoising/$i/align_noref_"$i"_resize.fasta
	
	# vsearch will search each dereplicated and aligned sequence in the subsetted "filterlen" dataset, recaltulating
	# the correct abundance (sizeout) and printing a new otu table. The dbmatched will create a new file with all the 
	# unique seqs and the correct abundance. However, vsearch cannot compute an exact similarity between aligned 
	# sequences and non aligned, so the gaps are removed beforehand with seqkit. The same program will be used
	# after to replace the headers of old aligned fasta file with the new correct headers
	
	sed -i 's/#OTU ID//' denoising/$i/ASV_counts_group"$i".tab
	# remove the "#OTU ID" data item, as before
	
	vsearch --sortbysize denoising/$i/align_noref_"$i"_resize.fasta \
				--output denoising/$i/align_noref_"$i"_resize_sort.fasta
	
	# sort the seqs by decreasing abundance
	
	seqkit replace -p '(.+)$' -r '{kv}' -I 1 -k \
		<(paste <(grep '^>' denoising/$i/align_noref_"$i"_resize_sort.fasta |\
					sed 's/>//' |\
					cut -d ";" -f1) \
				<(grep '^>' denoising/$i/align_noref_"$i"_resize_sort.fasta |\
					sed 's/>//')) \
		denoising/$i/align_noref_"$i".fasta \
		-o denoising/$i/align_noref_"$i"_corrected.fasta
	
	# seqkit will replace (in this case add) the size abundance in the old dereplicated and aligned sequences with the new
	# abundance obtained by the previous vsearch command. As the order of the sequences might be different, seqkit is
	# preferred as it uses regular expressions and a "key" file that specifies the match pattern to search for. This is
	# provided by the same new file of dereplicated sequences, using grep and including the name without and with the "size"
	# field. A new file is produced, with the aligned sequences and corrected headers

	vsearch --sortbysize <(sed 's/-/N/g' denoising/$i/align_noref_"$i"_corrected.fasta) \
				--output denoising/$i/align_noref_"$i"_corrected_sort.fasta
	
	# the file is then sorted again based on the new abundance values replacing the gap character with N
	# couldn't be done before as the previous vsearch command was searching exact, identical sequences

	seqkit seq -w 0 denoising/$i/align_noref_"$i"_corrected_sort.fasta |\
	sed 's/N/-/g' > denoising/$i/align_noref_"$i"_corrected_sort_oneline.fasta
	
	# the fasta file is transformed to single line and Ns are translated as gaps again
	
	# Entropy calculation, DnoisE
	
	# Here, -x indicate the position of the codon of the first base of each sequence (which, thanks to the refinement done
	# earlier with MEGAX, we know is the first). Python 3.8 is used, but any other version that support the programs
	# and its modules can be used (https://github.com/adriantich/DnoisE).

	python3.8 ~/bin/DnoisE/src/DnoisE.py -g \
				--fasta_input denoising/$i/align_noref_"$i"_corrected_sort_oneline.fasta \
				-x 1 \
				--csv_output denoising/$i/$i

	# After calculating the entropy of each sequence, we can run DnoisE. Entropy values are stored in a variable.
	# The program is run using the default criterion ("join by the lesser [abundance ratio / beta(d)]") as mentioned in -j
	# The alpha is set to 5 as default
	
	entropy=$(cat denoising/$i/"$i"_entropy_values.csv | tail -n 1 | cut -d "," -f 4-)

	python3.8 ~/bin/DnoisE/src/DnoisE.py -a 5 -c $cr \
		--fasta_input denoising/$i/align_noref_"$i"_corrected_sort_oneline.fasta \
		-x 1 \
		-e $entropy \
		--csv_output denoising/$i/denoised_$i \
		--fasta_output denoising/$i/denoised_$i \
		-j 1 \
		-y T
	
	# Next, format the csv ouput of DnoisE.py in order to resemble the output of swarm. This will allow us to calculate the
	# total abundance (keeping the info on the distribution between the samples) for each "mother" sequence which aggregated
	# other "daughter" sequences.

	for y in $(cat denoising/$i/denoised_"$i"_Adcorr_denoising_info.csv |\
				cut -d ',' -f 6 |\
				tail -n +2 |\
				sort |\
				uniq)
	do
		
		awk -F"," -v y="$y" '$6 == y {print $0}' denoising/$i/denoised_"$i"_Adcorr_denoising_info.csv |\
				 cut -d ',' -f 1 |\
				 tr '\n' ' ' |\
				 sed -e 's/^/\n'$y' /' >> denoising/$i/denoise_aggr
		
	done

	sed -i -e 's/$/\n/' denoising/$i/denoise_aggr
	
	cat denoising/$i/denoised_"$i"_Adcorr_denoising_info.csv |\
		cut -d ',' -f 6 |\
		tail -n +2 |\
		sort |\
		uniq |\
		grep -vw -f - <(awk -F"," '$6 == "" {print $0}' denoising/$i/denoised_"$i"_Adcorr_denoising_info.csv |\
						cut -d ',' -f1) >> denoising/$i/denoise_aggr
	
	sed -i '/^$/d' denoising/$i/denoise_aggr

	# aggregate each otu table based on "clustering" (in this case, aggregation by denoising) done by DnoisE.
	# This will be done using a slightly modified version of the "owi_recount_swarm.R" script by Owen S. Wangensteen
	# in GitHub Project Metabarpark 2017 (https://github.com/metabarpark/R_scripts_metabarpark)
	
	Rscript ../scripts/recount_swarm.R denoising/$i/denoise_aggr \
										denoising/$i/ASV_counts_group"$i".tab
	
done

seqruns_denoised=$(ls denoising/[1-9]*/denoised_[1-9]*_Adcorr_denoised* | tr '\n' ' ' | sed 's/ $/\n/')
# the filename of the fasta files of denoised sequences is stored in a variable separated by a space, indicating to the
# following cat command to concatenate those files

cat $seqruns_denoised > denoising/denoised_total.fasta

vsearch --derep_fulllength <(seqkit seq -g -w0 denoising/denoised_total.fasta) \
			--sizein \
			--sizeout \
			--output denoising/denoised_total_derep.fasta
# The concatenated denoised sequences are dereplicated again, taking the abundance info into account (--sizein) and
# giving the new abundance in the "size" field of the header (--sizeout) summing the abundances of the dereplicated
# sequences shared by samples of different runs.

seqruns_asvcount=$(ls -d denoising/*/ | sed 's/denoising\///' | sed 's/\///' | sort -g | tr '\n' ' ' | sed 's/ $/\n/')
# the same approach is performed for the otu tables of each different run.

Rscript ../scripts/aggregate_ASVcounts.R $seqruns_asvcount
# This script will, in addition to concatenating the otu tables, remove low abundance ASVs

seqkit seq -w0 denoising/denoised_total_derep.fasta |\
	grep -A1 --no-group-separator -f <(cat denoising/ASV_counts_total.tab |\
										tail +2 |\
										cut -f 1 |\
										sed -e 's/^/>/' -e 's/$/;/') > denoising/denoised_total_derep_abund.fasta
# the denoised and dereplicated sequences will be filtered to remove those with a low abundance as previously done

