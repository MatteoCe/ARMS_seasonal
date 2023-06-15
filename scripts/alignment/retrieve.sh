#!/bin/bash

set -u
set -e
set -x
#set -o pipefail

# this script will retrieve nucleotide sequences and informations from the NCBI of the accession numbers listed in the
# pre-compiled file "accn" in the scripts/alignment subfolder. The accession numbers correspond to mitogenomes of different
# metazoan phyla commonly found in benthic sessile metazoans. As we are only interested in the COX 1 region, the first
# step will list all the different namings of the regions of each mitogenome and print it on screen. The user can copy each
# useful naming of the region of interest (in this case coi, co1, COXI etc.) and write them on screen separated by single
# spaces. These will be saved in another text file (align/retrieve/genes), along with the start and end nucleotide position
# of each region, which will be used to download the exact COX 1 sequence.

# create new directory for the mitogenomes files
mkdir align/retrieve

# this for loop will search the accession numbers in the "accn" file and download the associated information regarding
# the gene name provided in the INSD  "Integrative Nucleotide Sequence Database". The namings of all regions associated
# to that accession number will be listed in the terminal, along with the number of times those namings occur.

for i in $(cat ../scripts/alignment/accn)

do

	efetch -db nucleotide -id "$i" -format gbc |\
	xtract -insd CDS gene >> align/retrieve/genes

done

echo
sort align/retrieve/genes -r -k 2,2 | awk -F"\t" '$2 != "" {print $0}' | cut -f2 | sort | uniq -c

# this function prompts the user for a region of interest (genes), reads that input, and then retrieves specific information
# from the NCBI. The user must enter the correct namings for the region of interest, separated by single spaces.
# As the mitogenomes in the "accn" file refer to more than COX 1 sequences, the exact coordinates of the starting and end
# positions of cytochrome sequences will be extracted written to the new file "coord", along with the accession number.

function genomes_interest() {

echo
read -p "Show region of interest" genes

	if [ "$genes" == "" ]; then
		genomes_interest
	else

		for i in $(cat ../scripts/alignment/accn)

		do

			filt="$(echo -n "$genes" |\
					tr ' ' '\n' |\
					sort -f |\
					uniq -i |\
					tr '\n' ' ' |\
					sed -e 's/ $//' -e 's/\ / -or INSDQualifier_value -equals /g')"
			efetch -db nucleotide -id "$i" -format gbc |\
			xtract -pattern INSDFeature -tab ":" -element INSDInterval_accession \
			-group INSDFeature -if INSDFeature_key -equals CDS -and INSDQualifier_name -equals gene \
			-and INSDQualifier_value -equals ''$filt'' \
			-unless INSDInterval_iscomp -or INSDFeature_operator -equals join \
			-or INSDQualifier_name -equals gene_synonym \
			-tab "" -element INSDInterval_from -lbl "-" -element INSDInterval_to >> align/retrieve/coord

		done

	fi


}

genomes_interest

# Accession number and coordinates will be used to download the extract the translation table code, which, according to
# the feature table adopted by the INSD, defines "the genetic code table used if other than the universal genetic code table"
# Later in the pipeline, thus not in this script, that information will be used by MACSE to properly translate the genomes
# sequences to the correct amino acid translations.

for i in $(cat align/retrieve/coord | grep ':' | cut -d ':' -f 1)

		do

			from=$(cat align/retrieve/coord | grep ':' | grep "$i" | cut -d ':' -f 2 | cut -d '-' -f1)
			to=$(cat align/retrieve/coord | grep ':' | grep "$i" | cut -d ':' -f 2 | cut -d '-' -f2)
			filt="$(echo -n "$genes" |\
					tr ' ' '\n' |\
					sort -f |\
					uniq -i |\
					tr '\n' ' ' |\
					sed -e 's/ $//' -e 's/\ / -or INSDQualifier_value -equals /g')"
			efetch -db nucleotide -id "$i" -format gbc |\
			xtract -pattern INSDFeature -tab ":" -element INSDInterval_accession \
			-group INSDFeature -if INSDFeature_key -equals CDS -and INSDInterval_from -equals ''$from'' -and INSDInterval_to -equals ''$to'' \
			-unless INSDQualifier_name -equals gene_synonym \
			-subset INSDQualifier -if INSDQualifier_name -equals transl_table -tab "|" -element INSDQualifier_value >> align/retrieve/accn_gencode

		done

# The entire mitogenome sequences are then retrieved and written to a fasta file
cat align/retrieve/coord | grep ':' | cut -d ':' -f 1 |\
efetch -db nucleotide -format fasta |\
sed '/^>/ s/ .*//' > align/retrieve/genome.fa

# Samtools will use the coordinates to cut the mitogenome sequences to the correct cytochrome regions
# (faidx command with -r option)
samtools faidx align/retrieve/genome.fa \
		-r <(cat align/retrieve/coord | grep ':') > align/retrieve/genome_coord.fa

# seqkit will ensure the sequences are one-line fasta (-w 0), all upper-case (-u) and without gaps (-g)
seqkit seq -ug -w0 align/retrieve/genome_coord.fa |\
	sed '/^>/s/\..*//g' > align/retrieve/clean_genome_coord.fasta
