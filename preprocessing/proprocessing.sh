#!/usr/bin/bash

# Date
DATE_NOW=$(date +'%d%m%Y')

# Set working directory
WORKING_DIRECTORY="/path/to/directory"

# Dataset plain text file
DATASET="/path/to/dataset/dataset.txt"

# Names of strains / sequences
STRAIN_NAMES="/path/to/strainnames/strains.csv"

# Raw read data
RAW_DATA="/path/to/raw_reads"

# Reference sequence
reference_D23580="/path/to/reference/reference.fa"


# GFF3 feature file for reference
gff3_annotated_D23580="/path/to/annotation/data.gff.gz"

# Chromosome 
CHROMOSOME='Chromosome'

# UDP-glucose pyrophosphorylase (galF)
GALF_=2183016

# gluconate-6-phosphate dehydrogenase (gnd)
GND_=2203520

# Reads file extension
READ_EXTENSION=".fastq.gz"

# Downstream analysis directory
DOWNSTREAM_ANALYSIS_DIR=${WORKING_DIRECTORY}/st313_analysis_${DATE_NOW}

# BAM sorted and indexed file directory
BAM_DIR_SORTED_INDEXED=$DOWNSTREAM_ANALYSIS_DIR/1_bam_sorted_indexed_

# BCF file directory
BCF_DIR=$DOWNSTREAM_ANALYSIS_DIR/2_bcf_

# VCF file directory
VCF_DIR=$DOWNSTREAM_ANALYSIS_DIR/3_vcf_

# VCF annotated directory
VCF_ANNOTATED_DIR=$DOWNSTREAM_ANALYSIS_DIR/4_vcf_annotated_

# Filtered variants directory
VARIANTS_DIR=$DOWNSTREAM_ANALYSIS_DIR/5_variants_

# Consensus O-Antigen sequences for all strains
OANTIGEN_FASTA_DIR=$DOWNSTREAM_ANALYSIS_DIR/6_oantigen_fasta_

# ORFs for all consensus sequences
OANTIGEN_ORFS_DIR=$DOWNSTREAM_ANALYSIS_DIR/7_oantigen_orfs_

# ORFs for all consensus sequences
OANTIGEN_ORFS_TRANS_DIR=$DOWNSTREAM_ANALYSIS_DIR/8_oantigen_orfs_trans_

# Protein translated grouped orthologous genes
ORFS__TRANS_ORTHOS_DIR=$DOWNSTREAM_ANALYSIS_DIR/9_oantigen_trans_orfs_orthologs_

# Aligned protein orthologous genes
ORFS__TRANS_ALIGNED_DIR=$DOWNSTREAM_ANALYSIS_DIR/10_oantigen_aligned_trans_orfs_orthologs_

# Aligned protein orthologous genes
ORFS_ORTHO_NUCL_DIR=$DOWNSTREAM_ANALYSIS_DIR/11_oantigen_orfs_orthologs_nucl_

# Codon aligned orthologous genes
ORFS__TRANS_CODON_ALIGNED_DIR=$DOWNSTREAM_ANALYSIS_DIR/12_oantigen_codon_aligned_trans_orfs_orthologs_

# Evolutionary analysis
DNDS_RATIO_DIR=$DOWNSTREAM_ANALYSIS_DIR/13_dnds_ratio_

# Final results
ANALYSIS_OUTPUT_DIR=$DOWNSTREAM_ANALYSIS_DIR/analysis_results_

# Directories array
pipeline_dirs=("${DOWNSTREAM_ANALYSIS_DIR}" "${BAM_DIR_SORTED_INDEXED}" "${BCF_DIR}" "${VCF_DIR}" "${VCF_ANNOTATED_DIR}" "${VARIANTS_DIR}" "${OANTIGEN_FASTA_DIR}" "${OANTIGEN_ORFS_DIR}" "${OANTIGEN_ORFS_TRANS_DIR}" "${ORFS__TRANS_ALIGNED_DIR}" "${ORFS_ORTHO_NUCL_DIR}" "${ORFS__TRANS_CODON_ALIGNED_DIR}" "${DNDS_RATIO_DIR}" "${ANALYSIS_OUTPUT_DIR}")

# Create directories
for dir in "${pipeline_dirs[@]}"; do
  rm -rf "$dir"
  mkdir -p "$dir"
done

# Loop through sequences in the target directory and for each sequence
# 1) Align each fasta file of contigs to reference genome with output SAM file
# 2) Convert SAM file to BAM file
# 3) Sort then index BAM file
# 4) Pile up BAM file with output BCF file
# 5) Call variants on pile up BCF file

# Initialize loop counter
counter=1

while IFS= read -r read_entry; do
	
	if [[ $read_entry =~ ^ERR.* || $read_entry =~ ^SRR.* ]]; then
		
		# Look up strain name from final dataset csv file based on accession number
		filename=$(awk -F, -v id="$read_entry" '$3 == id {print $1}' $STRAIN_NAMES)

		# BAM file
		BAM_FILE_SORTED_INDEXED="${BAM_DIR_SORTED_INDEXED}/${filename}_sorted_indexed.bam"

		bcf_="${BCF_DIR}/${filename}.bcf" 
		
		vcf_file="${VCF_DIR}/${filename}.vcf"

		gzvcf_file="${VCF_DIR}/${filename}.vcf.gz"

		# Variant effect prediction
		vep_output="${VCF_ANNOTATED_DIR}/${filename}.txt"

		# File Filtered alleles to O-Antigen region
		filtered_variant="${VARIANTS_DIR}/${filename}.txt"

		oag_region="${OANTIGEN_FASTA_DIR}/${filename}_oAg_region.fasta"

		orfs_="${OANTIGEN_ORFS_DIR}/${filename}_orfs.fasta" 

		orfs_trans_="${OANTIGEN_ORFS_TRANS_DIR}/${filename}_orfs_trans.faa" 

		if [[ ${read_entry} =~ ^ERR.* ]]; then 

			echo "...Processing strain # ${counter}: ${read_entry}"
			
			# Read 1					
			read1="${RAW_DATA}/${read_entry}_1${READ_EXTENSION}"

			# Read 2
			read2="${RAW_DATA}/${read_entry}_2${READ_EXTENSION}"

			if [[ -e $read1 && -e $read2 ]]; then 

				# ToDO

				# Processing paired reads
				minimap2 -ax sr $reference_D23580 $read1 $read2 | samtools sort -o $BAM_FILE_SORTED_INDEXED && samtools $BAM_FILE_SORTED_INDEXED

			fi

		elif [[ ${read_entry} =~ ^SRR ]]; then

			echo "...Processing strain # ${counter}: ${read_entry}"

			read="${RAW_DATA}/${read_entry}${READ_EXTENSION}"

			if [[ -e $read ]]; then 

				minimap2 -ax sr $reference_D23580 $read | samtools sort -o $BAM_FILE_SORTED_INDEXED && samtools $BAM_FILE_SORTED_INDEXED
			fi
			
		fi

		if [[ -e $BAM_FILE_SORTED_INDEXED ]]; then
			# Piling reads and calling variants
			bcftools mpileup -Ou -f $reference_D23580 $BAM_FILE_SORTED_INDEXED | bcftools call -O b --threads n -vc --ploidy 1 -p 0.05 | bcftools sort -Ob -o $bcf_  && bcftools index $bcf_


			# Filter variants to identified OAg region
			bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\n' $bcf_ -r "$CHROMOSOME:$GALF_-$GND_" -o $filtered_variant


			# Change BCF file to VCF - subset to O-Antigen region
			bcftools view $bcf_ -r "$CHROMOSOME:$GALF_-$GND_" -Ov -o $vcf_file


			# Bgzip variant file
			bgzip -c $vcf_file > $gzvcf_file

			# Index Bgzip vcf file
			bcftools index $gzvcf_file

			# Run VEP to predict and annotate variant effects on subset Bgzip vcf file
			/usr/bin/perl /opt/ensembl-vep-release-113/vep --input_file $gzvcf_file --output_file $vep_output --fasta $reference_D23580 --gff $gff3_annotated_D23580 --species bacteria --force_overwrite


			# Subset consensus whole genome to only O-Antigen region
			samtools faidx $reference_D23580 "$CHROMOSOME:$GALF_-$GND_" | bcftools consensus $gzvcf_file > $oag_region

			sed -i "s/^>\(.*\)/>${filename}/" $oag_region

			# Convert consensus sequence to nucleotide open reading frame
			getorf -sequence $oag_region -outseq $orfs_ -table 11 -find 2 -minsize 150 -reverse false

			# Translate open reading frame to protein sequence
			transeq -sequence $orfs_ -outseq $orfs_trans_ -table 11 -clean	

		fi

		counter=$((counter + 1))

	fi 

done < "$DATASET"

ulimit -n 300000

# Group orthologous genes
/opt/OrthoFinder/orthofinder -f ${OANTIGEN_ORFS_TRANS_DIR} -o ${ORFS__TRANS_ORTHOS_DIR} -og

# Align each group of orthologous genes
find ${ORFS__TRANS_ORTHOS_DIR} -type d -name "Orthogroup_Sequences" | while read -r orthodir; do

  for ortho_file in "$orthodir"/*; do

	ortho_name="${ortho_file##*/}" && ortho_name="${ortho_name%.*}"

	mafft --auto $ortho_file > "${ORFS__TRANS_ALIGNED_DIR}/${ortho_name}.fa"

	sed -i 's/_[0-9]\+$//' "${ORFS__TRANS_ALIGNED_DIR}/${ortho_name}.fa"

	ortho_nucl="${ORFS_ORTHO_NUCL_DIR}/${ortho_name}_nucl.fa"

	# Clean or create the output file
	touch "${ortho_nucl}" 

	# codon aligned 
	codon_aligned="${ORFS__TRANS_CODON_ALIGNED_DIR}/${ortho_name}_codon_aligned.fa"

	# Loop through each header in input orthofile and get the nucleotide sequence from orf file
	grep '^>' "${ORFS__TRANS_ALIGNED_DIR}/${ortho_name}.fa" | while read -r header; do

		id="${header#>}"

		# Extract file prefix (adjust this logic based on your header structure)
		file_prefix=$(echo "$id" | cut -d'_' -f1)  # change -f1-2 if needed

		# Find the file whose name starts with the prefix
		fasta_file=$(find "$OANTIGEN_ORFS_DIR" -type f -name "${file_prefix}_orfs.fasta" | head -n 1)

		if [[ -z "$fasta_file" ]]; then
			echo "No matching file found for $file_prefix"
			continue
		fi

	
		awk -v id="$id" 'BEGIN{RS=">";ORS=""} NR>1{h=substr($0,1,index($0,"\n")-1); sub(/ .*/,"",h); s=substr($0,index($0,"\n")+1); if(h==id) print ">" h "\n" s "\n"}' "$fasta_file" >> "$ortho_nucl"

	done

	# Codon align using Pal2nal
	/opt/pal2nal.v14/pal2nal.pl $ortho_file $ortho_nucl -output fasta > ${ORFS__TRANS_CODON_ALIGNED_DIR}/${ortho_name}_codon.fasta

	sed -i 's/^>\([^_]*\)_.*/>\1/' ${ORFS__TRANS_CODON_ALIGNED_DIR}/${ortho_name}_codon.fasta


  done

done



# Concatenate codon orthologs
/usr/bin/python3.10 /opt/AMAS/amas/AMAS.py concat -i ${ORFS__TRANS_CODON_ALIGNED_DIR}/*codon.fasta -f fasta -d dna -u phylip -t ${DNDS_RATIO_DIR}/codon_concatenated.phy

# Newick tree
iqtree2 -s ${DNDS_RATIO_DIR}/codon_concatenated.phy -m GTR+G -nt AUTO -pre ${DNDS_RATIO_DIR}/