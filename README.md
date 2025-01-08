# Epibase2025

### Download from SRA

```bash
conda activate /usr/local/src/conda_envs_binf

# Download with prefetch
cat sra_list.txt | xargs -n 1 -P 92 prefetch

# Split into _1 and _2 with fasterq-dump
parallel -j 92 fasterq-dump --split-files {} ::: $(cat sra_list.txt)

# old command, much slower! 
# cat sra_list.txt | xargs -n 1 -P 92 fasterq-dump -p

### Python script to change SRA names

#! usr/bin/python3

# rename_sra_files.py
# John Quensen
# 16 July 2019

# Creates a dictionary from a tab-delimited file of SRA ID's and sample names.
# Uses the dictionary to rename fastq files from SRA sequence IDs to sample names.
# Run the script from the directory containing the fastq sequences and the 
# tab-delimited file of SRA fastq file ID's and sample names.
# Files are assumed to be paired reads with names of the form:
# SRR8648700_1.fastq
# SRR8648700_2.fastq
# Lines in the tab-delimited file should look like (no header line):
# SRR8648702	112-20-1
# SRR8648701	112-20-2
# SRR8648700	112-30-1
# SRR8648699	112-30-2

# Usage is: python3 rename_sra_files.py sra_sample_names_file.tsv

import sys
import os

def main(args):
	name_file = sys.argv[1]
	hash={}
	for line in open(name_file, 'r').readlines():
		file_name, sample_name = line.strip().split()
		if(file_name not in hash.keys()):
			hash[file_name.strip()] = [sample_name.strip()]
		else:
			print("Warning: Duplicate file names..")
			sys.exit()
	for old_file_name in os.listdir():
		if old_file_name.endswith(".fastq"):
			sra = old_file_name.split("_")[0]
			suf = old_file_name.split("_")[1]
			sam = hash.get(sra)[0]
			new_file_name = "".join(str(sam) + "_" + str(suf))
			os.rename(old_file_name, new_file_name)

if __name__ == "__main__":
        usage = "python3 rename_sra_files.py sra_sample_names_file.tsv"
        if len(sys.argv) != 2 :
            print("Incorrect number of arguments.\nUsage is : ", usage)
            sys.exit()
        main(sys.argv[1:])
```


## Rename the SRA files

```bash
python rename_sra_files.py sra_sample_names_file.tsv
```

## Replace _1 and _2 with _R1 and _R2

```bash
rename "_1.fastq" "_R1.fastq" *.fastq
rename "_2.fastq" "_R2.fastq" *.fastq
```

# Combine reads from genome skim + seq capture for the same accessions

```bash
#!/bin/bash

# Directory where your fastq.gz files are located
fastq_dir="path/to/your/fastq/files"

# Move to the fastq directory
cd "$fastq_dir"


# Find all unique prefixes (everything before _S or _R fields)
for prefix in $(ls *_R1.fastq.gz | sed 's/_S[0-9].*//g' | sort -u); do
  echo "Processing prefix: $prefix"

  # Get all R1 and R2 files associated with the prefix
  r1_files=(${prefix}*_R1*.fastq.gz)
  r2_files=(${prefix}*_R2*.fastq.gz)

  if [ ${#r1_files[@]} -gt 0 ] && [ ${#r2_files[@]} -gt 0 ]; then
    # Find the longest filename (to preserve genus and species names)
    longest_r1=$(printf "%s\n" "${r1_files[@]}" | awk '{ print length, $0 }' | sort -nr | head -n1 | cut -d" " -f2)
    longest_r2=$(printf "%s\n" "${r2_files[@]}" | awk '{ print length, $0 }' | sort -nr | head -n1 | cut -d" " -f2)

    # Extract the base name (without the R1 or R2)
    base_name_r1=$(echo "$longest_r1" | sed 's/_R1.*//')
    base_name_r2=$(echo "$longest_r2" | sed 's/_R2.*//')

    # Use the longer of the two base names (R1 or R2)
    if [ ${#base_name_r1} -ge ${#base_name_r2} ]; then
      final_base_name="$base_name_r1"
    else
      final_base_name="$base_name_r2"
    fi

    # Combine R1 files
    echo "Combining R1 files for $prefix into ${final_base_name}_combined_R1.fastq.gz"
    cat "${r1_files[@]}" > "${final_base_name}_combined_R1.fastq.gz"

    # Combine R2 files
    echo "Combining R2 files for $prefix into ${final_base_name}_combined_R2.fastq.gz"
    cat "${r2_files[@]}" > "${final_base_name}_combined_R2.fastq.gz"
  else
    echo "Missing R1 or R2 files for $prefix!"
  fi
done

echo "All files processed successfully."
```

# FastP

```bash
## run fastp

for f1 in *_R1.fastq.gz
	do
        f2=${f1%%_R1.fastq.gz}"_R2.fastq.gz"
        fastp -i $f1 -I $f2 -w 16 --trim_poly_g --trim_poly_x -l 25 --cut_right -o "../fastp_polyg/fastp-$f1" -O "../fastp_polyg/fastp-$f2"
	done


# Captus
```bash
conda activate /home/cfb0001/.local/share/mamba/envs/captus

captus_assembly assemble -r fastp_polyg -t 48

extract -a 02_assemblies -n Angiosperms353

captus align -e 03_extractions
```

# Filter alignments

```bash
# python script to remove taxa represented for < 50 loci
# fasta_filter_debug.py
# make executable sudo chmod 777 
# usage: python3 fasta_filter_debug.py

import os
from collections import defaultdict

# Define input and output directories
input_dir = "/Data/cbarrett/004_2024_MH_EpiBase_seqcap/05_NR_5x_aligned/03_trimmed/05_naive/01_coding_NUC/02_NT"
output_dir = "/Data/cbarrett/004_2024_MH_EpiBase_seqcap/05_NR_5x_aligned_filtered"
os.makedirs(output_dir, exist_ok=True)

### # List of specific taxa to remove -- this part isn't working properly!
### taxa_to_remove = [
###     "fastp-US34_Aphyllorchis_montana_merged",
###     "fastp-US30_Uleiorchis_ulei_merged",
###     "fastp-Pseudovanilla_foliata",
###     "fastp-US31_Stereosandra_javanica_merged",
###     "fastp-US63_Elleanthus_aurantiacus_S46",
###     "fastp-Microchilus_plantagineus",
###     "fastp-Coelogyne_porrecta",
###     "fastp-Erythrorchis_altissima",
###     "fastp-Dactylostalix-ringensoides_S41",
###     "fastp-Yoania-prainii_S14",
###     "fastp-Yoania-japonica_S12",
###     "fastp-Mesadenus_lucayanus",
###     "fastp-Coelogyne_entomophobia"
### ]

# Print all taxa in the taxa_to_remove list for verification
print("\nList of taxa to remove:")
for taxon in taxa_to_remove:
    print(taxon)

# Function to parse fasta headers and remove extra text
def clean_header(header):
    return header.split()[0].strip()  # Strip any extra spaces around the taxon name

# Initialize a dictionary to count taxa occurrences
taxa_counts = defaultdict(int)

# Step 1: Clean headers and count taxa occurrences across all files
print("\nProcessing files and counting taxa...")
for filename in os.listdir(input_dir):
    if filename.endswith(".fna"):
        with open(os.path.join(input_dir, filename), 'r') as infile:
            output_path = os.path.join(output_dir, filename)
            with open(output_path, 'w') as outfile:
                for line in infile:
                    if line.startswith(">"):
                        taxon = clean_header(line.strip())
                        taxa_counts[taxon] += 1
                        outfile.write(f"{line}")
                    else:
                        outfile.write(line)

# Step 2: Filter taxa that occur in 50 or more alignments
taxa_to_keep = {taxon for taxon, count in taxa_counts.items() if count >= 50}

if not taxa_to_keep:
    print("\nNo taxa found in 50 or more alignments. Exiting.")
    exit(1)

# Print out all taxa that will be retained after filtering
print("\nTaxa to be kept (occurring in 50 or more alignments):")
for taxon in taxa_to_keep:
    print(taxon)

# Step 3: Remove specific taxa and keep only the filtered ones
print("\nFiltering taxa and saving output...")
for filename in os.listdir(output_dir):
    if filename.endswith(".fna"):
        temp_lines = []
        with open(os.path.join(output_dir, filename), 'r') as infile:
            lines = infile.readlines()
        
        write_taxon = False
        for line in lines:
            if line.startswith(">"):
                taxon = clean_header(line.strip())
                print(f"\nProcessing taxon: {taxon}")
                
                # Check if the taxon should be removed
                if taxon in taxa_to_keep and taxon not in taxa_to_remove:
                    print(f"Retaining taxon: {taxon}")
                    write_taxon = True
                else:
                    print(f"Removing taxon: {taxon}")
                    write_taxon = False
            if write_taxon:
                temp_lines.append(line)
        
        # Overwrite the file with filtered content
        with open(os.path.join(output_dir, filename), 'w') as outfile:
            outfile.writelines(temp_lines)

# Step 4: Output the list of removed and retained taxa
removed_taxa = [taxon for taxon in taxa_counts if taxon not in taxa_to_keep or taxon in taxa_to_remove]
retained_taxa = [taxon for taxon in taxa_to_keep if taxon not in taxa_to_remove]

print("\nTaxa removed:")
for taxon in removed_taxa:
    print(taxon)

print("\nTaxa retained:")
for taxon in retained_taxa:
    print(taxon)

print(f"\nFiltered output has been saved to {output_dir}")

# Create a list of taxa to further remove from all alignments:

nano taxa_to_remove.txt

fastp-US34_Aphyllorchis_montana_merged
fastp-US30_Uleiorchis_ulei_merged
fastp-Pseudovanilla_foliata
fastp-US31_Stereosandra_javanica_merged
fastp-US63_Elleanthus_aurantiacus_S46
fastp-Microchilus_plantagineus
fastp-Coelogyne_porrecta
fastp-Erythrorchis_altissima
fastp-Dactylostalix-ringensoides_S41
fastp-Yoania-prainii_S14
fastp-Yoania-japonica_S12
fastp-Mesadenus_lucayanus
fastp-Coelogyne_entomophobia
fastp-Calypso-bulbosa-Asia_S37
fastp-Calypso-bulbosa-var-americana_S4
fastp-Coelogyne_fimbriata
fastp-Coelogyne_sulphurea
fastp-Corallorhiza-bulbosa_S19
fastp-Corallorhiza_ekmanii
fastp-Corallorhiza-involuta_S1
fastp-Corallorhiza-macrantha_S21
fastp-Corallorhiza-maculata-var-maculata2_S36
fastp-Corallorhiza-maculata-var-mexicana_S17
fastp-Corallorhiza-mertensiana_S30
fastp-Corallorhiza-odontorhiza-Mexico_S47
fastp-Corallorhiza-odontorhiza-var-odontorhiza2_S38
fastp-Corallorhiza-striata-var-vreelandii_S7
fastp-Corallorhiza-striata-Sierra-Nevada_S9
fastp-Corallorhiza-wisteriana-Western-US_S28
fastp-Cremastra-saprophytica_S22
fastp-Cremastra-unguiculata_S43
fastp-Cremastra-variabilis_S24
fastp-Dendrobium_draconis
fastp-Dendrobium_ellipsophyllum
fastp-Dendrobium_toressae
fastp-Ephippianthus-sawadanus_S31
fastp-Govenia-superba_S6
fastp-Govenia_utriculata
fastp-Habenaria_arenaria
fastp-Maxillaria_platypetala
fastp-Oreorchis-coreana_S32
fastp-Oreorchis-erythrochrysea_S33
fastp-Oreorchis-fargesii_S34
fastp-Oreorchis-indica1_S35
fastp-Oreorchis-indica2_S46
fastp-Oreorchis-patens_S29
fastp-Pelexia_pterygantha
fastp-Vanilla_aff._odorata
fastp-Vanilla_dressleri
fastp-Vanilla_hartii


# Get taxon names from a newick file
grep -oP '[\w\.\-]+' bestmix_min100.treefile | sort | uniq > taxon_names.txt


# ctrl + X, Y, enter

# Remove taxa from alignments with seqkit

conda activate /usr/local/src/conda_envs/binf

mkdir filtered

for file in *.fna.fasta; do
    seqkit grep -v -f taxa_to_remove.txt "$file" -o "$file.filtered"
done
```bash



