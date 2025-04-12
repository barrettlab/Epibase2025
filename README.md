### All analyses were run in Linux or R
```bash
Torvalds L, Hamano J. 2024. The Linux kernel. Version 6.8. Available from: https://kernel.org
R Core Team. 2025. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.​
ChaptGPT and Github Copilot were used for debugging and code development in Inix and R:
OpenAI. (2025). ChatGPT 4o (April 8 version) [Large language model]. https://chat.openai.com
```

[Github Copilot](https://github.com/features/copilot?ef_id=_k_CjwKCAjwktO_BhBrEiwAV70jXs7vOsmIA9Wy5KM0DqvmR0kMfaYvsrDk3mYBan4xL4QmP4eD8pSiAhoCgbQQAvD_BwE_k_&OCID=AIDcmmb150vbv1_SEM__k_CjwKCAjwktO_BhBrEiwAV70jXs7vOsmIA9Wy5KM0DqvmR0kMfaYvsrDk3mYBan4xL4QmP4eD8pSiAhoCgbQQAvD_BwE_k_&gad_source=1&gclid=CjwKCAjwktO_BhBrEiwAV70jXs7vOsmIA9Wy5KM0DqvmR0kMfaYvsrDk3mYBan4xL4QmP4eD8pSiAhoCgbQQAvD_BwE)


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
```

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

# Clean up alignments with TrimAL
```bash
parallel trimal -in {} -out trimal_min100{.}.fasta -htmlout trimal_min100/{.}_trim.html -automated1 ::: *.fna.fasta
```

# Concatenated data
```bash
# Pull all files into geneious and keep only those with > 50 taxa
# Concatenate sequences in Geneious, and migrate alignments back to server

# Run iqtree2 under various models on concatenated data

# GTR+gamma
iqtree2 -s min100_concat.fasta -o fastp-Apostasia_nuda -m GTR+G -T 92 -B 1000 --prefix gtrg_min100

# MixtureFinder 
nohup iqtree2 -s min100_concat.fasta -o fastp-Apostasia_nuda -m MIX+MF -t BIONJ -T 92 --prefix bestmix_min100 &
iqtree2 -s min100_concat.fasta -o fastp-Apostasia_nuda -m "(MIX-XXXX}} -B 1000 -T 92 --prefix bestmix_min100_boot

# Ghost with H2,4,6 rate mixture classes (heterotachy), use BIC to choose among them
nohup iqtree2 -s filtered_prank_full_concat.fasta -o fastp-Apostasia_nuda -m GTR+FO*H2 -T 92 --prefix H2_min100 &
nohup iqtree2 -s filtered_prank_full_concat.fasta -o fastp-Apostasia_nuda -m GTR+FO*H4 -T 92 --prefix H4_min100 &
nohup iqtree2 -s filtered_prank_full_concat.fasta -o fastp-Apostasia_nuda -m GTR+FO*H6 -T 92 --prefix H6_min100 &

# Choose the best GHOST model and run bootstrapping
iqtree2 -s min100_concat.fasta -o fastp-Apostasia_nuda -m GTR+FO*HXXXX -T 92 -B 1000 --prefix HXXXX_min100_boot
```

# PRANK

```bash
# Run prank on all filtered alignments

#!/bin/bash

# Create the output directory if it doesn't exist
mkdir -p prank

# Find all .filtered files in the "filtered" directory and process them in parallel
find filtered/ -name "*.filtered" | parallel -j 80 ' \
    base=$(basename {} .filtered); \
    prank -d={} -o=prank/${base}_prank -F'
```


# Astral polytomy test

```bash
java -jar /usr/local/src/ASTRAL/astral.5.7.8.jar -i gts.tre -q st.tre -o output.tre --branch-annotate 10

## gts.tre = gene trees
## st.tre = wastral species tree
## --branch-annotate 10 = specify p-values for polytomy test, where p < 0.05 = reject polytomy (Ho) and p > 0.95 = fail to reject Ho

## plot the output in R, with red for significant polytomy and blue rejecting polytomy

# Load necessary libraries
library(ggtree)
library(treeio)
library(ggplot2)

# Load the tree with branch annotations
tree_file <- "output.tre"
tree <- read.tree(tree_file)

# Extract branch annotations (e.g., polytomy test p-values)
# Assuming branch annotations are stored as node labels
p_values <- as.numeric(tree$node.label)

# Define colors based on p-value thresholds
branch_colors <- ifelse(p_values < 0.05, "blue", 
                        ifelse(p_values > 0.95, "red", "gray"))

# Create a data frame for ggtree visualization
branch_data <- data.frame(node = (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode),
                          p_value = p_values,
                          color = branch_colors)

# Plot the tree with ggtree
p <- ggtree(tree) %<+% branch_data +
  geom_text2(aes(subset = !is.na(p_value), label = round(p_value, 3)),
             vjust = -0.5, hjust = 0.5, nudge_x = -0.02, size = 3, color = "black") +  # Add p-value labels and nudge them left
  geom_tiplab(size = 3) +  # Show tip labels
  geom_tree(aes(color = color)) +  # Color branches
  scale_color_identity() +  # Use predefined colors
  theme_tree2()

# Display the plot
print(p)

# Save the plot
ggsave("polytomy_test_tree.pdf", plot = p, width = 10, height = 20)
```

# MixtureFinder loop over all alignments
```bash
### Run MixtureFinder for 1-3 mixture classes as a loop over all alignments, choose best-fit mixture model

# Define directories
input_dir="prank"  # Directory with alignment files
output_dir="mix_loop_out"     # Directory to save IQ-TREE outputs
mkdir -p "$output_dir"            # Create output directory if it doesn't exist

# Initialize a results file to store the best model and BIC scores
bic_file="$output_dir/best_model_bic_summary.txt"
echo -e "Alignment\tBest_Model\tBIC_Score" > "$bic_file"  # Header for BIC scores summary

# Function to run ModelFinder for each alignment and capture BIC score
run_modelfinder_bic() {
    alignment="$1"
    base=$(basename "$alignment" .fas)
    
	# run for 1-4 mixture classes
    iqtree2 -s "$alignment" -m MIX+MFP -cmax 4 -B 1000 --prefix "$output_dir/${base}_model_selection" -T 2 > /dev/null 2>&1

	best_model=$(grep "Best-fit model according to BIC:" "$output_dir/${base}_model_selection.iqtree" | awk '{print $6}')
    bic_score=$(grep "BIC score:" "$output_dir/${base}_model_selection.iqtree" | awk '{print $3}')
    echo -e "${base}\t${best_model}\t${bic_score}" >> "$bic_file"
}

# Export functions and variables for GNU parallel
export -f run_modelfinder_bic
export output_dir model_set bic_file

# Step #2: Run ModelFinder in parallel for each alignment
echo "Starting ModelFinder for each alignment..."
find "$input_dir" -name "*.fas" | parallel -j 45 run_modelfinder_bic

# Wait for all ModelFinder processes to complete
wait  # Ensures all parallel jobs are completed before moving to the next steps
echo "ModelFinder analysis complete for all alignments. Summary of best models and BIC scores saved to $bic_file."

### Sum the BIC scores and # parameters across alignments
### Pull values and create a table: alignment, NP, BIC, logL


#!/bin/bash

# Define output file for the results table
output_file="free_parameters_bic_logl_summary.tsv"
echo -e "Filename\tFree_Parameters\tLog_Likelihood\tBIC_Score" > "$output_file"

# Loop through each .iqtree file in the directory
for file in mix_loop_out/*.iqtree; do
  echo "Processing file: $file"  # Debugging output to track file processing

  # Initialize variables to store extracted data
  free_params="NA"
  log_likelihood="NA"
  bic_score="NA"

  # Parse the .iqtree file line by line
  while IFS= read -r line; do
    # Extract the number of free parameters
    if [[ "$line" =~ "Number of free parameters" ]]; then
      free_params=$(echo "$line" | awk '{print $NF}')
    # Extract the log-likelihood (capture the negative or positive number only)
    elif [[ "$line" =~ "Log-likelihood of the tree:" ]]; then
      log_likelihood=$(echo "$line" | awk '{for(i=1;i<=NF;i++) if ($i ~ /^-?[0-9]+\.[0-9]+$/) print $i}')
    # Extract the BIC score
    elif [[ "$line" =~ "Bayesian information criterion (BIC) score:" ]]; then
      bic_score=$(echo "$line" | awk '{print $NF}')
    fi
  done < "$file"

  # Check if any fields are still NA, indicating a parsing failure
  if [[ "$free_params" == "NA" || "$log_likelihood" == "NA" || "$bic_score" == "NA" ]]; then
    echo "Warning: Missing data for $file. Manual review may be needed."  # Debugging output
  fi

  # Append extracted data to the output file
  echo -e "$(basename "$file")\t$free_params\t$log_likelihood\t$bic_score" >> "$output_file"
done
```

# Plotting QC, QD, QI, Branch Lengths (Quartet Sampling); run GLM and GAM models

```bash

setwd("E:/0001_research/0001_manuscripts/2024_A353_Epi_Base_MH/2024_10_15_filtered_analyses")

library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)

library(ape)
library(dplyr)
library(ggplot2)

# Load the Newick file containing the four trees
trees <- read.tree("annotated_qs.tre")

# Extract the trees in order
phylo_tree <- trees[[1]]  # Phylogram with branch lengths
qc_tree <- trees[[2]]     # QC tree
qd_tree <- trees[[3]]     # QD tree
qi_tree <- trees[[4]]     # QI tree

# Helper function to extract internal node metadata
extract_internal_metadata <- function(tree, label_prefix) {
  # Get internal node indices
  internal_nodes <- (Ntip(tree) + 1):max(tree$edge[, 2])
  
  # Extract metadata from node labels
  metadata <- sapply(tree$node.label, function(label) {
    if (!is.na(label) && grepl(paste0(label_prefix, "="), label)) {
      as.numeric(sub(paste0(".*", label_prefix, "="), "", label))
    } else {
      NA
    }
  })
  
  # Return metadata corresponding to internal nodes
  metadata[internal_nodes - Ntip(tree)]
}

# Extract internal branch lengths from the phylogram
internal_nodes <- (Ntip(phylo_tree) + 1):max(phylo_tree$edge[, 2])
branch_lengths <- phylo_tree$edge.length[match(internal_nodes, phylo_tree$edge[, 2])]

# Extract QC, QD, and QI values
qc_values <- extract_internal_metadata(qc_tree, "qc")
qd_values <- extract_internal_metadata(qd_tree, "qd")
qi_values <- extract_internal_metadata(qi_tree, "qi")

# Combine into a data frame
data <- data.frame(
  BranchLength = branch_lengths,
  QC = qc_values,
  QD = qd_values,
  QI = qi_values
)

# Print the data frame to check
print(data)

#######################

# Plotting QC vs. QI, highlighting EDE monophyly and relationship nodes in plot

library(ggplot2)
library(dplyr)
library(ggrepel)

# Define special node groups
red_nodes <- c(208:214, 221:222, 316)
blue_nodes <- c(215, 305, 307, 317, 320, 325, 328)

# Prepare the data
data <- data %>%
  mutate(
    qi = as.numeric(as.character(qi)),
    qc = as.numeric(as.character(qc)),
    group = case_when(
      node %in% red_nodes ~ "red_square",
      node %in% blue_nodes ~ "blue_triangle",
      TRUE ~ "gray_circle"
    )
  )

# Plot
ggplot(data, aes(x = qi, y = qc)) +
  # LOESS smoothing curve with CI ribbon
  geom_smooth(method = "loess", color = "black", fill = "lightgray",
              linewidth = 1, alpha = 0.4, se = TRUE) +
  
  # Plot points with different shapes/fills
  geom_point(aes(shape = group, fill = group),
             size = 3.5, stroke = 0.3, color = "black", alpha = 0.9) +
  
  # Add non-overlapping node labels for red/blue points
  geom_text_repel(data = data %>% filter(group != "gray_circle"),
                  aes(label = node),
                  size = 3.5,
                  fontface = "bold",
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  segment.color = "gray50",
                  segment.size = 0.3,
                  box.padding = 0.3,
                  point.padding = 0.2) +
  
  # Manual shape and fill definitions
  scale_shape_manual(values = c(
    gray_circle = 21,
    red_square = 22,
    blue_triangle = 24
  )) +
  scale_fill_manual(values = c(
    gray_circle = "gray60",
    red_square = "red",
    blue_triangle = "blue"
  )) +
  
  # Axis styling
  scale_x_continuous(name = "qi", breaks = scales::pretty_breaks(n = 6),
                     expand = expansion(mult = 0.05)) +
  scale_y_continuous(name = "qc", breaks = scales::pretty_breaks(n = 6),
                     expand = expansion(mult = 0.05)) +

  # Theme tweaks
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "none"
  )

##################################

# Box/jitter plot of BL x EDE/non-EDE

# Combine EDE categories into one
data_long <- data_long %>%
    mutate(HighlightGroupCombined = case_when(
        HighlightGroup %in% c("EDE Tribe Relationships", "EDE Tribe Monophyly") ~ "EDE Tribes",
        TRUE ~ "Other"
    ))

# Plot box/whisker/jitter plot with horizontal orientation
ggplot(data_long, aes(x = HighlightGroupCombined, y = BranchLength, fill = HighlightGroupCombined)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Boxplot without outliers
    geom_jitter(aes(color = HighlightGroup), size = 3, width = 0.2) + # Jitter points
    scale_fill_manual(values = c("EDE Tribes" = "blue", "Other" = "orange")) +
    scale_color_manual(values = c(
        "EDE Tribe Relationships" = "red",
        "EDE Tribe Monophyly" = "blue",
        "Other" = "orange"
    )) +
    labs(title = "Branch Lengths for Highlighted Groups",
         y = "Group",  # Switched x and y labels
         x = "Branch Length") +
    theme_minimal() +
    theme(legend.position = "right") + # Move legend for clarity
    coord_flip()  # Force horizontal orientation


##################################################



## Generalized linear model incorporating BranchLength interactions w QD and QC
# Convert HighlightGroup to a dummy variable

data_long <- data_long %>%
    mutate(EDE_Group = ifelse(HighlightGroup %in% c("EDE Tribe Monophyly", "EDE Tribe Relationships"), 1, 0))

data_long_wide <- data_long %>%
    pivot_wider(names_from = Metric, values_from = Value)



# Fit the GLM
glm_model <- glm(
    QC ~ QI * BranchLength + QD * BranchLength + EDE_Group,
    data = data_long_wide,
    family = gaussian()  # Assuming QC is continuous and normally distributed
)

# Summarize the model
summary(glm_model)

Call:
glm(formula = QC ~ QI * BranchLength + QD * BranchLength + EDE_Group, 
    family = gaussian(), data = data_long_wide)

Deviance Residuals: 
     Min        1Q    Median        3Q       Max  
-1.45390  -0.05805   0.10634   0.16496   0.72616  

Coefficients:
                Estimate Std. Error t value Pr(>|t|)    
(Intercept)     -2.50189    0.39810  -6.285 2.30e-09 ***
QI               3.33154    0.41358   8.055 9.46e-14 ***
BranchLength    -2.86131   16.85842  -0.170    0.865    
QD              -0.07701    0.14376  -0.536    0.593    
EDE_Group        0.05979    0.08981   0.666    0.506    
QI:BranchLength  3.79785   17.75292   0.214    0.831    
BranchLength:QD  3.27343    9.00195   0.364    0.717    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for gaussian family taken to be 0.1152682)

    Null deviance: 42.179  on 191  degrees of freedom
Residual deviance: 21.325  on 185  degrees of freedom
  (2 observations deleted due to missingness)
AIC: 138.93

Number of Fisher Scoring iterations: 2

## GAM (Generalized additive model, non-linear relationships)

library(mgcv)

# Fit the GAM model
qc_glm <- glm(QC ~ BranchLength * QD * QI, data = data_long_wide)
qc_gam <- gam(QC ~ s(BranchLength) + s(QD) + s(QI) + ti(BranchLength,QD,QI), data = data_long_wide)
# qc_glm_ede <- glm(QC ~ EDE_Group * BranchLength * QD * QI, data = data_long_wide)
# qc_gam_ede <- gam(QC ~ EDE_Group + s(BranchLength) + s(QD) + s(QI) + ti(EDE_Group,BranchLength,QD,QI), data = data_long_wide)

# Summarize the model
summary(qc_glm)
summary(qc_gam)
# summary(qc_glm_ede)
# summary(qc_gam_ede)

AIC(qc_glm,qc_gam)

Call:
glm(formula = QC ~ BranchLength * QD * QI, data = data_long_wide)

Deviance Residuals: 
     Min        1Q    Median        3Q       Max  
-1.43034  -0.05959   0.10066   0.12406   0.74344  

Coefficients:
                    Estimate Std. Error t value Pr(>|t|)    
(Intercept)          -3.6022     0.5012  -7.188 1.59e-11 ***
BranchLength         -7.5942    19.8615  -0.382   0.7026    
QD                    2.1307     1.1229   1.898   0.0593 .  
QI                    4.4770     0.5197   8.614 3.14e-15 ***
BranchLength:QD     159.0021    83.7190   1.899   0.0591 .  
BranchLength:QI       8.8413    20.8768   0.423   0.6724    
QD:QI                -2.4444     1.2457  -1.962   0.0512 .  
BranchLength:QD:QI -168.7619    92.4982  -1.824   0.0697 .  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for gaussian family taken to be 0.104247)

    Null deviance: 42.179  on 191  degrees of freedom
Residual deviance: 19.181  on 184  degrees of freedom
  (2 observations deleted due to missingness)
AIC: 120.59

Number of Fisher Scoring iterations: 2

> summary(qc_gam)

Family: gaussian 
Link function: identity 

Formula:
QC ~ s(BranchLength) + s(QD) + s(QI) + ti(BranchLength, QD, QI)

Parametric coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)  0.61349    0.02008   30.56   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Approximate significance of smooth terms:
                         edf Ref.df      F  p-value    
s(BranchLength)        1.000  1.000  0.054   0.8162    
s(QD)                  6.543  7.636  5.509 6.05e-06 ***
s(QI)                  4.821  5.871 13.726  < 2e-16 ***
ti(BranchLength,QD,QI) 1.000  1.000  4.756   0.0305 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-sq.(adj) =  0.665   Deviance explained = 68.9%
GCV = 0.079883  Scale est. = 0.073906  n = 192
> # summary(qc_glm_ede)
> # summary(qc_gam_ede)
> 
> AIC(qc_glm,qc_gam)
             df       AIC
qc_glm  9.00000 120.59037
qc_gam 15.36395  60.51905
 
 # Fit the GLM/GAM models
 qc_glm <- glm(QC ~ BranchLength * QD * QI, data = data_long_wide)
 qc_gam <- gam(QC ~ s(BranchLength) + s(QD) + s(QI) + ti(BranchLength,QI) + ti(BranchLength,QD) + ti(BranchLength,QD,QI), data = data_long_wide)
 summary(qc_glm)
 summary(qc_gam)
 AIC(qc_glm,qc_gam)

 # Summarize the model
 summary(qc_glm)

Call:
glm(formula = QC ~ BranchLength * QD * QI, data = data_long_wide)

Deviance Residuals: 
     Min        1Q    Median        3Q       Max  
-1.43034  -0.05959   0.10066   0.12406   0.74344  

Coefficients:
                    Estimate Std. Error t value Pr(>|t|)    
(Intercept)          -3.6022     0.5012  -7.188 1.59e-11 ***
BranchLength         -7.5942    19.8615  -0.382   0.7026    
QD                    2.1307     1.1229   1.898   0.0593 .  
QI                    4.4770     0.5197   8.614 3.14e-15 ***
BranchLength:QD     159.0021    83.7190   1.899   0.0591 .  
BranchLength:QI       8.8413    20.8768   0.423   0.6724    
QD:QI                -2.4444     1.2457  -1.962   0.0512 .  
BranchLength:QD:QI -168.7619    92.4982  -1.824   0.0697 .  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for gaussian family taken to be 0.104247)

    Null deviance: 42.179  on 191  degrees of freedom
Residual deviance: 19.181  on 184  degrees of freedom
  (2 observations deleted due to missingness)
AIC: 120.59

Number of Fisher Scoring iterations: 2

 summary(qc_gam)

Family: gaussian 
Link function: identity 

Formula:
QC ~ s(BranchLength) + s(QD) + s(QI) + ti(BranchLength, QD, QI)

Parametric coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)  0.61349    0.02008   30.56   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Approximate significance of smooth terms:
                         edf Ref.df      F  p-value    
s(BranchLength)        1.000  1.000  0.054   0.8162    
s(QD)                  6.543  7.636  5.509 6.05e-06 ***
s(QI)                  4.821  5.871 13.726  < 2e-16 ***
ti(BranchLength,QD,QI) 1.000  1.000  4.756   0.0305 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-sq.(adj) =  0.665   Deviance explained = 68.9%
GCV = 0.079883  Scale est. = 0.073906  n = 192

 
 AIC(qc_glm,qc_gam)
             df       AIC
qc_glm  9.00000 120.59037
qc_gam 15.36395  60.51905
```

# Plotting Quartet Scores on Tree

```bash
### Load libraries

if (!require(ape)) install.packages("ape", dependencies = TRUE); library(ape)
if (!require(phytools)) install.packages("phytools", dependencies = TRUE); library(phytools)
if (!require(gtools)) install.packages("gtools", dependencies = TRUE); library(gtools)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")
if (!require(cowplot)) install.packages("cowplot", dependencies = TRUE); library(cowplot)
if (!require(ggimage)) install.packages("ggimage", dependencies = TRUE); library(ggimage)
if (!require(treeio)) install.packages("treeio", dependencies = TRUE); library(treeio)
if (!require(gridExtra)) install.packages("gridExtra", dependencies = TRUE); library(gridExtra)

# Or if already installed
library(ape)
library(phytools)
library(gtools)
library(ggtree)
library(cowplot)
library(ggimage)
library(treeio)
library(gridExtra)
library(dplyr)

### Read in ASTRAL formatted tree with annotations (-u 2 option in ASTRAL for verbose output)
# AS_0all <- read.astral("prank_filtered_GT_verbose_rerooted.tre")
AS_0all <- read.astral("mixall_pranks_verbose_rooted.tre")
p_0all <- ggtree(AS_0all@phylo, ladderize=T, branch.length = "none") + geom_tiplab(size=3.5, hjust= -0.05) + xlim_tree(70)
write(AS_0all@phylo$tip.label, file= "CDS_alM_r_g_c0all_trees_BP10_SpeciesTree_annotQ_rooted2_tips.txt")

### Build a data frame for plotting q-values (3)
### Cite Sidonie's R code: https://github.com/sidonieB/scripts/blob/master/plot_Astral_trees_v2.R

QS_calc <- function(AS_A) {
  Q1 <- as.numeric(AS_A@data$q1) * 100
  Q <- as.data.frame(Q1)
  Q$Q2 <- as.numeric(AS_A@data$q2) * 100
  Q$Q3 <- as.numeric(AS_A@data$q3) * 100
  Q$node <- AS_A@data$node
  return(Q)
}


# Function to build a data frame for plotting f-values plus uninformative GTs (like phyparts)
QS2_calc <- function(AS_A) {
  Q1 <- as.numeric(AS_A@data$f1 / 318)
  Q <- as.data.frame(Q1)
  Q$Q2 <- as.numeric(AS_A@data$f2 / 318)
  Q$Q3 <- as.numeric(AS_A@data$f3 / 318)
  Q$Q4 <- (1 - (Q$Q1 + Q$Q2 + Q$Q3)) 
  Q$node <- AS_A@data$node
  return(Q)
}


# Run the function
Q_0all <- QS2_calc(AS_0all)

# Build node pies
pies_0all <- nodepie(Q_0all, cols=1:4, color=c(Q1='blue', Q2='red', Q3='orange', Q4='gray')) #, alpha=.6) # change style as needed, alpha is for transparency

# Create ggtree object with pie charts
p2_0all <- inset(p_0all, pies_0all, width=0.03, height=0.03)#, hjust=-.6) # change size if pies too big/small (better check first how it looks in the exported file)
# plot figure
p2_0all


### Add f1 (concordant) and f2+3 (discordant) gene tree numbers per node
p3 <- p2_0all %<+% Q_0all # adds a data frame to a ggtree object (operator = %<+%)

p3$data$Q1 # includes NA values for tips!

### Convert proportions to gene trees, rounded to integer

# Concordant GTs
p3$data <- p3$data %>%
	mutate(Q1_transformed = round(Q1 * 318))

# Sum of two other discordant quartet in gene trees
p3$data <- p3$data %>%
	mutate(Q23_transformed = round((Q2 + Q3) * 318))

### PLOT! May have to adjust sizes and "nudge" parameters

p4 <- p3 + 
    geom_nodelab(aes(label = Q1_transformed), 
                 nudge_x = -0.5, 
                 nudge_y = 0.25, 
                 size = 3)  +
    
    geom_nodelab(aes(label = Q23_transformed), 
                 nudge_x = -0.5,
                 nudge_y = -0.25,
                 size = 3)
p4
```

# Topology tests

## Searches based on 3 constraints: wastral (coalescent), zhang2023, and PE2024

```bash 
iqtree2 -s filtered_prank_concat.phy -m GTR+F+R7 -g st_nolpp_nobl.tre --prefix wastral_constraint -T 45 ; 
iqtree2 -s filtered_prank_concat.phy -m GTR+F+R7 -g zhang23_constraint_tree.tre --prefix zhang23_constraint -T 45 ; 
iqtree2 -s filtered_prank_concat.phy -m GTR+F+R7 -g pe24_constraint_tree.tre --prefix pe24_constraint -T 45 ;
```

## Concatenate output trees
```bash
cat unconstrained_concat.treefile wastral_constraint.treefile zhang23_constraint.treefile pe24_constraint.treefile > allconstraints.tre
```

## Run tests
```bash
 iqtree2 -s filtered_prank_concat.phy -m GTR+F+R7 -z allconstraints.tre -n 0 -zb 1000 -au -T 55 --prefix allconstraints

See allconstraints.trees for trees with branch lengths.

Tree      logL    deltaL  bp-RELL    p-KH     p-SH       c-ELW       p-AU
-------------------------------------------------------------------------
  1 -1736917.597       0   0.901 +  0.914 +      1 +       0.9 +    0.955 + 
  2 -1736999.485  81.888  0.0009 - 0.00145 - 0.00595 -  0.000882 -  0.00381 - 
  3 -1736985.653  68.056  0.0155 - 0.0171 - 0.0292 -    0.0156 -    0.024 - 
  4 -1736937.312  19.715  0.0826 + 0.0857 +  0.418 +    0.0835 +   0.0942 + 

deltaL  : logL difference from the maximal logl in the set.
bp-RELL : bootstrap proportion using RELL method (Kishino et al. 1990).
p-KH    : p-value of one sided Kishino-Hasegawa test (1989).
p-SH    : p-value of Shimodaira-Hasegawa test (2000).
c-ELW   : Expected Likelihood Weight (Strimmer & Rambaut 2002).
p-AU    : p-value of approximately unbiased (AU) test (Shimodaira, 2002).

1 unconstrained concatenated
2 wastral A353 coalescent
3 Zhang 2023 transcriptomes coalescent
4 Perez-Escobar A353 coalescent 

```

## Use quartet scores for constraint trees/gene trees in ASTRAL v.5.7.8

```bash
# Quartet scores among constraint trees in astral v.5.7.8
mkdir test_data

java -jar /usr/local/src/ASTRAL/astral.5.7.8.jar -i gts.tre -q unconstrained_concat.treefile  -o test_data/unconstrained-scored.tre  2> test_data/unconstrained-scored.log
java -jar /usr/local/src/ASTRAL/astral.5.7.8.jar -i gts.tre -q wastral_constraint.treefile    -o test_data/wastral-scored.tre  2> test_data/wastral-scored.log
java -jar /usr/local/src/ASTRAL/astral.5.7.8.jar -i gts.tre -q zhang23_constraint.treefile    -o test_data/zhang23-scored.tre  2> test_data/zhang23-scored.log
java -jar /usr/local/src/ASTRAL/astral.5.7.8.jar -i gts.tre -q pe24_constraint.treefile       -o test_data/pe24-scored.tre  2> test_data/pe24-scored.log
```

## Calculate RF distances among all trees
```bash
iqtree2 -rf_all allconstraints.tre --prefix rf_distances

###   4 4
###   Tree0       0 8 8 4		# unconstrained concatenated
###   Tree1       8 0 10 10		# wastral A353 coalescent
###   Tree2       8 10 0 8		# Zhang 2023 transcriptomes coalescent 
###   Tree3       4 10 8 0		# Perez-Escobar A353 coalescent
```

## Plot R-F distances using MDS in R

```{r}
# Load necessary libraries
library(ggplot2)

# Define the Robinson-Foulds distance matrix
rf_matrix <- matrix(c(
  0,  8,  8,  4,
  8,  0, 10, 10,
  8, 10,  0,  8,
  4, 10,  8,  0
), nrow=4, byrow=TRUE)

# Define tree names
tree_names <- c("Unconstrained Concatenated", "wASTRAL A353 Coalescent", 
                "Zhang 2023 Transcriptomes Coalescent", "Perez-Escobar A353 Coalescent")

# Perform classical MDS (Multidimensional Scaling)
mds <- cmdscale(rf_matrix, k=2)  # k=2 for two-dimensional scaling

# Convert to dataframe for plotting
mds_df <- as.data.frame(mds)
mds_df$Tree <- tree_names  # Add tree names as a column

# Create MDS plot using ggplot2
ggplot(mds_df, aes(x=V1, y=V2, label=Tree)) +
  geom_point(size=4, color="blue") +  # Add points
  geom_text(vjust=-0.5, hjust=0.5, size=4) +  # Add tree labels
  theme_minimal() +
  labs(title="MDS Plot of Robinson-Foulds Distances",
       x="MDS Dimension 1", 
       y="MDS Dimension 2")
```
## Plot gene tree representation on species tree estimate in R

```{r}
# Load required libraries
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Read the species tree
species_tree <- read.tree("st.tre")

# Read the gene trees
gene_trees <- read.tree("gts.tre")

# Count how many times each taxon appears across all gene trees
taxon_counts <- table(unlist(lapply(gene_trees, function(tree) tree$tip.label)))

# Convert taxon counts into a dataframe for plotting
taxon_data <- data.frame(
  taxon = names(taxon_counts),
  count = as.numeric(taxon_counts)
)

# Ensure taxa in species tree are included, filling missing ones with 0
taxon_data <- taxon_data %>%
  right_join(data.frame(taxon = species_tree$tip.label), by = "taxon") %>%
  tidyr::replace_na(list(count = 0))

# Plot species tree with extra space for bars
p_tree <- ggtree(species_tree) + 
  geom_tiplab(size = 3, align = TRUE, linetype = "dotted") +  # Align labels with dotted guide lines
  xlim(0, max(taxon_data$count) + 5)  # Extend x-axis to make space for bars

# Extract tree data to get y-axis positions for taxa
tree_data <- p_tree$data

# Merge taxon count data with tree node positions
taxon_data$y <- match(taxon_data$taxon, tree_data$label)  # Match taxon names to tree y-positions

# Define an offset value to shift the bars right
offset <- max(tree_data$x) + 2  # Shift bars to the right

# Create bar plot of missing data (shifted right)
p_bars <- ggplot(taxon_data, aes(x = count + offset, y = factor(y, levels = rev(y)))) +  
  geom_bar(stat = "identity", fill = "blue") +
  labs(x = "Gene Tree Count", y = NULL) +  
  theme_minimal() +
  theme(axis.text.y = element_blank(),  # Hide y-axis labels (tree already has them)
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  xlim(offset, offset + max(taxon_data$count) + 2)  # Ensure bars fit in allocated space

# Combine tree and bar plots side-by-side
p_tree + p_bars + plot_layout(widths = c(2, 1))  # Adjust layout widths
```
