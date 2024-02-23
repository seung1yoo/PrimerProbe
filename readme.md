# Introduction
This script is designed to identify the primer set utilized for microbiome data obtained from the Sequence Read Archive (SRA). The initial step in microbiome data analysis involves removing primer sequences. These primers are typically used to amplify desired regions of the 16S rRNA, and various primers may be employed for each experiment. To analyze data from multiple groups or institutions, it's essential to ascertain which primers were used for each dataset. Therefore, we coded a simple
script for this purpose.

# Basic approch
The basic approach involves counting the occurrences of known primer sequences in the top 1000 reads of the FASTQ files and returning the most frequently observed primer set. However, a more advanced approach would entail performing a BLAST search against a 16S sequence database with the top 1000 sequences from the FASTQ files to estimate the target region. This idea is planned for implementation in the future.
