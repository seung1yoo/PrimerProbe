# Introduction
This script is designed to identify the primer set utilized for microbiome data obtained from the Sequence Read Archive (SRA). The initial step in microbiome data analysis involves removing primer sequences. These primers are typically used to amplify desired regions of the 16S rRNA, and various primers may be employed for each experiment. To analyze data from multiple groups or institutions, it's essential to ascertain which primers were used for each dataset. Therefore, we coded a simple
script for this purpose.

# Basic approch
The basic approach involves counting the occurrences of known primer sequences in the top 1000 reads of the FASTQ files and returning the most frequently observed primer set. However, a more advanced approach would entail performing a BLAST search against a 16S sequence database with the top 1000 sequences from the FASTQ files to estimate the target region. This idea is planned for implementation in the future."

# Primer information

1. Primer set targeting the V3-V4 region:
   - Primer IDs: 341F, 805R
   - Forward primer: 341F (5'-CCTACGGGNGGCWGCAG-3')
   - Reverse primer: 805R (5'-GACTACHVGGGTATCTAATCC-3')
   - Typically used in 16S metagenomics primer sets provided by companies like Illumina and Thermo Fisher Scientific.

2. Primer set targeting the V4 region:
   - Primer IDs: 515F, 806R
   - Forward primer: 515F (5'-GTGCCAGCMGCCGCGGTAA-3')
   - Reverse primer: 806R (5'-GGACTACHVGGGTWTCTAAT-3')
   - Included in primer sets provided by companies like Illumina and Thermo Fisher Scientific.

3. Primer set targeting the V1-V3 region:
   - Primer IDs: 27F, 534R
   - Forward primer: 27F (5'-AGAGTTTGATCMTGGCTCAG-3')
   - Reverse primer: 534R (5'-CTGCTGCCTYCCGTA-3')
   - Widely used in various genomic studies, particularly in microbial diversity analysis.

4. Primer set targeting the V4-V5 region:
   - Primer IDs: 515F, 926R
   - Forward primer: 515F (5'-GTGYCAGCMGCCGCGGTAA-3')
   - Reverse primer: 926R (5'-GGACTACNVGGGTHTCTAAT-3')
   - Utilized in various ecological studies and environmental sampling.

5. Primer set targeting the V5-V6 region:
   - Primer IDs: 784F, 1061R
   - Forward primer: 784F (5'-GYACWCACCWGAGCTG-3')
   - Reverse primer: 1061R (5'-CCGTCACCTTGTTACGACTT-3')
   - Particularly used in studying microbial diversity in soil and aquatic ecosystems.

6. Primer set targeting the V6-V8 region:
   - Primer IDs: 968F, 1401R
   - Forward primer: 968F (5'-AACGCGAAGAACCTTAC-3')
   - Reverse primer: 1401R (5'-CGGTGTGTACAAGACCC-3')
   - Frequently employed in microbial ecology studies and community profiling.

This information was gathered through ChatGPT.
