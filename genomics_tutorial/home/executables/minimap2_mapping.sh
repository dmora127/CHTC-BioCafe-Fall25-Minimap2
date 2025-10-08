#!/bin/bash
# Use minimap2 to map the basecalled reads to the reference genome
./minimap2 -ax map-ont "$1" "$2" > "mapped_${2}_reads_to_genome.sam"

# Use samtools to sort our mapped reads BAM, required for downstream analysis
samtools sort "mapped_${2}_reads_to_genome.sam" -o "mapped_${2}_reads_to_genome_sam_sorted.bam"