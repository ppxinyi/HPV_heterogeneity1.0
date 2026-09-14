# HPV Integration Structure Reconstruction Pipeline

This repository contains the code used to reconstruct HPV integration structures from ONT long-read sequencing data.

The pipeline identifies human–viral junctions from aligned BAM files, extracts junction-supporting reads, and reconstructs local HPV integration structures through junction-guided assembly.

---

# Overview

## Input
- ONT long-read BAM file
- Reads aligned against combined Human + HPV reference genome

Example:

```bash
sample.sorted.bam
```


# Requirements

## Python
- Python >= 3.9

## Install environment

```bash
conda create -n hpv_assembly python=3.10 -y
conda activate hpv_assembly
conda install -c bioconda pysam -y
conda install pandas numpy networkx matplotlib -y
conda activate hpv_assembly
```
---

# Input BAM File

The BAM file should:

- Be sorted and indexed
- Contain supplementary alignments (`SA` tags)
- Be aligned against a combined Human + HPV reference genome

Example reference:

```bash
hg38_HPV16.fa
```

---

# Run example

Before running the pipeline, users need to modify the input and output paths in the Python script.

Example:

```bash
./assembly.sh \
    -b test.bam \
    -s sampleID \
    -c all(9,10,11 any chr want to run) \
    -o output_folder_name
```

---

## Important Note

This pipeline is still under active development, and additional features, automation steps, and parameter optimization will be released in future updates.

At the current stage, some reconstructed HPV integration structures still involve manual inspection and selection by comparing reconstructed tree structures back to BAM file alignments and read-level visualizations to ensure structural accuracy.

The current version is primarily intended for method development, benchmarking, and visualization of complex HPV integration events.
