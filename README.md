## Overview
Oxford Nanopore's (ONT) PCR-free Rapid Barcoding Kits (RBK) offer significant advantages, such as speed and cost efficiency, as well as read lengths of several KB's. However, despite the recent improvements of the technology, single base errors can still be a concern, particularly in the context of bacterial methylation sites. PCR-based library preparation protocols overcome this limitation, however long-range PCR comes with its own set of challenges such as limited read lengths. 
Here we present an experimental workflow combining PCR-free and PCR-based ONT reads, developed with the aim of matching accuracy of Illumina - ONT hybrid approaches. We successfully tested the workflow on a limited dataset consisting of 
20 <em>Corynebacterium diphtheriae</em> and 20 Vancomycin-resistant <em>Enterococcus faecium</em> (VRE) isolated collected during outbreaks.

**Note:** NanoporeHybridKP_v0.1.py is a prototype script, specifically designed to correct single base substitutions. The thresholds used for defining ambiguous assembly positions have not been optimized, and the workflow has only been tested with the two species mentioned above. Please use with caution. 
## Dependencies
- minimap2 (version 2.24)
- pysamstats (version 1.1.2)

## Input Data
- PCR-free Oxford Nanopore reads (tested with RBK)
- PCR-based Oxford Nanopore reads (RPB)
- PCR-free Oxford Nanopore assemblies (e.g. generated with Flye)

## Workflow

### Step 1: Generate pileup files for RBK reads
```bash
# Map PCR-free reads to the PCR-free assembly
conda activate minimap2
minimap2 -ax map-ont --secondary=no rbk_consensus.fasta rbk_reads.fastq -t 8 > rbk_rbk.sam
conda deactivate

conda activate samtools
samtools view -bS rbk_rbk.sam | samtools sort - > rbk_rbk.bam
samtools index rbk_rbk.bam && rm rbk_rbk.sam
conda deactivate

# Generate pysamstats pileup
conda activate pysamstats
pysamstats --type variation_strand --truncate --pad --fasta rbk_consensus.fasta rbk_rbk.bam 1> rbk_rbk.variation.txt
conda deactivate
```

### Step 2: Generate pileup files for RPB reads
```bash
# Map PCR-based reads to the PCR-free assembly
conda activate minimap2
minimap2 -ax map-ont --secondary=no rbk_consensus.fasta rpb_reads.fastq -t 8 > rbk_rpb.sam
conda deactivate

conda activate samtools
samtools view -bS rbk_rpb.sam | samtools sort - > rbk_rpb.bam
samtools index rbk_rpb.bam && rm rbk_rpb.sam
conda deactivate

# Generate pysamstats pileup
conda activate pysamstats
pysamstats --type variation_strand --truncate --pad --fasta rbk_consensus.fasta rbk_rpb.sam 1> rbk_rpb.variation.txt
conda deactivate
```

### Step 3: Run assembly correction
```bash
.venv/bin/activate
python nanoporehybridkp_v0.1.py -r rbk_rbk.variation.txt -p rbk_rpb.variation.txt -o out/corrected.fasta
```