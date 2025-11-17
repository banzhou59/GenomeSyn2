# [GenomeSyn-II: A Comparative Genomics Framework Integrating Synteny Visualization](https://cbi.gxu.edu.cn/zwzhou/GenomeSyn/GenomeSyn2-v1.0.0.tar.gz)

[![GitHub release](https://img.shields.io/github/v/release/banzhou59/GenomeSyn2?display_name=tag&sort=semver)](https://github.com/banzhou59/GenomeSyn2/releases)
[![License](https://img.shields.io/badge/license-GPL-blue.svg)](LICENSE)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/genomesyn2/README.html)
[![Bioconda Downloads](https://anaconda.org/bioconda/genomesyn2/badges/downloads.svg)](https://anaconda.org/bioconda/genomesyn2)

<table style="width:100%; max-width:100%; border-collapse: collapse;">
    <tr>
      <th>Author</th>
      <th>E-mail</th>
    </tr>
    <tr>
      <td><a href="https://github.com/banzhou59">Zu-Wen Zhou</a></td>
      <td><code>784012725@qq.com</code></td>
    </tr>
    <tr>
      <td><a href="https://lst.gxu.edu.cn/info/1077/2334.htm">Ling-Ling Chen</a></td>
      <td><code>jmsong@swu.edu.cn</code></td>
    </tr>
    <tr>
      <td><a href="http://agronomy.swu.edu.cn/info/1081/5164.htm">Jia-Ming Song</a></td>
      <td><code>llchen@gxu.edu.cn</code></td>
    </tr>
</table>

## Contents
- [Introduction](#Introduction)
- [Installation](#Installation)
- [Usage](#Usage)
  - [Test](#Test-installation)
  - [Inputs](#Inputs)
    - [1. Quick start](#1-Quick-Start-Guide-for-Running-GenomeSyn2)
    - [2. Drawing genome synteny diagrams and annotation information](#2-drawing-genome-synteny-diagrams-and-annotation-information)
    - [3. Local gene structure view](#3-local-gene-structure-view)
    - [4. Multi-parental origins contribution view](#5-multi-parental-origins-contribution-view)
      - [a. Calculate SNP density and SNP identity from a VCF file](#a-calculate-snp-density-and-snp-identity-from-a-vcf-file-to-visualize-multi-parental-origin-contributions)
      - [b. Plot multi-parental origins contribution](#b-based-on-snp-density-and-snp-identity-statistics-plot-the-multi-parental-origins-contribution)
  - [Outputs](#Outputs)
  - [Configuration File Structure](#Configuration-File-Structure)
- [Citation](#citation)

## Introduction

The tool integrates multi-genome synteny relationships with diverse statistical and annotation features, enabling systematic exploration of genome rearrangement patterns, tracing of genomic region origins, and quantitative comparison of parental contributions, thus providing an efficient and intuitive visualization platform for comparative genomics and population genetics research.

<div align="center"><img src="images/Figure1.png" alt="GenomeSyn2 Figure1" width="800"/></div>

## Installation

GenomeSyn2 provides three options for installing the required dependencies:

1.Install GenomeSyn2 directly from Bioconda:

```
conda create -n GenomeSyn2 -c bioconda genomesyn2

conda install -c bioconda genomesyn2
```

2.conda installation environment.yml

```
wget https://github.com/banzhou59/GenomeSyn2/releases/download/v1.0.0/GenomeSyn2-1.0.0.tar.gz

tar -zxvf GenomeSyn2-1.0.0.tar.gz

cd ./GenomeSyn2/

conda env create -f environment.yml

conda activate GenomeSyn2

# To make this change permanent, add the line to your ~/.bashrc file and then run source ~/.bashrc to apply it immediately.

export PATH=/your_path/GenomeSyn2/bin:$PATH
```

3.Install dependencies manually

You may manually install all required software listed in the documentation.

***GenomeSyn2 Environment Requirements***

GenomeSyn2 requires a Linux environment with both Perl and Python installed, along with several bioinformatics tools and libraries.
All dependencies can be easily installed using a Conda environment defined in the provided environment.yml file.

***1. Perl and perl modules***

**[Perl](https://www.perl.org)** ≥ 5.32 — the main language used for running GenomeSyn2 scripts

**[perl-bioperl-core](https://anaconda.org/bioconda/perl-bioperl-core)** — provides BioPerl functionalities for sequence and annotation processing

**[perl-svg](https://anaconda.org/bioconda/perl-svg)** — supports generation of scalable vector graphics (SVG) output

***2. Python and python modules***

**[Python](https://www.python.org)** ≥ 3.8 — required for auxiliary data visualization and conversion tools

**[cairosvg](https://cairosvg.org/)** — used for converting SVG files to PNG or PDF formats

***3. Other bioinformatics tools***

**[MUMmer4](https://anaconda.org/bioconda/mummer4)** — for whole-genome alignment and synteny detection

[minimap2](https://anaconda.org/bioconda/minimap2) — for fast and accurate sequence alignment

**[gffread](https://anaconda.org/bioconda/gffread)** — for extracting transcript and protein sequences from GFF/GTF files

**[Seqkit](https://anaconda.org/bioconda/seqkit)** — for efficient FASTA/FASTQ file manipulation

**[blast](https://anaconda.org/bioconda/blast)** — for sequence similarity searches

**[DIAMOND](https://anaconda.org/bioconda/diamond)** — for fast protein alignment

**[mmseqs2](https://anaconda.org/bioconda/mmseqs2)** — for large-scale sequence clustering and homology search

```
# Clone the github repository for GenomeSyn2.

git clone https://github.com/banzhou59/GenomeSyn2.git

cd ./GenomeSyn2/bin/

chmod +x *.pl *.sh

# or

wget https://github.com/banzhou59/GenomeSyn2/releases/download/v1.0.0/GenomeSyn2-1.0.0.tar.gz

tar -zxvf GenomeSyn2-1.0.0.tar.gz

# To make this change permanent, add the line to your ~/.bashrc file and then run source ~/.bashrc to apply it immediately.

export PATH=/your_path/GenomeSyn2/bin:$PATH
```

## Usage

## Test Installation

You can use the following commands to check whether GenomeSyn2 was installed successfully:

```
GenomeSyn2 --help

GenomeSyn2 --version
```
You may then download the example dataset to perform a test run:
```
# Download the test data of GenomeSyn2

wget https://cbi.gxu.edu.cn/zwzhou/GenomeSyn/GenomeSyn2_example_data.zip
```

## 1. Quick Start Guide for Running GenomeSyn2
a) Genome alignment
The following commands demonstrate how to run GenomeSyn2 using different genome alignment tools:

```
   GenomeSyn2 --align mummer --genome ./genome_path/ --outdir ./mummer/ --thread 30 > GS2.mummer.log

   GenomeSyn2 --align minimap2 --genome ./genome_path/ --outdir ./minimap2/ --thread 30 > GS2.minimap2.log
```

b) Protein alignment
The following commands demonstrate how to run GenomeSyn2 for protein alignment using different tools:

```
   GenomeSyn2 --align blastp --genome ./genome_path/ --gene ./gene_data/ --outdir ./blastp/ --thread 30 > GS2.blastp.log

   GenomeSyn2 --align mmseqs --genome ./genome_path/ --gene ./gene_data/ --outdir ./mmseqs/ --thread 30 > GS2.mmseqs.log

   GenomeSyn2 --align diamond --genome ./genome_path/ --gene ./gene_data/ --outdir ./diamond/ --thread 30 > GS2.diamond.log
```

Important note on input files:
Before running these commands, make sure that the genome FASTA files and the corresponding gene annotation files are renamed and sorted in numerical order, starting with a number followed by a dot (e.g., 1., 2., 3., etc.). The numbering of the genome files must exactly match the numbering of their corresponding annotation files. GenomeSyn2 will perform pairwise alignments based on this numerical order, ensuring that each genome is correctly matched with its annotation file.
Example directory structure:
```
genome_path/
├── 1.MH63RS3.fasta
├── 2.T.mark.fasta
└── 3.Y.mark.fasta

gene_data/
├── 1.MH63.gene.gff3
├── 2.T.gene.gff3
└── 3.Y.gene.gff3
```

## 2. Drawing genome synteny diagrams and annotation information
```
GenomeSyn2 --conf ? > anno.conf
GenomeSyn2 --anno ? >> anno.conf
GenomeSyn2 --conf anno.conf

# less anno.conf
------------------------------------------------------------------------------------------------------------
[genome_info]
# gonomes_filetype = (fasta/bed)
# Type of genome description (fasta or bed)
gonomes_filetype = bed
# List of genome chromosome sizes or fasta files
gonomes_list = chr_length.info.tsv

[synteny_info]
# line_type = (curve/line)
# Style for connecting syntenic blocks. (curve or line)
line_type = curve
# File containing synteny information between genomes.
synteny_list = synteny.info.tsv

[save_info]
# figure_type = (svg/pdf/png)
# File format for saving figures. (svg, pdf, png)
figure_type = pdf
# savefig1 / savefig2: Output filenames for figure 1 and figure 2.
savefig1 = GenomeSyn2.figure1.pdf
savefig2 = GenomeSyn2.figure2.pdf

[centromere_info]
centromere_list = centromere.info.tsv

[telomere_info]
telomere_list = telomere.info.tsv
telomere_color = #441680
opacity = 100%

[anno_info]
anno_number = [1,2,3,4,5,6,7]
anno_name = [PAV,SNP,TE,GC Content,Gypsy,Copia,Gene density]
anno_color = ['#5FB6DE','#0000FF','#3774B9','#000000','#00FF00','#F5F57A','#368F5C']
anno_type = [rectangle,barplot,barplot,lineplot,lineplot,lineplot,heatmap]
anno_position = [top,top,bottom,top,bottom,bottom,middle]
anno_height = [5,5,5,5,5,5,5]
min_max_value = [normal,auto,normal,0.4:0.5,normal,normal,normal]
anno_window = [none,none,100000,none,100000,100000,100000]
opacity = [50%,100%,100%,100%,100%,100%,100%]
file_type = [bed,bed,gff3,bed,gff3,gff3,gff3]
filter_type = [none,none,none,none,none,none,gene]
anno_list = [PAV.info.tsv,SNP.info.tsv,TE.info.tsv,GC.info.tsv,Gypsy.info.tsv,Copia.info.tsv,gene.info.tsv]
------------------------------------------------------------------------------------------------------------
```


## 3. Local gene structure view:
```
GenomeSyn2 --conf local.conf

# less local.conf
------------------------------------------
[genome_info]
gonomes_filetype = bed
gonomes_list = chr_length.info.tsv

[synteny_info]
line_type = curve
synteny_list = synteny.info.tsv

[show_region]
# region = (genome_Name:ChrID:start-end)
region = MH63:Chr10:24,850,000-24,885,000
gene_list = gene.info.tsv

[save_info]
figure_type = pdf
savefig1 = GenomeSyn2.figure1.pdf
------------------------------------------
```


## 4. Multi-parental origins contribution view:

## a) Calculate SNP density and SNP identity from a VCF file to visualize multi-parental origin contributions:

```
GenomeSyn2 --type identity --vcf ./parents.progeny.snps.genotype.Chr01.vcf --bin 50000 > GS2.vcf.log
```

The VCF file can be generated from either resequencing data or whole-genome assemblies. For detailed instructions on how to identify SNPs and produce the required VCF file, please refer to the [VCF_Preparation_Guide.md](docs/VCF_Preparation_Guide.md) document.

## b) Based on SNP density and SNP identity statistics, plot the multi-parental origins contribution:

```
GenomeSyn2 --type identity --identity ./SNP_identity.50Kb.bed --density ./SNP_density.50Kb.bed > GS2.vcf.log
```

## Outputs
🔹1. Outputs of Genome/Protein Alignment Mode (--align <mummer|minimap2|blastp|mmseqs|diamond>):
```
GenomeSyn2 --align <mummer|minimap2> --genome ./genome_path/ --outdir ./mummer/ --thread 30 > GS2.align.log
GenomeSyn2 --align <blastp|mmseqs|diamond> --genome ./genome_path/ --gene ./gene_data/ --outdir ./mummer/ --thread 30 > GS2.align.log
```
- 📁 **<outdir_name>** - Specifies the output directory or folder name.
    - 📁 **fa_bed** - Directory containing processed genome files that record chromosome length information.
    - 📁 **<align_name>** - Directory containing the results of alignments (genome or protein) generated using the --align option.
- 📄 **chr_length.info.tsv** - File recording the paths of chromosome length files for each genome, along with their corresponding variety names and assigned plotting colors. Each file contains the chromosome lengths of the respective genome.
- 📄 **genomes.info.tsv** - File recording the paths to genome files, variety names, and plotting colors.
- 📄 **synteny.info.tsv** - File recording the paths to the alignment results between the corresponding genomes, generated by tools such as MUMmer, Minimap2, BLASTp, MMseqs2, or Diamond.
- 📄 **total.conf** - The configuration file used for the current GenomeSyn2 run, containing plotting parameters, annotation settings, and selected regions for display.
- 📕 **GenomeSyn2.figure1.pdf** - Single-chromosome synteny block view.
- 📕 **GenomeSyn2.figure2.pdf** - Multi-chromosome synteny block view.


🔹2. Outputs of SNP Identity and Density Mode (--type <identity|density|unite>):
```
GenomeSyn2 --type identity --vcf ./parents.progeny.snps.genotype.Chr01.vcf --bin <bin_size> > GS2.vcf.log
```
- 📄 **SNP_identity.<bin_size>.bed** - This BED file reports the number of SNPs with identical genotypes shared among all samples within each genomic bin. The genomic bin size is defined by the --bin parameter.
- 📄 **SNP_density.<bin_size>.bed** - This BED file summarizes the SNP counts per sample within each genomic bin. The bin size is specified by the --bin parameter.
- 📕 **GenomeSyn.<bin_size>.pdf** - Multi-parental origins contribution view.

## Configuration File Structure

Please refer to [Configuration_File.README.md](docs/Configuration_File.README.md) for details on the configuration file format.

## Citation

Zhou, Z., Yu, Z., Huang, X., Liu, J., Guo, Y., Chen, L., Song, J., **2022**. GenomeSyn: a bioinformatics tool for visualizing genome synteny and structural variations. ***J. Genet. Genomics*** 49, 1174-1176. [https://doi.org/10.1016/j.jgg.2022.03.013](https://doi.org/10.1016/j.jgg.2022.03.013)




















