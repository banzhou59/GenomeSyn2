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

<div align="center"><img src="images/Figure1.png" alt="GenomeSyn2 Figure1" width="800"/></div>

## Contents
- [Introduction](#Introduction)
- [Installation](#Installation)
- [Usage](#Usage)
  - [Test](#Test)
  - [Inputs](#Inputs) 
  - [Outputs](#Outputs)
  - [Proof curation](#Proof-annotation)
- [All available options](#All-available-options)
- [Update history](#Update-history)
- [Citation](#Citation)

## Introduction

The tool integrates multi-genome synteny relationships with diverse statistical and annotation features, enabling systematic exploration of genome rearrangement patterns, tracing of genomic region origins, and quantitative comparison of parental contributions, thus providing an efficient and intuitive visualization platform for comparative genomics and population genetics research.


## Installation

GenomeSyn2 provides three options for installing the required dependencies:

1. Install GenomeSyn2 directly from Bioconda:

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

3. Install dependencies manually

You may manually install all required software listed in the documentation. This option is suitable for users who prefer full control over their environment.

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

**[MUMmer](https://anaconda.org/bioconda/mummer)** — for whole-genome alignment and synteny detection

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



### Test Installation

If installed successfully, you can check the version with:

```
GenomeSyn2 --help

# Download the test data of GenomeSyn2

wget https://cbi.gxu.edu.cn/zwzhou/GenomeSyn/GenomeSyn2_example_data.zip
```

# Quick start
1) The following commands demonstrate how to run GenomeSyn2 using different genome alignment software.

```
   GenomeSyn2 --align mummer --genome ./genome_path/ --outdir ./mummer/ --thread 30 > GS2.mummer.log

   GenomeSyn2 --align minimap2 --genome ./genome_path/ --outdir ./minimap2/ --thread 30 > GS2.minimap2.log
```

2) The following commands demonstrate how to run GenomeSyn2 for protein alignment using different alignment tools.

```
   GenomeSyn2 --align blastp --genome ./genome_path/ --gff ./gene_data/ --outdir ./blastp/ --thread 30 > GS2.blastp.log

   GenomeSyn2 --align mmseqs --genome ./genome_path/ --gff ./gene_data/ --outdir ./mmseqs/ --thread 30 > GS2.mmseqs.log

   GenomeSyn2 --align diamond --genome ./genome_path/ --gff ./gene_data/ --outdir ./diamond/ --thread 30 > GS2.diamond.log
```

3) Calculate SNP density and SNP identity from a VCF file to visualize multi-parental origin contributions:

```
   GenomeSyn2 --type identity --vcf ./parents.progeny.snps.genotype.Chr01.vcf --bin 50000 > GS2.vcf.log
```

4) Based on SNP density and SNP identity statistics, plot the multi-parental origins contribution:

```
   GenomeSyn2 --type identity --identity ./SNP_identity.50Kb.bed --density ./SNP_density.50Kb.bed > GS2.vcf.log
```


You can run GenomeSyn2 with a custom configuration file:
```
GenomeSyn2 --conf total.conf
```

# Configuration File Structure

Please refer to [Configuration_File.README.md](Configuration_File.README.md) for details on the configuration file format.


















