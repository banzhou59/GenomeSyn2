# [GenomeSyn-II: An Extensible Synteny Visualization Tool Supporting Multi-Genome Comparison and Flexible Annotation Display](http://cbi.gxu.edu.cn/GenomeSyn2/)

[![PyPI version](https://img.shields.io/pypi/v/GenomeSyn2.svg)](https://pypi.org/project/GenomeSyn2/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/genomesyn2/README.html)
[![Build Status](https://github.com/yourusername/GenomeSyn2/workflows/build/badge.svg)](https://github.com/yourusername/GenomeSyn2/actions)
[![Downloads](https://pepy.tech/badge/GenomeSyn2)](https://pepy.tech/project/GenomeSyn2)

The tool integrates multi-genome synteny relationships with diverse statistical and annotation features, enabling systematic exploration of genome rearrangement patterns, tracing of genomic region origins, and quantitative comparison of parental contributions, thus providing an efficient and intuitive visualization platform for comparative genomics and population genetics research.

|Author Info|                                                                  |
|-----------|------------------------------------------------------------------|
|  Authors  | [Zu-Wen Zhou](https://github.com/banzhou59)                      |
|           | [Ling-Ling Chen](https://lst.gxu.edu.cn/info/1077/2334.htm)      |
|           | [Jia-Ming Song](http://agronomy.swu.edu.cn/info/1081/5164.htm)   |
|   Email   | jmsong@swu.edu.cn, llchen@gxu.edu.cn, 784012725@qq.com           |

## Dependencies

GenomeSyn2 Environment Requirements

GenomeSyn2 requires a Linux environment with both Perl and Python installed, along with several bioinformatics tools and libraries.
All dependencies can be easily installed using a Conda environment defined in the provided environment.yml file.

```
conda env create -f environment.yml

conda activate GenomeSyn2
```

Perl and perl modules

Perl(https://www.perl.org) ≥ 5.32 — the main language used for running GenomeSyn2 scripts

[perl-bioperl-core](https://anaconda.org/bioconda/perl-bioperl-core) — provides BioPerl functionalities for sequence and annotation processing

[perl-svg](https://anaconda.org/bioconda/perl-svg) — supports generation of scalable vector graphics (SVG) output

Python(https://www.python.org) and python modules

Python ≥ 3.8 — required for auxiliary data visualization and conversion tools

[cairosvg](https://cairosvg.org/) — used for converting SVG files to PNG or PDF formats

External bioinformatics tools

[mummer](https://anaconda.org/bioconda/mummer) — for whole-genome alignment and synteny detection

[minimap2](https://anaconda.org/bioconda/minimap2) — for fast and accurate sequence alignment

[gffread](https://anaconda.org/bioconda/gffread) — for extracting transcript and protein sequences from GFF/GTF files

[seqkit](https://anaconda.org/bioconda/seqkit) — for efficient FASTA/FASTQ file manipulation

[blast](https://anaconda.org/bioconda/blast) — for sequence similarity searches

[diamond](https://anaconda.org/bioconda/diamond) — for fast protein alignment

[mmseqs2](https://anaconda.org/bioconda/mmseqs2) — for large-scale sequence clustering and homology search

After activating the Conda environment install GenomeSyn2 using one of the following options.

## Installation

### Installation options

1. Use pip to install the latest development version directly from this repo.

```
pip install git+https://github.com/banzhou59/GenomeSyn2.git
```

2. Install latest release from PyPi.

```
pip install GenomeSyn2
```

### Test Installation

If installed successfully, you can check the version with:

```
GenomeSyn2 --help
```

# Quick start

GenomeSyn2.pl --align mummer --genome ./genome_path/ --outdir ./mummer/ --thread 30 > GenomeSyn2.mummer.log

GenomeSyn2.pl --align minimap2 --genome ./genome_path/ --outdir ./minimap2/ --thread 30 > GenomeSyn2.minimap2.log

GenomeSyn2.pl --align blastp --genome ./genome_path/ --gff ./gene_data/ --outdir ./blastp/ --thread 30 > GenomeSyn2.blastp.log

GenomeSyn2.pl --align mmseqs --genome ./genome_path/ --gff ./gene_data/ --outdir ./mmseqs/ --thread 30 > GenomeSyn2.mmseqs.log

GenomeSyn2.pl --align diamond --genome ./genome_path/ --gff ./gene_data/ --outdir ./diamond/ --thread 30 > GenomeSyn2.diamond.log
























