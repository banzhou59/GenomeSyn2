# [GenomeSyn-II: A Comparative Genomics Framework Integrating Synteny Visualization](https://cbi.gxu.edu.cn/zwzhou/GenomeSyn/GenomeSyn2-v1.0.0.tar.gz)

[![GitHub release](https://img.shields.io/github/v/release/banzhou59/GenomeSyn2?label=release)](https://github.com/banzhou59/GenomeSyn2/releases)
[![License](https://img.shields.io/badge/license-GPL-blue.svg)](LICENSE)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/genomesyn2/README.html)
[![Bioconda Downloads](https://anaconda.org/bioconda/genomesyn2/badges/downloads.svg)](https://anaconda.org/bioconda/genomesyn2)

The tool integrates multi-genome synteny relationships with diverse statistical and annotation features, enabling systematic exploration of genome rearrangement patterns, tracing of genomic region origins, and quantitative comparison of parental contributions, thus providing an efficient and intuitive visualization platform for comparative genomics and population genetics research.

|    Author   |  E-mail                                                                |
|-------------|------------------------------------------------------------------------|
| [Zu-Wen Zhou](https://github.com/banzhou59)                 | `784012725@qq.com`     |
| [Ling-Ling Chen](https://lst.gxu.edu.cn/info/1077/2334.htm) | `jmsong@swu.edu.cn`    |
| [Jia-Ming Song](http://agronomy.swu.edu.cn/info/1081/5164.htm) | `llchen@gxu.edu.cn` |

## Dependencies

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


After installing the dependencies with Conda and activating the Conda environment, install GenomeSyn2 using **one** of the following options:

```
conda env create -f environment.yml

conda activate GenomeSyn2
```

## Installation

### Installation options

1. Use git to clone the project repository:
GenomeSyn2 does not require compilation and can be run directly.
To enable running GenomeSyn2.pl from any directory, it is recommended to add it to the environment variable $PATH.

```
git clone https://github.com/banzhou59/GenomeSyn2.git

cd ./GenomeSyn2/bin/

chmod +x *.pl *.sh

# To make this change permanent, add the line to your ~/.bashrc file and then run source ~/.bashrc to apply it immediately.

export PATH=/your_path/GenomeSyn2/bin:$PATH
```

2. Download the GenomeSyn2.tar.gz

```
wget https://cbi.gxu.edu.cn/zwzhou/GenomeSyn/GenomeSyn2-1.0.0.tar.gz

tar -zxvf GenomeSyn2.tar.gz

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

# Configuration file format

```
less total.conf

[genome_info]
# gonomes_filetype = (fasta/bed)
gonomes_filetype = bed
gonomes_list = chr_length.info.tsv

[synteny_info]
# line_type = (curve/line)
line_type = curve
synteny_list = synteny.info.tsv

[save_info]
# figure_type = (svg/pdf)
figure_type = pdf
savefig1 = GenomeSyn2.figure1.pdf
savefig2 = GenomeSyn2.figure2.pdf

[centromere_info]
centromere_list=centromere.info.tsv

[telomere_info]
telomere_list=telomere.info.tsv
telomere_color=#441680
opacity=100%

[show_region]
# region = (genome_Name:ChrID:start-end)
region = MH63:Chr10:24,850,000-24,885,000
gene_list = gene.info.tsv

[anno_info:no]
anno_number=[1,2,3,4,5,6,7]
anno_name=[PAV,SNP,TE,GC Content,Gypsy,Copia,Gene density]
anno_color=['#5FB6DE','#0000FF','#3774B9','#000000','#00FF00','#F5F57A','#368F5C']
anno_type=[rectangle,barplot,barplot,lineplot,lineplot,lineplot,heatmap]
anno_position=[top,top,bottom,top,bottom,bottom,middle]
anno_height=[5,5,5,5,5,5,5]
min_max_value=[normal,auto,normal,0.4:0.5,normal,normal,normal]
anno_window=[none,none,100000,none,100000,100000,100000]
opacity=[50%,100%,100%,100%,100%,100%,100%]
file_type=[bed,bed,gff3,bed,gff3,gff3,gff3]
filter_type=[none,none,none,none,none,none,gene]
anno_list=[PAV.info.tsv,SNP.info.tsv,TE.info.tsv,GC.info.tsv,Gypsy.info.tsv,Copia.info.tsv,gene.info.tsv]
```














