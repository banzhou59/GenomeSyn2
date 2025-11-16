## For detailed instructions on how to identify SNPs and produce the required VCF file.

Accurately deciphering the genetic composition and parental origins underlying progeny genotypes is a fundamental step for understanding the genetic mechanisms of complex traits and improving breeding materials. This is particularly important in multiparental populations (MPPs), such as nested association mapping (NAM) populations, multiparent advanced generation intercross (MAGIC) populations, and natural populations, where it allows the inference of the origins and contribution proportions of genomic segments from different genetic backgrounds.

#### Using NGS(Next-Generation Sequencing) data

To generate the VCF file, first select a reference genome. Then, use BWA to align the resequencing reads of both parents and progeny to the reference genome. After alignment, use BCFtools to perform SNP calling on the mapped data.

```
bwa index Reference.fasta

bwa mem -t 30 Reference.fasta Parent1_R1.fastq.gz Parent1_R2.fastq.gz | samtools sort -o Parent1.sort.bam

samtools index Parent1.sort.bam

bwa mem -t 30 Reference.fasta Parent2_R1.fastq.gz Parent2_R2.fastq.gz | samtools sort -o Parent2.sort.bam

samtools index Parent2.sort.bam

bwa mem -t 30 Reference.fasta Progeny_R1.fastq.gz Progeny_R2.fastq.gz | samtools sort -o Progeny.sort.bam

samtools index Progeny.sort.bam

bcftools mpileup -Ou -f Reference.fasta -a FORMAT/AD,FORMAT/DP --threads 30 Parent1.sort.bam Parent2.sort.bam Progeny.sort.bam | bcftools call -mv --ploidy 1 -Oz --threads 30 -o Output1.raw.vcf.gz

bcftools view -v snps Output1.raw.vcf.gz -Oz -o Output1.raw.snps.vcf.gz

bcftools index Output1.raw.snps.vcf.gz

bcftools filter -s LowQual -e 'QUAL<30 || INFO/DP<5' Output1.raw.snps.vcf.gz -Oz -o Output1.raw.snps.filtered.vcf.gz

bcftools index Output1.raw.snps.filtered.vcf.gz

bcftools view -f PASS Output1.raw.snps.filtered.vcf.gz -Oz -o Output1.snps.pass.vcf.gz

bcftools index Output1.snps.pass.vcf.gz

gunzip Output1.snps.pass.vcf.gz
```

#### Using Genome data

Select the progeny genome as the reference, and use Minimap2 to align both the parental and progeny genomes to the reference. Structural variants are then jointly called across multiple samples using BCFtools, followed by filtering based on depth, quality, and SNP sites. When using the progeny genome as the reference, it is important to align the progeny genome against itself. This self-alignment is necessary for subsequent SNP calling to accurately identify the number of SNPs with full identity.

```
minimap2 -t 30 -ax asm5 Progeny.fasta Parent1.fasta | samtools sort -o Progeny_vs_Parent1.sort.bam

samtools index Progeny_vs_Progeny.sort.bam

minimap2 -t 30 -ax asm5 Progeny.fasta Parent2.fasta | samtools sort -o Progeny_vs_Parent2.sort.bam

samtools index Progeny_vs_Parent1.sort.bam

# Self-alignment of the progeny genome

minimap2 -t 30 -ax asm5 Progeny.fasta Progeny.fasta | samtools sort -o Progeny_vs_Progeny.sort.bam

samtools index Progeny_vs_Parent2.sort.bam

bcftools mpileup -Ou -f Reference.fasta -a FORMAT/AD,FORMAT/DP --threads 30  Progeny_vs_Parent1.sort.bam Progeny_vs_Parent2.sort.bam Progeny_vs_Progeny.sort.bam | bcftools call -mv --ploidy 1 -Oz --threads 30 -o Output2.raw.vcf.gz

bcftools view -v snps Output2.raw.vcf.gz -Oz -o Output2.raw.snps.vcf.gz

bcftools index Output2.raw.snps.vcf.gz

bcftools filter -s LowQual -e 'QUAL<30 || INFO/DP<5' Output2.raw.snps.vcf.gz -Oz -o Output2.raw.snps.filtered.vcf.gz

bcftools index Output2.raw.snps.filtered.vcf.gz

bcftools view -f PASS Output2.raw.snps.filtered.vcf.gz -Oz -o Output2.snps.pass.vcf.gz

bcftools index Output2.snps.pass.vcf.gz

gunzip Output2.snps.pass.vcf.gz
```


## VCF File Format for GenomeSyn2

The VCF file used in GenomeSyn2 should follow the standard VCF format and include the following requirements:
### 1. Chromosome lengths
All chromosomes present in the VCF must have their lengths specified in the header using the ##contig=<ID=chr,length=xxxx> format. This information is mandatory for GenomeSyn2 to correctly calculate SNP density and identity across the genome.
### 2. Sample colors (optional)
Colors for each sample can be specified in the header using the ##color=<Sample=SampleName,color="#RRGGBB"> format. If a color is not provided, GenomeSyn2 will assign a default color automatically.
### 3. SNP genotype data
The VCF header line must include all sample names (varieties) after the FORMAT column. The body contains SNP genotype data in standard VCF format.
```
# less ./vcf/parents.progeny.snps.genotype.Chr01.vcf
##contig=<ID=Chr01,length=45027022>
##contig=<ID=Chr02,length=37301368>
##contig=<ID=Chr03,length=39893253>
##contig=<ID=Chr04,length=37319239>
##contig=<ID=Chr05,length=31307418>
##contig=<ID=Chr06,length=31921180>
##contig=<ID=Chr07,length=30877072>
##contig=<ID=Chr08,length=30492302>
##contig=<ID=Chr09,length=24892599>
##contig=<ID=Chr10,length=25690566>
##contig=<ID=Chr11,length=34100580>
##contig=<ID=Chr12,length=26942889>
##color=<Sample=HHZ,color="#39A5D6">
##color=<Sample=Kasalath,color="#43A98C">
##color=<Sample=N1-10-10A,color="#B8D891">
#CHROM POS	 ID	REF ALT	QUAL    FILTER INFO FORMAT HHZ Kasalath N1-10-10A
Chr01  5741  .  G   T   192.885	PASS   *    *      0/0 1/1      0/0
Chr01  5766  .  T   C   199.58	PASS   *    *      0/0 1/1      0/0
Chr01  55497 .  A   C   111.883 PASS   *    *      0/0 1/1      1/1
...
```


```
# less ./vcf/NAM.SNP_identity.50Kb.bed
#Chr  Start  End    HHZ Kasalath N1-10-10A
Chr01  1      50000  19  0        19
Chr01  50001  100000 40  2        42
Chr01  100001 150000 89  8        97
...

# less ./vcf/NAM.SNP_density.50Kb.bed
#Chr  Start   End    HHZ  Kasalath N1-10-10A
Chr01  1      50000  1    236      1
Chr01  50001  100000 5    425      2
Chr01  100001 150000 13   390      9
...
```




