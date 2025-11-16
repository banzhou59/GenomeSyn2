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




