# variants-pipelines

To create a variant calling pipeline that takes in fastq and outputs classified variant based in ACMG framework.

Here is the general approach:

1. Quality control: Perform quality control of the raw sequencing data to ensure that it is of high quality and meets the desired standards. This can be done using tools such as FastQC or MultiQC.

2. Read alignment: Align the sequencing reads to the reference genome using an alignment tool such as BWA or Bowtie2.

3. Variant calling: Call variants from the aligned reads using a variant caller such as GATK or FreeBayes. This step will produce a VCF file containing all the variants.

4. Variant annotation: Annotate the variants with additional information such as their impact on gene function, population frequency, and clinical significance using a variant annotation tool such as ANNOVAR or VEP.

5. Variant classification: Classify the variants according to the ACMG guidelines. This involves assessing the pathogenicity of the variants based on criteria such as their impact on protein function, their frequency in the general population, and their association with known disease-causing mutations. This step can be performed manually or using software such as InterVar or VarSome.

And the pipeline for the above approach using GATK could be:
```
# Index the reference genome
bwa index reference_genome.fasta

# Align reads to the reference genome
bwa mem -t 4 reference_genome.fasta sample_1.fastq sample_2.fastq | samtools sort -o aligned_reads.bam

# Mark duplicates and index the BAM file
gatk MarkDuplicates -I aligned_reads.bam -O marked_duplicates.bam -M marked_dup_metrics.txt
samtools index marked_duplicates.bam

# Call variants using GATK HaplotypeCaller
gatk HaplotypeCaller -R reference_genome.fasta -I marked_duplicates.bam -O raw_variants.vcf

# Annotate variants using ANNOVAR
annovar/table_annovar.pl raw_variants.vcf annovar/humandb/ -buildver hg38 -out annotated_variants -remove -protocol refGene,dbnsfp35c -operation g,f -nastring .

# Classify variants using InterVar
python intervar/InterVar.py -i annotated_variants.hg38_multianno.txt -o classified_variants.txt -b hg38 --refversion hg38 --buildver hg38 --input_type avinput
```

And the Non-GATK approach could be:
```
# Index the reference genome
bwa index reference.fasta

# Align the reads to the reference genome
bwa mem -t 8 reference.fasta sample_R1.fastq.gz sample_R2.fastq.gz > sample.sam

# Convert the SAM file to BAM and sort it
samtools view -bS sample.sam | samtools sort -@ 8 -o sample.sorted.bam

# Index the sorted BAM file
samtools index sample.sorted.bam

# Call variants using Freebayes
freebayes -f reference.fasta sample.sorted.bam > sample.vcf

# Annotate the variants with ANNOVAR or other annotation tools based on ACMG guidelines.
```