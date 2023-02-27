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

And the example outputs could be something like:
```
Variant Classification Report for Patient X

Variant ID: 123456
Gene Name: BRCA1
Variant Type: Missense
Allele Frequency: 0.001
Functional Impact Score: 0.98
In-silico Prediction: Damaging (SIFT), Likely Pathogenic (PolyPhen)

ACMG Classification:
- PM2: This variant is absent in population databases with a frequency of less than 0.01.
- PP3: Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact, etc.).
- PP4: Patient's phenotype or family history is highly specific for a disease with a single genetic etiology, and this variant has been identified in multiple patients with the same phenotype.

Overall Classification: Likely Pathogenic (class 4)

Comments:
This missense variant in BRCA1 has a functional impact score of 0.98 and is predicted to be damaging by SIFT and likely pathogenic by PolyPhen. It is not present in population databases with a frequency greater than 0.01, and multiple lines of computational evidence support its deleterious effect on the gene. The patient's family history is highly specific for breast cancer, and this variant has been identified in multiple patients with the same phenotype. Therefore, this variant is classified as Likely Pathogenic (class 4) according to ACMG guidelines.
```

Or in this format:
```
Variant Interpretation and Reporting
We identified a total of X variants in the patient's exome sequencing data. After applying our filtering criteria based on frequency, pathogenicity, and ACMG guidelines, we have narrowed down to Y variants of potential clinical significance.

Based on our analysis, we have identified the following variants that are likely responsible for the patient's symptoms:

Gene: ABCA4
Variant: c.5461-10T>C
Zygosity: Homozygous
ACMG classification: Likely Pathogenic
Interpretation: This variant has been previously reported in association with Stargardt disease, which is consistent with the patient's symptoms.

Gene: BRCA2
Variant: c.68-7T>A
Zygosity: Heterozygous
ACMG classification: Likely Benign
Interpretation: This variant is predicted to be benign based on its frequency in the general population and the fact that it is not located in a known functional domain of the protein.

We recommend that the patient undergo further confirmatory testing, such as Sanger sequencing or targeted gene panel testing, to confirm these findings. Additionally, genetic counseling should be offered to the patient and their family members to discuss the implications of these results and any potential treatment options or management strategies.
```

### How to execute the scripts
```
# Install the db
$ pyensembl install --release 97 --species homo_sapiens
# Build the image
$ cd variants-pipelines
$ docker-compose run --rm console bash
root@744a84ca81be:/# python3 --version
Python 3.6.9
root@744a84ca81be:/# pytest --version
pytest 7.0.1
root@744a84ca81be:/# ls app
Dockerfile  LICENSE  README.md  docker-compose.yml  inputs  outputs  requirements  src  tests
root@744a84ca81be:/# pytest app/src/test_interpret_variants.py 

```


### Reference data source
1. [mayo-test.vcf](https://bioinformaticstools.mayo.edu/research/vcf-miner-sample-vcfs/)