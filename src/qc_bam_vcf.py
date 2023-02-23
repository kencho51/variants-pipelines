import pysam
import vcf

# Define paths to input data
ref_fasta = "/inputs/reference/reference.fasta"
input_bam = "/path/to/input.bam"
input_vcf = "/path/to/input.vcf"

# Load reference genome
ref = pysam.FastaFile(ref_fasta)

# Load input BAM file and evaluate alignment quality
bam = pysam.AlignmentFile(input_bam, "rb")
for read in bam:
    # Check mapping quality
    if read.mapping_quality < 30:
        print(f"Low mapping quality read found: {read.query_name}")
    # Check alignment position
    if read.reference_start < 0 or read.reference_end > len(ref):
        print(f"Out-of-bounds alignment found: {read.query_name}")
    # Check read length
    if len(read.query_sequence) < 50:
        print(f"Short read found: {read.query_name}")

# Load input VCF file and evaluate variant quality
vcf_reader = vcf.Reader(filename=input_vcf)
for record in vcf_reader:
    # Check genotype quality
    if record.QUAL < 50:
        print(f"Low-quality genotype found: {record.CHROM}:{record.POS}")
    # Check variant type
    if record.var_type == "indel":
        print(f"Indel found: {record.CHROM}:{record.POS}")

# Clean up resources
bam.close()