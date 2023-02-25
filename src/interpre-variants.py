import pysam
import pyensembl

# Define paths to input data
ref_fasta = "/path/to/reference.fasta"
input_vcf = "/path/to/input.vcf"

# Load reference genome
ref = pysam.FastaFile(ref_fasta)

# Load Ensembl genome annotation
ensembl = pyensembl.EnsemblRelease(release=98)

# Load input VCF file and interpret variants
with open(input_vcf, "r") as f:
    for line in f:
        if line.startswith("#"):
            # Skip VCF header
            continue
        fields = line.strip().split("\t")
        chrom = fields[0]
        pos = int(fields[1])
        ref_allele = fields[3]
        alt_allele = fields[4]
        # Look up gene name and impact of variant
        overlapping_genes = ensembl.genes_at_locus(chrom, pos)
        for gene in overlapping_genes:
            for transcript in gene.transcripts:
                if transcript.overlaps(pos, pos+len(ref_allele)):
                    effects = transcript.effects(ref_allele, alt_allele, pos)
                    for effect in effects:
                        print(f"{chrom}:{pos} {ref_allele}>{alt_allele} in {gene.gene_name} ({effect.effect_type})")

# Clean up resources
ref.close()