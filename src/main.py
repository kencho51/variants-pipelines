from interpret_variants import interpret_variants

ref_fasta = "app/tests/data/reference.fasta"
input_vcf = "app/tests/data/input.vcf"

results = interpret_variants(ref_fasta, input_vcf)

for result in results:
    print(result)