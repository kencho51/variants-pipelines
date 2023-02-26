import io
import pytest
from unittest.mock import patch

from mymodule import process_vcf

def test_process_vcf():
    # Define input data
    ref_fasta = "tests/data/reference.fasta"
    input_vcf = "tests/data/input.vcf"

    # Define expected output
    expected_output = ["chr1:1000 A>G in gene1 (missense_variant)",
                       "chr2:2000 T>C in gene2 (synonymous_variant)"]

    # Mock Ensembl genome annotation object
    ensembl_mock = MagicMock()
    gene1_mock = MagicMock()
    gene1_mock.gene_name = "gene1"
    gene1_mock.transcripts = [MagicMock(overlaps=lambda x, y: True, effects=lambda x, y, z: [MagicMock(effect_type="missense_variant")])]
    gene2_mock = MagicMock()
    gene2_mock.gene_name = "gene2"
    gene2_mock.transcripts = [MagicMock(overlaps=lambda x, y: True, effects=lambda x, y, z: [MagicMock(effect_type="synonymous_variant")])]
    ensembl_mock.genes_at_locus.return_value = [gene1_mock, gene2_mock]

    # Mock pysam.FastaFile object
    ref_mock = MagicMock()

    # Mock open() function to return VCF data
    vcf_data = "#CHROM\tPOS\tREF\tALT\nchr1\t1000\tA\tG\nchr2\t2000\tT\tC\n"
    open_mock = mock_open(read_data=vcf_data)

    # Execute function
    with patch("pysam.FastaFile", return_value=ref_mock), \
         patch("pyensembl.EnsemblRelease", return_value=ensembl_mock), \
         patch("builtins.open", open_mock):
        output = io.StringIO()
        with redirect_stdout(output):
            process_vcf(ref_fasta, input_vcf)

    # Check if expected output is produced
    assert output.getvalue().strip().split("\n") == expected_output
