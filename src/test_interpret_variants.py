import io
import pytest
import pysam
import pyensembl
from unittest.mock import patch

# Import the function you want to test from the script
from interpret_variants import interpret_variants

def test_interpret_variants():
    # Define input data
    ref_fasta = "tests/data/reference.fasta"
    input_vcf = "tests/data/input.vcf"
    mock_input = f"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\nchr1\t1000\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n".encode('utf-8')

    # Define expected output
    expected_output = ["chr1:1000 A>G in gene1 (missense_variant)",
                       "chr2:2000 T>C in gene2 (synonymous_variant)"]

    # Patch the necessary objects
    with patch("builtins.open", return_value=io.StringIO(mock_input.decode('utf-8'))):
        with patch.object(pysam, "FastaFile", return_value=None) as mock_ref:
            with patch.object(pyensembl, "EnsemblRelease", return_value=None) as mock_ensembl:
                with patch("builtins.print") as mock_print:
                    # Call the function with the mock input data
                    interpret_variants(ref_fasta, input_vcf)

                    # Assert that the mock objects were called as expected
                    mock_ref.assert_called_once_with(ref_fasta)
                    mock_ensembl.assert_called_once_with(release=98)
                    mock_print.assert_called_once_with(expected_output)
