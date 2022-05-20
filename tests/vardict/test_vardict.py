#!/usr/bin/env python
import pytest  # type: ignore
import os 
from typer.testing import CliRunner
from pdb import set_trace as bp
from pv import app

runner = CliRunner()
vardict_single_calls = [
        ['vardict', 'single', 'filter', '--inputVcf', 'tests/vardict/data/single_test.vcf' , '--tsampleName', 'Myeloid200-1', '-ad','1', '-o', 'tests/vardict/data/single'], 

]

vardict_matched = [
    ['vardict', 'case-control', 'filter', '--inputVcf', 'tests/vardict/data/case_control_test.vcf' , '--tsampleName', 'C-C1V52M-L001-d', '-ad','1' , '-o', 'tests/vardict/data/two']
]

@pytest.mark.parametrize("call", vardict_single_calls)
def test_single(call):
    # vardict_tests['single']
    result = runner.invoke(app, call)
    assert result.exit_code == 0
    assert '' in result.stdout
    assert os.path.exists("tests/vardict/data/single/single_test_complex_STDfilter.vcf") == True
    assert os.path.exists("tests/vardict/data/single/single_test_STDfilter.vcf") == True
    assert os.path.exists("tests/vardict/data/single/single_test_STDfilter.txt") == True
    os.remove("tests/vardict/data/single/single_test_complex_STDfilter.vcf")
    os.remove("tests/vardict/data/single/single_test_STDfilter.vcf")
    os.remove("tests/vardict/data/single/single_test_STDfilter.txt")

@pytest.mark.parametrize("call", vardict_matched)
def test_two(call):
    result = runner.invoke(app, call)
    assert result.exit_code == 0
    assert '' in result.stdout
    assert os.path.exists("tests/vardict/data/two/case_control_test_complex_STDfilter.vcf") == True
    assert os.path.exists("tests/vardict/data/two/case_control_test_STDfilter.vcf") == True
    assert os.path.exists("tests/vardict/data/two/case_control_test_STDfilter.txt") == True
    os.remove("tests/vardict/data/two/case_control_test_complex_STDfilter.vcf")
    os.remove("tests/vardict/data/two/case_control_test_STDfilter.vcf")
    os.remove("tests/vardict/data/two/case_control_test_STDfilter.txt")


