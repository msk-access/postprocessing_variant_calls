import pytest  # type: ignore
import os
from typer.testing import CliRunner
from pdb import set_trace as bp
from postprocessing_variant_calls.main import app

runner = CliRunner()
maf_concat_files = [
    [
        "maf",
        "concat",
        "-f",
        "tests/data/maf/concat/maf2.maf",
        "-f",
        "tests/data/maf/concat/maf1.maf",
        "-o",
        "tests/data/maf/concat/output_maf.maf",
        "-h",
        "postprocessing_variant_calls/resources/maf_concat/header.txt",
    ]
]

maf_concat_paths = [
    [
        "maf",
        "concat",
        "-p",
        "tests/data/maf/concat/paths.txt",
        "-o",
        "tests/data/maf/concat/output_maf.maf",
        "-h",
        "postprocessing_variant_calls/resources/maf_concat/header.txt",
    ]
]

maf_annotate_maf_by_bed = [["maf", "annotate", "mafbybed", "--help"]]

maf_annotate_maf_by_bed = [["maf", "annotate", "mafbybed", "--help"]]

maf_subset = [
    [
        "maf",
        "subset",
        "-m",
        "tests/data/maf/subset/example_input.maf",
        "--ids",
        "postprocessing_variant_calls/resources/maf_subset/example_subset_ids.txt",
        "-o",
        "tests/data/maf/subset/example_output.maf",
    ]
]

maf_filter = [
    ["maf", "filter", "--help"],
    ["maf", "filter", "cmo_ch", "--help"],
    ["maf", "filter", "hotspot", "--help"],
    ["maf", "filter", "mappable", "--help"],
    ["maf", "filter", "non_common_variant", "--help"],
    ["maf", "filter", "non_hotspot", "--help"],
    ["maf", "filter", "not_complex", "--help"],
]


maf_tag = [
    ["maf", "tag", "--help"],
    ["maf", "tag", "cmo_ch", "--help"],
    ["maf", "tag", "common_variant", "--help"],
    ["maf", "tag", "germline_status", "--help"],
    ["maf", "tag", "prevalence_in_cosmicDB", "--help"],
    ["maf", "tag", "truncating_mut_in_TSG", "--help"],
    ["maf", "tag", "traceback", "--help"],
    ["maf", "tag", "by_variant_classification", "--help"],
]


# Data Driven Tests
@pytest.mark.parametrize("call", maf_concat_files)
def test_concat_files(call):
    result = runner.invoke(app, call)
    result.stdout
    assert result.exit_code == 0
    assert "" in result.stdout
    assert os.path.exists("tests/data/maf/concat/output_maf.maf") == True
    os.remove("tests/data/maf/concat/output_maf.maf")


@pytest.mark.parametrize("call", maf_concat_paths)
def test_concat_paths(call):
    result = runner.invoke(app, call)
    assert result.exit_code == 0
    assert "" in result.stdout
    assert os.path.exists("tests/data/maf/concat/output_maf.maf") == True
    os.remove("tests/data/maf/concat/output_maf.maf")


@pytest.mark.parametrize("call", maf_subset)
def test_concat_files(call):
    result = runner.invoke(app, call)
    result.stdout
    assert result.exit_code == 0
    assert "" in result.stdout
    assert os.path.exists("tests/data/maf/subset/example_output.maf") == True
    os.remove("tests/data/maf/subset/example_output.maf")


# Command Test
@pytest.mark.parametrize("call", maf_annotate_maf_by_bed)
def test_annotate_mafbybed(call):
    result = runner.invoke(app, call)
    assert result.exit_code == 0


@pytest.mark.parametrize("call", maf_tag)
def test_maf_tag(call):
    result = runner.invoke(app, call)
    assert result.exit_code == 0


@pytest.mark.parametrize("call", maf_filter)
def test_maf_filter(call):
    result = runner.invoke(app, call)
    assert result.exit_code == 0
