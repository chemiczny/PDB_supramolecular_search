from os.path import dirname, join

import pytest

DATA_DIR = join(dirname(__file__), "data")
CIF_1BP0 = join(DATA_DIR, "1bp0.cif")


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    # The package reads config.json and writes scratch/log files relative to the cwd.
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_primitive_cif2dict_matches_biopython(workdir):
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from supramolecular_search.primitive_cif2dict import primitiveCif2Dict

    keys = ["_refine.ls_d_res_high", "_exptl.method"]
    parsed = primitiveCif2Dict(CIF_1BP0, keys).result
    reference = MMCIF2Dict(CIF_1BP0)

    for key in keys:
        assert parsed[key] == reference[key]


def test_anion_templates_are_packaged(workdir):
    from supramolecular_search.anion_recogniser import getAllTemplates

    templates = getAllTemplates()

    assert {"O", "N", "C", "S", "F", "CL", "BR", "I"} <= set(templates)


def test_find_supramolecular_1bp0(workdir):
    from supramolecular_search import cif_analyser

    cif_analyser.findSupramolecular((CIF_1BP0, "1BP0", "test"))

    expected_rows = {
        "anionPi": 20,
        "planarAnionPi": 15,
        "linearAnionPi": 0,
        "cationPi": 20,
        "piPi": 16,
        "anionCation": 3,
        "hBonds": 21,
        "metalLigand": 2,
        "methylPi": 88,
    }
    scratch = workdir / cif_analyser.config["scratch"]
    for log_name, rows in expected_rows.items():
        log_file = scratch / f"{log_name}test.log"
        lines = log_file.read_text().splitlines() if log_file.exists() else []
        assert len(lines) == rows, log_name
