import copy
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from cof_vision import cofdraw_raspa as adapter
from cof_vision.cofdraw_raspa_worker import _cell_parameters
from cof_vision.cofdraw_science import import_cif

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture
def structure():
    return import_cif((REPO/"data/crystal-models/dia-imine.cif").read_text())


def prepare(tmp_path, structure, **overrides):
    options = dict(directory=tmp_path, structure=structure, gases=["CO2"], mole_fractions=[1.],
                   temperature_K=298., pressure_Pa=1e5, cycles=100, init_cycles=20,
                   seed=91357, charge_mode="zero", mode="gcmc")
    options.update(overrides)
    return adapter.prepare_run(**options)


@pytest.mark.parametrize("change,match", [
    ({"gases":["water"]}, "supported gases"),
    ({"gases":["CO2","CO2"],"mole_fractions":[.5,.5]}, "distinct"),
    ({"mole_fractions":[.8]}, "sum to one"),
    ({"mole_fractions":[float("nan")]}, "finite"),
    ({"cycles":50}, "cycles"),
    ({"seed":True}, "finite"),
    ({"temperature_K":18}, "temperature_K"),
    ({"charge_mode":"automatic"}, "explicit"),
    ({"charge_mode":"cif"}, "every framework atom"),
    ({"fugacity_coefficients":[0.]}, "fugacity coefficient"),
    ({"mode":"gcmc_fake"}, "mode"),
])
def test_input_constraints_before_any_job_write(tmp_path, structure, change, match):
    with pytest.raises(ValueError, match=match):
        prepare(tmp_path/"job", structure, **change)
    assert not (tmp_path/"job").exists()


def test_charges_are_explicit_and_neutral(tmp_path, structure):
    for atom in structure["atoms"]:
        atom["charge"] = 0.
    structure["atoms"][0]["charge"] = .2
    with pytest.raises(ValueError, match="neutral framework"):
        prepare(tmp_path/"bad",structure,charge_mode="cif")
    structure["atoms"][1]["charge"] = -.2
    manifest=prepare(tmp_path/"good",structure,charge_mode="cif")
    saved=json.loads((tmp_path/"good/engine-input.json").read_text())
    assert saved["structure"]["atoms"][0]["charge"] == .2
    assert manifest["framework_net_charge_e"] == 0


def test_periodic_padding_uses_cell_heights_and_preserves_input(tmp_path):
    structure={"cell":[[8.,0.,0.],[3.,9.,0.],[0.,0.,5.]],"atoms":[{"element":"C","fractional":[.2,.3,.4]}]}
    original=copy.deepcopy(structure)
    result=prepare(tmp_path,structure)
    assert structure == original
    assert min(result["supercell_heights_angstrom"]) >= 24
    assert result["simulation_atom_count"] == math.prod(result["repeats"])
    assert result["repeats"] == [4,3,5]


def test_mixture_preserves_physical_partial_fugacities(tmp_path,structure):
    result=prepare(tmp_path,structure,gases=["CO2","N2"],mole_fractions=[.15,.85],fugacity_coefficients=[.9,1.1])
    assert result["partial_fugacities_Pa"] == pytest.approx([13500.,93500.])
    assert result["internal_fugacity_coefficients"] == pytest.approx([.27,1.87])
    assert result["mole_fractions"] == [.15,.85]


def test_collisions_and_unknown_framework_parameters_rejected(tmp_path):
    base={"cell":[[30.,0,0],[0,30.,0],[0,0,30.]],"atoms":[{"element":"C","fractional":[0,0,0]},{"element":"C","fractional":[.01,0,0]}]}
    with pytest.raises(ValueError,match="short contacts"):
        prepare(tmp_path/"collision",base)
    base["atoms"]=[{"element":"U","fractional":[0,0,0]}]
    with pytest.raises(ValueError,match="No reviewed"):
        prepare(tmp_path/"unsupported",base)


def test_engine_resource_status_and_input_tamper(tmp_path,structure):
    assert adapter.engine_status()["available"]
    manifest=prepare(tmp_path,structure)
    resources={row["path"]:row["sha256"] for row in manifest["resource_files"]}
    assert set(resources)=={"force-field-profile.json","RASPA2-LICENSE.txt","RASPA3-LICENSE.txt"}
    assert resources["force-field-profile.json"]==manifest["profile_sha256"]
    assert (tmp_path/"force-field-profile.json").read_bytes()==(adapter.resource_directory()/"profile.json").read_bytes()
    for name,digest in resources.items():
        assert adapter._sha(tmp_path/name)==digest
    saved_input=json.loads((tmp_path/"engine-input.json").read_text())
    assert saved_input["manifest"]["resource_files"]==manifest["resource_files"]
    (tmp_path/"engine-input.json").write_text("{}")
    with pytest.raises(ValueError,match="changed"):
        adapter.parse_run(tmp_path,manifest)


def test_native_cell_angle_units():
    raspa=pytest.importorskip("raspalib")
    cell=np.array([[20.,0,0],[7.,18.,0],[2.,4.,25.]])
    native=raspa.SimulationBox(*_cell_parameters(cell))
    assert native.volume == pytest.approx(np.linalg.det(cell),rel=1e-12)


def test_native_ideal_single_and_binary_fugacity_with_uncertainty(tmp_path):
    report=tmp_path/"ideal.json"
    process=subprocess.run([sys.executable,str(Path(__file__).with_name("raspa_analytic_probe.py")),str(report)],cwd=tmp_path,capture_output=True,timeout=45)
    assert process.returncode==0,process.stderr.decode(errors="replace")
    results=json.loads(report.read_text())
    for experiment in results:
        for row in experiment:
            # Native block uncertainty must cover a known physical equation;
            # there are no LJ forces or fitted adsorption values in this test.
            assert abs(row["mean_molecules"]-row["ideal_gas_PV_kT"]) <= row["confidence95"]


def run_job(tmp_path,structure,**kwargs):
    manifest=prepare(tmp_path,structure,**kwargs)
    process=subprocess.run(adapter.command_for(tmp_path),capture_output=True,timeout=45)
    assert process.returncode==0, process.stderr.decode(errors="replace")
    return adapter.parse_run(tmp_path,manifest)


def test_real_five_gas_widom_is_reproducible(tmp_path,structure):
    options=dict(gases=["CO2","N2","CH4","H2","Ar"],mole_fractions=[.2]*5,mode="widom")
    first=run_job(tmp_path/"a",structure,**options)
    second=run_job(tmp_path/"b",structure,**options)
    assert first["components"] == second["components"]
    assert len(first["components"]) == 5
    output_files={row["path"]:row["sha256"] for row in first["output_files"]}
    for record in first["manifest"]["resource_files"]:
        assert output_files[record["path"]]==record["sha256"]
    for row in first["components"]:
        assert 0 < row["henry_mol_kg_Pa"] < 1
        assert row["henry_confidence95_mol_kg_Pa"] >= 0
        assert "loading_mol_kg" not in row
    assert first["simulation_volume_angstrom3"] == pytest.approx(structure["cell_volume_angstrom3"]*math.prod(first["manifest"]["repeats"]))


@pytest.mark.parametrize("gases,fractions", [(["CO2","N2"],[.15,.85]),(["CO2","CH4"],[.5,.5]),(["H2"],[1.]),(["Ar"],[1.])])
def test_real_gcmc_and_mixtures(tmp_path,structure,gases,fractions):
    result=run_job(tmp_path,structure,gases=gases,mole_fractions=fractions,cycles=200)
    assert result["engine_version"] == "3.1.0"
    assert [r["gas"] for r in result["components"]] == gases
    assert all(r["loading_mol_kg"] >= 0 for r in result["components"])
    assert all("henry_mol_kg_Pa" not in r for r in result["components"])
    assert (tmp_path/"loading-statistics.txt").is_file()


def test_real_nonzero_cif_charges(tmp_path,structure):
    for atom in structure["atoms"]:
        atom["charge"]=0.
    structure["atoms"][0]["charge"]=.1
    structure["atoms"][1]["charge"]=-.1
    result=run_job(tmp_path/"charged",structure,charge_mode="cif",mode="widom")
    zero=run_job(tmp_path/"zero",structure,charge_mode="zero",mode="widom")
    assert result["manifest"]["charge_mode"] == "cif"
    assert result["components"][0]["henry_mol_kg_Pa"] > 0
    assert result["components"][0]["henry_mol_kg_Pa"] != zero["components"][0]["henry_mol_kg_Pa"]
