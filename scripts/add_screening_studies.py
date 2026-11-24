"""Join completed adsorption studies to exact CIFs in a COFdraw explorer snapshot."""
import argparse
import csv
import gzip
import json
import math
import os
import statistics
import tempfile
from collections import defaultdict
from pathlib import Path


def records(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        yield from csv.DictReader(stream)


def enrich(snapshot, tables, native):
    by_id = {r["id"]: r for r in snapshot["records"]}
    metrics = snapshot["metrics"]

    def add(item, key, label, raw):
        if raw in (None, ""):
            return
        rid = item.get("parent") or item["id"]
        row = by_id[rid]
        if row.get("sha256") != (item.get("cif_sha256") or item["sha256"]):
            raise ValueError(f"CIF mismatch for {rid}")
        value = float(raw)
        if not math.isfinite(value) or key in row["values"]:
            raise ValueError(f"Invalid or duplicate {key} for {rid}")
        row["values"][key] = value
        metrics[key] = label

    for r in records(tables / "isotherms.csv"):
        if int(r["n_seeds"]) == 3:
            bar = float(r["pressure_bar"])
            if bar not in (1, 5.8, 10, 20, 35, 50, 65):
                raise ValueError("Unexpected methane isotherm pressure")
            add(r, "ch4_q_" + f"{bar:g}".replace(".", "_"), f"CH₄ uptake at {bar:g} bar (mmol/g)", r["qCH4_mmol_g"])
    for r in records(tables / "novel-storage.csv"):
        add(r, "additional_delivery", "Additional panel · CH₄ delivery (mmol/g)", r["delivery_mmol_g"])
        add(r, "additional_delivery_vv", "Additional panel · CH₄ delivery (cm³ STP/cm³)", r["delivery_vv"])
    for r in records(tables / "temperature-summary.csv"):
        gas = r["gas"]
        if gas not in ("CH4", "H2"):
            raise ValueError("Unexpected Widom guest")
        for t in (278, 298, 318):
            add(r, f"henry_{gas}_{t}", f"{gas} Henry coefficient at {t} K (mol/kg/Pa)", r[f"KH_{t}K_mol_kg_Pa"])
        if r["heat_status"] == "resolved" and int(r["observations"]) == 9:
            add(r, f"heat_{gas}", f"{gas} apparent dilute adsorption heat (kJ/mol)", r["apparent_heat_kJ_mol"])
    for r in records(tables / "dilute-selectivity.csv"):
        t = int(r["temperature_K"])
        if t not in (278, 298, 318):
            raise ValueError("Unexpected Widom temperature")
        if int(r["paired_resolved_seeds"]) == 3:
            add(r, f"henry_ratio_{t}", f"CH₄/H₂ Henry ratio at {t} K · triplicate", r["dilute_CH4_H2"])

    profiles = {"baseline": "02af8c3ce9e3ccc8b8754e0ad78eb952bbb0bdb014db0a4322b1f7dea1b43eb7",
                "uff-lb": "8e86953aa9a926445f49f688b952b487ff7999b09d3deb71b53f47ee50d9ac26"}
    seeds = (120071, 120121, 120181)
    paired = defaultdict(dict)
    for r in records(native):
        rid, profile = r["id"], r["profile"]
        pressure, seed = float(r["bar"]), int(r["seed"])
        if profiles.get(profile) != r["profile_sha256"] or seed not in seeds:
            raise ValueError("Force-field run protocol mismatch")
        if pressure == 35 and rid in {"COF-1","COF-5","COF-6","COF-8","COF-10","COF-102","COF-103"}:
            continue
        group = paired[rid, r["cif_sha256"]]
        point = (profile, seed, pressure)
        q = float(r["q_mol_kg"])
        if point in group or pressure not in (5.8, 65) or not math.isfinite(q) or by_id[rid].get("sha256") != r["cif_sha256"]:
            raise ValueError(f"Invalid force-field observation for {rid}")
        group[point] = q
    for (rid, cif), group in paired.items():
        if len(group) == 12:
            ref = {"id": rid, "cif_sha256": cif}
            for profile in profiles:
                mean = statistics.mean(group[profile,s,65] - group[profile,s,5.8] for s in seeds)
                add(ref, "ff_" + profile.replace("-", "_"), f"{profile} host · matched CH₄ delivery (mmol/g)", mean)
    snapshot["methods"].update(
        isotherm="Absolute methane uptake at 298 K; three seeds per pressure, 10,000 initialization and 25,000 production cycles.",
        additional_delivery="Separate 298 K endpoint panel: three seeds, 2,500 initialization and 5,000 production cycles per pressure.",
        thermal="Rigid neutral host. 15,000 Widom cycles per independent gas/temperature/seed at 278, 298 and 318 K. Three seeds per published coefficient. Apparent dilute heats use resolved three-temperature fits.",
        paired_force_fields="Absolute methane delivery, 65 minus 5.8 bar, at 298 K in the same rigid cell. Three paired seeds and both host profiles required. 20,000 initialization and 50,000 production cycles, identical methane, shifted 12 Å cutoff and Lorentz–Berthelot mixing.",
        host_profile_sha256=profiles)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("input", "tables", "native", "output"):
        parser.add_argument("--" + key, type=Path, required=True)
    args = parser.parse_args()
    with gzip.open(args.input, "rt", encoding="utf8") as stream:
        snapshot = json.load(stream)
    enrich(snapshot, args.tables, args.native)
    with tempfile.NamedTemporaryFile(dir=args.output.parent, delete=False) as temp:
        path = Path(temp.name)
    try:
        with gzip.open(path, "wt", encoding="utf8", compresslevel=5) as stream:
            json.dump(snapshot, stream, ensure_ascii=False, allow_nan=False, separators=(",",":"))
        os.replace(path, args.output)
    finally:
        path.unlink(missing_ok=True)
