"""Generate and independently check single-net TAPM/linker design models."""
import hashlib
import json

from generate_crystal_templates import ROOT, DEST, pcb, Structure, CifWriter
from generate_matched_crystal_templates import validate_imine_cell


def audit_substitutions(designs):
    from cof_vision.cofdraw_science import import_cif
    from cof_vision.cofdraw_crystal_linkers import fluorinated_model
    manifest = json.loads((DEST / "manifest.json").read_bytes())
    results = {}
    for row in manifest["reported_models"] + manifest["curated_models"] + manifest["templates"] + designs:
        if row.get("network_dimension") != 3 and row["id"] != "dia-imine":
            continue
        parent = import_cif((DEST / row["path"]).read_text("utf-8"))
        parent["provenance"]["model_name"] = row["model_name"]
        try:
            variant = fluorinated_model(parent, row["id"])
            count = len(variant["provenance"]["substituted_atom_indices"])
        except ValueError:
            count = 0
        results[row["id"]] = {"parent_sha256": row["sha256"], "substituted_sites": count}
    (DEST / "aryl-substitutions.json").write_text(json.dumps(results, indent=2)+"\n", encoding="utf-8", newline="\n")


def main():
    catalog = json.loads((ROOT / "data/framework-drawing/catalog/current.json").read_bytes())
    monomers = {m["id"]: m for m in catalog["monomer_presets"]}
    directory = DEST / "diamond"
    directory.mkdir(exist_ok=True)
    entries = []
    for identifier, precursor, core in [
        ("bpda", "bpda_cho", "BPNY_CHO_H"),
        ("tfta", "tfta", "BENZ_CHO_F_F"),
        ("dhta", "dhta", "BENZ_CHO_OH"),
        ("dmta", "dmtp", "BENZ_CHO_OMe"),
        ("dbta", "cv_53ad3ea36de7ff89", "BENZ_CHO_Br"),
    ]:
        name = f"D4_TTPM_NH2_H-L2_{core}-DIA_A-1"
        framework = pcb.Framework(name, save_bb=False, log_level="error", out_path=str(ROOT / "build/model-generator"))
        structure = Structure(framework.cellMatrix, framework.atom_types, framework.atom_pos, coords_are_cartesian=True)
        text = "# COFdraw single-net unrelaxed imine design.\n" + str(CifWriter(structure, symprec=None))
        audit = validate_imine_cell(text, monomers[precursor]["smiles"], monomers["tapm"]["smiles"])
        assert audit["inferred_periodic_dimensions"] == [3]
        path = directory / f"{identifier}.cif"
        path.write_text(text, encoding="utf-8", newline="\n")
        label = monomers[precursor]["label"].split(" · ")[0].split(" — ")[0]
        entries.append(dict(id="design-dia-"+identifier, model_name=f"TAPM + {label}",
            label=f"TAPM + {label} · single-net design", path=path.relative_to(DEST).as_posix(),
            sha256=hashlib.sha256(path.read_bytes()).hexdigest(), network_dimension=3,
            monomer_ids=["tapm", precursor], generator=name, linker_label=label,
            description="Single-net dia imine starting model from pyCOFBuilder. Explicit atoms and attachment graph checked. Cell and packing are unrelaxed.", audit=audit))
    (directory / "models.json").write_text(json.dumps(entries, indent=2, ensure_ascii=False)+"\n", encoding="utf-8", newline="\n")
    audit_substitutions(entries)


if __name__ == "__main__":
    main()
