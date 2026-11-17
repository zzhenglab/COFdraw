"""Repeatable local drawing/export timing, never an OCR accuracy benchmark."""
import argparse
import hashlib
import inspect
import json
from pathlib import Path
import statistics
import sys
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    sys.path.insert(0, str(args.source_root.resolve()))
    from rdkit import rdBase
    from cof_vision.framework_drawing import VERSION, build_patch
    from cof_vision.cofdraw_chemistry import prepare_template
    from cof_vision.cofdraw_core import render_styled, export_chemical_formats

    args.output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for index in range(args.repeats):
        start = time.perf_counter()
        prepared = prepare_template(
            "COc1cc(C=O)c(OC)cc1C=O",
            "Nc1ccc(-c2cc(-c3ccc(N)cc3)cc(-c3ccc(N)cc3)c2)cc1",
            "imine", "DMTP", "TAPB", 12, "hcb_3_2",
        )
        recipe = {"schema_version": VERSION, "recipe_id": f"dmtp-{index}",
                  "layout_refinement": "rigid_block_translation_v1", **prepared}
        prepare_s = time.perf_counter() - start
        drawing_start = time.perf_counter()
        patch = build_patch(recipe) if "prepared_patch" in inspect.signature(render_styled).parameters else None
        keywords = {"prepared_patch": patch} if patch is not None else {}
        audit = render_styled(recipe, args.output_dir, 1000, "standard", **keywords)
        drawing_s = time.perf_counter() - drawing_start
        export_start = time.perf_counter()
        keywords = {"source_mol": patch[0]} if patch is not None and "source_mol" in inspect.signature(export_chemical_formats).parameters else {}
        formats = export_chemical_formats(args.output_dir, recipe["recipe_id"], **keywords)
        export_s = time.perf_counter() - export_start
        record = {"iteration": index + 1, "prepare_seconds": prepare_s,
                  "drawing_seconds": drawing_s, "export_seconds": export_s,
                  "total_seconds": time.perf_counter() - start, "formats": formats,
                  "atoms": audit["atoms"], "ports": len(audit["continuation_ports"]),
                  "fragment_serialization_sha256": audit["fragment_serialization_sha256"],
                  "png_sha256": hashlib.sha256((args.output_dir / f"{recipe['recipe_id']}.png").read_bytes()).hexdigest(),
                  "svg_sha256": hashlib.sha256((args.output_dir / f"{recipe['recipe_id']}.svg").read_bytes()).hexdigest()}
        records.append(record)
        print(json.dumps(record), flush=True)
    sources = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
               for p in (args.source_root / "cof_vision").glob("*.py")}
    report = {"scope": "Local 240-atom DMTP-TAPB hcb_3_2 finite patch, 1000px Standard, prepare + draw + SDF/CDXML export; excludes app startup and browser download.",
              "rdkit": rdBase.rdkitVersion, "records": records, "source_sha256": sources,
              "median_total_seconds": statistics.median(r["total_seconds"] for r in records),
              "model_inference": False, "training": False}
    (args.output_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
