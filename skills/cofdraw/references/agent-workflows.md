# Query, calculate, and return inspectable results

The running API is authoritative. `GET /api/drawing/science/capabilities` returns current operation paths and actual JSON request schemas. The helper saves this response:

```shell
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/contracts science capabilities
```

## Retrieve actual local chemistry records

`POST /api/drawing/science/query` accepts:

```json
{"query":"Pyurea COF","record_types":["source"],"limit":10,"offset":0}
```

Save it as `query.json`, then run:

```shell
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/query science query --request query.json
```

Available record types are `monomer`, `linkage`, `topology`, `source`, `crystal`, `template`, `style`, `preset`, and `capability`. Omit the filter to search all types. An empty query lists records, with deterministic pagination. Query length is at most 500 characters; limit is 1–50 and offset 0–10000. Names, punctuation variants, Greek alpha/beta/gamma and DOI URLs are normalized. Every meaningful keyword must match; this is local record retrieval, not general language understanding or web search.

Each result has a stable namespaced `id`, original `record_id`, `record_type`, label, support status, provenance/hashes, qualifications and exact record data. `match.fields` identifies matching fields. `next_offset` continues a list. Source records also supply a GET path for the actual hash-bound recipe. Check `errors`: an invalid source/model binding is excluded rather than silently trusted.

Examples:

- `{"query":"10.1039/d1qm00231g","record_types":["source"]}` retrieves the source-qualified Pyurea entry and its monomers.
- `{"query":"boronate ester","record_types":["linkage"]}` finds reference chemistry whose automatic constructor is unsupported. A manually admitted boronate source graph can still be renderable.
- `{"query":"hcb","record_types":["topology","template","crystal"]}` keeps a reference net, finite drawing templates and independent periodic models distinct.
- `{"query":"PXRD","record_types":["capability"]}` returns the implemented calculation's request schema and limitations.

Answer from these records. Do not fill missing material properties with invented values, promote a reference-only transform to executable chemistry, or treat missing results as proof that a compound does not exist. Abbreviations remain paper-scoped; inspect full names, exact SMILES and DOI qualifiers when names collide.

## A periodic-structure inspection workflow

1. Query `crystal` records or import the user's explicit CIF through `science model`. Keep its original source/provenance. Built-in templates are separate unrelaxed reticular designs.
2. Read the returned `structure` object from `result.json`. Submit that exact object to `/science/inspect`, optionally selecting `bond_scale` and `tolerance_angstrom`. Inspect formula, cell metrics, density, contact warnings, cutoff sensitivity and network dimensionality.
3. For a calculated powder pattern, submit the same structure to `/science/pxrd` with explicit wavelength, angular range, step, Gaussian FWHM and isotropic B value. Keep CSV/JSON and label the profile as ideal calculated PXRD.
4. If asked for pore geometry, submit the explicit structure to `/science/porosity`. Carry sample settings, radii, sampling intervals and connectivity limitations with the values.
5. Return linked files and a concise comparison of actual outputs, preserving units and unsuccessful checks. Never describe a numerical or graph check as proof of synthetic feasibility, stability or experimental identity.

For example, after an actual model job in `outputs/model`, create inspect/PXRD request files without manually retyping coordinates:

```python
import json
from pathlib import Path

result = json.loads(Path("outputs/model/result.json").read_text(encoding="utf-8"))
structure = result["structure"]
Path("inspect.json").write_text(json.dumps({"structure": structure,
    "bond_scale": 1.15, "tolerance_angstrom": 0.15}), encoding="utf-8")
Path("pxrd.json").write_text(json.dumps({"structure": structure,
    "wavelength": 1.5406, "two_theta_min": 2, "two_theta_max": 50,
    "step": 0.02, "fwhm": 0.15, "b_iso": 0}), encoding="utf-8")
```

```shell
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/inspect science inspect --request inspect.json
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/pxrd science pxrd --request pxrd.json
```

## Compare supplied structures

`POST /api/drawing/science/screen` takes 1–8 entries, each `{"name":"descriptive label","structure":{...}}`, plus `include_porosity` (default false), `probe_radius` (Å, default 1.82), and `samples` (default 4000). Labels identify supplied inputs; they do not establish literature identity. Use the same conventions and parameters across candidates. The helper exposes this as `science screen --request screen.json` and saves returned CSV/JSON.

Inspect per-entry errors and warnings before ranking anything. Formula, listed-atom density, distance-inferred contact dimensionality and optional geometric porosity do not alone establish adsorption performance, catalytic activity, stability, selectivity or synthesis feasibility. Missing hydrogens, solvent, disorder or incomplete cell content change mass-normalized results. BET requires separately supplied adsorption data, not this screening table.

## Query saved calculations back

Scientific reports that provide `report_id` are stored under the service's local runtime data directory. They include operation, creation time, engine version, input parameters and input SHA-256; individual methods add their structure hash, units and qualifications. Use:

- `GET /api/drawing/science/reports?limit=20` for recent report metadata.
- `GET /api/drawing/science/reports/{report_id}` for the full saved report.
- A report's returned `downloads` paths for actual CSV, JSON or structure files.

These are local paths on the same loopback service, not public sharing URLs. A report read does not rerun the calculation. Keep IDs exactly as returned; do not guess IDs or filenames. Catalog queries and drawing measurements may return results without report IDs; preserve their `result.json` directly.

The helper creates a fresh isolated runtime for each call. Its `result.json`, `structure.json` when returned, `artifacts.json`, and downloaded files remain available after it stops. To query an earlier run through HTTP, start a headless service using that run's existing `runtime-data` directory and read its status URL, as described in `local-api.md`. Do not start a fresh empty runtime and expect another runtime's report history. Agents can also read the saved JSON files directly.

For a user-facing handoff, link the actual report and exports, state the concrete finding, and retain the model/source status plus material limitations. For example: “This supplied model has a 2D distance-inferred contact network at the reported cutoff; its ideal PXRD CSV and the full input-bound JSON report are attached.” Use the calculated values only after inspecting the result.
