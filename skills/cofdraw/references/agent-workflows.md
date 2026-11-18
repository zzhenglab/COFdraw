# Query, calculate, and return inspectable results

The running API is authoritative. `GET /api/drawing/science/capabilities` returns current operation paths and actual JSON request schemas. The helper saves this response:

```shell
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/contracts science capabilities
```

Use `science query` with `record_types:["literature"]` to retrieve the bundled Yaghi review ledger by title or DOI. Each paper carries actually inspected pages, checked monomer IDs, unresolved candidates and review limits. There are 100 selected index records; this is not an exhaustive main-paper/SI audit. Monomer query support includes current runtime roles separately from the preserved source catalog record. A functional-group match does not validate a reported source reaction or periodic crystal.

## Describe a drawing before rendering

`GET /api/drawing/describe/catalog` supplies built-in description examples and exact selectable source, monomer, linkage, template, style and palette IDs. The helper exposes it as `describe-catalog`:

```shell
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/description-options describe-catalog
```

For local interpretation, save a request such as `{"text":"Draw Pyurea with building-block highlights and a soft palette","mode":"rules"}` in `description.json`:

```shell
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/description-plan describe --request description.json
```

The request goes to `POST /api/drawing/describe`. Its fields are `text` (1–4000 characters), `mode` (`rules` by default, or `llm`), optional `selection`, and optional `llm` settings. This call returns a plan in `result.json`; it does not draw. Read `status`, evidence, assumptions and unresolved items. A ready plan provides a source recipe ID or a bounded prepare request plus presentation choices. Continue with the existing source/prepare/render workflow only after resolving missing or conflicting constraints. Do not treat selected topology as experimental identity.

`selection` accepts only `source_recipe_id`, `monomer_a_id`, `monomer_b_id`, `linkage`, `template_id`, `style_id` and `palette`, all taken from the actual catalog. Source and new design chemistry are mutually exclusive. Manual adjustments are revalidated locally with `mode: "rules"`, with no provider call. Recognized chemistry/geometry conflicts and unsupported operations remain blocking. With explicit manual or LLM choices, unmatched residual words are surfaced as review warnings; a valid parameter set does not guarantee complete semantic fidelity. Compare the returned choices and warnings with the user's intent before continuing.

Optional LLM mode is for an explicitly requested provider call. It needs `llm: {"base_url":"https://provider.example/v1","model":"user-chosen-model","api_key":"..."}` and no manual `selection`. The provider sees the description and minimal catalog labels/IDs; it cannot submit SMILES, recipes, code or new chemistry. Its proposed IDs go through the deterministic planner, and any unresolved provider item prevents drawing readiness. HTTPS is required except for exact localhost HTTP hosts; redirects are blocked. Keep credentials in the desktop field or an in-memory request, never in a shared report or code sample. The helper does not copy a supplied request file or its credentials into output files and returns generic description HTTP errors to avoid credential echoes.

The application stores no API configuration and discovers no environment keys. Do not turn a rules request into a paid/API request automatically, or retry provider failures silently. For compatibility and credential details, read `docs/llm-api.md` in the repository. No provider call is needed for built-in examples, explicit ID selections, rendering or scientific calculations.

For conversational follow-ups, direct clients use `POST /api/drawing/describe/chat` with the new `text` (1–3000 characters), explicit mode and optional `context` returned as `conversation.context` by the previous turn. `operation` is `drawing`, `style_only`, `help` or `reset`. The existing helper `describe` command remains the single-turn endpoint; it does not implicitly start a conversation. Read [conversation contracts](../../../docs/conversation-methods.md) for a complete example.

Carry only catalog IDs, bounded unresolved constraints, recognized site/instance requirements and the optional opaque local revision token. Do not send a whole recipe, custom SMILES or coordinates as context. A `drawing` plan can prepare/render only when ready. A `style_only` plan has no prepare request: confirm its `context_token` still matches the local active recipe and apply the returned style options to that same recipe. This preserves custom monomers and manual coordinates. A stale plan must be analyzed again. `help` opens allowlisted local tools; it supplies no calculation values. Reset is explicit and does not erase the drawing.

History and credentials remain session-only. A provider receives the submitted text and catalog descriptors/selected IDs, not previous chat text, pending local notes, project files or local molecular data. Review the current card before applying and report success only after the actual renderer returns. Unsupported chemistry and property targets remain unresolved across follow-ups; changing color does not fulfill them.

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

If the user is viewing the current finite drawing, first use its actual atoms instead of asking for a CIF unnecessarily. Call `science patch` with the current recipe or registered source ID, then pass the returned `viewer` to `finite-inspect` or `finite-powder` as `model`. For example, after saving a real patch response:

```python
import json
from pathlib import Path

patch = json.loads(Path("outputs/patch/result.json").read_text(encoding="utf-8"))
Path("finite-check.json").write_text(json.dumps({"model": patch["viewer"]}), encoding="utf-8")
Path("finite-scattering.json").write_text(json.dumps({"model": patch["viewer"],
    "wavelength": 1.5406, "two_theta_min": 2, "two_theta_max": 50,
    "step": 0.1}), encoding="utf-8")
```

```shell
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/finite-check science finite-inspect --request finite-check.json
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/finite-scattering science finite-powder --request finite-scattering.json
```

These results describe displayed atoms and finite Debye scattering, with continuation sites omitted and missing H unfilled. They cannot supply density, PLD or hkl peaks. For a periodic model associated with a design, submit the unchanged `prepared.json` recipe to `periodic-options`; inspect the offered IDs and match qualifications. If one is chosen, send the same recipe in `model.recipe_match` so the service revalidates the precursor/linkage/template binding. The resulting cell has separate unrelaxed coordinates. Do not construct a box around a finite patch or substitute an unrelated template when no match exists.

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

To compare periodic stacking, call `science stacking` separately for AA, AB, ABC and slipped using the same supplied single planar layer, spacing and slip convention. Then use each returned `structure` for PXRD with identical settings. The layer is screened for robust 2D contact connectivity and thickness; 3D or multilayer cells may correctly be rejected. Twisted finite-layer illustrations use `patch` plus `finite-powder`, not periodic stacking. UI comparison profiles normalize each curve independently; retain absolute/relative output fields and do not infer phase fractions from their peak heights.

## Compare supplied structures

`POST /api/drawing/science/screen` takes 1–8 entries, each `{"name":"descriptive label","structure":{...}}`, plus `include_porosity` (default false), `probe_radius` (Å, default 1.82), and `samples` (default 4000). Labels identify supplied inputs; they do not establish literature identity. Use the same conventions and parameters across candidates. The helper exposes this as `science screen --request screen.json` and saves returned CSV/JSON.

Inspect per-entry errors and warnings before ranking anything. Formula, listed-atom density, distance-inferred contact dimensionality and optional geometric porosity do not alone establish adsorption performance, catalytic activity, stability, selectivity or synthesis feasibility. Missing hydrogens, solvent, disorder or incomplete cell content change mass-normalized results. BET requires separately supplied adsorption data, not this screening table.

## Molecular accessibility and adsorption

For precursor complexity screening, `science synthesis-library` returns the actual finite-monomer RDKit SA Scores and qualifications. Rendered drawings carry `audit.precursor_synthetic_accessibility` for their declared precursors. Compare the stated 1–10 heuristic only; do not interpret it as a synthesis probability or framework feasibility score.

For transport geometry, use `science geometry --request geometry.json` with an explicit periodic structure and common probe/grid settings. Retain PLD/LCD bounds and `connectivity.status`; an unresolved zero cannot rank a structure as nonporous. For physical adsorption, use the separate [asynchronous RASPA workflow](science-api.md#local-adsorption-jobs), with actual model/charge/force-field choices and convergence checks. The helper saves the returned ZIP of simulation inputs and output as well as JSON/CSV; a terminated partial job must be labeled accordingly.

## Query saved calculations back

Scientific reports that provide `report_id` are stored under the service's local runtime data directory. They include operation, creation time, engine version, input parameters and input SHA-256; individual methods add their structure hash, units and qualifications. Use:

- `GET /api/drawing/science/reports?limit=20` for recent report metadata.
- `GET /api/drawing/science/reports/{report_id}` for the full saved report.
- A report's returned `downloads` paths for actual CSV, JSON or structure files.

These are local paths on the same loopback service, not public sharing URLs. A report read does not rerun the calculation. Keep IDs exactly as returned; do not guess IDs or filenames. Catalog queries and drawing measurements may return results without report IDs; preserve their `result.json` directly.

The helper creates a fresh isolated runtime for each call. Its `result.json`, `structure.json` when returned, `artifacts.json`, and downloaded files remain available after it stops. To query an earlier run through HTTP, start a headless service using that run's existing `runtime-data` directory and read its status URL, as described in `local-api.md`. Do not start a fresh empty runtime and expect another runtime's report history. Agents can also read the saved JSON files directly.

For a user-facing handoff, link the actual report and exports, state the concrete finding, and retain the model/source status plus material limitations. For example: “This supplied model has a 2D distance-inferred contact network at the reported cutoff; its ideal PXRD CSV and the full input-bound JSON report are attached.” Use the calculated values only after inspecting the result.
