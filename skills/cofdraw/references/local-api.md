# Local runtime and API

## Runtime

From a source checkout with Python 3.11 or newer, install the project's dependencies in a virtual environment; follow the repository README for the supported platform. Invoke the helper with that environment's Python. Source mode runs `scripts/launch_cofdraw.py` using the same interpreter. For the tested Windows dependency versions, use `requirements-cofdraw-desktop.txt`.

For an existing Windows release, extract the entire archive and pass the full `COFdraw.exe` path to `--exe`. The helper needs only standard-library Python in this mode. Native desktop viewing needs WebView2; headless API use does not open a webview.

By default, the helper reads template data from `<repo>/data` or the executable's sibling `data`. Override this with `--data-dir` when using another reviewed portable catalog. It copies `framework-drawing` into `<output>/runtime-data` and excludes generated `exports` and log files. It never changes the template catalog.

Direct launch, if a persistent service is useful:

```shell
python scripts/launch_cofdraw.py --headless --data-dir outputs/runtime-data --status-file outputs/server.json
```

First copy the reviewed `data/framework-drawing` into `outputs/runtime-data/framework-drawing` and create `outputs`. The status file reports the chosen loopback URL and process ID. Wait for `GET /api/desktop` to succeed; the initial status file alone does not establish readiness. Stop only the process you started when finished. The API binds to `127.0.0.1` and rejects cross-origin writes; it is a local application, not a hosted multiuser service.

## Endpoint contract

| Request | Result |
| --- | --- |
| `GET /api/desktop` | Version, offline flag, runtime data directory |
| `GET /api/drawing/library` | Source recipes, styles, monomers, templates, linkages and registry errors |
| `GET /api/drawing/sources/{recipe_id}` | Hash-bound source recipe and qualified metadata |
| `POST /api/drawing/monomer` | Validated monomer graph, SMILES, MOL block and depiction |
| `POST /api/drawing/prepare` | Editable finite recipe for explicit monomer inputs |
| `POST /api/drawing/render` | Audit and download URLs for actual generated files |
| `GET /api/drawing/exports/{drawing_id}/{extension}` | A generated SVG, PNG, MOL, SDF, CDXML or JSON file, when available |

Read valid IDs and supported status from `library` rather than maintaining a separate invented catalog. Requests reject unknown top-level fields. HTTP 422 carries the invalid-input reason.

### Source drawing

Render request:

```json
{"source_recipe_id":"tppa-pore-v1","style_id":"standard","size":1200}
```

`render` requires **exactly one** of `source_recipe_id` and `recipe`. Size is 500–2400 pixels. Optional fields are `style_seed` (0–2147483647) and `style_options` (object). Use the controls exposed by the library and the UI as the style-options reference. Retain the seed and resolved options for random styles.

### Custom design

Write this as `design.json` to create an explicitly chosen 3c+2c imine patch; the template is a design choice, not experimental topology:

```json
{
  "monomer_a":"O=Cc1cc(C=O)cc(C=O)c1",
  "monomer_b":"Nc1ccc(N)cc1",
  "label_a":"TFB",
  "label_b":"Pa",
  "linkage":"imine",
  "template_id":"hcb_3_2",
  "instances":12
}
```

```shell
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/design design --request design.json --style standard --size 1200
```

The helper sends this JSON to `/prepare`, saves `prepared.json`, then renders the returned `recipe`. Monomer inputs are finite, connected SMILES. Linkage support is bounded to the implemented transforms and compatible reactive groups. A reference-only linkage or unsupported template must not be presented as a working constructor.

For an edited recipe, prepare a full render request containing `{"recipe": {...}, "style_id": "standard", "size": 1200}` and pass it to the helper's `render --request render.json` subcommand. The API strips inherited source identity from all editable recipes. Render the unchanged registered source by ID when source identity should remain attached.

## Output and checks

Each helper invocation requires a new output directory. `library.json` preserves catalog selection context. A render creates `result.json`, `drawing.svg`, `drawing.png`, `drawing.mol`, and any returned optional formats, including the JSON recipe/audit bundle. `artifacts.json` records file sizes and SHA-256 hashes. `server-output.log` and `runtime-data/cofdraw.log` are local diagnostics.

The response `audit` records the graph scope, drawing status, monomers, template, style, and chemical export diagnostics. Check this before describing an output. CDXML export may be unavailable; only advertise formats present in `downloads`. A successful API response is a software check, not independent verification of the chemistry.
