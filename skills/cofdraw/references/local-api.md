# Local runtime and API

## Runtime

From a source checkout with Python 3.11 or newer, install the project's dependencies in a virtual environment; follow the repository README for the supported platform. Invoke the helper with that environment's Python. Source mode runs `scripts/launch_cofdraw.py` using the same interpreter. For the tested Windows dependency versions, use `requirements-cofdraw-desktop.txt`.

For an existing Windows release, extract the entire archive and pass the full `COFdraw.exe` path to `--exe`. The helper needs only standard-library Python in this mode. Native desktop viewing needs WebView2; headless API use does not open a webview.

By default, the helper reads template data from `<repo>/data` or the executable's sibling `data`. Override this with `--data-dir` when using another reviewed portable catalog. It copies `framework-drawing` and, when present, `crystal-models` into `<output>/runtime-data` and excludes generated `exports` and log files. It never changes the template catalog.

Direct launch, if a persistent service is useful:

```shell
python scripts/launch_cofdraw.py --headless --data-dir outputs/runtime-data --status-file outputs/server.json
```

First copy the reviewed `data/framework-drawing` into `outputs/runtime-data/framework-drawing`, copy `data/crystal-models` if using periodic templates, and create `outputs`. The status file reports the chosen loopback URL and process ID. Wait for `GET /api/desktop` to succeed; the initial status file alone does not establish readiness. Stop only the process you started when finished. The API binds to `127.0.0.1` and rejects cross-origin writes; it is a local application, not a hosted multiuser service.

## Endpoint contract

| Request | Result |
| --- | --- |
| `GET /api/desktop` | Version, offline flag, runtime data directory |
| `GET /api/drawing/library` | Source recipes, styles, monomers, templates, linkages and registry errors |
| `GET /api/drawing/sources/{recipe_id}` | Hash-bound source recipe and qualified metadata |
| `POST /api/drawing/monomer` | Validated monomer graph, SMILES, MOL block and depiction |
| `POST /api/drawing/prepare` | Editable finite recipe for explicit monomer inputs |
| `POST /api/drawing/render` | Audit and download URLs for actual generated files |
| `POST /api/drawing/layout` | Explicitly requested cleanup revision and diagnostics; inspect the returned recipe |
| `GET /api/drawing/exports/{drawing_id}/{extension}` | A generated SVG, PNG, MOL, SDF, CDXML, JSON, XYZ, CIF or PDB file, when available |

Read valid IDs and supported status from `library` rather than maintaining a separate invented catalog. Requests reject unknown top-level fields. HTTP 422 carries the invalid-input reason.

### Source drawing

Render request:

```json
{"source_recipe_id":"tppa-pore-v1","style_id":"standard","size":1200}
```

`render` requires **exactly one** of `source_recipe_id` and `recipe`. Size is 500–2400 pixels. Optional fields are `style_seed` (0–2147483647) and `style_options` (object). Use the controls exposed by the library and the UI as the style-options reference. Retain the seed and resolved options for random styles.

Use `style_id: "source"` for source-informed colors, continuation marks, fills and supported abbreviations. Source coordinates remain fixed in all normal rendering. For example, `pyurea-source-layout-v1` retains its compact pyrene-block arrangement and `ucof-source-layout-v1` retains the source-informed angular arrangement with red urethane links. Read the recipe's coordinate qualifications; neither is claimed to be exact pixel registration. Source style is a vector reconstruction, not a copy of the publication image.

### Coordinate editing

Recipes may use component fragments and instances, or an `explicit_graph`. For the latter, atom `index` and `xy: [x,y]` use Cartesian y-up drawing coordinates. Edit `explicit_graph.atoms[i].xy` while retaining all atom/bond, hydrogen, charge, isotope, stereo and origin fields. For component recipes with `atom_coordinates`, positions are keyed by `instance` and `fragment_atom_index`. Preserve these coordinates rather than changing only the coarse instance positions.

The audit's `atom_origin` maps output atoms back to these keys; direct graphs also include `explicit_atom_index`. Return edited recipes through `/render` as designed variants. `POST /layout` accepts the usual recipe/source selection plus `operation: "clean"`; it is an explicit geometry change, not a prerequisite for source rendering. Automatic cleanup is supported for assembled component recipes; direct `explicit_graph` recipes retain their supplied coordinates and reject this operation. Review the returned diagnostics before accepting a revision, then render it as an edited recipe. The helper exposes this as `layout --request layout.json` and saves the response in `result.json`.

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

### Geometry review

`half_space_stack` places shaded half-size van der Waals spheres behind the flat chemical drawing. These back layers illustrate stacking from the drawing coordinates. For coordinate-based views, use a periodic model or imported CIF in the 3D viewer. The `element_compact` style uses jewel-tone element colors and compact hydrogen labels. Set `style_options.hydrogen_labels` to `show` or `compact` independently. Compact labels omit H on neutral, nonaromatic N, O, P and S where stereochemistry permits. The chemical export retains the same atoms and bonds.

New `/prepare` designs use `bond_angle_pose_v2`: proper rotations and translations of rigid monomer fragments, with alternative planar connection-side choices only where stereochemistry is unspecified. Ring shapes and internal atom distances remain fixed. Explicitly stereospecified fragments use the conservative spacing path and request review. This is a depiction operation; a small angle error does not establish a physically stable or synthesizable COF.

Read `audit.layout_refinement.quality`, `final_max_junction_angle_deviation_degrees`, and `max_join_length_fraction_error`. The `regular` result requires junction deviation at most 15 degrees, join-length error at most 12%, and no detected atom overlaps or proper nonincident bond crossings. `audit.layout_diagnostics` contains the individual angle measurements, bond ratios, overlap pairs and crossing pairs. These diagnostics do not cover every possible label collision or all forms of segment overlap. Report a `needs_review` result, inspect the image, and use manual adjustment or a different construction where appropriate.

When comparing a reconstruction with a supplied original, bind the original crop and graph identities first. Use corresponding atoms or monomer landmarks and a global proper rotation, uniform scale and translation; do not deform individual rings to improve pixel overlap. Inspect both the overlay and chemistry. Image ink-distance scores can help rank candidates but also measure labels and unrelated marks, so they are not proof of atomwise registration. Retain explicit selected coordinates and qualify source-informed layouts separately from exact source registration. The bundled TpPa and TPB-DMTP revisions follow this approach; unchanged source requests use those coordinates directly.

Each helper invocation requires a new output directory. `library.json` preserves catalog selection context. A render creates `result.json`, `drawing.svg`, `drawing.png`, `drawing.mol`, and any returned optional formats, including the JSON recipe/audit bundle. `artifacts.json` records file sizes and SHA-256 hashes. `server-output.log` and `runtime-data/cofdraw.log` are local diagnostics.

The response `audit` records the graph scope, drawing status, monomers, template, style, and chemical export diagnostics. Check this before describing an output. CDXML export may be unavailable; only advertise formats present in `downloads`. A successful API response is a software check, not independent verification of the chemistry.
## Single-precursor and mixed-linker designs

Read `assembly_examples` from the library. Each example contains its input precursor IDs, sequence and a DOI where available. POST its `monomer_ids`, `linker_sequence` and optional `pore_layout` to `/api/drawing/prepare/assembly`. Run the returned recipe through the standard render endpoint.

For BFBAEB-COF use `{"monomer_ids":["bfbaeb"],"pore_view":"triangle"}` or choose `"hexagon"`. For SIOC-COF-1 use `{"monomer_ids":["tape","ta","bpda_cho"],"pore_layout":"kgm","linker_sequence":"BCBCBC"}`. Use `BBB` or `CCC` for either triangular environment in that same arrangement. For an ordinary mixed-linker design, use the library's `mixed_pore_patterns` for a four- or six-sided pore and omit `pore_layout`. Each sequence specifies a local chemical arrangement. Assigning a bulk occupancy or an experimental pore-size distribution requires independent evidence.

To use a new precursor, POST its `smiles` and `label` to `/api/drawing/assembly/precursor`. Include the returned ID in `monomer_ids` and the original input in `custom_precursors`. Read `assembly_roles` to check its suitability. Use `assembly.precursors` for the feed molecules when a reaction divides one precursor into several reacted fragments.

To reopen a compact COFex record, POST `{"text":"<complete COFex>"}` to `/api/drawing/cofex/resolve`. Inspect the returned choices. Use the supplied source ID, preparation request, assembly request or crystal project for the chosen interpretation. Preserve the distinction between a generated design and a cited reference phase. The shared Framework canvas accepts pasted COFex and can export a reaction figure with the displayed precursors, a user-editable name, and a measured pore arrow. Project JSON stores the figure name and arrangement.

Choose a larger connected view from the render result's `repeat_options.extents`. POST the same `recipe` or `source_recipe_id` to `/api/drawing/repeat` with the selected option's `repeats` (or its numeric `value` when omitted), `extent_mode` (defaults to `pores`), optional `color_pores` and rendering options. For the compact three-pore cluster, use `{"repeats":3,"extent_mode":"triad"}` when offered. Use the returned `downloads` for the displayed graph and read `summary.pore_count`, `pore_counts` or `pore_sequences` for the represented faces. Extent 1 uses the original render result. A source without a reviewed expansion retains its witnessed drawing.

The binary `kgm_4_2` template accepts ETTA with a linear aromatic dialdehyde for an imine design. Its drawing extents include the small triangular opening and a hexagon sharing edges with six triangles. For a repeated drawing with `summary.measurement_face`, POST an empty object to `/api/drawing/repeat/{drawing_id}/measure`. Read `measurement.focal_face` for the measured family, `measurement.diameter_nm` for the calibrated in-plane opening, and `audit.drawing_style.layout` for its placement in the visible SVG. The structural choices follow the dual-pore and mixed-linker precedents. [Binary Kagome](https://doi.org/10.1021/ja5092936) · [Mixed-linker Kagome](https://doi.org/10.1021/jacs.6b01244) · [Self-condensed Kagome](https://doi.org/10.1021/jacs.9b03463)


GET `/chemical-space.json` for the offline explorer snapshot. Each point contains a linkage family `f`, an unordered precursor-ID list `ids`, indexed [material name, DOI] references `refs`, and 2D t-SNE coordinates `xy`. Match all precursor IDs and the linkage. Each family's eligible binary pairs encode `min(id1,id2) * stride + max(id1,id2)`. Eligibility includes chains. An empty evidence set means no match in the screened index. The independent CIF property table records model IDs, source DOI, CIF SHA-256, density (g/cm³), geometric void percentage, sampled cavity diameter (Å), and volume sample count.
For the family screening panel, GET `/chemical-screen.json` for the row schema, reacted repeat descriptors, graph identifiers, modeled CIFs and grid results. Join candidate codes to `/chemical-space.json` within the same linkage family. GET `/api/drawing/space/blocks/{id}.svg` for a precursor depiction. POST `/api/drawing/space/export/{family}` and download its returned ZIP URL for batch work. Read `methods.json` and the status columns before using numeric fields. Uncomputed electronic and adsorption targets are initialized to zero and marked pending.


### Spatial imine drawings

The library's templates with `dimension: 3` and `shape_pair` describe finite molecular circuits from the selected parent net. Inspect both precursor roles and the site geometry before preparing a recipe. For example, a prismatic HFPTP-F aldehyde and a linear Bpy-55 diamine support the `acs_6_2` design. Use the extents returned by the render response before requesting repeat exports. Some quadrangular and ditopic pairs have both sheet and spatial options. Keep the chosen parent net, experimental phase assignment, and availability of an atomistic periodic model distinct.

For a spatial recipe, `POST /api/drawing/science/patch` can return `kind: spatial_building_block_net`. Its labeled vertices are building blocks in relative coordinates. Physical pore, density and diffraction calculations require an atomistic cell. Drawing measurements on unfolded spatial circuits return `quantity_label: Projected circuit opening`.
