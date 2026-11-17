# Text-to-drawing planning

The Describe planner turns supported wording into a **reviewable finite drawing request**. It does not infer an experimental COF identity from arbitrary prose, invent a monomer SMILES, choose a proven periodic topology, predict a synthesis, or produce a crystal structure. The final drawing still uses COFdraw's existing preparation/rendering code.

The local planner is deterministic and works offline. An optional language-model proposal can select catalog IDs; those IDs and the actual chemical preparation are still checked locally. A model proposal does not provide chemical evidence.

## Backend interface

```python
from cof_vision.cofdraw_describe import describe_catalog, describe_cof, validate_selection

# library is the exact current GET /api/drawing/library result.
choices = describe_catalog(library)
plan = describe_cof("Draw TFB and TAPB in pastel colors", library)
```

`describe_catalog` exposes sources, monomers, linkages, templates, styles, palettes, and demonstrations. Option rows have `id`, `label`, and `aliases`. Source rows retain qualifiers, hash, registration status, and precursor IDs. Monomer rows expose recognized compatible roles; they do not expose generated SMILES. Unsupported linkage records remain visible with `supported: false`.

`describe_cof(text, library, selection=None, *, allow_unmatched_words=False)` accepts 1–4000 characters. It returns:

- `status`: `ready`, `needs_input`, or `unsupported`.
- `summary`, `scope`, `warnings`, `unresolved`, and explicit review `assumptions`.
- `selection`: the seven source/monomer/linkage/template/style/palette IDs, with unused fields null.
- `matches`: matched catalog IDs and their text spans in the normalized input.
- `prepare_request`: an actual `PrepareDrawing` payload only for a ready new design.
- `source_recipe_id`: a registered recipe ID only for a ready source drawing.
- `style_id` and `style_options`, including the selected palette when supplied.
- `alternatives`: working demonstration prompts when the request needs clarification.

For a ready source plan, `prepare_request` is null. For a ready design, `source_recipe_id` is null. Unresolved plans expose neither actionable request. Callers must check `status` before creating a drawing, and preserve returned source qualifiers and style warnings.

## Names, roles, and conservative defaults

Names come from the active catalog: IDs, labels, full names, abbreviations, and the existing source precursor literature aliases. Unicode dashes, beta, subscripts, and whitespace are normalized. Longest-name matching prevents `TpPa` from taking precedence over `TpPa-SO3H`, or `COF-BP` over `COF-BP-Ni2`. Source aliases such as Pyurea and FCOF1 resolve to their specific registered finite drawing, with its recorded qualifications.

Short monomer abbreviations require the catalog's deliberate capitalization or an explicit monomer/plus-sign context. For example, `Pa` must not be inferred just from an incidental lowercase word. Colliding aliases such as `BD` retain all candidates until the user chooses one; the catalog contains distinct monomers with that abbreviation.

Two matched monomers are ordered into the reaction roles required by the implemented transform. A unique conventional supported linkage can be proposed, with an assumption recorded. Reduced-imine products are never the automatic default. When several transforms match, an existing unambiguous catalog starting-pair preference can be proposed; otherwise the linkage remains unresolved. A starting-pair preference creates a new design, not a source-identical material.

A unique compatible finite template can also be proposed for review. This is a drawing-layout choice, **not periodic-topology inference from monomer degree**. Explicitly requested incompatible templates fail. Tetrahedral nodes do not silently become planar nets. Every ready design calls the real `prepare_template` using the selected catalog monomer SMILES; role metadata alone cannot create a ready result.

Standard chemical style is the default. A pastel/soft palette request without an explicit style selects pastel regions behind dark chemical bonds for legibility. Other palette requests select colored building blocks so the color is visible. Both choices record an assumption. Explicit colored bonds retain the colored-bond style, while an explicit black-bond style stays black and warns that a palette may have no visible effect. Layers, tilted sheets, spheres, and cartoons remain qualified illustrations.

## Review selections and free wording

Selections allow only `source_recipe_id`, `monomer_a_id`, `monomer_b_id`, `linkage`, `template_id`, `style_id`, and `palette`, with IDs from the current catalog. Null values mean no selected value. Source and changed-design selections cannot be mixed. Unknown IDs, arbitrary SMILES, recipe objects, code, and additional fields cannot bypass the checks.

By default, unrecognized words or names keep a request unresolved. This prevents an unknown second monomer or extra constraint from being silently dropped. Vague requests such as “draw a COF” return examples and require choices. Negation and replacement language require a positive explicit request rather than a guessed interpretation.

Explicit di/tri/tetra/hexa-topic or numbered connected amine/aldehyde/hydrazide requirements are checked against the selected monomers' recognized functional groups. Point-group labels such as C3 are not interpreted as connectivity. Numeric instance/unit counts must match the finite template; a registered source does not silently change its incidence to satisfy a new count. A physical target such as a 2 nm pore cannot be promised by a drawing layout. Such targets, 3-D atomic structures, and XYZ/PDB coordinate generation are routed back to the explicit-input scientific tools. Drawing exports such as SVG and PNG remain supported.

For a language-model proposal or a user-reviewed nonempty selection, a caller may explicitly set `allow_unmatched_words=True`. Only residual-word handling changes: unmatched wording appears in a prominent warning requiring verification. Recognized source, monomer, linkage, and geometry conflicts, negation, unsupported 3D/property requests, ID allowlists, and actual preparation validation remain hard checks. The flag alone, without selected IDs, does not invent a choice.

This reviewed mode **does not guarantee semantic fidelity** for wording the local parser cannot understand. It can validate the selected chemistry and drawing request; the user must verify that the chosen identities and parameters reflect the intended free description. An unknown name is not independently identified merely because a proposed catalog selection is chemically compatible.

Explicitly reviewed style and palette changes may override the appearance wording; the override is recorded in `assumptions`. A selected chemical identity or template cannot override a contradictory recognized text constraint. Edit the description to change that instruction.

## Examples and verification

The built-in demonstrations include pastel TFB/TAPB imine design, qualified source TpPa, Pyurea highlights, F-COF-1 highlights and source colors, an acylhydrazone design, a square finite design, an unsupported boronate construction, and a vague request that needs a choice.

Focused tests exercise real catalog aliases and source hashes, longest-name precedence, collisions, short-name boundaries, reversed reaction roles, style effects, explicit assumptions, source/design separation, unknown names, negation, chemistry/template conflicts, unsupported requests, and attempted selection bypasses. Ready demonstration designs pass the same request schema and chemical/template validator as the API and build actual finite molecular graphs. Tests also verify that reviewed free wording keeps its warning and cannot relax recognized chemical constraints.
