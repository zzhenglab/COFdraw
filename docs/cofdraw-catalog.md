# Drawing catalog and source identity

`scripts/export_cofdraw_catalog.py` exports the portable catalog used by the
desktop drawing application. It reads the local `cofex.sqlite` and the four
hash-bound source recipes, without modifying either. The catalog contains six
precursor presets, four source-associated starting presets and twenty researched
linkage families. Only imine and beta-ketoenamine have bounded automatic monomer
transforms in the current drawing API.

Each source entry supplies the **full COFex notation**, including precursor
SMILES and the `COFex.unk.cat0;name` suffix, alongside its compact key. The full
notation is copied from the installed reference only after matching the DOI,
literature name and exact canonical precursor set. Duplicate matching records
are retained as a list. The notation describes a precursor-set reference;
it does not serialize a complete periodic framework. `unk` remains unknown.
`cat0` is explicitly qualified as the reference convention, without an
independent claim that this particular source proves noninterpenetration.

Names are paper-specific. The TFB-TAB drawing uses TAB for the structure commonly
called TAPB. The TPB-DMTP-COF paper labels its precursor drawings TAPB and DMTA.
The catalog retains these source labels and their matched graphs rather than
substituting an abbreviation globally. SCOF and TpPa-SO3H are retained as aliases.

Selecting a precursor preset and calling `/prepare` creates a **new design**.
Loading its `source_recipe_id` preserves the audited source incidence. A change
to monomers, linkage, ports or connectivity must clear the source identity.
Neither an alternating pore template nor its visual regularity determines a
unique topology, stacking, interpenetration or synthetic feasibility.

`source_styles` records observed colors, glyph conventions, abbreviations and
continuation marks separately from runtime style controls. Font substitutions,
sampled raster colors and generic layouts are approximate. In particular,
SO3H/HO3S and OCH3/H3CO orientations must retain their graph-bound attachment
atoms. Hidden atoms inside an abbreviation are not observed pixel coordinates.
Dashed, wavy and crossbar boundaries are continuation marks, not stereochemistry.

The export contains DOI/page/hash evidence but no private file paths, source PDF
text or publication images. The exporter refuses to overwrite an existing
catalog; produce a new version when the source registry or reference changes.
