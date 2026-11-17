"""Bounded local HTTP interface for explicit COF pore drawing recipes."""
from __future__ import annotations

import hashlib
import json
import math
import re
import threading
import uuid
from copy import deepcopy
from pathlib import Path
from typing import Literal

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel, ConfigDict, Field, model_validator


class PrepareDrawing(BaseModel):
    model_config = ConfigDict(extra="forbid")
    linkage: str = Field(default="imine", max_length=80)
    monomer_a: str = Field(min_length=1, max_length=4000)
    monomer_b: str = Field(min_length=1, max_length=4000)
    label_a: str = Field(default="A", max_length=100)
    label_b: str = Field(default="B", max_length=100)
    instances: Literal[4, 6, 8, 12] = 6
    template_id: str | None = Field(default=None, max_length=80)


class RenderDrawing(BaseModel):
    model_config = ConfigDict(extra="forbid")
    recipe: dict | None = None
    source_recipe_id: str | None = Field(default=None, max_length=150)
    size: int = Field(default=1600, ge=500, le=2400)
    style_id: str = Field(default="standard", max_length=80)
    style_seed: int | None = Field(default=None, ge=0, le=2147483647)
    style_options: dict = Field(default_factory=dict)

    @model_validator(mode="after")
    def exclusive_source(self):
        if (self.recipe is None) == (self.source_recipe_id is None):
            raise ValueError("Supply either an editable recipe or a registered source recipe ID")
        return self


class MonomerDrawing(BaseModel):
    model_config = ConfigDict(extra="forbid")
    smiles: str | None = Field(default=None, max_length=10000)
    molblock: str | None = Field(default=None, max_length=500000)
    graph: dict | None = None
    clean_layout: bool = False


def validate_recipe(recipe: dict) -> None:
    """Reject costly or degenerate graphs before the renderer allocates them."""
    from rdkit import Chem
    from .framework_drawing import VERSION

    if len(json.dumps(recipe, allow_nan=False)) > 200_000:
        raise ValueError("Recipe exceeds 200 kB")
    if recipe.get("schema_version") != VERSION:
        raise ValueError("Unsupported drawing recipe version")
    components = recipe.get("components")
    instances = recipe.get("instances")
    connections = recipe.get("connections")
    if not isinstance(components, dict) or not 1 <= len(components) <= 16:
        raise ValueError("Provide 1–16 components")
    if not isinstance(instances, list) or not 3 <= len(instances) <= 24:
        raise ValueError("Provide 3–24 ordered pore instances")
    if not isinstance(connections, list) or not 3 <= len(connections) <= 64:
        raise ValueError("Provide 3–64 explicit port connections")
    counts, ports = {}, {}
    for name, component in components.items():
        if not isinstance(name, str) or not name or len(name) > 100 or not isinstance(component, dict):
            raise ValueError("Invalid component name or definition")
        smiles = component.get("fragment_smiles")
        if not isinstance(smiles, str) or len(smiles) > 10000:
            raise ValueError("Each component needs a finite fragment SMILES")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or len(Chem.GetMolFrags(mol)) != 1:
            raise ValueError("Component fragment must be one valid connected graph")
        counts[name] = mol.GetNumAtoms()
        ports[name] = {a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetAtomicNum() == 0}
    total = 0
    positions, ids = [], set()
    for item in instances:
        if not isinstance(item, dict):
            raise ValueError("Each pore instance must be an object")
        name, component = item.get("id"), item.get("component")
        if not isinstance(name, str) or not name or len(name) > 100 or name in ids or component not in components:
            raise ValueError("Instances need unique IDs and existing components")
        ids.add(name)
        total += counts[component]
        incoming, outgoing = item.get("incoming_port"), item.get("outgoing_port")
        if incoming == outgoing or incoming not in ports[component] or outgoing not in ports[component]:
            raise ValueError("Choose two different existing ports for every pore instance")
        xy = item.get("layout_position")
        if not isinstance(xy, list) or len(xy) != 2 or any(not isinstance(x, (int, float)) or not math.isfinite(x) or abs(x) > 200 for x in xy):
            raise ValueError("Layout coordinates must be finite numbers between −200 and 200")
        positions.append(xy)
    if total > 1000:
        raise ValueError("Drawing limit is 1,000 atoms including continuation ports")
    for index, xy in enumerate(positions):
        if math.dist(xy, positions[index - 1]) < 0.05:
            raise ValueError("Adjacent pore instances need different positions")
    for edge in connections:
        if not isinstance(edge, dict):
            raise ValueError("Each connection must be an object")
        if edge.get("bond_order") not in (1, 2, 3):
            raise ValueError("Connection bond order must be 1, 2 or 3")
        for end in (edge.get("a"), edge.get("b")):
            if not isinstance(end, dict) or end.get("instance") not in ids or not isinstance(end.get("port"), int):
                raise ValueError("Every connection must identify an existing instance and integer port")


def create_drawing_router(data_root: Path) -> APIRouter:
    from .framework_drawing import LINKAGE_LIBRARY, VERSION, linkage_fragment, render_recipe, build_patch
    from .cofdraw_core import STYLES, MONOMERS, PRESETS, read_catalog, normalize_monomer, render_styled, export_chemical_formats
    from .cofdraw_chemistry import TEMPLATES, library_linkages, decorated_monomers, prepare_template, refresh_design_identity, design_signature
    from .cofdraw_styles import PALETTES, FONTS

    router = APIRouter(prefix="/api/drawing", tags=["COF drawing"])
    root = data_root.resolve()
    registry = root / "framework-drawing" / "registry.json"
    outputs = root / "framework-drawing" / "exports"
    lock = threading.Lock()

    def source_records():
        if not registry.is_file():
            return []
        return json.loads(registry.read_text(encoding="utf-8"))["recipes"]

    def source_recipe(recipe_id):
        matches = [r for r in source_records() if r.get("recipe_id") == recipe_id]
        if len(matches) != 1:
            raise ValueError("Source recipe is missing or ambiguous")
        record = matches[0]
        path = (registry.parent / record["path"]).resolve()
        if not path.is_relative_to(root) or not path.is_file() or path.suffix != ".json":
            raise ValueError("Source recipe is outside this data workspace")
        if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError("Source recipe hash mismatch")
        value = json.loads(path.read_text(encoding="utf-8"))
        if value.get("recipe_id") != recipe_id or value.get("cofex_id") != record.get("cofex_id"):
            raise ValueError("Registered source identity mismatch")
        return value, record["sha256"]

    def source_metadata(recipe, digest):
        from .framework_drawing import canonical
        value = read_catalog(root).get("source_metadata", {}).get(recipe["recipe_id"], {})
        if not value:
            return {}
        if value.get("cofex_id") != recipe.get("cofex_id") or value.get("source_recipe_sha256") != digest:
            raise ValueError("Catalog metadata does not match this source recipe revision")
        notation = value.get("cofex_notation", "")
        if " COFex." not in notation:
            raise ValueError("Catalog source is missing its full SMILES-based COFex notation")
        expected = {canonical(c["monomer_smiles"]) for c in recipe["components"].values()}
        named = {canonical(c["smiles"]) for c in value.get("monomers", [])}
        encoded = {canonical(s) for s in notation.split(" COFex.", 1)[0].split(".")}
        if expected != named or expected != encoded:
            raise ValueError("Catalog names/SMILES/COFex precursor sets do not agree")
        return value

    @router.get("/library")
    def library():
        sources, errors = [], []
        try:
            catalog = read_catalog(root)
        except (ValueError, OSError) as exc:
            catalog = {}; errors.append(str(exc))
        for row in source_records():
            try:
                value, digest = source_recipe(row["recipe_id"])
                metadata = source_metadata(value, digest)
                sources.append({"recipe_id": value["recipe_id"], "label": row.get("label") or value.get("label") or value["recipe_id"],
                                "cofex_id": value.get("cofex_id"), "sha256": digest, **metadata})
            except (ValueError, KeyError, OSError) as exc:
                errors.append(str(exc))
        return {"schema_version": VERSION, "linkages": library_linkages(catalog), "source_recipes": sources,
                "templates": TEMPLATES, "topology_catalog": catalog.get("topology_catalog", []),
                "style_controls": {"palettes":list(PALETTES),"fonts":FONTS,"edge_styles":["wavy_cross","dashed","crossbar","bond_stub"],
                                   "abbreviation_modes":["explicit","condensed","OMe"],"layer_count":[2,5],"layer_opacity":[.05,.7],
                                   "fill_opacity":[.05,.6],"font_scale":[.6,1.6],"bond_scale":[.5,2.5]},
                "styles": STYLES, "monomer_presets": decorated_monomers(catalog.get("monomer_presets", MONOMERS)),
                "cof_presets": catalog.get("cof_presets", PRESETS), "catalog_version": catalog.get("catalog_version", "builtin-1"),
                "style_evidence": catalog.get("style_evidence", catalog.get("source_styles", [])),
                "source_registry_errors": errors, "max_atoms": 1000,
                "presets": [{"instances": n, "label": f"{n // 2} A + {n // 2} B alternating pore"} for n in (6, 12)],
                "scope": "Finite chemical pore patches with explicit attachment ports; preset connectivity is a design choice, not a topology prediction."}

    @router.get("/sources/{recipe_id}")
    def source(recipe_id: str):
        try:
            value, digest = source_recipe(recipe_id)
            metadata = source_metadata(value, digest)
            return {"recipe": value, "source_sha256": digest, "status": "source_bound_provisional", "metadata": metadata}
        except (ValueError, KeyError, OSError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post("/monomer")
    def monomer(request: MonomerDrawing):
        try:
            return normalize_monomer(**request.model_dump())
        except (ValueError, KeyError, TypeError, RuntimeError, IndexError, OverflowError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post("/prepare")
    def prepare(request: PrepareDrawing):
        try:
            prepared = prepare_template(**request.model_dump())
            recipe = {"schema_version": VERSION, "recipe_id": "designed-pore", **prepared,
                      "topology": None, "topology_basis": "explicit_user_chosen_template_not_experimental",
                      "layout_refinement": "rigid_block_translation_v1",
                      "cofex_id": None, "cofex_notation": None, "cofex_string": None, "designed_variant": True, "ground_truth": False, "review_status": "user_design",
                      "design_note": "Reactive ports were chosen by their 2D angular fit to the selected template. Inspect incidence; a template does not establish experimental topology or synthesis feasibility."}
            recipe["design_template_signature"] = design_signature(recipe)
            validate_recipe(recipe)
            return {"recipe": recipe, "status": "designed_variant", "selected_template": recipe.get("selected_template"),
                    "compatibility": recipe.get("compatibility"), "limits": LINKAGE_LIBRARY[request.linkage]["limits"]}
        except (ValueError, KeyError, TypeError, RuntimeError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post("/render")
    def render(request: RenderDrawing):
        try:
            source_hash = None
            if request.source_recipe_id:
                recipe, source_hash = source_recipe(request.source_recipe_id)
                status = "source_bound_provisional"
                metadata = source_metadata(recipe, source_hash)
            else:
                recipe = deepcopy(request.recipe)
                # User edits must never inherit the identity or truth status of a source recipe.
                recipe["derived_from_recipe_id"] = recipe.get("recipe_id")
                recipe.update(cofex_id=None, cofex_notation=None, cofex_string=None, parent_cofex_id=None, source=None, reference=None, topology=None, topology_evidence=None,
                              designed_variant=True, ground_truth=False, review_status="user_design")
                validate_recipe(recipe)
                refresh_design_identity(recipe)
                status = "designed_variant"
                metadata = {}
            validate_recipe(recipe)
            drawing_id = uuid.uuid4().hex
            original_id = recipe.get("recipe_id")
            recipe["recipe_id"] = drawing_id  # filenames are always server-generated
            directory = outputs / drawing_id
            with lock:
                prepared_patch = build_patch(recipe)
                audit = render_styled(recipe, directory, request.size, request.style_id, request.style_seed, request.style_options,
                                      prepared_patch=prepared_patch)
                # Export support is explicit; the underlying MOL remains available.
                extra_formats = export_chemical_formats(directory, drawing_id, source_mol=prepared_patch[0])
            export_report = directory / f"{drawing_id}.chemical-export.json"
            if export_report.is_file():
                audit["chemical_exports"] = json.loads(export_report.read_text(encoding="utf-8"))
            audit["graph_scope"] = "finite_source_pore_patch_with_continuation_ports" if source_hash else "finite_designed_pore_patch_with_continuation_ports"
            audit.update(drawing_status=status, registered_recipe_sha256=source_hash, original_recipe_id=original_id)
            audit.update(design_cofex_notation=recipe.get("design_cofex_notation"), design_topology=recipe.get("design_topology"),
                         design_identity_qualification=recipe.get("design_identity_qualification"), selected_template=recipe.get("selected_template"))
            audit.update(cofex_string=metadata.get("cofex_notation"), cofex_notation=metadata.get("cofex_notation"),
                         literature_metadata=metadata, monomers=[{"component": name, "label": c.get("label", name), "smiles": c.get("monomer_smiles")} for name, c in recipe["components"].items()],
                         chemical_export_scope="Editable finite patch with continuation atoms. SVG styles are separate presentation; no complete periodic unit cell is implied.")
            urls = {ext: f"/api/drawing/exports/{drawing_id}/{ext}" for ext in ("svg", "png", "mol", "json", *extra_formats)}
            bundle = {"recipe": recipe, "audit": audit, "downloads": urls}
            (directory / f"{drawing_id}.json").write_text(json.dumps(bundle, indent=2) + "\n", encoding="utf-8")
            return {"drawing_id": drawing_id, "status": status, "audit": audit, "downloads": urls}
        except (ValueError, KeyError, TypeError, RuntimeError, OSError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.get("/exports/{drawing_id}/{extension}")
    def artifact(drawing_id: str, extension: Literal["svg", "png", "mol", "json", "sdf", "cdxml"]):
        if not re.fullmatch(r"[0-9a-f]{32}", drawing_id):
            raise HTTPException(404, "Drawing not found")
        path = (outputs / drawing_id / f"{drawing_id}.{extension}").resolve()
        if not path.is_relative_to(outputs.resolve()) or not path.is_file():
            raise HTTPException(404, "Drawing not found")
        types = {"svg": "image/svg+xml", "png": "image/png", "mol": "chemical/x-mdl-molfile", "json": "application/json", "sdf": "chemical/x-mdl-sdfile", "cdxml": "chemical/x-cdxml"}
        return FileResponse(path, media_type=types[extension], filename=f"cof-pore-{drawing_id[:8]}.{extension}")

    return router
