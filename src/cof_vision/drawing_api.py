"""Bounded local HTTP interface for explicit COF pore drawing recipes."""
from __future__ import annotations

import hashlib
import json
import math
import re
import threading
import uuid
from datetime import datetime, timezone
from copy import deepcopy
from pathlib import Path
from typing import Annotated, Literal

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse, Response
from pydantic import BaseModel, ConfigDict, Field, SecretStr, model_validator


class PrepareDrawing(BaseModel):
    model_config = ConfigDict(extra="forbid")
    linkage: str = Field(default="imine", max_length=80)
    monomer_a: str = Field(min_length=1, max_length=4000)
    monomer_b: str = Field(min_length=1, max_length=4000)
    label_a: str = Field(default="A", max_length=100)
    label_b: str = Field(default="B", max_length=100)
    instances: Literal[4, 6, 8, 12] = 6
    template_id: str | None = Field(default=None, max_length=80)


class CofexInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    text: str = Field(min_length=1, max_length=20000)


class AssemblyPrecursor(BaseModel):
    model_config = ConfigDict(extra="forbid")
    smiles: str = Field(min_length=1, max_length=4000)
    label: str = Field(default="Custom precursor", min_length=1, max_length=100)


class AssemblyDrawing(BaseModel):
    model_config = ConfigDict(extra="forbid")
    monomer_ids: list[Annotated[str, Field(min_length=1,max_length=150)]] = Field(min_length=1,max_length=3)
    linker_sequence: str = Field(default="",max_length=6,pattern="^[BC]*$")
    pore_view: Literal["triangle", "hexagon"] | None = None
    pore_layout: Literal["kgm"] | None = None
    custom_precursors: list[AssemblyPrecursor] = Field(default_factory=list, max_length=3)


class RenderDrawing(BaseModel):
    model_config = ConfigDict(extra="forbid")
    recipe: dict | None = None
    source_recipe_id: str | None = Field(default=None, max_length=150)
    size: int = Field(default=1600, ge=500, le=2400)
    style_id: str = Field(default="standard", max_length=80)
    style_seed: int | None = Field(default=None, ge=0, le=2147483647)
    recent_styles: list[str] = Field(default_factory=list, max_length=4)
    style_options: dict = Field(default_factory=dict)

    @model_validator(mode="after")
    def exclusive_source(self):
        if (self.recipe is None) == (self.source_recipe_id is None):
            raise ValueError("Supply either an editable recipe or a registered source recipe ID")
        return self


class RepeatDrawing(RenderDrawing):
    repeats: int = Field(ge=1, le=4)
    extent_mode: Literal['pores', 'cells', 'triad', 'star', 'triangle'] = 'pores'
    color_pores: bool = False


class MonomerDrawing(BaseModel):
    model_config = ConfigDict(extra="forbid")
    smiles: str | None = Field(default=None, max_length=10000)
    molblock: str | None = Field(default=None, max_length=500000)
    graph: dict | None = None
    clean_layout: bool = False


class TemplateSheet(BaseModel):
    cdxml: str = Field(min_length=1, max_length=500000)


class BinaryTemplateSheet(BaseModel):
    content_base64: str = Field(min_length=1, max_length=666668)


class LayoutDrawing(RenderDrawing):
    operation: Literal['clean'] = 'clean'


class ModelDrawing(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False)
    template_id: str | None = Field(default=None, max_length=80)
    cif_text: str | None = Field(default=None, max_length=2_000_000)
    repeats: list[Annotated[int, Field(ge=1, le=4)]] = Field(default=[1,1,1], min_length=3, max_length=3)
    recipe_match: RenderDrawing | None = None
    linker_edit: Literal['aryl_f'] | None = None


class CrystalDrawingView(BaseModel):
    model_config = ConfigDict(extra="forbid")
    arm_spread: int = Field(0, ge=-24, le=24)
    tilt: int = Field(0, ge=-45, le=45)


class DisplayCells(BaseModel):
    model_config = ConfigDict(extra='forbid')
    repeats: list[Annotated[int, Field(ge=1, le=4, strict=True)]] = Field(min_length=3, max_length=3)


class Patch3D(RenderDrawing):
    layers: int = Field(default=3, ge=1, le=6)
    layer_spacing: float = Field(default=3.4, ge=2.5, le=12, allow_inf_nan=False)
    stacking: Literal['AA','AB','ABC','slipped','twisted'] = 'AA'
    slip: list[Annotated[float, Field(ge=-20,le=20,allow_inf_nan=False)]] | None = Field(default=None,min_length=2,max_length=2)
    twist_degrees: float = Field(default=15, ge=-180, le=180, allow_inf_nan=False)


class PeriodicFromDrawing(RenderDrawing):
    layer_spacing: float | None = Field(default=None, ge=2.5, le=12, allow_inf_nan=False)
    repeats: list[Annotated[int, Field(ge=1, le=4)]] = Field(default=[1,1,1], min_length=3, max_length=3)


class DiamondDesign(BaseModel):
    model_config = ConfigDict(extra="forbid")
    monomer_a: str = Field(min_length=1, max_length=4000)
    monomer_b: str = Field(min_length=1, max_length=4000)
    label_a: str = Field(default="Monomer A", max_length=180)
    label_b: str = Field(default="Monomer B", max_length=180)
    linkage: Literal["imine","secondary_amine","acylhydrazone"] = "imine"
    topology: Literal["dia","lon","4h","6h"] = "dia"
    network_count: int = Field(default=1,ge=1,le=7,strict=True)
    repeats: list[Annotated[int, Field(ge=1,le=4)]] = Field(default=[1,1,1], min_length=3, max_length=3)


class NetworkHypothesis(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False)
    structure: dict
    offsets: list[list[Annotated[float, Field(ge=-1, le=1, allow_inf_nan=False)]]] = Field(min_length=1, max_length=4)
    repeats: list[Annotated[int, Field(ge=1, le=4)]] = Field(default=[1,1,1], min_length=3, max_length=3)


class MonomerAnalysis(BaseModel):
    model_config = ConfigDict(extra='forbid')
    smiles: str = Field(min_length=1, max_length=4000)
    conformers: int = Field(default=4, ge=1, le=8)


class SolubilityScreen(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)
    smiles: str | None = Field(default=None, min_length=1, max_length=10000)
    monomer_id: str | None = Field(default=None, min_length=1, max_length=80)
    solvent_ids: list[Annotated[str, Field(min_length=1, max_length=24)]] | None = Field(default=None, min_length=1, max_length=64)
    solvent_ml: float = Field(default=1.0, ge=0.01, le=1000)

    @model_validator(mode="after")
    def exclusive_monomer(self):
        if (self.smiles is None) == (self.monomer_id is None):
            raise ValueError("Supply a monomer ID or one SMILES")
        return self


class PoreAnalysis(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False)
    structure: dict
    probe_radius: float = Field(default=1.82, ge=0, le=3)
    samples: int = Field(default=12000, ge=1000, le=30000)


class GeometryAnalysis(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False, strict=True)
    structure: dict
    probe_radius: float = Field(default=1.82, ge=0, le=5)
    grid_spacing: float = Field(default=.6, ge=.1, le=3)
    surface_samples: int = Field(default=64, ge=16, le=256)
    refinement: bool = False
    include_map: bool = False

    @model_validator(mode='before')
    @classmethod
    def finite_controls(cls, value):
        # JSON parsers may accept nonstandard Infinity/NaN. Do not let them
        # enter FastAPI's validation-error JSON, which cannot serialize them.
        if isinstance(value, dict):
            for key in ('probe_radius', 'grid_spacing', 'surface_samples', 'refinement'):
                if isinstance(value.get(key), float) and not math.isfinite(value[key]):
                    raise HTTPException(422, f'{key} must be finite')
        return value


class BETAnalysis(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False)
    relative_pressures: list[float] = Field(min_length=3, max_length=2000)
    uptake_mmol_g: list[float] = Field(min_length=3, max_length=2000)
    fit_min: float = Field(default=.05, gt=0, lt=1)
    fit_max: float = Field(default=.30, gt=0, lt=1)
    cross_section_nm2: float = Field(default=.162, gt=0, le=1)


class StructureInspection(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False)
    structure: dict
    bond_scale: float = Field(default=1.15, ge=.9, le=1.35)
    tolerance_angstrom: float = Field(default=.15, ge=0, le=.3)


class PowderDrawing(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False)
    structure: dict
    wavelength: float = Field(default=1.5406, ge=.3, le=3)
    two_theta_min: float = Field(default=2, ge=.5, le=100)
    two_theta_max: float = Field(default=50, ge=1, le=150)
    step: float = Field(default=.02, ge=.01, le=.2)
    fwhm: float = Field(default=.15, ge=.02, le=2)
    b_iso: float = Field(default=0, ge=0, le=20)


class PeakComparison(PowderDrawing):
    observed_two_theta: list[Annotated[float, Field(gt=0, lt=180)]] = Field(min_length=1, max_length=128)
    tolerance_degrees: float = Field(default=.15, gt=0, le=1)
    min_relative_intensity: float = Field(default=1, ge=0, le=100)


class ReagentBalance(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False)
    monomer_a: str = Field(min_length=1, max_length=4000)
    monomer_b: str = Field(min_length=1, max_length=4000)
    linkage: str = Field(min_length=1, max_length=100)
    mass_a_mg: float = Field(ge=.000001, le=1_000_000)
    purity_a_percent: float = Field(default=100, ge=.000001, le=100)
    purity_b_percent: float = Field(default=100, ge=.000001, le=100)
    b_site_equivalents: float = Field(default=1, ge=.000001, le=100)


class ScreenEntry(BaseModel):
    model_config = ConfigDict(extra='forbid')
    name: str = Field(default='Structure', max_length=100)
    structure: dict


class FiniteInspection(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False)
    model: dict


class FinitePowder(FiniteInspection):
    wavelength: float = Field(default=1.5406, ge=.3, le=3)
    two_theta_min: float = Field(default=2, ge=.5, le=100)
    two_theta_max: float = Field(default=50, ge=1, le=150)
    step: float = Field(default=.1, ge=.01, le=.2)


class PeriodicStacking(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False)
    structure: dict
    stacking: Literal['AA', 'AB', 'ABC', 'slipped'] = 'AA'
    spacing: float = Field(default=3.4, ge=2.5, le=12)
    slip: list[Annotated[float, Field(ge=-20, le=20, allow_inf_nan=False)]] = Field(default=[1.5,0], min_length=2, max_length=2)
    corrugated_layer: bool = False


class ScreenBatch(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False)
    entries: list[ScreenEntry] = Field(min_length=1, max_length=8)
    include_porosity: bool = False
    probe_radius: float = Field(default=1.82, ge=0, le=3)
    samples: int = Field(default=4000, ge=1000, le=12000)


class KnowledgeQuery(BaseModel):
    model_config = ConfigDict(extra='forbid')
    query: str = Field(default='', max_length=500)
    record_types: list[str] | None = Field(default=None, max_length=10)
    limit: int = Field(default=20, ge=1, le=50)
    offset: int = Field(default=0, ge=0, le=10000)


class DescriptionSelection(BaseModel):
    model_config = ConfigDict(extra='forbid')
    source_recipe_id: str | None = Field(default=None, max_length=150)
    monomer_a_id: str | None = Field(default=None, max_length=150)
    monomer_b_id: str | None = Field(default=None, max_length=150)
    linkage: str | None = Field(default=None, max_length=80)
    template_id: str | None = Field(default=None, max_length=80)
    style_id: str | None = Field(default=None, max_length=80)
    palette: str | None = Field(default=None, max_length=80)


class DescriptionLLM(BaseModel):
    model_config = ConfigDict(extra='forbid')
    provider: Literal['compatible', 'openai', 'anthropic'] = 'compatible'
    base_url: str | None = Field(default=None, max_length=1000)
    model: str = Field(min_length=1, max_length=150)
    api_key: SecretStr = Field(default=SecretStr(''), exclude=True, repr=False)


class DescribeDrawing(BaseModel):
    model_config = ConfigDict(extra='forbid')
    text: str = Field(min_length=1, max_length=4000)
    mode: Literal['rules', 'llm'] = 'rules'
    selection: DescriptionSelection | None = None
    llm: DescriptionLLM | None = None


class ConversationContext(BaseModel):
    model_config = ConfigDict(extra='forbid')
    selection: DescriptionSelection = Field(default_factory=DescriptionSelection)
    pending: list[Annotated[str, Field(min_length=1, max_length=700)]] = Field(default_factory=list, max_length=12)
    requirements: list[Annotated[str, Field(min_length=1, max_length=80)]] = Field(default_factory=list, max_length=8)
    origin: Literal['plan', 'current_drawing'] = 'plan'
    token: str | None = Field(default=None, pattern=r'^[A-Za-z0-9_-]{1,80}$')


class DescribeConversation(DescribeDrawing):
    text: str = Field(min_length=1, max_length=3000)
    context: ConversationContext | None = None


def validate_recipe(recipe: dict) -> None:
    """Reject costly or degenerate graphs before the renderer allocates them."""
    from rdkit import Chem
    from .framework_drawing import VERSION

    if 'explicit_graph' in recipe:
        from .cofdraw_explicit import validate_explicit_recipe
        validate_explicit_recipe(recipe)
        return

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


def create_drawing_router(data_root: Path, *, output_root: Path | None = None) -> APIRouter:
    from .framework_drawing import LINKAGE_LIBRARY, VERSION, linkage_fragment, render_recipe, build_patch
    from .cofdraw_core import STYLES, PRESETS, catalog_monomers, read_catalog, normalize_monomer, render_styled, export_chemical_formats
    from .cofdraw_chemistry import TEMPLATES, library_linkages, decorated_monomers, prepare_template, refresh_design_identity, design_signature
    from .cofdraw_styles import PALETTES, FONTS, ELEMENT_SCHEMES, style_configuration
    from rdkit import Chem

    router = APIRouter(prefix="/api/drawing", tags=["COF drawing"])
    root = data_root.resolve()
    registry = root / "framework-drawing" / "registry.json"
    outputs = (output_root or root) / "framework-drawing" / "exports"
    lock = threading.Lock()
    from .cofdraw_adsorption import create_adsorption_router
    router.include_router(create_adsorption_router(outputs))

    @router.get('/space/screen-results')
    def local_screen_results():
        path = outputs/'research-screen.json.gz'
        if not path.is_file():
            return Response(status_code=204, headers={'Cache-Control': 'no-store'})
        return FileResponse(path, media_type='application/json',
                            headers={'Content-Encoding': 'gzip', 'Cache-Control': 'no-store'})

    @router.get('/space/blocks/{block}.svg')
    def space_block(block: int):
        from .cofdraw_family_screen import block_svg
        try:
            return Response(block_svg(block),media_type='image/svg+xml',headers={'Cache-Control':'private, max-age=3600'})
        except KeyError as exc:
            raise HTTPException(404,'Unknown building block') from exc

    @router.post('/space/export/{family}')
    def space_export(family: str):
        from .cofdraw_family_screen import export_family, screen
        if family not in screen()['families']:
            raise HTTPException(404,'No screened pairs in this family')
        key=uuid.uuid4().hex
        directory=outputs/key
        directory.mkdir(parents=True)
        export_family(family,directory/(key+'.zip'))
        return {'url':f'/api/drawing/exports/{key}/zip'}

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
        if not notation and recipe.get('explicit_graph'):
            return value
        if " COFex." not in notation:
            raise ValueError("Catalog source is missing its full SMILES-based COFex notation")
        expected = {canonical(c["monomer_smiles"]) for c in recipe["components"].values()}
        named = {canonical(c["smiles"]) for c in value.get("monomers", [])}
        encoded = {canonical(s) for s in notation.split(" COFex.", 1)[0].split(".")}
        if expected != named or expected != encoded:
            raise ValueError("Catalog names/SMILES/COFex precursor sets do not agree")
        return value

    @router.post("/cofex/resolve")
    def resolve_cofex(request: CofexInput):
        from .cofdraw_cofex import parse_record, phase_record
        try:
            text = request.text.strip()
            if not text:
                raise ValueError('Paste a COFex record first.')
            legacy = {}
            if text.startswith('{'):
                legacy = json.loads(text)
                if not isinstance(legacy,dict):
                    raise ValueError('Paste a COFex notation or a copied COFex record.')
                entry=legacy.get('drawing_cofex')
                text = entry['text'] if isinstance(entry,dict) and isinstance(entry.get('text'),str) else next((value for field in ('drawing_cofex','cofex_notation','design_cofex_notation','cofex_id')
                             if isinstance(value:=legacy.get(field),str) and value), '')
                if not text:
                    raise ValueError('The copied record contains no assigned COFex notation.')
            catalog = read_catalog(root)
            matches = []
            def matches_source(identity, notation):
                return text == identity or (notation and text in (notation, notation.split(';',1)[0]))
            for row in source_records():
                key = row['recipe_id']
                meta = catalog.get('source_metadata', {}).get(key, {})
                if matches_source(row.get('cofex_id'), meta.get('cofex_notation')):
                    value, digest = source_recipe(key)
                    meta = source_metadata(value,digest)
                    if matches_source(value.get('cofex_id'),meta.get('cofex_notation')):
                        matches.append(dict(source_recipe_id=key, label=meta.get('literature_name') or key))
            monomers = catalog_monomers(catalog)
            manifest = crystal_manifest()
            for item in manifest.get('reported_models',[]) + manifest.get('curated_models',[]):
                phase = phase_record(item,monomers)
                if phase and text in (phase['text'], phase['text'].split(';',1)[0]):
                    matches.append(dict(label=item['model_name']+' · published phase',crystal_project={
                        'template_id':item['id'],'version':2,'reference_only':True}))
            choices = matches
            if not matches:
                choices = parse_record(text, monomers)
            elif ' COFex.' in text and ';' not in text:
                # Multiple phases may use the same precursors. Also offer an
                # editable design when this compact notation has a builder.
                try:
                    choices.extend(parse_record(text, monomers))
                except ValueError:
                    pass
            if legacy.get('design_cofex_notation') and not legacy.get('drawing_cofex'):
                template = legacy.get('selected_template')
                linkage = legacy.get('linkage')
                choices = [p for p in choices if ('prepare_request' not in p or
                    ((not linkage or p['prepare_request']['linkage']==linkage) and
                     (not isinstance(template,dict) or p['prepare_request']['template_id']==template.get('id'))))]
                if not choices:
                    raise ValueError('Saved linkage and connectivity disagree with the COFex notation.')
            return {'choices': choices}
        except (ValueError, KeyError, RuntimeError, TypeError, RecursionError) as exc:
            raise HTTPException(422,str(exc)) from exc

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
        from .cofdraw_synthesis import synthetic_accessibility
        from .cofdraw_describe import NAMED_PAIRS, NAMED_DESIGNS
        from .cofdraw_mixed import MIXED_EXAMPLES, mixed_pore_patterns
        from .cofdraw_self import self_examples
        monomers = decorated_monomers(catalog_monomers(catalog))
        monomers = [{**m, 'synthetic_accessibility': synthetic_accessibility(m.get('smiles', ''))} for m in monomers]
        from .cofdraw_starters import topology_catalog
        return {"schema_version": VERSION, "linkages": library_linkages(catalog), "source_recipes": sources,
                "templates": TEMPLATES, "topology_catalog": topology_catalog(catalog.get("topology_catalog", [])),
                "style_controls": {"hydrogen_labels":["show","compact"], "element_schemes": ["style_default", *ELEMENT_SCHEMES], "palettes":list(PALETTES),"palette_colors":PALETTES,"fonts":FONTS,"edge_styles":["wavy_cross","dashed","crossbar","bond_stub"],
                                   "element_colors":{scheme:{Chem.GetPeriodicTable().GetElementSymbol(z):color for z,color in colors.items()} for scheme,colors in ELEMENT_SCHEMES.items()},
                                   "defaults":{s['id']:style_configuration(s['id'])[2] for s in STYLES},
                                   "abbreviation_modes":["explicit","condensed","OMe"],"layer_count":[2,5],"layer_opacity":[.05,.7],
                                   "fill_opacity":[.05,.6],"font_scale":[.6,1.6],"bond_scale":[.5,2.5]},
                "styles": STYLES, "monomer_presets": monomers,
                "assembly_examples": [*self_examples(), *MIXED_EXAMPLES],
                "mixed_pore_patterns": mixed_pore_patterns(),
                "cof_presets": catalog.get("cof_presets", PRESETS), "catalog_version": catalog.get("catalog_version", "builtin-1"),
                "named_pairs": [{"name":name,"a":a,"b":b,"linkage":link,"template_id":template,
                                 "doi":NAMED_DESIGNS[name][1]} for name,(a,b,link,template) in NAMED_PAIRS.items()],
                "style_evidence": catalog.get("style_evidence", catalog.get("source_styles", [])),
                "source_registry_errors": errors, "max_atoms": 1000,
                "presets": [{"instances": n, "label": f"{n // 2} A + {n // 2} B alternating pore"} for n in (6, 12)],
                "scope": "Finite chemical pore patches with explicit attachment ports; preset connectivity is a design choice, not a topology prediction."}

    @router.get('/describe/catalog')
    def description_catalog():
        from .cofdraw_describe import describe_catalog
        return describe_catalog(library())

    @router.post('/describe')
    def describe_drawing(request: DescribeDrawing):
        from .cofdraw_describe import describe_cof, describe_catalog
        text = request.text.strip()
        if not text:
            raise HTTPException(422, 'Describe a COF or choose a built-in example.')
        if request.mode == 'rules' and request.llm is not None:
            raise HTTPException(422, 'API settings are used only in API mode.')
        if request.mode == 'llm' and (request.llm is None or request.selection is not None):
            raise HTTPException(422, 'API mode needs API settings and no manual selection. Review selections locally in rules mode.')
        catalog = library()
        selection = request.selection.model_dump(exclude_none=True) if request.selection else None
        proposal = None
        if request.mode == 'llm':
            from .cofdraw_llm import select_with_llm
            try:
                proposal = select_with_llm(text, describe_catalog(catalog), base_url=request.llm.base_url, provider=request.llm.provider,
                                           model=request.llm.model,
                                           api_key=request.llm.api_key.get_secret_value())
                selection = proposal['selection']
            except (ValueError, RuntimeError) as exc:
                # Provider details and request credentials never enter stored reports.
                raise HTTPException(422, str(exc)) from None
        try:
            result = describe_cof(text, catalog, selection=selection,
                                  allow_unmatched_words=request.mode == 'llm' or bool(selection))
        except (ValueError, KeyError, TypeError, RuntimeError) as exc:
            raise HTTPException(422, str(exc)) from None
        result['interpretation'] = {'mode': request.mode}
        if proposal is not None:
            result['interpretation'].update(model=request.llm.model, explanation=proposal.get('explanation', ''))
            unresolved = proposal.get('unresolved', [])
            if unresolved:
                result['unresolved'] = list(dict.fromkeys(result.get('unresolved', []) + unresolved))
                if result['status'] == 'ready':
                    result['status'] = 'needs_input'
                result['prepare_request'] = None
                result['source_recipe_id'] = None
        # Description plans are transient. Keys, provider responses and request
        # bodies are not written to the scientific report store or project files.
        return result

    @router.post('/describe/chat')
    def describe_conversation(request: DescribeConversation):
        from .cofdraw_conversation import checked_context, describe_turn, model_turn_text, _help
        from .cofdraw_describe import describe_catalog
        text = request.text.strip()
        if not text:
            raise HTTPException(422, 'Write a message or choose an example.')
        if request.mode == 'rules' and request.llm is not None:
            raise HTTPException(422, 'API settings are used only in API mode.')
        if request.mode == 'llm' and (request.llm is None or request.selection is not None):
            raise HTTPException(422, 'API mode needs API settings and no manual selection.')
        current_library = library()
        context = request.context.model_dump(exclude_none=True) if request.context else None
        selection = request.selection.model_dump(exclude_none=True) if request.selection else None
        proposal = None
        try:
            context = checked_context(context, current_library)
            # Never echo a known credential through a user message, context,
            # exception, or the transient conversation response.
            secret = request.llm.api_key.get_secret_value() if request.llm else ''
            if secret and secret in json.dumps({'text': text, 'context': context}, ensure_ascii=False):
                raise ValueError('Keep API credentials in the password field, separate from chat messages.')
            local_help = _help(text) is not None or bool(re.fullmatch(r'(?:start over|start again|new conversation|reset conversation)[.! ]*', text, re.I))
            if request.mode == 'llm' and not local_help:
                from .cofdraw_llm import select_with_llm
                proposal = select_with_llm(model_turn_text(text, context, current_library), describe_catalog(current_library),
                    base_url=request.llm.base_url, provider=request.llm.provider, model=request.llm.model, api_key=secret)
                selection = proposal['selection']
            result = describe_turn(text, current_library, context=context, selection=selection,
                allow_unmatched_words=bool(selection) and (request.mode == 'llm' or request.selection is not None))
        except (ValueError, KeyError, TypeError, RuntimeError) as exc:
            raise HTTPException(422, str(exc)) from None
        result['interpretation'] = {'mode': 'rules' if local_help else request.mode}
        if proposal is not None:
            result['interpretation'].update(model=request.llm.model, explanation=proposal.get('explanation', ''))
            if proposal.get('unresolved'):
                result['unresolved'] = list(dict.fromkeys(result['unresolved'] + proposal['unresolved']))
                result['status'] = 'needs_input' if result['status'] == 'ready' else result['status']
                result['prepare_request'] = None; result['source_recipe_id'] = None
                pending = result['conversation']['context']['pending']
                result['conversation']['context']['pending'] = list(dict.fromkeys(pending + proposal['unresolved']))[:12]
        # No conversation/session store; the browser retains this in memory only.
        return result

    @router.get("/sources/{recipe_id}")
    def source(recipe_id: str):
        try:
            value, digest = source_recipe(recipe_id)
            metadata = source_metadata(value, digest)
            rule=crystal_manifest().get('source_layer_designs',{}).get(recipe_id)
            if rule and rule['source_sha256']==digest:
                metadata={**metadata,'periodic_design_inputs':{k:rule[k] for k in ('precursor_components','linkage','template')}}
            return {"recipe": value, "source_sha256": digest, "status": "source_bound_provisional", "metadata": metadata}
        except (ValueError, KeyError, OSError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post("/monomer")
    def monomer(request: MonomerDrawing):
        try:
            from .cofdraw_synthesis import synthetic_accessibility
            result = normalize_monomer(**request.model_dump())
            result['synthetic_accessibility'] = synthetic_accessibility(result['smiles'])
            return result
        except (ValueError, KeyError, TypeError, RuntimeError, IndexError, OverflowError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.get("/templates/installed")
    def installed_template_collections():
        from .cofdraw_template_sheet import installed_sheets
        return {"collections": list(installed_sheets())}

    @router.get("/templates/installed/{name}")
    def installed_templates(name: str):
        from .cofdraw_template_sheet import installed_sheets, read_binary_templates
        source = installed_sheets().get(name)
        if source is None:
            raise HTTPException(404, 'This ChemDraw collection is unavailable.')
        try:
            return read_binary_templates(source.read_bytes())
        except (ValueError, RuntimeError, OSError) as exc:
            raise HTTPException(422, 'The installed molecular palette could not be read.') from exc

    @router.post("/templates/binary")
    def import_binary_templates(request: BinaryTemplateSheet):
        import base64
        from .cofdraw_template_sheet import read_binary_templates
        try:
            return read_binary_templates(base64.b64decode(request.content_base64, validate=True))
        except (ValueError, RuntimeError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post("/templates/cdxml")
    def import_template_sheet(request: TemplateSheet):
        from .cofdraw_template_sheet import read_templates
        try:
            return read_templates(request.cdxml)
        except (ValueError, RuntimeError) as exc:
            raise HTTPException(422, str(exc)) from exc

    def custom_assembly_row(precursor):
        from rdkit import Chem
        mol = Chem.MolFromSmiles(precursor.smiles)
        if mol is None or not 1 <= mol.GetNumAtoms() <= 250 or len(Chem.GetMolFrags(mol)) != 1:
            raise ValueError("Choose one connected precursor with up to 250 atoms.")
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        key = 'custom-' + hashlib.sha256(smiles.encode()).hexdigest()[:24]
        row = decorated_monomers([{'id':key, 'label':precursor.label, 'abbreviation':precursor.label, 'smiles':smiles}])[0]
        if not row.get('assembly_roles'):
            raise ValueError("The precursor could not be checked for assembly.")
        return row

    @router.post("/assembly/precursor")
    def inspect_assembly_precursor(request: AssemblyPrecursor):
        try:
            return custom_assembly_row(request)
        except (ValueError, RuntimeError, IndexError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post("/prepare/assembly")
    def prepare_assembly(request: AssemblyDrawing):
        from .cofdraw_assembly import assemble
        try:
            rows = {r['id']:r for r in decorated_monomers(
                [r for r in catalog_monomers(read_catalog(root)) if r['id'] in request.monomer_ids])}
            for precursor in request.custom_precursors:
                row = custom_assembly_row(precursor)
                if row['id'] in rows or row['id'] not in request.monomer_ids:
                    raise ValueError("Each custom precursor must match one selected building block.")
                rows[row['id']] = row
            if any(k not in rows for k in request.monomer_ids):
                raise ValueError('Choose a recognized library entry.')
            if len(request.monomer_ids)==1 and request.linker_sequence:
                raise ValueError('A self-condensing design has no mixed-linker sequence.')
            recipe = assemble([rows[k] for k in request.monomer_ids], request.linker_sequence,
                              pore_view=request.pore_view, pore_layout=request.pore_layout)
            recipe['assembly']['precursors'] = [{'id':k,'label':rows[k]['abbreviation'], 'smiles':rows[k]['canonical_smiles']} for k in request.monomer_ids]
            if request.custom_precursors:
                recipe['assembly']['custom_precursors'] = [row.model_dump() for row in request.custom_precursors]
            recipe.update(schema_version=VERSION,recipe_id='designed-assembly',designed_variant=True,ground_truth=False,
                          review_status='user_design',layout_refinement='bond_angle_pose_v2',design_note=recipe['assembly']['scope'])
            recipe['design_template_signature']=design_signature(recipe)
            refresh_design_identity(recipe)
            validate_recipe(recipe)
            return {'recipe':recipe,'status':'designed_variant','limits':recipe['design_note']}
        except (ValueError, KeyError, TypeError, RuntimeError) as exc:
            raise HTTPException(422,str(exc)) from exc

    @router.post("/prepare")
    def prepare(request: PrepareDrawing):
        try:
            prepared = prepare_template(**request.model_dump())
            recipe = {"schema_version": VERSION, "recipe_id": "designed-pore", **prepared,
                      "topology": None, "topology_basis": "explicit_user_chosen_template_not_experimental",
                      "layout_refinement": "bond_angle_pose_v2",
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
                from .cofdraw_explicit import clear_source_registration
                clear_source_registration(recipe)
                validate_recipe(recipe)
                if not recipe.get('explicit_graph'):
                    refresh_design_identity(recipe)
                status = "designed_variant"
                metadata = {}
            validate_recipe(recipe)
            extent_options = drawing_extents(recipe, request.source_recipe_id, source_hash)
            drawing_id = uuid.uuid4().hex
            original_id = recipe.get("recipe_id")
            recipe["recipe_id"] = drawing_id  # filenames are always server-generated
            directory = outputs / drawing_id
            with lock:
                prepared_patch = build_patch(recipe)
                audit = render_styled(recipe, directory, request.size, request.style_id, request.style_seed, request.style_options,
                                      prepared_patch=prepared_patch, recent_styles=request.recent_styles)
                # Export support is explicit; the underlying MOL remains available.
                extra_formats = export_chemical_formats(directory, drawing_id, source_mol=prepared_patch[0])
            export_report = directory / f"{drawing_id}.chemical-export.json"
            if export_report.is_file():
                audit["chemical_exports"] = json.loads(export_report.read_text(encoding="utf-8"))
            audit["graph_scope"] = "finite_source_pore_patch_with_continuation_ports" if source_hash else "finite_designed_pore_patch_with_continuation_ports"
            audit.update(drawing_status=status, registered_recipe_sha256=source_hash, original_recipe_id=original_id)
            from .cofdraw_synthesis import precursor_scores
            audit['precursor_synthetic_accessibility'] = precursor_scores(recipe)
            audit.update(design_cofex_notation=recipe.get("design_cofex_notation"), design_topology=recipe.get("design_topology"),
                         design_identity_qualification=recipe.get("design_identity_qualification"), selected_template=recipe.get("selected_template"))
            audit.update(cofex_string=metadata.get("cofex_notation"), cofex_notation=metadata.get("cofex_notation"),
                         literature_metadata=metadata, monomers=[{"component": name, "label": c.get("label", name), "smiles": c.get("monomer_smiles")} for name, c in recipe.get("components",{}).items()],
                         chemical_export_scope="Editable finite patch with continuation atoms. SVG styles are separate presentation; no complete periodic unit cell is implied.")
            from .cofdraw_cofex import recipe_record
            audit['drawing_cofex'] = recipe_record(recipe, metadata, source=bool(source_hash))
            urls = {ext: f"/api/drawing/exports/{drawing_id}/{ext}" for ext in ("svg", "png", "mol", "json", *extra_formats)}
            bundle = {"recipe": recipe, "audit": audit, "downloads": urls}
            (directory / f"{drawing_id}.json").write_text(json.dumps(bundle, indent=2) + "\n", encoding="utf-8")
            return {"drawing_id": drawing_id, "status": status, "audit": audit, "downloads": urls,
                    "repeat_options": extent_options}
        except (ValueError, KeyError, TypeError, RuntimeError, OSError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.get("/exports/{drawing_id}/{extension}")
    def artifact(drawing_id: str, extension: Literal["svg", "png", "mol", "json", "sdf", "cdxml", "cif", "xyz", "pdb", "csv", "zip"]):
        if not re.fullmatch(r"[0-9a-f]{32}", drawing_id):
            raise HTTPException(404, "Drawing not found")
        path = (outputs / drawing_id / f"{drawing_id}.{extension}").resolve()
        if not path.is_relative_to(outputs.resolve()) or not path.is_file():
            raise HTTPException(404, "Drawing not found")
        types = {"svg": "image/svg+xml", "png": "image/png", "mol": "chemical/x-mdl-molfile", "json": "application/json", "sdf": "chemical/x-mdl-sdfile", "cdxml": "chemical/x-cdxml", "cif":"chemical/x-cif", "xyz":"chemical/x-xyz", "pdb":"chemical/x-pdb"}
        types['csv'] = 'text/csv; charset=utf-8'
        types['zip'] = 'application/zip'
        return FileResponse(path, media_type=types[extension], filename=f"cofdraw-{drawing_id[:8]}.{extension}")

    def selected_recipe(request):
        value = source_recipe(request.source_recipe_id)[0] if request.source_recipe_id else deepcopy(request.recipe)
        if not request.source_recipe_id:
            value.update(cofex_id=None, cofex_notation=None, cofex_string=None,
                         parent_cofex_id=None, source=None, reference=None,
                         topology=None, topology_evidence=None, designed_variant=True,
                         ground_truth=False, review_status='user_design')
            from .cofdraw_explicit import clear_source_registration
            clear_source_registration(value)
        validate_recipe(value)
        return value

    def drawing_extents(recipe, source_id, digest):
        from .cofdraw_repeat import repeat_options
        rule = crystal_manifest().get('source_layer_designs', {}).get(source_id) if source_id else None
        return repeat_options(recipe, source_rule=rule, source_sha256=digest)

    @router.post('/repeat-options')
    def available_repeats(request: RenderDrawing):
        try:
            recipe = selected_recipe(request)
            digest = source_recipe(request.source_recipe_id)[1] if request.source_recipe_id else None
            return drawing_extents(recipe, request.source_recipe_id, digest)
        except (ValueError, KeyError, TypeError, RuntimeError, OSError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post('/repeat')
    def repeat_drawing(request: RepeatDrawing):
        from .cofdraw_repeat import write_repeat
        try:
            if request.source_recipe_id:
                recipe, digest = source_recipe(request.source_recipe_id)
                rule = crystal_manifest().get('source_layer_designs', {}).get(request.source_recipe_id)
            else:
                recipe, digest, rule = selected_recipe(request), None, None
            key = uuid.uuid4().hex
            with lock:
                return write_repeat(recipe, outputs / key, key, request, source_rule=rule, source_sha256=digest)
        except (ValueError, KeyError, TypeError, RuntimeError, OSError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post('/repeat/{drawing_id}/measure')
    def measure_repeat(drawing_id: str):
        from .cofdraw_measure import measure_aperture, qualify_circuit_span
        if not re.fullmatch(r'[0-9a-f]{32}',drawing_id):
            raise HTTPException(404,'Drawing not found')
        path=(outputs/drawing_id/(drawing_id+'.json')).resolve()
        if not path.is_relative_to(outputs.resolve()) or not path.is_file():
            raise HTTPException(404,'Drawing not found')
        record=json.loads(path.read_bytes())
        audit=record.get('audit',{})
        if audit.get('graph_scope')!='finite_designed_connected_repeat' or 'drawing_aperture' not in audit:
            raise HTTPException(422,record.get('summary',{}).get('measurement_reason') or 'Choose a drawing with a resolved pore boundary.')
        layout=audit['drawing_style']['layout']
        try:
            result=measure_aperture(layout['atom_coordinates'],audit['drawing_aperture'])
        except ValueError as exc:
            raise HTTPException(422,str(exc)) from exc
        result['focal_face']=record['summary']['measurement_face']
        qualify_circuit_span(result,record['summary'].get('template'))
        result['method']='Inscribed-circle sampling within the selected molecular pore rim, using all explicit atom sites and their van der Waals radii.'
        return {'measurement':result,'audit':{'drawing_style':{'layout':layout}},'drawing_id':drawing_id}

    def scientific_result(result, operation='calculation', inputs=None):
        from . import __version__
        result=deepcopy(result)
        contents=result.pop('exports',{})
        identifier=uuid.uuid4().hex
        directory=outputs/identifier;directory.mkdir(parents=True)
        urls={}
        for extension,content in contents.items():
            if extension not in {'cif','xyz','pdb','sdf','mol','csv'}:continue
            if not isinstance(content,str) or len(content)>20_000_000:raise ValueError('Structure export exceeds the size limit')
            (directory/f'{identifier}.{extension}').write_text(content,encoding='utf-8',newline='')
            urls[extension]=f'/api/drawing/exports/{identifier}/{extension}'
        urls['json']=f'/api/drawing/exports/{identifier}/json'
        result['downloads']=urls
        result.update(report_id=identifier, operation=operation, engine_version=__version__,
                      report_schema_version='cofdraw-science-1',
                      created_at_utc=datetime.now(timezone.utc).isoformat())
        if inputs is not None:
            result['inputs']=deepcopy(inputs)
            result['input_sha256']=hashlib.sha256(json.dumps(inputs,sort_keys=True,allow_nan=False).encode()).hexdigest()
        report_text=json.dumps(result,indent=2,allow_nan=False)
        if len(report_text.encode('utf-8'))>40_000_000:raise ValueError('Scientific report exceeds the 40 MB limit')
        temporary=directory/f'{identifier}.json.tmp'
        temporary.write_text(report_text,encoding='utf-8')
        temporary.replace(directory/f'{identifier}.json')
        return result

    @router.post('/layout')
    def layout(request: LayoutDrawing):
        try:
            from .cofdraw_layout import clean_layout
            return clean_layout(selected_recipe(request))
        except (ValueError,KeyError,TypeError,RuntimeError) as exc:raise HTTPException(422,str(exc)) from exc

    def crystal_manifest():
        path=root/'crystal-models/manifest.json'
        catalog=json.loads(path.read_text(encoding='utf-8')) if path.is_file() else {'templates':[]}
        designs=path.parent/'diamond/models.json'
        catalog['linker_models']=json.loads(designs.read_bytes()) if designs.is_file() else []
        series=path.parent/'linker-series.json'
        rows={r['id']:r for group in ('templates','reported_models','curated_models','linker_models') for r in catalog.get(group,[])}
        if series.is_file():
            for family in json.loads(series.read_bytes()):
                choices={key:label for key,label in family['choices'].items() if key in rows}
                if family.get('include_diamond_designs'):
                    choices.update({r['id']:r['linker_label']+' · single-net design' for r in catalog['linker_models']})
                for key in choices:
                    rows[key]['linker_series']={'shared_node':family['shared_node'],'choices':choices,'sources':family['sources']}
        audit=path.parent/'aryl-substitutions.json'
        if audit.is_file():
            substitutions=json.loads(audit.read_bytes())
            for group in ('templates','reported_models','curated_models','linker_models'):
                for row in catalog.get(group,[]):
                    match=substitutions.get(row['id'],{})
                    row['aryl_f_sites']=match.get('substituted_sites',0) if match.get('parent_sha256')==row['sha256'] else 0
        from .cofdraw_starters import crystal_starters
        catalog['starter_designs'] = crystal_starters(catalog)
        return catalog

    def load_crystal_template(item, manifest):
        from .cofdraw_science import import_cif
        directory = (root/'crystal-models').resolve()
        path = (directory/item['path']).resolve()
        if not path.is_relative_to(directory) or hashlib.sha256(path.read_bytes()).hexdigest()!=item['sha256']:
            raise ValueError('Crystal template integrity check failed')
        structure = import_cif(path.read_text(encoding='utf-8'))
        structure['provenance'].update(template=item['id'], generator=manifest['generator'],
            model_name=item['model_name'], doi=manifest['citation_doi'], qualification=item['description'])
        if 'citation_doi' in item:
            structure['provenance'].update(kind='curated_crystal' if 'dataset_id' in item else 'reported_crystal',generator=item['generator'],
                doi=item['citation_doi'],network_dimension=item['network_dimension'])
            for field in ('data_doi','license','attribution','precursor_match','dataset_id','upstream_url'):
                if field in item:
                    structure['provenance'][field]=item[field]
        if 'linker_label' in item:
            structure['provenance'].update(kind='reticular_design', network_dimension=3,
                generator=item['generator'], monomer_ids=item['monomer_ids'])
        if item.get('cofex_phase'):
            from .cofdraw_cofex import phase_record
            record = phase_record(item,catalog_monomers(read_catalog(root)))
            if record:
                structure['provenance']['drawing_cofex'] = record
        return structure

    def crystal_exports(item,structure,repeats):
        from .cofdraw_science import structure_exports
        from .cofdraw_network import inspect_periodic_network, label_periodic_networks
        if item is None:
            # Imported CIFs can exceed the contact-search budget while remaining
            # usable coordinates. Preserve ordinary import in that case.
            try:
                network=inspect_periodic_network(structure)
            except ValueError as exc:
                result=structure_exports(structure,tuple(repeats))
                result.setdefault('warnings',[]).append(f'Connectivity illustration unavailable: {exc}')
                return result
        else:
            network=inspect_periodic_network(structure)
        result=structure_exports(structure,tuple(repeats),periodic_bonds=network['edges'])
        result['display_edges'] = network['edges']
        result['color_periodic_networks'] = 3 in network['dimensions']
        if result['color_periodic_networks']:
            label_periodic_networks(result['viewer'],structure,repeats,network['edges'])
            from .cofdraw_crystal_figures import adjacent_circuits, chemical_drawing, crystal_details, junction_patch
            result['crystal_details'] = crystal_details(structure, network['edges'])
            if item is not None or structure['provenance'].get('kind') in ('linker_design','monomer_dia_design'):
                result['chemical_drawing_svg'] = chemical_drawing(structure, network['edges'])
                star = junction_patch(structure, network['edges'])
                if star:
                    try:
                        result['chemical_junction_svg'] = chemical_drawing(structure,network['edges'],circuit=star)
                        result['junction_patch'] = star
                    except ValueError as exc:
                        result.setdefault('warnings',[]).append(str(exc))
                circuit = result['crystal_details'].get('ring')
                if circuit:
                    result['chemical_circuit_svg'] = chemical_drawing(structure, network['edges'], circuit=circuit)
                    joined=adjacent_circuits(circuit,structure['cell'])
                    if joined:
                        result['chemical_three_circuits_svg']=chemical_drawing(structure,network['edges'],circuit=joined)
                    result['circuit_extents']=[1,2,3,4,7] if joined else [1]
                    try:
                        result['chemical_framework_svg'] = chemical_drawing(structure, network['edges'], circuit=circuit, scheme=True)
                    except ValueError as exc:
                        result.setdefault('warnings',[]).append(f'Conventional circuit layout unavailable: {exc}')
                elif item is not None:
                    from .cofdraw_crystal_circuits import reference_framework
                    patches,topology=reference_framework(structure,network['edges'])
                    if patches:
                        result['reference_circuit_patches']={str(n):patch for n,patch in patches.items()}
                        result['circuit_extents']=sorted(patches)
                        result['chemical_circuit_svg']=chemical_drawing(structure,network['edges'],circuit=patches[1])
                        result['chemical_framework_svg']=result.get('chemical_junction_svg') or chemical_drawing(structure,network['edges'],circuit=patches[1],scheme=True)
                    if topology:
                        result['crystal_details']['periodic_topology']=topology
        if item is None:
            return result
        result['scope']=item['description']
        credit='; '.join(str(item[k]) for k in ('model_name','citation_doi','data_doi','attribution','license') if k in item)
        result['exports']['cif']=result['exports']['cif'].replace('data_COFdraw\n','data_COFdraw\n# '+credit+'\n',1)
        return result

    def reported_result(item,manifest,request):
        result=crystal_exports(item,load_crystal_template(item,manifest),request.repeats)
        result.update(available=True,cell_validation={'reported_model':item['id'],'source_sha256':item['sha256']})
        return scientific_result(result,'build-periodic',request.model_dump())

    @router.get('/science/library')
    def scientific_library():
        data=crystal_manifest()
        return {'templates':data['templates'],'reported_models':data.get('reported_models',[]),'linker_models':data['linker_models'],'starter_designs':data['starter_designs'],
                'curated_models':data.get('curated_models',[]),'methods':['RDKit ETKDG / MMFF or UFF','Periodic probe geometry','BET fit from supplied isotherm']}

    @router.get('/science/synthesis-library')
    def synthesis_library():
        from .cofdraw_synthesis import catalog_scores
        return catalog_scores(catalog_monomers(read_catalog(root)))

    @router.post('/science/periodic-options')
    def matched_periodic_options(request: RenderDrawing):
        try:
            from .cofdraw_periodic_options import periodic_options
            manifest=crystal_manifest()
            digest=source_recipe(request.source_recipe_id)[1] if request.source_recipe_id else None
            with lock:return periodic_options(selected_recipe(request),manifest,source_sha256=digest)
        except (ValueError,KeyError,TypeError,RuntimeError,OSError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.post('/science/model')
    def model(request: ModelDrawing):
        try:
            from .cofdraw_science import import_cif, structure_exports
            if bool(request.template_id)==bool(request.cif_text):raise ValueError('Choose one model template or import one CIF')
            if request.linker_edit and (not request.template_id or request.recipe_match):
                raise ValueError('Choose a catalog crystal to start a separate linker design.')
            if request.template_id:
                manifest=crystal_manifest()
                matches=[x for x in manifest['templates']+manifest.get('reported_models',[])+manifest.get('curated_models',[])+manifest['linker_models'] if x['id']==request.template_id]
                if len(matches)!=1:raise ValueError('Unknown crystal template')
                item=matches[0];structure=load_crystal_template(item,manifest)
                if request.recipe_match:
                    from .cofdraw_periodic_options import periodic_options
                    match_id=request.recipe_match.source_recipe_id
                    options=periodic_options(selected_recipe(request.recipe_match),manifest,
                                             source_sha256=source_recipe(match_id)[1] if match_id else None)
                    if not any(o['id']==item['id'] for o in options['options']):raise ValueError('This periodic model does not match the current precursor pair, linkage and unchanged connectivity template')
                    structure['provenance']['drawing_match']={**item['precursor_match'],'scope':options['scope']}
            else:
                if request.recipe_match:raise ValueError('Recipe matching applies only to validated built-in periodic designs')
                structure=import_cif(request.cif_text)
            with lock:
                if request.linker_edit:
                    from .cofdraw_crystal_linkers import fluorinated_model
                    structure=fluorinated_model(structure,request.template_id)
                    result=crystal_exports(None,structure,request.repeats)
                    result['scope']=structure['provenance']['qualification']
                    result['exports']['cif']=result['exports']['cif'].replace('data_COFdraw\n',
                        'data_COFdraw\n# Unrelaxed aryl-F design from '+request.template_id+'; inherited cell.\n',1)
                elif request.template_id:
                    result=crystal_exports(item,structure,request.repeats)
                else:
                    result=crystal_exports(None,structure,request.repeats)
                    result['scope']='Periodic atomic model with the supplied CIF geometry.'
            return scientific_result(result, 'model', request.model_dump())
        except (ValueError,KeyError,TypeError,RuntimeError,OSError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.post('/science/build-periodic')
    def build_periodic(request: PeriodicFromDrawing):
        try:
            from .cofdraw_periodic_build import build_periodic_layer
            from .cofdraw_periodic_options import periodic_options, validated_design
            from .cofdraw_science import structure_exports

            recipe = selected_recipe(request)
            manifest = crystal_manifest()
            with lock:
                from .cofdraw_spatial import SPATIAL
                key=(recipe.get('selected_template') or {}).get('id')
                if key in SPATIAL or key=='branched_4_3':
                    return {'available':False,'reason':'This spatial design has a building-block net. Atomic coordinates and a validated cell are required for a periodic CIF.'}
                source_digest=source_recipe(request.source_recipe_id)[1] if request.source_recipe_id else None
                source_rule=manifest.get('source_layer_designs',{}).get(request.source_recipe_id)
                if request.source_recipe_id and request.layer_spacing is None and not source_rule:
                    exact=next((r for r in manifest.get('reported_models',[])
                                if r.get('source_recipe_id')==request.source_recipe_id
                                and r.get('source_recipe_sha256')==source_digest),None)
                    if exact:
                        return reported_result(exact,manifest,request)
                try:
                    design = validated_design(recipe,source_rule=source_rule,source_sha256=source_digest)
                    options = periodic_options(recipe, manifest, design=design)
                except (ValueError, KeyError, TypeError) as exc:
                    return {'available': False, 'reason': str(exc)}
                basis=design[0].get('source_drawing_basis')
                matches = {o['id'] for o in options['options']}
                reported = next((entry for entry in manifest.get('reported_models', [])
                    if entry.get('precursor_match',{}).get('linkage') == design[0]['linkage']
                    and design[0]['selected_template']['id'] in entry['precursor_match']['drawing_template_ids']
                    and sorted(entry['precursor_match']['canonical_smiles']) == design[1]), None)
                if reported and request.layer_spacing is None and not basis:
                    return reported_result(reported,manifest,request)
                item = next((entry for entry in manifest['templates']
                             if entry['id'] in matches and entry['model_name'].endswith('-AA')), None)
                if item and request.layer_spacing is None and not basis:
                    structure = load_crystal_template(item, manifest)
                    structure['provenance']['drawing_match'] = item['precursor_match']
                    structure['provenance']['spacing_angstrom'] = structure['cell'][2][2]
                    validation = {'matched_template': item['id'], 'template_cif_sha256': item['sha256']}
                else:
                    try:
                        structure, validation = build_periodic_layer(design[0], spacing=request.layer_spacing)
                    except ValueError as exc:
                        return {'available': False, 'reason': str(exc)}
                if basis:
                    structure['provenance']['source_drawing_basis']=basis
                    validation['source_recipe_sha256']=source_digest
                result = structure_exports(structure, repeats=tuple(request.repeats), periodic_bonds=validation.get("periodic_bonds"))
                if validation.get("periodic_bonds"):
                    from .cofdraw_crystal_figures import precursor_scenes
                    result["crystal_details"]=precursor_scenes(structure,validation["periodic_bonds"])
            result.update(available=True, cell_validation=validation,
                scope='Unrelaxed periodic starting model of the selected precursor chemistry and planar connectivity. AA layer spacing is a model parameter.')
            provenance = structure['provenance']
            if basis:
                result['scope']=' '.join(filter(None,[result['scope'],basis['assumption'],basis.get('orientation_assumption')]))
                if 'reference_layer_spacing_angstrom' in basis:
                    result['scope']+=(f" The powder reference spacing is approximately {basis['reference_layer_spacing_angstrom']:.2f} Å. "
                                      f"This AA model uses {provenance['spacing_angstrom']:.2f} Å after interlayer contact screening.")
            result['exports']['cif'] = result['exports']['cif'].replace('data_COFdraw\n', (
                'data_COFdraw\n# COFdraw unrelaxed AA starting model.\n'
                + (f"# Source drawing {basis['recipe_id']} ({source_digest}); chosen periodic net and AA stacking.\n" if basis else '')
                + f"# Template {design[0]['selected_template']['id']}; layer spacing {provenance['spacing_angstrom']:.4f} angstrom.\n"
                + (f"# Fragment geometry: {provenance['fragment_geometry']}.\n" if 'fragment_geometry' in provenance else '')
                + (f"# Meso phenyl tilt {provenance['phenyl_tilt_degrees']} degrees.\n"
                   if provenance.get('phenyl_tilt_degrees') else '')), 1)
            return scientific_result(result, 'build-periodic', request.model_dump())
        except (ValueError, KeyError, TypeError, RuntimeError, OSError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post('/science/build-reticular')
    @router.post('/science/build-dia')
    def build_diamond_design(request: DiamondDesign):
        try:
            from .cofdraw_dia import build_dia
            with lock:
                structure,audit = build_dia(request.monomer_a, request.monomer_b, request.linkage,
                                            request.label_a, request.label_b, request.network_count, request.topology)
                result = crystal_exports(None, structure, request.repeats)
            result.update(available=True,cell_validation=audit,scope=structure["provenance"]["qualification"])
            result["exports"]["cif"] = result["exports"]["cif"].replace("data_COFdraw\n",
                f"data_COFdraw\n# Unrelaxed {request.topology} design from reacted precursors, {request.network_count} network(s).\n",1)
            return scientific_result(result,"build-dia",request.model_dump())
        except (ValueError,KeyError,TypeError,RuntimeError) as exc:
            raise HTTPException(422,str(exc)) from exc

    @router.post('/science/patch')
    def patch_model(request: Patch3D):
        try:
            from .cofdraw_science import patch_3d
            with lock:result=patch_3d(selected_recipe(request),layer_spacing=request.layer_spacing,layers=request.layers,stacking=request.stacking,slip=request.slip,twist_degrees=request.twist_degrees)
            result.setdefault('scope','Finite illustrative layers from scaled drawing coordinates; not a relaxed periodic crystal. X sites mark continuation boundaries.')
            return scientific_result(result, 'patch', request.model_dump())
        except (ValueError,KeyError,TypeError,RuntimeError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.post('/science/networks')
    def network_model(request: NetworkHypothesis):
        try:
            if len(json.dumps(request.structure)) > 2_000_000:raise ValueError('Structure input limit is 2 MB')
            from .cofdraw_arrangements import network_hypothesis
            with lock:result=network_hypothesis(request.structure,request.offsets,repeats=tuple(request.repeats))
            return scientific_result(result, 'networks', request.model_dump())
        except (ValueError,KeyError,TypeError,RuntimeError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.post('/science/monomer')
    def analyze_monomer(request: MonomerAnalysis):
        try:
            from .cofdraw_science import monomer_analysis
            with lock:result=monomer_analysis(request.smiles,conformers=request.conformers)
            from .cofdraw_synthesis import synthetic_accessibility
            result['synthetic_accessibility'] = synthetic_accessibility(request.smiles)
            return scientific_result(result, 'monomer', request.model_dump())
        except (ValueError,KeyError,TypeError,RuntimeError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.post('/science/solubility')
    def analyze_solubility(request: SolubilityScreen):
        from .cofdraw_solubility import screen_solubility
        from .framework_drawing import canonical
        try:
            catalog = read_catalog(root)
            monomers = {row['id']: row for row in catalog_monomers(catalog)}
            smiles = request.smiles
            if request.monomer_id is not None:
                entry = monomers.get(request.monomer_id)
                if entry is None:
                    raise ValueError('Unknown monomer ID')
                smiles = entry['smiles']
            loadings = [
                {**record, 'canonical_smiles': canonical(monomers[record['monomer_id']]['smiles'])}
                for record in catalog.get('recent_monomer_update', {}).get('preparation_loadings', [])
            ]
            result = screen_solubility(smiles, root / 'solubility',
                                       solvent_ids=request.solvent_ids, solvent_ml=request.solvent_ml,
                                       preparation_loadings=loadings)
            return scientific_result(result, 'solubility', request.model_dump())
        except OSError as exc:
            raise HTTPException(503, 'The offline solubility reference data could not be read') from exc
        except (ValueError, KeyError, TypeError, RuntimeError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post('/science/porosity')
    def analyze_porosity(request: PoreAnalysis):
        try:
            if len(json.dumps(request.structure))>2_000_000:raise ValueError('Structure input limit is 2 MB')
            from .cofdraw_science import pore_analysis
            with lock:result=pore_analysis(request.structure,probe_radius=request.probe_radius,samples=request.samples)
            return scientific_result(result, 'porosity', request.model_dump())
        except (ValueError,KeyError,TypeError,RuntimeError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.post('/science/geometry')
    def analyze_geometry(request: GeometryAnalysis):
        try:
            if len(json.dumps(request.structure)) > 2_000_000:
                raise ValueError('Structure input limit is 2 MB')
            from .cofdraw_porosity import geometric_porosity
            with lock:
                result = geometric_porosity(**request.model_dump())
            inputs=request.model_dump()
            if 'include_map' not in request.model_fields_set:
                inputs.pop('include_map')
            return scientific_result(result, 'geometry', inputs)
        except (ValueError, KeyError, TypeError, RuntimeError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post('/science/measure')
    def measure(request: RenderDrawing):
        try:
            from .cofdraw_measure import measure_drawing
            with lock:result=measure_drawing(selected_recipe(request))
            return scientific_result(result,'measure',request.model_dump())
        except (ValueError,KeyError,TypeError,RuntimeError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.post('/science/bet')
    def analyze_bet(request: BETAnalysis):
        try:
            from .cofdraw_science import fit_bet
            return scientific_result(fit_bet(**request.model_dump()), 'bet', request.model_dump())
        except (ValueError,KeyError,TypeError,RuntimeError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.get('/science/capabilities')
    def scientific_capabilities():
        from . import __version__
        from .cofdraw_knowledge import capability_inventory
        from .cofdraw_adsorption import AdsorptionRequest
        from .cofdraw_solubility import DATASET_DOI, SOLVENT_NAMES
        models={'model':ModelDrawing,'patch':Patch3D,'networks':NetworkHypothesis,'monomer':MonomerAnalysis,
                'solubility':SolubilityScreen,'porosity':PoreAnalysis,'geometry':GeometryAnalysis,'bet':BETAnalysis,'inspect':StructureInspection,
                'pxrd':PowderDrawing,'pxrd-peaks':PeakComparison,'reagent-balance':ReagentBalance,'screen':ScreenBatch,'query':KnowledgeQuery,
                'finite-inspect':FiniteInspection,'finite-powder':FinitePowder,'stacking':PeriodicStacking,'periodic-options':RenderDrawing,
                'measure':RenderDrawing,'build-periodic':PeriodicFromDrawing}
        return {'engine_version':__version__, 'schema_version':'cofdraw-capabilities-1',
                'operations':{name:{'method':'POST','path':f'/api/drawing/science/{name}',
                                    'request_schema':model.model_json_schema()} for name,model in models.items()},
                'workflows':capability_inventory(),
                'synthesis_library':{'method':'GET','path':'/api/drawing/science/synthesis-library'},
                'solubility':{'solvents':SOLVENT_NAMES,'dataset_doi':DATASET_DOI},
                'adsorption':{'engine':'/api/drawing/science/adsorption/engine',
                              'start':{'method':'POST','path':'/api/drawing/science/adsorption/jobs','request_schema':AdsorptionRequest.model_json_schema()},
                              'poll':'/api/drawing/science/adsorption/jobs/{job_id}',
                              'cancel':{'method':'POST','path':'/api/drawing/science/adsorption/jobs/{job_id}/cancel'},
                              'scope':'Asynchronous local RASPA runs; save input/output bundle and inspect convergence.'},
                'description':{'method':'POST','path':'/api/drawing/describe',
                               'catalog':'/api/drawing/describe/catalog',
                               'request_schema':DescribeDrawing.model_json_schema(),
                               'conversation':{'method':'POST','path':'/api/drawing/describe/chat',
                                               'request_schema':DescribeConversation.model_json_schema(),
                                               'scope':'Transient catalog-ID context; reviewed drawing plans, bound local appearance intents and qualified tool guidance.'}},
                'reports':{'list':'/api/drawing/science/reports','read':'/api/drawing/science/reports/{report_id}'},
                'scope':'Local catalog evidence and bounded calculations. Read each report method, units, provenance and warnings.'}

    @router.post('/science/query')
    def query_information(request: KnowledgeQuery):
        try:
            from .cofdraw_knowledge import query_catalog
            return query_catalog(root, **request.model_dump())
        except (ValueError,KeyError,TypeError,OSError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.get('/science/reports')
    def recent_reports(limit: int=20):
        if not 1<=limit<=100:raise HTTPException(422,'Choose a report limit between 1 and 100')
        reports=[]
        for directory in sorted(outputs.iterdir(),key=lambda p:p.stat().st_mtime,reverse=True)[:500] if outputs.exists() else []:
            data=stored_report(directory.name)
            if data:
                reports.append({key:data.get(key) for key in ['report_id','operation','kind','created_at_utc','engine_version','input_sha256','structure_sha256','downloads']})
                if len(reports)>=limit:break
        return {'reports':reports,'limit':limit,'scope':'Recent scientific reports in this local data directory; at most 500 output folders inspected.'}

    def stored_report(report_id):
        if not re.fullmatch(r'[0-9a-f]{32}',report_id):return None
        path=(outputs/report_id/f'{report_id}.json').resolve()
        try:
            if not path.is_relative_to(outputs.resolve()) or not path.is_file() or path.stat().st_size>40_000_000:return None
            data=json.loads(path.read_text(encoding='utf-8'))
            return data if isinstance(data,dict) and data.get('report_id')==report_id else None
        except (OSError,ValueError,UnicodeError):return None

    @router.get('/science/reports/{report_id}')
    def read_report(report_id: str):
        data=stored_report(report_id)
        if data is None:raise HTTPException(404,'Scientific report not found or unreadable')
        return data

    @router.post('/science/reports/{report_id}/circuits/{count}')
    def crystal_circuits(report_id: str, count: int, view: CrystalDrawingView = CrystalDrawingView()):
        from .cofdraw_crystal_figures import adjacent_circuits,chemical_drawing
        report=stored_report(report_id)
        if report is None:raise HTTPException(404,'Structure report was not found')
        ring=report.get('crystal_details',{}).get('ring')
        patch=report.get('reference_circuit_patches',{}).get(str(count))
        star=report.get('junction_patch') if count==1 else None
        if not ring and not patch and not star:raise HTTPException(422,'No connected circuit of this extent was traced in this crystal')
        try:
            structure=report['structure']
            edges=report['display_edges']
            provenance=structure.get('provenance',{})
            record=provenance.get('drawing_cofex',{})
            restored_order=False
            if record.get('scope')=='source' and record.get('linkage')=='imine' and not record.get('precursor_smiles'):
                # Older saved reports retain sorted COFex but lack A/B order.
                from .cofdraw_starters import PHASES
                from .framework_drawing import canonical
                phase=PHASES.get(provenance.get('template'))
                if phase:
                    by_id={r['id']:r['smiles'] for r in catalog_monomers(read_catalog(root))}
                    if all(k in by_id for k in phase[1]):
                        feeds=[canonical(by_id[k]) for k in phase[1]]
                        if '.'.join(sorted(feeds))==record.get('text','').partition(' ')[0]:
                            record={**record,'precursor_smiles':feeds}
                            structure={**structure,'provenance':{**provenance,'drawing_cofex':record}}
                            restored_order=True
            if report.get('cell_validation',{}).get('periodic_bonds'):
                structure={**structure,'chemical_bonds':report['cell_validation']['periodic_bonds']}
            junction=report.get('chemical_junction_svg') if star else None
            if star and (view.arm_spread or restored_order):
                junction=chemical_drawing(structure,edges,circuit=star,arm_spread=view.arm_spread)
            template=patch or (adjacent_circuits(ring,structure['cell'],count) if ring else None)
            if not template:
                if junction:return {'count':1,'junction_svg':junction,'junction_loops':[]}
                raise ValueError('These circuits do not share an assigned periodic edge')
            if view.tilt:
                import numpy as np
                normal=np.asarray(template['direction'])
                tangent=np.cross(normal,np.eye(3)[np.argmin(np.abs(normal))])
                angle=math.radians(view.tilt)
                template={**template,'direction':(normal*math.cos(angle)+tangent/np.linalg.norm(tangent)*math.sin(angle)).tolist()}
            projected=chemical_drawing(structure,edges,circuit=template)
            scheme=None
            if template.get('backbone') or patch:
                # A reference scheme has its own unobstructed chemical view.
                # Apply a user's tilt about that fitted view, not a newly
                # chosen camera after tilting the raw CIF projection.
                scheme_patch={**patch,'scheme_tilt':view.tilt} if patch and view.tilt and not patch.get('backbone') else template
                try:scheme=chemical_drawing(structure,edges,circuit=scheme_patch,scheme=True)
                except ValueError:pass
            loops=template.get('junction_loops')
            if loops is None:
                loops=[[[f+t for f,t in zip(structure['atoms'][i]['fractional'],shift)] for i,*shift in r['backbone']]
                       for r in template.get('rings',[template])]
            return {'count':count,'sites':len(template['site_images']),'scheme_svg':scheme,'projection_svg':projected,
                    'junction_svg':junction,'junction_loops':loops,'direction':template['direction']}
        except (ValueError, KeyError, TypeError, RuntimeError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.post('/science/reports/{report_id}/view')
    def repeat_view(report_id: str, request: DisplayCells):
        from .cofdraw_science import periodic_view
        report = stored_report(report_id)
        if report is None:
            raise HTTPException(404, 'Structure report was not found')
        if not report.get('structure', {}).get('periodic') or not report.get('viewer', {}).get('periodic'):
            raise HTTPException(422, 'Display cells require a periodic structure')
        edges = report.get('display_edges', report.get('cell_validation', {}).get('periodic_bonds'))
        try:
            viewer = periodic_view(report['structure'], request.repeats, periodic_bonds=edges,
                                   color_networks=bool(report.get('color_periodic_networks')))
            return {'viewer': viewer, 'repeats': request.repeats}
        except (ValueError, KeyError, TypeError, RuntimeError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post('/science/inspect')
    def inspect_structure(request: StructureInspection):
        try:
            if len(json.dumps(request.structure))>2_000_000:raise ValueError('Structure input limit is 2 MB')
            from .cofdraw_characterization import structure_summary
            from .cofdraw_network import inspect_periodic_network
            with lock:
                result={'kind':'periodic_structure_inspection','summary':structure_summary(request.structure),
                        'network':inspect_periodic_network(request.structure,bond_scale=request.bond_scale,tolerance_angstrom=request.tolerance_angstrom)}
            return scientific_result(result,'inspect',request.model_dump())
        except (ValueError,KeyError,TypeError,RuntimeError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.post('/science/finite-inspect')
    def inspect_current_finite(request: FiniteInspection):
        try:
            if len(json.dumps(request.model))>2_000_000:raise ValueError('Finite model input limit is 2 MB')
            from .cofdraw_current_model import inspect_finite
            with lock:result=inspect_finite(request.model)
            return scientific_result(result,'finite-inspect',request.model_dump())
        except (ValueError,KeyError,TypeError,RuntimeError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.post('/science/finite-powder')
    def finite_diffraction(request: FinitePowder):
        try:
            if len(json.dumps(request.model))>2_000_000:raise ValueError('Finite model input limit is 2 MB')
            from .cofdraw_current_model import finite_powder
            with lock:result=finite_powder(**request.model_dump())
            return scientific_result(result,'finite-powder',request.model_dump())
        except (ValueError,KeyError,TypeError,RuntimeError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.post('/science/stacking')
    def stack_periodic_model(request: PeriodicStacking):
        try:
            if len(json.dumps(request.structure))>2_000_000:raise ValueError('Structure input limit is 2 MB')
            from .cofdraw_stacking import periodic_stacking
            with lock:result=periodic_stacking(**request.model_dump())
            return scientific_result(result,'stacking',request.model_dump())
        except (ValueError,KeyError,TypeError,RuntimeError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.post('/science/pxrd')
    def diffraction(request: PowderDrawing):
        try:
            if len(json.dumps(request.structure))>2_000_000:raise ValueError('Structure input limit is 2 MB')
            from .cofdraw_characterization import powder_xrd
            with lock:result=powder_xrd(**request.model_dump())
            return scientific_result(result,'pxrd',request.model_dump())
        except (ValueError,KeyError,TypeError,RuntimeError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.post('/science/pxrd-peaks')
    def peak_comparison(request: PeakComparison):
        from .cofdraw_characterization import compare_powder_peaks
        try:
            if len(json.dumps(request.structure)) > 2_000_000:
                raise ValueError('Structure input limit is 2 MB')
            with lock:
                result = compare_powder_peaks(**request.model_dump())
            return scientific_result(result, 'pxrd-peaks', request.model_dump())
        except (ValueError, KeyError, TypeError, RuntimeError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post('/science/reagent-balance')
    def weigh_out(request: ReagentBalance):
        from .cofdraw_bench import reagent_balance
        try:
            return scientific_result(reagent_balance(**request.model_dump()), 'reagent-balance', request.model_dump())
        except (ValueError, KeyError, TypeError, RuntimeError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post('/science/screen')
    def screen(request: ScreenBatch):
        try:
            from .cofdraw_screening import screen_structures
            with lock:result=screen_structures(**request.model_dump())
            return scientific_result(result,'screen',request.model_dump())
        except (ValueError,KeyError,TypeError,RuntimeError) as exc:raise HTTPException(422,str(exc)) from exc

    return router
