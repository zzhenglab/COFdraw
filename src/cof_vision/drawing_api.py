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
from fastapi.responses import FileResponse
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


class LayoutDrawing(RenderDrawing):
    operation: Literal['clean'] = 'clean'


class ModelDrawing(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False)
    template_id: str | None = Field(default=None, max_length=80)
    cif_text: str | None = Field(default=None, max_length=2_000_000)
    repeats: list[Annotated[int, Field(ge=1, le=4)]] = Field(default=[1,1,1], min_length=3, max_length=3)
    recipe_match: RenderDrawing | None = None


class Patch3D(RenderDrawing):
    layers: int = Field(default=3, ge=1, le=6)
    layer_spacing: float = Field(default=3.4, ge=2.5, le=12, allow_inf_nan=False)
    stacking: Literal['AA','AB','ABC','slipped','twisted'] = 'AA'
    slip: list[Annotated[float, Field(ge=-20,le=20,allow_inf_nan=False)]] | None = Field(default=None,min_length=2,max_length=2)
    twist_degrees: float = Field(default=15, ge=-180, le=180, allow_inf_nan=False)


class PeriodicFromDrawing(RenderDrawing):
    layer_spacing: float | None = Field(default=None, ge=2.5, le=12, allow_inf_nan=False)
    repeats: list[Annotated[int, Field(ge=1, le=4)]] = Field(default=[1,1,1], min_length=3, max_length=3)


class NetworkHypothesis(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False)
    structure: dict
    offsets: list[list[Annotated[float, Field(ge=-1, le=1, allow_inf_nan=False)]]] = Field(min_length=1, max_length=4)
    repeats: list[Annotated[int, Field(ge=1, le=4)]] = Field(default=[1,1,1], min_length=3, max_length=3)


class MonomerAnalysis(BaseModel):
    model_config = ConfigDict(extra='forbid')
    smiles: str = Field(min_length=1, max_length=4000)
    conformers: int = Field(default=4, ge=1, le=8)


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
    base_url: str = Field(default='https://api.openai.com/v1', max_length=1000)
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
    from .cofdraw_adsorption import create_adsorption_router
    router.include_router(create_adsorption_router(outputs))

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
        monomers = decorated_monomers(catalog.get('monomer_presets', MONOMERS))
        monomers = [{**m, 'synthetic_accessibility': synthetic_accessibility(m.get('smiles', ''))} for m in monomers]
        return {"schema_version": VERSION, "linkages": library_linkages(catalog), "source_recipes": sources,
                "templates": TEMPLATES, "topology_catalog": catalog.get("topology_catalog", []),
                "style_controls": {"palettes":list(PALETTES),"palette_colors":PALETTES,"fonts":FONTS,"edge_styles":["wavy_cross","dashed","crossbar","bond_stub"],
                                   "abbreviation_modes":["explicit","condensed","OMe"],"layer_count":[2,5],"layer_opacity":[.05,.7],
                                   "fill_opacity":[.05,.6],"font_scale":[.6,1.6],"bond_scale":[.5,2.5]},
                "styles": STYLES, "monomer_presets": monomers,
                "cof_presets": catalog.get("cof_presets", PRESETS), "catalog_version": catalog.get("catalog_version", "builtin-1"),
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
                proposal = select_with_llm(text, describe_catalog(catalog), base_url=request.llm.base_url,
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
                    base_url=request.llm.base_url, model=request.llm.model, api_key=secret)
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
            from .cofdraw_synthesis import precursor_scores
            audit['precursor_synthetic_accessibility'] = precursor_scores(recipe)
            audit.update(design_cofex_notation=recipe.get("design_cofex_notation"), design_topology=recipe.get("design_topology"),
                         design_identity_qualification=recipe.get("design_identity_qualification"), selected_template=recipe.get("selected_template"))
            audit.update(cofex_string=metadata.get("cofex_notation"), cofex_notation=metadata.get("cofex_notation"),
                         literature_metadata=metadata, monomers=[{"component": name, "label": c.get("label", name), "smiles": c.get("monomer_smiles")} for name, c in recipe.get("components",{}).items()],
                         chemical_export_scope="Editable finite patch with continuation atoms. SVG styles are separate presentation; no complete periodic unit cell is implied.")
            urls = {ext: f"/api/drawing/exports/{drawing_id}/{ext}" for ext in ("svg", "png", "mol", "json", *extra_formats)}
            bundle = {"recipe": recipe, "audit": audit, "downloads": urls}
            (directory / f"{drawing_id}.json").write_text(json.dumps(bundle, indent=2) + "\n", encoding="utf-8")
            return {"drawing_id": drawing_id, "status": status, "audit": audit, "downloads": urls}
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
        return json.loads(path.read_text(encoding='utf-8')) if path.is_file() else {'templates':[]}

    def load_crystal_template(item, manifest):
        from .cofdraw_science import import_cif
        directory = (root/'crystal-models').resolve()
        path = (directory/item['path']).resolve()
        if not path.is_relative_to(directory) or hashlib.sha256(path.read_bytes()).hexdigest()!=item['sha256']:
            raise ValueError('Crystal template integrity check failed')
        structure = import_cif(path.read_text(encoding='utf-8'))
        structure['provenance'].update(template=item['id'], generator=manifest['generator'],
            model_name=item['model_name'], doi=manifest['citation_doi'], qualification=item['description'])
        return structure

    @router.get('/science/library')
    def scientific_library():
        data=crystal_manifest()
        return {'templates':data['templates'],'methods':['RDKit ETKDG / MMFF or UFF','Periodic probe geometry','BET fit from supplied isotherm']}

    @router.get('/science/synthesis-library')
    def synthesis_library():
        from .cofdraw_synthesis import catalog_scores
        return catalog_scores(read_catalog(root).get('monomer_presets', MONOMERS))

    @router.post('/science/periodic-options')
    def matched_periodic_options(request: RenderDrawing):
        try:
            from .cofdraw_periodic_options import periodic_options
            manifest=crystal_manifest()
            with lock:return periodic_options(selected_recipe(request),manifest)
        except (ValueError,KeyError,TypeError,RuntimeError,OSError) as exc:raise HTTPException(422,str(exc)) from exc

    @router.post('/science/model')
    def model(request: ModelDrawing):
        try:
            from .cofdraw_science import import_cif, structure_exports
            if bool(request.template_id)==bool(request.cif_text):raise ValueError('Choose one model template or import one CIF')
            if request.template_id:
                manifest=crystal_manifest()
                matches=[x for x in manifest['templates'] if x['id']==request.template_id]
                if len(matches)!=1:raise ValueError('Unknown crystal template')
                item=matches[0];structure=load_crystal_template(item,manifest)
                if request.recipe_match:
                    from .cofdraw_periodic_options import periodic_options
                    options=periodic_options(selected_recipe(request.recipe_match),manifest)
                    if not any(o['id']==item['id'] for o in options['options']):raise ValueError('This periodic model does not match the current precursor pair, linkage and unchanged connectivity template')
                    structure['provenance']['drawing_match']={**item['precursor_match'],'scope':options['scope']}
            else:
                if request.recipe_match:raise ValueError('Recipe matching applies only to validated built-in periodic designs')
                structure=import_cif(request.cif_text)
            with lock:result=structure_exports(structure,repeats=tuple(request.repeats))
            result['scope']='Periodic atomic model with an explicit cell. Built-in templates are unrelaxed reticular designs, not experimentally assigned crystals; imported CIFs retain supplied geometry.'
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
                try:
                    design = validated_design(recipe)
                    options = periodic_options(recipe, manifest, design=design)
                except (ValueError, KeyError, TypeError) as exc:
                    return {'available': False, 'reason': str(exc)}
                matches = {o['id'] for o in options['options']}
                item = next((entry for entry in manifest['templates']
                             if entry['id'] in matches and entry['model_name'].endswith('-AA')), None)
                if item and request.layer_spacing is None:
                    structure = load_crystal_template(item, manifest)
                    structure['provenance']['drawing_match'] = item['precursor_match']
                    structure['provenance']['spacing_angstrom'] = structure['cell'][2][2]
                    validation = {'matched_template': item['id'], 'template_cif_sha256': item['sha256']}
                else:
                    try:
                        structure, validation = build_periodic_layer(design[0], spacing=request.layer_spacing)
                    except ValueError as exc:
                        return {'available': False, 'reason': str(exc)}
                result = structure_exports(structure, repeats=tuple(request.repeats))
            result.update(available=True, cell_validation=validation,
                scope='Unrelaxed periodic starting model of the selected precursor chemistry and planar connectivity. AA layer spacing is a model parameter.')
            provenance = structure['provenance']
            result['exports']['cif'] = result['exports']['cif'].replace('data_COFdraw\n', (
                'data_COFdraw\n# COFdraw unrelaxed AA starting model.\n'
                f"# Template {design[0]['selected_template']['id']}; layer spacing {provenance['spacing_angstrom']:.4f} angstrom.\n"
                + (f"# Meso phenyl tilt {provenance['phenyl_tilt_degrees']} degrees.\n"
                   if provenance.get('phenyl_tilt_degrees') else '')), 1)
            return scientific_result(result, 'build-periodic', request.model_dump())
        except (ValueError, KeyError, TypeError, RuntimeError, OSError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post('/science/patch')
    def patch_model(request: Patch3D):
        try:
            from .cofdraw_science import patch_3d
            with lock:result=patch_3d(selected_recipe(request),layer_spacing=request.layer_spacing,layers=request.layers,stacking=request.stacking,slip=request.slip,twist_degrees=request.twist_degrees)
            result['scope']='Finite illustrative layers from scaled drawing coordinates; not a relaxed periodic crystal. X sites mark continuation boundaries.'
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
            return scientific_result(result, 'geometry', request.model_dump())
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
        models={'model':ModelDrawing,'patch':Patch3D,'networks':NetworkHypothesis,'monomer':MonomerAnalysis,
                'porosity':PoreAnalysis,'geometry':GeometryAnalysis,'bet':BETAnalysis,'inspect':StructureInspection,
                'pxrd':PowderDrawing,'screen':ScreenBatch,'query':KnowledgeQuery,
                'finite-inspect':FiniteInspection,'finite-powder':FinitePowder,'stacking':PeriodicStacking,'periodic-options':RenderDrawing,
                'measure':RenderDrawing,'build-periodic':PeriodicFromDrawing}
        return {'engine_version':__version__, 'schema_version':'cofdraw-capabilities-1',
                'operations':{name:{'method':'POST','path':f'/api/drawing/science/{name}',
                                    'request_schema':model.model_json_schema()} for name,model in models.items()},
                'workflows':capability_inventory(),
                'synthesis_library':{'method':'GET','path':'/api/drawing/science/synthesis-library'},
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

    @router.post('/science/screen')
    def screen(request: ScreenBatch):
        try:
            from .cofdraw_screening import screen_structures
            with lock:result=screen_structures(**request.model_dump())
            return scientific_result(result,'screen',request.model_dump())
        except (ValueError,KeyError,TypeError,RuntimeError) as exc:raise HTTPException(422,str(exc)) from exc

    return router
