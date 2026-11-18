"""Bounded, deterministic retrieval from local COFdraw records and capabilities.

This is catalog retrieval, not a language model, literature search, property
predictor or chemical substructure search. Returned catalog payloads are exact.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re
import unicodedata


RECORD_TYPES = ("monomer", "linkage", "topology", "source", "crystal", "template", "style", "preset", "capability", "literature")
_STOP_WORDS = {"find", "show", "me", "please", "what", "which", "the", "a", "an", "about", "for", "with", "of", "and", "is", "are", "in", "do", "does", "can", "available"}
_PLURALS = {k + "s": k for k in RECORD_TYPES}
_PLURALS.update({"capabilities": "capability", "topologies": "topology"})


def normalize_query(value):
    value = str(value).casefold().replace("β", "beta").replace("α", "alpha").replace("γ", "gamma")
    value = re.sub(r"https?://(?:dx\.)?doi\.org/", "", value)
    value = re.sub(r"\bdoi\s*:\s*", "", value)
    value = unicodedata.normalize("NFKD", value)
    value = "".join(c for c in value if not unicodedata.combining(c))
    return " ".join(re.findall(r"[a-z0-9]+", value))


def _hash(raw):
    return hashlib.sha256(raw).hexdigest()


def _read_json(path, max_bytes=4_000_000):
    if not path.is_file() or path.stat().st_size > max_bytes:
        raise ValueError("Missing or oversized local catalog file")
    raw = path.read_bytes()
    return json.loads(raw.decode("utf-8")), _hash(raw)


def _qualifications(row):
    result = []
    for key in ("qualification", "qualifications", "scope", "limitations", "runtime_scope", "representation_scope",
                "identity_qualification", "compatibility_scope", "name_scope", "template_scope", "prepare_scope",
                "node_symmetry_scope", "reactive_site_count_scope", "limits", "description"):
        values = row.get(key, [])
        if isinstance(values, str):
            values = [values]
        if isinstance(values, list):
            for value in values:
                if isinstance(value, str) and value not in result:
                    result.append(value)
    return result


def _record(kind, identifier, label, data, support, provenance, qualifications=()):
    return {"id": f"{kind}:{identifier}", "record_type": kind, "record_id": identifier,
            "label": label, "support": support, "provenance": provenance,
            "qualifications": list(qualifications), "data": deepcopy(data)}


def _bound_file(base, relative, expected_hash, max_bytes):
    if not isinstance(relative, str) or not isinstance(expected_hash, str):
        raise ValueError("Missing file binding")
    path = (base / relative).resolve()
    if not path.is_relative_to(base.resolve()) or not path.is_file() or path.stat().st_size > max_bytes:
        raise ValueError("File binding escapes the local bundle or exceeds its size limit")
    raw = path.read_bytes()
    if _hash(raw) != expected_hash:
        raise ValueError("Bundled file hash mismatch")
    return raw


def capability_inventory():
    """Describe implemented API operations, taking request schemas from models."""
    from . import drawing_api
    from .cofdraw_adsorption import AdsorptionRequest
    definitions = [
        ("cof-conversation", "/api/drawing/describe/chat", "DescribeConversation", "Discuss and revise a catalog-grounded COF drawing plan",
         "Local bounded follow-up intent and catalog IDs; explicit optional provider settings. Review before apply; never infer a periodic crystal or experimental identity from prose.",
         {"text":"Use the blue and coral palette", "mode":"rules"}, ["json"]),
        ("matching-periodic-design", "/api/drawing/science/periodic-options", "RenderDrawing", "Find periodic designs matching the current precursors and connectivity",
         "Exact canonical precursors, linkage, product fragments and unchanged chosen connectivity template. Separate unrelaxed coordinates, never a conversion of a finite drawing into an assigned crystal.", None, ["json"]),
        ("build-periodic-cell", "/api/drawing/science/build-periodic", "PeriodicFromDrawing", "Build and save a periodic CIF for the current drawing",
         "Uses a matched AA template or a fitted hexagonal or square cell with checked periodic bonds. Returns unrelaxed starting coordinates and validation or an explicit unavailable reason.", None, ["cif","xyz","pdb","json"]),
        ("finite-inspection", "/api/drawing/science/finite-inspect", "FiniteInspection", "Check atoms and contacts in the displayed finite structure without a CIF",
         "Displayed real atoms and supplied bonds only; continuation sites omitted, missing hydrogen not invented. No density, PLD or periodic dimensionality.", None, ["json"]),
        ("finite-scattering", "/api/drawing/science/finite-powder", "FinitePowder", "Debye powder scattering from displayed finite atoms",
         "All-pairs finite-cluster scattering with neutral-atom form factors, bounded atom/pair-point work. No periodic hkl or Bragg assignment; illustrative drawing coordinates are not relaxed crystals.", None, ["csv","json"]),
        ("periodic-stacking", "/api/drawing/science/stacking", "PeriodicStacking", "Compare AA, AB, ABC and slipped periodic layer hypotheses",
         "One screened periodic layer with contact and network guards. Explicit corrugated-layer mode checks the output bond graph. General 3D networks, multilayer input and arbitrary twist are rejected.", None, ["cif","xyz","pdb","json"]),
        ("describe-cof", "/api/drawing/describe", "DescribeDrawing", "Describe a COF and review a catalog-grounded drawing plan",
         "Offline name and word matching by default; optional explicitly configured LLM selects catalog IDs. Chemistry is validated locally. Ambiguity and unsupported requests require review; prose does not determine a unique crystal or establish a source identity.",
         {"text": "Draw TFB and TAPB with an imine linkage in pastel colors", "mode": "rules"}, ["json"]),
        ("local-query", "/api/drawing/science/query", "KnowledgeQuery", "Query local COF records and calculation capabilities",
         "Deterministic local keyword/name/DOI retrieval; no web search, property inference or fabricated answers.",
         {"query": "Pyurea", "record_types": ["source"], "limit": 10}, ["json"]),
        ("draw-source-or-design", "/api/drawing/render", "RenderDrawing", "Render source or designed finite chemical patches",
         "PNG/SVG and editable finite chemistry. Source qualifiers and saved coordinates remain distinct from user edits.",
         {"source_recipe_id": "pyurea-source-layout-v1", "style_id": "source", "size": 1200}, ["png", "svg", "mol", "sdf", "cdxml", "json"]),
        ("prepare-design", "/api/drawing/prepare", "PrepareDrawing", "Prepare a finite design from two monomers",
         "Only installed linkage transforms and compatible explicit templates; no synthesis feasibility or experimental topology prediction.",
         {"monomer_a": "O=Cc1cc(C=O)cc(C=O)c1", "monomer_b": "Nc1ccc(N)cc1", "linkage": "imine", "template_id": "hcb_3_2", "instances": 12}, ["json"]),
        ("monomer-conformers", "/api/drawing/science/monomer", "MonomerAnalysis", "Molecular descriptors and force-field conformers",
         "Complete finite molecules; seeded ETKDG and MMFF/UFF. Empirical charges and force-field energies are not DFT, band gaps or crystal stability.",
         {"smiles": "Nc1ccc(N)cc1", "conformers": 2}, ["xyz", "sdf", "json"]),
        ("monomer-solubility", "/api/drawing/science/solubility", "SolubilityScreen", "Screen neat-solvent solubility near room temperature",
         "Exact equilibrium records, bounded analog estimates and separate reported preparation loads. Mass is relative to starting solvent volume.",
         {"monomer_id": "ta", "solvent_ids": ["dmf", "dioxane", "ethanol"], "solvent_ml": 2}, ["csv", "json"]),
        ("periodic-model", "/api/drawing/science/model", "ModelDrawing", "Import periodic CIF or select an unrelaxed model",
         "Explicit cell and coordinates; supercell export does not optimize or experimentally validate the structure.",
         {"template_id": "hcb-imine-aa", "repeats": [2, 2, 1]}, ["cif", "xyz", "pdb", "json"]),
        ("finite-stack", "/api/drawing/science/patch", "Patch3D", "Illustrative layers from a finite drawing",
         "Nonperiodic drawing-derived coordinates, not a relaxed crystal; continuation sites remain X. AB requires explicit slip.",
         {"source_recipe_id": "ucof-source-layout-v1", "layers": 2, "layer_spacing": 3.4, "stacking": "AB", "slip": [1.7, 0]}, ["xyz", "json"]),
        ("network-hypothesis", "/api/drawing/science/networks", "NetworkHypothesis", "Compare explicitly translated periodic-network copies",
         "User-chosen fractional offsets with collision checks. Copy count is not an assigned interpenetration degree, cat value or stability result.",
         None, ["cif", "xyz", "pdb", "json"]),
        ("periodic-porosity", "/api/drawing/science/porosity", "PoreAnalysis", "Periodic hard-sphere pore geometry",
         "Requires an explicit periodic structure. Sampled volumes/surfaces and voxel connectivity are not BET area or certified transport accessibility.",
         None, ["json"]),
        ("periodic-geometry", "/api/drawing/science/geometry", "GeometryAnalysis", "Deterministic periodic geometry with pore-diameter bounds",
         "Explicit periodic hard-sphere model: grid lower/upper bounds for PLD/LCD and transport-connected probe-center volume; sampled probe-inflated surface connected by certified segments. Coarse grids can leave transport unresolved. Not Zeo++/Voronoi, adsorption, BET, helium porosimetry or experimental pore size.",
         None, ["json"]),
        ("local-adsorption", "/api/drawing/science/adsorption/jobs", "AdsorptionRequest", "Start a local RASPA3 adsorption simulation job",
         "Asynchronous rigid-framework classical simulation: separate Widom Henry coefficients, GCMC loadings and mixture selectivity, optional finite-difference zero-loading heat. Explicit model/charge/force-field/fugacity assumptions; defaults do not establish convergence. Poll the returned job_id until terminal; cancellation preserves completed stages. No finite-loading heat curve, chemisorption or synthesis prediction.",
         None, ["json", "csv", "zip"]),
        ("bet-fit", "/api/drawing/science/bet", "BETAnalysis", "Fit supplied adsorption-branch data",
         "Requires user-supplied measured data and a selected interval; no BET value is inferred from a drawing or geometric surface.",
         None, ["json"]),
        ("drawing-opening", "/api/drawing/science/measure", "RenderDrawing", "Estimate an opening in a 2D drawing",
         "Conventional bond scaling and sampled projected circle; not experimental pore size or pore-limiting diameter.",
         {"source_recipe_id": "pyurea-source-layout-v1"}, ["json"]),
        ("layout-revision", "/api/drawing/layout", "LayoutDrawing", "Explicit rigid-block coordinate cleanup",
         "Assembled component recipes only; direct explicit source graphs retain supplied coordinates. An edited layout loses active source registration.",
         {"source_recipe_id": "tppa-pore-v1", "operation": "clean"}, ["json"]),
        ("structure-inspection", "/api/drawing/science/inspect", "StructureInspection", "Cell, composition and periodic contact inspection",
         "Explicit periodic cell: formula, density, element-site densities and distance-inferred network contacts. Cutoff-dependent dimensionality is not bond order, formal valence, RCSR net assignment or active-site accessibility.",
         None, ["json"]),
        ("powder-diffraction", "/api/drawing/science/pxrd", "PowderDrawing", "Ideal powder X-ray diffraction PXRD pattern",
         "Explicit periodic structure and wavelength; neutral-atom form factors, reciprocal metric, Lorentz/polarization correction and chosen Gaussian broadening. No experimental phase identification, refinement, preferred orientation or instrumental-profile prediction.",
         None, ["csv", "json"]),
        ("structure-screen", "/api/drawing/science/screen", "ScreenBatch", "Compare supplied periodic structures in a bounded batch",
         "Up to eight named structures, with optional geometric porosity. Screening compares supplied models; no adsorption, synthesis feasibility, stability or experimental identity is predicted.",
         None, ["csv", "json"]),
    ]
    result = []
    for identifier, endpoint, model_name, label, scope, example, formats in definitions:
        schema = (AdsorptionRequest if model_name == "AdsorptionRequest" else getattr(drawing_api, model_name)).model_json_schema()
        result.append(_record("capability", identifier, label,
            {"method": "POST", "endpoint": endpoint, "request_schema": schema,
             "example_request": example, "possible_formats": formats, "scope": scope,
             "formats_qualification": "Only files actually listed in a result's downloads are available."},
            {"status": "implemented", "scope": scope}, {"kind": "local_runtime_contract", "request_model": model_name}, [scope]))
    for identifier, endpoint, label, scope in [
        ("adsorption-engine", "/api/drawing/science/adsorption/engine", "Inspect local adsorption engine availability",
         "Reports the pinned local RASPA3 CPU package/profile and supported gases without starting a simulation; availability does not establish numerical convergence or experimental validity."),
        ("precursor-sa-library", "/api/drawing/science/synthesis-library", "Offline finite-monomer SA Score catalog",
         "RDKit Contrib SA Score for eligible complete finite catalog monomers: 1 easier to 10 harder. A molecular heuristic, not ASKCOS, retrosynthesis, reaction support or COF yield/feasibility. Render audits separately summarize only declared finite precursors."),
        ("api-contracts", "/api/drawing/science/capabilities", "Read actual JSON request schemas",
         "Read the running service's supported operations and exact request schemas before constructing calculations."),
        ("saved-reports", "/api/drawing/science/reports", "List recent local machine-readable calculation reports",
         "Reports belong to this runtime data directory. Listing historical results does not rerun or independently validate them."),
        ("saved-report", "/api/drawing/science/reports/{report_id}", "Retrieve one saved calculation report",
         "Use a report_id returned by this service; retain operation, engine version, input hash, parameters, units and qualifications.")]:
        result.append(_record("capability", identifier, label, {"method": "GET", "endpoint": endpoint,
            "request_schema": None, "scope": scope}, {"status": "implemented", "scope": scope},
            {"kind": "local_runtime_contract"}, [scope]))
    return result


def _catalog_records(data_root):
    from .cofdraw_chemistry import TEMPLATES, library_linkages, inspect_monomer
    from .cofdraw_core import catalog_monomers
    from .cofdraw_styles import STYLES

    root = Path(data_root).resolve()
    errors, records, catalog = [], [], {}
    catalog_info = {"version": None, "sha256": None, "registry_sha256": None, "crystal_manifest_sha256": None}
    try:
        catalog, catalog_sha = _read_json(root / "framework-drawing/catalog/current.json")
        if not isinstance(catalog, dict):
            raise ValueError("Drawing catalog must be an object")
        catalog_info.update(version=catalog.get("catalog_version"), sha256=catalog_sha)
    except (ValueError, OSError) as exc:
        catalog = {}
        errors.append({"record_type": "catalog", "reason": str(exc)})
    catalog_provenance = {"kind": "bundled_catalog", "path": "framework-drawing/catalog/current.json", "sha256": catalog_info["sha256"]}
    coverage_binding=catalog.get('literature_coverage')
    if coverage_binding:
        try:
            coverage=json.loads(_bound_file(root/'framework-drawing',coverage_binding['path'],coverage_binding['sha256'],4_000_000).decode('utf-8'))
            for paper in coverage['papers']:
                records.append(_record('literature',paper['review_id'],paper['title'],paper,
                    {'status':paper['review_status'],'ground_truth':False,'full_main_and_si_audit':False},
                    {**catalog_provenance,'coverage_sha256':coverage_binding['sha256'],'source_dois':[paper['doi']] if paper.get('doi') else []},
                    [coverage['scope'],paper.get('review_scope',''),*coverage['limitations']]))
        except (ValueError,OSError,KeyError) as exc:
            errors.append({'record_type':'literature','reason':str(exc)})
    for row in catalog_monomers(catalog) if catalog else []:
        runtime=inspect_monomer(row['smiles'])
        records.append(_record("monomer", row["id"], row.get("full_name", row.get("label", row["id"])), row,
            {"status": "catalogued_monomer", "automatic_prepare_supported": row.get('automatic_prepare_supported') is not False and any(r['reactive_sites']>=2 for r in runtime['compatible_roles']),
             "runtime_compatible_roles":runtime['compatible_roles'],"runtime_reactive_groups":runtime['reactive_groups'],
             "scope": "Library presence does not establish a supported reaction, periodic net or synthesis feasibility."},
            {**catalog_provenance, "source_dois": row.get("source_dois", []), "source_evidence": row.get("source_evidence", [])}, _qualifications(row)))
    for identifier, runtime in library_linkages(catalog).items():
        row = next((r for r in catalog.get("linkage_catalog", []) if r["id"] == identifier), None)
        records.append(_record("linkage", identifier, runtime.get("label", identifier), {"catalog": row, "runtime": runtime},
            {"status": "bounded_automatic_transform" if runtime["supported"] else "reference_only",
             "automatic_prepare_supported": bool(runtime["supported"]), "scope": runtime.get("limits")},
            {**catalog_provenance, "runtime_authority": "cofdraw_chemistry.library_linkages"}, _qualifications(runtime)))
    for row in catalog.get("topology_catalog", []):
        actual_templates = [t["id"] for t in TEMPLATES if t["id"] in row.get("compatible_template_ids", [])]
        records.append(_record("topology", row["id"], row.get("label", row["id"]), row,
            {"status": "finite_template_available" if actual_templates else "reference_only",
             "template_ids": actual_templates, "periodic_net_assignment": False}, catalog_provenance, _qualifications(row)))
    for row in catalog.get("cof_presets", []):
        records.append(_record("preset", row["id"], row.get("label", row["id"]), row,
            {"status": "design_starting_inputs", "scope": "Preparing a preset creates a design; load its source ID separately to preserve source incidence."},
            catalog_provenance, _qualifications(row)))
    for kind, rows, scope, authority in [
        ("template", TEMPLATES, "Finite user-chosen connectivity; no experimental net inference.", "cofdraw_chemistry.TEMPLATES"),
        ("style", STYLES, "Presentation only; no chemical graph or source identity inference.", "cofdraw_styles.STYLES")]:
        for row in rows:
            records.append(_record(kind, row["id"], row.get("label", row["id"]), row,
                {"status": "implemented", "scope": scope}, {"kind": "local_runtime_registry", "authority": authority}, _qualifications(row) + [scope]))
    try:
        registry, registry_sha = _read_json(root / "framework-drawing/registry.json")
        catalog_info["registry_sha256"] = registry_sha
        if catalog.get("source_registry_sha256") not in (None, registry_sha):
            raise ValueError("Catalog registry hash mismatch; source records excluded")
        for binding in registry.get("recipes", []):
            identifier = binding.get("recipe_id", "unknown")
            try:
                raw = _bound_file(root / "framework-drawing", binding["path"], binding["sha256"], 1_500_000)
                recipe = json.loads(raw.decode("utf-8"))
                metadata = catalog.get("source_metadata", {}).get(identifier, {})
                if (recipe.get("recipe_id") != identifier or metadata.get("source_recipe_sha256") != binding["sha256"]
                        or recipe.get("cofex_id") != metadata.get("cofex_id")):
                    raise ValueError("Source identity/metadata binding mismatch")
                records.append(_record("source", identifier, metadata.get("display_name", metadata.get("literature_name", binding.get("label", identifier))),
                    {"registry": binding, "metadata": metadata},
                    {"status": "source_bound_provisional", "render_supported": True,
                     "source_registered": bool(metadata.get("source_registered", False)), "ground_truth": False},
                    {**catalog_provenance, "registry_sha256": registry_sha, "recipe_sha256": binding["sha256"],
                     "source_dois": metadata.get("source_dois", []), "review_evidence": metadata.get("review_evidence", [])},
                    _qualifications(metadata)))
                records[-1]["retrieval"] = {"method": "GET", "path": "/api/drawing/sources/" + identifier}
            except (ValueError, OSError, KeyError) as exc:
                errors.append({"record_type": "source", "record_id": identifier, "reason": str(exc)})
    except (ValueError, OSError, KeyError) as exc:
        errors.append({"record_type": "source", "reason": str(exc)})
    try:
        manifest, manifest_sha = _read_json(root / "crystal-models/manifest.json")
        catalog_info["crystal_manifest_sha256"] = manifest_sha
        for row in manifest.get("templates", []):
            try:
                _bound_file(root / "crystal-models", row["path"], row["sha256"], 2_000_000)
                records.append(_record("crystal", row["id"], row.get("label", row["id"]), row,
                    {"status": "unrelaxed_periodic_design", "model_import_supported": True, "experimentally_validated": False},
                    {"kind": "bundled_model", "manifest_sha256": manifest_sha, "model_sha256": row["sha256"],
                     "generator": manifest.get("generator"), "citation_doi": manifest.get("citation_doi")},
                    _qualifications(row) + [manifest.get("qualification", "Unrelaxed design; no experimental identity asserted.")] ))
            except (ValueError, OSError, KeyError) as exc:
                errors.append({"record_type": "crystal", "record_id": row.get("id"), "reason": str(exc)})
    except (ValueError, OSError, KeyError) as exc:
        errors.append({"record_type": "crystal", "reason": str(exc)})
    records.extend(capability_inventory())
    return records, catalog_info, errors


def _search_fields(value, path="", depth=0):
    if depth > 12:
        return []
    if isinstance(value, str):
        return [(path, normalize_query(value))]
    if isinstance(value, dict):
        return [(key, text) for name, child in value.items() if name not in {"request_schema", "sha256", "recipe_sha256", "registry_sha256", "manifest_sha256", "model_sha256"}
                for key, text in _search_fields(child, f"{path}.{name}".strip("."), depth + 1)]
    if isinstance(value, list):
        return [(key, text) for i, child in enumerate(value) for key, text in _search_fields(child, f"{path}[{i}]", depth + 1)]
    return []


def query_catalog(data_root, query="", record_types=None, limit=20, offset=0):
    """Return exact matching local records, deterministic ranking and pagination."""
    if not isinstance(query, str) or len(query) > 500:
        raise ValueError("Query must be a string of at most 500 characters")
    if type(limit) is not int or not 1 <= limit <= 50 or type(offset) is not int or not 0 <= offset <= 10000:
        raise ValueError("Use limit 1–50 and offset 0–10000")
    if record_types is not None and (not isinstance(record_types, list) or len(record_types) > 10
            or any(kind not in RECORD_TYPES for kind in record_types)):
        raise ValueError("Unknown record type; use the documented record_types list")
    allowed = set(RECORD_TYPES if record_types is None else record_types)
    records, catalog, errors = _catalog_records(data_root)
    normalized = normalize_query(query)
    tokens = list(dict.fromkeys(_PLURALS.get(t, t) for t in normalized.split() if t not in _STOP_WORDS))
    matched = []
    for record in records:
        if record["record_type"] not in allowed:
            continue
        fields = _search_fields(record)
        # Every keyword must be present; retrieval never invents a semantic answer.
        if not all(any(token in text.split() for _, text in fields) for token in tokens):
            continue
        hit_fields = [key for key, text in fields if tokens and any(t in text.split() for t in tokens)]
        primary = [normalize_query(record[k]) for k in ("id", "record_id", "label")]
        score = 100 if normalized and normalized in primary else 40 if normalized and any(normalized in p for p in primary) else 0
        score += sum(8 for token in tokens if any(token in p.split() for p in primary))
        score += min(len(hit_fields), 20)
        record["match"] = {"score": score, "fields": hit_fields[:24]}
        matched.append(record)
    matched.sort(key=lambda r: (-r["match"]["score"], r["record_type"], r["record_id"]))
    selected = matched[offset:offset + limit]
    return {"schema_version": "cofdraw-local-knowledge-1", "query": query, "normalized_query": normalized,
        "keywords": tokens, "match_mode": "all_normalized_keywords", "available_record_types": list(RECORD_TYPES),
        "record_types": [kind for kind in RECORD_TYPES if kind in allowed], "matched_count": len(matched),
        "returned_count": len(selected), "limit": limit, "offset": offset,
        "next_offset": offset + limit if offset + limit < len(matched) else None,
        "records": selected, "catalog": catalog, "errors": errors,
        "scope": "Exact local catalog/runtime records; no external literature search, semantic inference or property prediction.",
        "no_match_meaning": "No matching bundled record is not evidence that a material, reaction or property does not exist."}
