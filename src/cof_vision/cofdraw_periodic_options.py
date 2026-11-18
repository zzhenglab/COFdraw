"""Exact precursor/template matching to separately validated periodic designs."""
from __future__ import annotations

from copy import deepcopy
from rdkit import Chem


def _graph(smiles):
    molecule=Chem.MolFromSmiles(smiles or '')
    if molecule is None or molecule.GetNumAtoms()==0:
        raise ValueError('A complete declared precursor graph is required')
    return Chem.MolToSmiles(molecule,isomericSmiles=True)


def validated_design(recipe):
    """Recover the declared chemistry and confirm its finite template incidence."""
    from .cofdraw_chemistry import design_signature, prepare_template
    from .framework_drawing import linkage_fragment

    selected=recipe.get('selected_template') or {}
    if recipe.get('explicit_graph') is not None:
        raise ValueError('This explicit atom graph needs a declared periodic connectivity or a corresponding CIF.')
    if not selected.get('id') or recipe.get('design_template_signature')!=design_signature(recipe):
        raise ValueError('This drawing has no unchanged supported connectivity template.')
    components=list(recipe['components'].values())
    if len(components)!=2:
        raise ValueError('Periodic construction requires two complete declared precursors.')
    precursors=sorted(_graph(c.get('monomer_smiles')) for c in components)
    for component in components:
        expected=linkage_fragment(component['monomer_smiles'],recipe['linkage'],component['role'])
        if _graph(expected['fragment_smiles'])!=_graph(component['fragment_smiles']):
            raise ValueError('Edited product fragments do not match the declared precursor transform')
    fresh=prepare_template(recipe['components']['A']['monomer_smiles'],recipe['components']['B']['monomer_smiles'],
                           recipe['linkage'],'A','B',6,selected['id'])

    def normalized_signature(value):
        value=deepcopy(value)
        for component in value['components'].values():
            for field in ('monomer_smiles','fragment_smiles'):
                component[field]=_graph(component[field])
        return design_signature(value)

    if normalized_signature(fresh)!=normalized_signature(recipe):
        raise ValueError('Drawing incidence differs from the selected connectivity template')
    return fresh, precursors


def periodic_options(recipe, manifest, *, design=None):
    result={'options':[], 'preserves_drawing_coordinates':False,
      'scope':'Matching precursors and chosen connectivity identify separate unrelaxed periodic design hypotheses. They do not convert the drawing coordinates into a crystal or establish experimental identity.'}
    try:
        fresh,precursors=design if design is not None else validated_design(recipe)
    except (ValueError,KeyError,TypeError) as exc:
        result['reason']=str(exc)
        return result
    selected=fresh['selected_template']['id']
    for item in manifest.get('templates',[]):
        match=item.get('precursor_match') or {}
        if (match.get('linkage')==recipe.get('linkage')
                and selected in match.get('drawing_template_ids',[])
                and sorted(_graph(s) for s in match.get('canonical_smiles',[]))==precursors):
            result['options'].append({'id':item['id'],'label':item['label'],'atoms':item['atoms'],
              'formula':item['formula'],'model_sha256':item['sha256'],
              'match':match,'qualification':item['description']})
    result['precursor_canonical_smiles']=precursors
    result['drawing_template_id']=selected
    result['reason']=('Choose a matching periodic design, then inspect its geometry before calculating pore properties.' if result['options'] else
                      'No validated periodic design for this exact precursor pair, linkage and connectivity is bundled yet. Finite checks/scattering remain available; import a corresponding CIF for periodic properties.')
    return result
