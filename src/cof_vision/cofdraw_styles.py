"""Chemical, layered and schematic depictions of one unchanged finite graph."""
from __future__ import annotations
from copy import deepcopy
import math
import random
import re
import xml.etree.ElementTree as ET
import numpy as np

PALETTES={"classic":["#d82c2c","#264bd1","#388245","#873a91"],
          "teal_orange":["#168c89","#e68b35","#955ea5","#c6536c"],
          "purple_gold":["#8155a3","#d39a24","#3e8ea4","#8ca64b"],
          "pastel":["#cf7795","#6c9cc5","#89ae74","#c8a25b"],
          "ink":["#20252b"]*4}
PALETTES.update({
    'ocean_coral': ['#167d9a', '#ee796b', '#6796b5', '#cfb86d'],
    'forest_berry': ['#2d8068', '#b8567b', '#a6b96b', '#75658d'],
    'cobalt_citrus': ['#315bb8', '#e6ac31', '#29a092', '#cc6771'],
    'rose_slate': ['#c46e90', '#516b83', '#ceaa78', '#78a79c'],
    'lavender_mint': ['#9473c2', '#61a99a', '#d88f9c', '#d6b55e'],
    'earth': ['#a96743', '#6e8c63', '#7795aa', '#b68c98'],
    'colorblind': ['#0072b2', '#e69f00', '#009e73', '#cc79a7'],
    'sunlight': ['#e9b82c', '#32a4b0', '#9a77b9', '#df8066'],
    'lagoon': ['#236d87', '#62aeaa', '#5575b0', '#e4b682'],
    'orchid': ['#6d548e', '#c3739f', '#7999bd', '#a9b67d'],
    'ember': ['#a8464a', '#dd9768', '#665f91', '#bca570'],
    'jade': ['#32795f', '#a2b97d', '#617faa', '#d28983'],
})
FONTS=["Arial","Times New Roman","Cambria","Calibri","Verdana","Courier New","Segoe UI"]
STYLES=[
    {"id":"standard","label":"Standard chemical drawing","description":"Conventional black bonds and explicit continuation points.","category":"chemical"},
    {"id":"foreground_pore","label":"Foreground pore","description":"Trace the first pore in the building-block palette and fade neighbors beyond its boundary.","category":"chemical"},
    {"id":"source","label":"Source drawing style","description":"Retain recorded colors, halos, and faded regions from the source-associated drawing.","category":"chemical"},
    {"id":"mono_dashed","label":"Fine black · dashed edges","description":"Fine lines and dashed continuation bonds.","category":"chemical"},
    {"id":"mono_wave","label":"Black · wavy edges","description":"Transverse wavy continuation marks.","category":"chemical"},
    {"id":"red_blue","label":"Colored building blocks","description":"Component colors with wavy continuation cuts.","category":"chemical"},
    {"id":"serif","label":"Serif chemical labels","description":"Serif labels, fine lines and transverse edge cuts.","category":"chemical"},
    {"id":"compact_labels","label":"Condensed substituents","description":"OCH3, OEt and SO3H aliases bound to their exact atoms.","category":"chemical"},
    {"id":"hetero_colors","label":"Colored heteroatoms","description":"Blue nitrogen and red oxygen labels on black bonds.","category":"chemical"},
    {"id":"element_compact","label":"Element colors · compact labels","description":"Jewel-tone heteroatoms with inferable functional-group H labels hidden.","category":"chemical"},
    {"id":"node_highlight","label":"Building-block highlights","description":"Translucent colored envelopes behind full chemical structures.","category":"chemical"},
    {"id":"ring_fill","label":"Filled aromatic rings","description":"Soft colored fills inside small aromatic rings.","category":"chemical"},
    {"id":"pore_fill","label":"Filled pore","description":"A translucent pore interior beneath a chemical outline.","category":"chemical"},
    {"id":"pastel_blocks","label":"Pastel block regions","description":"Rounded pastel regions identify constituent building blocks.","category":"chemical"},
    {"id":"bold_color","label":"Bold colored structure","description":"Heavier bonds and bold labels with a selectable palette.","category":"chemical"},
    {"id":"stacked_eclipsed","label":"Transparent stacked layers","description":"Aligned decorative copies behind the front chemical layer.","category":"layered"},
    {"id":"stacked_slipped","label":"Offset stacked layers","description":"Illustrative translated layers with adjustable transparency.","category":"layered"},
    {"id":"perspective","label":"Tilted sheet projection","description":"An oblique drawing projection of the same finite patch.","category":"layered"},
    {"id":"cartoon_nodes","label":"Cartoon nodes and struts","description":"Geometric building-block symbols with labeled connections.","category":"schematic"},
    {"id":"cartoon_overlay","label":"Chemical + cartoon regions","description":"Full chemical drawing over soft building-block symbols.","category":"schematic"},
    {"id":"glow_outline","label":"Soft molecular glow","description":"Soft colored light follows the chemical bonds while preserving the sharp drawing.","category":"chemical"},
    {"id":"molecular_spheres","label":"Molecular spheres illustration","description":"Shaded atoms and rods on the existing 2D coordinates; an illustrative molecular projection.","category":"layered"},
    {"id":"gradient_stack","label":"Gradient layers","description":"A sharp chemical drawing over tilted sheets with depth and color gradients.","category":"layered",
     "flat_foreground":True,"design_reference":"https://www.sciencedirect.com/science/article/pii/S2451929423003662"},
    {"id":"sphere_cutaway","label":"Spheres and bonds","description":"Shaded atoms fade into a fully readable chemical drawing.","category":"layered",
     "flat_foreground":True,"design_reference":"https://doi.org/10.1021/jacs.5c01990"},
    {"id":"pore_perspective","label":"Through the pore","description":"Look into a channel of translucent, shaded atoms beneath the chemical drawing.","category":"layered","flat_foreground":True},
    {"id":"half_space_stack","label":"Half-space-filling layers","description":"A sharp chemical drawing over translucent layers of shaded atoms at half van der Waals radius.","category":"layered","flat_foreground":True},
]
for row in STYLES:
    row["preview_kind"]=row["id"]
    row["representation_scope"]=("Illustrative layers/projection; not a crystal structure or CIF." if row["category"]=="layered" else
                                  "Building-block schematic; consult monomers and chemical exports for atom-level structure." if row["category"]=="schematic" else
                                  "Chemical depiction; presentation does not change connectivity.")

# Pair colors and line weights with the presentation rather than drawing each
# setting independently. The saved resolved style is still fully editable.
ELEMENT_SCHEMES = {
    "conventional": {7:"#245bd1", 8:"#d82c2c", 9:"#268a4b", 16:"#b18d0c", 15:"#c05c22", 5:"#bd598c", 17:"#268a4b", 35:"#93552b", 53:"#823c9c"},
    "deep_ink": {7:"#244487", 8:"#ad3540", 9:"#207458", 16:"#946c09", 15:"#a95522", 5:"#9e557e", 17:"#437b38", 35:"#86432b", 53:"#704582"},
    "jewel": {7:"#5944a6", 8:"#c64059", 9:"#10857d", 16:"#ab7213", 15:"#c36225", 5:"#9e528c", 17:"#317544", 35:"#98452f", 53:"#7f3c80"},
}
RANDOM_STYLES = (
    ('hetero_colors', ('classic', 'ink')),
    ('gradient_stack', ('lagoon', 'orchid', 'ember', 'jade')),
    ('half_space_stack', ('ocean_coral', 'lavender_mint', 'purple_gold')),
    ('sphere_cutaway', ('lagoon', 'orchid', 'jade', 'rose_slate')),
    ('pastel_blocks', ('lavender_mint', 'pastel', 'forest_berry', 'ember')),
    ('ring_fill', ('cobalt_citrus', 'orchid', 'jade', 'ocean_coral')),
    ('red_blue', ('colorblind', 'rose_slate', 'earth', 'lagoon')),
    ('glow_outline', ('sunlight', 'ocean_coral', 'lavender_mint')),
)

STYLE_FAMILIES = {
    **dict.fromkeys(('standard', 'hetero_colors', 'red_blue', 'serif', 'bold_color',
                     'mono_dashed', 'compact_labels', 'element_compact'), 'lines'),
    **dict.fromkeys(('gradient_stack', 'half_space_stack', 'pore_perspective', 'stacked_eclipsed', 'stacked_slipped'), 'layers'),
    **dict.fromkeys(('node_highlight', 'pastel_blocks', 'cartoon_overlay'), 'block regions'),
    **dict.fromkeys(('molecular_spheres', 'sphere_cutaway'), 'spheres'),
}


def style_configuration(style_id, seed=None, options=None, recent_styles=()):
    options=dict(options or {})
    allowed={"palette","font_family","font_scale","bond_scale","layer_count","layer_opacity","fill_opacity","edge_style","transparent_background","abbreviation_mode",
             "element_scheme", "hydrogen_labels", "custom_palette", "highlight_width", "glow_strength", "layer_spacing", "projection_tilt", "sphere_scale", "boundary_cut"}
    if set(options)-allowed: raise ValueError("Unknown style options: "+", ".join(sorted(set(options)-allowed)))
    if style_id=="random":
        seed=random.SystemRandom().randrange(2**31) if seed is None else seed
        rng=random.Random(seed)
        recent = {STYLE_FAMILIES.get(s,s) for s in recent_styles[-4:]}
        choices = [row for row in RANDOM_STYLES if STYLE_FAMILIES.get(row[0], row[0]) not in recent]
        style_id, palettes = rng.choice(choices)
        defaults={"palette":rng.choice(palettes),"font_scale":rng.choice([1.,1.05,1.1]),
                  "bond_scale":rng.choice([.95,1.1,1.25]),"fill_opacity":rng.choice([.16,.21,.27])}
        if style_id in ('gradient_stack','half_space_stack'):
            defaults.update(layer_count=rng.choice([3,4,5]), layer_spacing=rng.choice([.038,.044,.05]),
                            layer_opacity=rng.choice([.46,.55,.64]), projection_tilt=rng.choice([37,44,50]))
        elif style_id=='sphere_cutaway':
            defaults.update(sphere_scale=rng.choice([.9,1.,1.12]))
        defaults["element_scheme"] = rng.choice(list(ELEMENT_SCHEMES))
        defaults.update(options); options=defaults
    if style_id not in {s["id"] for s in STYLES}: raise ValueError("Unknown drawing style")
    defaults={"palette":"classic","font_family":"Times New Roman" if style_id=="serif" else "Arial","font_scale":1.,"bond_scale":1.,
              "layer_count":3,"layer_opacity":.24,"fill_opacity":.18,"edge_style":"dashed" if style_id=="mono_dashed" else "crossbar" if style_id=="serif" else "wavy_cross",
              "transparent_background":False,"abbreviation_mode":"condensed" if style_id=="compact_labels" else "explicit"}
    defaults.update(element_scheme="style_default", hydrogen_labels="show", boundary_cut="linkage", custom_palette=None, highlight_width=1., glow_strength=.55,
                    layer_spacing=.04, projection_tilt=35., sphere_scale=1.)
    if style_id in ('pastel_blocks', 'cartoon_overlay'):
        defaults['palette'] = 'lavender_mint'
    elif style_id == 'glow_outline':
        defaults['palette'] = 'sunlight'
    elif style_id in ('stacked_eclipsed', 'stacked_slipped', 'perspective'):
        defaults['palette'] = 'ocean_coral'
        defaults['layer_opacity'] = .38
    elif style_id in ('gradient_stack','half_space_stack'):
        defaults.update(palette='lagoon', layer_count=5, layer_spacing=.044, layer_opacity=.56, projection_tilt=43.)
    elif style_id == 'pore_perspective':
        defaults.update(palette='cobalt_citrus', layer_count=5, layer_spacing=.04, layer_opacity=.64,
                        projection_tilt=12., transparent_background=True)
    elif style_id == 'element_compact':
        defaults.update(element_scheme='jewel',hydrogen_labels='compact')
    elif style_id == 'sphere_cutaway':
        defaults.update(palette='orchid', sphere_scale=1.)
    defaults.update(options)
    if defaults["boundary_cut"] not in ('linkage', 'aryl', 'hub'):
        raise ValueError("Choose linkage, aryl or hub for the figure boundary")
    if defaults["hydrogen_labels"] not in ("show","compact"):
        raise ValueError("Choose show or compact for functional-group H labels")
    if defaults["element_scheme"] not in ("style_default", *ELEMENT_SCHEMES):
        raise ValueError("Choose a supported atom color scheme")
    if defaults["palette"] not in PALETTES or defaults["font_family"] not in FONTS: raise ValueError("Choose a supported palette and font")
    if defaults["edge_style"] not in {"wavy_cross","dashed","crossbar","bond_stub"}: raise ValueError("Unknown continuation mark")
    if defaults["abbreviation_mode"] not in {"explicit","condensed","OMe",'source_aldehyde'}: raise ValueError("Unknown abbreviation style")
    custom = defaults['custom_palette']
    if custom is not None and (not isinstance(custom, list) or not 2 <= len(custom) <= 8 or
                               any(not isinstance(c, str) or not re.fullmatch(r'#[0-9a-fA-F]{6}', c) for c in custom)):
        raise ValueError('custom_palette must contain 2–8 hexadecimal #RRGGBB colors')
    for key,lo,hi in (("font_scale",.6,1.6),("bond_scale",.5,2.5),("layer_opacity",.05,.7),("fill_opacity",.05,.6),
                      ('highlight_width',.5,2.5),('glow_strength',0.,1.),('layer_spacing',.015,.09),
                      ('projection_tilt',0.,65.),('sphere_scale',.5,1.6)):
        value=defaults[key]
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not lo<=value<=hi:
            raise ValueError(f"{key} must be between {lo} and {hi}")
    if type(defaults["layer_count"]) is not int or not 2<=defaults["layer_count"]<=5: raise ValueError("Choose 2–5 illustrative layers")
    if type(defaults["transparent_background"]) is not bool: raise ValueError("Transparency must be true or false")
    return style_id,seed,defaults


def resolved_palette(options):
    return options.get('custom_palette') or PALETTES[options['palette']]

def hull(points):
    points=sorted(set(tuple(map(float,p)) for p in points))
    if len(points)<3:return points
    def cross(o,a,b):return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    lower=[];upper=[]
    for p in points:
        while len(lower)>=2 and cross(lower[-2],lower[-1],p)<=0:lower.pop()
        lower.append(p)
    for p in reversed(points):
        while len(upper)>=2 and cross(upper[-2],upper[-1],p)<=0:upper.pop()
        upper.append(p)
    return lower[:-1]+upper[:-1]

def _segments(mol, xy, *, real_only=False):
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if real_only and not (mol.GetAtomWithIdx(a).GetAtomicNum() and mol.GetAtomWithIdx(b).GetAtomicNum()):
            continue
        yield f'M {xy[a][0]:.3f},{xy[a][1]:.3f} L {xy[b][0]:.3f},{xy[b][1]:.3f}'


def _mix(first, second, fraction):
    channels = [round((1-fraction)*int(first[i:i+2],16)+fraction*int(second[i:i+2],16)) for i in (1,3,5)]
    return '#'+''.join(f'{value:02x}' for value in channels)


def _gradient_layers(tag, defs, root, mol, xy, palette, bond, size, options, sphere_sheet=None):
    paths = ' '.join(_segments(mol, xy, real_only=True))
    if not paths:
        return
    count = options['layer_count']-1
    center = (xy.min(0)+xy.max(0))/2
    tilt = math.radians(options['projection_tilt'])
    vertical = .96-.31*(1-math.cos(tilt))
    skew = -.18*options['projection_tilt']
    spacing = min(options['layer_spacing']*size, size*.20/count)
    # Keep the front chemical layer clear where an underlying sheet crosses it.
    mask = ET.SubElement(defs,tag('mask'),{'id':'front-bond-clearance','maskUnits':'userSpaceOnUse',
                                          'x':'0','y':'0','width':str(size),'height':str(size)})
    ET.SubElement(mask,tag('rect'),{'width':str(size),'height':str(size),'fill':'white'})
    ET.SubElement(mask,tag('path'),{'d':' '.join(_segments(mol,xy)), 'fill':'none','stroke':'black',
                                     'stroke-width':str(bond*.42),'stroke-linecap':'round'})
    sheets = ET.SubElement(root,tag('g'),{'data-layer':'gradient-stack','mask':'url(#front-bond-clearance)'})
    ring_faces = ' '.join('M '+' L '.join(f'{xy[i][0]:.3f},{xy[i][1]:.3f}' for i in ring)+' Z'
                          for ring in mol.GetRingInfo().AtomRings() if 3<=len(ring)<=8)
    for k in range(count,0,-1):
        weight = k/count
        name = f'depth-ramp-{k}'
        layer = ET.SubElement(sheets,tag('g'),{'data-illustrative-layer':str(k),
            'transform':f'translate({center[0]+k*spacing*.24},{center[1]+k*spacing}) '
                        f'skewX({skew}) scale(.97,{vertical}) translate({-center[0]},{-center[1]})',
            'opacity':str(options['layer_opacity']*(1-.24*weight))})
        if sphere_sheet is not None:
            ET.SubElement(layer,tag('use'),{'{http://www.w3.org/1999/xlink}href':'#'+sphere_sheet})
            continue
        gradient = ET.SubElement(defs,tag('linearGradient'),{'id':name,'gradientUnits':'userSpaceOnUse',
            'x1':str(size*.19),'y1':str(size*.23),'x2':str(size*.84),'y2':str(size*.83)})
        for stop, color in (('0%',_mix(palette[0], '#ffffff', .16+.22*weight)),
                            ('55%', _mix(palette[1], '#ffffff', .17)),
                            ('100%', _mix(palette[-1], '#ffffff', .20))):
            ET.SubElement(gradient,tag('stop'),{'offset':stop,'stop-color':color})
        base = {'d':paths,'fill':'none','stroke-linecap':'round','stroke-linejoin':'round'}
        ET.SubElement(layer,tag('path'),{**base,'stroke':_mix(palette[0],'#253342',.48),
            'stroke-width':str(bond*.68),'transform':f'translate(0,{bond*.086})'})
        ET.SubElement(layer,tag('path'),{'d':ring_faces,'fill':f'url(#{name})'})
        ET.SubElement(layer,tag('path'),{**base,'stroke':f'url(#{name})','stroke-width':str(bond*.65)})


def _sphere_layer(tag, defs, mol, xy, colors, audit, bond, options, *, cutaway=False, half_space=False, component_carbon=False):
    from rdkit import Chem
    atom_colors={6:'#6c8190',7:'#476fcd',8:'#db5d59',16:'#d3ac41',5:'#c99a7d',9:'#61b887',17:'#61b887'}
    atom_colors.update(ELEMENT_SCHEMES.get(options.get("element_scheme"), {}))
    for color in sorted(set(atom_colors.values())|set(colors.values())):
        gradient=ET.SubElement(defs,tag('radialGradient'),{'id':'sphere-'+color[1:],'cx':'29%','cy':'24%','r':'76%'})
        for stop, shade in (('0%','#ffffff'),('35%',color),('100%',_mix(color,'#22313b',.72))):
            ET.SubElement(gradient,tag('stop'),{'offset':stop,'stop-color':shade})
    spheres=ET.Element(tag('g'),{'data-layer':'illustrative-spheres','data-coordinate-scope':'Illustrative shaded atoms on drawing coordinates'})
    ET.SubElement(spheres,tag('path'),{'d':' '.join(_segments(mol,xy,real_only=cutaway or half_space)),
        'stroke':'#91a5ad','stroke-width':str(bond*(.27 if cutaway else .16)),'stroke-linecap':'round','fill':'none'})
    for i,atom in enumerate(mol.GetAtoms()):
        p=xy[i]; atomic=atom.GetAtomicNum()
        if not atomic and (cutaway or half_space):
            continue
        component=colors[audit['atom_origin'][i]['component']]
        color=component if (cutaway or component_carbon) and atomic==6 else atom_colors.get(atomic,component)
        radius=bond*(.43 if cutaway else .22 if atomic==6 else .26)*options['sphere_scale']
        if half_space:
            radius = .5 * Chem.GetPeriodicTable().GetRvdw(atomic) * bond / 1.4 * options['sphere_scale']
        attrs={'cx':str(p[0]),'cy':str(p[1]),'r':str(radius),'stroke':color,'stroke-width':'.5','data-atom':str(i)}
        if atomic:
            attrs['fill']='url(#sphere-'+color[1:]+')'
        else:
            attrs.update(fill='none',r=str(bond*.10),stroke='#6d7b83',**{'stroke-dasharray':'2,2'})
        ET.SubElement(spheres,tag('circle'),attrs)
    if cutaway:
        lo, hi = xy.min(0), xy.max(0)
        ramp = ET.SubElement(defs,tag('linearGradient'),{'id':'cutaway-fade','gradientUnits':'userSpaceOnUse',
            'x1':str(lo[0]),'y1':str(lo[1]),'x2':str(hi[0]),'y2':str(lo[1]+.39*(hi[1]-lo[1]))})
        for offset, opacity in (('0%','0.94'),('34%','0.80'),('67%','0'),('100%','0')):
            ET.SubElement(ramp,tag('stop'),{'offset':offset,'stop-color':'white','stop-opacity':opacity})
        mask=ET.SubElement(defs,tag('mask'),{'id':'sphere-cutaway-mask','maskUnits':'userSpaceOnUse',
            'x':'0','y':'0','width':'100%','height':'100%'})
        ET.SubElement(mask,tag('rect'),{'width':'100%','height':'100%','fill':'url(#cutaway-fade)'})
        spheres.set('mask','url(#sphere-cutaway-mask)')
    return spheres


def decorate_svg(svg,mol,audit,recipe,xy,style_id,options,size):
    ET.register_namespace('',"http://www.w3.org/2000/svg")
    tag=lambda name:"{http://www.w3.org/2000/svg}"+name
    original=ET.fromstring(svg)
    chemical=ET.Element(tag('g'),{"data-layer":"chemical-graph"})
    for child in original:
        if child.tag.endswith('defs'):continue
        if child.tag.endswith('rect') and child.get('width')=='100%':continue
        chemical.append(deepcopy(child))
    root=ET.Element(tag('svg'),{"width":str(size),"height":str(size),"viewBox":f"0 0 {size} {size}",
        'data-drawing-style':style_id,'data-graph-sha256':audit['fragment_serialization_sha256'],
        'data-boundary-cut':options['boundary_cut'], 'data-figure-atoms':str(mol.GetNumAtoms())})
    defs=ET.SubElement(root,tag('defs'))
    for child in original:
        if child.tag.endswith('defs'):
            for definition in child: defs.append(deepcopy(definition))
    if not options['transparent_background']:ET.SubElement(root,tag('rect'),{"width":"100%","height":"100%","fill":"white"})
    full=ET.Element(tag('g'))
    names=list(dict.fromkeys(o['instance'] for o in audit['atom_origin']))
    components=list(dict.fromkeys(o['component'] for o in audit['atom_origin']))
    palette=resolved_palette(options)
    colors={c:palette[i%len(palette)] for i,c in enumerate(components)}
    source_style=recipe.get('source_style',{}) if style_id=='source' else {}
    colors.update(source_style.get('component_colors',{}))
    fills=source_style.get('component_fills',{})
    lengths=[np.linalg.norm(np.asarray(xy)[b.GetBeginAtomIdx()]-np.asarray(xy)[b.GetEndAtomIdx()]) for b in mol.GetBonds()]
    bond_length=float(np.median(lengths)) if lengths else size*.04
    if style_id in ('gradient_stack','half_space_stack'):
        sheet_id = None
        if style_id == 'half_space_stack':
            sheet_id = 'half-space-sheet'
            sheet = _sphere_layer(tag,defs,mol,xy,colors,audit,bond_length,options,half_space=True)
            sheet.set('id',sheet_id)
            defs.append(sheet)
        _gradient_layers(tag,defs,root,mol,xy,palette,bond_length,size,options,sheet_id)
    if style_id == 'pore_perspective':
        sheet = _sphere_layer(tag, defs, mol, xy, colors, audit, bond_length, options,
                              half_space=True, component_carbon=True)
        sheet.set('id', 'perspective-atoms')
        defs.append(sheet)
        # Parallel sheets seen by a pinhole camera shrink about a common vanishing point.
        # Keep the front-plane coordinates intact for the chemical drawing and picking.
        center = (xy.min(0)+xy.max(0))/2
        tilt = math.radians(options['projection_tilt'])
        focus = center + size*math.sin(tilt)*np.array([.24, -.16])
        stack = ET.SubElement(root, tag('g'), {'data-layer':'pore-perspective'})
        count = options['layer_count']
        for depth in reversed(range(count)):
            factor = 1/(1+depth*options['layer_spacing']*2.8)
            opacity = options['layer_opacity']*(1-.55*depth/(count-1)) if depth else .23
            ET.SubElement(stack, tag('use'), {
                '{http://www.w3.org/1999/xlink}href':'#perspective-atoms',
                'data-illustrative-layer':str(depth), 'opacity':str(opacity),
                'transform':f'translate({focus[0]:.3f},{focus[1]:.3f}) scale({factor:.6f}) '
                            f'translate({-focus[0]:.3f},{-focus[1]:.3f})'})
    halo=source_style.get('halo')
    if style_id=='glow_outline' or halo:
        blur=ET.SubElement(defs,tag('filter'),{'id':'molecular-glow','x':'-30%','y':'-30%','width':'160%','height':'160%'})
        ET.SubElement(blur,tag('feGaussianBlur'),{'stdDeviation':str(bond_length*.13*options['highlight_width'])})
        glow=ET.SubElement(full,tag('g'),{'data-layer':'molecular-glow','filter':'url(#molecular-glow)',
                                        'opacity':str(halo.get('opacity',.55) if halo else options['glow_strength'])})
        # Trace bonds, rather than a convex disk, so light follows the chemical
        # silhouette and leaves the pore clear.
        for b in mol.GetBonds():
            a,c=b.GetBeginAtomIdx(),b.GetEndAtomIdx();p,q=xy[a],xy[c]
            color=halo.get('color','#f3d74d') if halo else colors[audit['atom_origin'][a]['component']]
            ET.SubElement(glow,tag('path'),{'d':f'M {p[0]},{p[1]} L {q[0]},{q[1]}','fill':'none','stroke':color,
                'stroke-width':str(bond_length*(halo.get('radius_bonds',.22)*2 if halo else .52)*options['highlight_width']),
                'stroke-linecap':'round','stroke-linejoin':'round'})
    centers={}; owner={}
    for name in names:
        indices=[i for i,o in enumerate(audit['atom_origin']) if o['instance']==name and mol.GetAtomWithIdx(i).GetAtomicNum()]
        if not indices:
            continue  # A standalone continuation port is not a building block.
        centers[name]=np.asarray(xy)[indices].mean(0); owner[name]=audit['atom_origin'][indices[0]]['component']
        points=np.asarray(xy)[indices]
        if style_id in ('node_highlight','pastel_blocks','cartoon_overlay') or owner[name] in fills:
            if style_id=='pastel_blocks':
                # Rounded envelopes follow the building block orientation;
                # axis-aligned boxes falsely suggest rectangular monomers.
                outline=hull(points)
                ET.SubElement(full,tag('polygon'),{'points':' '.join(f'{x:.3f},{y:.3f}' for x,y in outline),
                    'fill':colors[owner[name]],'stroke':colors[owner[name]],'stroke-width':str(bond_length*.68*options['highlight_width']),
                    'stroke-linejoin':'round','opacity':str(options['fill_opacity']), 'data-highlight-instance':name})
            elif style_id=='cartoon_overlay':
                radius=float(max(np.linalg.norm(points-centers[name],axis=1)))*1.1
                ET.SubElement(full,tag('circle'),{"cx":str(centers[name][0]),"cy":str(centers[name][1]),"r":str(radius),"fill":colors[owner[name]],"opacity":str(options['fill_opacity']),"stroke":colors[owner[name]],"stroke-width":str(size*.003)})
            else:
                outline=hull(points)
                fill=fills.get(owner[name],colors[owner[name]])
                ET.SubElement(full,tag('polygon'),{"points":" ".join(f"{x:.3f},{y:.3f}" for x,y in outline),"fill":fill,"opacity":str(options['fill_opacity']),"stroke":fill,"stroke-width":str(bond_length*.45*options['highlight_width']),"stroke-linejoin":"round",'data-highlight-instance':name})
    if style_id=='ring_fill':
        for ring in mol.GetRingInfo().AtomRings():
            if not 3<=len(ring)<=8 or not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):continue
            pts=np.asarray(xy)[list(ring)]; center=pts.mean(0); pts=center+.92*(pts-center)
            ET.SubElement(full,tag('polygon'),{"points":" ".join(f"{x:.3f},{y:.3f}" for x,y in pts),"fill":colors[audit['atom_origin'][ring[0]]['component']],"opacity":str(options['fill_opacity'])})
    face_colors = recipe.get('repeat_face_colors',[])
    faces = recipe.get('repeat_faces')
    if ((style_id=='pore_fill' or face_colors) and len(centers)>=3 and
            (faces or recipe.get('layout_kind') not in ('open_chain','net_patch'))):
        if faces is None:
            ids = [i['id'] for i in recipe.get('instances',[]) if i['id'] in centers]
            polygons = [[centers[i] for i in ids] if len(ids)>=3 else hull(list(centers.values()))]
        else:
            polygons = [[centers[i] for i in face] for face in faces]
        for i,vertices in enumerate(polygons):
            pts=np.asarray(vertices)
            center=pts.mean(0);pts=center+.82*(pts-center)
            fill=face_colors[i] if face_colors else palette[1]
            ET.SubElement(full,tag('polygon'),{"points":" ".join(f"{x:.3f},{y:.3f}" for x,y in pts),
                "fill":fill,"opacity":str(options['fill_opacity']),"data-pore-fill":"true"})
    if style_id=='cartoon_nodes':
        connections={(edge['a']['instance'],edge['b']['instance']) for edge in recipe.get('connections',[])}
        if not connections:
            connections={(audit['atom_origin'][b.GetBeginAtomIdx()]['instance'],audit['atom_origin'][b.GetEndAtomIdx()]['instance']) for b in mol.GetBonds()}
        for first,second in sorted(connections):
            if first==second or first not in centers or second not in centers:continue
            a,b=centers[first],centers[second]
            ET.SubElement(full,tag('line'),{"x1":str(a[0]),"y1":str(a[1]),"x2":str(b[0]),"y2":str(b[1]),"stroke":"#45545d","stroke-width":str(size*.011),"stroke-linecap":"round"})
        for name,center in centers.items():
            component=recipe.get('components',{}).get(owner[name],{'label':owner[name]}); count=component.get('connectivity',2);radius=size*.03
            if count==2:
                ET.SubElement(full,tag('circle'),{"cx":str(center[0]),"cy":str(center[1]),"r":str(radius*.75),"fill":colors[owner[name]]})
            else:
                pts=[center+radius*np.array([math.cos(-math.pi/2+2*math.pi*i/count),math.sin(-math.pi/2+2*math.pi*i/count)]) for i in range(min(count,12))]
                ET.SubElement(full,tag('polygon'),{"points":" ".join(f"{x},{y}" for x,y in pts),"fill":colors[owner[name]]})
            label=ET.SubElement(full,tag('text'),{"x":str(center[0]),"y":str(center[1]+radius+size*.019),"text-anchor":"middle","font-family":options['font_family'],"font-size":str(size*.018),"fill":"#283c40"})
            label.text=str(component.get('label',owner[name]))[:24]
    elif style_id in ('molecular_spheres','sphere_cutaway'):
        cutaway=style_id=='sphere_cutaway'
        full.append(_sphere_layer(tag,defs,mol,xy,colors,audit,bond_length,options,cutaway=cutaway))
        if cutaway:
            full.append(chemical)
    else:full.append(chemical)
    if style_id in ('stacked_eclipsed','stacked_slipped','perspective'):
        count=options['layer_count'] if style_id!='perspective' else 1
        tilt=math.radians(options['projection_tilt']);vertical=max(.34,math.cos(tilt))
        spacing=size*options['layer_spacing'];total=(count-1)*spacing
        fit=min(.90,(size*.84)/(size*.9+total*1.6))
        for i in reversed(range(count)):
            layer=deepcopy(full);offset=(i-(count-1)/2)*spacing;slip=offset*.8 if style_id=='stacked_slipped' else 0
            layer.set('transform',f"translate({slip},{-offset}) translate({size/2},{size/2}) rotate(-10) skewX(-18) scale({fit},{fit*vertical}) translate({-size/2},{-size/2})")
            layer.set('data-illustrative-layer',str(i))
            if i:layer.set('opacity',str(options['layer_opacity']*(1-.12*(i-1))))
            root.append(layer)
    else:root.append(full)
    return ET.tostring(root,encoding='unicode')
