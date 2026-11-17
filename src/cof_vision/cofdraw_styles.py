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
FONTS=["Arial","Times New Roman","Cambria","Calibri","Verdana","Courier New","Segoe UI"]
STYLES=[
    {"id":"standard","label":"Standard chemical drawing","description":"Conventional black bonds and explicit continuation points.","category":"chemical"},
    {"id":"mono_dashed","label":"Fine black · dashed edges","description":"Fine lines and dashed continuation bonds.","category":"chemical"},
    {"id":"mono_wave","label":"Black · wavy edges","description":"Transverse wavy continuation marks.","category":"chemical"},
    {"id":"red_blue","label":"Colored building blocks","description":"Component colors with wavy continuation cuts.","category":"chemical"},
    {"id":"serif","label":"Serif chemical labels","description":"Serif labels, fine lines and transverse edge cuts.","category":"chemical"},
    {"id":"compact_labels","label":"Condensed substituents","description":"OCH3, OEt and SO3H aliases bound to their exact atoms.","category":"chemical"},
    {"id":"hetero_colors","label":"Colored heteroatoms","description":"Blue nitrogen and red oxygen labels on black bonds.","category":"chemical"},
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
]
for row in STYLES:
    row["preview_kind"]=row["id"]
    row["representation_scope"]=("Illustrative layers/projection; not a crystal structure or CIF." if row["category"]=="layered" else
                                  "Building-block schematic; consult monomers and chemical exports for atom-level structure." if row["category"]=="schematic" else
                                  "Chemical depiction; presentation does not change connectivity.")

def style_configuration(style_id, seed=None, options=None):
    options=dict(options or {})
    allowed={"palette","font_family","font_scale","bond_scale","layer_count","layer_opacity","fill_opacity","edge_style","transparent_background","abbreviation_mode"}
    if set(options)-allowed: raise ValueError("Unknown style options: "+", ".join(sorted(set(options)-allowed)))
    if style_id=="random":
        seed=random.SystemRandom().randrange(2**31) if seed is None else seed
        rng=random.Random(seed); style_id=rng.choice(STYLES[1:])["id"]
        defaults={"palette":rng.choice(list(PALETTES)),"font_family":rng.choice(FONTS),"font_scale":rng.choice([.9,1,1.1]),
                  "bond_scale":rng.choice([.8,1,1.3]),"fill_opacity":rng.choice([.12,.2,.3]),"layer_opacity":rng.choice([.16,.25,.35]),
                  "layer_count":rng.choice([2,3,4]),"edge_style":rng.choice(["wavy_cross","dashed","crossbar"])}
        defaults.update(options); options=defaults
    if style_id not in {s["id"] for s in STYLES}: raise ValueError("Unknown drawing style")
    defaults={"palette":"classic","font_family":"Times New Roman" if style_id=="serif" else "Arial","font_scale":1.,"bond_scale":1.,
              "layer_count":3,"layer_opacity":.24,"fill_opacity":.18,"edge_style":"dashed" if style_id=="mono_dashed" else "crossbar" if style_id=="serif" else "wavy_cross",
              "transparent_background":False,"abbreviation_mode":"condensed" if style_id=="compact_labels" else "explicit"}
    defaults.update(options)
    if defaults["palette"] not in PALETTES or defaults["font_family"] not in FONTS: raise ValueError("Choose a supported palette and font")
    if defaults["edge_style"] not in {"wavy_cross","dashed","crossbar","bond_stub"}: raise ValueError("Unknown continuation mark")
    if defaults["abbreviation_mode"] not in {"explicit","condensed","OMe"}: raise ValueError("Unknown abbreviation style")
    for key,lo,hi in (("font_scale",.6,1.6),("bond_scale",.5,2.5),("layer_opacity",.05,.7),("fill_opacity",.05,.6)):
        value=defaults[key]
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not lo<=value<=hi:
            raise ValueError(f"{key} must be between {lo} and {hi}")
    if type(defaults["layer_count"]) is not int or not 2<=defaults["layer_count"]<=5: raise ValueError("Choose 2–5 illustrative layers")
    if type(defaults["transparent_background"]) is not bool: raise ValueError("Transparency must be true or false")
    return style_id,seed,defaults

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

def decorate_svg(svg,mol,audit,recipe,xy,style_id,options,size):
    ET.register_namespace('',"http://www.w3.org/2000/svg")
    tag=lambda name:"{http://www.w3.org/2000/svg}"+name
    original=ET.fromstring(svg)
    chemical=ET.Element(tag('g'),{"data-layer":"chemical-graph"})
    for child in original:
        if child.tag.endswith('rect') and child.get('width')=='100%':continue
        chemical.append(deepcopy(child))
    root=ET.Element(tag('svg'),{"width":str(size),"height":str(size),"viewBox":f"0 0 {size} {size}"})
    if not options['transparent_background']:ET.SubElement(root,tag('rect'),{"width":"100%","height":"100%","fill":"white"})
    full=ET.Element(tag('g'))
    names=list(dict.fromkeys(o['instance'] for o in audit['atom_origin']))
    components=list(dict.fromkeys(o['component'] for o in audit['atom_origin']))
    colors={c:PALETTES[options['palette']][i%4] for i,c in enumerate(components)}
    centers={}; owner={}
    for name in names:
        indices=[i for i,o in enumerate(audit['atom_origin']) if o['instance']==name and mol.GetAtomWithIdx(i).GetAtomicNum()]
        centers[name]=np.asarray(xy)[indices].mean(0); owner[name]=audit['atom_origin'][indices[0]]['component']
        points=np.asarray(xy)[indices]
        if style_id in ('node_highlight','pastel_blocks','cartoon_overlay'):
            if style_id=='pastel_blocks':
                lo=points.min(0)-size*.009;hi=points.max(0)+size*.009
                ET.SubElement(full,tag('rect'),{"x":str(lo[0]),"y":str(lo[1]),"width":str(hi[0]-lo[0]),"height":str(hi[1]-lo[1]),"rx":str(size*.018),"fill":colors[owner[name]],"opacity":str(options['fill_opacity'])})
            elif style_id=='cartoon_overlay':
                radius=float(max(np.linalg.norm(points-centers[name],axis=1)))*1.1
                ET.SubElement(full,tag('circle'),{"cx":str(centers[name][0]),"cy":str(centers[name][1]),"r":str(radius),"fill":colors[owner[name]],"opacity":str(options['fill_opacity']),"stroke":colors[owner[name]],"stroke-width":str(size*.003)})
            else:
                outline=hull(centers[name]+1.14*(points-centers[name]))
                ET.SubElement(full,tag('polygon'),{"points":" ".join(f"{x:.3f},{y:.3f}" for x,y in outline),"fill":colors[owner[name]],"opacity":str(options['fill_opacity']),"stroke":colors[owner[name]],"stroke-width":str(size*.008),"stroke-linejoin":"round"})
    if style_id=='ring_fill':
        for ring in mol.GetRingInfo().AtomRings():
            if not 3<=len(ring)<=8 or not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):continue
            pts=np.asarray(xy)[list(ring)]; center=pts.mean(0); pts=center+.92*(pts-center)
            ET.SubElement(full,tag('polygon'),{"points":" ".join(f"{x:.3f},{y:.3f}" for x,y in pts),"fill":colors[audit['atom_origin'][ring[0]]['component']],"opacity":str(options['fill_opacity'])})
    if style_id=='pore_fill' and recipe.get('layout_kind')!='open_chain':
        pts=np.asarray([centers[i['id']] for i in recipe['instances']]);center=pts.mean(0);pts=center+.82*(pts-center)
        ET.SubElement(full,tag('polygon'),{"points":" ".join(f"{x:.3f},{y:.3f}" for x,y in pts),"fill":PALETTES[options['palette']][1],"opacity":str(options['fill_opacity'])})
    if style_id=='cartoon_nodes':
        for edge in recipe['connections']:
            a,b=centers[edge['a']['instance']],centers[edge['b']['instance']]
            ET.SubElement(full,tag('line'),{"x1":str(a[0]),"y1":str(a[1]),"x2":str(b[0]),"y2":str(b[1]),"stroke":"#45545d","stroke-width":str(size*.011),"stroke-linecap":"round"})
        for name,center in centers.items():
            component=recipe['components'][owner[name]]; count=component.get('connectivity',3);radius=size*.03
            if count==2:
                ET.SubElement(full,tag('circle'),{"cx":str(center[0]),"cy":str(center[1]),"r":str(radius*.75),"fill":colors[owner[name]]})
            else:
                pts=[center+radius*np.array([math.cos(-math.pi/2+2*math.pi*i/count),math.sin(-math.pi/2+2*math.pi*i/count)]) for i in range(min(count,12))]
                ET.SubElement(full,tag('polygon'),{"points":" ".join(f"{x},{y}" for x,y in pts),"fill":colors[owner[name]]})
            label=ET.SubElement(full,tag('text'),{"x":str(center[0]),"y":str(center[1]+radius+size*.019),"text-anchor":"middle","font-family":options['font_family'],"font-size":str(size*.018),"fill":"#283c40"})
            label.text=str(component.get('label',owner[name]))[:24]
    else:full.append(chemical)
    if style_id in ('stacked_eclipsed','stacked_slipped','perspective'):
        count=options['layer_count'] if style_id!='perspective' else 1
        for i in reversed(range(count)):
            layer=deepcopy(full);offset=i*size*.027;slip=i*size*.025 if style_id=='stacked_slipped' else 0
            layer.set('transform',f"translate({offset+slip},{-offset}) translate({size/2},{size/2}) rotate(-12) skewX(-24) scale(.75,.58) translate({-size/2},{-size/2})")
            if i:layer.set('opacity',str(options['layer_opacity']))
            root.append(layer)
    else:root.append(full)
    return ET.tostring(root,encoding='unicode')
