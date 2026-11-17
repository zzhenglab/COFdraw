"""Source-coordinate chemical redraws, separate from graph recognition.

Explicit source anchors and an audited graph drive SVG primitives. The source
pixels are not traced into the output. This is a depiction transcription,
not independent OCR or a new real-image example.
"""
from __future__ import annotations
import hashlib,html,json,math,re
from pathlib import Path
import numpy as np
from rdkit import Chem
from .depiction_rings import symmetric_depiction_rings

def _key(value):
    return value if isinstance(value,str) else f"{value['instance']}:{value['fragment_atom_index']}"

def _tuple(value):
    label=_key(value);instance,index=label.rsplit(':',1);return instance,int(index)

def _label_svg(label, mode='tspan'):
    if mode=='unicode':return html.escape(label.translate(str.maketrans('0123456789','₀₁₂₃₄₅₆₇₈₉')))
    return re.sub(r'(\d+)',r'<tspan font-size="70%" baseline-shift="sub">\1</tspan>',html.escape(label))


def graph_coordinates(mol, audit, annotation):
    positions = {}
    for row in annotation['atom_coordinates']:
        key = (row['instance'], row['fragment_atom_index'])
        if key in positions:
            raise ValueError('Duplicate source atom coordinate')
        positions[key] = np.asarray(row['xy'],dtype=float)
    for group in annotation.get('style',{}).get('collapsed_groups',[]):
        anchor=_tuple(group['anchor'])
        for row in group['hidden_atoms']:
            positions[_tuple(row)]=positions[anchor]
    coords=[]
    for origin in audit['atom_origin']:
        key=(origin['instance'],origin['fragment_atom_index'])
        if key not in positions:
            raise ValueError(f'Missing source coordinate {key}')
        coords.append(positions[key])
    if len(coords)!=mol.GetNumAtoms() or not np.isfinite(coords).all():
        raise ValueError('Invalid source coordinates')
    return np.asarray(coords)


def render_source_svg(mol, audit, annotation, output):
    """Source-sized vector rendering with chemical bonds, labels and ports.

    Aromatic resonance is rendered in one valid Kekule form. Coordinates alone
    do not introduce stereochemical labels. Collapsed display groups need an
    explicit graph atom list; graph exports remain uncollapsed.
    """
    coords=graph_coordinates(mol,audit,annotation)
    style=annotation.get('style',{});w,h=annotation['canvas_size']
    width=float(style.get('bond_width',1.6));font=float(style.get('atom_font_size',14))
    colors=style.get('component_colors',{});overrides=style.get('atom_labels',{})
    key=lambda idx:f"{audit['atom_origin'][idx]['instance']}:{audit['atom_origin'][idx]['fragment_atom_index']}"
    color=lambda idx:style.get('atom_colors',{}).get(key(idx),colors.get(audit['atom_origin'][idx]['component'],'#202020'))
    drawmol=Chem.Mol(mol);Chem.Kekulize(drawmol,clearAromaticFlags=True)
    symmetric_depiction_rings(mol)
    rings=list(mol.GetRingInfo().AtomRings())
    hidden=set();group_labels={};group_options={}
    aromatic_choices={frozenset(pair) for pair in style.get('aromatic_double_bonds',[])}
    aromatic_scope=set().union(*aromatic_choices) if aromatic_choices else set()
    if aromatic_choices:
        lookup={key(i):i for i in range(mol.GetNumAtoms())}
        for pair in aromatic_choices:
            if len(pair)!=2 or not pair.issubset(lookup):raise ValueError('Invalid aromatic drawing bond')
            first,second=[lookup[k] for k in pair];bond=mol.GetBondBetweenAtoms(first,second)
            if bond is None or not bond.GetIsAromatic():raise ValueError('Depiction resonance override is restricted to aromatic graph bonds')
        if any(sum(k in p for p in aromatic_choices)!=1 for k in aromatic_scope):
            raise ValueError('Aromatic depiction bonds must form a nonoverlapping matching')
    for group in style.get('collapsed_groups',[]):
        anchor=_key(group['anchor']);group_labels[anchor]=group['label'];group_options[anchor]=group
        hidden.update(_key(x) for x in group['hidden_atoms'])
    body=[]
    def line(p,q,col,extra=''):
        body.append(f'<path d="M {p[0]:.3f},{p[1]:.3f} L {q[0]:.3f},{q[1]:.3f}" fill="none" stroke="{col}" stroke-width="{width}" stroke-linecap="round" {extra}/>')
    for bond in drawmol.GetBonds():
        a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx();p,q=coords[a],coords[b]
        if key(a) in hidden or key(b) in hidden:continue
        delta=q-p;length=np.linalg.norm(delta)
        if length<.2:raise ValueError('Coincident bonded source atoms')
        tangent=delta/length;normal=np.array([-tangent[1],tangent[0]])
        acol,bcol=color(a),color(b)
        bond_override=style.get('bond_colors',{}).get('|'.join(sorted([key(a),key(b)])))
        if bond_override:acol=bcol=bond_override
        isport=mol.GetAtomWithIdx(a).GetAtomicNum()==0 or mol.GetAtomWithIdx(b).GetAtomicNum()==0
        if isport:
            port=a if mol.GetAtomWithIdx(a).GetAtomicNum()==0 else b
            start=q if port==a else p;end=p if port==a else q
            kind=style.get('port_styles',{}).get(key(port),style.get('port_style','dashed'))
            port_color=style.get('port_colors',{}).get(key(port),color(b if port==a else a))
            if kind=='hidden':
                # The graph port can be inferred beyond a visible terminal
                # atom even when the source draws no continuation stroke.
                pass
            elif kind=='bond_stub':
                # Unmarked open imine/alkyne ends retain the visible bond
                # multiplicity; a double stub must not become a single dash.
                order=mol.GetBondBetweenAtoms(a,b).GetBondTypeAsDouble()
                if order not in (1,2,3):raise ValueError('Bond stub requires an explicit bond order')
                d=end-start;n=np.array([-d[1],d[0]])/np.linalg.norm(d)
                shifts=[np.zeros(2)] if order==1 else [-n*float(style.get('double_offset',1.5)),n*float(style.get('double_offset',1.5))] if order==2 else [-n*float(style.get('triple_offset',3)),np.zeros(2),n*float(style.get('triple_offset',3))]
                for shift in shifts:line(start+shift,end+shift,port_color)
            elif kind=='crossbar':
                # A transverse end mark is a continuation symbol, not a
                # second chemical bond or a finite substituent.
                order=mol.GetBondBetweenAtoms(a,b).GetBondTypeAsDouble()
                if order not in (1,2,3):raise ValueError('Crossbar continuation requires an explicit bond order')
                d=end-start;n=np.array([-d[1],d[0]])/np.linalg.norm(d)
                center=start+float(style.get('crossbar_fraction',1.))*d
                span=float(style.get('crossbar_span',font*.9))
                if not math.isfinite(span) or span<=0:raise ValueError('Invalid crossbar span')
                shifts=[np.zeros(2)] if order==1 else [-n*float(style.get('double_offset',1.5)),n*float(style.get('double_offset',1.5))] if order==2 else [-n*float(style.get('triple_offset',3)),np.zeros(2),n*float(style.get('triple_offset',3))]
                for shift in shifts:line(start+shift,end+shift,port_color)
                line(center-n*span/2,center+n*span/2,style.get('crossbar_color',port_color))
            elif kind=='wavy_cross':
                # Literature often draws a short wave ACROSS an ordinary bond.
                # This is a continuation marker, not an unspecified stereo bond.
                d=end-start;l=np.linalg.norm(d);t=d/l;n=np.array([-t[1],t[0]])
                order=mol.GetBondBetweenAtoms(a,b).GetBondTypeAsDouble()
                if order not in (1,2,3):raise ValueError('Wavy continuation requires an explicit bond order')
                shifts=[np.zeros(2)] if order==1 else [-n*float(style.get('double_offset',1.5)),n*float(style.get('double_offset',1.5))] if order==2 else [-n*float(style.get('triple_offset',3)),np.zeros(2),n*float(style.get('triple_offset',3))]
                for shift in shifts:line(start+shift,end+shift,port_color)
                center=start+float(style.get('wave_cross_fraction',.72))*d
                span=float(style.get('wave_cross_span',font*.9))
                amplitude=float(style.get('wave_amplitude',1.5))
                points=[center+(u-.5)*span*n+amplitude*math.sin(u*math.pi*6)*t for u in np.linspace(0,1,41)]
                path='M '+' L '.join(f'{v[0]:.3f},{v[1]:.3f}' for v in points)
                body.append(f'<path d="{path}" fill="none" stroke="{port_color}" stroke-width="{width}"/>')
            elif kind=='wavy':
                d=end-start;l=np.linalg.norm(d);n=np.array([-d[1],d[0]])/l
                points=[start+t*d+n*float(style.get('wave_amplitude',1.5))*math.sin(t*math.pi*6) for t in np.linspace(0,1,41)]
                path='M '+' L '.join(f'{v[0]:.3f},{v[1]:.3f}' for v in points)
                body.append(f'<path d="{path}" fill="none" stroke="{port_color}" stroke-width="{width}"/>')
            else:
                pattern=style.get('port_dasharray','1.5,3')
                if kind=='dashed':
                    # Explicit segments survive SVG-to-PDF rasterizers that
                    # ignore stroke-dasharray on presentation attributes.
                    values=[float(x) for x in (re.split(r'[,\s]+',pattern.strip()) if isinstance(pattern,str) else pattern)]
                    if not values or any(not math.isfinite(v) or v<=0 for v in values):raise ValueError('Invalid continuation dash pattern')
                    if len(values)%2:values*=2
                    d=end-start;l=np.linalg.norm(d);t=d/l;offset=0.;j=0
                    while offset<l:
                        nxt=min(l,offset+values[j%len(values)])
                        if j%2==0:line(start+offset*t,start+nxt*t,port_color)
                        offset=nxt;j+=1
                else:line(start,end,port_color)
            continue
        order=bond.GetBondTypeAsDouble()
        if key(a) in aromatic_scope and key(b) in aromatic_scope and mol.GetBondBetweenAtoms(a,b).GetIsAromatic():
            order=2 if frozenset([key(a),key(b)]) in aromatic_choices else 1
        if order==1:
            if acol==bcol:line(p,q,acol)
            else:
                mid=(p+q)/2;line(p,mid,acol);line(mid,q,bcol)
        elif order==2:
            local=[ring for ring in rings if a in ring and b in ring and len(ring)<=8]
            if local:
                ring=min(local,key=len);center=coords[list(ring)].mean(axis=0)
                if np.dot(center-(p+q)/2,normal)<0:normal=-normal
                off=float(style.get('ring_double_offset',3.4))*normal
                trim=float(style.get('ring_double_trim',.18))
                line(p,q,acol);line(p+off+tangent*length*trim,q+off-tangent*length*trim,acol)
            else:
                off=normal*float(style.get('double_offset',1.5))
                for shift in [off,-off]:
                    if acol==bcol:line(p+shift,q+shift,acol)
                    else:
                        mid=(p+q)/2;line(p+shift,mid+shift,acol);line(mid+shift,q+shift,bcol)
        elif order==3:
            offset=float(style.get('triple_offset',3))
            for shift in [-normal*offset,normal*offset,np.zeros(2)]:line(p+shift,q+shift,acol)
        else:raise ValueError('Unsupported source drawing bond order')
    for atom in mol.GetAtoms():
        idx=atom.GetIdx();symbol=atom.GetSymbol()
        if key(idx) in hidden:continue
        if atom.GetAtomicNum()==0:continue
        if atom.GetAtomicNum()==6 and key(idx) not in group_labels and key(idx) not in overrides:continue
        hydrogens=atom.GetTotalNumHs()
        hydrogen_suffix=('H'+(str(hydrogens) if hydrogens>1 else '')) if hydrogens else ''
        label=group_labels.get(key(idx),overrides.get(key(idx),symbol+hydrogen_suffix))
        if atom.GetFormalCharge():label+=('+' if atom.GetFormalCharge()>0 else '-')
        options=group_options.get(key(idx),{})
        fs=options.get('font_size',style.get('atom_font_sizes',{}).get(key(idx),font))
        weight=options.get('font_weight',style.get('atom_font_weights',{}).get(key(idx),style.get('font_weight','normal')))
        x,y=options.get('label_xy',coords[idx]);text_width=max(fs*.58,len(label)*fs*.56)
        body.append(f'<rect x="{x-text_width/2-1:.3f}" y="{y-fs*.52:.3f}" width="{text_width+2:.3f}" height="{fs*1.04:.3f}" fill="white"/>')
        family=options.get('font_family',style.get('font_family','Arial'));label_color=options.get('color',color(idx))
        # Use an explicit alphabetic baseline: PDF SVG renderers may ignore
        # dominant-baseline, letting bonds show through the top of atom glyphs.
        baseline=y+fs*.35
        body.append(f'<text x="{x:.3f}" y="{baseline:.3f}" text-anchor="middle" font-family="{html.escape(family)}" font-style="{"italic" if options.get("italic") else "normal"}" font-weight="{html.escape(str(weight))}" font-size="{fs}" fill="{label_color}">{_label_svg(label,style.get("subscript_style","tspan"))}</text>')
    for label in style.get('extra_labels',[]):
        x,y=label['xy'];fs=label.get('font_size',font)
        weight=label.get('font_weight',style.get('font_weight','normal'))
        body.append(f'<text x="{x}" y="{y}" text-anchor="middle" dominant-baseline="central" font-family="{html.escape(label.get("font_family",style.get("font_family","Arial")))}" font-weight="{html.escape(str(weight))}" font-size="{fs}" fill="{label.get("color","#202020")}">{html.escape(label["text"])}</text>')
    if annotation.get('center_label'):
        label=annotation['center_label'];x,y=label['xy']
        body.append(f'<text x="{x}" y="{y}" text-anchor="middle" dominant-baseline="central" font-family="Times New Roman" font-size="{label.get("font_size",34)}" font-weight="bold">{html.escape(label["text"])}</text>')
    svg=f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}"><rect width="100%" height="100%" fill="white"/>'+''.join(body)+'</svg>'
    Path(output).write_text(svg,encoding='utf-8')
    return {'renderer':'source-chemical-svg-1','graph_atoms':mol.GetNumAtoms(),'source_anchored_atoms':len(coords)-len(hidden),
            'graph_atoms_hidden_by_source_abbreviation':len(hidden),
            'graph_serialization_sha256':audit['fragment_serialization_sha256'],
            'svg_sha256':hashlib.sha256(svg.encode()).hexdigest(),'chemical_graph_changed':False,
            'mode':'manual_source_depiction_transcription_not_model_prediction',
            'stereochemistry':'No stereo inferred from these source coordinates. Separate graph assignments retained.'}


def overlay_metrics(source, reconstruction, *, excluded_rectangles=(), tolerance_pixels=2, pixel_threshold=160, ink_measure='min_channel'):
    """Ink registration after source-defined masks, never a chemistry metric."""
    from PIL import Image
    from scipy.ndimage import distance_transform_edt
    original=np.asarray(Image.open(source).convert('RGB'))
    rendered=np.asarray(Image.open(reconstruction).convert('RGB'))
    if original.shape!=rendered.shape:raise ValueError('Overlay images need equal source pixel dimensions')
    mask=np.ones(original.shape[:2],bool)
    for x0,y0,x1,y1 in excluded_rectangles:mask[y0:y1,x0:x1]=False
    if ink_measure=='min_channel':
        channel1,channel2=original.min(axis=2),rendered.min(axis=2)
    elif ink_measure=='luminance':
        weights=np.array([.2126,.7152,.0722]);channel1=original@weights;channel2=rendered@weights
    else:raise ValueError('Unknown ink measure')
    ink1=(channel1<pixel_threshold)&mask;ink2=(channel2<pixel_threshold)&mask
    intersection=(ink1&ink2).sum();union=(ink1|ink2).sum()
    precision=(distance_transform_edt(~ink1)[ink2]<=tolerance_pixels).mean() if ink2.any() else 0
    recall=(distance_transform_edt(~ink2)[ink1]<=tolerance_pixels).mean() if ink1.any() else 0
    f1=2*precision*recall/(precision+recall) if precision+recall else 0
    return {'ink_iou':float(intersection/union if union else 0),'tolerant_ink_precision':float(precision),
            'tolerant_ink_recall':float(recall),'tolerant_ink_f1':float(f1),'tolerance_pixels':tolerance_pixels,
            'excluded_rectangles':list(excluded_rectangles),'pixel_threshold':pixel_threshold,'ink_measure':ink_measure,
            'not_chemical_accuracy':True,'fit_uses_source_pixels':True}
