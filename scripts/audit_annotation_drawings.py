"""Build reviewed COF targets from precursor graphs and compare their bounded drawings.

The builder gets feed SMILES, a linkage and a proposed net. Source atom graphs
are used only to compare finished products. The output includes unresolved rows.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
from functools import lru_cache
import hashlib
import json
import math
from statistics import median
import os
from pathlib import Path
import re
import subprocess
import sys

from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from cof_vision.cofdraw_boundary import aryl_boundary
from cof_vision.cofdraw_assembly import assemble
from cof_vision.cofdraw_chemistry import decorated_monomers, design_signature, prepare_template
from cof_vision.cofdraw_cofex import parse_record, recipe_record
from cof_vision.cofdraw_repeat import repeat_patch
from cof_vision.cofdraw_self_kagome import hub_boundary
from cof_vision.depiction_rings import symmetric_depiction_rings
from cof_vision.drawing_api import validate_recipe
from cof_vision.framework_drawing import VERSION, build_patch, matches_linkage_family


LINKAGES = {
    'imine': 'imine', 'aldimine': 'imine', 'aryl c=n imines':'imine', 'imine as drawn':'imine',
    'phenolic imine':'imine', 'enol–imine':'imine', 'acylhydrazone': 'acylhydrazone',
    'oxalyl bis(acylhydrazone)':'acylhydrazone',
    'β-ketoenamine': 'beta_ketoenamine', 'beta-ketoenamine': 'beta_ketoenamine', 'beta_ketoenamine':'beta_ketoenamine',
    'azine': 'azine', 'secondary amine': 'secondary_amine', 'nitrone': 'nitrone',
    'cyanovinylene': 'cyano_vinylene_sp2_carbon', 'cyano-vinylene': 'cyano_vinylene_sp2_carbon',
    'α-cyanovinylene':'cyano_vinylene_sp2_carbon', 'cyano-substituted vinylene': 'cyano_vinylene_sp2_carbon',
    'vinylene': 'vinylene_methyltriazine', 'trans-vinylene': 'vinylene_methyltriazine',
    'amide': 'amide', 'imide': 'imide', 'cyclic imide': 'imide',
    'boronate ester':'boronate_ester', 'boronate-ester':'boronate_ester', 'cyclic boronate ester':'boronate_ester',
    'boronate_ester':'boronate_ester', 'boroxine':'boroxine', 'triazine':'triazine',
}
RECORD = re.compile(r'(\S+)\s+COFex\.([a-zA-Z0-9_-]+)\.cat(\d+)(?:;.*)?$')


@lru_cache(maxsize=8192)
def precursor(smiles):
    if not isinstance(smiles,str) or not smiles or '.' in smiles:
        return ''
    mol=Chem.MolFromSmiles(smiles)
    if mol is None or any(not a.GetAtomicNum() or a.GetAtomMapNum() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return ''
    return Chem.MolToSmiles(mol)


def length_spread(lengths):
    """Percent of skeletal bonds outside 65% to 135% of the drawing's median."""
    if len(lengths)<3 or not all(math.isfinite(d) for d in lengths):
        return ''
    scale=median(lengths)
    if scale<=0:
        return ''
    return round(100*sum(d<.65*scale or d>1.35*scale for d in lengths)/len(lengths),2)


DRAWN_ELEMENTS={'B','C','N','O','F','Si','P','S','Cl','Br','I'}


def drawn_bond_spread(model):
    if not model.GetNumConformers():
        return ''
    positions=model.GetConformer().GetPositions()
    if any(abs(z)>.001 for z in positions[:,2]):
        return ''
    xy=positions[:,:2]
    return length_spread([math.dist(xy[b.GetBeginAtomIdx()],xy[b.GetEndAtomIdx()]) for b in model.GetBonds()
                         if {b.GetBeginAtom().GetSymbol(),b.GetEndAtom().GetSymbol()}<=DRAWN_ELEMENTS])


def verified_file(record,root,suffix):
    if not isinstance(record,dict) or not isinstance(record.get('path'),str) or not record.get('sha256'):
        return None
    try:
        path=Path(record['path']).resolve()
        if not path.is_relative_to(root) or path.suffix.lower()!=suffix or not path.is_file() or path.stat().st_size>12_000_000:
            return None
        data=path.read_bytes()
    except (OSError, ValueError, RuntimeError):
        return None
    return data if hashlib.sha256(data).hexdigest()==record['sha256'] else None


def reference_spread(patch,root):
    if root is None or not isinstance(patch,dict):
        return ''
    # Stored annotation geometry is independent of the builder. Read only
    # hash-identified local drawing assets and never supply them to prepare.
    records=[patch.get('molfile'),patch.get('mol')]
    if not any(records):
        graph=verified_file(patch.get('graph'),root,'.json')
        if graph:
            try:
                data=json.loads(graph)
            except (UnicodeDecodeError, ValueError):
                return ''
            if isinstance(data,dict):
                records=[data.get('molfile'),data.get('mol')]
    for record in records:
        block=verified_file(record,root,'.mol')
        if not block:
            continue
        lines=block.decode('utf-8-sig',errors='replace').splitlines()
        if len(lines)<5 or 'V2000' not in lines[3] or '3D' in lines[1]:
            continue
        try:
            n,m=int(lines[3][:3]),int(lines[3][3:6])
            atoms=[(float(line[:10]),float(line[10:20]),float(line[20:30]),line[31:34].strip()) for line in lines[4:4+n]]
            bonds=[(int(line[:3])-1,int(line[3:6])-1) for line in lines[4+n:4+n+m]]
            if len(atoms)!=n or len(bonds)!=m or any(abs(x[2])>.001 for x in atoms) or any(min(i,j)<0 or max(i,j)>=n for i,j in bonds):
                continue
            return length_spread([math.dist(atoms[i][:2],atoms[j][:2]) for i,j in bonds if {atoms[i][3],atoms[j][3]}<=DRAWN_ELEMENTS])
        except (ValueError,IndexError):
            continue
    return ''


def inventory_targets(path,reference_root=None):
    with path.open(encoding='utf-8-sig') as handle:
        for line in handle:
            row=json.loads(line)
            if row['task'] not in ('framework','framework_formation'):
                continue
            for index,target in enumerate(row['labels'].get('frameworks') or [{}],1):
                members=target.get('monomers') or []
                feeds=[precursor(m.get('smiles')) if isinstance(m,dict) else '' for m in members]
                complete=bool(feeds) and all(feeds) and all(target.get(k) is not False for k in
                    ('monomer_graphs_complete','monomer_list_complete','precursor_identities_complete'))
                raw=target.get('linkage') or ''
                declared=str(raw.get('type','') if isinstance(raw,dict) else raw).strip()
                label=declared.lower()
                linkage=LINKAGES.get(label) or LINKAGES.get(re.split(r'\s+\(',label,1)[0],'')
                reported=str(target.get('topology') or '').strip()
                match=re.fullmatch(r'(?:[23]d )?([a-z]{3})(?:\s*;.*)?',reported.lower())
                net=match[1] if match else ''
                source=field=''
                conflict=False
                for k in ('cofex_notation','cofex_id','proposed_cofex_id'):
                    value=target.get(k)
                    match=RECORD.fullmatch(value) if isinstance(value,str) else None
                    if match:
                        source,field=value,k
                        other=[precursor(s) for s in match[1].split('.')]
                        conflict=bool(complete and sorted(set(other))!=sorted(set(feeds)))
                        if not conflict and match[2]!='unk' and not net:
                            net=match[2].lower()
                        break
                smiles='.'.join(sorted(set(filter(None,feeds))))
                query=smiles+' COFex.'+(net or 'unk')+'.cat0' if complete and not conflict else ''
                patch=target.get('connected_patch')
                bounded=(patch.get('bounded_smiles') or patch.get('smiles') or '') if isinstance(patch,dict) else ''
                yield dict(annotation=row['id'],target=index,doi=row['doi'],task=row['task'],name=target.get('name') or '',
                           source_linkage=declared,linkage=linkage,source_net=net,reported_topology=reported,
                           reported_catenation=target.get('catenation',''),precursor_smiles=smiles,complete=complete,
                           source_record=source,source_record_field=field,identity_conflict=conflict,
                           single_network_design_query=query,annotated_graph=bounded,
                           annotated_mol_bond_outliers_pct=reference_spread(patch,reference_root))


def same_graph(model,target):
    """Exact labelled isomorphism by joint color refinement and individualization."""
    if model.GetNumAtoms()!=target.GetNumAtoms() or model.GetNumBonds()!=target.GetNumBonds():
        return False
    def graph(mol):
        labels=[(a.GetAtomicNum(),a.GetIsotope(),a.GetFormalCharge(),a.GetTotalNumHs(),
                 a.GetNumRadicalElectrons(),a.GetIsAromatic()) for a in mol.GetAtoms()]
        edges=[{} for _ in labels]
        for bond in mol.GetBonds():
            i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
            edges[i][j]=edges[j][i]=(str(bond.GetBondType()),bond.GetIsAromatic())
        return labels,edges
    left,le=graph(model)
    right,re=graph(target)
    if Counter(left)!=Counter(right):
        return False
    def compare(first,second):
        while True:
            colors={}
            def refine(labels,edges):
                result=[]
                for atom,neighbors in enumerate(edges):
                    signature=(labels[atom],tuple(sorted((kind,labels[j]) for j,kind in neighbors.items())))
                    result.append(colors.setdefault(signature,len(colors)))
                return result
            previous=len(set(first))
            first,second=refine(first,le),refine(second,re)
            counts=Counter(first)
            if counts!=Counter(second):
                return False
            if len(counts)==previous:
                break
        if len(counts)==len(first):
            lookup={color:i for i,color in enumerate(second)}
            mapping=[lookup[color] for color in first]
            return all(re[mapping[i]].get(mapping[j])==kind
                       for i,neighbors in enumerate(le) for j,kind in neighbors.items())
        color=min((n,c) for c,n in counts.items() if n>1)[1]
        first[first.index(color)]=len(counts)
        for j,c in enumerate(second):
            if c==color:
                candidate=second.copy()
                candidate[j]=len(counts)
                if compare(first,candidate):
                    return True
        return False
    return compare(left,right)


@lru_cache(maxsize=128)
def target_graph(smiles,atoms):
    if not isinstance(smiles,str) or not smiles:
        return None
    mol=Chem.MolFromSmiles(smiles,sanitize=False)
    if mol is None:
        return None
    if mol.GetNumAtoms()!=atoms:
        return 'different atom count'
    mol.UpdatePropertyCache(strict=False)
    symmetric_depiction_rings(mol)
    if atoms < 150:
        if Chem.SanitizeMol(mol,catchErrors=True):
            return None
    else:
        # Full aromaticity perception can enumerate exponentially many fused
        # macrocycles. Score the explicit bonds and local five/six-member rings.
        mol.UpdatePropertyCache(strict=True)
        Chem.SetAromaticity(mol,model=Chem.AromaticityModel.AROMATICITY_SIMPLE)
        Chem.AssignRadicals(mol)
    return mol.ToBinary()


def local_aromatic_graph(molecule):
    # RDKit's global aromaticity model can mark an entire small COF pore
    # aromatic, including genuine amide junctions. Apply identical local
    # ring perception to both graphs for a separate representation check.
    mol=Chem.Mol(molecule)
    mol.UpdatePropertyCache(strict=True)
    symmetric_depiction_rings(mol)
    Chem.Kekulize(mol,clearAromaticFlags=True)
    Chem.SetAromaticity(mol,Chem.AromaticityModel.AROMATICITY_SIMPLE)
    Chem.AssignRadicals(mol)
    mol.UpdatePropertyCache(strict=True)
    return mol


def compare_bounded(model, smiles):
    target=target_graph(smiles,model.GetNumAtoms())
    if target is None:
        return 'No parseable bounded target'
    if isinstance(target,bytes):
        if same_graph(model,Chem.Mol(target)):
            return 'Exact bounded connectivity, stereo unscored'
        original=Chem.MolFromSmiles(smiles,sanitize=False)
        if same_graph(local_aromatic_graph(model),local_aromatic_graph(original)):
            return 'Matching connectivity under local aromaticity, stereo unscored'
    return 'Different bounded graph or depiction extent'


@lru_cache(maxsize=128)
def construct(smiles,link,net,star,exchange='',hub_cut=False):
    choices=parse_record(exchange or smiles+' COFex.'+(net or 'unk')+'.cat0',[])
    errors=[]
    for plan in sorted(choices,key=lambda p:'crystal_project' in p):
        try:
            if 'prepare_request' in plan:
                request=plan['prepare_request']
                if not matches_linkage_family(link,request['linkage']):
                    continue
                recipe=prepare_template(**request)
            elif 'assembly_request' in plan:
                request=plan['assembly_request']
                rows=decorated_monomers([dict(**m,id='custom-'+hashlib.sha256(m['smiles'].encode()).hexdigest()[:24],
                    abbreviation=m['label']) for m in request['custom_precursors']])
                by_id={r['id']:r for r in rows}
                ordered=[by_id[k] for k in request['monomer_ids']]
                recipe=assemble(ordered,request.get('linker_sequence') or '',request.get('pore_view'),request.get('pore_layout'))
                recipe['assembly']['precursors']=[dict(smiles=r['smiles'],label=r['abbreviation'],id=r['id']) for r in ordered]
                if not matches_linkage_family(link,recipe['linkage']):
                    continue
            elif 'crystal_project' in plan:
                project=plan['crystal_project']
                if not matches_linkage_family(link,project['linkage']):
                    continue
                from cof_vision.cofdraw_dia import build_dia
                from cof_vision.cofdraw_crystal_figures import junction_patch,chemical_drawing
                a,b=project['monomers']
                cell,review=build_dia(a['smiles'],b['smiles'],project['linkage'],a['label'],b['label'],project['network_count'],net)
                edges=review['periodic_bonds']
                if '<svg' not in chemical_drawing(cell,edges,circuit=junction_patch(cell,edges)):
                    raise ValueError('Molecular projection unavailable')
                record=cell['provenance']['drawing_cofex']
                return record,None,'Periodic design with molecular projection'
            else:
                continue
            recipe.update(schema_version=VERSION,recipe_id='precursor-audit',designed_variant=True)
            recipe['design_template_signature']=design_signature(recipe)
            validate_recipe(recipe)
            if star:
                mode='star' if len(smiles.split('.'))==1 else 'pores'
                mol,_,_,summary,_=repeat_patch(recipe,2,extent_mode=mode)
                if summary['pore_count']!=7:
                    raise ValueError('Seven-pore neighborhood unavailable')
                extent='Hexagon and six adjacent triangles'
                if hub_cut:
                    mol=hub_boundary(mol,recipe)
                    if mol is None:
                        raise ValueError('No valid outer node boundary for this star')
                    extent+=' with the outer arms cut at their hubs'
            else:
                mol,_=build_patch(recipe)
                extent='Default molecular patch'
            return recipe_record(recipe),mol.ToBinary(),extent
        except (ValueError,KeyError,IndexError,RuntimeError) as error:
            errors.append(str(error))
    raise ValueError(errors[0] if errors else 'No matching precursor reaction and net construction.')


def with_status(row,status,detail=''):
    return {**{k:v for k,v in row.items() if k!='annotated_graph'},'status':status,'detail':detail,
            'design_cofex':'','editable_record':'','atoms':'','drawn_bond_outliers_pct':'','tested_extent':'','target_comparison':'','aryl_cut_comparison':'','star_comparison':'','matched_star_extent':''}


def check(row):
    if not row['complete']:
        return with_status(row,'Incomplete or unverified precursors')
    if row['identity_conflict']:
        return with_status(row,'Source COFex and feed identities disagree')
    if not row['linkage']:
        return with_status(row,'Linkage needs review or a construction rule')
    star=row['source_net']=='kgm' and bool(re.search(r'(six|6) triangles',row['reported_topology'],re.I))
    try:
        record,binary,extent=construct(row['precursor_smiles'],row['linkage'],row['source_net'],star)
    except (ValueError,RuntimeError,KeyError,IndexError) as error:
        return with_status(row,'Construction unavailable',str(error))
    result=with_status(row,'Built from precursors')
    result.update(design_cofex=record['compact'],editable_record=record['text'],tested_extent=extent)
    if binary is None:
        result['target_comparison']='Periodic projection, bounded graph unscored'
    else:
        model=Chem.Mol(binary)
        result['atoms']=model.GetNumAtoms()
        result['drawn_bond_outliers_pct']=drawn_bond_spread(model)
        try:
            result['target_comparison']=compare_bounded(model,row['annotated_graph'])
            if result['target_comparison']=='Different bounded graph or depiction extent':
                clipped=aryl_boundary(model)
                if clipped is not None:
                    comparison=compare_bounded(clipped,row['annotated_graph'])
                    if comparison.startswith('Exact'):
                        result['aryl_cut_comparison']='Exact bounded connectivity at aryl cuts, stereo unscored'
                    elif comparison.startswith('Matching'):
                        result['aryl_cut_comparison']=comparison
                if row['source_net']=='kgm' and '.' not in row['precursor_smiles'] and not result['aryl_cut_comparison']:
                    for hub_cut in (False, True):
                        try:
                            _,binary,scope=construct(row['precursor_smiles'],row['linkage'],'kgm',True,hub_cut=hub_cut)
                        except (ValueError, KeyError, IndexError, RuntimeError):
                            continue
                        comparison=compare_bounded(Chem.Mol(binary),row['annotated_graph'])
                        if comparison.startswith(('Exact','Matching')):
                            result.update(star_comparison=comparison,matched_star_extent=scope)
                            break
        except (ValueError,RuntimeError) as error:
            result.update(target_comparison='Bounded comparison unavailable',detail=str(error))
    return result


def limit_worker_memory():
    """Cap native allocations of a single audit worker; retain the parent process."""
    limit = 1536 * 1024**2
    if sys.platform.startswith('linux'):
        import resource
        pages = int(Path('/proc/self/statm').read_text().split()[0])
        cap = pages * os.sysconf('SC_PAGE_SIZE') + limit
        _, hard = resource.getrlimit(resource.RLIMIT_AS)
        if hard != resource.RLIM_INFINITY:
            cap = min(cap, hard)
        resource.setrlimit(resource.RLIMIT_AS, (cap, hard))
    elif os.name == 'nt':
        import ctypes as c
        from ctypes import wintypes as w
        class BasicLimits(c.Structure):
            _fields_ = [('process_time',c.c_int64),('job_time',c.c_int64),('flags',w.DWORD),
                        ('min_ws',c.c_size_t),('max_ws',c.c_size_t),('active',w.DWORD),
                        ('affinity',c.c_size_t),('priority',w.DWORD),('scheduling',w.DWORD)]
        class Limits(c.Structure):
            _fields_ = [('basic',BasicLimits),('io',c.c_uint64*6),('process_memory',c.c_size_t),
                        ('job_memory',c.c_size_t),('peak_process',c.c_size_t),('peak_job',c.c_size_t)]
        kernel = c.WinDLL('kernel32', use_last_error=True)
        kernel.CreateJobObjectW.argtypes = [c.c_void_p, w.LPCWSTR]
        kernel.CreateJobObjectW.restype = w.HANDLE
        kernel.SetInformationJobObject.argtypes = [w.HANDLE,c.c_int,c.c_void_p,w.DWORD]
        kernel.AssignProcessToJobObject.argtypes = [w.HANDLE,w.HANDLE]
        kernel.GetCurrentProcess.restype = w.HANDLE
        job = kernel.CreateJobObjectW(None,None)
        limits = Limits()
        limits.basic.flags = 0x100  # JOB_OBJECT_LIMIT_PROCESS_MEMORY
        limits.process_memory = limit
        if not job or not kernel.SetInformationJobObject(job,9,c.byref(limits),c.sizeof(limits)):
            raise c.WinError(c.get_last_error())
        if not kernel.AssignProcessToJobObject(job,kernel.GetCurrentProcess()):
            raise c.WinError(c.get_last_error())
        # The worker owns this handle until it exits.
        return job


def check_batch(rows,timeout):
    executable=sys.executable
    env=os.environ.copy()
    if os.name=='nt':
        # Killing the venv launcher alone would leave its native child running.
        executable=sys._base_executable
        env['PYTHONPATH']=os.pathsep.join(filter(None,sys.path))
    try:
        proc=subprocess.run([executable,__file__,'--batch'],input=json.dumps(rows).encode(),stdout=subprocess.PIPE,
                            timeout=timeout,env=env,check=False)
        output=proc.stdout
        failure = f'Isolated worker exited with code {proc.returncode}' if proc.returncode else ''
    except subprocess.TimeoutExpired as error:
        output=error.stdout or b''
        failure=''
    completed=[json.loads(line) for line in output.splitlines() if line.endswith(b'}')]
    remaining=rows[len(completed):]
    if remaining and len(rows)==1:
        return [with_status(rows[0],'Worker could not complete' if failure else 'Calculation time limit',
                            failure or f'Exceeded {timeout} seconds for one target')]
    if remaining:
        completed.extend(check_batch(remaining[:1],timeout))
        if remaining[1:]:
            completed.extend(check_batch(remaining[1:],timeout))
    return completed


def main():
    RDLogger.DisableLog('rdApp.*')
    if sys.argv[1:]==['--batch']:
        _job=limit_worker_memory()
        for row in json.load(sys.stdin.buffer):
            print(json.dumps(check(row)),flush=True)
        return
    if sys.argv[1:2]==['draw']:
        parser=argparse.ArgumentParser(description='Export a finite drawing from a precursor COFex record.')
        parser.add_argument('record',help='Compact or complete COFex with precursor SMILES')
        parser.add_argument('output',type=Path,help='Destination MOL file with the drawing coordinates')
        parser.add_argument('--linkage',required=True,help='Reaction rule ID, e.g. imine or cyano_vinylene_sp2_carbon')
        parser.add_argument('--star',action='store_true',help='A Kagome hexagon with its six neighboring triangles')
        cut=parser.add_mutually_exclusive_group()
        cut.add_argument('--aryl-cut',action='store_true',help='Place continuation marks at outer aryl carbons where possible')
        cut.add_argument('--hub-cut',action='store_true',help='Trim outward arms to their braced hubs in a self-condensed Kagome star')
        args=parser.parse_args(sys.argv[2:])
        _job=limit_worker_memory()
        match=RECORD.fullmatch(args.record.strip())
        if not match:
            parser.error('A precursor COFex record is required.')
        if args.hub_cut and not args.star:
            parser.error('--hub-cut requires --star')
        record,binary,extent=construct(match[1],args.linkage,match[2],args.star,args.record,hub_cut=args.hub_cut)
        if not binary:
            parser.error('This construction is a periodic model. Use its 3D workspace for crystal exports.')
        mol=Chem.Mol(binary)
        if args.aryl_cut:
            clipped=aryl_boundary(mol)
            if clipped is not None:
                mol=clipped
        args.output.write_text(Chem.MolToMolBlock(mol,kekulize=False),encoding='utf8')
        print(json.dumps(dict(cofex=record['compact'],extent=extent,atoms=mol.GetNumAtoms(),output=str(args.output))))
        return
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--reference-root',type=Path,help='Directory containing hash-verified annotated MOL/graph files for optional drawing-quality comparison.')
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--timeout',type=int,default=30,help='Seconds before retrying a target in isolation.')
    args=parser.parse_args()
    if args.workers<1 or args.timeout<1:
        parser.error('Workers and timeout must be positive.')
    targets=list(inventory_targets(args.manifest,args.reference_root.resolve() if args.reference_root else None))
    if not targets:
        parser.error('Manifest contains no framework or formation targets.')
    pending=[]
    results={}
    def key(row):
        return row['annotation'],row['target']
    for row in targets:
        if row['complete'] and not row['identity_conflict'] and row['linkage']:
            pending.append(row)
        else:
            results[key(row)]=check(row)
    pending.sort(key=lambda r:(r['precursor_smiles'],r['linkage'],r['source_net']))
    batches=[pending[i:i+24] for i in range(0,len(pending),24)]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for batch in pool.map(lambda items:check_batch(items,args.timeout),batches):
            results.update((key(r),r) for r in batch)
    counts=Counter()
    with args.output.open('w',encoding='utf-8-sig',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(results[key(targets[0])]))
        writer.writeheader()
        for target in targets:
            row=results[key(target)]
            writer.writerow(row)
            counts[row['status']]+=1
    print(json.dumps(dict(annotations=len({t['annotation'] for t in targets}),targets=len(targets),statuses=counts),indent=2))


if __name__=='__main__':
    main()
