"""Original vector schemes for the linkage gallery. Stars represent framework continuations."""
import math
import re
from xml.sax.saxutils import escape

from rdkit import Chem
from rdkit.Chem import rdAbbreviations, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

CHO = '*c1ccc(C=O)cc1'
AMINE = '*c1ccc(N)cc1'
BORON = '*c1ccc(B(O)O)cc1'
CATECHOL = '*c1ccc(O)c(O)c1'
AMINOPHENOL = '*c1ccc(O)c(N)c1'
SALICYL = '*c1cc(O)c(C=O)cc1'
HYDRAZIDE = '*c1ccc(C(=O)NN)cc1'
IMINE = '*c1ccc(C=Nc2ccc(*)cc2)cc1'
HALIDE = '*c1ccc(Br)cc1'

# Example reactants, product and description. Broad families display one specified example.
SCHEMES = {
    'imine': ([CHO, AMINE], IMINE, 'Aryl aldehyde and primary arylamine'),
    'beta-ketoenamine': ([SALICYL, AMINE], '*C1=CC(=O)C(=CNc2ccc(*)cc2)C=C1', 'Condensation followed by keto–enamine tautomerization'),
    'hydrazone': ([CHO, HYDRAZIDE], '*c1ccc(C=NNC(=O)c2ccc(*)cc2)cc1', 'Aldehyde and acylhydrazide example'),
    'azine': ([CHO, 'NN'], '*c1ccc(C=NN=Cc2ccc(*)cc2)cc1', 'Two aldehyde groups condense with hydrazine'),
    'boronate-ester': ([BORON, CATECHOL], '*c1ccc(B2Oc3cc(*)ccc3O2)cc1', 'Boronic acid and vicinal diol'),
    'boroxine': ([BORON], '*B1OB(*)OB(c2ccc(*)cc2)O1', 'Three boronic acid groups form a boroxine ring'),
    'triazine': (['*c1ccc(C#N)cc1'], '*c1ccc(-c2nc(-c3ccc(*)cc3)nc(-c3ccc(*)cc3)n2)cc1', 'Aryl nitrile cyclotrimerization'),
    'c=c': ([CHO, '*c1ccc(CC#N)cc1'], '*c1ccc(C=C(C#N)c2ccc(*)cc2)cc1', 'Aldehyde and arylacetonitrile · cyano-vinylene example'),
    'imide': (['*c1ccc2c(c1)C(=O)OC2=O', AMINE], '*c1ccc2c(c1)C(=O)N(c1ccc(*)cc1)C2=O', 'Cyclic anhydride and primary arylamine'),
    'benzoxazole': ([CHO, AMINOPHENOL], '*c1ccc(-c2nc3cc(*)ccc3o2)cc1', 'Aldehyde and ortho-aminophenol · oxidative cyclization'),
    'benzothiazole': ([CHO, '*c1ccc(S)c(N)c1'], '*c1ccc(-c2nc3cc(*)ccc3s2)cc1', 'Aldehyde and ortho-aminothiophenol · oxidative cyclization'),
    'amide': ([IMINE], '*c1ccc(C(=O)Nc2ccc(*)cc2)cc1', 'Oxidation of an imine to an amide'),
    'reduced-imine': ([IMINE], '*c1ccc(CNc2ccc(*)cc2)cc1', 'Reduction of an imine to a secondary amine'),
    'phenazine': (['*c1ccc2c(c1)C(=O)C(=O)C=C2', '*c1ccc(N)c(N)c1'], '*c1ccc2c(c1)ccc1nc3cc(*)ccc3nc12', 'Vicinal diketone and ortho-diamine'),
    'azo': ([AMINE, '*c1ccc(N=O)cc1'], '*c1ccc(N=Nc2ccc(*)cc2)cc1', 'Arylamine and nitrosoarene coupling'),
    'nitrone': ([CHO, '*c1ccc(NO)cc1'], '*c1ccc(C=[N+]([O-])c2ccc(*)cc2)cc1', 'Aldehyde and N-arylhydroxylamine'),
    'keto-hydrazide': ([SALICYL, HYDRAZIDE], '*C1=CC(=O)C(=CNNC(=O)c2ccc(*)cc2)C=C1', 'Hydroxyaldehyde and acylhydrazide · keto tautomer'),
    'c-n': ([HALIDE, AMINE], '*c1ccc(Nc2ccc(*)cc2)cc1', 'Aryl C–N coupling example'),
    'c-o': ([HALIDE, '*c1ccc(O)cc1'], '*c1ccc(Oc2ccc(*)cc2)cc1', 'Aryl ether coupling example'),
    'c-s': ([HALIDE, '*c1ccc(S)cc1'], '*c1ccc(Sc2ccc(*)cc2)cc1', 'Aryl thioether coupling example'),
    'dioxin': ([CATECHOL, 'N#Cc1c(F)c(F)c(C#N)c(F)c1F'], '*c1ccc2c(c1)Oc1c(C#N)c(F)c(F)c(C#N)c1O2', 'Vicinal diol and tetrafluoroterephthalonitrile'),
    'squaraine': ([AMINE, 'O=C1C(=O)C(O)=C1O'], '*c1ccc([NH+]=C2C(=O)C(Nc3ccc(*)cc3)=C2[O-])cc1', 'Squaric acid and arylamines · one resonance form'),
    'benzobisoxazole': ([CHO, 'Nc1c(O)cc(N)c(O)c1'], '*c1ccc(-c2nc3cc4oc(-c5ccc(*)cc5)nc4cc3o2)cc1', 'Aldehydes and diaminohydroquinone · double oxidative cyclization'),
    'quinoline': (['*c1ccc(N)c(C=O)c1', '*c1ccc(C(=O)C)cc1'], '*c1ccc2nc(-c3ccc(*)cc3)ccc2c1', 'Friedländer annulation of an ortho-aminoaldehyde and an aryl methyl ketone'),
    'imidazole': ([CHO, '*c1ccc(C(=O)C(=O)c2ccc(*)cc2)cc1'], '*c1ccc(-c2nc(-c3ccc(*)cc3)c(-c3ccc(*)cc3)[nH]2)cc1', 'Aldehyde, vicinal diketone and ammonium acetate'),
    'pyridine': ([CHO, '*c1ccc(C(=O)C)cc1'], '*c1ccc(-c2cc(-c3ccc(*)cc3)nc(-c3ccc(*)cc3)c2)cc1', 'Aldehyde and two aryl methyl ketone groups with ammonium acetate'),
}


def fragment(smiles, x, y, width, height, face=0):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'Invalid illustrated fragment {smiles}')
    # Retain the full product graph. Condense conventional end groups on the precursor line.
    if face:
        mol = rdAbbreviations.CondenseMolAbbreviations(mol, rdAbbreviations.ParseAbbreviations(
            'CHO *C=O CHO OHC\nCN *C#N CN NC'), .45)
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    coords = conf.GetPositions()
    ports = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==0 and not a.HasProp('atomLabel')]
    if len(ports)==2 or (len(ports)==1 and face):
        start = coords[ports[0],:2]
        end = coords[ports[1],:2] if len(ports)==2 else coords[:,:2].mean(axis=0)
        angle = -math.atan2(end[1]-start[1],end[0]-start[0]) + (math.pi if face<0 else 0)
        for i,(px,py,_) in enumerate(coords):
            conf.SetAtomPosition(i,(px*math.cos(angle)-py*math.sin(angle),px*math.sin(angle)+py*math.cos(angle),0))
    drawer = rdMolDraw2D.MolDraw2DSVG(width,height)
    opt = drawer.drawOptions()
    opt.useBWAtomPalette()
    opt.clearBackground=False
    opt.bondLineWidth=1.7
    opt.fixedBondLength=17.5
    opt.padding=.055
    for i in ports:
        opt.atomLabels[i]=' '
    rdMolDraw2D.PrepareAndDrawMolecule(drawer,mol)
    drawer.FinishDrawing()
    source=drawer.GetDrawingText()
    body=source[source.index('>', source.index('<svg'))+1:source.rindex('</svg>')]
    bond_ids={mol.GetAtomWithIdx(i).GetBonds()[0].GetIdx() for i in ports}
    def boundary(m):
        tag=m.group()
        if int(m.group(1)) in bond_ids:
            tag=tag.replace("style='","style='stroke-dasharray:3.2,2.5;")
        return tag
    body=re.sub(r"<path class='bond-(\d+) [^>]*>",boundary,body)
    return f'<g transform="translate({x} {y})">{body}</g>'


def illustration(key, name, product):
    reactants,product,note=SCHEMES.get(key,([],product,'Representative product motif'))
    width,height=340,216
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img"><title>{escape(name)}</title>']
    if reactants:
        if len(reactants)==2:
            parts.extend([fragment(reactants[0],8,0,152,72,1),fragment(reactants[1],180,0,152,72,-1)])
            parts.append('<path d="M164 36h12m-6 -6v12 M84 69v5h172v-5" fill="none" stroke="#525552" stroke-width="1"/>')
        else:
            parts.append(fragment(reactants[0],25,0,290,73,1))
        parts.append('<path d="M170 75v22m-4 -6 4 7 4 -7" fill="none" stroke="#343733" stroke-width="1.4"/>')
        if key in ('imidazole','pyridine'):
            parts.append('<text x="180" y="91" font-family="Arial,sans-serif" font-size="12">NH₄OAc</text>')
        parts.append(fragment(product,6,100,328,114))
    else:
        parts.append(fragment(product,12,14,316,183))
    return ''.join(parts)+'</svg>',note
