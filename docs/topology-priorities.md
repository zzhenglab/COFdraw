# COF topology priorities for the monomer picker

Among 1,242 records in a CoRE COF survey, 770 of the 1,099 two-dimensional entries were assigned hcb and 260 were assigned sql. Of 143 three-dimensional entries, 79 were assigned dia. These are counts of structure records. A material with multiple reported stackings or interpenetration states can contribute more than once. The counts below come from the authors’ reported-materials column and cover selected nets. [Published survey](https://doi.org/10.1021/acsomega.5c10086) · [Authors’ tabulation, Tables 1 and 2](https://doi.org/10.26434/chemrxiv-2025-085wt)

“Connections” below means directions in which the reacted unit joins the framework. A linear two-connected spacer can extend an edge. Condensation can also create a new junction, such as a three-connected boroxine ring. For each pair the picker needs both a supported linkage and compatible directions at the participating sites.

| Net | Dimension | Records | Representative node and partner geometry | Monomer picker |
| --- | --- | ---: | --- | --- |
| hcb | 2D | 770 | Planar trigonal 3 + linear 2, or two planar trigonal 3 units | Available |
| sql | 2D | 260 | Planar four-arm node + linear 2 or another planar four-arm node. Square or rhombic pores | Available |
| dia | 3D | 79 | Tetrahedral 4 + linear 2. Two complementary tetrahedral 4 units can also make a diamond net | Molecular 4 + 2 cage. Periodic builder supports qualified 4 + 2 and 4 + 4 pairs |
| kgm | 2D | 41 | Planar 4 + linear 2. The four arms ideally form alternating 60° and 120° angles around a vertex | Available alongside sql |
| ctn | 3D | 8 | Tetrahedral 4 + trigonal planar 3 | Available as a finite circuit drawing |
| pts | 3D | 7 | Tetrahedral 4 + square-planar 4 | Available |
| stp | 3D | 7 | Triangular-prismatic 6 + square-planar 4 | Available as a finite circuit drawing |
| fes | 2D | 6 | Three-connected vertices shared by a four-membered and two eight-membered net rings. Mapping a real monomer to these junctions needs a worked example | Candidate for later |
| bor | 3D | 5 | Tetrahedral 4 + trigonal planar 3, joined with different periodic connectivity from ctn | Available as a finite circuit drawing |
| hxl | 2D | 5 | Planar hexagonal 6 + linear 2 | Available |
| kgd | 2D | 5 | Planar hexagonal 6 + planar trigonal 3 | Available as a finite sheet drawing |
| bex | 2D | 5 | Reported route uses triangular tritopic and rectangular tetratopic monomers with partial site occupation | Model site occupation explicitly |
| ffc | 3D | 5 | Trigonal 3 + planar four-arm units in a spatial arrangement with two inequivalent four-connected sites | Candidate for later |
| tbo | 3D | 4 | Trigonal 3 + square-planar 4, connected across different planes | Available as a finite circuit drawing |

## Original priorities and current coverage

1. **ctn and bor.** Both are implemented for supported linkages and geometries. COF-105 and COF-108 give reference structures with the same local geometry but different nets. Some original examples use boronate chemistry, which also requires a molecular linkage rule. [El-Kaderi et al., Science](https://yaghi.berkeley.edu/pdfPublications/2007designSynFull.pdf)
2. **dia and lon.** The periodic builder supports both backbones, including eligible pairs of different tetrahedral monomers. LZU-111 is a reference for the alternative lon backbone. Show one network first and treat interpenetration as a separate structural choice. [Ma et al., Science, Figs. 3–4](https://doi.org/10.1126/science.aat7679)
3. **kgd and stp.** Both are implemented. The six reactive directions surround one plane in kgd and occupy two triangular faces of a prism in stp. COF-360 and TUS-64 provide molecular examples for the two drawings. [COF-360](https://doi.org/10.1021/jacs.9b05626) · [TUS-64](https://onlinelibrary.wiley.com/doi/full/10.1002/ange.202300172)
4. **bex and ffc.** These remain candidates for partial occupation and differently oriented planar units. The picker already offers tbo. A 3 + 4 reagent pair alone cannot choose the dimension. [Sub-stoichiometric bex COFs](https://doi.org/10.1038/s41467-019-10574-6) · [Interpenetrated ffc COFs](https://doi.org/10.1039/C9SC04882K)

As a cross-check, counting the 315 distinct COF IDs in the published deCOFpose plotting CSV gives 232 hcb, 39 sql and 30 dia entries. Within dia, the table has 28 records with a 2 + 4 fragment pair and two with 4 + 4. This smaller dataset selects for successful fragmentation and electronic-structure analysis, so the counts above use the larger survey. [Dataset and CSV](https://zenodo.org/records/14855145) · [deCOFpose study](https://doi.org/10.1039/D5QM00727E)
