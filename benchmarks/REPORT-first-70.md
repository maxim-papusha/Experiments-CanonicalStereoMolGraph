# Benchmark: RDKit vs StereoMolGraph stereoisomer enumeration

**Status:** preliminary, 70/300 CIP Validation Suite molecules.  
**Source:** [data/smiles.txt](../data/smiles.txt), 300 records.  
**Baseline script:** [benchmarks/compare_stereoisomers.py](compare_stereoisomers.py).  
**First run:** https://github.com/maxim-papusha/Experiments-CanonicalStereoMolGraph/actions/runs/43314294057

## Methods

The same SMILES are used in both engines. For each molecule, analyze both
`specified` (retain original stereochemistry) and `unassigned` (call
`Chem.RemoveStereochemistry` first) modes.

**RDKit:** `Chem.FindPotentialStereo`, then
`EnumerateStereoisomers(onlyUnassigned=True, unique=True, maxIsomers=0,
tryEmbedding=False)`. RDKit 2026.09.2 in the GitHub Actions runner.

**StereoMolGraph:** explicitly add hydrogens before using
`StereoMolGraph.from_rdmol(stereo_complete=False,resonance=False)`. The
adapter creates many **possible geometry descriptors** for atoms/bonds.
To isolate enumeration from site-nomination differences, remove existing
descriptors not in the RDKit `FindPotentialStereo` site list; never add a
missing descriptor. The experimental branch's `detect_stereo_elements`
and `enumerate_stereo_assignments(unique=True)` are then called.

**Crucial caveat:** This is a **matched-site enumeration benchmark**,
*not* an independent comparison of potential stereocenter detection.
It compares counts only; matching counts do not prove matching molecular
structures or physical validity. RDKit additionally performs
`Chem.AssignStereochemistry(cleanIt=True, force=True)` for each candidate
and deduplicates via canonical isomeric SMILES. StereoMolGraph deduplicates
via its own graph hash and stereo-aware isomorphism.

## Summary (first 70 molecules, two modes)

| Outcome | Case count |
|---|---:|
| Compared | 135 |
| Equal unique-isomer counts | 127 |
| Different counts | 8 |
| Skipped: missing SMG descriptor | 1 |
| Skipped: >7 RDKit potential sites | 4 |
| Total considered | 140 |

All 8 unequal comparisons came from **4 source molecules**, each in both
the specified and unassigned modes.

### Explicit discrepancies (unassigned mode)

| ID | RDKit | StereoMolGraph | Main explanation |
|---|---:|---:|---|
| VS006 | 1 | 3 | Bridgehead stereo in a symmetric cage is not a set of freely independent configurations |
| VS008 | 1 | 2 | Similarly constrained polycyclic/bridged case; RDKit cleans ineffective stereotags |
| VS042 | 4 | 2 | P center with constitutionally distinct oxo and negatively charged O substituents becomes spuriously symmetric in SMG |
| VS044 | 4 | 2 | Same P-center issue |

In `specified` mode, VS006 and VS008 remain (1 vs 3) and (1 vs 2);
VS042 and VS044 become (2 vs 1).

### Why VS042 and VS044 lose stereoisomers

VS042: `CC1=C(C=CC=C1)/N=N/P(=O)([O-])OC`  
VS044: `COP(=O)(/N=N/C1=CC(=CC=C1)F)[O-]`

RDKit recognizes a tetrahedral phosphorus center plus an azo E/Z bond.
After removing specified stereo: 2 x 2 = 4 variants.

The current SMG `rdmol2graph.py` writes only atom symbols and connectivity
for the molecular graph, not formal charge, isotope, or bond-order labels.
The `P=O` oxygen and `P-O^-` oxygen are thus indistinguishable as
connected O atoms. In the graph isomorphism used for deduplication,
switching the P configuration can be treated as an equivalent
permutation. Consequently, **2** outputs are retained instead of **4**.

This is **not primarily an enumeration-loop problem**. The graph's
chemical attributes and the equality/isomorphism invariants need to
distinguish constitutionally distinct atoms and bonds.

### Why VS006 and VS008 overenumerate

VS006: `C1[C@H]2C[C@H]3C[C@@H]1C[C@@H](C2)C3`  
VS008: `C=1[C@@H]2C=C[C@H](C1)C=C2`

RDKit's parsing and stereo cleanup yield a nonstereo canonical SMILES
for these heavily constrained ring systems. `FindPotentialStereo` may
still report individual candidate sites, but
`EnumerateStereoisomers` cleans stereochemical tags when producing
the completed structures. In SMG, locally represented tetrahedral
orientations are combined without a global ring-stereo feasibility
normalization. This produces **3** and **2** formal graph configurations,
respectively.

Note that `tryEmbedding=False` was used; the difference therefore
already appears during RDKit stereochemical perception/cleanup and
canonicalization, without 3D embedding. It should not be described
exclusively as an embedding/energy restriction.

### Candidate detection mismatch: VS014

VS014: `[O-][N+](C=1C=CC(=CC1)[S@@](OCC)=O)=O`

After removing stereochemistry RDKit finds the sulfoxide S center;
StereoMolGraph's adapter does not create the corresponding undefined
descriptor, because its untagged inference path primarily recognizes
four-coordinate SP3 atoms and the sulfur has three explicit neighbors
plus a lone pair. The benchmark therefore **skips** the comparison
rather than secretly creating a site.

This is an upstream candidate-nomination issue, not a deficiency of the
input-descriptor-only enumeration algorithm.

### Descriptor overgeneration at ingestion

For 68 unassigned-mode inputs evaluated far enough to construct SMG:
- **66** had *more* raw SMG descriptors than RDKit potential stereo sites.
- Average raw SMG descriptor count was **10.31**, versus **1.76**
  RDKit candidate sites in the same 68 inputs.
- The maximum raw SMG descriptor count in this subset was **46**.

Thus, using every raw descriptor from `from_rdmol` as an independent
stereogenic site would dramatically overgenerate assignments. The
benchmark deliberately filters raw descriptors to isolate the
enumeration comparison.

## Prioritized remedies

1. **Preserve constitutional labels** (bond order, formal charge,
   isotope and relevant atom attributes) in SMG import, color refinement
   and graph-isomorphism feasibility before trusting stereoisomer counts.
2. **Specify an explicit candidate input contract.** An unknown parity
   means an eligible geometry site, not necessarily a proven stereocenter.
3. **Add stereo cleanup/normalization for cyclic constraints and inactive
   sites**, with exact tests. Do not confuse graph-isomorphism differences
   with physically distinct stereoisomers.
4. **Test complete isomer sets** on matched chemical representations,
   not just their cardinalities.
5. Keep an exhaustive unpruned assignment generator as the reference
   oracle before introducing McKay-style automorphism pruning.

## Scope limits

The CIP validation dataset contains more stereochemical types than
either of these automatic enumeration pipelines necessarily supports,
including helical, atropisomeric and other specialized configurations.
No general correctness claim follows from agreement with RDKit alone.
