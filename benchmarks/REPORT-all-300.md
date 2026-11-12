# RDKit vs StereoMolGraph: 300-molecule stereoisomer comparison

**Dataset:** [data/smiles.txt](../data/smiles.txt), 300 entries from the CIP Validation Suite.  
**Experiment branch:** stereoisomer-comparison.  
**StereoMolGraph branch:** stereoisomer (experimental baseline; no McKay pruning).  
**RDKit version:** 2026.09.2, in the GitHub Actions runner.  
**Successful full run:** [GitHub Actions #43314861870](https://github.com/maxim-papusha/Experiments-CanonicalStereoMolGraph/actions/runs/43314861870).  
**Reproduction script:** [compare_stereoisomers.py](compare_stereoisomers.py).  
**All 600 case-level outcomes:** [results-300.csv](results-300.csv).

## 1. Controlled comparison: what is measured?

The same source SMILES is used in two modes:

- **specified:** retain stereochemical assignments in the source SMILES,
  enumerate only still-unassigned RDKit stereogenic sites;
- **unassigned:** remove the stereochemistry from a copy of the RDKit
  molecule, then enumerate the resulting unspecified sites.

RDKit uses \`FindPotentialStereo\` and
\`EnumerateStereoisomers(onlyUnassigned=True, unique=True, maxIsomers=0,
tryEmbedding=False)\`.

StereoMolGraph inputs are constructed using
\`StereoMolGraph.from_rdmol(Chem.AddHs(mol), stereo_complete=False,
resonance=False)\`. Its adapter includes many geometry descriptors that
are not independent stereogenic variables. To make a controlled
**same-site enumeration** comparison, input descriptors not corresponding
to RDKit's potential tetrahedral atom or double-bond stereo sites are
removed *before enumeration*. Missing sites are reported as a failure
to construct a comparable input; they are **not added**.

The library's \`detect_stereo_elements\` and
\`enumerate_stereo_assignments(unique=True)\` operate only on descriptors
already present in this filtered input. Candidate generation is
the Cartesian product of the supplied local domains; completed graph
isomorphism is used for deduplication.

**This is not a standalone validation of StereoMolGraph's potential
stereogenic-unit detection.** RDKit defines the matched-site candidate
set by design. Agreement here does not establish chemical correctness
or identity of the **sets** of stereoisomers: the current benchmark
compares counts only.

## 2. Conditions and overall results

The full input contains 300 SMILES. Each is examined in both modes,
for **600 cases**. Individual cases were skipped if:

- there were more than 7 RDKit potential stereo sites;
- there were more than 50 heavy atoms;
- a RDKit-listed site had no corresponding pre-existing StereoMolGraph
  descriptor;
- the per-case computation exceeded 3 seconds.

The generator also has a 128-formal-assignments cap (no cases in this
run were excluded by that condition).

| Status | Case count |
|---|---:|
| Both engines compared | **521** |
| Matching unique-output counts | **463** |
| Different unique-output counts | **58** |
| Missing stereo descriptors in SMG input | 9 |
| More than 7 RDKit sites | 40 |
| More than 50 heavy atoms | 22 |
| Timed out (3 s per case) | 8 |
| **Total** | **600** |

Among comparable cases, **463/521 = 88.9%** matched in count.
However, many specified-mode cases have a trivial 1-versus-1
agreement, so this is not a meaningful global accuracy percentage.

The modes analyzed separately:

| Mode | Compared | Same count | Different count |
|---|---:|---:|---:|
| specified | 266 | 245 | 21 |
| unassigned | 255 | 218 | 37 |

In the unassigned mode, 203 compared cases had more than one distinct
RDKit stereoisomer, and **180 of those 203** matched the SMG count.

The 58 different comparisons represent **41 distinct molecules**.
These are not all the same kind of problem.

### Overcomplete raw adapter descriptors

Among 262 unassigned-mode cases that reached the SMG-import stage:

- 237 had more raw SMG geometry descriptors than potential sites
  reported by RDKit;
- the mean raw descriptor count was **9.68**;
- the mean RDKit potential-site count was **2.15**;
- the highest raw descriptor count was **46**.

This does **not** imply that all excess descriptors are bugs:
some encode planar/ring geometry or hybridization without being
stereoisomer-enumeration variables. It does demonstrate why simply
enumerating every \`parity=None\` descriptor produced by the current
adapter is inappropriate without an explicit candidate-site contract.

## 3. Diagnostic examples and likely causes

### A. Constitutional information missing from SMG graphs

**VS042:** \`CC1=C(C=CC=C1)/N=N/P(=O)([O-])OC\`  
**VS044:** \`COP(=O)(/N=N/C1=CC(=CC=C1)F)[O-]\`

| Molecule | Mode | RDKit | SMG |
|---|---|---:|---:|
| VS042 | specified | 2 | 1 |
| VS042 | unassigned | 4 | 2 |
| VS044 | specified | 2 | 1 |
| VS044 | unassigned | 4 | 2 |

RDKit identifies a stereogenic P center and an E/Z azo bond.
The P substituents include both \`P=O\` and \`P–O^-\` oxygen atoms.
\`rdmol2graph.py\` currently constructs connectivity using atom types
and edges, not the atom charge and bond-order attributes needed
to distinguish those oxygen environments. The stereo-aware graph
isomorphism then can merge formally different P configurations.

This is primarily an issue with the **constitutional graph model used
for isomorphism**, not with the Cartesian-product iterator.

The same representation problem is unambiguously exposed by the
isotope-containing series:

- **VS170, VS175–VS183, VS186–VS187:** in unassigned mode,
  RDKit reports 2 configurations and SMG reports 1.
- **VS119:** RDKit reports 4 and SMG 3 in unassigned mode.

These structures contain isotopic distinctions such as
\`[2H]\`, \`[3H]\`, \`[13CH2]\`, \`[17O]\`, \`[18O]\`,
or isotope-labelled iodine. The SMG importer does not preserve isotope
numbers as atom labels, which introduces artificial symmetry.

**Priority:** preserve the relevant chemical atom attributes (isotope,
formal charge, radical state, etc.) and bond attributes (bond order,
aromaticity, directional constraints as appropriate) in the graph
invariants and isomorphism test. Merely copying attributes without
incorporating them into equality/isomorphism is insufficient.

### B. Ineffective configurations in rigid cyclic frameworks

**VS006:** \`C1[C@H]2C[C@H]3C[C@@H]1C[C@@H](C2)C3\`  
**VS008:** \`C=1[C@@H]2C=C[C@H](C1)C=C2\`

Both modes:

| Molecule | RDKit | SMG |
|---|---:|---:|
| VS006 | 1 | 3 |
| VS008 | 1 | 2 |

With RDKit, these compounds are serialized to SMILES without the
source's tetrahedral stereo tags. RDKit's enumeration algorithm
calls \`Chem.AssignStereochemistry(cleanIt=True, force=True)\` for
each complete candidate and deduplicates isomeric canonical SMILES.
This removes ineffective formal stereotags in these systems.

The SMG algorithm retains and compares explicitly represented
tetrahedral configurations at the candidate sites. It does not
yet implement RDKit-like *stereogenic-unit normalization* in a cyclic
environment.

Because \`tryEmbedding=False\`, the discrepancy occurs even without
3D conformer generation. No claim is made that the RDKit result
alone proves a particular physically accessible configuration count.

Additional cyclic or constrained examples worth investigating:

- VS188–VS190 (constrained exocyclic/cyclic E/Z descriptors);
- VS246–VS248 (coupled alkene and ring stereochemistry);
- VS280 and VS286 (cage systems);
- VS293 and VS296–VS299 (rigid cyclic configurations).

### C. Conditional stereogenicity, symmetry and complex rings

Some nontrivial disagreements require a separate exact
stereo-aware and site-anchored analysis:

| ID | Mode | RDKit | SMG | Context |
|---|---|---:|---:|---|
| VS173 | unassigned | 8 | 16 | Coupled polycyclic and multiple E/Z units |
| VS215 | unassigned | 4 | 2 | Heterocyclic fused polycyclic system |
| VS217 | unassigned | 3 | 2 | Fused polycyclic system |
| VS279 | specified | 1 | 2 | Multiple pseudoasymmetric/related stereocenters |
| VS280 | unassigned | 1 | 8 | Highly symmetric cage |
| VS293 | unassigned | 6 | 8 | Several coupled cyclic E/Z units |
| VS294 | unassigned | 6 | 7 | Polycyclic oxygen-rich system |
| VS295 | unassigned | 6 | 7 | Related oxygen-rich system |

These differences **are not yet proven to be errors in a particular
algorithm**. Some may be a result of missing constitutional labels,
others of cyclic stereo cleanup, and others of distinct symmetry
conventions. Equality of isomer sets must be checked explicitly
before drawing conclusions. The VS294/VS295 pair is especially
interesting for conditional stereogenicity in rings.

### D. Candidate site absent from the SMG importer

**VS014** has sulfoxide sulfur with a lone pair. After removing
stereochemistry, RDKit still identifies an unassigned tetrahedral
stereogenic unit, but the graph adapter does not create an
undefined S stereo descriptor (its untagged tetrahedral heuristic
primarily recognizes four-coordinate atoms).

Other missing-descriptor cases: VS074, VS128, VS129, VS138, VS147, VS237.
There were 9 such **mode-specific comparisons** in the full run.

The enumeration API was explicitly designed to take its complete
candidate set as an **input**: this is therefore an upstream
stereo-site-nomination/adapter limitation, not a reason for the
enumerator to add new descriptors.

## 4. Differences between the two generation algorithms

### RDKit (for this benchmark)

1. Detect potential stereo sites using RDKit stereoperception.
2. Flip the local tetrahedral and bond-stereo states.
3. Reassign/clean stereochemistry after creating each candidate
   (\`AssignStereochemistry(cleanIt=True, force=True)\`).
4. Canonicalize isomeric SMILES and suppress repeated strings.
5. Optionally embed 3D structures; **disabled here**.

### StereoMolGraph experimental branch

1. Use only the already present descriptors (no new sites).
2. For each undefined descriptor, call its \`get_isomers()\` method
   to construct the local configuration domain.
3. Enumerate the Cartesian product of local state choices.
4. Replace states at existing sites and compare complete stereo
   graphs using SMG hashing plus stereo-aware isomorphism.
5. **No** global inactive-site cleanup, geometry-feasibility
   certification or McKay canonical augmentation yet.

This explains why the engines can agree on candidate sites and local
domains yet produce different numbers of distinct outputs.

## 5. Recommended development order

1. **Fix representation fidelity first.** Extend the graph's
   constitutional attributes used in colors and VF2 feasibility to
   include isotopes, formal charges and bond orders (with explicit
   semantics for resonance and aromaticity). Re-run VS042, VS044 and
   VS170–VS187 to ensure lost distinctions are recovered.
2. **Define and test exact inactivity/conditional stereogenicity** for
   completed configurations. Do not greedily delete a site whose
   ligands share initial colors. Start with VS006, VS008 and the
   cyclic series VS188–VS190.
3. **Disentangle structural isomorphism and geometric feasibility**:
   representation-only graph equivalence is not identical to
   physical realizability or configurational stability.
4. **Compare isomer identities rather than only numbers**, using
   a common chemically complete constitution and matching
   stereodescriptors through atom maps; canonical SMILES alone
   does not cover every non-tetrahedral/axial class.
5. **Improve supplied descriptor quality** separately from the
   enumerator; solve lone-pair-containing input cases without
   altering the contract that enumeration never inserts sites.
6. **Only then add group-action/McKay-style pruning.** Exact pruning
   must reproduce the unpruned reference's equivalence classes on
   the agreed representation. Pruning does not repair a faulty
   isomorphism relation or the definition of realizable stereo.

## 6. Reproducibility and caveats

Run in this branch:

\`\`\`bash
python benchmarks/compare_stereoisomers.py \
  --dataset data/smiles.txt \
  --limit 0 \
  --max-sites 7 \
  --max-product 128 \
  --max-heavy-atoms 50 \
  --timeout-per-case 3 \
  --output results/comparison.json
\`\`\`

The workflow installs the live \`StereoMolGraph/stereoisomer\` branch at
execution time, not a pinned commit. For publication-quality work,
pin both the RDKit build and library commit and archive the full JSON.
The included CSV records the **case-level output reported by the
successful run**, while the GitHub Actions run stores the full JSON as
a downloadable artifact.

The CIP Validation Suite was developed primarily to validate
stereochemical descriptors, **not as an independently certified
enumeration reference** for all represented chemistry. RDKit does
not support all helical, axial, planar and coordination stereochemical
classes through the enumeration API. Some "matches" mean that both
systems simply returned the input without recognizing a specialized
stereochemical element.
