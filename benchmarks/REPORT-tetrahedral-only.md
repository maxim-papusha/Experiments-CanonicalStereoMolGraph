# Tetrahedral-only stereoisomer enumeration: focused benchmark

**Status:** completed baseline comparison; no fixes to the enumeration algorithm
are implied by count disagreements alone.

- Input: [300 CIP Validation Suite molecules](../data/smiles.txt)
- RDKit: **2026.09.2**
- StereoMolGraph: experimental \`stereoisomer\` branch
- Code: [compare_tetrahedral_only.py](compare_tetrahedral_only.py)
- [Full successful workflow run](https://github.com/maxim-papusha/Experiments-CanonicalStereoMolGraph/actions/runs/43316529675)
- [Machine-readable 15 discrepant cases](tetrahedral-discrepancies.csv)

## Scope

The filter keeps **only isotope-free molecules whose bonds are all single**.
Aromatic, double, triple and dative bonds are excluded, including carbonyl
and aromatic rings. This is a stricter restriction than excluding E/Z alone,
chosen to eliminate constitutional-bond-order confounders from the current
StereoMolGraph representation. It does **not** exclude saturated rings and
bridged polycyclic systems.

All eligible stereo sites must be of RDKit type
\`Chem.StereoType.Atom_Tetrahedral\`.
StereoMolGraph is allowed to enumerate **only descriptors already on
the input graph**. To isolate the enumeration from automatic candidate
nomination, the benchmark filters the adapter-generated descriptor
sites to precisely those already found by RDKit; it never creates
a new descriptor. If an expected site is missing, that case is excluded.

Two modes are compared:
1. \`specified\`: preserve stereochemical assignments in the source SMILES
   *as RDKit interprets them* (RDKit can clean ineffective tags).
2. \`unassigned\`: strip stereo tags from the input then enumerate
   possible tetrahedral configurations.

Each mode compares \`EnumerateStereoisomers(unique=True)\` against
\`enumerate_stereo_assignments(unique=True)\`. Neither uses 3D
embedding or McKay-style automorphism pruning. Counts are not
independent proof of distinct chemical stereoisomers.

## Results

Two attempts per source molecule: **600 mode-specific cases**.

| Outcome | Number of cases |
|---|---:|
| Completed comparison | 152 |
| Matching counts | **137** |
| Different counts | **15** |
| Excluded: isotope label | 32 |
| Excluded: aromatic, double, triple or other non-single bond | 352 |
| Excluded: more than 50 heavy atoms | 16 |
| Excluded: more than 7 candidate sites | 36 |
| Excluded: no RDKit tetrahedral site | 3 |
| Excluded: missing descriptor in SMG input | 4 |
| Per-case timeout (3 s) | 5 |
| Total | 600 |

**All 15 different-count cases belong to 10 ring-containing
molecules.** No acyclic molecule in the successfully compared subset
showed a count mismatch.

In the unassigned mode, there were **69 comparisons producing more than
one RDKit stereoisomer**; 67 of these matched the StereoMolGraph count.
This is a useful targeted result, but a larger test set and exact
cross-structure matching are still required.

### All discrepant tetrahedral-only molecules

| VS ID | Mode(s) differing | RDKit count | SMG count | Interpretation |
|---|---|---:|---:|---|
| VS006 | both | 1 | 3 | Rigid bridged cage, dependent configurations |
| VS207 | specified | 1 | 2 | Multiple ring stereosites |
| VS279 | specified | 1 | 2 | Symmetric rings and possible pseudoasymmetry |
| VS280 | both | 1 | 8 | Highly constrained, symmetric cage |
| VS286 | specified | 1 | 2 | Spiro/bridged cyclic constraints |
| VS294 | unassigned | 6 | 7 | Oxygen-containing bridged cyclic system |
| VS295 | unassigned | 6 | 7 | Related bridged cyclic system |
| VS296 | both | 1 | 2 | Rigid cyclic geometry |
| VS297 | both | 1 | 2 | Rigid cyclic geometry |
| VS298 | both | 1 | 2 | Bicyclic tetrahedral configuration |

The \`specified\` and \`unassigned\` modes may have different
formal search domains. Detailed values for every discrepant
mode are in [tetrahedral-discrepancies.csv](tetrahedral-discrepancies.csv).

The SMG output is a graph-theoretic set of independent **local
tetrahedral descriptor assignments**, filtered by its existing
stereo-graph isomorphism. RDKit assigns/cleans tetrahedral stereo
based on its molecular stereoperception and canonical isomeric
SMILES. Rigid-ring constraints, inactive descriptors and conditional
stereogenicity can therefore produce different counts, without
either count alone establishing the correct physical isomer count.

## What this benchmark tells us about the next algorithm

The next useful research target is **stereogenicity under fixed and
partial assignments of tetrahedral descriptors in cyclic graphs**.

For an existing site $s$, let $D_s$ denote its two local tetrahedral
configurations. At a complete assignment $\sigma_{\neg s}$ of all
other sites, compare

$$
G(\sigma_{\neg s},s=+1)\stackrel{?}{\cong}_s
G(\sigma_{\neg s},s=-1).
$$

The rooted isomorphism fixes the central atom of $s$ while preserving
already specified stereo elsewhere. The site is stereogenic in that
context only if the two rooted graphs differ. An initial equality
of refinement colors must **not** be used to permanently discard the
site, because remote assignments can activate it.

The priority sequence is:

1. **VS006:** smallest particularly useful cyclic regression example;
   inspect which local tetrahedral configurations are actually
   distinguishable and how RDKit interprets its original tags.
2. **VS294 and VS295:** test conditional stereogenicity and exact
   graph equivalence; the 6-versus-7 disagreement is not explained
   by isotope, bond order or E/Z because all have been filtered.
3. **VS296–VS298:** test whether cyclic constraints couple local
   tetrahedral configurations.
4. **VS280:** advanced highly symmetric-cage case, useful once
   simpler cases work.

Do not add new sites inside the stereoisomer generator.
Any detection or cleaning must instead classify the **existing**
descriptors as active, inactive or conditional and only transform
configurations already represented by those sites.

### Further caution

The present benchmark deliberately uses RDKit's potential
stereocenter identification as a matching-site input filter and so
**does not benchmark independent tetrahedral site detection**.
It also measures counts, not identity of all enumerated isomers.
A formally correct stereo-isomer set will require comparing the
completed structures under a common representation, and chemical
feasibility for highly constrained bridges may require separate
geometry checks. No McKay-style orbit pruning should be introduced
before the unpruned reference is validated on these systems.
