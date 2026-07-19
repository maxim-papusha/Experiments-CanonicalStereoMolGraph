# Experiments: Canonical StereoMolGraph

Supporting code and notebooks for canonical atom enumeration of
StereoMolGraph / MolGraph structures.

Repository: https://github.com/maxim-papusha/Experiments-CanonicalStereoMolGraph

StereoMolGraph (canon branch): https://github.com/maxim-papusha/StereoMolGraph/tree/canon

## Contents

| Path | Description |
|------|-------------|
| `src/canonical_smg/canon_minimal.py` | Theory-faithful DFS + \(\phi\) enumeration (iterative queue) |
| `notebooks/test_canon_minimal.ipynb` | Shuffle / stability tests for `canon_minimal` |
| `notebooks/test_canon_topo.ipynb` | Tests for the tuned experimental canon in StereoMolGraph |
| `notebooks/atom_invariants.ipynb` | Atom-invariant exploration |
| `data/smiles.txt` | CIP Validation Suite SMILES list |

## Environment

Python 3.11+.

```bash
cd experiments/SI_CanonicalSMG   # or clone this repo
pip install -e .
```

This installs notebook dependencies and **StereoMolGraph from the `canon` branch**:

```text
stereomolgraph @ git+https://github.com/maxim-papusha/StereoMolGraph.git@canon
```

For a local editable StereoMolGraph checkout instead:

```bash
pip install -e ../../packages/StereoMolGraph
pip install -e .
```

(ensure that checkout is on the `canon` branch).

## Quick check

```python
from rdkit import Chem
from stereomolgraph import StereoMolGraph
from canonical_smg import canon_atom_num

g = StereoMolGraph.from_rdmol(
    Chem.AddHs(Chem.MolFromSmiles("C[C@H](O)N")),
    stereo_complete=True,
)
print(dict(sorted(canon_atom_num(g).items(), key=lambda x: x[1])))
```

## License

MIT — see `LICENSE`.
