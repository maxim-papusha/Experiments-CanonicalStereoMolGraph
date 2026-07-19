"""Minimal stereo-aware canonical atom enumeration.

Implements the DFS + comparison function φ from the Canonization theory
section (writing/Canonization/2_Theory.tex), with iterative queue-based
backtracking in the style of ``canon.py`` (no recursion):

- Partial map M : atoms → {1, …, |A|}
- Frontier search with a stack of tie-sets (queue)
- φ = (ℓ_a, bonds, stereodescriptors), lexicographically minimized
- Unenumerated atoms use sentinel |A| in M̂
- Stereodescriptor key M̂★(s) = min_{π ∈ P} π · M̂(s)
"""

from __future__ import annotations

import itertools
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from stereomolgraph import AtomId, MolGraph, StereoMolGraph
from stereomolgraph.stereodescriptors import Stereo

CanonNum = int
PhiKey = tuple[Any, ...]
PhiSeq = tuple[PhiKey, ...]


# ---------------------------------------------------------------------------
# Comparison function φ (theory)
# ---------------------------------------------------------------------------


def _m_hat(
    atom: AtomId | None,
    mapping: Mapping[AtomId, CanonNum],
    n_atoms: int,
) -> int:
    """M̂(a): assigned index, or sentinel |A| if unenumerated.

    Lone-pair placeholders (``None`` in an LS) use |A|+1 so they sort after
    every real atom index and the unenumerated sentinel.
    """
    if atom is None:
        return n_atoms + 1
    return mapping[atom] if atom in mapping else n_atoms


def _stereo_atom_tuple(stereo: Stereo) -> tuple[AtomId | None, ...]:
    """Atom order for comparison; parity -1 → inverted geometry."""
    if stereo.parity == -1 and hasattr(stereo, "_inverted_atoms"):
        return tuple(stereo._inverted_atoms())  # type: ignore[attr-defined]
    return tuple(stereo.atoms)


def _permutation_index_sets(stereo: Stereo) -> Sequence[tuple[int, ...]]:
    """Index permutations defining P (or all perms if parity is undefined)."""
    n = len(stereo.atoms)
    if stereo.parity is None:
        return tuple(itertools.permutations(range(n)))
    return tuple(tuple(p) for p in stereo.PERMUTATION_GROUP)


def _m_star(
    stereo: Stereo,
    mapping: Mapping[AtomId, CanonNum],
    n_atoms: int,
) -> tuple[int, ...]:
    """M̂★(s) = min_{π ∈ P} π · M̂(s)."""
    atoms = _stereo_atom_tuple(stereo)
    filled = tuple(_m_hat(a, mapping, n_atoms) for a in atoms)
    best: tuple[int, ...] | None = None
    for perm in _permutation_index_sets(stereo):
        image = tuple(filled[i] for i in perm)
        if best is None or image < best:
            best = image
    assert best is not None
    return best


def _ell_S(stereo: Stereo) -> tuple[Any, ...]:
    """Stereodescriptor label ℓ_S(s): type only (orientation is in M̂★)."""
    return (stereo.__class__.__name__,)


def _ell_b(_a: AtomId, _b: AtomId) -> tuple[Any, ...]:
    """Bond label ℓ_b. Constant for plain MolGraph / SMG (no reaction tags)."""
    return ()


def phi(
    atom: AtomId,
    *,
    mapping: Mapping[AtomId, CanonNum],
    atom_labels: Mapping[AtomId, Any],
    neighbors: Mapping[AtomId, frozenset[AtomId] | set[AtomId]],
    stereo_of_atom: Mapping[AtomId, Sequence[Stereo]],
    n_atoms: int,
) -> PhiKey:
    """Comparison key φ(G_S, M, a) from the theory section."""
    atom_key = atom_labels[atom]

    bond_key = tuple(
        sorted(
            (_ell_b(atom, nbr), _m_hat(nbr, mapping, n_atoms))
            for nbr in neighbors[atom]
        )
    )

    stereo_key = tuple(
        sorted(
            (_ell_S(s), _m_star(s, mapping, n_atoms))
            for s in stereo_of_atom.get(atom, ())
        )
    )

    return (atom_key, bond_key, stereo_key)


def _worse_than(seq_a: Sequence[PhiKey], seq_b: Sequence[PhiKey]) -> bool:
    """True if seq_a is lexicographically worse (larger) than seq_b.

    Empty ``seq_b`` means no best yet, so nothing is worse than it.
    Compared only on the shared prefix length (like ``canon.py``).
    """
    if not seq_b:
        return False
    return tuple(seq_a) > tuple(seq_b)


# ---------------------------------------------------------------------------
# Iterative search state (queue / backtracking, as in canon.py)
# ---------------------------------------------------------------------------


@dataclass
class _Parameters:
    """Immutable graph data for the search."""

    atoms: frozenset[AtomId]
    n_atoms: int
    nbrhds: Mapping[AtomId, frozenset[AtomId] | set[AtomId]]
    atom_labels: Mapping[AtomId, Any]
    stereo_of_atom: Mapping[AtomId, list[Stereo]]


@dataclass
class _CanonState:
    """Mutable search state.

    ``queue`` is a stack of remaining candidate sets at each depth (FILO).
    Python dicts preserve insertion order; ``popitem()`` undoes the last map.
    """

    best_labels: list[PhiKey] = field(default_factory=list)
    best_mapping: dict[AtomId, CanonNum] = field(default_factory=dict)

    current_labels: list[PhiKey] = field(default_factory=list)
    current_mapping: dict[AtomId, CanonNum] = field(default_factory=dict)

    frontier: set[AtomId] = field(default_factory=set)
    queue: list[set[AtomId]] = field(default_factory=list)


def _candidate_phi(atom_id: AtomId, params: _Parameters, state: _CanonState) -> PhiKey:
    return phi(
        atom_id,
        mapping=state.current_mapping,
        atom_labels=params.atom_labels,
        neighbors=params.nbrhds,
        stereo_of_atom=params.stereo_of_atom,
        n_atoms=params.n_atoms,
    )


def _update_queue(params: _Parameters, state: _CanonState) -> None:
    """Score the frontier, push the min-φ tie set, prune dominated branches.

    An empty tie set means the branch is finished (or pruned). A complete
    mapping with a better Φ updates ``best_*``.
    """
    if state.frontier:
        frontier = state.frontier
    else:
        # New connected component, or first atom
        frontier = set(params.atoms) - state.current_mapping.keys()

    candidate_labels = {
        atom_id: _candidate_phi(atom_id, params, state) for atom_id in frontier
    }

    atoms: set[AtomId] = set()
    best_step_label: PhiKey = ()

    if candidate_labels:
        best_step_label = min(candidate_labels.values())
        potential_new_labels = state.current_labels + [best_step_label]

        # Prune if this prefix is already worse than the best complete Φ
        if not _worse_than(
            potential_new_labels,
            state.best_labels[: len(potential_new_labels)],
        ):
            atoms = {
                atom_id
                for atom_id, lab in candidate_labels.items()
                if lab == best_step_label
            }

    state.queue.append(atoms)
    state.current_labels.append(best_step_label)

    if (
        not atoms
        and params.n_atoms == len(state.current_mapping)
        and (
            not state.best_labels
            or not _worse_than(state.current_labels, state.best_labels)
        )
    ):
        # Prefer strictly better Φ; on equal Φ keep the first found
        if not state.best_labels or tuple(state.current_labels) < tuple(
            state.best_labels
        ):
            state.best_labels = state.current_labels.copy()
            state.best_mapping = state.current_mapping.copy()


def _initialize(g: MolGraph) -> _Parameters:
    atom_labels = {a: t for a, t in zip(g.atoms, g.atom_types)}
    stereo_of_atom: dict[AtomId, list[Stereo]] = defaultdict(list)
    if isinstance(g, StereoMolGraph):
        for s in g.stereo.values():
            for atom in s.atoms:
                if atom is not None:
                    stereo_of_atom[atom].append(s)

    return _Parameters(
        atoms=frozenset(g.atoms),
        n_atoms=len(g.atoms),
        nbrhds=g.neighbors,
        atom_labels=atom_labels,
        stereo_of_atom=stereo_of_atom,
    )


def _run(g: MolGraph) -> tuple[dict[AtomId, CanonNum], list[PhiKey]]:
    """Iterative DFS with queue of tie-sets (same control flow as ``canon.py``)."""
    params = _initialize(g)
    if params.n_atoms == 0:
        return {}, []

    state = _CanonState()
    # next enumeration index is 1, 2, …, |A| (theory)
    next_num: CanonNum = 1

    _update_queue(params, state)

    while state.queue:
        # ---- step forward ----
        if state.queue[-1]:
            next_atom = state.queue[-1].pop()

            state.current_mapping[next_atom] = next_num
            next_num += 1

            # incremental frontier update
            state.frontier.discard(next_atom)
            state.frontier |= set(params.nbrhds[next_atom]) - state.current_mapping.keys()

            _update_queue(params, state)

        # ---- step back ----
        else:
            state.queue.pop()
            state.current_labels.pop()
            if state.queue:
                atom_id, _ = state.current_mapping.popitem()
                next_num -= 1

                not_frontier_anymore = {
                    nbr
                    for nbr in params.nbrhds[atom_id]
                    if not any(
                        nbr2 in state.current_mapping for nbr2 in params.nbrhds[nbr]
                    )
                }
                state.frontier -= not_frontier_anymore
                if set(params.nbrhds[atom_id]).intersection(state.current_mapping.keys()):
                    state.frontier.add(atom_id)

    return state.best_mapping, state.best_labels


def canon_atom_num(g: MolGraph) -> Mapping[AtomId, CanonNum]:
    """Return a canonical atom enumeration M_C : A → {1, …, |A|}.

    Independent of the input atom order for isomorphic labeled graphs under
    the minimal φ of the theory section.
    """
    mapping, _labels = _run(g)
    return mapping


def canon_phi_sequence(g: MolGraph) -> PhiSeq:
    """Return the lexicographically minimal comparison sequence Φ."""
    _mapping, labels = _run(g)
    return tuple(labels)
