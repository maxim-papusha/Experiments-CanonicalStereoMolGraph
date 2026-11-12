"""Compare RDKit and StereoMolGraph stereochemical enumeration on CIP SMILES.

Methodological controls:
- Source SMILES: data/smiles.txt from this repository.
- SMG input sites: ONLY descriptors created by StereoMolGraph.from_rdmol().
- For a *matched-site* comparison, delete parser-created sites NOT listed by
  RDKit's FindPotentialStereo(). Never add a missing descriptor.
- Report missing/extra site counts separately: candidate perception and
  assignment enumeration are different problems.
- Strip all input stereo for "unassigned", keep input stereo for "specified".
- Only count assignments; it does not establish matching molecular structures.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
import traceback
from collections import Counter
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions,
)

from stereomolgraph import StereoMolGraph
from stereomolgraph.experimental.stereoisomer_draft import (
    detect_stereo_elements,
    enumerate_stereo_assignments,
)


def load_smiles(path: Path):
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            continue
        yield {"smiles": fields[0], "id": fields[1], "metadata": fields[2:]}


def stereo_keys(rd_mol: Chem.Mol):
    """Index RDKit potential stereo candidates before adding hydrogens."""
    infos = list(Chem.FindPotentialStereo(rd_mol))
    keys = set()
    unsupported = Counter()
    for info in infos:
        if info.type == Chem.StereoType.Atom_Tetrahedral:
            keys.add(("atom", int(info.centeredOn)))
        elif info.type == Chem.StereoType.Bond_Double:
            bond = rd_mol.GetBondWithIdx(int(info.centeredOn))
            keys.add((
                "bond",
                frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())),
            ))
        else:
            unsupported[str(info.type)] += 1
    return keys, dict(unsupported), infos


def smg_keys(graph: StereoMolGraph):
    return (
        {("atom", int(a)) for a in graph.atom_stereo}
        | {("bond", frozenset(b)) for b in graph.bond_stereo}
    )


def ensure_matched_only_existing(graph: StereoMolGraph, desired: set):
    """Delete extra candidates introduced by import; never add a descriptor."""
    raw = smg_keys(graph)
    for a in tuple(graph.atom_stereo):
        if ("atom", a) not in desired:
            graph.delete_atom_stereo(a)
    for b in tuple(graph.bond_stereo):
        if ("bond", b) not in desired:
            graph.delete_bond_stereo(b)
    return raw


def enum_rdkit_count(mol: Chem.Mol):
    opts = StereoEnumerationOptions(
        onlyUnassigned=True,
        unique=True,
        maxIsomers=0,
        tryEmbedding=False,
    )
    return sum(1 for _ in EnumerateStereoisomers(mol, options=opts))


def analyze_case(record, mode: str, max_sites: int, max_product: int):
    start = time.monotonic()
    result = {
        "id": record["id"],
        "mode": mode,
        "smiles": record["smiles"],
        "metadata": record["metadata"],
        "status": "unknown",
    }
    mol = Chem.MolFromSmiles(record["smiles"])
    if mol is None:
        result["status"] = "rdkit_parse_failure"
        return result

    result["heavy_atoms"] = mol.GetNumHeavyAtoms()
    if mode == "unassigned":
        Chem.RemoveStereochemistry(mol)

    desired, unsupported, info = stereo_keys(mol)
    result["rdkit_sites"] = len(desired)
    result["rdkit_stereo_types_unhandled"] = unsupported
    result["rdkit_specified"] = sum(
        x.specified == Chem.StereoSpecified.Specified for x in info
    )

    if len(desired) > max_sites:
        result["status"] = "skip_too_many_rdkit_sites"
        return result

    # Candidate source is explicitly the existing SMG adapter.
    explicit = Chem.AddHs(Chem.Mol(mol))
    graph = StereoMolGraph.from_rdmol(
        explicit, stereo_complete=False, resonance=False
    )
    raw = ensure_matched_only_existing(graph, desired)
    missing = desired - raw
    extra = raw - desired
    result["smg_raw_sites"] = len(raw)
    result["smg_extra_sites"] = len(extra)
    result["smg_missing_sites"] = len(missing)
    result["smg_missing_keys"] = [str(k) for k in sorted(missing, key=repr)]
    result["smg_extra_keys"] = [str(k) for k in sorted(extra, key=repr)]

    analysis = detect_stereo_elements(graph)
    result["smg_status_counts"] = dict(
        Counter(site.classification for site in analysis.sites)
    )
    result["smg_domain_sizes"] = [len(site.states) for site in analysis.sites]
    result["smg_influence_edges"] = sum(
        len(edges) for edges in analysis.dependencies.values()
    )
    result["smg_total_local_assignments"] = math.prod(
        len(site.states) for site in analysis.sites
    )

    if missing:
        result["status"] = "skip_missing_descriptors"
    elif result["smg_total_local_assignments"] > max_product:
        result["status"] = "skip_too_many_assignments"
    else:
        result["rdkit_unique"] = enum_rdkit_count(mol)
        result["smg_unique"] = sum(
            1 for _ in enumerate_stereo_assignments(graph, analysis, unique=True)
        )
        result["smg_raw_count"] = result["smg_total_local_assignments"]
        result["status"] = (
            "match" if result["rdkit_unique"] == result["smg_unique"]
            else "DIFFER"
        )

    result["duration_s"] = round(time.monotonic() - start, 3)
    return result


def run(args):
    records = list(load_smiles(args.dataset))
    if args.ids:
        ids = set(args.ids.split(","))
        selected = [r for r in records if r["id"] in ids]
    else:
        selected = records[:args.limit] if args.limit else records
    rows = []
    for r in selected:
        for mode in args.modes:
            try:
                row = analyze_case(r, mode, args.max_sites, args.max_product)
            except Exception as e:
                row = {
                    "id": r["id"],
                    "mode": mode,
                    "smiles": r["smiles"],
                    "status": "ERROR",
                    "error": f"{type(e).__name__}: {e}",
                    "traceback": traceback.format_exc(limit=3),
                }
            rows.append(row)
            print(
                f"{row['id']} {mode:10s} {row['status']:28s} "
                f"rd-sites={row.get('rdkit_sites','?')} "
                f"smg-raw={row.get('smg_raw_sites','?')} "
                f"missing={row.get('smg_missing_sites','?')} "
                f"rd-unique={row.get('rdkit_unique','-')} "
                f"smg-unique={row.get('smg_unique','-')} "
                f"error={row.get('error','')}",
                flush=True,
            )
    counts = Counter(r["status"] for r in rows)
    summary = {
        "dataset": str(args.dataset),
        "total_dataset_rows": len(records),
        "selected_molecules": len(selected),
        "row_count": len(rows),
        "rdkit_version": rdBase.rdkitVersion,
        "counts_by_status": dict(counts),
        "compared": counts["match"] + counts["DIFFER"],
        "sampled_ids": [r["id"] for r in selected],
        "filters": {
            "max_sites": args.max_sites,
            "max_product": args.max_product,
            "modes": args.modes,
        },
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("SUMMARY_JSON=" + json.dumps({
        k: v for k, v in summary.items() if k != "rows"
    }, sort_keys=True), flush=True)
    print("DIFFERENCES_JSON=" + json.dumps([
        {k: r.get(k) for k in (
            "id", "mode", "rdkit_unique", "smg_unique", "rdkit_sites",
            "smg_raw_sites", "smg_missing_sites", "smg_extra_sites",
            "smg_total_local_assignments",
        )}
        for r in rows if r["status"] == "DIFFER"
    ]), flush=True)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, default=Path("data/smiles.txt"))
    parser.add_argument("--limit", type=int, default=70)
    parser.add_argument("--ids", default="")
    parser.add_argument("--max-sites", type=int, default=7)
    parser.add_argument("--max-product", type=int, default=256)
    parser.add_argument(
        "--modes", nargs="+", default=["specified", "unassigned"],
        choices=["specified", "unassigned"]
    )
    parser.add_argument(
        "--output", type=Path, default=Path("results/comparison.json")
    )
    sys.exit(run(parser.parse_args()))
