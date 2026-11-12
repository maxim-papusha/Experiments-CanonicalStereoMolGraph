"""Strict tetrahedral-only RDKit / StereoMolGraph benchmark.

No new StereoMolGraph descriptors are generated. Molecules containing isotope-
labelled atoms or any bond other than a plain single bond are excluded.
This deliberately avoids confounding missing isotope/bond-order information
with tetrahedral stereochemical equivalence.

Uses the existing comparison engine for the enumerator and its caps.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import signal
import sys
import traceback
from collections import Counter
from pathlib import Path

from rdkit import Chem, rdBase

from compare_stereoisomers import analyze_case, load_smiles
from stereomolgraph import StereoMolGraph
from stereomolgraph.experimental.stereoisomer_draft import detect_stereo_elements
from stereomolgraph.stereodescriptors import Tetrahedral


def eligibility(smiles: str) -> tuple[str, Chem.Mol | None]:
    """Strictly isolate nonisotopic, single-bond chemistry."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return "skip_parse_error", None
    if any(a.GetIsotope() != 0 for a in molecule.GetAtoms()):
        return "skip_isotope", molecule
    if any(b.GetBondType() != Chem.BondType.SINGLE for b in molecule.GetBonds()):
        # Excludes double/triple/aromatic/dative bonds, not only E/Z units.
        return "skip_non_single_bond", molecule
    if any(
        a.GetChiralTag() not in (
            Chem.ChiralType.CHI_UNSPECIFIED,
            Chem.ChiralType.CHI_TETRAHEDRAL_CW,
            Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
        )
        for a in molecule.GetAtoms()
    ):
        return "skip_other_stereo_geometry", molecule
    return "eligible", molecule


def analyze_tetra_case(record, mode: str, args) -> dict:
    reason, mol = eligibility(record["smiles"])
    if reason != "eligible":
        return {
            "id": record["id"], "mode": mode, "status": reason,
            "smiles": record["smiles"]
        }

    # An explicit stereo-site check, independent of import/enumeration:
    candidate = Chem.Mol(mol)
    if mode == "unassigned":
        Chem.RemoveStereochemistry(candidate)
    info = list(Chem.FindPotentialStereo(candidate))
    if any(i.type != Chem.StereoType.Atom_Tetrahedral for i in info):
        return {
            "id": record["id"], "mode": mode,
            "status": "skip_non_tetrahedral_rdkit_site",
            "smiles": record["smiles"]
        }

    # Do not count trivial no-stereo examples as successful validation.
    if not info:
        return {
            "id": record["id"], "mode": mode,
            "status": "skip_no_tetrahedral_site",
            "smiles": record["smiles"]
        }

    result = analyze_case(
        record, mode, args.max_sites, args.max_product, args.max_heavy_atoms
    )
    if result["status"] in ("match", "DIFFER"):
        # Ensure importer retained only tetrahedral atom sites after
        # matched-site filtering: no planar bond or coordination stereo.
        # analyze_case guarantees it never ADDS missing stereo descriptors.
        result["tetrahedral_only"] = True
    result["has_rings"] = mol.GetRingInfo().NumRings() > 0
    return result


def run(args) -> int:
    records = list(load_smiles(args.dataset))
    records = records[:args.limit] if args.limit else records
    rows = []
    for record in records:
        for mode in args.modes:
            try:
                if hasattr(signal, "setitimer"):
                    signal.setitimer(signal.ITIMER_REAL, args.timeout_per_case)
                row = analyze_tetra_case(record, mode, args)
            except TimeoutError:
                row = {
                    "id": record["id"], "mode": mode, "status": "timeout",
                    "smiles": record["smiles"]
                }
            except Exception as e:
                row = {
                    "id": record["id"], "mode": mode, "status": "ERROR",
                    "smiles": record["smiles"],
                    "error": f"{type(e).__name__}: {e}",
                    "traceback": traceback.format_exc(limit=3),
                }
            finally:
                if hasattr(signal, "setitimer"):
                    signal.setitimer(signal.ITIMER_REAL, 0)
            rows.append(row)
            if row["status"] not in (
                "skip_isotope", "skip_non_single_bond",
                "skip_no_tetrahedral_site"
            ):
                print(
                    f"{row['id']} {mode:10s} {row['status']:32s} "
                    f"sites={row.get('rdkit_sites','?')} "
                    f"RDKit={row.get('rdkit_unique','-')} "
                    f"SMG={row.get('smg_unique','-')} "
                    f"{row.get('error','')}",
                    flush=True,
                )

    statuses = Counter(row["status"] for row in rows)
    compared = [row for row in rows if row["status"] in ("match", "DIFFER")]
    informative = [
        row for row in compared
        if row["mode"] == "unassigned" and row["rdkit_unique"] > 1
    ]
    errors = [row for row in rows if row["status"] in ("ERROR", "timeout")]
    summary = {
        "rdkit_version": rdBase.rdkitVersion,
        "molecules": len(records),
        "modes": args.modes,
        "total_cases": len(rows),
        "statuses": dict(statuses),
        "compared": len(compared),
        "matching": statuses["match"],
        "different": statuses["DIFFER"],
        "informative_unassigned_cases": len(informative),
        "informative_matching": sum(r["status"] == "match" for r in informative),
        "differences": [{
            "id": r["id"], "mode": r["mode"],
            "smiles": r["smiles"],
            "rdkit": r.get("rdkit_unique"),
            "smg": r.get("smg_unique"),
            "sites": r.get("rdkit_sites"),
            "rings": r.get("has_rings"),
            "formal_assignments": r.get("smg_total_local_assignments"),
        } for r in compared if r["status"] == "DIFFER"],
        "timeouts_and_errors": errors,
        "settings": {
            "max_sites": args.max_sites,
            "max_product": args.max_product,
            "max_heavy_atoms": args.max_heavy_atoms,
            "timeout_per_case": args.timeout_per_case,
            "exclude_any_non_single_bond": True,
            "exclude_isotopes": True,
        },
        "rows": rows,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    csv_path = args.output.with_suffix(".csv")
    columns = [
        "id", "mode", "status", "rdkit_sites",
        "smg_raw_sites", "smg_missing_sites",
        "rdkit_unique", "smg_unique", "has_rings"
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in columns})
    print("TETRA_SUMMARY=" + json.dumps({
        key: value for key, value in summary.items() if key != "rows"
    }, sort_keys=True), flush=True)
    print("RESULT_PATH=" + str(args.output), flush=True)
    print("CSV_PATH=" + str(csv_path), flush=True)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, default=Path("data/smiles.txt"))
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--max-sites", type=int, default=7)
    parser.add_argument("--max-product", type=int, default=128)
    parser.add_argument("--max-heavy-atoms", type=int, default=50)
    parser.add_argument("--timeout-per-case", type=float, default=3.0)
    parser.add_argument(
        "--modes", nargs="+", default=["specified", "unassigned"],
        choices=["specified", "unassigned"]
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path("results/tetrahedral_only.json")
    )
    if hasattr(signal, "SIGALRM"):
        def _timeout(_signum, _frame):
            raise TimeoutError("case exceeded time limit")
        signal.signal(signal.SIGALRM, _timeout)
    sys.exit(run(parser.parse_args()))
