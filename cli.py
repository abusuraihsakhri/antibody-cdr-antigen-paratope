#!/usr/bin/env python3
"""
Command-Line Interface for Antibody CDR-H3 Modeling, Paratope Contacts & Developability.

Domain: Computational Structural Biology & Immunoinformatics
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any, Dict, List, Optional, Tuple

import antibody_cdr_antigen_paratope as acap


def format_table(headers: List[str], rows: List[List[Any]]) -> str:
    """Format a clean ascii/unicode table for terminal output."""
    if not rows:
        return "(no records)"
    str_rows = [[str(cell) for cell in row] for row in rows]
    col_widths = [len(h) for h in headers]
    for row in str_rows:
        for i, cell in enumerate(row):
            if i < len(col_widths):
                col_widths[i] = max(col_widths[i], len(cell))
            else:
                col_widths.append(len(cell))

    header_line = " | ".join(h.ljust(col_widths[i]) for i, h in enumerate(headers))
    sep_line = "-+-".join("-" * col_widths[i] for i in range(len(col_widths)))
    data_lines = [
        " | ".join(cell.ljust(col_widths[i]) for i, cell in enumerate(row))
        for row in str_rows
    ]
    return f"{header_line}\n{sep_line}\n" + "\n".join(data_lines)


def cmd_scan(args: argparse.Namespace) -> None:
    if args.json_data:
        try:
            raw_data = json.loads(args.json_data)
        except json.JSONDecodeError as e:
            print(f"Error parsing JSON data: {e}", file=sys.stderr)
            sys.exit(1)
        # Expect list of dicts with {"chain": "H", "id": 33, "aa": "Y", "bsa": 95.2}
        residue_bsa = {}
        for item in raw_data:
            key = (item["chain"], int(item["id"]), item["aa"])
            residue_bsa[key] = float(item["bsa"])
    elif args.file:
        with open(args.file, "r", encoding="utf-8") as f:
            raw_data = json.load(f)
        residue_bsa = {}
        for item in raw_data:
            key = (item["chain"], int(item["id"]), item["aa"])
            residue_bsa[key] = float(item["bsa"])
    else:
        # Default sample dataset
        residue_bsa = {
            ("H", 31, "S"): 78.5,
            ("H", 33, "Y"): 95.2,
            ("H", 52, "W"): 30.4,
            ("H", 100, "D"): 115.0,
            ("H", 101, "Y"): 88.9,
            ("H", 102, "F"): 66.3,
            ("L", 32, "Y"): 55.0,
            ("L", 91, "S"): 20.0,
            ("L", 96, "W"): 70.4,
        }

    results = acap.perform_alanine_scan(residue_bsa)

    if args.json:
        out = [
            {
                "chain": r.chain,
                "residue_id": r.residue_id,
                "residue_name": r.residue_name,
                "bsa_a2": r.bsa_a2,
                "ddg_kcal_mol": r.ddg_kcal_mol,
                "classification": r.classification,
            }
            for r in results
        ]
        print(json.dumps(out, indent=2))
        return

    print("=" * 75)
    print("  PARATOPE ALANINE SCANNING & HOT-SPOT PREDICTION")
    print("  Threshold: ddG >= 1.5 kcal/mol (Clackson & Wells)")
    print("=" * 75)
    headers = ["Chain", "Pos", "Res", "BSA (A^2)", "ddG (kcal/mol)", "Classification"]
    rows = [
        [r.chain, r.residue_id, r.residue_name, f"{r.bsa_a2:.1f}", f"{r.ddg_kcal_mol:.3f}", r.classification]
        for r in results
    ]
    print(format_table(headers, rows))


def cmd_decompose(args: argparse.Namespace) -> None:
    if args.json_data:
        raw_data = json.loads(args.json_data)
        residue_bsa = {(item["chain"], int(item["id"]), item["aa"]): float(item["bsa"]) for item in raw_data}
    else:
        residue_bsa = {
            ("H", 31, "S"): 78.5,
            ("H", 33, "Y"): 95.2,
            ("H", 58, "Y"): 42.0,
            ("H", 105, "D"): 115.0,
            ("H", 106, "Y"): 88.9,
            ("H", 108, "F"): 66.3,
            ("L", 32, "Y"): 55.0,
            ("L", 58, "S"): 20.0,
            ("L", 108, "W"): 70.4,
        }

    decomps = acap.decompose_cdr_loops(residue_bsa, numbering_scheme=args.scheme)

    if args.json:
        out = [
            {
                "loop": d.loop,
                "total_bsa_a2": d.total_bsa_a2,
                "percentage_of_paratope": d.percentage_of_paratope,
                "hotspots_count": d.hotspots_count,
                "critical_hotspots_count": d.critical_hotspots_count,
            }
            for d in decomps
        ]
        print(json.dumps(out, indent=2))
        return

    print("=" * 75)
    print(f"  CDR LOOP BURIED SURFACE AREA DECOMPOSITION ({args.scheme.upper()})")
    print("=" * 75)
    headers = ["CDR Loop", "Total BSA (A^2)", "Paratope %", "Hotspots", "Critical Hotspots"]
    rows = [
        [d.loop, f"{d.total_bsa_a2:.1f}", f"{d.percentage_of_paratope:.1f}%", d.hotspots_count, d.critical_hotspots_count]
        for d in decomps
    ]
    print(format_table(headers, rows))


def cmd_cdrh3(args: argparse.Namespace) -> None:
    profile = acap.profile_cdrh3_loop(
        cdrh3_sequence=args.sequence,
        vh_residue_94=args.r94,
        vh_residue_101=args.r101,
        vh_residue_103=args.r103,
    )

    if args.json:
        out = {
            "sequence": profile.sequence,
            "length": profile.length,
            "length_category": profile.length_category,
            "base_conformation": profile.base_conformation,
            "torso_salt_bridge": profile.torso_salt_bridge,
            "aromatic_fraction": profile.aromatic_fraction,
            "gly_ser_fraction": profile.gly_ser_fraction,
            "net_charge_ph74": profile.net_charge_ph74,
        }
        print(json.dumps(out, indent=2))
        return

    print("=" * 70)
    print("  CDR-H3 LOOP STRUCTURAL & CONFORMATIONAL PROFILE")
    print("=" * 70)
    print(f"  Sequence:            {profile.sequence}")
    print(f"  Loop Length:         {profile.length} residues ({profile.length_category})")
    print(f"  Base Conformation:   {profile.base_conformation} (Shirai-Morea Rules)")
    print(f"  Torso Salt Bridge:   {'Present (H94/H101)' if profile.torso_salt_bridge else 'Absent'}")
    print(f"  Aromatic Fraction:   {profile.aromatic_fraction * 100:.1f}%")
    print(f"  Gly/Ser Flexibility: {profile.gly_ser_fraction * 100:.1f}%")
    print(f"  Net Charge (pH 7.4): {profile.net_charge_ph74:+.2f}")
    print("=" * 70)


def cmd_developability(args: argparse.Namespace) -> None:
    cdr_dict = None
    if args.cdrs:
        try:
            cdr_dict = json.loads(args.cdrs)
        except json.JSONDecodeError:
            pass

    report = acap.generate_developability_report(
        antibody_id=args.antibody_id,
        vh_sequence=args.vh,
        vl_sequence=args.vl,
        cdr_sequences=cdr_dict,
    )

    if args.json:
        out = {
            "antibody_id": report.antibody_id,
            "isoelectric_point": report.isoelectric_point,
            "net_charge_ph74": report.net_charge_ph74,
            "net_charge_ph60": report.net_charge_ph60,
            "aggregation_risk": report.aggregation_risk,
            "hydrophobic_patches": report.hydrophobic_patches,
            "liabilities_count": report.liabilities_count,
            "chemical_liabilities": [
                {
                    "type": l.liability_type,
                    "location": l.location,
                    "motif": l.motif,
                    "pos": l.start_pos,
                    "risk": l.risk_level,
                }
                for l in report.chemical_liabilities
            ],
            "framework_humanness_vh": report.framework_humanness_vh,
            "framework_humanness_vl": report.framework_humanness_vl,
            "overall_humanness": report.overall_humanness,
            "immunogenicity_tier": report.immunogenicity_tier,
            "estimated_tm_celsius": report.estimated_tm_celsius,
            "composite_developability_score": report.composite_developability_score,
            "developability_tier": report.developability_tier,
        }
        print(json.dumps(out, indent=2))
        return

    print("=" * 75)
    print(f"  DEVELOPABILITY & LIABILITY AUDIT: {report.antibody_id}")
    print("=" * 75)
    print(f"  Composite Score:      {report.composite_developability_score}/100 [{report.developability_tier}]")
    print(f"  Isoelectric Point pI: {report.isoelectric_point}")
    print(f"  Net Charge (pH 7.4):  {report.net_charge_ph74:+.2f} (pH 6.0: {report.net_charge_ph60:+.2f})")
    print(f"  Estimated Tm:         {report.estimated_tm_celsius} °C")
    print(f"  Aggregation Risk:     {report.aggregation_risk}")
    print(f"  Overall Humanness:    {report.overall_humanness}% [{report.immunogenicity_tier}]")
    print(f"    - VH Humanness:     {report.framework_humanness_vh}%")
    print(f"    - VL Humanness:     {report.framework_humanness_vl}%")
    print(f"  Chemical Liabilities: {report.liabilities_count} detected")

    if report.chemical_liabilities:
        headers = ["Type", "Location", "Motif", "Pos", "Risk"]
        rows = [
            [l.liability_type, l.location, l.motif, l.start_pos, l.risk_level]
            for l in report.chemical_liabilities
        ]
        print("\n" + format_table(headers, rows))
    print("=" * 75)


def cmd_interactive() -> None:
    print("\n--- Paratope & Developability Interactive Console ---")
    print("Select analysis mode:")
    print("1. CDR-H3 Loop Profiling & Conformation")
    print("2. Paratope Alanine Scan")
    print("3. Full Developability & Liability Screening")

    choice = input("\nSelect option [1-3]: ").strip()
    if choice == "1":
        seq = input("Enter CDR-H3 amino acid sequence: ").strip().upper()
        if not seq:
            seq = "ARDGYYYGMDV"
        r94 = input("Residue at VH-94 [default: R]: ").strip().upper() or "R"
        r101 = input("Residue at VH-101 [default: D]: ").strip().upper() or "D"
        r103 = input("Residue at VH-103 [default: W]: ").strip().upper() or "W"
        prof = acap.profile_cdrh3_loop(seq, r94, r101, r103)
        print(f"\nResult: Conformation={prof.base_conformation}, Length={prof.length} ({prof.length_category}), Aromatics={prof.aromatic_fraction*100:.1f}%, Charge={prof.net_charge_ph74:+.2f}")
    elif choice == "2":
        print("\nUsing standard benchmark paratope BSA set:")
        bsa_sample = {("H", 33, "Y"): 95.0, ("H", 100, "D"): 115.0, ("H", 101, "Y"): 88.0, ("L", 96, "W"): 70.0}
        scans = acap.perform_alanine_scan(bsa_sample)
        for s in scans:
            print(f"Chain {s.chain} Res {s.residue_id}{s.residue_name}: ddG={s.ddg_kcal_mol:.2f} kcal/mol -> {s.classification}")
    elif choice == "3":
        ab_id = input("Antibody ID: ").strip() or "mAb-Test"
        vh = input("VH Sequence: ").strip() or "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDYSGMDYWGQGTLVTVSS"
        vl = input("VL Sequence: ").strip() or "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPYTFGQGTKVEIK"
        rep = acap.generate_developability_report(ab_id, vh, vl)
        print(f"\nResult: Score={rep.composite_developability_score}/100 [{rep.developability_tier}], pI={rep.isoelectric_point}, Humanness={rep.overall_humanness}%")
    else:
        print("Invalid selection.")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="antibody-cdr-antigen-paratope",
        description="Antibody CDR-H3 Modeling, Paratope-Epitope Contact Analysis & Developability Engine",
    )
    subparsers = parser.add_subparsers(dest="command")

    # scan
    p_scan = subparsers.add_parser("scan", help="Paratope alanine scanning & hotspot calling")
    p_scan.add_argument("--json-data", help="JSON string containing residue BSA records")
    p_scan.add_argument("--file", help="Path to JSON file with BSA records")
    p_scan.add_argument("--json", action="store_true", help="Output machine-readable JSON")

    # decompose
    p_dec = subparsers.add_parser("decompose", help="Decompose BSA into CDR loop contributions")
    p_dec.add_argument("--json-data", help="JSON string of BSA records")
    p_dec.add_argument("--scheme", choices=["IMGT", "Kabat"], default="IMGT", help="Numbering scheme")
    p_dec.add_argument("--json", action="store_true", help="Output machine-readable JSON")

    # cdrh3
    p_h3 = subparsers.add_parser("cdrh3", help="Profile CDR-H3 loop conformation and properties")
    p_h3.add_argument("--sequence", required=True, help="CDR-H3 amino acid sequence")
    p_h3.add_argument("--r94", default="R", help="Residue at position 94 (default: R)")
    p_h3.add_argument("--r101", default="D", help="Residue at position 101 (default: D)")
    p_h3.add_argument("--r103", default="W", help="Residue at position 103 (default: W)")
    p_h3.add_argument("--json", action="store_true", help="Output machine-readable JSON")

    # developability
    p_dev = subparsers.add_parser("developability", help="Full developability & liability scoring")
    p_dev.add_argument("--antibody-id", default="mAb-001", help="Identifier for antibody")
    p_dev.add_argument("--vh", required=True, help="VH variable domain amino acid sequence")
    p_dev.add_argument("--vl", required=True, help="VL variable domain amino acid sequence")
    p_dev.add_argument("--cdrs", help="Optional JSON dict of CDR sequences")
    p_dev.add_argument("--json", action="store_true", help="Output machine-readable JSON")

    # interactive
    subparsers.add_parser("interactive", help="Interactive command line walkthrough")

    args = parser.parse_args(argv)

    if args.command == "scan":
        cmd_scan(args)
    elif args.command == "decompose":
        cmd_decompose(args)
    elif args.command == "cdrh3":
        cmd_cdrh3(args)
    elif args.command == "developability":
        cmd_developability(args)
    elif args.command == "interactive":
        cmd_interactive()
    else:
        parser.print_help()
        return 0

    return 0


if __name__ == "__main__":
    sys.exit(main())
