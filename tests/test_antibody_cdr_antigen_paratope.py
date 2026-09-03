#!/usr/bin/env python3
"""
Comprehensive Unit Test Suite for Antibody CDR-H3 Modeling & Paratope Analytics.
"""

import unittest
import math
from antibody_cdr_antigen_paratope import (
    calculate_ddg_from_bsa,
    classify_hotspot,
    perform_alanine_scan,
    decompose_cdr_loops,
    analyze_interface_contacts,
    profile_cdrh3_loop,
    scan_hydrophobic_patches,
    find_chemical_liabilities,
    compute_net_charge,
    compute_isoelectric_point,
    compute_sequence_identity,
    generate_developability_report,
    HotSpotResult,
    LoopDecomposition,
    InterfaceContactSummary,
    ChemicalLiability,
    CDRH3StructureProfile,
    DevelopabilityReport,
)


class TestAntibodyCdrParatope(unittest.TestCase):

    # 1. Delta-Delta-G and Hotspot Thresholds
    def test_calculate_ddg_standard(self):
        # 100 A^2 BSA * 22 cal/(mol*A^2) = 2200 cal/mol = 2.2 kcal/mol
        ddg = calculate_ddg_from_bsa(100.0)
        self.assertAlmostEqual(ddg, 2.2, places=3)

    def test_calculate_ddg_zero(self):
        ddg = calculate_ddg_from_bsa(0.0)
        self.assertEqual(ddg, 0.0)

    def test_calculate_ddg_negative_raises(self):
        with self.assertRaises(ValueError):
            calculate_ddg_from_bsa(-10.0)

    def test_classify_hotspot_tiers(self):
        self.assertEqual(classify_hotspot(3.0), "Critical Hotspot")
        self.assertEqual(classify_hotspot(2.5), "Critical Hotspot")
        self.assertEqual(classify_hotspot(1.8), "Hotspot")
        self.assertEqual(classify_hotspot(1.5), "Hotspot")
        self.assertEqual(classify_hotspot(0.8), "Moderate")
        self.assertEqual(classify_hotspot(0.3), "Neutral")

    # 2. Alanine Scanning
    def test_alanine_scan_sorting_and_classification(self):
        bsa_data = {
            ("H", 33, "Y"): 120.0,  # ddG = 2.64 -> Critical Hotspot
            ("H", 52, "W"): 80.0,   # ddG = 1.76 -> Hotspot
            ("L", 32, "A"): 30.0,   # ddG = 0.66 -> Moderate
            ("L", 50, "G"): 10.0,   # ddG = 0.22 -> Neutral
        }
        results = perform_alanine_scan(bsa_data)
        self.assertEqual(len(results), 4)
        # Should be sorted descending by ddG
        self.assertEqual(results[0].residue_id, 33)
        self.assertEqual(results[0].classification, "Critical Hotspot")
        self.assertEqual(results[1].residue_id, 52)
        self.assertEqual(results[1].classification, "Hotspot")
        self.assertEqual(results[2].residue_id, 32)
        self.assertEqual(results[2].classification, "Moderate")
        self.assertEqual(results[3].residue_id, 50)
        self.assertEqual(results[3].classification, "Neutral")

    def test_alanine_scan_empty_raises(self):
        with self.assertRaises(ValueError):
            perform_alanine_scan({})

    def test_alanine_scan_negative_bsa_raises(self):
        with self.assertRaises(ValueError):
            perform_alanine_scan({("H", 31, "S"): -5.0})

    # 3. CDR Loop Decomposition
    def test_decompose_cdr_loops_imgt(self):
        bsa_data = {
            ("H", 30, "S"): 50.0,   # CDR-H1 (27-38)
            ("H", 33, "Y"): 100.0,  # CDR-H1 (27-38) -> Hotspot
            ("H", 60, "T"): 40.0,   # CDR-H2 (56-65)
            ("H", 110, "W"): 130.0, # CDR-H3 (105-117) -> Critical Hotspot
            ("L", 30, "N"): 30.0,   # CDR-L1 (27-38)
            ("L", 110, "Y"): 50.0,  # CDR-L3 (105-117)
        }
        decomps = decompose_cdr_loops(bsa_data, numbering_scheme="IMGT")
        loop_map = {d.loop: d for d in decomps}

        self.assertAlmostEqual(loop_map["CDR-H1"].total_bsa_a2, 150.0)
        self.assertEqual(loop_map["CDR-H1"].hotspots_count, 1)
        self.assertAlmostEqual(loop_map["CDR-H3"].total_bsa_a2, 130.0)
        self.assertEqual(loop_map["CDR-H3"].critical_hotspots_count, 1)
        self.assertAlmostEqual(loop_map["CDR-H2"].total_bsa_a2, 40.0)
        self.assertEqual(loop_map["CDR-L2"].total_bsa_a2, 0.0)

    def test_decompose_cdr_loops_kabat(self):
        bsa_data = {
            ("H", 33, "Y"): 80.0,   # CDR-H1 (31-35)
            ("H", 55, "G"): 45.0,   # CDR-H2 (50-65)
            ("H", 98, "D"): 90.0,   # CDR-H3 (95-102)
        }
        decomps = decompose_cdr_loops(bsa_data, numbering_scheme="Kabat")
        loop_map = {d.loop: d for d in decomps}
        self.assertAlmostEqual(loop_map["CDR-H1"].total_bsa_a2, 80.0)
        self.assertAlmostEqual(loop_map["CDR-H2"].total_bsa_a2, 45.0)
        self.assertAlmostEqual(loop_map["CDR-H3"].total_bsa_a2, 90.0)

    # 4. Interface Contact Network
    def test_interface_contact_summary(self):
        hbonds = [("H33", "A105"), ("H100", "A108")]
        salt_bridges = [("H101", "A45")]
        pi_stacking = [("H33", "A80")]
        hydrophobic = [("H52", "A22"), ("L96", "A24")]

        summary = analyze_interface_contacts(hbonds, salt_bridges, pi_stacking, hydrophobic)
        self.assertEqual(summary.hydrogen_bonds, 2)
        self.assertEqual(summary.salt_bridges, 1)
        self.assertEqual(summary.pi_stacking, 1)
        self.assertEqual(summary.hydrophobic_contacts, 2)
        self.assertEqual(summary.total_interactions, 6)
        # Energy = 2*1.5 + 1*3.0 + 1*2.0 + 2*0.8 = 3.0 + 3.0 + 2.0 + 1.6 = 9.6 kcal/mol
        self.assertAlmostEqual(summary.estimated_interaction_energy_kcal, 9.6)

    # 5. CDR-H3 Loop Profiling
    def test_cdrh3_profiling_kinked_base(self):
        # Kinked: R94 + D101 + W103
        profile = profile_cdrh3_loop(
            cdrh3_sequence="ARDGYYYGMDV",
            vh_residue_94="R",
            vh_residue_101="D",
            vh_residue_103="W",
        )
        self.assertEqual(profile.base_conformation, "Kinked")
        self.assertTrue(profile.torso_salt_bridge)
        self.assertEqual(profile.length, 11)
        self.assertEqual(profile.length_category, "Medium")
        self.assertGreater(profile.aromatic_fraction, 0.2)

    def test_cdrh3_profiling_extended_base(self):
        # Extended: A94 (no basic) + D101
        profile = profile_cdrh3_loop(
            cdrh3_sequence="AGGYGMDV",
            vh_residue_94="A",
            vh_residue_101="D",
            vh_residue_103="W",
        )
        self.assertEqual(profile.base_conformation, "Extended")
        self.assertFalse(profile.torso_salt_bridge)
        self.assertEqual(profile.length_category, "Short")

    def test_cdrh3_profiling_invalid_amino_acid(self):
        with self.assertRaises(ValueError):
            profile_cdrh3_loop("ARDGYY123MDV")

    # 6. Hydrophobic Patches
    def test_hydrophobic_patch_detection(self):
        # Hydrophobic stretch: VILFIVA
        seq = "SSGVIILFIVAASS"
        patches = scan_hydrophobic_patches(seq, window=7, threshold=2.0)
        self.assertTrue(len(patches) > 0)
        # Hydrophilic stretch
        polar_seq = "SSSDDDDEEEERRR"
        patches_polar = scan_hydrophobic_patches(polar_seq, window=7, threshold=1.0)
        self.assertEqual(len(patches_polar), 0)

    # 7. Chemical Liabilities
    def test_find_chemical_liabilities_motifs(self):
        seqs = {
            "CDR-H1": "SYAMSW",     # Oxidation (M, W)
            "CDR-H2": "AISGSGGSTYY",
            "CDR-H3": "ARNGDGYW",   # Deamidation (NG), Isomerization (DG), Oxidation (W)
            "VL": "NTSPDP",         # Glycosylation (NTS), Acid-labile (DP)
        }
        liabs = find_chemical_liabilities(seqs)
        types = [l.liability_type for l in liabs]
        self.assertIn("Deamidation", types)
        self.assertIn("Isomerization", types)
        self.assertIn("Oxidation", types)
        self.assertIn("Glycosylation", types)
        self.assertIn("Acid Labile", types)

    def test_unpaired_cysteines(self):
        seqs = {"CDR-H3": "ARDCYYCGMDV"}  # Even count (2 Cys) -> No unpaired alert
        liabs1 = find_chemical_liabilities(seqs)
        self.assertNotIn("Unpaired Cys", [l.liability_type for l in liabs1])

        seqs_odd = {"CDR-H3": "ARDCYYGMDV"}  # Odd count (1 Cys) -> Unpaired alert
        liabs2 = find_chemical_liabilities(seqs_odd)
        self.assertIn("Unpaired Cys", [l.liability_type for l in liabs2])

    # 8. Net Charge and Isoelectric Point (pI)
    def test_net_charge_basic_peptide(self):
        # Sequence with 3 Lysines and 2 Arginines
        seq = "KKKRR"
        charge = compute_net_charge(seq, ph=7.4)
        self.assertGreater(charge, 4.0)

    def test_net_charge_acidic_peptide(self):
        # Sequence with 4 Aspartate and 3 Glutamate
        seq = "DDDDEEE"
        charge = compute_net_charge(seq, ph=7.4)
        self.assertLess(charge, -6.0)

    def test_isoelectric_point_calculation(self):
        # Basic peptide should have high pI (>9)
        pI_basic = compute_isoelectric_point("KKKKRRRR")
        self.assertGreater(pI_basic, 10.0)

        # Acidic peptide should have low pI (<4)
        pI_acidic = compute_isoelectric_point("DDDDEEEE")
        self.assertLess(pI_acidic, 4.0)

    # 9. Sequence Identity and Humanness
    def test_sequence_identity_exact(self):
        seq = "EVQLVESGGGLVQPGGSLRLSCAAS"
        self.assertEqual(compute_sequence_identity(seq, seq), 100.0)

    def test_sequence_identity_mismatch(self):
        seq1 = "AAAAA"
        seq2 = "AAAAT"
        self.assertEqual(compute_sequence_identity(seq1, seq2), 80.0)

    # 10. Full Developability Report
    def test_generate_developability_report_clean(self):
        vh = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDYSGMDYWGQGTLVTVSS"
        vl = "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPYTFGQGTKVEIK"
        cdrs = {
            "CDR-H1": "SYAMS",
            "CDR-H2": "AISGSGGSTYYADSVKG",
            "CDR-H3": "DYSGMDY",
            "CDR-L1": "RASQSISSYLN",
            "CDR-L2": "AASSLQS",
            "CDR-L3": "QQSYSTPYT",
        }
        report = generate_developability_report("mAb-001", vh, vl, cdrs)
        self.assertEqual(report.antibody_id, "mAb-001")
        self.assertGreater(report.composite_developability_score, 60.0)
        self.assertIn(report.developability_tier, ["High Developability", "Moderate Developability"])
        self.assertGreater(report.overall_humanness, 70.0)
        self.assertGreater(report.estimated_tm_celsius, 60.0)

    def test_developability_empty_sequence_raises(self):
        with self.assertRaises(ValueError):
            generate_developability_report("mAb-Bad", "", "")

    # 11. Additional Edge & Corner Cases
    def test_ultra_long_cdrh3_loop(self):
        # 24 aa CDR-H3
        seq = "ARDGYYYGMDVARDGYYYGMDVAA"
        profile = profile_cdrh3_loop(seq)
        self.assertEqual(profile.length_category, "Ultra-long")
        self.assertEqual(profile.length, 24)

    def test_extreme_ph_net_charge(self):
        seq = "ACDEFGHIKLMNPQRSTVWY"
        charge_ph1 = compute_net_charge(seq, ph=1.0)
        charge_ph14 = compute_net_charge(seq, ph=14.0)
        # At pH 1, positive groups charged, negative groups protonated -> net positive
        self.assertGreater(charge_ph1, 0.0)
        # At pH 14, positive groups deprotonated, negative groups ionized -> net negative
        self.assertLess(charge_ph14, 0.0)

    def test_single_amino_acid_pi(self):
        pi_k = compute_isoelectric_point("K")
        self.assertGreater(pi_k, 9.0)
        pi_d = compute_isoelectric_point("D")
        self.assertLess(pi_d, 4.0)

    def test_empty_interface_contacts(self):
        summary = analyze_interface_contacts([], [], [], [])
        self.assertEqual(summary.total_interactions, 0)
        self.assertEqual(summary.estimated_interaction_energy_kcal, 0.0)

    def test_cli_execution_smoke(self):
        from cli import main
        import io
        from unittest.mock import patch
        with patch("sys.stdout", new=io.StringIO()) as fake_out:
            code = main(["cdrh3", "--sequence", "ARDGYYYGMDV", "--json"])
            self.assertEqual(code, 0)
            self.assertIn("Kinked", fake_out.getvalue())

    def test_cli_batch(self):
        from cli import main
        import os
        import tempfile
        sample_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "sample.csv")
        with tempfile.TemporaryDirectory() as tmpdir:
            out_file = os.path.join(tmpdir, "out_batch.csv")
            code = main(["batch", "-i", sample_path, "-o", out_file])
            self.assertEqual(code, 0)
            self.assertTrue(os.path.exists(out_file))


if __name__ == "__main__":
    unittest.main()

