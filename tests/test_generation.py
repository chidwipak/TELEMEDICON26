"""Unit tests for constrained analogue generation (RDKit, no data/weights needed)."""
import unittest

from rdkit import Chem, DataStructs

from _common import REQUIRED_COLUMNS, TOY_FS_PRECEDENT, TOY_H_PRECEDENT, toy_db
from cliffforge import generate as G


class TestSmilesHandling(unittest.TestCase):
    def test_flat_smiles_is_canonical_and_stereo_free(self):
        self.assertEqual(G.flat_smiles("C[C@H](N)C(=O)O"), G.flat_smiles("C[C@@H](N)C(=O)O"))
        self.assertEqual(G.flat_smiles("OC(=O)C(N)C"), G.flat_smiles("CC(N)C(=O)O"))

    def test_invalid_smiles_returns_none(self):
        self.assertIsNone(G.flat_smiles("not_a_smiles"))
        self.assertIsNone(G.flat_smiles("C1CC"))  # unclosed ring

    def test_attach_returns_sanitized_molecule(self):
        m = G.attach("c1ccccc1[*:1]", "Cl[*:1]")
        self.assertIsNotNone(m)
        self.assertEqual(Chem.MolToSmiles(m), "Clc1ccccc1")
        self.assertEqual(Chem.MolToSmiles(G.attach("c1ccccc1[*:1]", "H")), "c1ccccc1")

    def test_attach_rejects_garbage(self):
        self.assertIsNone(G.attach("not_smiles", "Cl[*:1]"))


class TestEditSites(unittest.TestCase):
    def test_cuts_report_sidechain_size_and_environment(self):
        mol = Chem.MolFromSmiles("COC(=O)c1ccc(F)cc1")
        cuts = G.cuts(mol, max_side=12, min_scaf_frac=0.6)
        sides = {side for _, side, _, _ in cuts}
        self.assertIn("F[*:1]", sides)           # aromatic fluorine as removable sidechain
        self.assertIn("H", sides)                # aromatic C-H sites are editable
        for _, side, n_side, _ in cuts:
            self.assertEqual(n_side == 0, side == "H")
            self.assertLessEqual(n_side, 12)


class TestToyGeneration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.db = toy_db(TOY_FS_PRECEDENT)
        cls.parent = "COC(=O)c1ccc(F)cc1"
        cls.df = G.generate(cls.parent, cls.db, known={}, seed_id="toy")
        old = G.MATCH_H_SITES
        G.MATCH_H_SITES = True          # opt-in switch: see the KNOWN LIMITATION note in generate.py
        try:
            cls.db_h = toy_db(TOY_H_PRECEDENT)
            cls.df_h = G.generate("COC(=O)c1ccccc1", cls.db_h, known={}, seed_id="toyH")
        finally:
            G.MATCH_H_SITES = old

    def test_observed_transformation_is_applied(self):
        self.assertEqual(list(self.df.candidate_smiles), ["COC(=O)c1ccc(Cl)cc1"])
        row = self.df.iloc[0]
        self.assertEqual(row.transformation, "F[*:1]>>Cl[*:1]")
        self.assertEqual(int(row.n_precedent_pairs), 3)
        self.assertEqual(row.filter_status, "pass")

    def test_output_schema(self):
        for df in (self.df, self.df_h):
            missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
            self.assertEqual(missing, [])
            self.assertTrue(df.valid_rdkit.dtype == bool)
            self.assertTrue(df.parent_tanimoto.between(0, 1).all())

    def test_parent_candidate_relationship_and_edit_metadata(self):
        for df, parent in ((self.df, self.parent), (self.df_h, "COC(=O)c1ccccc1")):
            pm = Chem.MolFromSmiles(parent)
            pfp = G.fp(pm)
            for r in df[df.valid_rdkit].itertuples():
                cm = Chem.MolFromSmiles(r.candidate_smiles)
                self.assertIsNotNone(cm)                                         # sanitizes
                self.assertNotEqual(Chem.MolToSmiles(cm, isomericSmiles=False),  # never the parent itself
                                    Chem.MolToSmiles(pm, isomericSmiles=False))
                self.assertEqual(r.edit_size_heavy_atoms, r.removed_heavy_atoms + r.added_heavy_atoms)
                self.assertEqual(r.heavy_atom_delta, r.candidate_heavy_atoms - r.parent_heavy_atoms)
                self.assertEqual(r.heavy_atom_delta, r.added_heavy_atoms - r.removed_heavy_atoms)
                self.assertAlmostEqual(r.parent_tanimoto, DataStructs.TanimotoSimilarity(pfp, G.fp(cm)), places=9)
                scaffold = Chem.AdjustQueryProperties(Chem.MolFromSmiles(r.scaffold_smiles.replace("[*:1]", "*")))
                self.assertTrue(cm.HasSubstructMatch(scaffold))                  # shared scaffold is preserved

    def test_deduplication(self):
        df = self.df_h
        self.assertGreater(len(df), len(df[df.unique_candidate]))                # equivalent ring positions -> duplicates
        uniq = df[df.unique_candidate].candidate_smiles
        self.assertEqual(len(uniq), uniq.nunique())
        dups = df[~df.unique_candidate]
        self.assertTrue(dups.filter_reason.str.contains("duplicate|same_as_parent").all())
        self.assertTrue((dups.filter_status == "rejected").all())

    def test_funnel_is_monotone(self):
        f = G.funnel(self.df_h)
        self.assertGreaterEqual(f["generated_edits"], f["valid_rdkit"])
        self.assertGreaterEqual(f["valid_rdkit"], f["unique_non_parent"])
        self.assertGreaterEqual(f["unique_non_parent"], f["pass_all_filters"])

    def test_precedent_filter_rejects_unsupported_edits(self):
        df = G.generate(self.parent, self.db, known={}, seed_id="toy", min_precedent=99)
        self.assertTrue((df.filter_status == "rejected").all())
        self.assertTrue(df.filter_reason.str.contains("precedent<99").all())

    def test_transform_db_contains_only_supplied_chemistry(self):
        sides = set(self.db.n_side)
        for forbidden in ("Br[*:1]", "O[*:1]"):
            self.assertNotIn(forbidden, sides)
        self.assertIn("Cl[*:1]", sides)

    def test_default_does_not_match_hydrogen_sites_known_limitation(self):
        """Documents the behaviour of the reported runs: H<->X transformations are not matched by default."""
        self.assertFalse(G.MATCH_H_SITES)
        db = toy_db(TOY_H_PRECEDENT)
        df = G.generate("COC(=O)c1ccccc1", db, known={}, seed_id="toyH0")
        self.assertEqual(len(df), 0)

    def test_known_compound_annotation_carries_no_labels(self):
        known = {"COC(=O)c1ccc(Cl)cc1": "test"}
        df = G.generate(self.parent, self.db, known=known, seed_id="toy")
        self.assertEqual(df.iloc[0].known_in_moonshot, "test")
        self.assertFalse(any("pIC50" in c and "measured" in c for c in df.columns))


if __name__ == "__main__":
    unittest.main()
