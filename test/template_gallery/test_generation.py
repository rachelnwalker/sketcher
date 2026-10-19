"""Geometry and sidechain invariants for the gallery generator."""

import importlib.util
from pathlib import Path
import unittest

try:
    import numpy as np
    from rdkit import Chem, Geometry
except ImportError as exc:
    raise unittest.SkipTest(
        'Requires the native RDKit gallery runtime') from exc

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    'gallery_generate', ROOT / 'tools/template_gallery/generate.py')
gallery = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gallery)


class GalleryGenerationTest(unittest.TestCase):

    def test_variants_are_distinct_and_preserve_core_bonds(self):
        base = Chem.MolFromSmiles('C1CC2CCC1C2')
        variants, limitation = gallery.make_variants(base, 7, 5)
        self.assertIsNone(limitation)
        self.assertEqual(len(variants), 5)
        self.assertEqual(len({gallery.canonical(v['mol']) for v in variants}),
                         5)
        for variant in variants:
            mol = variant['mol']
            Chem.SanitizeMol(mol)
            for bond in base.GetBonds():
                new_bond = mol.GetBondBetweenAtoms(bond.GetBeginAtomIdx(),
                                                   bond.GetEndAtomIdx())
                self.assertEqual(new_bond.GetBondType(), bond.GetBondType())
            self.assertTrue(
                all(i >= base.GetNumAtoms() for i in variant['sidechainAtoms']))
        self.assertEqual(variants[-1]['substituent'], 'double substitution')

    def test_source_coordinates_cannot_seed_generated_pair(self):
        base = Chem.MolFromSmiles('C1CC2CCC1C2')
        gallery.rdDepictor.Compute2DCoords(base)
        variant = gallery.make_variants(base, 2, 2)[0][0]
        variant['mol'] = Chem.Mol(base)
        first = gallery.depict_case(variant, base.GetNumAtoms())
        conf = variant['mol'].GetConformer()
        for idx in range(conf.GetNumAtoms()):
            p = conf.GetAtomPosition(idx)
            conf.SetAtomPosition(idx, Geometry.Point3D(p.x + 1000, p.y - 700,
                                                       0))
        second = gallery.depict_case(variant, base.GetNumAtoms())
        self.assertEqual(first['enabled']['svg'], second['enabled']['svg'])
        self.assertEqual(first['disabled']['svg'], second['disabled']['svg'])

    def test_alignment_uses_rotation_and_translation_without_reflection(self):
        a = Chem.MolFromSmiles('CCCO')
        gallery.rdDepictor.Compute2DCoords(a)
        b = Chem.Mol(a)
        original = np.array(a.GetConformer().GetPositions())[:, :2]
        rotated = original @ np.array([[0, -1], [1, 0]]) + [100, 30]
        for idx, (x, y) in enumerate(rotated):
            b.GetConformer().SetAtomPosition(
                idx, Geometry.Point3D(float(x), float(y), 0))
        rmsd, _ = gallery.align_pair(a, b, a.GetNumAtoms())
        self.assertLess(rmsd, 1e-8)
        for idx, (x, y) in enumerate(original):
            b.GetConformer().SetAtomPosition(
                idx, Geometry.Point3D(float(-x), float(y), 0))
        rmsd, _ = gallery.align_pair(a, b, a.GetNumAtoms())
        self.assertGreater(rmsd, 0.1)

    def test_crossings_require_interior_intersection(self):
        self.assertTrue(gallery.crossing((0, 0), (1, 1), (0, 1), (1, 0)))
        self.assertFalse(gallery.crossing((0, 0), (1, 0), (1, 0), (2, 1)))
        self.assertFalse(gallery.crossing((0, 0), (1, 0), (0.5, 0), (2, 0)))


if __name__ == '__main__':
    unittest.main()
