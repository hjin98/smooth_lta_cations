import copy
import random
import unittest
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice
from pymatgen.core.operations import SymmOp

from lta_cation_population_grouped import (
    CandidateVariant,
    PhysicalSite,
    PlacedIon,
    build_site_symmetry_group,
    canonical_site_occupation_signature,
    site_occupation_signature,
    transform_site_signature,
    generate_unique_lta_cation_batch_from_specs,
)


def make_site(site_id, frac, variant_offset=0.01):
    frac = np.asarray(frac, dtype=float)
    variants = [
        CandidateVariant(frac=frac.copy(), site_label="NA3", ring_class="4", site_id=site_id, variant_id=0),
        CandidateVariant(frac=(frac + np.array([variant_offset, 0.0, 0.0])) % 1.0, site_label="NA3", ring_class="4", site_id=site_id, variant_id=1),
    ]
    return PhysicalSite(site_id=site_id, ring_class="4", site_label="NA3", variants=variants)


def ion(element, site, variant_id=0):
    v = site.variants[variant_id]
    return PlacedIon(element, v.frac.copy(), site.site_label, site.ring_class, site.site_id, variant_id)


class UniqueBatchTests(unittest.TestCase):
    def setUp(self):
        self.lattice = Lattice.cubic(10.0)
        # Four sites related by quarter-turns around the fractional origin.
        sites = [
            make_site("s0", [0.25, 0.25, 0.0]),
            make_site("s1", [0.75, 0.25, 0.0]),
            make_site("s2", [0.75, 0.75, 0.0]),
            make_site("s3", [0.25, 0.75, 0.0]),
        ]
        self.by_ring = {"4": sites, "6": [], "8": []}
        identity = SymmOp.from_xyz_str("x,y,z")
        rot90 = SymmOp.from_rotation_and_translation(
            [[0, -1, 0], [1, 0, 0], [0, 0, 1]], [0, 0, 0]
        )
        self.group = build_site_symmetry_group(
            self.lattice,
            self.by_ring,
            [identity, rot90],
            mapping_tolerance_ang=0.2,
        )

    def test_generator_closure_finds_four_rotations(self):
        self.assertEqual(self.group.order, 4)

    def test_rotated_populations_have_same_canonical_signature(self):
        sites = self.by_ring["4"]
        adjacent_a = [ion("Li", sites[0]), ion("Li", sites[1])]
        adjacent_b = [ion("Li", sites[1]), ion("Li", sites[2])]
        opposite = [ion("Li", sites[0]), ion("Li", sites[2])]
        ca = canonical_site_occupation_signature(self.by_ring, adjacent_a, self.group)
        cb = canonical_site_occupation_signature(self.by_ring, adjacent_b, self.group)
        co = canonical_site_occupation_signature(self.by_ring, opposite, self.group)
        self.assertEqual(ca, cb)
        self.assertNotEqual(ca, co)

    def test_variant_choice_does_not_change_signature(self):
        sites = self.by_ring["4"]
        a = [ion("Na", sites[0], 0)]
        b = [ion("Na", sites[0], 1)]
        self.assertEqual(site_occupation_signature(self.by_ring, a), site_occupation_signature(self.by_ring, b))
        self.assertEqual(
            canonical_site_occupation_signature(self.by_ring, a, self.group),
            canonical_site_occupation_signature(self.by_ring, b, self.group),
        )

    def test_species_are_part_of_equivalence_state(self):
        sites = self.by_ring["4"]
        li_na = [ion("Li", sites[0]), ion("Na", sites[1])]
        na_li = [ion("Na", sites[0]), ion("Li", sites[1])]
        # A quarter turn can move both labels, but cannot exchange chemical species.
        c1 = canonical_site_occupation_signature(self.by_ring, li_na, self.group)
        c2 = canonical_site_occupation_signature(self.by_ring, na_li, self.group)
        self.assertNotEqual(c1, c2)

    def test_signature_transform_uses_old_to_new_mapping(self):
        sig = ("Li", "empty", "Na", "empty")
        quarter_turn = next(p for p in self.group.permutations if p[0] == 1)
        transformed = transform_site_signature(sig, quarter_turn)
        self.assertEqual(transformed[1], "Li")
        self.assertEqual(transformed[3], "Na")


class UniqueBatchIntegrationTests(unittest.TestCase):
    def test_batch_api_returns_pairwise_distinct_canonical_signatures(self):
        cif = Path(__file__).with_name("synthetic_c4_lta.cif")
        batch = generate_unique_lta_cation_batch_from_specs(
            cif,
            batch_size=6,
            n4="Li:2,Na:2,K:0",
            n6="Li:0,Na:0,K:0",
            n8="Li:0,Na:0,K:0",
            seed=789,
            variant_group_cutoffs={"4": 0.05, "6": 0.05, "8": 0.05},
            max_attempts=100,
            _enforce_supported_scope=False,
        )
        self.assertEqual(len(batch.members), 6)
        self.assertEqual(len(set(batch.canonical_signatures)), 6)
        self.assertEqual(batch.symmetry_group.order, 4)

    def test_impossible_batch_fails_instead_of_emitting_duplicates(self):
        cif = Path(__file__).with_name("synthetic_c4_full_lta.cif")
        with self.assertRaisesRegex(RuntimeError, "Could generate only 1"):
            generate_unique_lta_cation_batch_from_specs(
                cif,
                batch_size=2,
                n4="Li:8,Na:0,K:0",
                n6="Li:0,Na:0,K:0",
                n8="Li:0,Na:0,K:0",
                seed=321,
                variant_group_cutoffs={"4": 0.05, "6": 0.05, "8": 0.05},
                max_attempts=8,
                _enforce_supported_scope=False,
            )


if __name__ == "__main__":
    unittest.main()
