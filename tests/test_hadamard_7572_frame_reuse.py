import json
from pathlib import Path
import unittest

import numpy as np

from experiments.scalar_lift import reuse_hadamard_7572_frame as reuse


class FrameReuseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=reuse.certificate()

    def test_exact_completed_recipes_match_saved_certificate(self):
        path=Path(__file__).resolve().parents[1]/'certificates/scalar-lift/hadamard-7572-frame-reuse.json'
        self.assertEqual(self.report,json.loads(path.read_text()))
        self.assertEqual(self.report['output_orders'],[17668,22716,32812])
        self.assertTrue(self.report['new_real_orders_claimed'])
        self.assertEqual(self.report['new_real_orders_claimed_to_best_of_authors_knowledge'],
                         [17668,22716,32812])
        self.assertTrue(self.report['novelty_convention']['new_order_claims_are_qualified_by_authors_knowledge'])
        for record in self.report['outputs']:
            self.assertEqual(record['new_order_status'],'new_to_best_of_authors_knowledge')
            self.assertTrue(record['new_Hadamard_order_claimed'])
            self.assertFalse(record['prior_construction_identified'])
            self.assertEqual(record['source_used_alphabet'],[0,1,3,4])
            self.assertEqual(record['output_gram_scale'],record['output_order'])
            self.assertTrue(record['complete_algebraic_endpoint_Gram_verified'])

    def test_literal_entry_blocks_and_strict_indices(self):
        endpoint=reuse.endpoint('catalogue-p7-q6')
        for exponent in (0,1,3,4):
            i,j=next(zip(*np.where(endpoint.exponents==exponent)))
            for a,b in ((0,0),(631,1889),(2523,2523)):
                expected=(1 if exponent<3 else -1)*int((endpoint.P if exponent%3==0 else endpoint.Q)[a,b])
                self.assertEqual(endpoint.entry(int(i)*2524+a,int(j)*2524+b),expected)
        for invalid in (True,False,1.0,-1,endpoint.order):
            with self.assertRaises(ValueError):endpoint.entry(invalid,0)
            with self.assertRaises(ValueError):endpoint.entry(0,invalid)
        self.assertEqual(endpoint.entry(np.int64(0),np.int64(0)),endpoint.entry(0,0))

    def test_exact_source_gram_rejects_changed_entry_and_nonflat_palette(self):
        endpoint=reuse.endpoint('catalogue-p7-q6')
        changed=endpoint.exponents.copy()
        changed[0,0]=(changed[0,0]+1)%6
        with self.assertRaises(ValueError):reuse.verify_source_gram(changed)
        changed[0,0]=2
        with self.assertRaises(ValueError):reuse.ReusedEndpoint('invalid',changed,endpoint.P,endpoint.Q)


if __name__=='__main__':unittest.main()
