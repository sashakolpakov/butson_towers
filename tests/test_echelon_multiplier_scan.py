"""Independent finite row-count controls and strategy-certificate replay."""
from itertools import product
import unittest

from experiments.scalar_lift import echelon_multiplier_scan as scan


class EchelonMultiplierTests(unittest.TestCase):
    def test_small_chart_rows_by_literal_signed_orbit_words(self):
        # Enumerate the actual eight signed coordinate values on every orbit,
        # rather than reusing either norm-equation enumeration in the module.
        for chart in scan.arithmetic_charts(4):
            p, m = chart['p'], chart['m']
            N, action = chart['nonzero_orbits'], chart['coordinate_action']
            realized = set()
            for word in product(range(8), repeat=N):
                sums = [1, 0, 0, 0]
                for state in word:
                    axis, sign = state//2, 1-2*(state % 2)
                    if action == 'C3' and axis:
                        for j in range(1, 4):
                            sums[j] += sign*(m//3)
                    else:
                        sums[axis] += sign*m
                if sum(x*x for x in sums) == p:
                    realized.add((sums[0], *sorted(map(abs, sums[1:]), reverse=True)))
            self.assertEqual(realized, {tuple(r['row']) for r in chart['necessary_DC_rows']})

    def test_five_dense_comparison_rows_are_replayed_and_costed(self):
        report = scan.build_report()
        rows = {r['p']: r for r in report['dense_Williamson_comparison']['rows']}
        self.assertEqual(set(rows), {463, 919, 1801, 1951, 2269})
        for p, m, N in [(463, 22, 21), (919, 54, 17), (1801, 72, 25),
                        (1951, 78, 25), (2269, 84, 27)]:
            self.assertEqual((rows[p]['m'], rows[p]['nonzero_orbits']), (m, N))
            self.assertEqual(rows[p]['unfiltered_axis_bits'], 2*N)
            self.assertEqual(int(rows[p]['unfiltered_axis_word_count']), 4**N)
        low = report['tiers'][0]
        self.assertEqual(low['source_listed_selected_orders'], [2676, 4460, 5892, 9260])
        self.assertFalse(report['interpretation']['matrix_search_performed'])


if __name__ == '__main__':
    unittest.main()
