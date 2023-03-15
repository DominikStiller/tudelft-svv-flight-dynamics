from unittest import TestCase

from fd.util import mean_not_none


class TestUtil(TestCase):
    def test_mean_not_none(self):
        self.assertAlmostEqual(mean_not_none([0, None, 1]), 0.5)
        self.assertAlmostEqual(mean_not_none([None, 3.4, 3.8, None]), 3.6)
