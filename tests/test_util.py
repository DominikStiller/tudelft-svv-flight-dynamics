from unittest import TestCase

import pandas as pd
from pandas._testing import assert_frame_equal

from fd.util import mean_not_none, get_closest, mean_not_nan_df


class TestUtil(TestCase):
    def test_get_closest(self):
        df = pd.DataFrame([[1], [2], [3], [4]], index=[0, 1.3, 4.5, 5.6])

        # Single rows
        self.assertEqual(get_closest(df, -3)[0], 1)
        self.assertEqual(get_closest(df, 0)[0], 1)
        self.assertEqual(get_closest(df, 0.5)[0], 2)
        self.assertEqual(get_closest(df, 1.3)[0], 2)
        self.assertEqual(get_closest(df, 4.51)[0], 4)
        self.assertEqual(get_closest(df, 100)[0], 4)

        # Multiple rows
        self.assertListEqual(list(get_closest(df, [0.5, 0.5, 0.5])[0]), [2, 2, 2])
        self.assertListEqual(list(get_closest(df, [0.5, 1.3, 4.51])[0]), [2, 2, 4])
        self.assertListEqual(list(get_closest(df, [-3, 100])[0]), [1, 4])

    def test_mean_not_none(self):
        self.assertAlmostEqual(mean_not_none([0, None, 1]), 0.5)
        self.assertAlmostEqual(mean_not_none([None, 3.4, 3.8, None]), 3.6)

    def test_mean_not_nan_df(self):
        df1 = pd.DataFrame({"a": [1, 2, pd.NA], "b": [4, 5, 6]})
        df2 = pd.DataFrame({"a": [7, pd.NA, pd.NA], "b": [pd.NA, 11, 12]})
        df3 = pd.DataFrame({"a": [13, 14, pd.NA], "b": [pd.NA, 17, pd.NA]})

        df_mean_expected = pd.DataFrame({"a": [7, 8, pd.NA], "b": [4, 11, 9]})

        df_mean = mean_not_nan_df([df1, df2, df3])

        assert_frame_equal(df_mean, df_mean_expected)
