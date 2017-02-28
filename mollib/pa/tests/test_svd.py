import unittest
import re

from mollib import Molecule
from mollib.pa.process_molecule import Process
from mollib.pa.data_readers import read_pa_string
from mollib.pa.svd import calc_pa_SVD
from mollib.pa.analysis import calc_statistics
from mollib.pa import settings


pna_rdc = """
3N-H -12.76
4N-H 10.38
5N-H 24.1
6N-H 32.01
7N-H 18.29
10N-H 10.28
11N-H 0.18
12N-H 21.05
13N-H 26.43
14N-H 22.59
15N-H -14.69
16N-H -11.95
17N-H -33.63
20N-H 21.15
21N-H -28.36
23N-H 31.37
25N-H 25.54
26N-H 27.83
27N-H 27.43
28N-H 31.04
29N-H 27.24
30N-H 32.26
31N-H 24.44
32N-H 29.29
33N-H 30.54
34N-H 27.52
35N-H 14.01
39N-H -16.96
40N-H 12.77
41N-H -36.71
42N-H 2.57
43N-H -6.2
44N-H 28.46
45N-H 23.34
46N-H 24.14
47N-H 30.2
48N-H 17.46
49N-H 24.27
50N-H 25.49
51N-H 12.97
52N-H 16.81
54N-H 15.95
55N-H 28.73
56N-H -11.9
58N-H -18.41
59N-H 6.02
60N-H -33.17
61N-H -42.13
62N-H -19.61
63N-H 20.24
64N-H -35.4
65N-H -24.24
66N-H -9.85
67N-H 19.67
68N-H 30.03
69N-H 5.25
70N-H -3.1
71N-H -24.93
72N-H -9.21
"""
#
# pna_rdc = """
# 2N-H 13.2
# 3N-H -8.78
# 4N-H -6.81
# 5N-H -3.19
# 6N-H -0.92
# 7N-H 4.62
# 8N-H -13.62
# 10N-H -5.57
# 11N-H 14.87
# 12N-H -0.59
# 13N-H 1.42
# 14N-H -3.66
# 15N-H -11.84
# 16N-H 0.30
# 17N-H -2.18
# 18N-H 17.18
# 20N-H 11.28
# 21N-H -3.94
# 22N-H -20.53
# 23N-H -1.11
# 25N-H 3.23
# 26N-H -0.39
# 27N-H -2.61
# 28N-H 3.01
# 29N-H -0.82
# 30N-H -1.01
# 31N-H -1.83
# 32N-H 3.81
# 33N-H -0.72
# 34N-H -2.31
# 35N-H 10.024
# 36N-H 20.59
# 39N-H -14.97
# 40N-H 0.42
# 41N-H -24.12
# 42N-H -9.53
# 43N-H -10.97
# 44N-H -1.69
# 45N-H 4.15
# 46N-H 0.97
# 47N-H -1.79
# 48N-H 20.029
# 49N-H 3.89
# 50N-H -1.84
# 51N-H 13.56
# 52N-H 22.73
# 54N-H 4.17
# 55N-H -0.022
# 56N-H -16.10
# 57N-H -22.51
# 58N-H -10.29
# 59N-H -9.23
# 60N-H -22.90
# 61N-H -16.95
# 62N-H -15.85
# 64N-H -19.87
# 65N-H -18.07
# 66N-H -13.20
# 67N-H -5.26
# 68N-H -1.75
# 69N-H -4.42
# 70N-H -12.02
# 71N-H -19.48
#
# 2H -4.7
# 3H -0.1
# 4H  2.6
# 5H -4.7
# 6H  1.6
# 7H -8.2
# 8H  3.0
# 10H -4.1
# 11H -2.0
# 12H -2.2
# 13H  6.1
# 14H -0.8
# 15H  8.0
# 16H  1.6
# 17H  0.4
# 18H -8.2
# 20H  0.6
# 21H -2.6
# 23H -0.1
# 25H -3.1
# 26H  3.9
# 27H -4.6
# 28H  2.7
# 29H -3.3
# 30H  1.0
# 31H  1.4
# 32H -2.5
# 33H  1.4
# 34H -0.7
# 35H -0.3
# 36H -5.2
# 39H  0.9
# 40H  1.6
# 41H  3.0
# 42H  6.4
# 43H  1.0
# 44H  4.8
# 45H -4.2
# 46H  4.1
# 47H -2.0
# 48H -3.9
# 49H  0.2
# 50H  3.5
# 51H -1.1
# 52H -6.0
# 54H -4.2
# 55H -6.2
# 56H  6.4
# 57H  2.8
# 58H  2.4
# 59H  4.0
# 60H  5.1
# 61H  6.9
# 62H  3.1
# 63H -5.1
# 64H  8.1
# 65H  2.0
# 66H  2.5
# 67H -1.2
# 68H -4.2
# 69H -1.5
# 70H  1.9
# 71H  6.5
# """


re_label_sort = re.compile(r'[A-Z]?\.?(\d+)(.+)')
def sort_key(label):
    m = re_label_sort.match(label)
    if m:
        groups = m.groups()
        return (groups[1], int(groups[0]))
    else:
        return None

class TestSVD(unittest.TestCase):

    def test_svd_ubiquitin_NH(self):
        """Test the SVD of ubiquitin NH RDCs"""

        mol = Molecule('2MJB')
        process = Process(mol)
        magnetic_interactions = process.process()

        data = read_pa_string(pna_rdc)

        data_pred, Saupe_components, stats = calc_pa_SVD(magnetic_interactions,
                                                         data)

        print("{:.1f}%".format(stats['Q']*100.))

        rss = 0.
        count = 0
        for label in sorted(data_pred, key=sort_key):
            if label in data:
                fmt = "{:<10} {:5.1f} {:5.1f}"
                print(fmt.format(label, data[label].value,
                                        data_pred[label].value))
            else:
                fmt = "{:<10} {:5.1f}"
                print(fmt.format(label, data_pred[label].value))

    def test_stats(self):
        """Test the calculation of the Q-factor."""
        mol = Molecule('2MJB')

        # Calculate the Q-factor first from the bond lengths
        process = Process(mol)
        settings.calculate_from_bonds = True
        magnetic_interactions = process.process()

        data = read_pa_string(pna_rdc)
        data_pred, Saupe_components, stats1 = calc_pa_SVD(magnetic_interactions,
                                                          data)

        # The fit Q-factor should be better than 10%
        self.assertLessEqual(stats1['Q'], 0.10)

        # Now calculate it from static DCCs. This Q-factor should be different
        # (yet still low)
        process = Process(mol)
        settings.calculate_from_bonds = False
        magnetic_interactions = process.process()

        data = read_pa_string(pna_rdc)
        data_pred, Saupe_components, stats2 = calc_pa_SVD(magnetic_interactions,
                                                          data)

        # The fit Q-factor should be better than 10%
        self.assertLessEqual(stats2['Q'], 0.10)

        # The residual sum squared should be different between calculating the
        # RDCs using bond lengths vs static values
        self.assertNotEqual(stats1['RSS'], stats2['RSS'])
