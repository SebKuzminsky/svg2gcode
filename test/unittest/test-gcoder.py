#!/usr/bin/env python3

import unittest
import os
import sys

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'svgpathtools'))
import svgpathtools

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import gcoder


class test_split_path_at_intersections(unittest.TestCase):

    def test_no_intersection(self):
        p = svgpathtools.path.Path()
        p.append(svgpathtools.path.Line((0+0j), (10+0j)))
        p.append(svgpathtools.path.Line((10+0j), (10+10j)))
        p.append(svgpathtools.path.Line((10+10j), (0+10j)))
        p.append(svgpathtools.path.Line((0+10j), (0+0j)))

        out = gcoder.split_path_at_intersections(p, debug=True)

        self.assertEqual(len(out), 1)
        self.assertEqual(svgpathtools.path.Path(*out[0]), p)


if __name__ == '__main__':
    unittest.main()
