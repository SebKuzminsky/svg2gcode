#!/usr/bin/env python2

import gcoder

gcoder.relative_arc_centers()

gcoder.g0(z=1)
gcoder.g0(x=1, y=0)

gcoder.speed(1000)
gcoder.feed(10)

gcoder.g1(z=0)

path = []
path.append(gcoder.line(x=4))
path.append(gcoder.arc_ccw(x=5, y=1, i=0, j=1))
path.append(gcoder.line(y=4))
path.append(gcoder.arc_ccw(x=4, y=5, i=-1, j=0))
path.append(gcoder.line(x=1))
path.append(gcoder.arc_ccw(x=0, y=4, i=0, j=-1))
path.append(gcoder.line(y=1))
path.append(gcoder.arc_ccw(x=1, y=0, i=1, j=0))

gcoder.z_path2(path=path, depth_of_cut=0.750, z_start=0, z_top_of_work=0, z_target=-1)

gcoder.g1(x=4, y=0.025, z=-0.975)
gcoder.g0(z=1)

gcoder.m2()

