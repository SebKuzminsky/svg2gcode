#!/usr/bin/env python2

import gcoder

gcoder.init()
gcoder.relative_arc_centers()

gcoder.speed(1000)
gcoder.set_feed_rate(10)
gcoder.spindle_on()

path = []
path.append(gcoder.line(x=4))
path.append(gcoder.arc_ccw(x=5, y=1, i=0, j=1))
path.append(gcoder.line(y=4))
path.append(gcoder.arc_ccw(x=4, y=5, i=-1, j=0))
path.append(gcoder.line(x=1))
path.append(gcoder.arc_ccw(x=0, y=4, i=0, j=-1))
path.append(gcoder.line(y=1))
path.append(gcoder.arc_ccw(x=1, y=0, i=1, j=0))

gcoder.comment("first")

gcoder.comment("first: getting into position")
gcoder.g0(z=1)
gcoder.g0(x=1, y=0)
gcoder.g0(z=0.1)
gcoder.g1(z=0)

gcoder.comment("first: z_path2")
gcoder.z_path2(path=path, depth_of_cut=0.750, z_target=-1)

gcoder.comment("first: getting out of the cut")
gcoder.g1(x=4, y=0.025, z=-0.975)
gcoder.g0(z=1)


gcoder.comment("second")

gcoder.comment("second: getting into position")
gcoder.g0(z=1)
gcoder.g0(x=1, y=0)
gcoder.g0(z=0.1)
gcoder.g1(z=0)

gcoder.comment("second: z_path2")
gcoder.z_path2(path=path, depth_of_cut=10, z_target=-2)

gcoder.comment("second: getting out of the cut")
gcoder.g1(x=4, y=0.025, z=-1.975)
gcoder.g0(z=1)

gcoder.m2()

