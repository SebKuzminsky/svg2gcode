#!/usr/bin/env python3

import gcoder


gcoder.init()

gcoder.set_feed_rate(1.23)
gcoder.g1(x=1)

path = []
path.append({'x':0, 'y':0})
path.append({'z':-1})
path.append({'x':2, 'y':0})
path.append({'x':2, 'y':2})
path.append({'x':0, 'y':2})
path.append({'x':0, 'y':0})
path.append({'z':1})
gcoder.g1(path=path)

gcoder.g83(x=1, y=2, z=-3, retract=0.1, delta=0.2)

gcoder.drill_hog(diameter=0.250, retract=0.1, delta=0.2, z_drill=-1, x0=0, y0=0, x1=-1, y1=2)
gcoder.drill_hog(diameter=0.250, retract=0.1, delta=0.2, z_drill=-1, x0=0, y0=0, x1=-1, y1=0.506, xy_finishing_allowance=0.001)

gcoder.comment('saw_square, with rapid_plunge not specified')
gcoder.saw_square(
    x_start = 0,
    y_start = 0,
    z_start = 0,
    x_end = 1,
    y_end = 0,
    z_end = -1,
    max_doc = 0.260
)

gcoder.comment('saw_square, with rapid_plunge False')
gcoder.saw_square(
    x_start = 0,
    y_start = 0,
    z_start = 0,
    x_end = 1,
    y_end = 0,
    z_end = -1,
    max_doc = 0.260,
    rapid_plunge = False
)

gcoder.comment('saw_square, with rapid_plunge False and final_retract False')
gcoder.saw_square(
    x_start = 0,
    y_start = 0,
    z_start = 0,
    x_end = 1,
    y_end = 0,
    z_end = -1,
    max_doc = 0.260,
    rapid_plunge = False,
    final_retract = False
)

gcoder.comment('z_path')
path = []
path.append({'x':0, 'y':0})
path.append({'x':2, 'y':0})
path.append({'x':2, 'y':2})
path.append({'x':0, 'y':2})
gcoder.z_path(path, 0.1, 1.0, 0.0, -0.250)

gcoder.done()
