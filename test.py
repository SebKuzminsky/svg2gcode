#!/usr/bin/env python2

import gcoder


gcoder.init()

gcoder.feed(1.23)
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

gcoder.comment('z_path')
path = []
path.append({'x':0, 'y':0})
path.append({'x':2, 'y':0})
path.append({'x':2, 'y':2})
path.append({'x':0, 'y':2})
gcoder.z_path(path, 0.1, 1.0, 0.0, -0.250)


gcoder.comment('z_path2')
path = []
path.append(gcoder.line(x=2, y=0))
path.append(gcoder.arc_ccw(x=3, y=1, i=2, j=1))
path.append(gcoder.line(x=3, y=2))
path.append(gcoder.arc_ccw(x=2, y=3, i=2, j=2))
path.append(gcoder.line(x=1, y=3))
path.append(gcoder.arc_ccw(x=0, y=2, i=1, j=2))
path.append(gcoder.line(x=0, y=1))
path.append(gcoder.arc_ccw(x=1, y=0, i=1, j=1))

gcoder.absolute_arc_centers()
gcoder.feed(10)
gcoder.z_path2(path, 0.1, 1.0, 0.0, -0.250)


gcoder.done()
