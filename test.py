#!/usr/bin/env python2

import gcoder


gcoder.init()

gcoder.feed(1.23)
gcoder.g1(x=1, a=123)

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

gcoder.done()
