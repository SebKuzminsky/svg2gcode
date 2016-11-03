#!/usr/bin/env python2

import gcoder

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
