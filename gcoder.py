
def feed(feed_rate_units_per_minute):
    print "F %.4f" % feed_rate_units_per_minute


def speed(spindle_rpm):
    print "S %.4f" % spindle_rpm


def g0(path=None, x=None, y=None, z=None, a=None, b=None, c=None, u=None, v=None, w=None):
    if path is not None:
        for waypoint in path:
            g0(**waypoint)
    else:
        print "G0",
        if x is not None: print "X%.4f" % x,
        if y is not None: print "Y%.4f" % y,
        if z is not None: print "Z%.4f" % z,
        if a is not None: print "A%.4f" % a,
        if b is not None: print "B%.4f" % b,
        if c is not None: print "C%.4f" % c,
        if u is not None: print "U%.4f" % u,
        if v is not None: print "V%.4f" % v,
        if w is not None: print "W%.4f" % w,
        print


def g1(path=None, x=None, y=None, z=None, a=None, b=None, c=None, u=None, v=None, w=None):
    if path is not None:
        for waypoint in path:
            g1(**waypoint)
    else:
        print "G1",
        if x is not None: print "X%.4f" % x,
        if y is not None: print "Y%.4f" % y,
        if z is not None: print "Z%.4f" % z,
        if a is not None: print "A%.4f" % a,
        if b is not None: print "B%.4f" % b,
        if c is not None: print "C%.4f" % c,
        if u is not None: print "U%.4f" % u,
        if v is not None: print "V%.4f" % v,
        if w is not None: print "W%.4f" % w,
        print


def g83(retract, delta, x=None, y=None, z=None):
    print "G83",
    if x is not None: print "X%.4f" % x,
    if y is not None: print "Y%.4f" % y,
    if z is not None: print "Z%.4f" % z,
    print "R%.4f" % retract,
    print "Q%.4f" % delta,
    print
