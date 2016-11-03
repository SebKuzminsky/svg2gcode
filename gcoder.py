import math


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


def drill_hog(diameter, retract, delta, z_drill, x0, y0, x1, y1, finishing_allowance=None):

    """Drills as many holes as will fit in a rectangular grid, within the
    rectangle defined by (x0, y0) and (x1, y1).  If finishing_tolerance
    is specified, all the holes will stay at least that far away from
    the edges of the specified rectangle."""

    radius = diameter/2.0

    min_x = min(x0, x1)
    max_x = max(x0, x1)

    min_y = min(y0, y1)
    max_y = max(y0, y1)

    if finishing_allowance != None:
        min_x = min_x + finishing_allowance
        max_x = max_x - finishing_allowance
        min_y = min_y + finishing_allowance
        max_y = max_y - finishing_allowance

    x_range = max_x - min_x
    y_range = max_y - min_y

    num_in_x = int(math.floor(x_range / diameter))
    num_in_y = int(math.floor(y_range / diameter))

    min_x = min_x + radius
    max_x = max_x - radius

    min_y = min_y + radius
    max_y = max_y - radius

    x_range = x_range - diameter
    y_range = y_range - diameter

    for x_index in range(0, num_in_x):
        if num_in_x > 1:
            x = min_x + ((x_index / float(num_in_x - 1)) * x_range)
        else:
            x = min_x + (x_range / 2.0)

        for y_index in range(0, num_in_y):
            if num_in_y > 1:
                y = min_y + ((y_index / float(num_in_y - 1)) * y_range)
            else:
                y = min_y + (y_range / 2.0)

            print "G83 ",
            print "X%.4f" % x,
            print "Y%.4f" % y,
            print "Z%.4f" % z_drill,
            print "R%.4f" % retract,
            print "Q%.4f" % delta,
            print

