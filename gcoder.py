import math


def init():
    print
    print "; init"
    print "G20          (inch)"
    print "G17          (xy plane)"
    print "G90          (absolute)"
    print "G91.1        (arc centers are relative to arc starting point)"
    print "G40          (disable cutter comp)"
    print "G54          (switch to coordinate system 1)"
    print "G94          (units/minute feed mode)"
    print "G99          (in canned cycles, retract to the Z coordinate specified by the R word)"
    print "G64 P0.0005  (enable path blending, but stay within 0.0005 of the programmed path)"
    print "G49          (turn off tool length compensation)"
    print "G80          (turn off canned cycles)"
    print


def comment(msg):
    if msg:
        print ";", msg
    else:
        print


def absolute():
    print "G90"


def absolute_arc_centers():
    print "G90.1"


def relative_arc_centers():
    print "G91.1"


def cancel_cutter_comp():
    g40()


def spindle_on():
    print "M3"


def spindle_off():
    print "M5"


def quill_up():
    absolute()
    cancel_cutter_comp()
    print "G53 G0 Z0"
    spindle_off()


def presentation_position():
    imperial()
    quill_up()

    # rapid to presentation position
    # table centered in X, all the way forward towards the user
    print "G53 G0 X9 Y12"


def m2():
    print
    print "M2"


def done():
    print
    print "; done"
    presentation_position()
    print "M2"


def imperial():
    print "G20"


def feed(feed_rate_units_per_minute):
    print "F %.4f" % feed_rate_units_per_minute


def speed(spindle_rpm):
    print "S %d" % spindle_rpm


def g0(path=None, x=None, y=None, z=None, a=None, b=None, c=None, u=None, v=None, w=None):
    if path is not None:
        print
        print "; g0 path"
        for waypoint in path:
            g0(**waypoint)
        print
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
        print
        print "; g1 path"
        for waypoint in path:
            g1(**waypoint)
        print
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


def g2(x=None, y=None, z=None, i=None, j=None, p=None):
    print "G2",
    if x is not None: print "X%.4f" % x,
    if y is not None: print "Y%.4f" % y,
    if z is not None: print "Z%.4f" % z,
    if i is not None: print "I%.4f" % i,
    if j is not None: print "J%.4f" % j,
    if p is not None: print "P%.4f" % p,
    print


def g3(x=None, y=None, z=None, i=None, j=None, p=None):
    print "G3",
    if x is not None: print "X%.4f" % x,
    if y is not None: print "Y%.4f" % y,
    if z is not None: print "Z%.4f" % z,
    if i is not None: print "I%.4f" % i,
    if j is not None: print "J%.4f" % j,
    if p is not None: print "P%.4f" % p,
    print


def g40():
    print "G40"


def g42_1(comp_diameter):
    print "G42.1 D%.4f" % comp_diameter


def g83(retract, delta, x=None, y=None, z=None):
    print "G83",
    if x is not None: print "X%.4f" % x,
    if y is not None: print "Y%.4f" % y,
    if z is not None: print "Z%.4f" % z,
    print "R%.4f" % retract,
    print "Q%.4f" % delta,
    print


def drill_hog(diameter, retract, delta, z_drill, x0, y0, x1, y1, xy_finishing_allowance=None, z_finishing_allowance=None):

    """Drills as many evenly spaced holes as will fit in a rectangular
    grid, within the rectangle defined by (x0, y0) and (x1, y1).
    The specified rectangle describes the material contour, the holes
    will be inset from the edges by the drill's radius.

    If finishing_tolerance is specified, all the holes will stay at least
    that far away from the material contour the specified rectangle.
    If z_finishing_allowance is specified the holes will end that far
    above the specified drill depth."""

    print
    print "; drill hog"

    radius = diameter/2.0

    min_x = min(x0, x1)
    max_x = max(x0, x1)

    min_y = min(y0, y1)
    max_y = max(y0, y1)

    if xy_finishing_allowance != None:
        min_x = min_x + xy_finishing_allowance
        max_x = max_x - xy_finishing_allowance
        min_y = min_y + xy_finishing_allowance
        max_y = max_y - xy_finishing_allowance

    if z_finishing_allowance != None:
        z_drill = z_drill + z_finishing_allowance

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

    print


def z_path(path, depth_of_cut, z_start, z_top_of_work, z_target):

    """This function traverses a path (a list of waypoints), cutting a
    little deeper on each pass.  The waypoints are (X, Y) coordinates.

    This function ramps down on the first of these cuts (to our depth
    of cut), then goes around at that depth.  Repeat until we're reach
    the bottom, then cut away the final ramp.

    The first move of this function is to the first waypoint in the path.
    When the function returns the controlled point is once again at the
    first waypoint in the path, all the way down at Z=z_target."""

    z = z_start

    if z > z_top_of_work:
        z = z_top_of_work
        g1(z = z)

    while z > z_target:
        z = z - depth_of_cut
        if z < z_target:
            z = z_target

        for waypoint in path:
            g1(x=waypoint['x'], y=waypoint['y'], z=z)

    # Cut away the last ramp we left behind.
    g1(**path[0])
