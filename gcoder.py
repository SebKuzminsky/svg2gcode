import math


class line(object):

    """The Line class represents a linear feed move (g1) to the specified
    endpoint.

    It can be used as an element in the list of moves passed to
    z_path2()."""

    def __init__(self, x=None, y=None, z=None):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        have_arg = False
        r = "Line("

        if self.x is not None:
            r += "x=%.4f" % self.x
            have_arg = True

        if self.y is not None:
            if have_arg:
                r += ", "
            r += "y=%.4f" % self.y
            have_arg = True

        if self.z is not None:
            if have_arg:
                r += ", "
            r += "z=%.4f" % self.z
            have_arg = True

        r += ")"
        return r


class arc(object):

    """arc() is a base class representing a circular feed move (g2 or g3)
    to the specified endpoint.  It is not suitable for use by itself, you
    should use one of its subclasses, arc_cw() or arc_ccw(), instead."""

    def __init__(self, x=None, y=None, z=None, i=None, j=None, p=None):
        self.x = x
        self.y = y
        self.z = z
        self.i = i
        self.j = j
        self.p = p

    def __str__(self):
        have_arg = False
        r = self.__class__.__name__ + "("

        if self.x is not None:
            r += "x=%.4f" % self.x
            have_arg = True

        if self.y is not None:
            if have_arg:
                r += ", "
            r += "y=%.4f" % self.y
            have_arg = True

        if self.z is not None:
            if have_arg:
                r += ", "
            r += "z=%.4f" % self.z
            have_arg = True

        if self.i is not None:
            if have_arg:
                r += ", "
            r += "i=%.4f" % self.i
            have_arg = True

        if self.j is not None:
            if have_arg:
                r += ", "
            r += "j=%.4f" % self.j
            have_arg = True

        if self.p is not None:
            if have_arg:
                r += ", "
            r += "p=%.4f" % self.p
            have_arg = True

        r += ")"
        return r


class arc_cw(arc):
    """The arc_cw() class represents a circular clockwise feed move (g2)
    to the specified endpoint.

    It can be used as an element in the list of moves passed to
    z_path2()."""

class arc_ccw(arc):
    """The arc_ccw() class represents a circular counter-clockwise feed
    move (g3) to the specified endpoint.

    It can be used as an element in the list of moves passed to
    z_path2()."""


def init():
    print
    print "; init"
    print "G20          (inch)"
    print "G17          (xy plane)"
    print "G90          (absolute)"
    print "G91.1        (arc centers are relative to arc starting point)"
    cutter_comp_off()
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


def spindle_on():
    print "M3"


def spindle_off():
    print "M5"


def quill_up():
    absolute()
    cutter_comp_off()
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
    """Clockwise arc feed."""
    if i is None and j is None:
        raise TypeError, "gcoder.g2() without i or j"
    print "G2",
    if x is not None: print "X%.4f" % x,
    if y is not None: print "Y%.4f" % y,
    if z is not None: print "Z%.4f" % z,
    if i is not None: print "I%.4f" % i,
    if j is not None: print "J%.4f" % j,
    if p is not None: print "P%.4f" % p,
    print


def g3(x=None, y=None, z=None, i=None, j=None, p=None):
    """Counter-clockwise arc feed."""
    if i is None and j is None:
        raise TypeError, "gcoder.g3() without i or j"
    print "G3",
    if x is not None: print "X%.4f" % x,
    if y is not None: print "Y%.4f" % y,
    if z is not None: print "Z%.4f" % z,
    if i is not None: print "I%.4f" % i,
    if j is not None: print "J%.4f" % j,
    if p is not None: print "P%.4f" % p,
    print


#
# Cutter compensation handling.
#

def cutter_comp_off():
    print "G40          (cutter comp off)"

def cancel_cutter_comp():
    print "; gcoder: calling program used obsolete cancel_cutter_comp() function, use cutter_comp_off() instead"
    cutter_comp_off()

def g40():
    print "; gcoder: calling program used obsolete g40() function, use cutter_comp_off() instead"
    cutter_comp_off()


def cutter_comp_left(**kwargs):

    """Enable cutter diameter compensation on the left side of the
    programmed path.

    When called with no argument, uses the diameter of the currently
    loaded tool (from the tool table).

    When called with the `diameter` argument, uses the specified diameter.

    When called with the `tool` argument (and without the `diameter`
    argument), uses the diameter of the specified tool number (from the
    tool table)."""

    if 'diameter' in kwargs:
        print "G41.1 D%.4f   (cutter comp left, diameter mode)" % kwargs['diameter']
    elif 'tool' in kwargs:
        print "G41 D%d   (cutter comp left, tool-number mode)" % kwargs['tool']
    else:
        print "G41   (cutter comp left, current tool)"


def cutter_comp_right(**kwargs):

    """Enable cutter diameter compensation on the right side of the
    programmed path.

    When called with no argument, uses the diameter of the currently
    loaded tool (from the tool table).

    When called with the `diameter` argument, uses the specified diameter.

    When called with the `tool` argument (and without the `diameter`
    argument), uses the diameter of the specified tool number (from the
    tool table)."""

    if 'diameter' in kwargs:
        print "G42.1 D%.4f   (cutter comp right, diameter mode)" % kwargs['diameter']
    elif 'tool' in kwargs:
        print "G42 D%d   (cutter comp right, tool-number mode)" % kwargs['tool']
    else:
        print "G42   (cutter comp right, current tool)"

def g42_1(comp_diameter):
    print "; gcoder: calling program used obsolete g42_1() function, use cutter_comp_right() instead"
    cutter_comp_right(diameter=comp_diameter)


def g81(retract, x=None, y=None, z=None):
    print "G81",
    if x is not None: print "X%.4f" % x,
    if y is not None: print "Y%.4f" % y,
    if z is not None: print "Z%.4f" % z,
    print "R%.4f" % retract,
    print


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
    The motion is this:

        1. Set Z to z_start.

        2. If z_top_of_work is below z_start, set Z level to z_top_of_work
           and feed down to Z (otherwise don't move the controlled point).

        3. Reduce Z by depth_of_cut, but not below z_target.

        4. Feed to each waypoint in path, in order starting with the
           first and ending with the last.

        5. After arriving at the last waypoint, if Z is not yet down to
           z_target: feed to the first waypoint while ramping down by
           depth_of_cut (but not below z_target), then go back to step
           4 for another trip around the path at this Z level.

        6. After reaching step 5 with Z at z_target, feed back to the
           first waypoint while keeping Z at the z_target level, thereby
           cutting away the ramp left by the previous iteration.


    The first move of this function is to the first waypoint in the
    path, at a Z level that's depth-of-cut below the lower of z_start
    and z_top_of_material (but not below z_target).  If you position
    the cutter above the *last* waypoint in the path, you'll get a nice
    consistent ramp down to the first waypoint each time around the path.

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


def z_path2(path, depth_of_cut, z_start, z_top_of_work, z_target):

    """This function traverses a path (a list of Line, ArcCW, and ArcCCW
    objects), cutting a little deeper on each pass.  The motion is this:

        1. Set Z to z_start.

        2. If z_top_of_work is below z_start, set Z level to z_top_of_work
           and feed down to Z (otherwise don't move the controlled point).

        3. Reduce Z by depth_of_cut, but not below z_target.

        4. Feed to each waypoint in path, in order starting with the
           first and ending with the last.

        5. After arriving at the last waypoint, if Z is not yet down to
           z_target: feed to the first waypoint while ramping down by
           depth_of_cut (but not below z_target), then go back to step
           4 for another trip around the path at this Z level.

        6. After reaching step 5 with Z at z_target, feed back to the
           first waypoint while keeping Z at the z_target level, thereby
           cutting away the ramp left by the previous iteration.


    The first move of this function is to the first waypoint in the
    path, at a Z level that's depth-of-cut below the lower of z_start
    and z_top_of_material (but not below z_target).  If you position
    the cutter above the *last* waypoint in the path, you'll get a nice
    consistent ramp down to the first waypoint each time around the path.

    When the function returns the controlled point is once again at the
    first waypoint in the path, all the way down at Z=z_target."""

    def handle_item(item):
        if type(item) is line:
            g1(x=item.x, y=item.y, z=z)
        elif type(item) is arc_cw:
            g2(x=item.x, y=item.y, z=z, i=item.i, j=item.j, p=item.p)
        elif type(item) is arc_ccw:
            g3(x=item.x, y=item.y, z=z, i=item.i, j=item.j, p=item.p)
        else:
            raise TypeError('z_path2() only accepts line(), arc_cw(), and arc_ccw() objects')

    z = z_start

    if z > z_top_of_work:
        z = z_top_of_work
        g1(z = z)

    while z > z_target:
        z = z - depth_of_cut
        if z < z_target:
            z = z_target

        for item in path:
            handle_item(item)

    # Cut away the last ramp we left behind.
    handle_item(path[0])


def helix_hole(x, y, z_retract, z_start, z_bottom, diameter, doc):

    """This function helix-mills a hole.  The motion is this:

        1. Rapid to Z=z_retract.

        2. Rapid to the vincinity of (X, Y).

        3. Rapid to Z=z_start.

        4. Helix down to Z=z_bottom, descending not more than Z=doc
           per revolution.

        5. One more full circle at Z=z_bottom, to flatten the floor.

        6. Feed to the center of the hole and up off the floor a
           little bit.

        7. Rapid up to Z=z_retract."""

    r = diameter / 2.0

    z_range = z_start - z_bottom
    full_circles = math.ceil(z_range / doc)

    absolute_arc_centers()

    # get in position for the cut
    g0(z=z_retract)
    g0(x=x+r, y=y)
    g0(z=z_start)

    # helix down, then flatten the bottom
    g2(x=x+r, y=y, z=z_bottom, p=full_circles, i=x, j=y)
    g2(x=x+r, y=y, z=z_bottom, i=x, j=y)

    # extract the tool from the work
    g1(x=x, y=y, z=z_bottom + 0.025)
    g0(z=z_retract)


def saw_square(x_start, y_start, z_start, x_end, y_end, z_end, max_doc, rapid_plunge=True, final_retract=True):

    """Cuts back and forth between (X=x_start, Y=y_start) and (X=x_end,
    Y=y_end), plunging Z down (rapid or feed) at the end of each pass.

    The actual depth of cut may be reduced a little from max_doc to
    achieve equal depth of cut on each pass, while minimizing the number
    of passes.

    Upon return the tool will be positioned at either (X=x_start,
    Y=y_start) or at (X=x_end, Y=y_end), and at either Z=z_start (if
    final_retract is True) or Z=z_end (if final_retract is False).

    Motion:

        Initial Motion:

            Rapid to X=x_start, Y=y_start.

            Spindle on.

            Rapid to Z=z_start.

        Cycle:

            Plunge Z down by actual_doc, but not below z_end (rapid or
            feed, determined by rapid_plunge).

            Feed to X=x_end, Y=y_end.

            If Z is at z_end, goto Done.

            Plunge Z down by actual_doc, but not below z_end (rapid or
            feed, determined by rapid_plunge).

            Feed to X=x_start, Y=y_start.

            If Z is at z_end, goto Done.

            Goto Cycle.

        Done:

            If final_retract is True, rapid to Z=z_start."""

    z_range = z_start - z_end
    num_passes = math.ceil(z_range / max_doc)
    doc = z_range / num_passes

    g0(x=x_start, y=y_start)

    spindle_on();

    z = z_start
    g0(z=z)

    while z > z_end:
        z = z - doc
        if z < z_end:
            z = z_end
        if rapid_plunge:
            g0(z=z)
        else:
            g1(z=z)
        g1(x=x_end, y=y_end)

        if z > z_end:
            z = z - doc
            if z < z_end:
                z = z_end
            if rapid_plunge:
                g0(z=z)
            else:
                g1(z=z)
            g1(x=x_start, y=y_start)

    if final_retract:
        g0(z=z_start)

