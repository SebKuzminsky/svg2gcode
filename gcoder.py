from __future__ import print_function

import cairosvg.parser
import math
import os
import re
import sys

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'svgpathtools'))
import svgpathtools


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


class svg():
    def __init__(self, svg_file):
        self.svg_file = svg_file
        self.height = 0.0

        self.cairo = cairosvg.parser.Tree(url=self.svg_file)

        m = re.match('([0-9.]+)([a-zA-Z]*)', self.cairo['height'])
        if m == None:
            raise SystemExit, "failed to parse SVG height: %s" % c['height']

        if len(m.groups()) == 1:
            self.height = float(m.group(1))
        elif len(m.groups()) == 2:
            print("FIXME: there's a string here ('%s') specifying the units, i hope it's 1:1 with the SVG user units, and i hope it's mm" % m.group(2), file=sys.stderr)
            self.height = float(m.group(1))
        else:
            raise SystemExit, "weird result from re"

        self.paths, self.attributes = svgpathtools.svg2paths(self.svg_file)


    def to_mm_x(self, x_mm):
        if type(x_mm) != float:
            raise SystemExit, 'non-float input'
        return x_mm


    def to_mm_y(self, y_mm):
        if type(y_mm) != float:
            raise SystemExit, 'non-float input'
        return self.height - y_mm


    def to_mm(self, xy):
        if type(xy) != complex:
            raise SystemExit, 'non-complex input'
        x = self.to_mm_x(xy.real)
        y = self.to_mm_y(xy.imag)
        return (x, y)


def split_path_at_intersections(path_list):

    """`path_list` is a list of connected path segments.  This function
    identifies each place where the path intersects intself, and splits
    each non-self-intersecting subset of the path into a separate
    path list.  This may involve splitting segments.

    Returns a list of path lists."""

    first_path = []
    second_path = []
    for i in range(len(path_list)):
        this_seg = path_list[i]
        for j in range(i+2, len(path_list)):
            if (i == 0) and (j == len(path_list) - 1):
                continue
            other_seg = path_list[j]

            print("intersecting:", file=sys.stderr)
            print("    this:", this_seg, file=sys.stderr)
            print("    other:", other_seg, file=sys.stderr)
            intersections = this_seg.intersect(other_seg)
            print("    intersections:", len(intersections), file=sys.stderr)

            # FIXME: deal with multiple intersections here
            for intersection in intersections:
                this_first_seg, this_second_seg = this_seg.split(intersection[0])
                other_first_seg, other_second_seg = other_seg.split(intersection[1])
                if type(this_seg) == svgpathtools.path.Arc:
                    print("split an arc:", this_seg, file=sys.stderr)
                    print("    t:", intersection[0], file=sys.stderr)
                    print("    ", this_first_seg, file=sys.stderr)
                    print("    ", this_second_seg, file=sys.stderr)
                if type(other_seg) == svgpathtools.path.Arc:
                    print("split an arc:", other_seg, file=sys.stderr)
                    print("    t:", intersection[1], file=sys.stderr)
                    print("    ", other_first_seg, file=sys.stderr)
                    print("    ", other_second_seg, file=sys.stderr)

                # FIXME: This fixup is bogus, but the two segments'
                # `t` parameters don't put the intersection at the
                # same point...
                other_first_seg.end = this_first_seg.end
                other_second_seg.start = other_first_seg.end

                assert(complex_close_enough(this_first_seg.end, this_second_seg.start))
                assert(complex_close_enough(this_first_seg.end, other_first_seg.end))
                assert(complex_close_enough(this_first_seg.end, other_second_seg.start))

                first_path.append(this_first_seg)
                first_path.append(other_second_seg)
                for k in range(j+1, len(path_list)):
                    first_path.append(path_list[k])

                second_path.append(this_second_seg)
                for k in range(i+1, j):
                    second_path.append(path_list[k])
                second_path.append(other_first_seg)

                first_paths = split_path_at_intersections(first_path)
                second_paths = split_path_at_intersections(second_path)

                return first_paths + second_paths

        # This_seg did not intersect any of the other segments in the
        # path list, so it goes in the first path.
        first_path.append(this_seg)

    # This path list did not intersect itself, so we return a list
    # containing just the input path.
    return [first_path]


def approximate_path_area(path):

    """Approximates the path area by converting each Arc to 1,000
    Lines."""

    assert(path.isclosed())
    tmp = svgpathtools.path.Path()
    for seg in path:
        if type(seg) == svgpathtools.path.Arc:
            for i in range(0, 1000):
                t0 = i/1000.0
                t1 = (i+1)/1000.0
                l = svgpathtools.path.Line(start=seg.point(t0), end=seg.point(t1))
                tmp.append(l)
        else:
            tmp.append(seg)
    return tmp.area()


def offset_path(path, offset_distance, steps=100):
    """Takes an svgpathtools.path.Path object, `path`, and a float
    distance, `offset_distance`, and returns the parallel offset curve
    (in the form of another svgpathtools.path.Path)."""

    # FIXME: Join adjacent non-intersecting line segments with arcs.

    # This only works on closed paths.
    print("input path:", file=sys.stderr)
    print(path, file=sys.stderr)
    assert(path.isclosed())


    #
    # First generate a list of Path elements (Lines and Arcs),
    # corresponding to the offset versions of the Path elements in the
    # input path.
    #

    offset_path_list = []
    for seg in path:
        if type(seg) == svgpathtools.path.Line:
            start = seg.point(0) + (offset_distance * seg.normal(0))
            end = seg.point(1) + (offset_distance * seg.normal(1))
            offset_path_list.append(svgpathtools.Line(start, end))

        elif type(seg) == svgpathtools.path.Arc and (seg.radius.real == seg.radius.imag):
            # Circular arcs remain arcs, elliptical arcs become linear
            # approximations below.
            # FIXME: this only works for concave arcs
            new_radius = seg.radius.real - offset_distance
            start = seg.point(0) + (offset_distance * seg.normal(0))
            end = seg.point(1) + (offset_distance * seg.normal(1))
            if new_radius > 0.002:
                radius = complex(seg.radius.real - offset_distance, seg.radius.imag - offset_distance)
                offset_arc = svgpathtools.path.Arc(
                    start = start,
                    end = end,
                    radius = radius,
                    rotation = seg.rotation,
                    large_arc = seg.large_arc,
                    sweep = seg.sweep
                )
                offset_path_list.append(offset_arc)
            elif new_radius > 0.0:
                # Offset Arc still exists but radius is smaller than
                # the minimum that LinuxCNC accepts, replace with a Line.
                offset_arc = svgpathtools.path.Line(start = start, end = end)
                offset_path_list.append(offset_arc)

        else:
            # Deal with any segment that's not a line or a circular arc.
            # This includes elliptic arcs and bezier curves.  Use linear
            # approximation.
            #
            # FIXME: Steps should probably be computed dynamically to make
            #     the length of the *offset* line segments manageable.
            points = []
            for k in range(steps):
                t = k / float(steps)
                normal = seg.normal(t)
                offset_vector = offset_distance * normal
                points.append(seg.point(t) + offset_vector)
            for k in range(len(points)-1):
                start = points[k]
                end = points[k+1]
                offset_path_list.append(svgpathtools.Line(start, end))


    #
    # Find all the places where one segment intersects the next, and
    # trim to the intersection.
    #

    for i in range(len(offset_path_list)):
        this_seg = offset_path_list[i]
        if (i+1) < len(offset_path_list):
            next_seg = offset_path_list[i+1]
        else:
            next_seg = offset_path_list[0]

        # FIXME: I'm not sure about this part.
        print("intersecting", file=sys.stderr)
        print("    this", this_seg, file=sys.stderr)
        print("    next", next_seg, file=sys.stderr)
        intersections = this_seg.intersect(next_seg)
        print("    intersections:", intersections, file=sys.stderr)
        for intersection in intersections:
            point = this_seg.point(intersection[0])
            if not complex_close_enough(point, this_seg.end):
                this_seg.end = this_seg.point(intersection[0])
                next_seg.start = this_seg.end


    #
    # Find all the places where adjacent segments do not end/start close
    # to each other, and join them with Arcs.
    #

    joined_offset_path_list = []
    for i in range(len(offset_path_list)):
        this_seg = offset_path_list[i]
        if (i+1) < len(offset_path_list):
            next_seg = offset_path_list[i+1]
        else:
            next_seg = offset_path_list[0]

        if complex_close_enough(this_seg.end, next_seg.start):
            joined_offset_path_list.append(this_seg)
            continue

        # FIXME: Choose values for `large_arc` and `sweep` correctly here.
        joining_arc = svgpathtools.path.Arc(
            start = this_seg.end,
            end = next_seg.start,
            radius = complex(offset_distance, offset_distance),
            rotation = 0,
            large_arc = False,
            sweep = True  # sweep means "clockwise"
        )
        joined_offset_path_list.append(this_seg)
        joined_offset_path_list.append(joining_arc)
        print("these segments don't join:", file=sys.stderr)
        print(this_seg, file=sys.stderr)
        print(next_seg, file=sys.stderr)
        print("adding joining arc:", file=sys.stderr)
        print(joining_arc, file=sys.stderr)

    offset_path_list = joined_offset_path_list


    #
    # Find the places where the path intersects itself, split into
    # multiple separate paths in those places.
    #

    offset_paths_list = split_path_at_intersections(offset_path_list)


    #
    # Smooth the path: adjacent segments whose start/end points are
    # "close enough" to each other are adjusted so they actually touch.
    #

    for path_list in offset_paths_list:
        for i in range(len(path_list)):
            this_seg = path_list[i]
            if (i+1) < len(path_list):
                next_seg = path_list[i+1]
            else:
                next_seg = path_list[0]
            if complex_close_enough(this_seg.end, next_seg.start):
                next_seg.start = this_seg.end
            else:
                print("discontinuity (seg %d and following):" % i, file=sys.stderr)
                print("    this_seg.end:", this_seg.end, file=sys.stderr)
                print("    next_seg.start:", next_seg.start, file=sys.stderr)


    #
    # Done!  Convert each path list to a Path object, sanity check,
    # and return them.
    #

    offset_paths = []

    path_area = approximate_path_area(path)

    for path_list in offset_paths_list:
        offset_path = svgpathtools.Path(*path_list)
        print("offset path:", file=sys.stderr)
        print(offset_path, file=sys.stderr)
        print("continuous?", offset_path.iscontinuous(), file=sys.stderr)
        print("", file=sys.stderr)
        assert(offset_path.isclosed())
        offset_path_area = approximate_path_area(offset_path)
        if path_area * offset_path_area > 0.0:
            # input path and offset path go in the same direction
            offset_paths.append(offset_path)

    return offset_paths


def path_to_gcode(svg, path):
    absolute_arc_centers()
    (x, y) = svg.to_mm(path[0].start)
    g0(z=10.000)
    g0(x=x, y=y)

    spindle_on()
    g1(z=0.000)

    for element in path:
        if type(element) == svgpathtools.path.Line:
            (start_x, start_y) = svg.to_mm(element.start)
            (end_x, end_y) = svg.to_mm(element.end)
            g1(x=end_x, y=end_y)
        elif type(element) == svgpathtools.path.Arc:
            # FIXME: g90.1 or g91.1?
            if element.radius.real != element.radius.imag:
                raise ValueError, "arc radii differ: %s", element
            (end_x, end_y) = svg.to_mm(element.end)
            (center_x, center_y) = svg.to_mm(element.center)
            if element.sweep:
                g2(x=end_x, y=end_y, i=center_x, j=center_y)
            else:
                g3(x=end_x, y=end_y, i=center_x, j=center_y)
        else:
            # Deal with any segment that's not a line or a circular arc,
            # this includes elliptic arcs and bezier curves.  Use linear
            # approximation.
            #
            # FIXME: The number of steps should probably be dynamically
            #     adjusted to make the length of the *offset* line
            #     segments manageable.
            steps = 1000
            for k in range(steps):
                t = k / float(steps)
                end = element.point(t)
                (end_x, end_y) = svg.to_mm(end)
                g1(x=end_x, y=end_y)

    g0(z=10.000)


# These keep track of where the most recent move left the controlled
# point, or None if the position is not known.
current_x = None
current_y = None
current_z = None
current_a = None
current_b = None
current_c = None
current_u = None
current_v = None
current_w = None


# When comparing floats, a difference of less than epsilon counts as no
# difference at all.
epsilon = 0.000000001

def close_enough(a, b):
    """Returns True if the two numbers `a` and `b` are within `epsilon`
    (1e-9) of each other, False if they're farther apart."""
    return abs(a - b) < epsilon

def complex_close_enough(a, b):
    """Returns True if the two complex numbers `a` and `b` are within
    `epsilon` (1e-9) of each other, False if they're farther apart."""
    diff = complex(a.real - b.real, a.imag - b.imag)
    mag = math.sqrt(pow(diff.real, 2) + pow(diff.imag, 2))
    if mag < epsilon:
        return True
    return False


def init():
    print()
    print("; init")
    print("G20          (inch)")
    print("G17          (xy plane)")
    print("G90          (absolute)")
    print("G91.1        (arc centers are relative to arc starting point)")
    cutter_comp_off()
    print("G54          (switch to coordinate system 1)")
    print("G94          (units/minute feed mode)")
    print("G99          (in canned cycles, retract to the Z coordinate specified by the R word)")
    print("G64 P0.0005  (enable path blending, but stay within 0.0005 of the programmed path)")
    print("G49          (turn off tool length compensation)")
    print("G80          (turn off canned cycles)")
    print()


def comment(msg):
    if msg:
        print(";", msg)
    else:
        print()


def absolute():
    print("G90")


def absolute_arc_centers():
    print("G90.1")


def relative_arc_centers():
    print("G91.1")


def spindle_on():
    print("M3")


def spindle_off():
    print("M5")


def quill_up():
    absolute()
    cutter_comp_off()
    print("G53 G0 Z0")
    current_z = None
    spindle_off()


def presentation_position():
    imperial()
    quill_up()

    # rapid to presentation position
    # table centered in X, all the way forward towards the user
    print("G53 G0 X9 Y12")
    current_x = None
    current_y = None


def m2():
    print()
    print("M2")


def done():
    print()
    print("; done")
    presentation_position()
    print("M2")


def imperial():
    print("G20")


def metric():
    print("G21")


def feed(feed_rate_units_per_minute):
    print("F %.4f" % feed_rate_units_per_minute)


def speed(spindle_rpm):
    print("S %d" % spindle_rpm)


# FIXME: g0(path) should be merged or replaced by z_path() somehow
def g0(path=None, x=None, y=None, z=None, a=None, b=None, c=None, u=None, v=None, w=None):
    global current_x
    global current_y
    global current_z
    global current_a
    global current_b
    global current_c
    global current_u
    global current_v
    global current_w

    if path is not None:
        print()
        print("; g0 path")
        for waypoint in path:
            g0(**waypoint)
        print()
    else:
        print("G0", end='')
        if x is not None:
            current_x = x
            print(" X%.4f" % x, end='')
        if y is not None:
            current_y = y
            print(" Y%.4f" % y, end='')
        if z is not None:
            current_z = z
            print(" Z%.4f" % z, end='')
        if a is not None:
            current_a = a
            print(" A%.4f" % a, end='')
        if b is not None:
            current_b = b
            print(" B%.4f" % b, end='')
        if c is not None:
            current_c = c
            print(" C%.4f" % c, end='')
        if u is not None:
            current_u = u
            print(" U%.4f" % u, end='')
        if v is not None:
            current_v = v
            print(" V%.4f" % v, end='')
        if w is not None:
            current_w = w
            print(" W%.4f" % w, end='')
        print()


def g1(path=None, x=None, y=None, z=None, a=None, b=None, c=None, u=None, v=None, w=None):
    global current_x
    global current_y
    global current_z
    global current_a
    global current_b
    global current_c
    global current_u
    global current_v
    global current_w

    if path is not None:
        print()
        print("; g1 path")
        for waypoint in path:
            g1(**waypoint)
        print()
    else:
        print("G1", end='')
        if x is not None:
            current_x = x
            print(" X%.4f" % x, end='')
        if y is not None:
            current_y = y
            print(" Y%.4f" % y, end='')
        if z is not None:
            current_z = z
            print(" Z%.4f" % z, end='')
        if a is not None:
            current_a = a
            print(" A%.4f" % a, end='')
        if b is not None:
            current_b = b
            print(" B%.4f" % b, end='')
        if c is not None:
            current_c = c
            print(" C%.4f" % c, end='')
        if u is not None:
            current_u = u
            print(" U%.4f" % u, end='')
        if v is not None:
            current_v = v
            print(" V%.4f" % v, end='')
        if w is not None:
            current_w = w
            print(" W%.4f" % w, end='')
        print()


def g2(x=None, y=None, z=None, i=None, j=None, p=None):
    global current_x
    global current_y
    global current_z

    """Clockwise arc feed."""
    if i is None and j is None:
        raise TypeError, "gcoder.g2() without i or j"
    print("G2", end='')
    if x is not None:
        current_x = x
        print(" X%.4f" % x, end='')
    if y is not None:
        current_y = y
        print(" Y%.4f" % y, end='')
    if z is not None:
        current_z = z
        print(" Z%.4f" % z, end='')
    if i is not None: print(" I%.4f" % i, end='')
    if j is not None: print(" J%.4f" % j, end='')
    if p is not None: print(" P%.4f" % p, end='')
    print()


def g3(x=None, y=None, z=None, i=None, j=None, p=None):
    global current_x
    global current_y
    global current_z

    """Counter-clockwise arc feed."""
    if i is None and j is None:
        raise TypeError, "gcoder.g3() without i or j"
    print("G3", end='')
    if x is not None:
        current_x = x
        print(" X%.4f" % x, end='')
    if y is not None:
        current_y = y
        print(" Y%.4f" % y, end='')
    if z is not None:
        current_z = z
        print(" Z%.4f" % z, end='')
    if i is not None: print(" I%.4f" % i, end='')
    if j is not None: print(" J%.4f" % j, end='')
    if p is not None: print(" P%.4f" % p, end='')
    print()


#
# Cutter compensation handling.
#

def cutter_comp_off():
    print("G40          (cutter comp off)")

def cancel_cutter_comp():
    print("; gcoder: calling program used obsolete cancel_cutter_comp() function, use cutter_comp_off() instead")
    cutter_comp_off()

def g40():
    print("; gcoder: calling program used obsolete g40() function, use cutter_comp_off() instead")
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
        print("G41.1 D%.4f   (cutter comp left, diameter mode)" % kwargs['diameter'])
    elif 'tool' in kwargs:
        print("G41 D%d   (cutter comp left, tool-number mode)" % kwargs['tool'])
    else:
        print("G41   (cutter comp left, current tool)")


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
        print("G42.1 D%.4f   (cutter comp right, diameter mode)" % kwargs['diameter'])
    elif 'tool' in kwargs:
        print("G42 D%d   (cutter comp right, tool-number mode)" % kwargs['tool'])
    else:
        print("G42   (cutter comp right, current tool)")

def g42_1(comp_diameter):
    print("; gcoder: calling program used obsolete g42_1() function, use cutter_comp_right() instead")
    cutter_comp_right(diameter=comp_diameter)


def g81(retract, x=None, y=None, z=None):
    global current_x
    global current_y
    global current_z

    print("G81", end='')
    if x is not None:
        current_x = x
        print(" X%.4f" % x, end='')
    if y is not None:
        current_y = y
        print(" Y%.4f" % y, end='')
    if z is not None:
        print(" Z%.4f" % z, end='')
    print(" R%.4f" % retract, end='')
    print()
    # FIXME: keep track of retract mode, set Z correctly here
    current_z = None


def g83(retract, delta, x=None, y=None, z=None):
    global current_x
    global current_y
    global current_z

    print("G83", end='')
    if x is not None:
        current_x = x
        print(" X%.4f" % x, end='')
    if y is not None:
        current_y = y
        print(" Y%.4f" % y, end='')
    if z is not None:
        print(" Z%.4f" % z, end='')
    print(" R%.4f" % retract, end='')
    print(" Q%.4f" % delta, end='')
    print()
    # FIXME: keep track of retract mode, set Z correctly here
    current_z = None


def drill_hog(diameter, retract, delta, z_drill, x0, y0, x1, y1, xy_finishing_allowance=None, z_finishing_allowance=None):

    """Drills as many evenly spaced holes as will fit in a rectangular
    grid, within the rectangle defined by (x0, y0) and (x1, y1).
    The specified rectangle describes the material contour, the holes
    will be inset from the edges by the drill's radius.

    If finishing_tolerance is specified, all the holes will stay at least
    that far away from the material contour the specified rectangle.
    If z_finishing_allowance is specified the holes will end that far
    above the specified drill depth."""

    print()
    print("; drill hog")

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

            g83(x=x, y=y, z=z_drill, delta=delta, retract=retract)

    print()


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


def z_path2(path, depth_of_cut, z_target):

    """This function traverses a path (a list of Line, ArcCW, and ArcCCW
    objects), cutting a little deeper on each pass.

    z_path2() has a local variable named "Z" that tracks the Z level
    of the current pass.  On entry to the function it is initialized
    to gcoder's current Z position.  The Z variable overrides any Z
    coordinates specified in the path.

    The motion is this:

        1. Initialize the local variable Z to gcoder's current Z
           coordinate.

        2. Reduce Z by depth_of_cut, but not below z_target.

        3. Feed to each waypoint in path, in order starting with the
           first and ending with the last.

        4. After arriving at the last waypoint, if Z is not yet down to
           z_target, go to step 2.

        5. After reaching step 4 with Z at z_target, feed back to the
           first waypoint while keeping Z at the z_target level, thereby
           cutting away the ramp left by the previous iteration.

    The first move of this function is to the first waypoint in the path,
    at a Z level that's depth-of-cut below the starting Z level (but not
    below z_target).  If you position the cutter at the *last* waypoint
    in the path, at the start-of-material Z level, you'll get a nice
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

    z = current_z

    # Shrink depth_of_cut so all passes are equally deep, instead of
    # letting the final pass be "whatever's left over".
    z_range = current_z - z_target
    num_passes = math.ceil(float(z_range) / depth_of_cut)
    depth_of_cut = z_range / num_passes

    while not close_enough(z, z_target):
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

    global epsilon

    z_range = z_start - z_end
    num_passes = math.ceil(z_range / max_doc)
    doc = z_range / num_passes

    g0(x=x_start, y=y_start)

    spindle_on();

    z = z_start
    g0(z=z)

    while not close_enough(z, z_end):
        z = z - doc
        if z < z_end:
            z = z_end
        if rapid_plunge:
            g0(z=z)
        else:
            g1(z=z)
        g1(x=x_end, y=y_end)

        if not close_enough(z, z_end):
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

