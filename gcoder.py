#
# gcoder.py - python library for writing g-code
#
# Copyright (C) 2018-2020 Sebastian Kuzminsky
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 


from __future__ import print_function

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
    GCODE_ORIGIN_IS_VIEWBOX_LOWER_LEFT=0
    GCODE_ORIGIN_IS_SVG_ORIGIN=1

    def __init__(self, svg_file, gcode_origin=GCODE_ORIGIN_IS_VIEWBOX_LOWER_LEFT):
        self.svg_file = svg_file
        self.gcode_origin = gcode_origin

        print("gcode_origin:", gcode_origin, file=sys.stderr)

        #
        # We need to convert from whatever units the SVG input file
        # (self) is in, to mm for the output gcode.
        #
        # The SVG spec lets SVG files program paths in arbitrary "user"
        # units.  The conversion from user units to mm is done as follows:
        #
        # The SVG file specifies its width and height in one of the
        # accepted "real-world" units, all of which can be converted
        # to mm with a straight-forward linear scale (see the viewport
        # handling below).
        #
        # The SVG file *may* specify a viewBox.
        #
        # If no viewBox is specified, the viewport units are used for
        # the user units.
        #
        # If a viewBox is specified, it provides the (X, Y) coordinates
        # (in user units) of the lower left corner of the viewport and
        # the width and height (again in user units) of the viewport.
        #
        # From all this we compute an x_scale and y_scale, such that:
        #
        # x_gcode(mm) = x_svg * x_scale
        # y_gcode(mm) = y_svg * y_scale
        #
        # (This ignores the fact that SVG specifies Y upside-down compared
        # to gcode, and that the viewBox can be offset from the origin,
        # see svg.x_to_mm() and svg.y_to_mm() below.)
        #

        self.doc = svgpathtools.Document(self.svg_file)

        self.results = self.doc.flatten_all_paths()
        self.paths = [result.path for result in self.results]

        self.svg_attributes = self.doc.root.attrib


        #
        # Deal with the viewport.
        #

        val, units, scale = self._parse_height_width(self.svg_attributes['height'])
        self.viewport_height = val
        self.viewport_y_units = units
        self.x_scale = scale

        val, units, scale = self._parse_height_width(self.svg_attributes['width'])
        self.viewport_width = val
        self.viewport_x_units = units
        self.y_scale = scale

        print("svg viewport:", file=sys.stderr)
        print("    width: %.3f%s" % (self.viewport_width, self.viewport_x_units), file=sys.stderr)
        print("    height: %.3f%s" % (self.viewport_height, self.viewport_y_units), file=sys.stderr)


        #
        # Deal with the viewBox.
        #

        if 'viewBox' in self.svg_attributes:
            print("svg viewBox:", self.svg_attributes['viewBox'], file=sys.stderr)
            (x_min, y_min, width, height) = re.split(',|(?: +(?:, *)?)', self.svg_attributes['viewBox'])
            self.viewBox_x = float(x_min)
            self.viewBox_y = float(y_min)
            self.viewBox_width = float(width)
            self.viewBox_height = float(height)

            self.x_scale *= self.viewport_width / self.viewBox_width
            self.y_scale *= self.viewport_height / self.viewBox_height

            print("    viewBox_x:", self.viewBox_x, file=sys.stderr)
            print("    viewBox_y:", self.viewBox_y, file=sys.stderr)
            print("    viewBox_width:", self.viewBox_width, file=sys.stderr)
            print("    viewBox_height:", self.viewBox_height, file=sys.stderr)
            print("    x_scale:", self.x_scale, file=sys.stderr)
            print("    y_scale:", self.y_scale, file=sys.stderr)

        else:
            self.viewBox_x = 0.0
            self.viewBox_y = 0.0
            self.viewBox_width = self.viewport_width
            self.viewBox_height = self.viewport_height


    def _parse_height_width(self, s):
        m = re.match('^([0-9.]+)([a-zA-Z]*)$', s)
        if m == None or len(m.groups()) != 2:
            raise SystemExit("failed to parse SVG viewport height/width: %s" % s)

        val = float(m.group(1))
        unit = m.group(2)

        units = {
            # "px" (or "no units") is 96 dpi: 1 inch/96 px * 25.4 mm/1 inch = 25.4/96 mm/px
            '': 25.4/96,
            'px': 25.4/96,

            # "pt" is 72 dpi: 1 inch/72 pt * 25.4 mm/1 inch = 25.4/72 mm/pt
            'pt': 25.4/72,

            # Units are Picas "pc", 6 dpi: 1 inch/6 pc * 25.4 mm/1 inch = 25.4/6 mm/pc
            'pc': 25.4/6,

            'cm': 10.0,
            'mm': 1.0,

            # Units are inches: 25.4 mm/1 inch
            'in': 25.4
        }

        if unit not in units:
            raise SystemExit("unknwn SVG viewport units: '%s'" % unit)

        scale = units[unit]

        return (val, unit, scale)


    def x_to_mm(self, x):
        if type(x) != float:
            raise SystemExit('non-float input')
        if self.gcode_origin == self.GCODE_ORIGIN_IS_SVG_ORIGIN:
            out = x * self.x_scale
        elif self.gcode_origin == self.GCODE_ORIGIN_IS_VIEWBOX_LOWER_LEFT:
            out = (x - self.viewBox_x) * self.x_scale
        return out


    def y_to_mm(self, y):
        if type(y) != float:
            raise SystemExit('non-float input')
        # Y is upside down in SVG.
        if self.gcode_origin == self.GCODE_ORIGIN_IS_SVG_ORIGIN:
            out = -y * self.y_scale
        elif self.gcode_origin == self.GCODE_ORIGIN_IS_VIEWBOX_LOWER_LEFT:
            out = (self.viewBox_height - (y - self.viewBox_y)) * self.y_scale
        return out


    def xy_to_mm(self, xy):
        if type(xy) != complex:
            raise SystemExit('non-complex input')
        x = self.x_to_mm(xy.real)
        y = self.y_to_mm(xy.imag)
        return (x, y)


def close_path(p):
    def join_segments(this_seg, next_seg):
        if this_seg.end != next_seg.start:
            if close_enough(this_seg.end, next_seg.start):
                avg = (this_seg.end + next_seg.start) / 2.0
                this_seg.end = avg
                next_seg.start = avg
            else:
                raise ValueError("segments are not even close to closed: %s, %s" % (this_seg, next_seg))
        return (this_seg, next_seg)

    for i in range(len(p)-1):
        this_seg = p[i]
        next_seg = p[i+1]
        (this_seg, next_seg) = join_segments(this_seg, next_seg)
        p[i] = this_seg
        p[i+1] = next_seg

    this_seg = p[-1]
    next_seg = p[0]
    (this_seg, next_seg) = join_segments(this_seg, next_seg)
    p[-1] = this_seg
    p[0] = next_seg

    p.closed = True

    return p


def split_path_at_intersections(path_list, debug=False):

    """`path_list` is a list of connected path segments, or a
    svgpathtools.path.Path() object.  This function identifies
    each place where the path intersects intself, and splits each
    non-self-intersecting subset of the path into a separate path list.
    This may involve splitting segments.

    Returns a list of path lists."""

    # If path_list is a Path object, convert it to a regular list (of
    # segments), because it's easier to work with.
    if type(path_list) == svgpathtools.path.Path:
        l = []
        for i in range(len(path_list)):
            l.append(path_list[i])
        path_list = l

    def find_earliest_intersection(path_list, this_seg_index):
        this_seg = path_list[this_seg_index]
        if debug: print("looking for earliest intersection of this seg(%d):" % this_seg_index, this_seg, file=sys.stderr)

        earliest_this_t = None
        earliest_other_seg_index = None
        earliest_other_t = None

        for other_seg_index in range(this_seg_index+2, len(path_list)):
            other_seg = path_list[other_seg_index]
            if debug: print("    other[%d]:" % other_seg_index, other_seg, file=sys.stderr)
            intersections = this_seg.intersect(other_seg)
            if len(intersections) == 0:
                continue
            if debug: print("        intersect!  %s" % intersections, file=sys.stderr)

            # The intersection that comes earliest in `this_seg` is
            # the interesting one, except that intersections at the
            # segments' endpoints don't count.
            for intersection in intersections:
                if close_enough(intersection[0], 0.0) or close_enough(intersection[0], 1.0):
                    if debug: print("            at end of this segment, ignoring", file=sys.stderr)
                    continue

                if intersection[0] > 1.0 or intersection[0] < 0.0:
                    if debug: print("            off the end of this segment?!  ignoring", file=sys.stderr)
                    continue

                if (earliest_this_t == None) or (intersection[0] < earliest_this_t):
                    if debug: print("        earliest!", file=sys.stderr)
                    earliest_this_t = intersection[0]
                    earliest_other_seg_index = other_seg_index
                    earliest_other_t = intersection[1]

        return earliest_this_t, earliest_other_seg_index, earliest_other_t


    if debug: print("splitting path:", file=sys.stderr)
    if debug: print("    ", path_list, file=sys.stderr)

    # This is a list of pairs.  Each pair represents a place where the
    # input path crosses itself.  The two members of the pair are the
    # indexes of the segments that end at the intersection point.
    intersections = []

    this_seg_index = 0
    while this_seg_index < len(path_list):
        this_seg = path_list[this_seg_index]

        this_t, other_seg_index, other_t = find_earliest_intersection(path_list, this_seg_index)
        if this_t == None:
            this_seg_index += 1
            continue
        other_seg = path_list[other_seg_index]

        # Found the next intersection.  Split the segments and note
        # the intersection.

        if debug: print("intersection:", file=sys.stderr)
        if debug: print("    this:", file=sys.stderr)
        if debug: print("        %d @ %f" % (this_seg_index, this_t), file=sys.stderr)
        if debug: print("        %s" % this_seg, file=sys.stderr)
        if debug: print("    other:", file=sys.stderr)
        if debug: print("        %d @ %f" % (other_seg_index, other_t), file=sys.stderr)
        if debug: print("        %s" % other_seg, file=sys.stderr)

        this_first_seg, this_second_seg = this_seg.split(this_t)
        other_first_seg, other_second_seg = other_seg.split(other_t)
        if debug: print("split this seg:", this_seg, file=sys.stderr)
        if debug: print("    t:", this_t, file=sys.stderr)
        if debug: print("    ", this_first_seg, file=sys.stderr)
        if debug: print("    ", this_second_seg, file=sys.stderr)
        if debug: print("split other seg:", other_seg, file=sys.stderr)
        if debug: print("    t:", other_t, file=sys.stderr)
        if debug: print("    ", other_first_seg, file=sys.stderr)
        if debug: print("    ", other_second_seg, file=sys.stderr)

        # FIXME: This fixup is bogus, but the two segments'
        # `t` parameters don't put the intersection at the
        # same point...
        other_first_seg.end = this_first_seg.end
        other_second_seg.start = other_first_seg.end

        assert(close_enough(this_first_seg.end, this_second_seg.start))
        assert(close_enough(this_first_seg.end, other_first_seg.end))
        assert(close_enough(this_first_seg.end, other_second_seg.start))

        assert(close_enough(this_first_seg.start, this_seg.start))
        assert(close_enough(this_second_seg.end, this_seg.end))

        assert(close_enough(other_first_seg.start, other_seg.start))
        assert(close_enough(other_second_seg.end, other_seg.end))

        # Replace the old (pre-split) this_seg with the first sub-segment.
        path_list[this_seg_index] = this_first_seg

        # Insert the second sub-segment after the first one.
        path_list.insert(this_seg_index+1, this_second_seg)

        # We inserted a segment before other_seg, so we increment
        # its index.
        other_seg_index += 1

        # Replace the old (pre-split) other_seg with the first sub-segment.
        path_list[other_seg_index] = other_first_seg

        # Insert the second sub-segment after the first one.
        path_list.insert(other_seg_index+1, other_second_seg)

        for i in range(len(intersections)):
            if debug: print("bumping intersection:", file=sys.stderr)
            if debug: print("    ", intersections[i], file=sys.stderr)
            if intersections[i][1] >= this_seg_index:
                intersections[i][1] += 1  # for this_seg that got split
            if intersections[i][1] >= other_seg_index:
                intersections[i][1] += 1  # for other_seg that got split
            if debug: print("    ", intersections[i], file=sys.stderr)

        # Add this new intersection we just made.
        i = [this_seg_index, other_seg_index]
        if debug: print("    new:", i, file=sys.stderr)
        intersections.append(i)

        # Look for intersections in the remainder of this_seg (the second
        # part of the split).
        this_seg_index += 1

    if len(intersections) > 0:
        if debug: print("found some intersections:", file=sys.stderr)
        for i in intersections:
            if debug: print("    ", i, file=sys.stderr)
            if debug: print("        ", path_list[i[0]], file=sys.stderr)
            if debug: print("        ", path_list[i[1]], file=sys.stderr)
    else:
        if debug: print("path does not self-intersect", file=sys.stderr)

    paths = []
    while True:
        if debug: print("starting a new path", file=sys.stderr)
        path = []
        # Start at the first unused segment
        seg_index = 0
        for seg_index in range(len(path_list)):
            if path_list[seg_index] != None:
                break

        while seg_index < len(path_list):
            if path_list[seg_index] == None:
                # Done with this path.
                break

            if debug: print("    adding segment %d:" % seg_index, path_list[seg_index], file=sys.stderr)
            path.append(path_list[seg_index])
            path_list[seg_index] = None

            i = None
            for i in intersections:
                if seg_index == i[0] or seg_index == i[1]:
                    break
            if debug: print("i:", i, file=sys.stderr)
            if debug: print("seg_index:", seg_index, file=sys.stderr)
            if (i is not None) and (i[0] == seg_index):
                # This segment is the first entrance to an intersection,
                # take the second exit.
                if debug: print("    intersection!", file=sys.stderr)
                seg_index = i[1] + 1
            elif (i is not None) and (i[1] == seg_index):
                # This segment is the second entrance to an intersection,
                # take the first exit.
                if debug: print("    intersection!", file=sys.stderr)
                seg_index = i[0] + 1
            else:
                # This segment doesn't end in an intersection, just go
                # to the next one.
                seg_index += 1

        if path == []:
            break

        paths.append(path)

    return paths


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


def offset_paths(path, offset_distance, steps=100, debug=False):
    """Takes an svgpathtools.path.Path object, `path`, and a float
    distance, `offset_distance`, and returns the parallel offset curves
    (in the form of a list of svgpathtools.path.Path objects)."""


    def is_enclosed(path, check_paths):

        """`path` is an svgpathtools.path.Path object, `check_paths`
        is a list of svgpath.path.Path objects.  This function returns
        True if `path` lies inside any of the paths in `check_paths`,
        and returns False if it lies outside all of them."""

        seg = path[0]
        point = seg.point(0.5)

        for i in range(len(check_paths)):
            test_path = check_paths[i]
            if path == test_path:
                continue
            # find outside_point, which lies outside other_path
            (xmin, xmax, ymin, ymax) = test_path.bbox()
            outside_point = complex(xmax+100, ymax+100)
            if svgpathtools.path_encloses_pt(point, outside_point, test_path):
                if debug: print("point is within path", i, file=sys.stderr)
                return True
        return False


    def intersect(this_seg, next_seg, intersection):
        this_point = this_seg.point(intersection[0])
        next_point = next_seg.point(intersection[1])

        if debug:
            print(f"    this_seg: {this_seg}", file=sys.stderr)
            print(f"    next_seg: {next_seg}", file=sys.stderr)
            print(f"    intersection: {intersection} {this_point} {next_point}", file=sys.stderr)

        if close_enough(intersection[0], 0.0) and close_enough(intersection[1], 1.0):
            # Start of `this_seg` touches end of `next_seg`, that's ok.
            point = (this_point + next_point) / 2
            this_seg.start = point
            next_seg.end = point
            # If you change an Arc you have to re-parameterize it.
            if type(this_seg) is svgpathtools.path.Arc:
                this_seg._parameterize()
            if type(next_seg) is svgpathtools.path.Arc:
                next_seg._parameterize()

        elif close_enough(intersection[0], 1.0) and close_enough(intersection[1], 0.0):
            # End of `this_seg` touches start of `next_seg`, that's ok.
            point = (this_point + next_point) / 2
            this_seg.end = point
            next_seg.start = point
            # If you change an Arc you have to re-parameterize it.
            if type(this_seg) is svgpathtools.path.Arc:
                this_seg._parameterize()
            if type(next_seg) is svgpathtools.path.Arc:
                next_seg._parameterize()

        else:
            # Trim the end off `this_seg` and the start off `next_seg`
            # so they meet at their intersection.
            this_seg = this_seg.cropped(0.0, intersection[0])
            next_seg = next_seg.cropped(intersection[1], 1.0)


    def remove_connected_zero_length_segments(path_list):
        if debug: print("removing too-short connected segments...", file=sys.stderr)
        new_path_list = []
        for i in range(len(path_list)):
            prev_seg = path_list[i-1]
            this_seg = path_list[i]
            next_seg = path_list[(i+1) % len(path_list)]
            if close_enough(prev_seg.end, this_seg.start) and close_enough(this_seg.end, next_seg.start) and this_seg.length() < epsilon:
                if debug: print(f"removing {this_seg.length()} long segment: {this_seg}", file=sys.stderr)
                midpoint = (prev_seg.end + next_seg.start) / 2
                prev_seg.end = midpoint
                next_seg.start = midpoint
                continue
            new_path_list.append(this_seg)
        return new_path_list


    # This only works on closed paths.
    if debug: print("input path:", file=sys.stderr)
    if debug: print(path, file=sys.stderr)
    if debug: print("offset:", offset_distance, file=sys.stderr)
    assert(path.isclosed())


    #
    # First generate a list of Path elements (Lines and Arcs),
    # corresponding to the offset versions of the Path elements in the
    # input path.
    #

    if debug: print("generating offset segments...", file=sys.stderr)

    offset_path_list = []
    for seg in path:
        if type(seg) == svgpathtools.path.Line:
            if close_enough(seg.point(0), seg.point(1)):
                if debug: print("    skipping zero-length line segment", file=sys.stderr)
                continue
            start = seg.point(0) + (offset_distance * seg.normal(0))
            end = seg.point(1) + (offset_distance * seg.normal(1))
            seg = svgpathtools.Line(start, end)
            offset_path_list.append(seg)
            if debug: print("    %s" % offset_path_list[-1], file=sys.stderr)

        elif type(seg) == svgpathtools.path.Arc and (seg.radius.real == seg.radius.imag):
            # Circular arcs remain arcs, elliptical arcs become linear
            # approximations below.
            #
            # Polygons (input paths) are counter-clockwise.
            #
            # Positive offsets are to the inside of the polygon, negative
            # offsets are to the outside.
            #
            # If this arc is counter-clockwise (sweep == False),
            # *subtract* the `offset_distance` from its radius, so
            # insetting makes the arc smaller and outsetting makes
            # it larger.
            #
            # If this arc is clockwise (sweep == True), *add* the
            # `offset_distance` from its radius, so insetting makes the
            # arc larger and outsetting makes it smaller.
            #
            # If the radius of the offset arc is negative, use its
            # absolute value and invert the sweep.

            if seg.sweep == False:
                new_radius = seg.radius.real - offset_distance
            else:
                new_radius = seg.radius.real + offset_distance

            start = seg.start + (offset_distance * seg.normal(0))
            end = seg.end + (offset_distance * seg.normal(1))
            sweep = seg.sweep

            flipped = False
            if new_radius < 0.0:
                if debug: print("    inverting Arc!", file=sys.stderr)
                flipped = True
                new_radius = abs(new_radius)
                sweep = not sweep

            if new_radius > minimum_arc_radius:
                radius = complex(new_radius, new_radius)
                offset_arc = svgpathtools.path.Arc(
                    start = start,
                    end = end,
                    radius = radius,
                    rotation = seg.rotation,
                    large_arc = seg.large_arc,
                    sweep = sweep
                )
            else:
                # Offset Arc radius is smaller than the minimum that
                # LinuxCNC accepts, replace with a Line.
                if debug: print("    arc too small, replacing with a line", file=sys.stderr)
                if flipped:
                    old_start = start
                    start = end
                    end = old_start
                offset_arc = svgpathtools.path.Line(start = start, end = end)
            offset_path_list.append(offset_arc)
            if debug: print("    %s" % offset_path_list[-1], file=sys.stderr)

        else:
            # Deal with any segment that's not a line or a circular arc.
            # This includes elliptic arcs and bezier curves.  Use linear
            # approximation.
            #
            # FIXME: Steps should probably be computed dynamically to make
            #     the length of the *offset* line segments manageable.
            points = []
            for k in range(steps+1):
                t = k / float(steps)
                normal = seg.normal(t)
                offset_vector = offset_distance * normal
                points.append(seg.point(t) + offset_vector)
            for k in range(len(points)-1):
                start = points[k]
                end = points[k+1]
                seg = svgpathtools.Line(start, end)
                offset_path_list.append(seg)
            if debug: print("    (long list of short lines)", file=sys.stderr)

    offset_path_list = remove_connected_zero_length_segments(offset_path_list)


    #
    # Find all the places where one segment intersects the next, and
    # trim to the intersection.
    #

    if debug: print("trimming intersecting segments...", file=sys.stderr)

    if len(offset_path_list) == 2:
        i = 0
        this_seg = offset_path_list[i]

        next_i = 1
        next_seg = offset_path_list[next_i]

        if debug: print("intersecting 2-segment path", file=sys.stderr)
        if debug: print("    this", this_seg, file=sys.stderr)
        if debug: print("        length", this_seg.length(), file=sys.stderr)
        if debug: print("    next", next_seg, file=sys.stderr)
        if debug: print("        length", next_seg.length(), file=sys.stderr)
        intersections = this_seg.intersect(next_seg)
        if debug: print("    intersections:", intersections, file=sys.stderr)
        for intersection in intersections:
            intersect(this_seg, next_seg, intersection)
            offset_path_list[i] = this_seg
            offset_path_list[next_i] = next_seg

            if debug: print("    trimmed:", file=sys.stderr)
            if debug: print("        this", this_seg, file=sys.stderr)
            if debug: print("            length", this_seg.length(), file=sys.stderr)
            if debug: print("        next", next_seg, file=sys.stderr)
            if debug: print("            length", next_seg.length(), file=sys.stderr)

    else:
        for i in range(len(offset_path_list)):
            this_seg = offset_path_list[i]

            next_i = (i + 1) % len(offset_path_list)
            next_seg = offset_path_list[next_i]

            if debug: print("intersecting", file=sys.stderr)
            if debug: print("    this", this_seg, file=sys.stderr)
            if debug: print("        length", this_seg.length(), file=sys.stderr)
            if debug: print("    next", next_seg, file=sys.stderr)
            if debug: print("        length", next_seg.length(), file=sys.stderr)
            intersections = this_seg.intersect(next_seg)
            if debug: print("    intersections:", intersections, file=sys.stderr)
            if len(intersections) > 0:
                intersection = intersections[0]
                if debug:
                    this_point = this_seg.point(intersections[0][0])
                    next_point = next_seg.point(intersections[0][1])
                    print("    first intersection: {} {} {}".format(intersection, this_point, next_point), file=sys.stderr)
                # Trim the end off `this_seg` and the start off `next_seg`
                # so they meet at their intersection.
                this_seg = this_seg.cropped(0.0, intersection[0])
                next_seg = next_seg.cropped(intersection[1], 1.0)

                offset_path_list[i] = this_seg
                offset_path_list[next_i] = next_seg

                if debug: print("    trimmed:", file=sys.stderr)
                if debug: print("        this", this_seg, file=sys.stderr)
                if debug: print("            length", this_seg.length(), file=sys.stderr)
                if debug: print("        next", next_seg, file=sys.stderr)
                if debug: print("            length", next_seg.length(), file=sys.stderr)

    offset_path_list = remove_connected_zero_length_segments(offset_path_list)


    #
    # Find all the places where adjacent segments do not end/start close
    # to each other, and join them with Arcs.
    #

    if debug: print("joining non-connecting segments with arcs...", file=sys.stderr)

    joined_offset_path_list = []
    for i in range(len(offset_path_list)):
        this_seg = offset_path_list[i]
        if (i+1) < len(offset_path_list):
            next_seg = offset_path_list[i+1]
        else:
            next_seg = offset_path_list[0]

        if close_enough(this_seg.end, next_seg.start):
            joined_offset_path_list.append(this_seg)
            continue

        if debug: print("these segments don't touch end to end:", file=sys.stderr)
        if debug: print("    this", this_seg, file=sys.stderr)
        if debug: print("    next", next_seg, file=sys.stderr)
        if debug: print("    error: %s (%.7f)" % (this_seg.end-next_seg.start, abs(this_seg.end-next_seg.start)), file=sys.stderr)

        # FIXME: Choose values for `large_arc` and `sweep` correctly here.
        # I think the goal is to make the joining arc tangent to the segments it joins.
        # large_arc should always be False
        # sweep means "clockwise" (but +Y is down)
        if debug: print("determining joining arc:", file=sys.stderr)
        if debug: print("    this_seg ending normal:", this_seg.normal(1), file=sys.stderr)
        if debug: print("    next_seg starting normal:", next_seg.normal(0), file=sys.stderr)

        sweep_arc = svgpathtools.path.Arc(
            start = this_seg.end,
            end = next_seg.start,
            radius = complex(offset_distance, offset_distance),
            rotation = 0,
            large_arc = False,
            sweep = True
        )
        sweep_start_error = this_seg.normal(1) - sweep_arc.normal(0)
        sweep_end_error = next_seg.normal(0) - sweep_arc.normal(1)
        sweep_error = pow(abs(sweep_start_error), 2) + pow(abs(sweep_end_error), 2)
        if debug: print("    sweep arc starting normal:", sweep_arc.normal(0), file=sys.stderr)
        if debug: print("    sweep arc ending normal:", sweep_arc.normal(1), file=sys.stderr)
        if debug: print("    sweep starting error:", sweep_start_error, file=sys.stderr)
        if debug: print("    sweep end error:", sweep_end_error, file=sys.stderr)
        if debug: print("    sweep error:", sweep_error, file=sys.stderr)

        antisweep_arc = svgpathtools.path.Arc(
            start = this_seg.end,
            end = next_seg.start,
            radius = complex(offset_distance, offset_distance),
            rotation = 0,
            large_arc = False,
            sweep = False
        )
        antisweep_start_error = this_seg.normal(1) - antisweep_arc.normal(0)
        antisweep_end_error = next_seg.normal(0) - antisweep_arc.normal(1)
        antisweep_error = pow(abs(antisweep_start_error), 2) + pow(abs(antisweep_end_error), 2)
        if debug: print("    antisweep arc starting normal:", antisweep_arc.normal(0), file=sys.stderr)
        if debug: print("    antisweep arc ending normal:", antisweep_arc.normal(1), file=sys.stderr)
        if debug: print("    antisweep starting error:", antisweep_start_error, file=sys.stderr)
        if debug: print("    antisweep end error:", antisweep_end_error, file=sys.stderr)
        if debug: print("    antisweep error:", antisweep_error, file=sys.stderr)

        joining_arc = None
        if sweep_error < antisweep_error:
            if debug: print("joining arc is sweep", file=sys.stderr)
            joining_arc = sweep_arc
        else:
            if debug: print("joining arc is antisweep", file=sys.stderr)
            joining_arc = antisweep_arc

        if debug: print("joining arc:", file=sys.stderr)
        if debug: print(joining_arc, file=sys.stderr)
        if debug: print("    length:", joining_arc.length(), file=sys.stderr)
        if debug: print("    start-end distance:", joining_arc.start-joining_arc.end, file=sys.stderr)

        # FIXME: this is kind of arbitrary
        # FIXME: we should really just drop any segment that doesn't move
        #     the controlled point enough that the current coordinates
        #     change.  That's currently 0.0001mm, but should probably
        #     be configurable.
        joining_seg = joining_arc
        if joining_arc.length() < 1e-4:
            joining_seg = svgpathtools.path.Line(joining_arc.start, joining_arc.end)
            if debug: print("    too short!  replacing with a line:", joining_seg, file=sys.stderr)

        joined_offset_path_list.append(this_seg)
        joined_offset_path_list.append(joining_seg)

    offset_path_list = joined_offset_path_list


    #
    # Find the places where the path intersects itself, split into
    # multiple separate paths in those places.
    #

    if debug: print("splitting path at intersections...", file=sys.stderr)

    offset_paths_list = split_path_at_intersections(offset_path_list, debug=debug)

    new_offset_paths_list = []
    for path_list in offset_paths_list:
        new_path_list = remove_connected_zero_length_segments(path_list)
        new_offset_paths_list.append(new_path_list)

    offset_paths_list = new_offset_paths_list


    #
    # Smooth the path: adjacent segments whose start/end points are
    # "close enough" to each other are adjusted so they actually touch.
    #

    if debug: print("smoothing paths...", file=sys.stderr)

    for path_list in offset_paths_list:
        for i in range(len(path_list)):
            this_seg = path_list[i]
            next_seg = path_list[(i+1) % len(path_list)]
            midpoint = (this_seg.end + next_seg.start) / 2
            next_seg.start = this_seg.end = midpoint


    #
    # Convert each path list to a Path object and sanity check.
    #

    if debug: print("converting path lists to paths...", file=sys.stderr)

    offset_paths = []
    for path_list in offset_paths_list:
        offset_path = svgpathtools.Path(*path_list)
        if debug: print("offset path:", file=sys.stderr)
        if debug: print(offset_path, file=sys.stderr)
        assert(offset_path.isclosed())
        offset_paths.append(offset_path)


    #
    # The set of paths we got from split_path_at_intersections() has
    # zero or more 'true paths' that we actually want to return, plus
    # zero or more 'false paths' that should be discarded.
    #
    # When offsetting a path to the inside, the false paths will be
    # outside the true path and will wind in the opposite direction of
    # the input path.
    #
    # When offsetting a path to the outside, the false paths will be
    # inside the true paths, and will wind in the same direction as the
    # input path.
    #
    # [citation needed]
    #

    if debug: print("pruning false paths...", file=sys.stderr)

    path_area = approximate_path_area(path)
    if debug: print("input path area:", path_area, file=sys.stderr)

    keepers = []

    if offset_distance > 0:
        # The offset is positive (inwards), discard paths with opposite
        # direction from input path, and paths inside any other path.
        for offset_path in offset_paths:
            if debug: print("checking path:", offset_path, file=sys.stderr)
            offset_path_area = approximate_path_area(offset_path)
            if debug: print("offset path area:", offset_path_area, file=sys.stderr)
            if is_enclosed(offset_path, offset_paths):
                if debug: print("path is enclosed, dropping", file=sys.stderr)
            elif path_area * offset_path_area < 0.0:
                # Input path and offset path go in the opposite directions,
                # drop offset path.
                if debug: print("wrong direction, dropping", file=sys.stderr)
            else:
                if debug: print(f"good direction, keeping", file=sys.stderr)
                keepers.append(offset_path)

    else:
        # The offset is negative (outwards), discard paths that lie
        # inside any other path and have the same winding direction as
        # the input path.
        for offset_path in offset_paths:
            if debug: print("checking path:", offset_path, file=sys.stderr)
            if is_enclosed(offset_path, offset_paths):
                if debug: print("    enclosed", file=sys.stderr)
                # This path is enclosed, check the winding direction.
                offset_path_area = approximate_path_area(offset_path)
                if debug: print("offset path area:", offset_path_area, file=sys.stderr)
                if path_area * offset_path_area > 0.0:
                    if debug: print("    winding is the same as input, dropping", file=sys.stderr)
                    continue
                else:
                    if debug: print("    winding is opposite input", file=sys.stderr)
            else:
                if debug: print("    not enclosed", file=sys.stderr)
            if debug: print("    keeping", file=sys.stderr)
            keepers.append(offset_path)

    offset_paths = keepers

    return offset_paths



def path_segment_to_gcode(svg, segment, z=None):
    if type(segment) == svgpathtools.path.Line:
        (start_x, start_y) = svg.xy_to_mm(segment.start)
        (end_x, end_y) = svg.xy_to_mm(segment.end)
        g1(x=end_x, y=end_y, z=z)
    elif type(segment) is svgpathtools.path.Arc and (segment.radius.real == segment.radius.imag):
        # FIXME: g90.1 or g91.1?
        (end_x, end_y) = svg.xy_to_mm(segment.end)
        (center_x, center_y) = svg.xy_to_mm(segment.center)

        # In SVG sweep==True means clockwise and sweep==False means
        # counter-clockwise, but in svgpathtools it's the opposite.
        if segment.sweep:
            g2(x=end_x, y=end_y, z=z, i=center_x, j=center_y)
        else:
            g3(x=end_x, y=end_y, z=z, i=center_x, j=center_y)
    else:
        # Deal with any segment that's not a line or a circular arc,
        # this includes elliptic arcs and bezier curves.  Use linear
        # approximation.
        #
        # FIXME: The number of steps should probably be dynamically
        #     adjusted to make the length of the *offset* line
        #     segments manageable.
        if z is not None:
            raise ValueError("Z value specified for non-Line, non-circular-Arc segment: %s", segment)
        steps = 1000
        for k in range(steps+1):
            t = k / float(steps)
            end = segment.point(t)
            (end_x, end_y) = svg.xy_to_mm(end)
            g1(x=end_x, y=end_y)


def path_to_gcode(svg, path, z_traverse=10, z_approach=None, z_top_of_material=0, z_cut_depth=0, lead_in=True, lead_out=True, feed=None, plunge_feed=None, ramp_slope=None, max_depth_of_cut=None, work_holding_tabs=0, work_holding_tab_height=0.5, work_holding_tab_width=10.0, work_holding_tab_locations=[], debug=False):

    """Prints the G-code corresponding to the input `path`.

Arguments:

    `path`: The SVG path to emit g-code for.

    `z_traverse` (float, default 10.0): Z level for safe traversals/rapids
        across the work piece.  The lead-in move (if enabled) goes here
        before any other motion.  The lead-out move (if enabled) leaves
        the tool up here.

    `z_approach` (float, defaults to 0.5mm above z_top_of_material):
        Above this Z level we rapid, below it we feed.

    `z_top_of_material` (float, default 0.0): Z level where the work
        piece material starts.

    `z_cut_depth` (float, default 0.0): Z level to cut down to.

    `lead_in` (boolean, default True): Perform preliminary motion
        (before starting to cut the path) as described below.

    `lead_out` (boolean, default True): Perform final motion (after
        finishing cutting the path) as described below.

    `feed`: Feed rate to use when ramping into the cut and when cutting
        along the path.

    `plunge_feed`: Feed rate to use when plunging into the cut.

    `ramp_slope`: The Z depth per unit path length of the ramp slope.

    `max_depth_of_cut` (float, optional): Max axial depth of cut.  The cut
        starts at `z_top_of_material` and goes down to `z_cut_depth`;
        if this distance is greater than `max_depth_of_cut` then the
        cut will consist of multiple passes.  The axial depth-of-cut of
        each pass will be the same, and may be adjusted down from the
        specified max slightly, to cut to the specified depth in the
        smallest number of equal-depth passes.  If `max_depth_of_cut` is
        not specified, then the full depth will be cut in a single pass.

    `work_holding_tabs` (integer, optional, default 0): Number of
        work-holding tabs to add.

    `work_holding_tab_height` (float, optional, default 1.0): Height
        (above `z_cut_depth`) of the work holding tabs (if any).

    `work_holding_tab_width` (float, optional, default 10.0): Width of
        each work-holding tab (if any).  Note that this function does
        not take into account the cutter diameter, so the width of the
        actual material of the work-holding tabs will be one cutter
        diameter smaller than this number.

    `work_holding_tab_locations` (list of float, optional): Starting
        locations of the work holding tabs, specified as distances
        along the path.  The length of the list must be the same as
        `work_holding_tabs`.  If omitted, the first tab starts at 0.0
        and all the tabs are evenly spaced along the path.

    Motion:

        This function produces motion in three parts: Preliminary Motion,
        one or more Cutting Passes, and Final Motion.

        Preliminary Motion:

            This function begins with preliminary motion, controlled by the
            `lead_in` argument.  If `lead_in` is True:

                * rapid Z to the `z_traverse` level

                * rapid X and Y to the start of the first path segment

                * rapid Z to the `z_approach` level

            If `lead_in` is False:

                * feed X and Y to the start of the first segment at
                  `feed` feed rate

            At the completion of Preliminary Motion, the tool is
            positioned at the X and Y coordinates of the start of the
            first segment in the path, but at an unknown Z level.

            The Z level may be `z_approach` (if the `lead_in` argument
            is True) or Z may be unchanged from where it was on entry
            to the function (if the `lead_in` argument is False).

        Cutting Pass:

            The cut (from `z_top_of_material` to `z_cut_depth`) is made
            in the minimum number of equal-depth cutting passes.

            If `max_depth_of_cut` is None (not specified) the cut is
            made in a single full-depth pass.

            If `max_depth_of_cut` has a value, the cut is made in
            the smallest possible number of equal-depth passes, where
            the per-pass depth of cut is as large as possible without
            exceeding `max_depth_of_cut`.

            If `ramp_slope` is specified, then the combination of
            `ramp_slope`, `max_depth_of_cut` (or its default of "full
            depth"), and the length of the path may be overconstrained.
            If the specified ramp doesn't reach the target depth of cut
            before the end of the path, then the ramp slope is increased
            so the ramp reaches the target depth of cut exactly at the
            end of the path.

            A cutting pass consists of Entry Motion to get the tool into
            the work, followed by Main Motion cutting the path.

            Entry Motion:

                This gets the tool into the work.

                If `ramp_slope` has a value:

                    * Feed Z to the top of this pass (at `feed` feed
                      rate).

                    * Move X and Y along the path while ramping Z down to
                      the bottom of this pass at the specified slope rate.
                      The bottom of the pass is depth-of-cut below the
                      top of the pass.  The ramp may span segments, and
                      may terminate anywhere in any segment.  The segment
                      and location-within-the-segment where the ramp
                      finishes is noted for later cleanup.

                If `ramp_slope` is None:

                    * Feed Z to the bottom of this pass (at `plunge_feed`
                      feed rate).

                At the end of Entry Motion the tool is positioned
                somewhere on the path, with Z at the bottom of the pass.

            Main Motion:

                This cuts the path.

                Main Motion feeds along the path with Z fixed at the
                bottom of the pass, starting at whatever (X, Y) position
                Entry Motion left it and ending at the end of the path.

                If the Entry Motion used a ramp, that ramp is in front
                of the tool at the end of the pass.

        Final Motion:

            At the end of the last cutting pass, the tool is positioned
            at the (X, Y) starting point of the path, with Z at the
            `z_cut_depth` level.

            If the Entry Motion used ramp, then the final ramp remains
            in front of the tool and is cut away now, leaving (X, Y)
            at the end of the ramp move.

            If the `lead_out` argument is False, no Final Motion takes
            place and the tool is left where the Path Motion left it.

            If `lead_out` is True:

                * feed Z to `z_approach` level

                * rapid Z to `z_traverse` level

"""

    def cut_segment_skip_tabs(svg, seg_index, z_start, z_end, z_tab):
        """We're about to cut the specified segment, while moving the
        tool from `z_start` where it is now to `z_end` at the end of
        the segment, except it shouldn't cut into any tabs."""

        segment = path[seg_index]

        if segment_is_in_tab[seg_index]:
            # This segment is part of a tab.  If our Z
            # makes us interact with it, adjust things
            # to the tool goes over instead of through.
            if current_z < z_tab:
                # Feed up to clear the tab.
                g1(z=z_tab)
            if z_end < z_tab:
                # Adjust segment end Z from to clear tab.
                z_end = z_tab

        else:
            # Not in a tab.  If our Z is above where it should be
            # (because we just cleared a tab), plunge down then cut.
            if not close_enough(abs(current_z - z_start), 0.0):
                # Plunge to catch up to ramp.
                if plunge_feed is not None:
                    set_feed_rate(plunge_feed)
                g1(z=z_start)
                if plunge_feed is not None and feed is not None:
                    set_feed_rate(feed)

        # Inhibit Z coordinate if it's not needed.
        if close_enough(current_z, z_end):
            z_end = None
        path_segment_to_gcode(svg, segment, z=z_end)


    absolute_arc_centers()
    (x, y) = svg.xy_to_mm(path[0].start)

    if z_approach == None:
        z_approach = 0.5 + z_top_of_material

    # This is only used if work_holding_tabs > 0.
    z_tab = z_cut_depth + work_holding_tab_height

    path_length = path.length()
    if debug: print("path length: %f" % path_length, file=sys.stderr)
    length_between_tabs = path_length
    if work_holding_tabs > 0:
        # The user has requested N tabs.  There will be N sections of
        # the path between the tabs.
        length_between_tabs = (path_length - (work_holding_tabs * work_holding_tab_width)) / work_holding_tabs

    # Find the start and end of each work-holding tab.
    #
    # Split the segments there so that tabs start and end on a segment
    # boundary.
    #
    # Keep track of which segments are in tabs vs between tabs.

    # This array has an entry for each segment in `path`.  The entry is
    # True if the segment is entirely within a tab, False if the entry
    # is entirely outside all tabs.
    segment_is_in_tab = []

    for tab in range(work_holding_tabs):

        #
        # First deal with the segment where the tab *starts*.
        #

        start = tab * (work_holding_tab_width + length_between_tabs)
        if debug: print("default start of tab %d: %f" % (tab, start), file=sys.stderr)
        if len(work_holding_tab_locations) == work_holding_tabs:
            start = work_holding_tab_locations[tab]
            if debug: print("user-specified start of tab %d: %f" % (tab, start), file=sys.stderr)

        start_T = start / path_length
        if start_T > 1.0:
            raise ValueError("error: tab %d start at %f, past end of path at %f" % (tab, start, path_length))
        seg_index, t = path.T2t(start_T)

        if close_enough(t, 0.0):
            # The tab starts super close to the beginning of the segment,
            # keep the segment unsplit and count the whole thing as
            # "in the tab".

            # Everything from where the siit array ends, up to and
            # including the segment before this one is *before* this
            # tab starts.
            segment_is_in_tab += (seg_index - len(segment_is_in_tab)) * [False]

            # This segment segment is *in* this tab.
            segment_is_in_tab += [True]
            if debug: print("tab %d starts at %f (at start of seg %d, t %f)" % (tab, start, seg_index, t), file=sys.stderr)

        elif close_enough(t, 1.0):
            # The tab starts super close to the end of the segment,
            # keep the segment unsplit and count the whole thing as
            # "not in the tab".

            # Everything from where the siit array ends, up to and
            # including the segment before this one is *before* this
            # tab starts.
            segment_is_in_tab += (seg_index - len(segment_is_in_tab) + 1) * [False]

            if debug: print("tab %d starts at %f (at end of seg %d, t %f)" % (tab, start, seg_index, t), file=sys.stderr)

        else:
            # The tab starts somewhere in the middle of the segment,
            # split it.

            seg = path[seg_index]
            new_segments = seg.split(t)

            print("t=%f, 0.len=%f 1.len=%f" % (t, new_segments[0].length(), new_segments[1].length()), file=sys.stderr)

            path[seg_index] = new_segments[0]
            path.insert(seg_index+1, new_segments[1])

            # Everything from where the siit array ends, up to and including
            # the first part of the split segment is *before* this tab starts.
            segment_is_in_tab += (seg_index - len(segment_is_in_tab) + 1) * [False]

            if debug:
                print("tab %d starts at %f (in middle of seg %d, t %f)" % (tab, start, seg_index, t), file=sys.stderr)
                print("    original seg (%f): %s" % (seg.length(), seg), file=sys.stderr)
                print("    split seg 0 (%f): %s" % (new_segments[0].length(), new_segments[0]), file=sys.stderr)
                print("    split seg 1 (%f): %s" % (new_segments[1].length(), new_segments[1]), file=sys.stderr)

        #
        # Next deal with the segment where the tab *ends*.
        #

        end = start + work_holding_tab_width
        end_T = end / path_length
        if end_T > 1.0:
            raise ValueError("error: tab %d ends at %f, past end of path at %f" % (tab, end, path_length))
        seg_index, t = path.T2t(end_T)
        seg = path[seg_index]
        new_segments = seg.split(t)

        if close_enough(t, 0.0) or new_segments[0].length() <= 1e-2:
            # The tab end super close to the beginning of the segment,
            # keep the segment unsplit and count the whole thing as
            # "not in the tab".

            # Everything from where the siit array ends, up to and
            # including the segment before this one is *in* this tab.
            segment_is_in_tab += (seg_index - len(segment_is_in_tab)) * [True]

            # This segment segment is *not* in this tab.
            segment_is_in_tab += [False]

            if debug: print("tab %d ends at %f (at start of seg %d, t %f)" % (tab, end, seg_index, t), file=sys.stderr)

        elif close_enough(t, 1.0) or new_segments[1].length() <= 1e-2:
            # The tab ends super close to the end of the segment,
            # keep the segment unsplit and count the whole thing as
            # "in the tab".

            # Everything from where the siit array ends, up to and
            # including this segment is *in* this tab.
            segment_is_in_tab += (seg_index - len(segment_is_in_tab) + 1) * [True]

            if debug: print("tab %d ends at %f (at end of seg %d, t %f)" % (tab, end, seg_index, t), file=sys.stderr)

        else:
            # The tab ends somewhere in the middle of the segment,
            # split it.
            path[seg_index] = new_segments[0]
            path.insert(seg_index+1, new_segments[1])

            # Everything from where the siit array ends, up to and
            # including the first part of the split segment is *in*
            # this tab.
            segment_is_in_tab += (seg_index - len(segment_is_in_tab) + 1) * [True]

            if debug:
                print("tab %d ends at %f (in middle of seg %d, t %f)" % (tab, end, seg_index, t), file=sys.stderr)
                print("    original seg (%f): %s" % (seg.length(), seg), file=sys.stderr)
                print("    split seg 0 (%f): %s" % (new_segments[0].length(), new_segments[0]), file=sys.stderr)
                print("    split seg 1 (%f): %s" % (new_segments[1].length(), new_segments[1]), file=sys.stderr)

    # Everything from where the siit array ends, up to and including
    # where the path ends, is after the final tab.
    segment_is_in_tab += (len(path) - len(segment_is_in_tab)) * [False]


    #
    # Preliminary Motion
    #

    if lead_in == True:
        g0(z=z_traverse)
        g0(x=x, y=y)
        spindle_on()
        g0(z=z_approach)
    else:
        spindle_on()
        set_feed_rate(feed)
        g1(x=x, y=y)



    #
    # Cutting Passes
    #

    cut_depth = z_top_of_material - z_cut_depth

    if max_depth_of_cut == None:
        max_depth_of_cut = cut_depth

    num_passes_f = cut_depth / max_depth_of_cut
    num_passes = 0
    if close_enough(num_passes_f, round(num_passes_f)):
        num_passes = int(round(num_passes_f))
    else:
        num_passes = math.ceil(cut_depth / max_depth_of_cut)

    depth_of_cut = cut_depth / num_passes

    pass_num = 1
    while pass_num <= num_passes:
        z_top_of_pass = z_top_of_material - ((pass_num - 1) * depth_of_cut)
        z_bottom_of_pass = z_top_of_material - (pass_num * depth_of_cut)

        comment("path pass %d/%d, z from %f to %f" % (pass_num, num_passes, z_top_of_pass, z_bottom_of_pass))

        pass_num += 1


        #
        # Entry Motion
        #

        global current_z

        segment_after_ramp = 0

        if ramp_slope is not None:
            # Ramp in.

            # Start the ramp at the top of material
            if lead_in == True and feed is not None:
                set_feed_rate(feed)
            g1(z=z_top_of_pass)

            # If the path is so short that the ramp doesn't reach the
            # full requested depth of cut, then increase the ramp slope
            # so it does.
            # FIXME: maybe throw an exception here and force the user
            # to reconcile ramp_slope, depth_of_cut, and path length?
            full_path_ramp_slope = depth_of_cut / path.length()
            if full_path_ramp_slope > ramp_slope:
                print("WARNING: adjusting ramp slope from %.3f to %.3f to reach target depth of cut %.3f by the end of the path (length=%.3f)" % (ramp_slope, full_path_ramp_slope, depth_of_cut, path.length()), file=sys.stderr)
                ramp_slope = full_path_ramp_slope

            # During ramping, mix plunge_feed and feed according to
            # the slope.
            if feed is not None and plunge_feed is not None:
                ramp_feed = ((ramp_slope * plunge_feed) + feed) / (1.0 + ramp_slope)
                set_feed_rate(ramp_feed)
            elif plunge_feed is not None:
                set_feed_rate(plunge_feed)
            elif feed is not None:
                set_feed_rate(feed)

            # This variable shadows current_z, except that it does *not*
            # raise to skip over work-holding tabs.  It's the ramping
            # slope, ignoring tabs.
            z_ignoring_tabs = current_z

            # Ramp down, segment by segment, until we come to the bottom
            # of the pass or the end of the path.
            while (z_ignoring_tabs - 1e-6) > z_bottom_of_pass:
                i = 0
                while i < len(path):
                    segment = path[i]
                    length = segment.length()
                    z_at_segment_start = z_ignoring_tabs
                    z_at_segment_end = z_at_segment_start - (ramp_slope * length)
                    z_ignoring_tabs = z_at_segment_end

                    if (z_at_segment_end - 1e-6) > z_bottom_of_pass:
                        # Not bottoming out on this segment, it is in
                        # the middle of the ramp somewhere.
                        cut_segment_skip_tabs(svg, i, z_at_segment_start, z_at_segment_end, z_tab)
                        i += 1

                    else:
                        # Bottoming out somewhere on this segment!
                        # This finishes the ramp.  Cut to the end of the
                        # ramp, then break out of the ramping loop and let
                        # the normal "bottom of pass" loop finish the cut.
                        #
                        # Split the segment in two so the first part
                        # ends at the end of the ramp.  Take care to
                        # keep the segment_is_in_tab[] array up to date.
                        #
                        # In this segment, as t goes from 0.0 to 1.0:
                        # z(t) = z_at_segment_start - (ramp_slope * length * t)
                        #
                        # Setting z(t) equal to z_bottom_of_pass and
                        # solving for t gives us this equation:
                        t = (z_at_segment_start - z_bottom_of_pass) / (ramp_slope * length)
                        print("ramp bottom, t=%f" % (t), file=sys.stderr)

                        if close_enough(t, 0.0):
                            # Bottoming out right at the start of this segment, don't split.
                            segment_after_ramp = i;
                            if debug: print("ramp bottoms out at start of seg %d, t %f" % (i, t), file=sys.stderr)

                        elif close_enough(t, 1.0):
                            # Bottoming out right at the end of this segment, don't split.
                            cut_segment_skip_tabs(svg, i, z_at_segment_start, z_bottom_of_pass, z_tab)
                            segment_after_ramp = i+1;
                            if debug: print("ramp bottoms out at end of seg %d, t %f" % (i, t), file=sys.stderr)

                        else:
                            # Bottoming out in the middle of this segment, split.
                            new_segments = segment.split(t)
                            path[i] = new_segments[0]
                            path.insert(i+1, new_segments[1])
                            segment_is_in_tab.insert(i+1, segment_is_in_tab[i])
                            segment_after_ramp = i+1
                            cut_segment_skip_tabs(svg, i, z_at_segment_start, z_bottom_of_pass, z_tab)
                            if debug:
                                print("ramp bottoms out in middle of seg %d, t %f" % (i, t), file=sys.stderr)
                                print("    original seg (%f): %s" % (segment.length(), segment), file=sys.stderr)
                                print("    split seg 0 (%f): %s" % (new_segments[0].length(), new_segments[0]), file=sys.stderr)
                                print("    split seg 1 (%f): %s" % (new_segments[1].length(), new_segments[1]), file=sys.stderr)

                        # In any case we hit the bottom of the ramp,
                        # so we're done with Entry Motion and move on
                        # to Main Motion.
                        break

        else:
            # Plunge in.
            plunge_z = z_bottom_of_pass
            if segment_is_in_tab[0]:
                plunge_z = z_tab

            if current_z > plunge_z:
                if plunge_feed is not None:
                    set_feed_rate(plunge_feed)
                elif feed is not None:
                    set_feed_rate(feed)
                g1(z=plunge_z)


        #
        # Main Motion
        #

        if feed is not None:
            set_feed_rate(feed)

        for i in range(segment_after_ramp, len(path)):
            cut_segment_skip_tabs(svg, i, z_bottom_of_pass, z_bottom_of_pass, z_tab)

    # Done with all cutting passes.

    #
    # Final Motion
    #

    # Clean up the final ramp (if any).
    for i in range(0, segment_after_ramp):
        cut_segment_skip_tabs(svg, i, z_bottom_of_pass, z_bottom_of_pass, z_tab)

    if lead_out:
        g1(z=z_approach)
        g0(z=z_traverse)


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

current_feed = None


# Minimum arc radius accepted by LinuxCNC 2.7 (in mm).  An arc with
# a radius less than this will be rejected by the g-code interpreter,
# and the program will fail to run.  Any arcs smaller than this, replace
# with straight lines.
minimum_arc_radius = 0.00128

# When comparing floats, a difference of less than epsilon counts as no
# difference at all.
epsilon = 1e-5

def close_enough(a, b):
    """Returns True if the two numbers `a` and `b` are within `epsilon`
    (1e-6) of each other, False if they're farther apart."""
    return abs(a - b) < epsilon


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


def return_to_old_z():
    """After finishing a canned cycle, return the tool to the old Z
    coordinate (before starting the canned cycle)."""

    print("G98   ; after canned cycle, return to the old Z")


def return_to_r():
    """After finishing a canned cycle, return the tool to the R-word
    (retract level) of the canned cycle."""

    print("G99   ; after canned cycle, return to the R-word")


def spindle_on():
    print("M3")


def spindle_off():
    print("M5")


def path_blend(tolerance=None):
    print("G64 P%.4f (enable path blending with tolerance)" % tolerance)


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


def set_feed_rate(feed_rate_units_per_minute):
    print("F %.4f" % feed_rate_units_per_minute)


def speed(spindle_rpm):
    print("S %d" % spindle_rpm)


def coord_to_str(val=None):
    if val == None:
        return ""
    if close_enough(val, 0.0):
        # This avoids the "0.000"/"-0.000" confusion, mostly to make the tests consistent.
        val = 0.000
    return "%.4f" % val


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
            print(" X%s" % coord_to_str(x), end='')
        if y is not None:
            current_y = y
            print(" Y%s" % coord_to_str(y), end='')
        if z is not None:
            current_z = z
            print(" Z%s" % coord_to_str(z), end='')
        if a is not None:
            current_a = a
            print(" A%s" % coord_to_str(a), end='')
        if b is not None:
            current_b = b
            print(" B%s" % coord_to_str(b), end='')
        if c is not None:
            current_c = c
            print(" C%s" % coord_to_str(c), end='')
        if u is not None:
            current_u = u
            print(" U%s" % coord_to_str(u), end='')
        if v is not None:
            current_v = v
            print(" V%s" % coord_to_str(v), end='')
        if w is not None:
            current_w = w
            print(" W%s" % coord_to_str(w), end='')
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
            print(" X%s" % coord_to_str(x), end='')
        if y is not None:
            current_y = y
            print(" Y%s" % coord_to_str(y), end='')
        if z is not None:
            current_z = z
            print(" Z%s" % coord_to_str(z), end='')
        if a is not None:
            current_a = a
            print(" A%s" % coord_to_str(a), end='')
        if b is not None:
            current_b = b
            print(" B%s" % coord_to_str(b), end='')
        if c is not None:
            current_c = c
            print(" C%s" % coord_to_str(c), end='')
        if u is not None:
            current_u = u
            print(" U%s" % coord_to_str(u), end='')
        if v is not None:
            current_v = v
            print(" V%s" % coord_to_str(v), end='')
        if w is not None:
            current_w = w
            print(" W%s" % coord_to_str(w), end='')
        print()


def g2(x=None, y=None, z=None, i=None, j=None, p=None):
    global current_x
    global current_y
    global current_z

    """Clockwise arc feed."""
    if i is None and j is None:
        raise TypeError("gcoder.g2() without i or j")
    print("G2", end='')
    if x is not None:
        current_x = x
        print(" X%s" % coord_to_str(x), end='')
    if y is not None:
        current_y = y
        print(" Y%s" % coord_to_str(y), end='')
    if z is not None:
        current_z = z
        print(" Z%s" % coord_to_str(z), end='')
    if i is not None: print(" I%s" % coord_to_str(i), end='')
    if j is not None: print(" J%s" % coord_to_str(j), end='')
    if p is not None: print(" P%s" % coord_to_str(p), end='')
    print()


def g3(x=None, y=None, z=None, i=None, j=None, p=None):
    global current_x
    global current_y
    global current_z

    """Counter-clockwise arc feed."""
    if i is None and j is None:
        raise TypeError("gcoder.g3() without i or j")
    print("G3", end='')
    if x is not None:
        current_x = x
        print(" X%s" % coord_to_str(x), end='')
    if y is not None:
        current_y = y
        print(" Y%s" % coord_to_str(y), end='')
    if z is not None:
        current_z = z
        print(" Z%s" % coord_to_str(z), end='')
    if i is not None: print(" I%s" % coord_to_str(i), end='')
    if j is not None: print(" J%s" % coord_to_str(j), end='')
    if p is not None: print(" P%s" % coord_to_str(p), end='')
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
        print(" X%s" % coord_to_str(x), end='')
    if y is not None:
        current_y = y
        print(" Y%s" % coord_to_str(y), end='')
    if z is not None:
        print(" Z%s" % coord_to_str(z), end='')
    print(" R%s" % coord_to_str(retract), end='')
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
        print(" X%s" % coord_to_str(x), end='')
    if y is not None:
        current_y = y
        print(" Y%s" % coord_to_str(y), end='')
    if z is not None:
        print(" Z%s" % coord_to_str(z), end='')
    print(" R%s" % coord_to_str(retract), end='')
    print(" Q%s" % coord_to_str(delta), end='')
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

