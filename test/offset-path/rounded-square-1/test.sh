#!/bin/bash

../../../svg2gcode --offset 1 ../rounded-square-equal-radii.svg >| result.ngc
diff -u expected.ngc result.ngc || exit 1

rs274 -t /dev/null -v /dev/null -g result.ngc >| result.canon
diff -u expected.canon result.canon || exit 1
