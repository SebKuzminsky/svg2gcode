#!/bin/bash

do_test() {
    SVG_BASENAME="$1"
    OFFSET="$2"

    TEST_BASENAME=${SVG_BASENAME}.offset-${OFFSET}

    echo testing ${TEST_BASENAME}

    ../../svg2gcode --offset ${OFFSET} ${SVG_BASENAME}.svg >| ${TEST_BASENAME}.result.ngc 2>> ${TEST_BASENAME}.stderr
    diff -u ${TEST_BASENAME}.expected.ngc ${TEST_BASENAME}.result.ngc || exit 1

    rs274 -t /dev/null -v /dev/null -g ${TEST_BASENAME}.result.ngc >| ${TEST_BASENAME}.result.canon 2>> ${TEST_BASENAME}.stderr
    diff -u ${TEST_BASENAME}.expected.canon ${TEST_BASENAME}.result.canon || exit 1

    rm ${TEST_BASENAME}.result.ngc
    rm ${TEST_BASENAME}.result.canon
    rm ${TEST_BASENAME}.stderr
}


for SVG in *.svg; do
    SVG_BASENAME=$(basename ${SVG} .svg)
    for EXPECTED in $(ls -1 ${SVG_BASENAME}.offset-*.expected.ngc 2> /dev/null); do
        EXPECTED_BASENAME=$(basename ${EXPECTED} .expected.ngc)
        OFFSET=$(echo ${EXPECTED_BASENAME} | sed -re 's/^.*\.offset-(.*)$/\1/')
        do_test ${SVG_BASENAME} ${OFFSET}
    done
done
