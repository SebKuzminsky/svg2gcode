all: svg2gcode.1

%.1: %.adoc
	a2x -v -d manpage -f manpage $<
