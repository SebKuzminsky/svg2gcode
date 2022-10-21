 #!/usr/bin/env python
"""\
Simple g-code streaming script for grbl

Provided as an illustration of the basic communication interface
for grbl. When grbl has finished parsing the g-code block, it will
return an 'ok' or 'error' response. When the planner buffer is full,
grbl will not send a response until the planner buffer clears space.

G02/03 arcs are special exceptions, where they inject short line 
segments directly into the planner. So there may not be a response 
from grbl for the duration of the arc.

---------------------
The MIT License (MIT)

Copyright (c) 2012 Sungeun K. Jeon

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
---------------------
"""

import serial
import time
import os
import fnmatch

# Open grbl serial port
grblSerialPort = serial.Serial('/dev/ttyUSB0',115200)
log = open('log/failed_tests_grbl.log','w')
log.close()
for root, dirnames, filenames in os.walk('svg2gcode'):
    for filename in fnmatch.filter(filenames, 'expected.ngc'):
        gcodeFileName=os.path.join(root, filename)
        print(gcodeFileName)
        # Open g-code file
        gcode = open(gcodeFileName,'r')
        log = open('log/failed_tests_grbl.log','a')

        grblSerialPort.write(("\r\n\r\n").encode("utf8")) # Wake up grbl
        time.sleep(2)   # Wait for grbl to initialize 
        grblSerialPort.flushInput()  # Flush startup text in serial input

        # Stream g-code to grbl
        for line in gcode:
            lineStrip = line.strip() # Strip all EOL characters for consistency
            if lineStrip != 'M2':
                grblSerialPort.write((lineStrip + '\n').encode("utf8")) # Send g-code block to grbl
                if b'error' in (grblSerialPort.readline()).strip():
                    print('  Test Error: ' + lineStrip)
                    log.write(gcodeFileName + '\n')
                    break
                print('  Test OK: ' + lineStrip)
        # Close file and serial port
        gcode.close()

grblSerialPort.close()
log.close()
