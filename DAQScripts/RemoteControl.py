# Script to remotely control a LeCroy oscilloscope
# For each trigger two waveforms per channel are acquired:
#   wfa = detailed waveform (typ. 2000 points) centered at trigger time
#   wfb = full waveform (up to 1000 points)

import re
import numpy
import pyvisa
from pyvisa.resources import MessageBasedResource
import sys
import time

# Sampling interval in ns (normally: 25 ps @ 40 GS/s)
ns_sample = 0.025

# Number of points for "wfa"
points_a = '2000'
npa = int(points_a)

# Maximum number of points for "wfb"
points_b = '1000'

# Sparsification parameters for waveforms "a" and "b"
# assuming 25 ps/point, waveform "b" will be sampled at 10 ns/point
sparse_a = '1'
sparse_b = '400'

## Ask the number of triggers to be collected
ntrig = int(input('Enter the number of triggers: '))
print("Requested triggers: "+str(ntrig)+"\n")

# Ask the number of channels to be acquired
# TBD - default is all 4 channels

## Ask the trigger time position from the beginning of the waveform
## and compute the first datapoint for the finely sampled 'wfa'
ns_delay = int(input('Enter the trigger time position (ns): '))
first = int(ns_delay / ns_sample) - (npa/2)
first_s = str(first)
print("First data point (wfa) = "+first_s+"\n")

# Ask for choice between ASCII and binary data transfer
# AT THE MOMENT ONLY ascii IS POSSIBLE, binary is not tested
data_transfer = ''
#while data_transfer != 'ascii' and data_transfer != 'binary':
#    data_transfer = input('Choose between ascii and binary data transfer: ')
data_transfer = 'ascii'
    
## Last question: ask the filename for saving waveform data
filnam = input(('Enter the filename (with extension) for waveform saving: '))
if bool(filnam):
    print("Filename for waveforms: "+filnam+"\n")
    fout = open(filnam,'w')

orig_stdout = sys.stdout

rm = pyvisa.ResourceManager()

# Connect to the 'scope over Ethernet
scope = rm.open_resource('TCPIP0::193.205.66.32::inst0::INSTR',
                         resource_pyclass=MessageBasedResource)

# Set timeout in ms (1500 s for radioactive source op.) and clear buffers
scope.timeout = 1500000
scope.clear()

# Remove comm header from scope answers
scope.write("CHDR OFF")

print(scope.query("*IDN?"))

# Read Vdiv for Ch 1
scope.write("C1:VDIV?")
value = scope.read()
value = "Ch.1 vertical scale = " + value   
print (value)

# (Set) / Read Tdiv
#scope.write("TDIV 5E-4")
scope.write("TDIV?")
value = scope.read()
value = "Horizontal scale = " + value   
print (value)

# Read memory size
scope.write("MSIZ?")
value = scope.read()
value = "Memory size = " + value   
print (value)

# Read trigger selection
scope.write("TRIG_SELECT?")
value = scope.read()
print(value)

# Read trigger level
scope.write("TRLV?")
value = scope.read()
value = "Trigger level = " + value   
print (value)

# Read trigger slope
scope.write("TRIG_SLOPE?")
value = scope.read()
value = "Trigger slope = " + value   
print (value)

# Read trigger mode
scope.write("TRMD?")
value = scope.read()
value = "Trigger mode = " + value   
print (value)

# Set (and read) trigger mode to STOP - legacy
scope.write("TRMD STOP")
scope.write("TRMD?")
value = scope.read()
value = "New Trigger mode = " + value   
print (value)

### Loop on triggers

for i in range(ntrig):
    
    # Set trigger mode to SINGLE - VBS syntax or legacy syntax
    #command = "VBS 'app.Acquisition.triggermode = \"single\"' "
    #scope.write(command)
    scope.write("TRMD SINGLE")
    scope.write("TRMD?")
    value = scope.read()
    value = "New trigger mode = " + value   
    print (value)

    # Wait for next trigger
    scope.write("WAIT")
    scope.write("*OPC?")
    value = scope.read()
    
    print(" scope has captured a waveform\n")

    print('{0:10s} {1:5d}'.format("Trigger n.",i))
    if bool(filnam):
        sys.stdout = fout
        print('{0:10s} {1:5d}'.format("Trigger n.",i))
        sys.stdout = orig_stdout
    
    # Get waveform info (HORIZ_INTERVAL IS THE SAME FOR ALL CH's)
    #print(scope.query("C1:INSPECT? HORIZ_INTERVAL"))
    scope.write("C4:INSPECT? HORIZ_INTERVAL")
    line = scope.read()
    print('C4: '+line)
    m = re.search(r'\d+\.\d+([eE][+-]\d+)',line)
    horiz = float(m.group(0))
    
    #print(scope.query("C1:INSPECT? VERTICAL_GAIN"))
    print(scope.query("C4:INSPECT? VERTICAL_GAIN"))

    #print(scope.query("C1:INSPECT? TRIGGER_TIME"))
    scope.write("C4:INSPECT? TRIGGER_TIME")
    value = scope.read()
    print(value)
    if bool(filnam):
        sys.stdout = fout
        print(value)
        sys.stdout = orig_stdout
    
    #print(scope.query("C1:INSPECT? WAVE_ARRAY_1"))
    #print(scope.query("C1:INSPECT? WAVE_ARRAY_2"))
    print(scope.query("C4:INSPECT? WAVE_ARRAY_1"))
    #print(scope.query("C4:INSPECT? WAVE_ARRAY_2"))

    # Get full description for Ch 1 waveform (ASCII format)
    #print(scope.query("C1:INSPECT? WAVEDESC"))

    # Limit the number of data points acquired from the waveform:
    #   SP = sparsing parameter; NP = max. n. of data points;
    #   FP = first data point (counting from 0);
    #   SN = segment number

    # Loop on the channels

    for j in range(4):

        j1 = j+1
        chan = 'C'+str(j1)

        print("reading channel "+chan)
    
        scope.write(chan+":WFSU SP,"+sparse_a+",NP,"+points_a+",FP,"+first_s+",SN,1")
        
        if data_transfer == 'ascii':
            # Get selected data points -> waveform "a" (ASCII format)
            scope.write(chan+":INSPECT? DATA_ARRAY_1")
            wfdata = scope.read()
            if bool(filnam):
                sys.stdout = fout
                print('{0:20s}{1:s}{2:18s}{3:e}{4:10s}{5:2s}{6:s}'.format(
                "Waveform with SP of ",sparse_a,", Horiz_interval of ",horiz,
                ", Channel ",chan,":"))
                print(wfdata)
                sys.stdout = orig_stdout
            else:
                print(wfdata)
        elif data_transfer == 'binary':
            scope.write(chan+":WF? DAT1")
            dataword = scope.read()
            # skip first 16 bytes when CHDR SHORT
            wfbindata = numpy.fromstring(dataword, dtype = numpy.int16)
            print(wfbindata)
            wfldata = numpy.array(wfbindata, dtype = numpy.float)
            print(wfldata)

        if data_transfer == 'ascii':
            scope.write(chan+":WFSU SP,"+sparse_b+",NP,"+points_b+",FP,0,SN,1")
            # Get selected data points -> waveform "b" (ASCII format)
            scope.write(chan+":INSPECT? DATA_ARRAY_1")
            wfdata = scope.read()    
            if bool(filnam):
                sys.stdout = fout
                print('{0:20s}{1:s}{2:18s}{3:e}{4:10s}{5:2s}{6:s}'.format(
                "Waveform with SP of ",sparse_b,", Horiz_interval of ",horiz,
                ", Channel ",chan,":"))
                print(wfdata)
                sys.stdout = orig_stdout
            else:
                print(wfdata)

if bool(filnam):
    fout.close()

# Disconnect the oscilloscope
scope.close()
rm.close()
