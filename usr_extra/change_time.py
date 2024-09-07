
#import io
#import os
import struct
import sys
import re

if len(sys.argv) !=2 :
    print('\nChange the time inside *0.f0* file based on timestep')
    print('Usage: python3 ./%s avg<case_name>.nek5000 \n\n'%(sys.argv[0]))
    quit(0)

fnek5000 = sys.argv[1]


try:
    print('Reading '+fnek5000,end='')
    f = open(fnek5000)
except FileNotFoundError:
    print('%s file not found!'%fnek5000)
    quit(1)
else:
    with f:
        inlines = f.readlines()
        
        cname = inlines[0].split()[1].split('%')[0]
        regx = inlines[0].split()[1].split('%')[1]
        firsttimestep = int(inlines[1].split()[1])
        numtimesteps = int(inlines[2].split()[1])
        num_zero = int(regx[1])

        print('   [ %s %d %d %d ]'%(cname,num_zero,firsttimestep,numtimesteps))


for i in range(numtimesteps):
    ifile = i + firsttimestep
    fname = cname+'0'*num_zero+'.f%05d'%ifile
    print('  %s'%fname,end ="")

    with open(fname,'r+b') as f:
        inline = f.read(132).decode("utf-8")
    
        header = inline.split()
        time = float(header[7])
        timestep = int(header[8])
        newTime = timestep * 1.0
    
        newHeader = inline[:38] + '%20.13E'%newTime + inline[58:]
        print('   step=%d time=%g newTime=%g'%(timestep,time,newTime))
#       print(inline)
#       print(newHeader)
        f.seek(0)
        f.write(newHeader.encode("utf-8"))

