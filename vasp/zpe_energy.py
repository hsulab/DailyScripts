#!/usr/bin/env python
from math import exp
from math import log
import sys
import os

if len(sys.argv) == 1:
    T = 500.0
elif len(sys.argv) == 2:
    T = float(sys.argv[1])
else:
    print('too many inputs')

detail=True
if detail:
    os.system('zp |tee zpe')
else:
    os.system('zp > zpe')
strlist=list(open('zpe'))
wlist=[ float(str.split('2PiTHz')[1].split()[0]) for str in strlist if 'f/i' not in str ]
h = 6.63*10**(-34)
c = 3*10**10
KB = 1.38*10**(-23)
Nmol = 6.023*10**23
eV2J = 96485.310432
R = 8.314
hcw= [ h*c*x for x in wlist ]
hcw_kbt= [ x/(KB*T) for x in hcw ]

zp=sum([ x*Nmol/2 for x in hcw ]) #in J/mol
S=sum([ R*(x/(exp(x)-1)-log(1-exp(-x))) for x in hcw_kbt ]) #in J/mol/K
U=sum([ R*x/KB/(exp(x/KB/T)-1) for x in hcw ]) #in J/mol
TS=T*S

zp_ineV=zp/eV2J
TS_ineV=TS/eV2J
U_ineV=U/eV2J

total=zp_ineV+U_ineV-TS_ineV
if detail:
    print('Temperature = %s K\tZPE = %3.2f eV\tU = %3.2f eV\tTS = %3.2f eV\tZPE+U-TS = %3.2f eV\n'%(T,zp_ineV,U_ineV,TS_ineV,total))
    print(wlist)
else:
    print('%3.2f   %3.2f   %3.2f   %3.2f'%(zp_ineV,U_ineV,TS_ineV,total))
    print(wlist)
