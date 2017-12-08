#!/usr/bin/python

from scipy import fromfile
from pylab import *
from sys import argv

#c=loadtxt('exp-info',dtype='float',delimiter=' ')
#print c
NI,NJ,NK = float(argv[1]),float(argv[2]),float(argv[3])
fn = argv[4]

a=fromfile(fn)
b=a.reshape(-1,NK,NJ,NI)
figure(1)
contourf(b[-1,-1,:,:])
colorbar()
title('xy end')

figure(2)
contourf(b[-1,:,:,2])
title('yz end')
figure(3)
contourf(b[0,:,:,2])
title('yz start')
show()
