# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 19:11:00 2016

@author: Evyatar
"""

import time
from math import *
import numpy as np
from scipy import optimize
from Instance import *

 
def const1(x):
    return 1.0
def slope1(x):
    return x
   
csr = Instance(slope1,const1,1.0)
detStr = csr.optDetStr()
if detStr[0]==float("Inf"):
     print "The optimal deterministic strategy never switch and has a competitive ratio of",detStr[1]
elif detStr[1]==0.0:
    print "The optimal deterministic strategy start at buy state and has a competitive ratio of",detStr[1]
else:
    print "The optimal deterministic strategy switch at time",detStr[0],"and has a competitive ratio of",detStr[1]




