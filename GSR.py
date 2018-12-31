# -*- coding: utf-8 -*-
"""
Created on Fri May 13 09:56:44 2016

@author: Admin
"""
from math import *
for i in range(1,15):
    eps = 10**(-i)
    num1 = log(eps**(-2))/log(1+eps)
    num2 = num1/(eps**(-3))
    num3 = num1/(eps**(-2))
    num4 = num1/(eps**(-1))
    num5 = num4/log(1/eps)
    num6 = (num1/(log(1/eps)))*eps
    print(i,eps,num1,num2,num3,num4,num5,num6)
    
    
    
    