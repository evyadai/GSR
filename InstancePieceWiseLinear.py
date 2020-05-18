# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 19:11:00 2016

@author: Evyatar
"""

from math import *
import numpy as np

np.random.seed(1)

import scipy
from Instance import *
from Linear import *
from PieceWiseLinear import *

import matplotlib
import matplotlib.mlab as mlab
# insert "%matplotlib tk"
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle
from matplotlib import cm


class InstancePieceWiseLinear(Instance):

    def __init__(self, fr, fb, c):
        self.fr = fr
        self.fb = fb
        self.c = np.float64(c)
        self.times = fr.times
        self.norm_cost, self.norm_time = 1.0, 1.0
        if self.fr[1] == self.fb[1] == 1:
            self.norm = True
        else:
            self.norm = False

    def __str__(self):
        s = ""
        for t, lr, lb in zip(self.times, self.fr.lins, self.fb.lins):
            s = s + "Time>= " + str(t) + ":  Rent: " + str(lr) + " Buy: " + str(lb) + "\n"
        return s

    def intersection(self):
        eps_num = 1e-10
        for i, t in enumerate(self.times):
            int_t = div(self.fb.lins[i].offset - self.fr.lins[i].offset, self.fr.lins[i].slope - self.fb.lins[i].slope)
            if i < (len(self.times) - 1):
                if int_t >= t and int_t <= self.times[i + 1]:
                    return int_t

            elif int_t + eps_num >= t:
                return int_t

    def normalize(self):
        # we need to init new instance and return this
        intersect = self.intersection()
        if intersect is None:
            return
        val_inter = self.fr[intersect]

        # print intersect,val_inter
        # all slopes need to be updated (*intersections/val_intersect)
        # all times need to be updated (/intersection)
        # all offsets need to be updated (/val_inter)

        normTimes = [t / intersect for t in self.times]
        diagLinn = Linear(1.0, 0.0)
        normFr = [Linear(l.slope * intersect / val_inter, l.offset / val_inter) for l in self.fr.lins]
        normFb = [Linear(l.slope * intersect / val_inter, l.offset / val_inter) for l in self.fb.lins]
        for ind in range(len(normTimes)):
            next_t = normTimes[ind + 1] if ind < len(normTimes) - 1 else np.inf
            if 1. > normTimes[ind] and 1. <= next_t:
                ind_intersect = ind + 1

        if 1.0 not in self.times:
            normTimes.insert(ind_intersect, 1.0)
            normFr.insert(ind_intersect, normFr[ind_intersect - 1])
            normFb.insert(ind_intersect, normFb[ind_intersect - 1])

        normed = InstancePieceWiseLinear(PieceWiseLinear(normFr, normTimes), PieceWiseLinear(normFb, normTimes),
                                         self.c / val_inter)
        normed.norm = True
        normed.norm_time = intersect
        normed.norm_cost = val_inter
        return normed

    def plot(self, num_fig=1):
        t_max = self.times[-1] * 2.0
        tt = self.times[:]
        tt.append(t_max)
        cost_r = [self.fr[t] for t in tt]
        cost_b = [self.fb[t] for t in tt]
        xx = [t for t in tt]

        plt.figure(num_fig)
        plt.clf()
        plt.hold(True)
        plt.plot(xx, cost_r, 'r',linewidth=2)
        plt.plot(xx, cost_b, 'b',linewidth=2)
        plt.axis([0.0, xx[-1] * 1.0, 0.0, 1.1 * max(np.max(cost_r), np.max(cost_b))])
        fs = 36
        plt.xlabel("Time", fontsize=fs)
        plt.ylabel("Cost", fontsize=fs)
        # print cost_b,xx
        return plt

    def plot_q(self, eps, norm=True, num_fig=1):
        cr, q, T = self.approxOptStr(eps)
        Q = np.cumsum(q)
        if norm:
            xx = [t * self.norm_time for t in T]
            Q = Q * self.norm_cost
        else:
            xx = T[:]
        plt.figure(num_fig)
        plt.hold(True)
        plt.title("Instance")
        plt.plot(xx, Q, 'g')
        plt.legend(["Rent", "buy", "Probabilty"], loc=0)
        plt.xlabel("Time")
        plt.ylabel("Cost (Probabilty)")
        plt.text(xx[-2], self.norm_cost, "Probabilty 1.0")
        plt.text(xx[-2], 0.5 * self.norm_cost, "Probabilty 0.5")
        plt.show()

        print "plot_q", xx, Q

    def cr_z(self, z):
        if z == 0:
            return self.fb.lins[0].offset / self.fr.lins[0].offset
        elif z == np.inf:
            if self.fr.lins[-1].slope == self.fb.lins[-1].slope:
                return self.fr[np.inf] / self.fb[np.inf]
            return self.fr.lins[-1].slope / self.fb.lins[-1].slope
        else:
            raise Exception("Error in InstnacePieceWiseLinear.cr_z")

    def cr_0(self):
        return div(self.fb.lins[0].offset, self.fr.lins[0].offset)

    def cr_inf(self):
        return div(self.fr.lins[-1].slope, self.fb.lins[-1].slope) if self.fr.lins[-1].slope > 0 else self.fr(
            self.times[-1] + 10) / self.fb(self.times[-1] + 10)

    def cr_inf_overall(self):
        ts = self.times
        max_ts = np.max([self.fr(t)/self.fb(t) for t in ts])
        max_inf = self.cr_inf()
        return np.max([max_ts,max_inf])



    def zm(self):
        pass

    def cr_zm(self):
        return h(self.zm())

    # todo - g in linear g(z)=max fr/fb(y) (y>z)
    def g(self, z):
        pass

    # todo -h(z)= (fr(z)+c)/fb(z)
    def h(self, z):
        pass

    # we iterate on all the lines and find if there is inersection between h and g
    # if there is we find the exact intersection
    # otherwise we return inf
    def intersection_time(self):
        # print "inter"
        for ind, t in enumerate(self.times):
            if t >= 1:
                # at point 1 g<h and then g always up
                # we need to divide to 2 parts:
                # 1.h decrese or not up - then we need only to look on the start and end
                # 2.h increase - then we need to solve
                # TODO - try to checjk sopes or use ratio
                g1 = self.fr[t] / self.fb[t]
                h1 = (self.fr[t] + self.c) / self.fb[t]
                if ind == len(self.times) - 1:
                    next_t = np.inf
                    g2_n, g2_d = self.fr[next_t], self.fb[next_t]
                    if g2_d == np.inf:
                        g2 = div(self.fr.lins[ind].slope, self.fb.lins[ind].slope)
                    else:
                        g2 = div(g2_n, g2_d)

                    h2_n, h2_d = self.fr[next_t] + self.c, self.fb[next_t]
                    if h2_d == np.inf:
                        h2 = div(self.fr.lins[ind].slope, self.fb.lins[ind].slope)
                    else:
                        h2 = div(h2_n, h2_d)
                else:
                    next_t = self.times[ind + 1]
                    g2 = self.fr[next_t] / self.fb[next_t]
                    h2 = (self.fr[next_t] + self.c) / self.fb[next_t]
                # print "inter",ind,g1,g2,h1,h2
                # TODO - find intersection better (maybe the intersection function will do this)
                # if max(g1,g2) >=min(h1,h2):
                if True:
                    # we have here intersection - what is the solve?
                    # if g1<=g2 then the function is monotone - simple solve
                    if g1 < g2:
                        # we need to find zs:
                        roots_intersection = ratio_intersection(self.fr.lins[ind], self.fb.lins[ind],
                                                                self.fr.lins[ind] + self.c, self.fb.lins[ind])
                    if g1 >= g2:
                        # then the function is const
                        roots_intersection = ratio_intersection(Linear(0, g1), Linear(0, 1.0),
                                                                self.fr.lins[ind] + self.c, self.fb.lins[ind])
                    inds = np.where((roots_intersection >= t) & (roots_intersection <= next_t))[0]
                    if inds.size > 0:
                        return roots_intersection[inds[0]]
                else:
                    pass

        return float("inf")

    def middle_time_cr(self):
        z_s = self.intersection_time()
        # print "z_s at z_m",z_s
        # the competitve is always h(z) - we need to find the maximum
        cr_m = (self.fr[1] + self.c) / self.fb[1]
        z_m = 1.0
        # because at each interval h always goes
        for ind, t in enumerate(self.times):
            if t >= 1 and t <= z_s:
                ht = (self.fr(t) + self.c) / self.fb(t)
                if ht < cr_m:
                    z_m = t
                    cr_m = ht
                # print "iter at z_m",ind,t,ht,cr_m

        # at the end we need to check z_s or np.inf
        if t < z_s:
            t = z_s
            ht = self.cr_inf() if z_s == np.inf else (self.fr(t) + self.c) / self.fb(t)
            if ht < cr_m:
                z_m = t
                cr_m = ht
            # print "last iter at z_m",ht,cr_m,self.cr_inf() if z_s==np.inf else (self.fr(t)+self.c)/self.fb(t),z_s==np.inf,self.fr(t),self.c,self.fb(t),t
        # print "middle_time_cr",(z_m,cr_m)
        return (z_m, cr_m)

        """
                if ind==len(self.times)-1:
                    next_t=min(np.inf,z_s)
                    h2_n,h2_d = self.fr[next_t]+self.c,self.fb[next_t]
                    if h2_d==np.inf:
                        h2=div(self.fr.lins[ind].slope,self.fb.lins[ind].slope)
                    else:
                        h2=div(h2_n,h2_d)
                else:
                    next_t=min(self.times[ind+1],z_s)
                    h2 = (self.fr[next_t]+self.c)/self.fb[next_t]
                if cr_m > min(h1,h2):
                    cr_m = min(h1,h2)
                    if h1<=h2:
                        z_m=self.times[ind]
                    else:
                        z_m=next_t
        return (z_m,cr_m)
        
        """

    def middle_time(self):
        return self.middle_time_cr()[0]

    # todo
    def between(self, t, ind):
        if ind == len(self.times) - 1:
            if t >= self.times[ind]:
                return True
        elif (t >= self.times[ind]) and (t <= self.times[ind + 1]):
            return True
        return False

    def calc_t1(self, eps):
        # t1
        # f_{r}(z)+c =(1+eps)(f_{rb}(0) +c+f_b(z)
        # frb(t1)-eps*f_b(t1)=(1+eps)f_{rb}(0)+eps*c

        # t1 = min(np.float64(1.0),)

        for ind in range(len(self.times)):
            if self.times[ind] >= 1.0:
                return 1.0
            else:
                lin = self.fr.lins[ind] - self.fb.lins[ind] * (1.0 + eps)
                t1_candidate = lin.inverse((1 + eps) * self.f_rb(0.0) + eps * self.c)
                if (t1_candidate >= self.times[ind]) and (t1_candidate <= self.times[ind + 1]):
                    return min(np.float64(1.0), t1_candidate)

    # todo
    # fix situation when the ratio raise
    def calc_t_N_minus_1(self, eps):
        global lin
        # t_{N-1}\gets \min\{f_{r}^{-1}(c/ \epsilon) ,\min\{x:x\geq 1,\forall y\geq x,\frac{f_r(x)}{f_b(x)}\geq (1+\epsilon)\frac{f_r(y)}{f_b(y)} \}\} $

        if (self.c <= self.fr[1.0] * eps):
            return 1.0
        val1 = max(1.0, self.fr.inverse(self.c / eps))
        maxRatio = self.cr_inf()
        # print self.fr.lins[-1],val1,self.c/eps
        # print maxRatio/(1+eps),self.fr.inverse(self.c/eps)
        # print "val1",val1

        if maxRatio < np.inf:
            maxRatioLin = [maxLin(self.fr.lins[ind - 1], self.fb.lins[ind - 1], self.times[ind - 1], self.times[ind])[0]
                           for ind in range(1, len(self.times))]
            maxRatioLin.append(maxLin(self.fr.lins[-1], self.fb.lins[-1], self.times[-1], np.inf)[0])
            # print maxRatioLin,maxLin(self.fr.lins[-1],self.fb.lins[-1],self.times[-1],np.inf)
            # fr = maxRatio/(1+eps) *fb
            # print "max Ratios",maxRatioLin
            for ind in range(len(self.times)):
                if self.times[ind] < 1.0:
                    continue
                # TODO - if ratio is down we need to consider from (ind+1) else ind
                maxRatio = max(maxRatioLin[ind:])
                # print ind,maxRatio,self.times[ind],eps

                lin_cond = self.fr.lins[ind] - self.fb.lins[ind] * (maxRatio / (1 + eps))
                val2_candidate = lin_cond.inverse(0.0)  # this is the point where the ratio is fine
                if (val2_candidate < self.times[ind]) and (
                        self.fr[self.times[ind]] / self.fb[self.times[ind]] >= maxRatio / (1 + eps)):
                    # this case when the candidate is out of range (and we take the start of the range)
                    # print "oor",val2_candidate,self.times[ind],maxRatio,self.fr.lins[ind],self.fb.lins[ind]
                    val2_candidate = self.times[ind]
                    # print self.between(val2_candidate,ind)
                if self.between(val2_candidate,
                                ind):  # (val2_candidate>=self.times[ind]) and (val2_candidate<= self.times[ind+1]) :
                    val2_cand_ratio = self.fr[val2_candidate] / self.fb[val2_candidate]
                    if (val2_cand_ratio + 1e-5 >= (maxRatio / (1 + eps))):
                        # print "return",ind,min(val1,val2_candidate),val2_cand_ratio,maxRatio/(1+eps),val1,val2_candidate
                        # print self.fr[50]/self.fb[50]
                        return min(val1, val2_candidate)

        return val1

    def valid(self):
        if min([self.times[1:][j] - self.times[:-1][j] for j in range(len(self.times) - 1)]) < 0.:
            return False
        if min([min(self.fr.lins[i].slope, self.fb.lins[i].slope) for i in range(len(self.times))]) < 0:
            return False
        return True

