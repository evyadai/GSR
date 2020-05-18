from math import *
import numpy as np

np.random.seed(1)

import scipy
from Instance import *
from Linear import *
from PieceWiseLinear import *
from InstancePieceWiseLinear import *

import matplotlib
import matplotlib.mlab as mlab
# insert "%matplotlib tk"
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle
from matplotlib import cm


# check
"""
pr=PieceWiseLinear([Linear(1.0/3,2.0/3.0),Linear(1.0/3,2.0/3.0),Linear(4.0/3,0)],[0.0,1.0,1.5])
pb=PieceWiseLinear([Linear(1.0,0.0),Linear(1.0,0.0),Linear(1.0,0.0)],[0.0,1.0,4.0/3])
exmp2=InstancePieceWiseLinear(pr,pb,1.0)



r0,b0=35.0,45.0
t1,t2=6.0,10.0
ar1,ab1=6.0,6.0
ending=True
if ending:
    endT=6.0
    pr=PieceWiseLinear([Linear(r0,0.0),Linear(r0-t1*ar1,ar1),Linear(r0-t1*ar1,ar1),Linear(r0+endT*ar1,0.0)],[0.0,t1,t2,t1+endT])
    pb=PieceWiseLinear([Linear(b0,0.0),Linear(b0,0.0),Linear(b0-t2*ab1,ab1),Linear(b0-t2*ab1+ab1*(t1+endT),0.0)],[0.0,t1,t2,t1+endT])
else:
    pr=PieceWiseLinear([Linear(r0,0.0),Linear(r0-t1*ar1,ar1),Linear(r0-t1*ar1,ar1)],[0.0,t1,t2])
    pb=PieceWiseLinear([Linear(b0,0.0),Linear(b0,0.0),Linear(b0-t2*ab1,ab1)],[0.0,t1,t2])
cons=InstancePieceWiseLinear(pr,pb,(b0-r0)*1.35)
cons.plot(False,2)
cons.normalize()
print cons.calc_t_N_minus_1(0.01)
cons.plot(False,3)
#plt.figure(3)
#plt.plot(np.arange(1.0,2.5,0.01),[cons.fr[i]/cons.fb[i] for i in np.arange(1.0,2.5,0.01)],'g-')
cons.plot_q(0.01,True,2)



%[b,a1,break,a2]
%HOT - 8G=40sh +1G=6sh
%12G - 50sh + 1G=6sh
%first we deal only what happen after the first 8G
%50s->1s , 1.5G-1t
% fr=[0.8,6*1.5/50,4/1.5,6*1.5/50];
% fb=[1,0,4/1.5,6*1.5/50];
% c=1.0;
%result : Never switch


% fr=[40,6,4,6];
% fb=[50,0,4,6];
% c=50;

%4.2 - we have no tight
% fr=[0.25,0.75,5/3,0.2];
% fb=[1,0,5/3,0.1];
% c=1;

%4.3 - we have unoptimal tight
%Intersection (t=8/9,cost=1)
fr=[1/3,0.75,4/3,0];
fb=[1,0,4/3,0];
c=1;
"""


def UT_linear():
    global i, xx, gg, hh
    print "Unit test linear ski rental problems"
    eps = 0.1
    eps_num = 1e-10
    rr = np.arange(0.0, 1.02, 0.1)
    num_iter = 100
    for eps in [0.1]:
        print "epsilon", eps
        for i in range(num_iter):
            a1 = max(np.random.choice(rr), 0.05)
            a2 = 2
            while a2 >= a1:
                a2 = np.random.choice(rr)
            # a1,a2=1.0,0.2
            b1, b2 = 1 - a1, 1 - a2
            c = (b2 - b1) * np.random.uniform(1, 5)
            # a1,b1,a2,b2,c=0.5,0.5,0.0,1.0,0.684677189538

            # print ",".join(str(u) for u in [a1,b1,a2,b2,c])
            pr = PieceWiseLinear([Linear(a1, b1), Linear(a1, b1)], [0, 1])
            pb = PieceWiseLinear([Linear(a2, b2), Linear(a2, b2)], [0, 1])
            cons = InstancePieceWiseLinear(pr, pb, c)
            # test 1

            xx = np.arange(1.0, 20.0, 0.01)
            gg = np.array([cons.fr[x] / cons.fb[x] for x in xx])

            gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
            hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
            gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
            hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
            # plt = cons.plot()
            # plt.hold(True)
            # plt.plot(xx,gg,"g")
            # plt.plot(xx,hh,"y")
            assert (np.abs(cons.intersection() - 1.0) < eps_num)
            # print cons.calc_t1(eps)

            t1 = cons.calc_t1(eps)
            # print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
            assert (-eps_num < t1 < 1.0 + eps_num)
            assert ((cons.fr[t1] + c) <= eps_num + (1 + eps) * (cons.f_rb(0) + c + cons.fb[t1]))
            assert ((np.abs(t1 - 1.0) < eps_num) or np.abs(
                cons.fr[t1] + c - (1 + eps) * (cons.f_rb(0) + c + cons.fb[t1])) < eps_num)

            alp = cons.calc_t_N_minus_1(eps)
            # print alp
            assert (alp >= 1.0)
            # print alp,cons.rho_inf(),cons.fr(alp)
            flag1 = (cons.fr[alp] / cons.fb[alp]) <= eps_num + (1 + eps) * cons.cr_inf()
            # print alp,c,cons.fr[alp]
            flag2 = (c <= eps_num + eps * cons.fr[alp])
            assert (flag1 or flag2)

            beta = (b2 - b1) / c
            rho_prob_anylsis = np.exp(beta) / (np.exp(beta) - b2 + b1)

            # print cons.approxOptStr(eps)[0],rho_prob_anylsis
            ratio = cons.approxOptStr(eps)[0] / rho_prob_anylsis
            assert ((ratio <= eps_num + (1 + eps)) and (ratio >= -eps_num + 1 / (1 + eps)))

            # check det
            rho_det_anylsis = min(div(b2, b1), c + 1, div(a1, a2))
            t_det_anylsis = [0, 1, np.inf][np.argmin([div(b2, b1), c + 1, div(a1, a2)])]
            # print cons.intersection_time()
            # print cons.middle_time()
            # print cons.middle_time_cr()
            # print cons.cr_0(),cons.cr_inf()
            t_det_gen, rho_det_gen = cons.optDetStr()
            ratio = rho_det_gen / rho_det_anylsis
            # print ratio,rho_det_gen,rho_det_anylsis
            assert ((ratio <= eps_num + (1 + eps)) and (ratio >= -eps_num + 1 / (1 + eps)))

            if t_det_gen == np.inf and t_det_anylsis == np.inf:
                ratio = 0.0
            else:
                ratio = t_det_gen - t_det_anylsis

            assert ((ratio <= eps_num) and (ratio >= -eps_num))

            # we need to check cr_0,cr_inf

            # cons.plot(num_fig=i)
            # check norm - return new instance
            # check plot
            # print ""

    return True


def create_instance(br0, bb0, ar_s, ab_s, ts, c):
    lins_r, lins_b = [Linear(ar_s[0], br0)], [Linear(ab_s[0], bb0)]
    for ind, t in enumerate(ts[1:]):
        yr, yb = lins_r[-1][t], lins_b[-1][t]
        lins_r.append(Linear(ar_s[ind + 1], -ar_s[ind + 1] * t + yr))
        lins_b.append(Linear(ab_s[ind + 1], -ab_s[ind + 1] * t + yb))
    pr = PieceWiseLinear(lins_r, ts)
    pb = PieceWiseLinear(lins_b, ts)
    return InstancePieceWiseLinear(pr, pb, c)


def UT_PWL_5params():
    global cons
    print "Unit test piece-wise linear(5 parameters) ski rental problems"
    eps = 0.05
    eps_num = 1e-10
    rr = np.arange(0.0, 1.02, 0.1)
    num_iter = 100

    for i in range(num_iter):
        # print i
        ar0 = np.round(np.random.uniform(0.06, 1), 1)
        ab0 = np.round(np.random.uniform(0.0, ar0 - 0.07), 1)
        ar1 = np.round(np.random.uniform(0.2, 5.0), 1)
        ab1 = np.round(np.random.uniform(0.2, ar1 - 0.1), 1)
        c = max(ar0 - ab0, np.round(10 ** np.random.uniform(-0.3, 1.0), 1))

        # ar0,ab0,ar1,ab1,c =0.3,0.2,2.4,0.3,2.3
        # print ",".join([str(u) for u in [ar0,ab0,ar1,ab1,c]])

        cons = create_instance(1 - ar0, 1 - ab0, [ar0, ar1], [ab0, ab1], [0., 1.], c)
        plt = cons.plot(num_fig=1)
        plt.hold(True)
        # print cons.optDetStr()

        xx = np.arange(1.0, 20.0, 0.1)
        gg = np.array([cons.fr[x] / cons.fb[x] for x in xx])
        gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
        hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])

        # plt.plot(xx,np.array([cons.fr[x]/cons.fb[x] for x in xx]),"g")
        gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
        hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
        plt.plot(xx, gg, "g", xx, hh, "k")
        intr = cons.intersection()
        # print intr
        assert (np.abs(intr - 1.0) < eps_num)

        t_1 = cons.calc_t1(eps)
        # print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
        assert (-eps_num < t_1 < 1.0 + eps_num)
        assert ((cons.fr[t_1] + c) <= eps_num + (1 + eps) * (cons.f_rb(0) + c + cons.fb[t_1]))
        assert ((np.abs(t_1 - 1.0) < eps_num) or np.abs(
            cons.fr[t_1] + c - (1 + eps) * (cons.f_rb(0) + c + cons.fb[t_1])) < eps_num)

        alp = cons.calc_t_N_minus_1(eps)

        assert (alp >= 1.0)
        # print alp,cons.rho_inf(),cons.fr(alp)
        flag1 = (cons.fr[alp] / cons.fb[alp]) >= -eps_num + cons.cr_inf() / (1 + eps)
        # print alp,c,cons.fr[alp]
        flag2 = (c <= eps_num + eps * cons.fr[alp])
        # print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c

        assert (flag1 or flag2)

        cr0 = cons.cr_0()
        cr0_thr = div(1 - ab0, 1 - ar0)
        assert (np.abs(cr0_thr - cr0) < eps_num or (min(cr0, cr0_thr) == np.inf))

        cri = cons.cr_inf()

        cri_thr = div(ar1, ab1) if ar1 > 0 else cons.fr[2] / cons.fb[2]
        assert (np.abs(cri_thr - cri) < eps_num or (min(cri, cri_thr) == np.inf))
        # TODO - check det
        # we have cr_0 and cr_inf

        # fr = ar1(x-1)+1=ar1*x-ar1+1,fb=ab1(x-1)+1=ab1*x-ab1+1
        # (a1b2-a2b1)=ar1*(-ab1+1) - ab1*(-ar1+1)=ar1-ab1>0
        # g(z)= max(fr/fb,y<=z) = fr(z)/fb(z)
        # z_s = np.inf
        # min(h)
        # (a1b2-a2b1)=ar1*(-ab1+1) - ab1*(-ar1+1+c)=ar1-ab1-ab1*c
        hd = ar1 - ab1 - ab1 * c
        cr_1 = (cons.fr(1) + c) / cons.fb(1)
        if hd >= 0:
            zm = 1
            cr_zm = cr_1
        else:
            zm = None
            cr_zm = cri + eps  # because h is decreasing and cr_inf <= h(z) z->inf

        rho_det_anylsis = min(cr0, cr_zm, cri)
        t_det_anylsis = [0, zm, np.inf][np.argmin(np.array([cr0, cr_zm, cri]))]
        t_det_gen, rho_det_gen = cons.optDetStr()
        ratio = rho_det_gen / rho_det_anylsis
        assert ((ratio <= eps_num + (1 + eps)) and (ratio >= -eps_num + 1 / (1 + eps)))

        if t_det_gen == np.inf and t_det_anylsis == np.inf:
            ratio = 0.0
        else:
            ratio = t_det_gen - t_det_anylsis

        assert ((ratio <= eps_num) and (ratio >= -eps_num))

        tup = cons.approxOptStr(eps)
        ratio = tup[0] / rho_det_anylsis
        # print "ratios 5params",ratio,rho_det_anylsis,tup[0],t_det_anylsis
        assert (ratio <= eps_num + (1 + eps))
        # ((ratio<=eps_num+(1+eps)) and (ratio>= -eps_num+1/(1+eps)))
    return True


def UT_PWL_6params():
    global cons, gg, hh, xx
    print "Unit test piece-wise linear(6 parameters) ski rental problems"
    eps = 0.05
    eps_num = 1e-10
    rr = np.arange(0.0, 1.02, 0.1)
    num_iter = 100

    for i in range(num_iter):
        # print i
        ar0 = np.round(np.random.uniform(0.06, 1), 1)
        ab0 = np.round(np.random.uniform(0.0, ar0 - 0.07), 1)
        ar1 = np.round(np.random.uniform(0.2, 5.0), 1)
        ab1 = np.round(np.random.uniform(0.2, ar1 - 0.1), 1)
        c = max(ar0 - ab0, np.round(10 ** np.random.uniform(-0.3, 1.0), 1))
        t1 = np.round(np.random.uniform(1.2, 5.0), 1)

        # ar0,ab0,ar1,ab1,c,t1 = 0.5,0.3,0.2,0.2,0.8,1.6
        # print ",".join([str(u) for u in [ar0,ab0,ar1,ab1,c,t1]])

        cons = create_instance(1 - ar0, 1 - ab0, [ar0, ar0, ar1], [ab0, ab0, ab1], [0., 1., t1], c)
        cons = cons.normalize()
        plt = cons.plot()
        plt.hold(True)

        assert (np.max(np.abs(np.array(
            [cons.fr.lins[i](cons.times[i + 1]) - cons.fr.lins[i + 1](cons.times[i + 1]) for i in
             range(len(cons.times) - 1)]))) < eps)
        assert (np.max(np.abs(np.array(
            [cons.fb.lins[i](cons.times[i + 1]) - cons.fb.lins[i + 1](cons.times[i + 1]) for i in
             range(len(cons.times) - 1)]))) < eps)

        xx = np.arange(1.0, 50.0, 0.01)
        gg = np.array([cons.fr[x] / cons.fb[x] for x in xx])

        gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
        hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
        # plt.plot(xx,gg,"g",xx,hh,"k")
        # plt.plot(xx,np.array([cons.fr[x]/cons.fb[x] for x in xx]),"g")
        gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
        hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
        plt.plot(xx, gg, "g")
        plt.plot(xx, hh, "y")

        assert (np.abs(cons.intersection() - 1.0) < eps_num)

        t_1 = cons.calc_t1(eps)
        # print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
        assert (-eps_num < t_1 < 1.0 + eps_num)
        assert ((cons.fr[t_1] + c) <= eps_num + (1 + eps) * (cons.f_rb(0) + c + cons.fb[t_1]))
        assert ((np.abs(t_1 - 1.0) < eps_num) or np.abs(
            cons.fr[t_1] + c - (1 + eps) * (cons.f_rb(0) + c + cons.fb[t_1])) < eps_num)

        alp = cons.calc_t_N_minus_1(eps)

        assert (alp >= 1.0)
        # print alp,cons.rho_inf(),cons.fr(alp)
        flag1 = (cons.fr[alp] / cons.fb[alp]) >= -eps_num + cons.cr_inf() / (1 + eps)
        # print alp,c,cons.fr[alp]
        flag2 = (c <= eps_num + eps * cons.fr[alp])
        # print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c

        assert (flag1 or flag2)

        cr0 = cons.cr_0()
        cr0_thr = div(1 - ab0, 1 - ar0)
        assert (np.abs(cr0_thr - cr0) < eps_num or (min(cr0, cr0_thr) == np.inf))

        cri = cons.cr_inf()

        cri_thr = div(ar1, ab1) if ar1 > 0 else cons.fr[t1 + 10] / cons.fb[t1 + 10]
        # print cri_thr,cri
        assert (np.abs(cri_thr - cri) < eps_num or (min(cri, cri_thr) == np.inf))
        # TODO - check det
        # we have cr_0 and cr_inf

        # fr = ar0(x-1)+1=ar0*x-ar0+1,fb=ab0(x-1)+1=ab1*x-ab1+1
        # (a1b2-a2b1)=ar0*(-ab0+1) - ab0*(-ar0+1)=ar0-ab0>0
        # g(z)= max(fr/fb,y<=z) = fr(z)/fb(z)
        # 1<t<t1 h>g

        # fr = ar1(x-t1)+ar0(t1-1)+1=ar1*x+(ar0-ar1)*t1+1-ar0,fb=ab1*x+(ab0-ab1)*t1+1-ab0
        # (a1b2-a2b1)=ar1*((ab0-ab1)*t1+1-ab0) -ab1*((ar0-ar1)*t1+1-ar0)=ar1*((ab0)*t1+1-ab0) -ab1*((ar0)*t1+1-ar0)
        # g(z)= max(fr/fb,y<=z) = fr(z)/fb(z)
        # t1<t
        # if gd(t1)<0 =>g(>t1)=g(t1), else zs=inf

        # first we will check the slopes

        # we have 2 sections
        # 1<t<t1 g=fr/fb if gd down then g=g(1) else g=g
        # h = (fr+c)/fb
        """
        h=lambda t:((cons.fr(t)+c)/cons.fb(t))
        hd1 = (cons.fr.lins[1].slope*cons.fb.lins[1].offset)-(cons.fb.lins[1].slope*(cons.fr.lins[1].offset+c))

        gd1 = (cons.fr.lins[1].slope*cons.fb.lins[1].offset)-(cons.fb.lins[1].slope*cons.fr.lins[1].offset)
        if gd1<=0:
            g=lambda t:(cons.fr(1)/cons.fb(1))
        else:
            h=lambda t:(cons.fr(t)/cons.fb(t))

        if g(t1)>=max(h(1),h(t1)):
            pass
            #we have intersection here
        else:
            #the same for the next interval



        gd = ar1*((ab0)*t1+1-ab0) -ab1*((ar0)*t1+1-ar0)
        if gd>0:
            pass
            #does h monotone inreases ? 
            #hd = rd*b-(r+c)*bd=




        #z_s = np.inf
        #min(h)
        #(a1b2-a2b1)=ar1*(-ab1+1) - ab1*(-ar1+1+c)=ar1-ab1-ab1*c
        hd =ar1-ab1-ab1*c
        cr_1= (cons.fr(1)+c)/cons.fb(1)
        if hd>=0:
            zm=1
            cr_zm=cr_1
        else:
            zm=None
            cr_zm = cri+eps #because h is decreasing and cr_inf <= h(z) z->inf
        """

        # we will find cr_zm numeric
        zs_numer = xx[gg >= hh][0] if np.sum(gg >= hh) > 0 else np.inf
        zm_numer = xx[np.argmin(hh[xx <= zs_numer])] if (
            (cri > np.min(hh[xx <= zs_numer]) or zs_numer < np.inf)) else np.inf
        # cr_zm_numer = np.max(hh[xx<=zs_numer]) if ((cri>np.max(hh)or zs_numer<np.inf)) else cri
        zs = cons.intersection_time()
        zm, cr_zm = cons.middle_time_cr()
        print "zs", zs, zs_numer, zm, zm_numer, cri, cons.norm_time
        rho_det_anylsis = min(cr0, cr_zm, cri)
        t_det_anylsis = [0, zm, np.inf][np.argmin(np.array([cr0, cr_zm, cri]))]
        # print "det inline",rho_det_anylsis,t_det_anylsis,cr_zm,cr0,cri
        t_det_gen, rho_det_gen = cons.optDetStr()
        ratio = rho_det_gen / rho_det_anylsis
        assert ((ratio <= eps_num + (1 + eps)) and (ratio >= -eps_num + 1 / (1 + eps)))

        if t_det_gen == np.inf and t_det_anylsis == np.inf:
            ratio = 0.0
        else:
            ratio = t_det_gen - t_det_anylsis

        assert ((ratio <= eps_num) and (ratio >= -eps_num))
        # cons.approxOptStr(eps)
        tup = cons.approxOptStr(eps)
        ratio = tup[0] / rho_det_anylsis
        # print "ratios 6params",ratio,rho_det_anylsis,tup[0],t_det_anylsis
        # assert (ratio<= eps_num+(1+eps))
        # print "no assert of probalistic startegy at 6params"

    return True


def exmp_det():
    global cons, gg, hh, xx, xx_e, cr_xx_e, cons_org, c, td, br0, bb0, ar1, dict_res
    print "Examples of instances for learning about determistic strategy"
    eps = 0.05
    eps_num = 1e-10
    dict_res = {}
    bb0 = 50.0
    ar1 = 6.0
    genr = ((br0, c, td) for br0 in np.arange(10.0, 48.0, 2) for c in np.arange((bb0 - br0) + 0.2, bb0 - 0.1, 2) for td
            in np.arange((bb0 - br0) / ar1 + 0.1, 10.0, 0.5))
    for i, (br0, c, td) in enumerate(genr):
        """
        #hot
        br0=40
        bb0=50
        ts=[0,12,18]
        ar_s = [0,6,0]
        ab_s = [0,0,0]
        c = 10
        """
        """
        #P_mobile
        br0=35
        bb0=40
        ts=[0,15,20,25,30]
        ar_s = [0,5,5,0,0]
        ab_s = [0,0,5,5,0]
        c =40
        """

        assert ((td - eps_num) > ((bb0 - br0) / ar1)), "Bad TD " + ",".join(str(u) for u in [c, br0, td])
        # print ",".join(str(si) for si in [c,br0,td])
        assert (c > (bb0 - br0 - eps_num)), "Bad C " + ",".join(str(u) for u in [c, br0, td])

        ts = [0, 12, 12 + td]
        ar_s = [0, ar1, 0]
        ab_s = [0, 0, 0]

        cons_org = create_instance(br0, bb0, ar_s, ab_s, ts, c)
        cons = cons_org.normalize()
        # plt=cons.plot()
        # plt.hold(True)
        zs = cons.intersection_time()

        assert (np.max(np.abs(np.array(
            [cons.fr.lins[i](cons.times[i + 1]) - cons.fr.lins[i + 1](cons.times[i + 1]) for i in
             range(len(cons.times) - 1)]))) < eps)
        assert (np.max(np.abs(np.array(
            [cons.fb.lins[i](cons.times[i + 1]) - cons.fb.lins[i + 1](cons.times[i + 1]) for i in
             range(len(cons.times) - 1)]))) < eps)

        xx = np.arange(1.0, 20.0, 0.01)
        gg = np.array([cons.fr[x] / cons.fb[x] for x in xx])

        gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
        hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
        # plt.plot(xx,gg,"g",xx,hh,"k")
        # plt.plot(xx,np.array([cons.fr[x]/cons.fb[x] for x in xx]),"g")
        gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
        hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
        # plt.plot(xx,gg,"g")
        # plt.plot(xx,hh,"y")

        assert (np.abs(cons.intersection() - 1.0) < eps_num)

        t_1 = cons.calc_t1(eps)
        # print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
        assert (-eps_num < t_1 < 1.0 + eps_num)
        assert ((cons.fr[t_1] + c) <= eps_num + (1 + eps) * (cons.f_rb(0) + c + cons.fb[t_1]))
        # assert((np.abs(t_1-1.0)<eps_num) or np.abs(cons.fr[t_1]+c-(1+eps)*(cons.f_rb(0)+c+cons.fb[t_1]))<eps_num)
        # print "no assert of t_1"

        alp = cons.calc_t_N_minus_1(eps)

        assert (alp >= 1.0)
        # print alp,cons.rho_inf(),cons.fr(alp)
        flag1 = (cons.fr[alp] / cons.fb[alp]) >= -eps_num + cons.cr_inf() / (1 + eps)
        # print alp,c,cons.fr[alp]
        flag2 = (c <= eps_num + eps * cons.fr[alp])
        # print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c

        # assert(flag1 or flag2)

        cr0 = cons.cr_0()
        cr0_thr = div(bb0, br0)
        # assert(np.abs(cr0_thr-cr0)<eps_num or (min(cr0,cr0_thr)==np.inf ))
        # print "no assert of cr0"
        cri = cons.cr_inf()

        cri_thr = div(ar_s[-1], ab_s[-1]) if ar_s[-1] > 0 else cons.fr[20] / cons.fb[20]
        # print cri_thr,cri
        assert (np.abs(cri_thr - cri) < eps_num or (min(cri, cri_thr) == np.inf))

        zs_numer = xx[gg >= hh][0] if np.sum(gg >= hh) > 0 else np.inf
        zm_numer = xx[np.argmin(hh[xx <= zs_numer])] if (
            (cri > np.min(hh[xx <= zs_numer]) or zs_numer < np.inf)) else np.inf
        # cr_zm_numer = np.max(hh[xx<=zs_numer]) if ((cri>np.max(hh)or zs_numer<np.inf)) else cri
        zs = cons.intersection_time()
        zm, cr_zm = cons.middle_time_cr()
        # print "zs",zs,zs_numer,zm,zm_numer,cri,cons.norm_time
        rho_det_anylsis = min(cr0, cr_zm, cri)
        t_det_anylsis = [0, zm, np.inf][np.argmin(np.array([cr0, cr_zm, cri]))]
        # print "det inline",rho_det_anylsis,t_det_anylsis,cr_zm,cr0,cri
        t_det_gen, rho_det_gen = cons.optDetStr()
        print "det str", br0, c, td, rho_det_gen, t_det_gen

        # we will build based on XX the competitve ratio
        xx_e = np.arange(0.0, 20.0, 0.05)
        xx_e[-1] = np.inf
        cr_xx_e = [np.max([cons.cr(z, y) for y in xx_e]) for z in xx_e]
        # print "best cr and time",c,np.min(cr_xx_e),xx[np.argmin(cr_xx_e)],np.max([cons.cr(t_det_gen,y) for y in xx_e])
        diff = np.min(cr_xx_e) - np.max([cons.cr(t_det_gen, y) for y in xx_e])
        assert ((diff <= eps_num) and (diff >= -eps_num))
        ratio = rho_det_gen / rho_det_anylsis
        assert ((ratio <= eps_num + (1 + eps)) and (ratio >= -eps_num + 1 / (1 + eps)))

        if t_det_gen == np.inf and t_det_anylsis == np.inf:
            ratio = 0.0
        else:
            ratio = t_det_gen - t_det_anylsis

        assert ((ratio <= eps_num) and (ratio >= -eps_num))
        # cons.approxOptStr(eps)
        tup = cons.approxOptStr(eps)
        ratio = tup[0] / rho_det_anylsis
        # print "ratios 6params",ratio,rho_det_anylsis,tup[0],t_det_anylsis
        assert (ratio <= eps_num + (1 + eps)), "prob less good than det"
        # print "no assert of probalistic startegy at 6params"
        markersize = 20
        if False:
            if t_det_gen == 0:
                plt.plot(c, td, "b*", markersize=markersize)
            elif t_det_gen < np.inf:
                plt.plot(c, td, "gx", markersize=markersize)
            elif t_det_gen == np.inf:
                plt.plot(c, td, "ro", markersize=markersize)
            else:
                print "no det choice"
        """
        simArtists = [plt.Line2D((0,1),(0,0), color=c, marker=mark, linestyle='-') for c,mark in zip(["b","g","r"],["*","x","*"])]
        plt.legend(simArtists,["Strategy starts at buy","Strategy switch at intersection","Strategy never switch"])
        plt.ylabel("Slope of the rent option",fontsize=20)
        plt.xlabel("Switch cost",fontsize=20)
        plt.title("Best strategy as function of switch cost and slope",fontsize=30)
        """
        dict_res[c, td, br0] = (rho_det_gen, t_det_gen)
        if i % 100 == -1 % 100:
            print "saved", i
            pickle.dump(dict_res, open("dict_res_det_c_td_br0.pkl", "w+"))

    pickle.dump(dict_res, open("dict_res_det_c_td_br0.pkl", "w+"))

    return True


def exmp_det_for_middle_time_graph():
    global cons, gg, hh, xx, xx_e, cr_xx_e, cons_org, c, td, br0, bb0, ar1, dict_res
    print "Examples of instances for learning about determistic strategy"
    eps = 0.01
    eps_num = 1e-10
    dict_res = {}
    bb0 = 50.0
    ar1 = 5.0

    #g = max fr/fb (1+100t)/1
    #h (fr+c)/fb

    #when fr get more acelarate (slope bigger - we can see that g "catch" h)
    #cons_org = create_instance(0, 1, [1,2,0.1], [0,0,0.1], [0,1,2], 1) - zm finite but still h is 1
    cons_org = create_instance(0, 1, [1,2,0.6], [0,0,0.5], [0,1,1.3], 1)
    cons = cons_org.normalize()
    zm = cons.intersection_time()
    print(zm,cons.middle_time(),cons.optDetStr())
    tt=np.linspace(0,10,1000)
    fr_tt = np.array([cons.fr(t) for t in tt])
    fb_tt = np.array([cons.fb(t) for t in tt])
    plt.plot(tt,fr_tt,"r")
    plt.plot(tt,fb_tt,"b")
    fg_tt = fr_tt / np.clip(fb_tt, 0.01, None)
    fg_tt = [np.max(fg_tt[:(i + 1)]) for i in range(len(fr_tt))]
    plt.plot(tt, fg_tt , "g")
    plt.plot(tt,(fr_tt+1)/np.clip(fb_tt,0.01,None),"y")

    plt.legend(["Rent", "Buy", "g(z)", "h(z)"], loc=0, fontsize=20)
    plt.text(0.5, 4, "Switch cost:  1", fontsize=25)
    plt.text(5.8, 2, "middle time:  6.3", fontsize=25)
    plt.plot(6.3, 1.6, "k*",markersize=15)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)


    plt.show()

    return

    C1 = np.arange(2, 50, 1)
    T1 = np.arange(0.05, 20, 0.1)
    B1 = np.arange(2, 50, 0.1)
    genr = ((br0, c, td) for br0 in B1 for c in np.arange((bb0 - br0) + 0.2, bb0 - 0.1, 2) for td
            in np.arange((bb0 - br0) / ar1 + 0.1, 10.0, 0.2))
    for i, (br0, c, td) in enumerate(genr):


        assert ((td - eps_num) > ((bb0 - br0) / ar1)), "Bad TD " + ",".join(str(u) for u in [c, br0, td])
        # print ",".join(str(si) for si in [c,br0,td])
        assert (c > (bb0 - br0 - eps_num)), "Bad C " + ",".join(str(u) for u in [c, br0, td])

        ts = [0, 10, 10 + td]
        ar_s = [0, ar1, 0]
        ab_s = [0, 0, 0]

        cons_org = create_instance(br0, bb0, ar_s, ab_s, ts, c)
        cons = cons_org.normalize()
        # plt=cons.plot()
        # plt.hold(True)

        zm = cons.middle_time()
        if zm==np.inf:
            continue
        if zm==1.0:
            continue
        assert (np.max(np.abs(np.array(
            [cons.fr.lins[i](cons.times[i + 1]) - cons.fr.lins[i + 1](cons.times[i + 1]) for i in
             range(len(cons.times) - 1)]))) < eps)
        assert (np.max(np.abs(np.array(
            [cons.fb.lins[i](cons.times[i + 1]) - cons.fb.lins[i + 1](cons.times[i + 1]) for i in
             range(len(cons.times) - 1)]))) < eps)

        xx = np.arange(1.0, 20.0, 0.01)
        gg = np.array([cons.fr[x] / cons.fb[x] for x in xx])

        gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
        hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
        # plt.plot(xx,gg,"g",xx,hh,"k")
        # plt.plot(xx,np.array([cons.fr[x]/cons.fb[x] for x in xx]),"g")
        gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
        hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
        # plt.plot(xx,gg,"g")
        # plt.plot(xx,hh,"y")

        assert (np.abs(cons.intersection() - 1.0) < eps_num)

        t_1 = cons.calc_t1(eps)
        # print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
        assert (-eps_num < t_1 < 1.0 + eps_num)
        assert ((cons.fr[t_1] + c) <= eps_num + (1 + eps) * (cons.f_rb(0) + c + cons.fb[t_1]))
        # assert((np.abs(t_1-1.0)<eps_num) or np.abs(cons.fr[t_1]+c-(1+eps)*(cons.f_rb(0)+c+cons.fb[t_1]))<eps_num)
        # print "no assert of t_1"

        alp = cons.calc_t_N_minus_1(eps)

        assert (alp >= 1.0)
        # print alp,cons.rho_inf(),cons.fr(alp)
        flag1 = (cons.fr[alp] / cons.fb[alp]) >= -eps_num + cons.cr_inf() / (1 + eps)
        # print alp,c,cons.fr[alp]
        flag2 = (c <= eps_num + eps * cons.fr[alp])
        # print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c

        # assert(flag1 or flag2)

        cr0 = cons.cr_0()
        cr0_thr = div(bb0, br0)
        # assert(np.abs(cr0_thr-cr0)<eps_num or (min(cr0,cr0_thr)==np.inf ))
        # print "no assert of cr0"
        cri = cons.cr_inf()

        cri_thr = div(ar_s[-1], ab_s[-1]) if ar_s[-1] > 0 else cons.fr[20] / cons.fb[20]
        # print cri_thr,cri
        assert (np.abs(cri_thr - cri) < eps_num or (min(cri, cri_thr) == np.inf))

        zs_numer = xx[gg >= hh][0] if np.sum(gg >= hh) > 0 else np.inf
        zm_numer = xx[np.argmin(hh[xx <= zs_numer])] if (
            (cri > np.min(hh[xx <= zs_numer]) or zs_numer < np.inf)) else np.inf
        # cr_zm_numer = np.max(hh[xx<=zs_numer]) if ((cri>np.max(hh)or zs_numer<np.inf)) else cri
        zs = cons.intersection_time()
        zm, cr_zm = cons.middle_time_cr()
        # print "zs",zs,zs_numer,zm,zm_numer,cri,cons.norm_time
        rho_det_anylsis = min(cr0, cr_zm, cri)
        t_det_anylsis = [0, zm, np.inf][np.argmin(np.array([cr0, cr_zm, cri]))]
        # print "det inline",rho_det_anylsis,t_det_anylsis,cr_zm,cr0,cri
        t_det_gen, rho_det_gen = cons.optDetStr()
        print "det str", br0, c, td, rho_det_gen, t_det_gen

        # we will build based on XX the competitve ratio
        xx_e = np.arange(0.0, 20.0, 0.05)
        xx_e[-1] = np.inf
        cr_xx_e = [np.max([cons.cr(z, y) for y in xx_e]) for z in xx_e]
        # print "best cr and time",c,np.min(cr_xx_e),xx[np.argmin(cr_xx_e)],np.max([cons.cr(t_det_gen,y) for y in xx_e])
        diff = np.min(cr_xx_e) - np.max([cons.cr(t_det_gen, y) for y in xx_e])
        assert ((diff <= eps_num) and (diff >= -eps_num))
        ratio = rho_det_gen / rho_det_anylsis
        assert ((ratio <= eps_num + (1 + eps)) and (ratio >= -eps_num + 1 / (1 + eps)))

        if t_det_gen == np.inf and t_det_anylsis == np.inf:
            ratio = 0.0
        else:
            ratio = t_det_gen - t_det_anylsis

        assert ((ratio <= eps_num) and (ratio >= -eps_num))
        # cons.approxOptStr(eps)
        tup = cons.approxOptStr(eps)
        ratio = tup[0] / rho_det_anylsis
        # print "ratios 6params",ratio,rho_det_anylsis,tup[0],t_det_anylsis
        assert (ratio <= eps_num + (1 + eps)), "prob less good than det"
        # print "no assert of probalistic startegy at 6params"
        markersize = 20
        if False:
            if t_det_gen == 0:
                plt.plot(c, td, "b*", markersize=markersize)
            elif t_det_gen < np.inf:
                plt.plot(c, td, "gx", markersize=markersize)
            elif t_det_gen == np.inf:
                plt.plot(c, td, "ro", markersize=markersize)
            else:
                print "no det choice"
        """
        simArtists = [plt.Line2D((0,1),(0,0), color=c, marker=mark, linestyle='-') for c,mark in zip(["b","g","r"],["*","x","*"])]
        plt.legend(simArtists,["Strategy starts at buy","Strategy switch at intersection","Strategy never switch"])
        plt.ylabel("Slope of the rent option",fontsize=20)
        plt.xlabel("Switch cost",fontsize=20)
        plt.title("Best strategy as function of switch cost and slope",fontsize=30)
        """
        dict_res[c, td, br0] = (rho_det_gen, t_det_gen)
        if i % 100 == -1 % 100:
            print "saved", i
            pickle.dump(dict_res, open("dict_res_det_c_td_br0.pkl", "w+"))

    pickle.dump(dict_res, open("dict_res_det_c_td_br0.pkl", "w+"))

    return True

def exmp_prob():
    global cons, gg, hh, xx, cons_org, c, td, br0, bb0, ar1, dict_res, C, f, A_eq, b_ub, b_eq, q, T, t_det_anylsis, rho_det_anylsis, t_det_gen, rho_det_gen
    print "Examples of instances for learning about probabilties strategy"
    eps = 0.05
    eps_num = 1e-10
    num_iter = 10
    dict_res = {}
    bb0 = 50.0
    ar1 = 6.0
    genr = ((br0, c, td) for br0 in np.arange(10.0, 48.0, 2) for c in np.arange((bb0 - br0) + 0.2, bb0 - 0.1, 2) for td
            in np.arange((bb0 - br0) / ar1 + 0.1, 10.0, 0.5))
    for i, (br0, c, td) in enumerate(genr):
        # br0,c,td=10.0,40.2,6.76666666667
        assert (td - eps_num) > ((bb0 - br0) / ar1), "BAD TD" + ",".join(str(si) for si in [br0, c, td])

        assert (c > (bb0 - br0 - eps_num)), "BAD C" + ",".join(str(si) for si in [br0, c, td])

        ts = [0, 12, 12 + td]
        ar_s = [0, ar1, 0]
        ab_s = [0, 0, 0]

        cons_org = create_instance(br0, bb0, ar_s, ab_s, ts, c)
        cons = cons_org.normalize()
        # plt=cons.plot()
        # plt.hold(True)

        assert (np.max(np.abs(np.array(
            [cons.fr.lins[i](cons.times[i + 1]) - cons.fr.lins[i + 1](cons.times[i + 1]) for i in
             range(len(cons.times) - 1)]))) < eps), "BAD fr" + ",".join(str(si) for si in [br0, c, td])
        assert (np.max(np.abs(np.array(
            [cons.fb.lins[i](cons.times[i + 1]) - cons.fb.lins[i + 1](cons.times[i + 1]) for i in
             range(len(cons.times) - 1)]))) < eps), "BAD fb" + ",".join(str(si) for si in [br0, c, td])

        xx = np.arange(1.0, 20.0, 0.01)
        gg = np.array([cons.fr[x] / cons.fb[x] for x in xx])

        gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
        hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
        # plt.plot(xx,gg,"g",xx,hh,"k")
        # plt.plot(xx,np.array([cons.fr[x]/cons.fb[x] for x in xx]),"g")
        gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
        hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
        # plt.plot(xx,gg,"g")
        # plt.plot(xx,hh,"y")

        assert (np.abs(cons.intersection() - 1.0) < eps_num), "BAD intersection" + ",".join(
            str(si) for si in [br0, c, td])

        t_1 = cons.calc_t1(eps)
        # print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
        assert (-eps_num < t_1 < 1.0 + eps_num), "BAD t1_1 " + ",".join(str(si) for si in [br0, c, td])
        assert ((cons.fr[t_1] + cons.c) <= eps_num + (1 + eps) * (
                cons.f_rb(0) + cons.c + cons.fb[t_1])), "BAD t1_2 " + ",".join(str(si) for si in [br0, c, td])
        tol = 0.5
        assert ((np.abs(t_1 - 1.0) < eps_num) or np.abs(cons.fr[t_1] + cons.c - (1 + eps) * (
                cons.f_rb(0) + cons.c + cons.fb[t_1])) < eps_num), "BAD t1_3 " + ",".join(
            str(si) for si in [br0, c, td]) + "|" + str(t_1) + "|" + str(
            np.abs(cons.fr[t_1] + c - (1 + eps) * (cons.f_rb(0) + c + cons.fb[t_1]))) + "|" + str(
            cons.fr[1.0] + c) + "|" + str((1 + eps) * (cons.f_rb(0) + c + cons.fb[1.0]))
        # print "no assert of t_1"

        alp = cons.calc_t_N_minus_1(eps)

        assert (alp >= 1.0), "BAD alp" + ",".join(str(si) for si in [br0, c, td])
        # print alp,cons.rho_inf(),cons.fr(alp)
        flag1 = (cons.fr[alp] / cons.fb[alp]) >= -eps_num + cons.cr_inf() / (1 + eps)
        # print alp,c,cons.fr[alp]
        flag2 = (cons.c <= eps_num + eps * cons.fr[alp])
        # print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c

        assert (flag1 or flag2), "BAD flags" + ",".join(str(si) for si in [br0, c, td])

        cr0 = cons.cr_0()
        cr0_thr = div(bb0, br0)
        assert (np.abs(cr0_thr - cr0) < eps_num or (min(cr0, cr0_thr) == np.inf)), "BAD cr0" + ",".join(
            str(si) for si in [br0, c, td])

        cri = cons.cr_inf()

        cri_thr = div(ar_s[-1], ab_s[-1]) if ar_s[-1] > 0 else cons.fr[20] / cons.fb[20]
        # print cri_thr,cri
        assert (np.abs(cri_thr - cri) < eps_num or (min(cri, cri_thr) == np.inf)), "BAD cri" + ",".join(
            str(si) for si in [br0, c, td])

        zs_numer = xx[gg >= hh][0] if np.sum(gg >= hh) > 0 else np.inf
        zm_numer = xx[np.argmin(hh[xx <= zs_numer])] if (
            (cri > np.min(hh[xx <= zs_numer]) or zs_numer < np.inf)) else np.inf
        # cr_zm_numer = np.max(hh[xx<=zs_numer]) if ((cri>np.max(hh)or zs_numer<np.inf)) else cri
        zs = cons.intersection_time()
        zm, cr_zm = cons.middle_time_cr()
        # print "zs",zs,zs_numer,zm,zm_numer,cri,cons.norm_time
        rho_det_anylsis = min(cr0, cr_zm, cri)
        t_det_anylsis = [0, zm, np.inf][np.argmin(np.array([cr0, cr_zm, cri]))]
        # print "det inline",rho_det_anylsis,t_det_anylsis,cr_zm,cr0,cri
        t_det_gen, rho_det_gen = cons.optDetStr()
        print "prob str", c, td, br0

        # we will build based on XX the competitve ratio
        xx_e = np.arange(0.0, 20.0, 0.05)
        xx_e[-1] = np.inf
        cr_xx_e = [np.max([cons.cr(z, y) for y in xx_e]) for z in xx_e]
        # print "best cr and time",c,np.min(cr_xx_e),xx[np.argmin(cr_xx_e)],np.max([cons.cr(t_det_gen,y) for y in xx_e])
        diff = np.min(cr_xx_e) - np.max([cons.cr(t_det_gen, y) for y in xx_e])
        eps_diff = 0.02
        assert (np.abs(diff) <= eps_diff), "bad opt det" + ",".join(str(si) for si in [br0, c, td]) + "|diff= " + str(
            diff)
        ratio = rho_det_gen / rho_det_anylsis
        assert ((ratio <= eps_num + (1 + eps)) and (ratio >= -eps_num + 1 / (1 + eps))), "bad opt det" + ",".join(
            str(si) for si in [br0, c, td])

        if t_det_gen == np.inf and t_det_anylsis == np.inf:
            ratio = 0.0
        else:
            ratio = t_det_gen - t_det_anylsis

        assert ((ratio <= eps_num) and (ratio >= -eps_num)), "bad opt det ratio" + ",".join(
            str(si) for si in [br0, c, td])
        # cons.approxOptStr(eps)
        tup = cons.approxOptStr(eps)
        q = tup[1]
        T = tup[2]
        C, f, A_eq, b_ub, b_eq = tup[4]
        ratio = tup[0] / rho_det_anylsis
        # print "ratios 6params",ratio,rho_det_anylsis,tup[0],t_det_anylsis
        assert (ratio <= eps_num + (1 + eps)), "prob less good then det|ratio={}|det_cr={}|prob_cr={}".format(ratio,
                                                                                                              rho_det_anylsis,
                                                                                                              tup[
                                                                                                                  0]) + "|" + ",".join(
            str(si) for si in [br0, c, td])
        # TODO - make some random probs for seeing if the LP works good
        markersize = 20
        if False:
            if t_det_gen == 0:
                plt.plot(c, td, "b*", markersize=markersize)
            elif t_det_gen < np.inf:
                plt.plot(c, td, "gx", markersize=markersize)
            elif t_det_gen == np.inf:
                plt.plot(c, td, "ro", markersize=markersize)
            else:
                print "no det choice"
        # dict_res[c,td,ar1,br0,bb0]
        dict_res[c, td, br0] = (tup[0], q[0], q[-1], q, len(T))
        if i == -1:
            print "saved", i
            pickle.dump(dict_res, open("dict_res_prob_c_td_ar1_br0.pkl", "w+"))

    return True


def exmp_unifed():
    global cons, gg, hh, xx, cons_org, c, td, br0, bb0, ar1, dict_res, C, f, A_eq, b_ub, b_eq, q, T, t_det_anylsis, rho_det_anylsis, t_det_gen, rho_det_gen
    print "Examples of instances for learning about probabilties strategy"
    eps = 0.05
    eps_num = 1e-10
    dict_res = {}
    bb0 = 50.0
    ar1 = 5.0
    genr = ((br0, c, td) for br0 in np.arange(10.0, 48.0, 2) for c in np.arange((bb0 - br0) + 0.2, bb0 - 0.1, 2) for td
            in np.arange(0.5, 10.0, 0.5))
    for ind, (br0, c, td) in enumerate(genr):
        if ((td - eps_num) < ((bb0 - br0) / ar1)):
            dict_res_slice[br0] = [np.inf, 1, None, None, None]
            continue
        assert (td - eps_num) > ((bb0 - br0) / ar1), "BAD TD" + ",".join(str(si) for si in [br0, c, td])

        assert (c > (bb0 - br0 - eps_num)), "BAD C" + ",".join(str(si) for si in [br0, c, td])

        ts = [0, 10, 10 + td]
        ar_s = [0, ar1, 0]
        ab_s = [0, 0, 0]

        cons_org = create_instance(br0, bb0, ar_s, ab_s, ts, c)
        cons = cons_org.normalize()
        if ind == (-1):
            plt = cons.plot()
            plt.hold(True)
            plt.axis([0.0, 1.2, 0, 1.2])

            plt.annotate('Initial rent cost', xy=(0.3, 0.24), xytext=(0.3, 0.03),
                         arrowprops=dict(facecolor='black', shrink=1),
                         )
            plt.annotate('Time diffrnece', xy=(1.03, 0.6), xytext=(0.52, 0.6),
                         arrowprops=dict(facecolor='black', shrink=1),
                         )
            plt.show()
            break

        assert (np.max(np.abs(np.array(
            [cons.fr.lins[i](cons.times[i + 1]) - cons.fr.lins[i + 1](cons.times[i + 1]) for i in
             range(len(cons.times) - 1)]))) < eps), "BAD fr" + ",".join(str(si) for si in [br0, c, td])
        assert (np.max(np.abs(np.array(
            [cons.fb.lins[i](cons.times[i + 1]) - cons.fb.lins[i + 1](cons.times[i + 1]) for i in
             range(len(cons.times) - 1)]))) < eps), "BAD fb" + ",".join(str(si) for si in [br0, c, td])

        xx = np.arange(1.0, 20.0, 0.01)
        gg = np.array([cons.fr[x] / cons.fb[x] for x in xx])

        gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
        hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
        # plt.plot(xx,gg,"g",xx,hh,"k")
        # plt.plot(xx,np.array([cons.fr[x]/cons.fb[x] for x in xx]),"g")
        gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
        hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
        # plt.plot(xx,gg,"g")
        # plt.plot(xx,hh,"y")

        assert (np.abs(cons.intersection() - 1.0) < eps_num), "BAD intersection" + ",".join(
            str(si) for si in [br0, c, td])

        t_1 = cons.calc_t1(eps)
        # print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
        assert (-eps_num < t_1 < 1.0 + eps_num), "BAD t1_1 " + ",".join(str(si) for si in [br0, c, td])
        assert ((cons.fr[t_1] + cons.c) <= eps_num + (1 + eps) * (
                cons.f_rb(0) + cons.c + cons.fb[t_1])), "BAD t1_2 " + ",".join(str(si) for si in [br0, c, td])
        tol = 0.5
        assert ((np.abs(t_1 - 1.0) < eps_num) or np.abs(cons.fr[t_1] + cons.c - (1 + eps) * (
                cons.f_rb(0) + cons.c + cons.fb[t_1])) < eps_num), "BAD t1_3 " + ",".join(
            str(si) for si in [br0, c, td]) + "|" + str(t_1) + "|" + str(
            np.abs(cons.fr[t_1] + c - (1 + eps) * (cons.f_rb(0) + c + cons.fb[t_1]))) + "|" + str(
            cons.fr[1.0] + c) + "|" + str((1 + eps) * (cons.f_rb(0) + c + cons.fb[1.0]))
        # print "no assert of t_1"

        alp = cons.calc_t_N_minus_1(eps)

        assert (alp >= 1.0), "BAD alp" + ",".join(str(si) for si in [br0, c, td])
        # print alp,cons.rho_inf(),cons.fr(alp)
        flag1 = (cons.fr[alp] / cons.fb[alp]) >= -eps_num + cons.cr_inf() / (1 + eps)
        # print alp,c,cons.fr[alp]
        flag2 = (cons.c <= eps_num + eps * cons.fr[alp])
        # print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c

        assert (flag1 or flag2), "BAD flags" + ",".join(str(si) for si in [br0, c, td])

        cr0 = cons.cr_0()
        cr0_thr = div(bb0, br0)
        assert (np.abs(cr0_thr - cr0) < eps_num or (min(cr0, cr0_thr) == np.inf)), "BAD cr0" + ",".join(
            str(si) for si in [br0, c, td])

        cri = cons.cr_inf()

        cri_thr = div(ar_s[-1], ab_s[-1]) if ar_s[-1] > 0 else cons.fr[20] / cons.fb[20]
        # print cri_thr,cri
        assert (np.abs(cri_thr - cri) < eps_num or (min(cri, cri_thr) == np.inf)), "BAD cri" + ",".join(
            str(si) for si in [br0, c, td])

        zs_numer = xx[gg >= hh][0] if np.sum(gg >= hh) > 0 else np.inf
        zm_numer = xx[np.argmin(hh[xx <= zs_numer])] if (
            (cri > np.min(hh[xx <= zs_numer]) or zs_numer < np.inf)) else np.inf
        # cr_zm_numer = np.max(hh[xx<=zs_numer]) if ((cri>np.max(hh)or zs_numer<np.inf)) else cri
        zs = cons.intersection_time()
        zm, cr_zm = cons.middle_time_cr()
        # print "zs",zs,zs_numer,zm,zm_numer,cri,cons.norm_time
        rho_det_anylsis = min(cr0, cr_zm, cri)
        t_det_anylsis = [0, zm, np.inf][np.argmin(np.array([cr0, cr_zm, cri]))]
        # print "det inline",rho_det_anylsis,t_det_anylsis,cr_zm,cr0,cri
        t_det_gen, rho_det_gen = cons.optDetStr()
        # print "prob str", c, td, br0

        # we will build based on XX the competitve ratio
        xx_e = np.arange(0.0, 20.0, 0.05)
        xx_e[-1] = np.inf
        cr_xx_e = [np.max([cons.cr(z, y) for y in xx_e]) for z in xx_e]
        # print "best cr and time",c,np.min(cr_xx_e),xx[np.argmin(cr_xx_e)],np.max([cons.cr(t_det_gen,y) for y in xx_e])
        diff = np.min(cr_xx_e) - np.max([cons.cr(t_det_gen, y) for y in xx_e])
        eps_diff = 0.02
        assert (np.abs(diff) <= eps_diff), "bad opt det" + ",".join(str(si) for si in [br0, c, td]) + "|diff= " + str(
            diff)
        ratio = rho_det_gen / rho_det_anylsis
        assert ((ratio <= eps_num + (1 + eps)) and (ratio >= -eps_num + 1 / (1 + eps))), "bad opt det" + ",".join(
            str(si) for si in [br0, c, td])

        if t_det_gen == np.inf and t_det_anylsis == np.inf:
            ratio = 0.0
        else:
            ratio = t_det_gen - t_det_anylsis

        assert ((ratio <= eps_num) and (ratio >= -eps_num)), "bad opt det ratio" + ",".join(
            str(si) for si in [br0, c, td])
        # cons.approxOptStr(eps)
        tup = cons.approxOptStr(eps)
        rho_prob = tup[0]
        q = tup[1]
        T = tup[2]
        C, f, A_eq, b_ub, b_eq = tup[4]
        ratio = tup[0] / rho_det_anylsis
        # print "ratios 6params",ratio,rho_det_anylsis,tup[0],t_det_anylsis
        assert (ratio <= eps_num + (1 + eps)), "prob less good then det|ratio={}|det_cr={}|prob_cr={}".format(ratio,
                                                                                                              rho_det_anylsis,
                                                                                                              tup[
                                                                                                                  0]) + "|" + ",".join(
            str(si) for si in [br0, c, td])
        # TODO - make some random probs for seeing if the LP works good
        markersize = 20
        if False:
            if t_det_gen == 0:
                plt.plot(c, td, "b*", markersize=markersize)
            elif t_det_gen < np.inf:
                plt.plot(c, td, "gx", markersize=markersize)
            elif t_det_gen == np.inf:
                plt.plot(c, td, "ro", markersize=markersize)
            else:
                print "no det choice"
        # dict_res[c,td,ar1,br0,bb0]
        dict_res[br0, c, td] = (t_det_gen, rho_det_gen, rho_prob, q, T)
        if ind % 100 == 0:
            print "saved", ind
            pickle.dump(dict_res, open("dict_res_uniified_br0_c_td.pkl", "w+"))
        pickle.dump(dict_res, open("import dict_res_uniified_br0_c_td.pkl", "w+"))
    return True


def exmp_p6():
    # TODO
    # add switch cost to text
    # add middle time to text (which is optimal)
    global cons, gg, hh, xx, cons_org, c, td, br0, bb0, ar1, dict_res, cr_xx_e, xx_e
    print "Example of instance for determinstic strategy "
    eps = 0.05
    eps_num = 1e-10
    num_iter = 1

    # $a_r^1,a_b^1,a_r^2,a_b^2,a_r^3 = 0.5,0,2,1,1$

    # $t_1,t_2,c=2 ,5 ,1$}

    br0, bb0 = 0.5, 1
    ts = [0, 2, 5]
    ar_s = [0.5, 2, 2]
    ab_s = [0, 1, 1]
    # if ar_s[-1] = 2 then optimal is 2.0 but z_m=inf else z_m =6.15 but it is not optimal
    c = 1

    cons_org = create_instance(br0, bb0, ar_s, ab_s, ts, c)
    cons = cons_org.normalize()
    plt = cons.plot()
    plt.hold(True)

    assert (np.max(np.abs(np.array([cons.fr.lins[i](cons.times[i + 1]) - cons.fr.lins[i + 1](cons.times[i + 1]) for i in
                                    range(len(cons.times) - 1)]))) < eps), "BAD fr" + ",".join(
        str(si) for si in [br0, c, td])
    assert (np.max(np.abs(np.array([cons.fb.lins[i](cons.times[i + 1]) - cons.fb.lins[i + 1](cons.times[i + 1]) for i in
                                    range(len(cons.times) - 1)]))) < eps), "BAD fb" + ",".join(
        str(si) for si in [br0, c, td])

    xx = np.arange(1.0, 20.0, 0.01)
    gg = np.array([cons.fr[x] / cons.fb[x] for x in xx])

    gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
    hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
    # plt.plot(xx,gg,"g",xx,hh,"k")
    # plt.plot(xx,np.array([cons.fr[x]/cons.fb[x] for x in xx]),"g")
    gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
    hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
    plt.plot(xx, gg, "g")
    plt.plot(xx, hh, "y")

    plt.legend(["Rent", "Buy", "g(z)", "h(z)"], loc=0, fontsize=20)
    plt.text(0.5, 1.33, "Switch cost:  1", fontsize=20)
    plt.text(1.5, 2.7, "middle time:  2", fontsize=20)
    plt.plot(2, 2.5, "y*")
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    assert (np.abs(cons.intersection() - 1.0) < eps_num), "BAD intersection" + ",".join(str(si) for si in [br0, c, td])

    t_1 = cons.calc_t1(eps)
    # print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
    assert (-eps_num < t_1 < 1.0 + eps_num), "BAD t1_1 " + ",".join(str(si) for si in [br0, c, td])
    assert ((cons.fr[t_1] + cons.c) <= eps_num + (1 + eps) * (
            cons.f_rb(0) + cons.c + cons.fb[t_1])), "BAD t1_2 " + ",".join(str(si) for si in [br0, c, td])
    tol = 0.5
    assert ((np.abs(t_1 - 1.0) < eps_num) or np.abs(
        cons.fr[t_1] + cons.c - (1 + eps) * (cons.f_rb(0) + cons.c + cons.fb[t_1])) < eps_num), "BAD t1_3 " + ",".join(
        str(si) for si in [br0, c, td]) + "|" + str(t_1) + "|" + str(
        np.abs(cons.fr[t_1] + c - (1 + eps) * (cons.f_rb(0) + c + cons.fb[t_1]))) + "|" + str(
        cons.fr[1.0] + c) + "|" + str((1 + eps) * (cons.f_rb(0) + c + cons.fb[1.0]))
    # print "no assert of t_1"

    alp = cons.calc_t_N_minus_1(eps)

    assert (alp >= 1.0), "BAD alp" + ",".join(str(si) for si in [br0, c, td])
    # print alp,cons.rho_inf(),cons.fr(alp)
    flag1 = (cons.fr[alp] / cons.fb[alp]) >= -eps_num + cons.cr_inf() / (1 + eps)
    # print alp,c,cons.fr[alp]
    flag2 = (cons.c <= eps_num + eps * cons.fr[alp])
    # print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c

    assert (flag1 or flag2), "BAD flags" + ",".join(str(si) for si in [br0, c, td])

    cr0 = cons.cr_0()
    cr0_thr = div(bb0, br0)
    assert (np.abs(cr0_thr - cr0) < eps_num or (min(cr0, cr0_thr) == np.inf)), "BAD cr0" + ",".join(
        str(si) for si in [br0, c, td])

    cri = cons.cr_inf()

    cri_thr = div(ar_s[-1], ab_s[-1]) if ar_s[-1] > 0 else cons.fr[20] / cons.fb[20]
    # print cri_thr,cri
    assert (np.abs(cri_thr - cri) < eps_num or (min(cri, cri_thr) == np.inf)), "BAD cri" + ",".join(
        str(si) for si in [br0, c, td])

    zs_numer = xx[gg >= hh][0] if np.sum(gg >= hh) > 0 else np.inf
    g_zs_numer = gg[gg >= hh][0] if np.sum(gg >= hh) > 0 else np.inf
    zm_numer = xx[np.argmin(hh[xx <= zs_numer])] if (
        (cri > np.min(hh[xx <= zs_numer]) or zs_numer < np.inf)) else np.inf
    # cr_zm_numer = np.max(hh[xx<=zs_numer]) if ((cri>np.max(hh)or zs_numer<np.inf)) else cri
    zs = cons.intersection_time()
    zm, cr_zm = cons.middle_time_cr()
    print "zs", zs, zs_numer, zm, zm_numer, cri, cons.norm_time
    plt.plot(zs, g_zs_numer, "*")
    print zs, g_zs_numer
    rho_det_anylsis = min(cr0, cr_zm, cri)
    t_det_anylsis = [0, zm, np.inf][np.argmin(np.array([cr0, cr_zm, cri]))]
    # print "det inline",rho_det_anylsis,t_det_anylsis,cr_zm,cr0,cri
    t_det_gen, rho_det_gen = cons.optDetStr()
    print "det str", rho_det_gen, t_det_gen

    # we will build based on XX the competitve ratio
    xx_e = np.arange(0.0, 20.0, 0.05)
    xx_e[-1] = np.inf
    cr_xx_e = [np.max([cons.cr(z, y) for y in xx_e]) for z in xx_e]
    # print "best cr and time",c,np.min(cr_xx_e),xx[np.argmin(cr_xx_e)],np.max([cons.cr(t_det_gen,y) for y in xx_e])
    diff = np.min(cr_xx_e) - np.max([cons.cr(t_det_gen, y) for y in xx_e])
    diff = 0 if np.min(cr_xx_e) == np.inf and np.max([cons.cr(t_det_gen, y) for y in xx_e]) == np.inf else diff
    assert ((diff <= eps_num) and (diff >= -eps_num)), str(diff) + "|" + str(np.min(cr_xx_e)) + "|" + str(
        np.max([cons.cr(t_det_gen, y) for y in xx_e]))
    ratio = rho_det_gen / rho_det_anylsis
    assert ((ratio <= eps_num + (1 + eps)) and (ratio >= -eps_num + 1 / (1 + eps))), "bad opt det" + ",".join(
        str(si) for si in [br0, c, td])

    if t_det_gen == np.inf and t_det_anylsis == np.inf:
        ratio = 0.0
    else:
        ratio = t_det_gen - t_det_anylsis

    assert ((ratio <= eps_num) and (ratio >= -eps_num)), "bad opt det ratio" + ",".join(str(si) for si in [br0, c, td])
    # cons.approxOptStr(eps)
    tup = cons.approxOptStr(eps)
    ratio = tup[0] / rho_det_anylsis
    # print "ratios 6params",ratio,rho_det_anylsis,tup[0],t_det_anylsis
    # assert (ratio<= eps_num+(1+eps)),"prob less good then det"
    print "assert  PROB!!!"
    plt.axis([0, 5, 0, 8])
    plt.show()
    return True


# ar1,ab1 - control the parameters before 1
# ar2,ab2 - control the parametres after 1 unti t1
# ar3 - the slope of rent after t2 - can be 0 or ab2
# t1 - the time when rent decrease his slope to ab2
# t2 - the time when buy decreses his slope to 0
# c - swith cost

def UT_PWL_8params():
    global cons
    print "Unit test linear ski rental problems"
    eps = 0.05
    eps_num = 1e-10
    num_iter = 100

    for i in range(num_iter):
        print i
        ar1 = np.round(np.random.uniform(0.06, 1), 1)
        ab1 = np.round(np.random.uniform(0.0, ar1 - 0.07), 1)
        ar2 = np.round(np.random.uniform(0.2, 5.0), 1)
        ab2 = np.round(np.random.uniform(0.2, ar2 - 0.1), 1)
        c = min(ar1 - ab1, np.round(10 ** np.random.uniform(-0.3, 1.0), 1))
        ar3 = np.random.choice([0.0, ab2])
        t1 = np.round(np.random.uniform(1.2, 5.0), 1)
        # TODO - make sure t2 is before intersection
        t2 = np.round(np.random.uniform(t1 + 0.1, 5.0), 1)
        ab3 = ab2  # 0.0
        # ar1,ab1,ar2,ab2,ar3,t1,t2,c = 0.5, 0.3 ,0.2 ,0.2 ,0.0, 2.7 ,3.7 ,0.2
        # print ar1,ab1,ar2,ab2,ar3,t1,t2,c

        pr = PieceWiseLinear([Linear(ar1, 1 - ar1), Linear(ar2, 1 - ar2), Linear(ab2, -ab2 * t1 + ar2 * t1 + 1 - ar2),
                              Linear(ar3, -ar3 * t2 + ab2 * t2 - ab2 * t1 + ar2 * t1 + 1 - ar2)], [0.0, 1., t1, t2])
        pb = PieceWiseLinear([Linear(ab1, 1 - ab1), Linear(ab2, 1 - ab2), Linear(ab2, 1 - ab2),
                              Linear(ab3, -ab3 * t2 + ab2 * t2 + 1 - ab2)], [0.0, 1., t1, t2])
        cons = InstancePieceWiseLinear(pr, pb, c)

        # cons.plot(num_fig=i)
        # print cons.optDetStr()

        xx = np.arange(1.0, 20.0, 0.1)
        gg = np.array([cons.fr[x] / cons.fb[x] for x in xx])
        gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
        hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])
        # plt.plot(xx,gg,"g",xx,hh,"k")
        # plt.plot(xx,np.array([cons.fr[x]/cons.fb[x] for x in xx]),"g")
        gg = np.array([np.max(gg[:j + 1]) for j in range(gg.size)])
        hh = np.array([(cons.fr[x] + c) / cons.fb[x] for x in xx])

        assert (np.abs(cons.intersection() - 1.0) < eps_num)

        t_1 = cons.calc_t1(eps)
        # print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
        assert (-eps_num < t_1 < 1.0 + eps_num)
        assert ((cons.fr[t_1] + c) <= eps_num + (1 + eps) * (cons.f_rb(0) + c + cons.fb[t_1]))
        assert ((np.abs(t_1 - 1.0) < eps_num) or np.abs(
            cons.fr[t_1] + c - (1 + eps) * (cons.f_rb(0) + c + cons.fb[t_1])) < eps_num)

        alp = cons.calc_t_N_minus_1(eps)

        assert (alp >= 1.0)
        # print alp,cons.rho_inf(),cons.fr(alp)
        flag1 = (cons.fr[alp] / cons.fb[alp]) >= eps_num + cons.cr_inf() / (1 + eps)
        # print alp,c,cons.fr[alp]
        flag2 = (c <= eps_num + eps * cons.fr[alp])
        # print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c

        # assert(flag1 or flag2)

        # beta = (b2-b1)/c
        # rho_prob_anylsis = np.exp(beta)/(np.exp(beta)-b2+b1)

        # print cons.approxOptStr(eps)[0],rho_prob_anylsis
        # ratio = cons.approxOptStr(eps)[0]/rho_prob_anylsis
        # assert ((ratio<=eps_num+(1+eps)) and (ratio>= -eps_num+1/(1+eps)))

        # check det
        # rho_det_anylsis = min(div(b2,b1),c+1,div(a1,a2))
        # t_det_anylsis = [0,1,np.inf][np.argmin([div(b2,b1),c+1,div(a1,a2)])]

        # t_det_gen,rho_det_gen = cons.optDetStr()
        # ratio = rho_det_gen/rho_det_anylsis
        # assert((ratio<=eps_num+(1+eps)) and (ratio>= -eps_num+1/(1+eps)))

        # if t_det_gen==np.inf and t_det_anylsis==np.inf:
        #    ratio=0.0
        # else:
        #    ratio = t_det_gen-t_det_anylsis

        # assert((ratio<=eps_num) and (ratio>= -eps_num))

        cr0 = cons.cr_0()
        cr0_thr = div(1 - ab1, 1 - ar1)
        assert (np.abs(cr0_thr - cr0) < eps_num or (min(cr0, cr0_thr) == np.inf))

        cri = cons.cr_inf()

        cri_thr = np.inf if ar3 > 0 else cons.fr[50] / cons.fb[50]
        # assert(np.abs(cri_thr-cri)<eps_num or (min(cri,cri_thr)==np.inf ))
        # TODO - check det

        cons.approxOptStr(eps)

    return True


def exmp1_cellular():
    # first example of determinstic solution
    # HOT - 8G=40sh +1G=6sh
    # 12G - 50sh + 1G=6sh
    num_exmp = 10
    i = 0
    while i < num_exmp:
        r0, b0 = np.random.uniform(35.0, 40.0), np.random.uniform(45.0, 70.)
        t1, t2 = np.random.uniform(2.0, 8.0), np.random.uniform(9.0, 14.)
        ar1, ab1 = 6.0, 6.0
        factor_c = np.random.uniform(1.0, 2.0)
        ending = True
        if ending:
            endT = 6.0
            pr = PieceWiseLinear(
                [Linear(0.0, r0), Linear(ar1, r0 - t1 * ar1), Linear(ar1, r0 - t1 * ar1), Linear(0.0, r0 + endT * ar1)],
                [0.0, t1, t2, t1 + endT])
            pb = PieceWiseLinear([Linear(0.0, b0), Linear(0.0, b0), Linear(ab1, b0 - t2 * ab1),
                                  Linear(0.0, b0 - t2 * ab1 + ab1 * (t1 + endT))], [0.0, t1, t2, t1 + endT])
        else:
            pr = PieceWiseLinear([Linear(0.0, r0), Linear(ar1, r0 - t1 * ar1), Linear(ar1, r0 - t1 * ar1)],
                                 [0.0, t1, t2])
            pb = PieceWiseLinear([Linear(0.0, b0), Linear(0.0, b0), Linear(ab1, b0 - t2 * ab1)], [0.0, t1, t2])

        cons = InstancePieceWiseLinear(pr, pb, (b0 - r0) * factor_c)

        intr = cons.intersection()
        if intr == np.inf:
            continue
        if not cons.valid():
            continue

        cons.plot(num_fig=i + 1)

        # cons_norm = cons.normalize()
        # cons_norm.plot(num_fig=2)

        z_opt, _ = cons.optDetStr()
        z_plot = np.minimum(z_opt, plt.xticks()[0][-2])
        z_plot = np.maximum(z_plot, 1)
        plt.plot(i + 1);
        plt.plot([z_plot, z_plot], [0, plt.yticks()[0][-1]], "-k")
        if z_opt > intr and z_opt < np.inf:
            print i, True
        else:
            pass
            print False, intr, z_opt, z_plot
        i = i + 1
        # print [str(lin) for lin in cons.fb.lins]
        # print cons.times
        # print cons_norm.optDetStr()
    # print cons_norm.middle_time(),cons_norm.intersection_time(),cons_norm.optDetStr(),cons_norm.c
    # TODO - plot the results of the DET on the graph


def exmp2_cellular():
    # second  example of DET - here we will show intersing properities of
    # instances and their strategy - consult Boaz
    # HOT - 8G=40sh +1G=6sh
    # 12G - 50sh + 1G=6sh
    # try only 2 slopes it will suppose to work
    r0, b0 = 35.0, 45.0
    t1, t2 = 6.0, 10.0
    ar1, ab1 = 6.0, 6.0
    factor_c = 1.0
    ending = True
    if ending:
        endT = 6.0
        pr = PieceWiseLinear(
            [Linear(0.0, r0), Linear(ar1, r0 - t1 * ar1), Linear(ar1, r0 - t1 * ar1), Linear(0.0, r0 + endT * ar1)],
            [0.0, t1, t2, t1 + endT])
        pb = PieceWiseLinear([Linear(0.0, b0), Linear(0.0, b0), Linear(ab1, b0 - t2 * ab1),
                              Linear(0.0, b0 - t2 * ab1 + ab1 * (t1 + endT))], [0.0, t1, t2, t1 + endT])
    else:
        pr = PieceWiseLinear([Linear(r0, 0.0), Linear(r0 - t1 * ar1, ar1), Linear(r0 - t1 * ar1, ar1)], [0.0, t1, t2])
        pb = PieceWiseLinear([Linear(b0, 0.0), Linear(b0, 0.0), Linear(b0 - t2 * ab1, ab1)], [0.0, t1, t2])
    cons = InstancePieceWiseLinear(pr, pb, (b0 - r0) * factor_c)
    cons.plot(num_fig=1)

    cons_norm = cons.normalize()
    cons_norm.plot(num_fig=2)


def exmp3_cellular():
    # third  example of TIGHT (4.1) - we need to find example of instance without tight startegy
    # better with f_r unbounded (ending is false)

    # CODE function which find TIGHT and reports when they could not find one
    # FIND what is the optimal PROB
    # ANYLSIS why we could never find such strategy - HARD
    # we might use some points as milestone and show the unextinence of tight during bounds on CR
    # we can also show imposiblity for smooth startegy

    r0, b0 = 3
    t1, t2 = 6.0, 10.0
    ar1, ab1 = 6.0, 6.0
    factor_c = 1.0
    ending = False
    if ending:
        endT = 6.0
        pr = PieceWiseLinear(
            [Linear(0.0, r0), Linear(ar1, r0 - t1 * ar1), Linear(ar1, r0 - t1 * ar1), Linear(0.0, r0 + endT * ar1)],
            [0.0, t1, t2, t1 + endT])
        pb = PieceWiseLinear([Linear(0.0, b0), Linear(0.0, b0), Linear(ab1, b0 - t2 * ab1),
                              Linear(0.0, b0 - t2 * ab1 + ab1 * (t1 + endT))], [0.0, t1, t2, t1 + endT])
    else:
        # pr=PieceWiseLinear([Linear(r0,0.0),Linear(r0-t1*ar1,ar1),Linear(r0-t1*ar1,ar1)],[0.0,t1,t2])
        # pb=PieceWiseLinear([Linear(b0,0.0),Linear(b0,0.0),Linear(b0-t2*ab1,ab1)],[0.0,t1,t2])

        pr = PieceWiseLinear([Linear(0.0, r0), Linear(ar1, r0 - t1 * ar1), Linear(ar1, r0 - t1 * ar1)], [0.0, t1, t2])
        pb = PieceWiseLinear([Linear(0.0, b0), Linear(0.0, b0), Linear(ab1, b0 - t2 * ab1)], [0.0, t1, t2])

    cons = InstancePieceWiseLinear(pr, pb, (b0 - r0) * factor_c)
    cons.plot(num_fig=1)

    cons_norm = cons.normalize()
    cons_norm.plot(num_fig=2)


def exmp4_cellular():
    global cons, cr, cr_lp
    # forth  example of TIGHT (4.2) - show instace which is bound and has tight strategy
    # which are not optimal (HARD - maybe impossible)
    # we can maybe approximate the optimal PROB with EXP curve and then show tight is not optimally
    # Q - how could we be sure this is th eonly TIGHT ?
    eps = 0.01
    print "exmp4 - tight which not optimal (bound case)"
    eps_num = 1e-10
    num_iter = 1

    for i in range(num_iter):
        print i
        ar1 = np.round(np.random.uniform(0.06, 1), 1)
        ab1 = np.round(np.random.uniform(0.0, ar1 - 0.07), 1)
        ar2 = np.round(np.random.uniform(0.3, 5.0), 1)
        ab2 = np.round(np.random.uniform(0.1, ar2 - 0.1), 1)
        c = max(ar1 - ab1, np.round(10 ** np.random.uniform(-0.3, 1.0), 1))
        ar3 = 0.0  # ab2 #np.random.choice([0.0,ab2])
        t1 = np.round(np.random.uniform(1.2, 5.0), 1)
        # TODO - make sure t2 is before intersection
        t2 = np.round(np.random.uniform(t1 + 0.1, 5.0), 1)
        ab3 = 0.0  # ab2 #0.0
        ar1, ab1, ar2, ab2, ar3, ab3, t1, t2, c = 0.9, 0.0, 0.9, 0.0, 0.0, 0.0, 1.5, 2.0, 2.0
        ar1, ab1, ar2, ab2, ar3, ab3, ar4, t1, t2, c = 0.9, 0.0, 0.9, 0.0, 0.0, 0.0, 5.0, 1.1, 3.0, 1.0
        # print ar1,ab1,ar2,ab2,ar3,ab3,t1,t2,c

        pr = PieceWiseLinear([Linear(ar1, 1 - ar1), Linear(ar2, 1 - ar2), Linear(ar3, -ar3 * t1 + ar2 * t1 + 1 - ar2),
                              Linear(ar4, -ar4 * t2 + ar3 * t2 - ar3 * t1 + ar2 * t1 + 1 - ar2)], [0.0, 1., t1, t2])
        pb = PieceWiseLinear([Linear(ab1, 1 - ab1), Linear(ab2, 1 - ab2), Linear(ab2, 1 - ab2),
                              Linear(ab3, -ab3 * t2 + ab2 * t2 + 1 - ab2)], [0.0, 1., t1, t2])

        # example
        ar1, ab1, ar2, ab2, ar3, ab3, t1, t2 = 0.9, 0.0, 4.0, 0.0, 0.0, 0.0, 0.4, 0.55
        # m(x-x0)+y0=n => n = -mx0+y0
        pr = PieceWiseLinear([Linear(ar1, 1 - ar1), Linear(ar2, -ar2 * t1 + ar1 * t1 + 1 - ar1),
                              Linear(ar3, -ar3 * t2 + ar2 * t2 - ar2 * t1 + ar1 * t1 + 1 - ar1)], [0.0, t1, t2])
        pb = PieceWiseLinear([Linear(ab1, 1 - ab1), Linear(ab2, -ab2 * t1 + ab1 * t1 + 1 - ab1),
                              Linear(ab3, -ab3 * t2 + ab2 * t2 - ab2 * t1 + ab1 * t1 + 1 - ab1)], [0.0, t1, t2])

        cons = InstancePieceWiseLinear(pr, pb, c)
        cons = create_instance(0.1, 1.0, [0.9, 0.9, 0.0, 0.3, 0], [0.0, 0.0, 0.0, 0.2], [0.0, 1.0, 1.3, 2.2], 1.0)
        cons.plot(num_fig=i)
        cons_norm = cons.normalize()
        cons_norm.plot(num_fig=i)
        cr, q, T, cr_lp = cons_norm.approxOptStr(eps)

        plt.plot(T, cr + cr_lp, "k*-")
        print len(T), np.cumsum(q).shape
        plt.plot(T[:-1], np.cumsum(q)[:-1], "g*-")
        print "CR", cr
        # print T
    return True


def exmp5_cellular():
    # fifth  example of appendix - show parameters of instance for APPROX
    r0, b0 = 35.0, 45.0
    t1, t2 = 6.0, 10.0
    ar1, ab1 = 6.0, 6.0
    factor_c = 1.0
    ending = False
    if ending:
        endT = 6.0
        pr = PieceWiseLinear(
            [Linear(0.0, r0), Linear(ar1, r0 - t1 * ar1), Linear(ar1, r0 - t1 * ar1), Linear(0.0, r0 + endT * ar1)],
            [0.0, t1, t2, t1 + endT])
        pb = PieceWiseLinear([Linear(0.0, b0), Linear(0.0, b0), Linear(ab1, b0 - t2 * ab1),
                              Linear(0.0, b0 - t2 * ab1 + ab1 * (t1 + endT))], [0.0, t1, t2, t1 + endT])
    else:
        # pr=PieceWiseLinear([Linear(r0,0.0),Linear(r0-t1*ar1,ar1),Linear(r0-t1*ar1,ar1)],[0.0,t1,t2])
        # pb=PieceWiseLinear([Linear(b0,0.0),Linear(b0,0.0),Linear(b0-t2*ab1,ab1)],[0.0,t1,t2])

        pr = PieceWiseLinear([Linear(0.0, r0), Linear(ar1, r0 - t1 * ar1), Linear(ar1, r0 - t1 * ar1)], [0.0, t1, t2])
        pb = PieceWiseLinear([Linear(0.0, b0), Linear(0.0, b0), Linear(ab1, b0 - t2 * ab1)], [0.0, t1, t2])

    cons = InstancePieceWiseLinear(pr, pb, (b0 - r0) * factor_c)
    cons.plot(num_fig=1)

    cons_norm = cons.normalize()
    cons_norm.plot(num_fig=2)


def plot_results_det():
    global dict_res
    file_name = "dict_res_det_c_td_br0.pkl"
    dict_res = pickle.load(open(file_name, "r"))
    markersize = 18
    fontsize = 20

    dict_marker = {np.inf: "o", 1.0: "d", 0.0: "*"}
    fig = plt.figure(2)
    fig.clf()
    ax = Axes3D(fig)
    sc = ax.scatter(xs=[k[0] for k in dict_res.keys() if dict_res[k][1] == np.inf],
                    ys=[k[1] for k in dict_res.keys() if dict_res[k][1] == np.inf],
                    zs=[k[2] for k in dict_res.keys() if dict_res[k][1] == np.inf], marker="o",
                    c=[(dict_res[k][0]) for k in dict_res.keys() if dict_res[k][1] == np.inf], cmap='plasma', s=60)
    sc = ax.scatter(xs=[k[0] for k in dict_res.keys() if dict_res[k][1] == 1.0],
                    ys=[k[1] for k in dict_res.keys() if dict_res[k][1] == 1.0],
                    zs=[k[2] for k in dict_res.keys() if dict_res[k][1] == 1.0], marker="d",
                    c=[(dict_res[k][0]) for k in dict_res.keys() if dict_res[k][1] == 1.0], cmap='plasma', s=60)
    sc = ax.scatter(xs=[k[0] for k in dict_res.keys() if dict_res[k][1] == 0.0],
                    ys=[k[1] for k in dict_res.keys() if dict_res[k][1] == 0.0],
                    zs=[k[2] for k in dict_res.keys() if dict_res[k][1] == 0.0], marker="*",
                    c=[(dict_res[k][0]) for k in dict_res.keys() if dict_res[k][1] == 0.0], cmap='plasma', s=60)

    ax.set_xlabel("Switch cost", fontsize=fontsize)
    ax.set_ylabel("Time diffrenece", fontsize=fontsize)
    ax.set_zlabel("Initial rent price", fontsize=fontsize)
    ax.set_title("Determinstic strategies properties", fontsize=25)
    simArtists = [plt.Line2D((4, 2, 40), (4, 2.001, 40), color="k", marker=m, linestyle='-', markersize=15) for m in
                  ["*", "d", "o"]]

    ax.legend(simArtists,
              ["best startegy start with buy", "best startegy switch at middle time", "best strategy never buy"])
    # ax.legend(sc)
    fig.colorbar(sc)
    return
    for k in dict_res.keys():
        (c, td) = k
        rho_det_gen, t_det_gen = dict_res[k]
        if t_det_gen == 0:
            plt.plot(c, td, "b*", markersize=markersize)
        elif t_det_gen < np.inf:
            plt.plot(c, td, "gx", markersize=markersize)
        elif t_det_gen == np.inf:
            plt.plot(c, td, "ro", markersize=markersize)
        else:
            print "no det choice"
    plt.figure(2)
    plt.clf()
    plt.hold(True)
    plt.xlabel("Switch cost", fontsize=fontsize)
    plt.ylabel("Time diffrenece", fontsize=fontsize)
    colormap = plt.cm.gist_ncar
    # you can add cmap

    sc = plt.scatter([k[0] for k in dict_res.keys()], [k[1] for k in dict_res.keys()], s=100,
                     c=[(dict_res[k][0]) for k in dict_res.keys()], cmap='plasma')
    plt.colorbar()

    simArtists = [ax.Line3D((4, 2, 40), (4, 2.001, 40), color=c, marker=m, linestyle='-') for m in ["*", "d", "0"]]

    # ax.legend(simArtists,["best startegy start with buy","best startegy switch at middle time","best strategy never buy"])
    ax.legend(sc)

    plt.show()
    # 3D - https://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html
    # Axes3D.scatter(xs, ys, zs=0, zdir='z', s=20, c=None, depthshade=True, *args, **kwargs)

    return


def plot_results_prob():
    global dict_res
    file_name = "dict_res_prob_c_td_ar1_br0.pkl"
    dict_res = pickle.load(open(file_name, "r"))
    markersize = 18
    fontsize = 20

    ind_feature = 1  # 0-cr,1-q0,2-q_inf
    strs = ["Competetive ratio", "Start at buy probability", "Never switch probability"]
    for ind_feature in range(3):
        fig = plt.figure(3 + ind_feature)
        fig.clf()
        ax = Axes3D(fig)
        sc = ax.scatter([k[0] for k in dict_res.keys()], [k[1] for k in dict_res.keys()],
                        zs=[k[3] for k in dict_res.keys()], marker="o",
                        c=[(dict_res[k][ind_feature]) for k in dict_res.keys()], cmap='plasma', s=60)

        ax.set_xlabel("Switch cost", fontsize=fontsize)
        ax.set_ylabel("Time diffrenece", fontsize=fontsize)
        ax.set_zlabel("Initial rent price", fontsize=fontsize)
        ax.set_title("Probabilities strategies properties (" + strs[ind_feature] + ")", fontsize=30)
        # simArtists = [plt.Line2D((4,2,40),(4,2.001,40), color="k", marker=m, linestyle='-',markersize=15) for m in ["*","d","o"]]

        # ax.legend(simArtists,["best startegy start with buy","best startegy switch at middle time","best strategy never buy"])
        # ax.legend(sc)
        fig.colorbar(sc)
    return


def plot_results_det_2D():
    global dict_res_slice, dict_res
    # we need to create  surface by using meshgrid and  building Z (competetive ratio and partition)
    file_name = "dict_res_uniified_br0_c_td.pkl"
    dict_res = pickle.load(open(file_name, "r"))
    # dict_res[br0,c,td] = (t_det_gen, rho_det_gen,rho_prob,  q, T)
    # genr = ((br0, c, td) for br0 in np.arange(10.0, 38.0, 2) for c in np.arange((bb0 - br0) + 0.2, bb0 - 0.1, 2) for td
    #     in np.arange((bb0 - br0) / ar1 + 0.1, 10.0, 0.5))
    print np.max([dict_res[k][1] for k in dict_res.keys()])
    # print [sorted(list(set([k[i] for k in dict_res.keys()]))) for i in range(3)]
    targets = [np.arange(10, 39, 4), np.arange(4.2, 39, 4), np.arange(0.9, 10, 0.2)]
    starts = [10, 4.2, 0.9]
    ends = [39, 39, 10]
    tols_max = [2, 2, 0.1]
    tols_r = [4, 4, 1]

    ind_target = 0
    parm_strs = ["Initial rent cost", "Switch cost", "Time difference"]
    normalize = matplotlib.colors.Normalize(vmin=1, vmax=1.8)

    for ind, target in enumerate(np.arange(starts[ind_target], ends[ind_target], tols_r[ind_target])):
        rest = list(range(3))
        rest.remove(ind_target)
        dict_res_slice = {(k[rest[0]], k[rest[1]]): dict_res[k] for k in dict_res.keys() if
                          abs(k[ind_target] - target) < 0.005}
        xs = set([k[0] for k in dict_res_slice.keys()])
        ys = list(set([k[1] for k in dict_res_slice.keys()]))
        # (y - eps_num) > ((40 - x) / 5)

        if True:
            for x in xs:
                for y in ys:
                    # x=br0, y=td
                    xys = np.array(dict_res_slice.keys())
                    ds = np.min(np.sum((np.array([x, y]) - xys) ** 2, axis=1))
                    if ds > 0.1:
                        print x, y
                        if y <= ((40 - x) / 5):
                            print x, y
                            dict_res_slice[x, y] = [np.inf, 1, None, None, None]

        markersize = 18
        fontsize = 20

        dict_marker = {np.inf: "o", 1.0: "d", 0.0: "*"}
        plt.figure(10 + ind)
        plt.clf()
        for km in dict_marker.keys():
            plt.scatter(x=[k[0] for k in dict_res_slice.keys() if dict_res_slice[k][0] == km],
                        y=[k[1] for k in dict_res_slice.keys() if dict_res_slice[k][0] == km], marker=dict_marker[km],
                        c=[(dict_res_slice[k][1]) for k in dict_res_slice.keys() if dict_res_slice[k][0] == km],
                        cmap='plasma', s=60, norm=normalize, linewidth=0)
        plt.colorbar()
        plt.title(parm_strs[ind_target] + "=" + str(target), fontsize=fontsize)
        plt.xlabel(parm_strs[rest[0]], fontsize=fontsize)
        plt.ylabel(parm_strs[rest[1]], fontsize=fontsize)

        ax = plt.gca()

        # simArtists = [plt.Line2D((4, 2, 40), (4, 2.001, 40), color="k", marker=m, linestyle='-') for m in ["*", "d", "o"]]
        simArtists = [plt.Line2D([0], [0], color="k", marker=m, linestyle='-') for m in
                      ["*", "d", "o"]]

        ax.legend(simArtists,
                  ["best startegy start with buy", "best startegy switch at middle time", "best strategy never buy"])
    plt.show()
    return


def plot_results_det_1D():
    global dict_res_slice, dict_res
    # first we need to slice
    file_name = "dict_res_uniified_br0_c_td.pkl"
    dict_res = pickle.load(open(file_name, "r"))
    # dict_res[br0,c,td] = (t_det_gen, rho_det_gen,rho_prob,  q, T)
    # genr = ((br0, c, td) for br0 in np.arange(10.0, 38.0, 2) for c in np.arange((bb0 - br0) + 0.2, bb0 - 0.1, 2) for td
    #     in np.arange((bb0 - br0) / ar1 + 0.1, 10.0, 0.5))
    # print [sorted(list(set([k[i] for k in dict_res.keys()]))) for i in range(3)]
    targets = [np.arange(10, 39, 4), np.arange(4.2, 39, 4), np.arange(0.9, 10, 0.2)]
    starts = [10, 4.2, 0.9]
    ends = [39, 39, 10]
    tols_max = [2, 2, 0.1]
    tols_r = [4, 4, 1]

    parm_strs = ["Initial rent cost", "Switch cost", "Time difference"]

    target_c = 34.0
    target_td = 12.0
    dict_res_slice = {k[0]: dict_res[k] for k in dict_res.keys() if
                      (abs(k[1] - target_c) < 0.005) and (abs(k[2] - target_td) < 0.005)}

    eps = 0.05
    eps_num = 1e-10
    dict_res_slice = {}
    bb0 = 50.0
    ar1 = 5.0
    td = target_td
    c = target_c
    genr = ((br0) for br0 in np.arange(10.0, 48.0, 0.2) if (c > (bb0 - br0 - eps_num)))
    for ind, (br0) in enumerate(genr):
        if ((td - eps_num) < ((bb0 - br0) / ar1)):
            dict_res_slice[br0] = [np.inf, 1, None, None, None]
            continue

        assert (td - eps_num) > ((bb0 - br0) / ar1), "BAD TD" + ",".join(str(si) for si in [br0, c, td])

        assert (c > (bb0 - br0 - eps_num)), "BAD C" + ",".join(str(si) for si in [br0, c, td])

        ts = [0, 10, 10 + td]
        ar_s = [0, ar1, 0]
        ab_s = [0, 0, 0]

        cons_org = create_instance(br0, bb0, ar_s, ab_s, ts, c)
        cons = cons_org.normalize()
        assert (np.max(np.abs(np.array(
            [cons.fr.lins[i](cons.times[i + 1]) - cons.fr.lins[i + 1](cons.times[i + 1]) for i in
             range(len(cons.times) - 1)]))) < eps), "BAD fr" + ",".join(str(si) for si in [br0, c, td])
        assert (np.max(np.abs(np.array(
            [cons.fb.lins[i](cons.times[i + 1]) - cons.fb.lins[i + 1](cons.times[i + 1]) for i in
             range(len(cons.times) - 1)]))) < eps), "BAD fb" + ",".join(str(si) for si in [br0, c, td])

        assert (np.abs(cons.intersection() - 1.0) < eps_num), "BAD intersection" + ",".join(
            str(si) for si in [br0, c, td])

        t_1 = cons.calc_t1(eps)
        # print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
        assert (-eps_num < t_1 < 1.0 + eps_num), "BAD t1_1 " + ",".join(str(si) for si in [br0, c, td])
        assert ((cons.fr[t_1] + cons.c) <= eps_num + (1 + eps) * (
                cons.f_rb(0) + cons.c + cons.fb[t_1])), "BAD t1_2 " + ",".join(str(si) for si in [br0, c, td])
        tol = 0.5
        assert ((np.abs(t_1 - 1.0) < eps_num) or np.abs(cons.fr[t_1] + cons.c - (1 + eps) * (
                cons.f_rb(0) + cons.c + cons.fb[t_1])) < eps_num), "BAD t1_3 " + ",".join(
            str(si) for si in [br0, c, td]) + "|" + str(t_1) + "|" + str(
            np.abs(cons.fr[t_1] + c - (1 + eps) * (cons.f_rb(0) + c + cons.fb[t_1]))) + "|" + str(
            cons.fr[1.0] + c) + "|" + str((1 + eps) * (cons.f_rb(0) + c + cons.fb[1.0]))
        # print "no assert of t_1"

        alp = cons.calc_t_N_minus_1(eps)

        assert (alp >= 1.0), "BAD alp" + ",".join(str(si) for si in [br0, c, td])
        # print alp,cons.rho_inf(),cons.fr(alp)
        flag1 = (cons.fr[alp] / cons.fb[alp]) >= -eps_num + cons.cr_inf() / (1 + eps)
        # print alp,c,cons.fr[alp]
        flag2 = (cons.c <= eps_num + eps * cons.fr[alp])
        # print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c

        assert (flag1 or flag2), "BAD flags" + ",".join(str(si) for si in [br0, c, td])

        cr0 = cons.cr_0()
        cr0_thr = div(bb0, br0)
        assert (np.abs(cr0_thr - cr0) < eps_num or (min(cr0, cr0_thr) == np.inf)), "BAD cr0" + ",".join(
            str(si) for si in [br0, c, td])

        cri = cons.cr_inf()

        cri_thr = div(ar_s[-1], ab_s[-1]) if ar_s[-1] > 0 else cons.fr[20] / cons.fb[20]
        # print cri_thr,cri
        assert (np.abs(cri_thr - cri) < eps_num or (min(cri, cri_thr) == np.inf)), "BAD cri" + ",".join(
            str(si) for si in [br0, c, td])

        t_det_gen, rho_det_gen = cons.optDetStr()
        dict_res_slice[br0] = (t_det_gen, rho_det_gen, None, None, None)

    markersize = 18
    fontsize = 20
    normalize = matplotlib.colors.Normalize(vmin=1, vmax=1.8)
    dict_marker = {np.inf: "o", 1.0: "d", 0.0: "*"}
    plt.figure(10)
    plt.clf()
    dict_size = {"o": 100, "d": 100, "*": 130}
    for km in dict_marker.keys():
        pass
        """plt.scatter(x=[k for k in dict_res_slice.keys() if dict_res_slice[k][0] == km],
                        y=[dict_res_slice[k][1] for k in dict_res_slice.keys() if dict_res_slice[k][0] == km], marker=dict_marker[km],
                        c=[(dict_res_slice[k][1]) for k in dict_res_slice.keys() if dict_res_slice[k][0] == km],
                        cmap='plasma', s=dict_size[dict_marker[km]],norm=normalize,linewidth=0)"""
    xs = sorted(dict_res_slice.keys())
    plt.plot(xs, [dict_res_slice[x][1] for x in xs], "r-")
    print("xs: ",xs)
    print("ys: ",[dict_res_slice[x][1] for x in xs])
    print [
        "change from {0} to {1} at {2} ({3})".format(dict_res_slice[xs[i]][0], dict_res_slice[xs[i + 1]][0], xs[i], i)
        for i in range(len(xs) - 1) if dict_res_slice[xs[i]][0] != dict_res_slice[xs[i + 1]][0]]
    fs = 12
    plt.plot([24, 24], [0, 1.68], "k--")
    plt.text(20, 1.2, "Never switch", {"fontsize": fs})
    plt.plot([29.6, 29.6], [0, 1.68], "k--")
    plt.text(24.1, 1.4, "Switch at intersection time", {"fontsize": fs})
    plt.text(32, 1.0, "Start at buy", {"fontsize": fs})

    # plt.colorbar()
    plt.title(parm_strs[2] + "=" + str(target_td) + ", " + parm_strs[1] + "=" + str(target_c), fontsize=fontsize)
    plt.xlabel(parm_strs[0], fontsize=fontsize)
    plt.ylabel("Competitive ratio", fontsize=fontsize)
    ks = sorted(dict_res_slice.keys())
    ys = [dict_res_slice[k][1] for k in ks]
    # plt.plot(ks,[dict_res_slice[k][1] for k in ks],"k-")
    plt.axis([15, 48, 1, 1.8])
    ax = plt.gca()

    # simArtists = [plt.Line2D((4, 2, 40), (4, 2.001, 40), color="k", marker=m, linestyle='-') for m in ["*", "d", "o"]]
    simArtists = [plt.Line2D([0], [0], color="k", marker=m, linestyle='-') for m in
                  ["*", "d", "o"]]

    # ax.legend(simArtists,["best startegy start with buy","best startegy switch at middle time","best strategy never buy"])
    # plt.legend(simArtists,["best startegy start with buy","best startegy switch at middle time","best strategy never buy"],loc=2)
    plt.show()
    return

    ax.set_xlabel("Switch cost", fontsize=fontsize)
    ax.set_ylabel("Time diffrenece", fontsize=fontsize)
    ax.set_zlabel("Initial rent price", fontsize=fontsize)
    ax.set_title("Determinstic strategies properties", fontsize=25)
    simArtists = [plt.Line2D((4, 2, 40), (4, 2.001, 40), color="k", marker=m, linestyle='-', markersize=15) for m in
                  ["*", "d", "o"]]

    ax.legend(simArtists,
              ["best startegy start with buy", "best startegy switch at middle time", "best strategy never buy"])
    # ax.legend(sc)
    fig.colorbar(sc)
    return


def example_insatnce():
    c = 24.0
    td = 5.0
    br0 = 40
    eps = 0.05
    eps_num = 1e-10
    dict_res_slice = {}
    bb0 = 50.0
    ar1 = 5.0
    ind_target = 1
    ts = [0, 10, 10 + td]
    ar_s = [0, ar1, 0]
    ab_s = [0, 0, 0]

    cons_org = create_instance(br0, bb0, ar_s, ab_s, ts, c)
    cons = cons_org.normalize()
    plt = cons.plot()
    plt.hold(True)
    plt.axis([0.0, (ts[1] + td + 2) / cons.norm_time, 0.0, cons.fr(ts[-1]) + 0.2])
    plt.plot([(ts[1] + td) / cons.norm_time] * 2, [0.0, cons.fr(ts[-1]) + 0.2], 'k--')
    plt.plot([ts[1] / cons.norm_time] * 2, [0.0, cons.fr(ts[-1]) + 0.2], 'k--')

    #plt.annotate('Initial rent cost', xy=(0.3, cons.fr(0)), xytext=(0.3, 0.05),
     #            arrowprops=dict(facecolor='black', shrink=1), fontsize=22
     #            )
    plt.annotate("", xy=((ts[1] + td) / cons.norm_time, 0.7),
                 xytext=((ts[1] ) / cons.norm_time, 0.7),
                 arrowprops=dict(arrowstyle="<->",linewidth=2), fontsize=25
                 )
    plt.annotate("Time difference", xy=((ts[1] + td) / cons.norm_time, 0.715),
                 xytext=((ts[1] +1.0) / cons.norm_time, 0.72),
                  fontsize=28
                 )
    plt.annotate("Rent", xy=(0.5, 1.02),
                 xytext=(0.5, 1.02),
                  fontsize=28
                 )

    plt.annotate("Buy", xy=(0.5, 0.72),
                 xytext=(0.5, 0.72),
                  fontsize=28
                 )
    plt.show()
    return

    # choose params
    # plot
    # arrow initial rent cost
    # arrow time diffrence (+2 dasheed lines)
    # axis time,cost


def plot_results_det_2D_surf():
    global dict_res_slice, dict_res
    # first we need to slice

    # dict_res[br0,c,td] = (t_det_gen, rho_det_gen,rho_prob,  q, T)
    # genr = ((br0, c, td) for br0 in np.arange(10.0, 38.0, 2) for c in np.arange((bb0 - br0) + 0.2, bb0 - 0.1, 2) for td
    #     in np.arange((bb0 - br0) / ar1 + 0.1, 10.0, 0.5))
    # print [sorted(list(set([k[i] for k in dict_res.keys()]))) for i in range(3)]
    targets = [np.arange(10, 39, 4), np.arange(4.2, 39, 4), np.arange(0.9, 10, 0.2)]
    starts = [10, 4.2, 0.9]
    ends = [39, 39, 10]
    tols_max = [2, 2, 0.1]
    tols_r = [4, 4, 1]

    parm_strs = ["Initial rent cost", "Switch cost", "Time difference"]

    target_c = 24.0
    target_td = 12.0
    target_br0 = 40
    eps = 0.05
    eps_num = 1e-10
    dict_res_slice = {}
    bb0_init = 50.0
    ar1 = 5.0
    ind_target = 1
    tol_b, tol_t = 0.1, 0.01
    if ind_target == 0:
        br0 = target_br0
        C1 = np.arange(2, 50, 1)
        T1 = np.arange(0.05, 20, 0.1)
        C, T = np.meshgrid(C1, T1)
        Z = np.zeros_like(C)
    else:
        c = target_c
        B1 = np.arange(2, 50, tol_b)
        T1 = np.arange(0.05, 20, tol_t)
        B, T = np.meshgrid(B1, T1)
        Z = np.zeros_like(B)
        R = np.zeros_like(B)
    # TODO - no itersection
    tic = time.time()
    for ind_x in range(Z.shape[0]):
        for ind_y in range(Z.shape[1]):
            change_c = False
            if ind_target == 1:
                br0, td = B[ind_x, ind_y], T[ind_x, ind_y]
            else:
                c, td = C[ind_x, ind_y], T[ind_x, ind_y]
            if ((td - eps_num) < ((bb0_init - br0) / ar1)):
                # trival case - rent always less then buy
                Z[ind_x, ind_y] = 1.0
                R[ind_x, ind_y] = np.inf
                continue
            assert (td - eps_num) > ((bb0_init - br0) / ar1), "BAD TD" + ",".join(str(si) for si in [br0, c, td])

            if (c <= (bb0_init - br0 - eps_num)):
                # print "bad c", c, br0, td,"before"
                bb0 = br0 + c
                change_c = True
            else:
                change_c = False
                bb0 = bb0_init

            if (br0 > bb0):
                # trival case - rent always higher then buy
                Z[ind_x, ind_y] = 1.0
                R[ind_x, ind_y] = 0.0
                continue

            assert (c > (bb0 - br0 - eps_num)), "BAD C" + ",".join(str(si) for si in [br0, c, td])

            ts = [0, 10, 10 + td]
            ar_s = [0, ar1, 0]
            ab_s = [0, 0, 0]

            cons_org = create_instance(br0, bb0, ar_s, ab_s, ts, c)
            try:
                cons = cons_org.normalize()
                if cons is None:
                    print "None", br0, c, td, bb0
                    # plt = cons_org.plot()
                    # plt.show()
                    return
            except:
                print "no intersection", br0, bb0, ar_s, ab_s, ts, c, td
                # plt=cons_org.plot()
                # plt.show()
                return

            assert (np.max(np.abs(np.array(
                [cons.fr.lins[i](cons.times[i + 1]) - cons.fr.lins[i + 1](cons.times[i + 1]) for i in
                 range(len(cons.times) - 1)]))) < eps), "BAD fr" + ",".join(str(si) for si in [br0, c, td])
            assert (np.max(np.abs(np.array(
                [cons.fb.lins[i](cons.times[i + 1]) - cons.fb.lins[i + 1](cons.times[i + 1]) for i in
                 range(len(cons.times) - 1)]))) < eps), "BAD fb" + ",".join(str(si) for si in [br0, c, td])

            assert (np.abs(cons.intersection() - 1.0) < eps_num), "BAD intersection" + ",".join(
                str(si) for si in [br0, c, td])

            t_1 = cons.calc_t1(eps)
            # print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
            assert (-eps_num < t_1 < 1.0 + eps_num), "BAD t1_1 " + ",".join(str(si) for si in [br0, c, td])
            assert ((cons.fr[t_1] + cons.c) <= eps_num + (1 + eps) * (
                    cons.f_rb(0) + cons.c + cons.fb[t_1])), "BAD t1_2 " + ",".join(str(si) for si in [br0, c, td])
            tol = 0.5
            assert ((np.abs(t_1 - 1.0) < eps_num) or np.abs(cons.fr[t_1] + cons.c - (1 + eps) * (
                    cons.f_rb(0) + cons.c + cons.fb[t_1])) < eps_num), "BAD t1_3 " + ",".join(
                str(si) for si in [br0, c, td]) + "|" + str(t_1) + "|" + str(
                np.abs(cons.fr[t_1] + c - (1 + eps) * (cons.f_rb(0) + c + cons.fb[t_1]))) + "|" + str(
                cons.fr[1.0] + c) + "|" + str((1 + eps) * (cons.f_rb(0) + c + cons.fb[1.0]))
            # print "no assert of t_1"

            alp = cons.calc_t_N_minus_1(eps)

            assert (alp >= 1.0), "BAD alp" + ",".join(str(si) for si in [br0, c, td])
            # print alp,cons.rho_inf(),cons.fr(alp)
            flag1 = (cons.fr[alp] / cons.fb[alp]) >= -eps_num + cons.cr_inf() / (1 + eps)
            # print alp,c,cons.fr[alp]
            flag2 = (cons.c <= eps_num + eps * cons.fr[alp])
            # print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c

            assert (flag1 or flag2), "BAD flags" + ",".join(str(si) for si in [br0, c, td])

            cr0 = cons.cr_0()
            cr0_thr = div(bb0, br0)
            assert (np.abs(cr0_thr - cr0) < eps_num or (min(cr0, cr0_thr) == np.inf)), "BAD cr0" + ",".join(
                str(si) for si in [br0, c, td])

            cri = cons.cr_inf()

            cri_thr = div(ar_s[-1], ab_s[-1]) if ar_s[-1] > 0 else cons.fr[20] / cons.fb[20]
            # print cri_thr,cri
            assert (np.abs(cri_thr - cri) < eps_num or (min(cri, cri_thr) == np.inf)), "BAD cri" + ",".join(
                str(si) for si in [br0, c, td])

            t_det_gen, rho_det_gen = cons.optDetStr()
            Z[ind_x, ind_y] = rho_det_gen
            R[ind_x, ind_y] = t_det_gen

    print "time elpassed ", time.time() - tic
    fig = plt.figure(1)
    ax = fig.gca(projection='3d')
    ax = plt.gca()
    R[R == np.inf] = 2 * np.max(R[R != np.inf])
    target_plot = Z
    if ind_target == 0:
        pickle.dump((C, T, Z), open("BTZ.pkl", "wb"))
        #numpy.savetxt("foo.csv", a, delimiter=",")
        surf = ax.plot_surface(C, T, target_plot, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
    else:
        pickle.dump((B, T, Z, R), open("BTZ_ind_target1.pkl", "wb"))

        surf = ax.plot_surface(B, T, target_plot, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        np.savetxt("x_figure4.csv",B, delimiter=",")
        np.savetxt("y_figure4.csv", T, delimiter=",")
        np.savetxt("z_figure4.csv", Z, delimiter=",")


    fontsize = 25
    if ind_target == 1:
        plt.title(parm_strs[1] + "=" + str(target_c), fontsize=fontsize)
        plt.xlabel(parm_strs[0], fontsize=fontsize)
        plt.ylabel(parm_strs[2], fontsize=fontsize)
    else:
        plt.title(parm_strs[0] + "=" + str(target_br0), fontsize=fontsize)
        plt.xlabel(parm_strs[1], fontsize=fontsize)
        plt.ylabel(parm_strs[2], fontsize=fontsize)

    fig.colorbar(surf)

    fig2 = plt.figure(2)
    ax = fig2.gca(projection='3d')
    ax = plt.gca()

    target_plot = R
    if ind_target == 0:
        surf2 = ax.plot_surface(C, T, target_plot, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
    else:
        surf2 = ax.plot_surface(B, T, target_plot, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        np.savetxt("x_figure5.csv", B, delimiter=",")
        np.savetxt("y_figure5.csv", T, delimiter=",")
        np.savetxt("z_figure5.csv", R, delimiter=",")

    fontsize = 25
    if ind_target == 1:
        plt.title(parm_strs[1] + "=" + str(target_c), fontsize=fontsize)
        plt.xlabel(parm_strs[0], fontsize=fontsize)
        plt.ylabel(parm_strs[2], fontsize=fontsize)
    else:
        plt.title(parm_strs[0] + "=" + str(target_br0), fontsize=fontsize)
        plt.xlabel(parm_strs[1], fontsize=fontsize)
        plt.ylabel(parm_strs[2], fontsize=fontsize)

    fig.colorbar(surf2)


    plt.show()
    return


def plot_results_prob_2D_surf():
    # TODO (BPS) - add q_0,q_inf,q_middle
    # TODO (BPS) - make sure there is only 1 colorbar
    global dict_res_slice, dict_res
    # first we need to slice

    # dict_res[br0,c,td] = (t_det_gen, rho_det_gen,rho_prob,  q, T)
    # genr = ((br0, c, td) for br0 in np.arange(10.0, 38.0, 2) for c in np.arange((bb0 - br0) + 0.2, bb0 - 0.1, 2) for td
    #     in np.arange((bb0 - br0) / ar1 + 0.1, 10.0, 0.5))
    # print [sorted(list(set([k[i] for k in dict_res.keys()]))) for i in range(3)]
    targets = [np.arange(10, 39, 4), np.arange(4.2, 39, 4), np.arange(0.9, 10, 0.2)]
    starts = [10, 4.2, 0.9]
    ends = [39, 39, 10]
    tols_max = [2, 2, 0.1]
    tols_r = [4, 4, 1]

    parm_strs = ["Initial rent cost", "Switch cost", "Time difference"]

    target_c = 24.0
    target_td = 12.0
    target_br0 = 40
    eps = 0.05
    eps_num = 1e-10
    dict_res_slice = {}
    bb0_init = 50.0
    ar1 = 5.0
    ind_target = 1
    cnt = 0
    if ind_target == 0:
        br0 = target_br0
        C1 = np.arange(2, 50, 1)
        T1 = np.arange(0.05, 20, 0.1)
        C, T = np.meshgrid(C1, T1)
        Z = np.zeros_like(C)
    else:
        c = target_c
        B1 = np.arange(2, 50, 0.1)
        T1 = np.arange(0.05, 20, 0.1)
        B, T = np.meshgrid(B1, T1)
        Z = np.zeros_like(B)
        R0 = np.zeros_like(B)
        Rm = np.zeros_like(B)
        Ri = np.zeros_like(B)
    # TODO - no itersection
    tic = time.time()
    for ind_x in range(Z.shape[0]):
        for ind_y in range(Z.shape[1]):
            change_c = False
            if ind_target == 1:
                # print ind_x,ind_y
                br0, td = B[ind_x, ind_y], T[ind_x, ind_y]
            else:
                c, td = C[ind_x, ind_y], T[ind_x, ind_y]
            if ((td - eps_num) < ((bb0_init - br0) / ar1)):
                Z[ind_x, ind_y] = 1.0
                Ri[ind_x, ind_y] = 1
                continue
            assert (td - eps_num) > ((bb0_init - br0) / ar1), "BAD TD" + ",".join(str(si) for si in [br0, c, td])

            if (c <= (bb0_init - br0 - eps_num)):
                # print "bad c", c, br0, td,"before"
                bb0 = br0 + c
                change_c = True
            else:
                change_c = False
                bb0 = bb0_init

            if (br0 > bb0):
                Z[ind_x, ind_y] = 1.0
                R0[ind_x, ind_y] = 1.0
                continue
            cnt += 1
            assert (c > (bb0 - br0 - eps_num)), "BAD C" + ",".join(str(si) for si in [br0, c, td])

            ts = [0, 10, 10 + td]
            ar_s = [0, ar1, 0]
            ab_s = [0, 0, 0]

            cons_org = create_instance(br0, bb0, ar_s, ab_s, ts, c)
            try:
                cons = cons_org.normalize()
                if cons is None:
                    print "None", br0, c, td, bb0
                    # plt = cons_org.plot()
                    # plt.show()
                    return
            except:
                print "no intersection", br0, bb0, ar_s, ab_s, ts, c, td
                # plt=cons_org.plot()
                # plt.show()
                return

            assert (np.max(np.abs(np.array(
                [cons.fr.lins[i](cons.times[i + 1]) - cons.fr.lins[i + 1](cons.times[i + 1]) for i in
                 range(len(cons.times) - 1)]))) < eps), "BAD fr" + ",".join(str(si) for si in [br0, c, td])
            assert (np.max(np.abs(np.array(
                [cons.fb.lins[i](cons.times[i + 1]) - cons.fb.lins[i + 1](cons.times[i + 1]) for i in
                 range(len(cons.times) - 1)]))) < eps), "BAD fb" + ",".join(str(si) for si in [br0, c, td])

            assert (np.abs(cons.intersection() - 1.0) < eps_num), "BAD intersection" + ",".join(
                str(si) for si in [br0, c, td])

            t_1 = cons.calc_t1(eps)
            # print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
            assert (-eps_num < t_1 < 1.0 + eps_num), "BAD t1_1 " + ",".join(str(si) for si in [br0, c, td])
            assert ((cons.fr[t_1] + cons.c) <= eps_num + (1 + eps) * (
                    cons.f_rb(0) + cons.c + cons.fb[t_1])), "BAD t1_2 " + ",".join(str(si) for si in [br0, c, td])
            tol = 0.5
            assert ((np.abs(t_1 - 1.0) < eps_num) or np.abs(cons.fr[t_1] + cons.c - (1 + eps) * (
                    cons.f_rb(0) + cons.c + cons.fb[t_1])) < eps_num), "BAD t1_3 " + ",".join(
                str(si) for si in [br0, c, td]) + "|" + str(t_1) + "|" + str(
                np.abs(cons.fr[t_1] + c - (1 + eps) * (cons.f_rb(0) + c + cons.fb[t_1]))) + "|" + str(
                cons.fr[1.0] + c) + "|" + str((1 + eps) * (cons.f_rb(0) + c + cons.fb[1.0]))
            # print "no assert of t_1"

            alp = cons.calc_t_N_minus_1(eps)

            assert (alp >= 1.0), "BAD alp" + ",".join(str(si) for si in [br0, c, td])
            # print alp,cons.rho_inf(),cons.fr(alp)
            flag1 = (cons.fr[alp] / cons.fb[alp]) >= -eps_num + cons.cr_inf() / (1 + eps)
            # print alp,c,cons.fr[alp]
            flag2 = (cons.c <= eps_num + eps * cons.fr[alp])
            # print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c

            assert (flag1 or flag2), "BAD flags" + ",".join(str(si) for si in [br0, c, td])

            cr0 = cons.cr_0()
            cr0_thr = div(bb0, br0)
            assert (np.abs(cr0_thr - cr0) < eps_num or (min(cr0, cr0_thr) == np.inf)), "BAD cr0" + ",".join(
                str(si) for si in [br0, c, td])

            cri = cons.cr_inf()

            cri_thr = div(ar_s[-1], ab_s[-1]) if ar_s[-1] > 0 else cons.fr[20] / cons.fb[20]
            # print cri_thr,cri
            assert (np.abs(cri_thr - cri) < eps_num or (min(cri, cri_thr) == np.inf)), "BAD cri" + ",".join(
                str(si) for si in [br0, c, td])

            cr, q, T_cons, _, _ = cons.approxOptStr(eps)
            Z[ind_x, ind_y] = cr
            R0[ind_x, ind_y] = q[0] if T_cons[0] == 0 else 0.0
            Ri[ind_x, ind_y] = q[-1] if T_cons[-1] == np.inf else 0.0
            Rm[ind_x, ind_y] = 1 - R0[ind_x, ind_y] - Ri[ind_x, ind_y]

    print "time elpassed ", time.time() - tic, "counter", cnt
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # ax = plt.gca()
    fontsize = 25

    if ind_target == 0:
        pickle.dump((C, T, Z), open("BTZ.pkl", "wb"))
        # surf = ax.plot_surface(C, T, Z, cmap=cm.coolwarm,
        #             linewidth=0, antialiased=False)
    else:
        pickle.dump((B, T, Z, R0, Rm, Ri), open("BTZ.pkl", "wb"))

        np.savetxt("x_figure6.csv", B, delimiter=",")
        np.savetxt("y_figure6.csv", T, delimiter=",")
        np.savetxt("z_figure6.csv", Z, delimiter=",")

        np.savetxt("x_figure7.csv", B, delimiter=",")
        np.savetxt("y_figure7.csv", T, delimiter=",")
        np.savetxt("z_figure7.csv", R0, delimiter=",")

        np.savetxt("x_figure8.csv", B, delimiter=",")
        np.savetxt("y_figure8.csv", T, delimiter=",")
        np.savetxt("z_figure8.csv", Ri, delimiter=",")

        np.savetxt("x_figure_rm.csv", B, delimiter=",")
        np.savetxt("y_figure_rm.csv", T, delimiter=",")
        np.savetxt("z_figure_rm.csv", Rm, delimiter=",")



        z_title = [r"$q_0$", r"$q_m$", r"$q_{\infty}$"]
        for ind, O in enumerate([R0, Rm, Ri]):
            fig = plt.figure(ind + 1)
            ax = fig.gca(projection='3d')
            ax = plt.gca()
            surf = ax.plot_surface(B, T, O, cmap=cm.coolwarm,
                                   linewidth=0, antialiased=False)
            plt.title(parm_strs[1] + "=" + str(target_c), fontsize=fontsize)
            plt.xlabel(parm_strs[0], fontsize=fontsize)
            plt.ylabel(parm_strs[2], fontsize=fontsize)
            ax.set_zlabel(z_title[ind], fontsize=fontsize)
            # fig.colorbar(surf)

    if ind_target == 1:
        plt.title(parm_strs[1] + "=" + str(target_c), fontsize=fontsize)
        plt.xlabel(parm_strs[0], fontsize=fontsize)
        plt.ylabel(parm_strs[2], fontsize=fontsize)
    else:
        plt.title(parm_strs[0] + "=" + str(target_br0), fontsize=fontsize)
        plt.xlabel(parm_strs[1], fontsize=fontsize)
        plt.ylabel(parm_strs[2], fontsize=fontsize)

    fig.colorbar(surf)
    plt.show()
    return


def plot_results_prob_cdf():
    global dict_res_slice, dict_res
    # first we need to slice
    file_name = "dict_res_uniified_br0_c_td.pkl"
    dict_res = pickle.load(open(file_name, "r"))
    # dict_res[br0,c,td] = (t_det_gen, rho_det_gen,rho_prob,  q, T)
    # genr = ((br0, c, td) for br0 in np.arange(10.0, 38.0, 2) for c in np.arange((bb0 - br0) + 0.2, bb0 - 0.1, 2) for td
    #     in np.arange((bb0 - br0) / ar1 + 0.1, 10.0, 0.5))
    print [sorted(list(set([k[i] for k in dict_res.keys()]))) for i in range(3)]
    num_iter = 10
    inds = np.random.choice(len(dict_res), num_iter, False)
    ks = [dict_res.keys()[i] for i in inds]

    parm_strs = ["Initial rent cost", "Switch cost", "Time difference"]
    # dict_res_slice = {k: dict_res[k] for k in ks}

    br0, c, td = 12.0, 30.0, 10.0
    eps = 0.05
    bb0 = 50.0
    ar1 = 5.0

    ts = [0, 10, 10 + td]
    ar_s = [0, ar1, 0]
    ab_s = [0, 0, 0]

    cons_org = create_instance(br0, bb0, ar_s, ab_s, ts, c)
    plt = cons_org.plot()
    plt.show()

    cons = cons_org.normalize()

    tup = cons.approxOptStr(eps)
    rho_prob = tup[0]
    q = tup[1]
    T = tup[2]
    C, f, A_eq, b_ub, b_eq = tup[4]
    fontsize = 20
    plt.figure(1)
    plt.clf()
    plt.title("CDF of the best probabilistic strategy", fontsize=fontsize)
    plt.xlabel("Consumption (Time)", fontsize=fontsize)
    plt.ylabel("Probability", fontsize=fontsize)
    plt.plot(cons.norm_time * np.array(T[:-1]), np.cumsum(q[:-1]), "g-*")
    print("xs",cons.norm_time * np.array(T[:-1]))
    print("ys",np.cumsum(q[:-1]))
    np.savetxt("x_figure9.csv", (cons.norm_time * np.array(T[:-1])), delimiter=",")
    np.savetxt("y_figure9.csv", (np.cumsum(q[:-1])), delimiter=",")




    plt.show()
    return

    for ind, k in enumerate(ks):
        markersize = 18
        fontsize = 20
        plt.figure(10 + ind)
        plt.clf()
        plt.title("CDF of the best probabilistic strategy", fontsize=fontsize)
        plt.xlabel("Consumption (Time)", fontsize=fontsize)
        plt.ylabel("Probability", fontsize=fontsize)
        plt.plot(dict_res[k][4][:-1], np.cumsum(dict_res[k][3][:-1]), "g-*")

    ax = plt.gca()

    # simArtists = [plt.Line2D((4, 2, 40), (4, 2.001, 40), color="k", marker=m, linestyle='-') for m in ["*", "d", "o"]]
    simArtists = [plt.Line2D([0], [0], color="k", marker=m, linestyle='-') for m in
                  ["*", "d", "o"]]

    # ax.legend(simArtists,["best startegy start with buy","best startegy switch at middle time","best strategy never buy"])
    plt.legend(simArtists,
               ["best strategy start with buy", "best strategy switch at middle time", "best strategy never buy"],
               loc=2)
    plt.show()
    return


def plot_results_det_1D_cloud():
    global dict_res_slice, dict_res
    # first we need to slice
    file_name = "dict_res_uniified_cloud_br0_c_td.pkl"
    #dict_res = pickle.load(open(file_name, "r"))
    # dict_res[br0,c,td] = (t_det_gen, rho_det_gen,rho_prob,  q, T)
    # genr = ((br0, c, td) for br0 in np.arange(10.0, 38.0, 2) for c in np.arange((bb0 - br0) + 0.2, bb0 - 0.1, 2) for td
    #     in np.arange((bb0 - br0) / ar1 + 0.1, 10.0, 0.5))
    # print [sorted(list(set([k[i] for k in dict_res.keys()]))) for i in range(3)]


    ratio_yearly = 8.0
    ratio_3_yearly=15.0
    monthly_cost = 100.0
    if False:
        target_c = ratio_yearly*monthly_cost
        bb0 = ratio_yearly * monthly_cost
        plan_duration =12
    else:
        target_c = ratio_3_yearly*monthly_cost
        bb0 = ratio_3_yearly * monthly_cost
        plan_duration = 36


    eps = 0.01
    eps_num = 1e-10
    dict_res_slice = {}

    br0=0.0
    ar1 = monthly_cost
    months_in_year=12
    c = target_c
    genr = ((td) for td in np.arange(28.0, 29.0) )
    for ind, (td) in enumerate(genr):
        if ((td - eps_num) < ((bb0 - br0) / ar1)):
            dict_res_slice[br0] = [np.inf, 1, None, None, None]
            raise Exception("not supported")
            continue

        assert (td - eps_num) > ((bb0 - br0) / ar1), "BAD TD" + ",".join(str(si) for si in [br0, c, td])

        assert (c > (bb0 - br0 - eps_num)), "BAD C" + ",".join(str(si) for si in [br0, c, td])
        if td > plan_duration:
            ts = [0, plan_duration,td]
            ar_s = [ar1,ar1,0]
            ab_s = [0, ar1, 0]
        else:
            ts = [0, td]
            ar_s = [ar1,  0]
            ab_s = [0, 0]

        cons_org = create_instance(br0, bb0, ar_s, ab_s, ts, c)
        cons = cons_org.normalize()
        assert (np.max(np.abs(np.array(
            [cons.fr.lins[i](cons.times[i + 1]) - cons.fr.lins[i + 1](cons.times[i + 1]) for i in
             range(len(cons.times) - 1)]))) < eps), "BAD fr" + ",".join(str(si) for si in [br0, c, td])
        assert (np.max(np.abs(np.array(
            [cons.fb.lins[i](cons.times[i + 1]) - cons.fb.lins[i + 1](cons.times[i + 1]) for i in
             range(len(cons.times) - 1)]))) < eps), "BAD fb" + ",".join(str(si) for si in [br0, c, td])

        assert (np.abs(cons.intersection() - 1.0) < eps_num), "BAD intersection" + ",".join(
            str(si) for si in [br0, c, td])

        t_1 = cons.calc_t1(eps)
        # print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
        assert (-eps_num < t_1 < 1.0 + eps_num), "BAD t1_1 " + ",".join(str(si) for si in [br0, c, td])
        assert ((cons.fr[t_1] + cons.c) <= eps_num + (1 + eps) * (
                cons.f_rb(0) + cons.c + cons.fb[t_1])), "BAD t1_2 " + ",".join(str(si) for si in [br0, c, td])
        tol = 0.5
        assert ((np.abs(t_1 - 1.0) < eps_num) or np.abs(cons.fr[t_1] + cons.c - (1 + eps) * (
                cons.f_rb(0) + cons.c + cons.fb[t_1])) < eps_num), "BAD t1_3 " + ",".join(
            str(si) for si in [br0, c, td]) + "|" + str(t_1) + "|" + str(
            np.abs(cons.fr[t_1] + c - (1 + eps) * (cons.f_rb(0) + c + cons.fb[t_1]))) + "|" + str(
            cons.fr[1.0] + c) + "|" + str((1 + eps) * (cons.f_rb(0) + c + cons.fb[1.0]))
        # print "no assert of t_1"

        alp = cons.calc_t_N_minus_1(eps)

        assert (alp >= 1.0), "BAD alp" + ",".join(str(si) for si in [br0, c, td])
        # print alp,cons.rho_inf(),cons.fr(alp)
        flag1 = (cons.fr[alp] / cons.fb[alp]) >= -eps_num + cons.cr_inf() / (1 + eps)
        # print alp,c,cons.fr[alp]
        flag2 = (cons.c <= eps_num + eps * cons.fr[alp])
        # print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c

        assert (flag1 or flag2), "BAD flags" + ",".join(str(si) for si in [br0, c, td])

        cr0 = cons.cr_0()
        cr0_thr = div(bb0, br0)
        assert (np.abs(cr0_thr - cr0) < eps_num or (min(cr0, cr0_thr) == np.inf)), "BAD cr0" + ",".join(
            str(si) for si in [br0, c, td])

        cri = cons.cr_inf()

        cri_thr = div(ar_s[-1], ab_s[-1]) if ar_s[-1] > 0 else cons.fr[20] / cons.fb[20]
        # print cri_thr,cri
        assert (np.abs(cri_thr - cri) < eps_num or (min(cri, cri_thr) == np.inf)), "BAD cri" + ",".join(
            str(si) for si in [br0, c, td])

        t_det_gen, rho_det_gen = cons.optDetStr()
        cr, q, T_cons, _, _ = cons.approxOptStr(eps)
        dict_res_slice[br0] = (t_det_gen, rho_det_gen, None, None, None)
        print td,t_det_gen,rho_det_gen,cr,q[0],q[-1]
        print len(T_cons[:-1]),len(np.cumsum(q[:-1]))
        np.savetxt("fig11_time.csv", np.array(T_cons[:-1])*15, delimiter=",")
        np.savetxt("fig11_prob.csv", np.cumsum(q[:-1]), delimiter=",")

    return
    markersize = 18
    fontsize = 20
    normalize = matplotlib.colors.Normalize(vmin=1, vmax=1.8)
    dict_marker = {np.inf: "o", 1.0: "d", 0.0: "*"}
    plt.figure(10)
    plt.clf()
    dict_size = {"o": 100, "d": 100, "*": 130}
    for km in dict_marker.keys():
        pass
        """plt.scatter(x=[k for k in dict_res_slice.keys() if dict_res_slice[k][0] == km],
                        y=[dict_res_slice[k][1] for k in dict_res_slice.keys() if dict_res_slice[k][0] == km], marker=dict_marker[km],
                        c=[(dict_res_slice[k][1]) for k in dict_res_slice.keys() if dict_res_slice[k][0] == km],
                        cmap='plasma', s=dict_size[dict_marker[km]],norm=normalize,linewidth=0)"""
    xs = sorted(dict_res_slice.keys())
    plt.plot(xs, [dict_res_slice[x][1] for x in xs], "r-")
    print("xs: ",xs)
    print("ys: ",[dict_res_slice[x][1] for x in xs])
    print [
        "change from {0} to {1} at {2} ({3})".format(dict_res_slice[xs[i]][0], dict_res_slice[xs[i + 1]][0], xs[i], i)
        for i in range(len(xs) - 1) if dict_res_slice[xs[i]][0] != dict_res_slice[xs[i + 1]][0]]
    fs = 12
    plt.plot([24, 24], [0, 1.68], "k--")
    plt.text(20, 1.2, "Never switch", {"fontsize": fs})
    plt.plot([29.6, 29.6], [0, 1.68], "k--")
    plt.text(24.1, 1.4, "Switch at intersection time", {"fontsize": fs})
    plt.text(32, 1.0, "Start at buy", {"fontsize": fs})

    # plt.colorbar()
    plt.title(parm_strs[2] + "=" + str(target_td) + ", " + parm_strs[1] + "=" + str(target_c), fontsize=fontsize)
    plt.xlabel(parm_strs[0], fontsize=fontsize)
    plt.ylabel("Competitive ratio", fontsize=fontsize)
    ks = sorted(dict_res_slice.keys())
    ys = [dict_res_slice[k][1] for k in ks]
    # plt.plot(ks,[dict_res_slice[k][1] for k in ks],"k-")
    plt.axis([15, 48, 1, 1.8])
    ax = plt.gca()

    # simArtists = [plt.Line2D((4, 2, 40), (4, 2.001, 40), color="k", marker=m, linestyle='-') for m in ["*", "d", "o"]]
    simArtists = [plt.Line2D([0], [0], color="k", marker=m, linestyle='-') for m in
                  ["*", "d", "o"]]

    # ax.legend(simArtists,["best startegy start with buy","best startegy switch at middle time","best strategy never buy"])
    # plt.legend(simArtists,["best startegy start with buy","best startegy switch at middle time","best strategy never buy"],loc=2)
    plt.show()
    return

    ax.set_xlabel("Switch cost", fontsize=fontsize)
    ax.set_ylabel("Time diffrenece", fontsize=fontsize)
    ax.set_zlabel("Initial rent price", fontsize=fontsize)
    ax.set_title("Determinstic strategies properties", fontsize=25)
    simArtists = [plt.Line2D((4, 2, 40), (4, 2.001, 40), color="k", marker=m, linestyle='-', markersize=15) for m in
                  ["*", "d", "o"]]

    ax.legend(simArtists,
              ["best startegy start with buy", "best startegy switch at middle time", "best strategy never buy"])
    # ax.legend(sc)
    fig.colorbar(sc)
    return






flag, flag2 = None, None
# flag = UT_linear()
if not flag:
    print "Linear tests  has not been done"

# flag2 = UT_PWL_5params()
if not flag2:
    print "5 params tests  has not been done"

flag = flag and flag2
if flag:
    print "Good Test"
print ""

# plot_results_det_2D()


#plot_results_det_2D_surf()

#plot_results_det_1D()
#plot_results_prob_cdf()
# plot_results_prob()
# exmp_unifed()
#example_insatnce()

# bb0_init problem
#plot_results_prob_2D_surf()

# exmp_p6()
#exmp_det_for_middle_time_graph()
plot_results_det_1D_cloud()