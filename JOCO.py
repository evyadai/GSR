from sympy import *
a1,a2,b1,b2,c,r,T,t,l = symbols('a1 a2 b1 b2 c r T t l')
init_printing(use_unicode=True)

#P6-L33
R0 = (b2+a2/r*(1-exp(-r*T)))/(b1+a1/r*(1-exp(-r*T)))
res0 = diff(R0,T)
#print(res)

#P6-L39
Ri = (b1+a1/r*(1-exp(-r*T)))/(b2+a2/r*(1-exp(-r*T)))
res = diff(Ri,T)
#print(res)

#P7-L2
#if T<t then CR=1 so T>=t
R1t = (b1+a1/r*(1-exp(-r*t))+c*exp(-r*t)+a2/r*(exp(-r*t)-exp(-r*T)))/(b1+a1/r*(1-exp(-r*T)))
res = diff(R1t,T)
#print(res)

#P7-L18 (t>= T*)
#if T<T* then offline choose 1 and CR=1
#if T*<T<t then offline choose 2 and cost is 1 so we get
# Ri which increase at time (so adv choose t)
#if T>= t then offline choose 2 and the cost is:
R2t = (b1+a1/r*(1-exp(-r*t))+c*exp(-r*t)+a2/r*(exp(-r*t)-exp(-r*T)))/(b2+a2/r*(1-exp(-r*T)))
res = diff(R2t,T)
#print(res)

#P7-L21
R2 = (b1+a1/r*(1-exp(-r*t))+c*exp(-r*t))/(b2+a2/r*(1-exp(-r*t)))
res = diff(R2,t)

#P7-L21
ta = -1/r*log(1-((b2-b1)*r)/(a1-a2))
R2_ta = (b1+a1/r*(1-exp(-r*ta))+c*exp(-r*ta))/(b2+a2/r*(1-exp(-r*ta)))
"""
#P8-L23
#b2=0.8,b1=0.1,c=4,cr*=4.5
#lambda=1.01 - I_lambda:cr<4.545/
#if t is small(<1) and advrsary choose high value(>1)
CR1_high = (0.1+0.9*t+4+0.7*(T-t))/(0.8+0.2*T)
Tadv_high=1
CR1_Tadv_high = (0.1+0.9*t+4+0.7*(Tadv_high-t))/(0.8+0.2*Tadv_high)
#which is 0.2t+4.8 >4.8
#if t is small(<1) and advrsary choose low value(<1)
CR1_low = (0.1+0.9*t+4+0.7*(T-t))/(0.1+0.9*T)
Tadv=t
CR1_Tadv_low = (0.1+0.9*t+4)/(0.1+0.9*t)
#which is (0.9t+4.1)/(0.9t+0.1) => 5<CR<41
#if t is high(>1)and adv could choose before t
CR2_low =  (0.1+0.9*T)/(0.8+0.2*T)
Tadv=1
CR2_Tadv_low = 0.1+0.9*t #=> 

#if t is high(>1)and adv could choose after t
CR2_high =  (0.1+0.9*t+4+0.7*(T-t))/(0.8+0.2*T)
Tadv=t
CR2_Tadv_high = (0.1+0.9*t+4)/(0.8+0.2*t) # 4.5<CR<5.125 => t>62

#if t  is never then adv choose >1
CR3_high = (0.1+0.9*T)/(0.8+0.2*T) #TADV=Inf,CR=4.5

#now we know that T<1 then
CR_risk = 1
"""

"""
#P8-L23
#b2=0.8,b1=0.1,c=3,cr*=4
#lambda=1.1 - I_lambda:cr<4.4
#if t is small(<1) and advrsary choose high value(>1)
CR1_high = (0.1+0.9*t+3+0.7*(T-t))/(0.8+0.2*T)
Tadv_high=1
CR1_Tadv_high = (0.1+0.9*t+3+0.7*(Tadv_high-t))/(0.8+0.2*Tadv_high)
#which is 0.2t+3.8 >3.8
#if t is small(<1) and advrsary choose low value(<1)
CR1_low = (0.1+0.9*t+3+0.7*(T-t))/(0.1+0.9*T)
Tadv=t
CR1_Tadv_low = (0.1+0.9*t+3)/(0.1+0.9*t) 
#which is (0.9t+3.1)/(0.1t+0.9) => 4<CR<31 => 0.87<t<1

#if t is high(>1)and adv could choose before t
CR2_low =  (0.1+0.9*T)/(0.8+0.2*T)
Tadv=1
CR2_Tadv_low = 1 #=> 

#if t is high(>1)and adv could choose after t
CR2_high =  (0.1+0.9*t+3+0.7*(T-t))/(0.8+0.2*T)
Tadv=t
CR2_Tadv_high = (0.1+0.9*t+3)/(0.8+0.2*t) # 4<CR<4.5 => 1<t<20

#if t  is never then adv choose >1
CR3_high = (0.1+0.9*T)/(0.8+0.2*T) #TADV=Inf,CR=4.5

#I_lambda   0.87<t<20
#now we know that T<1 then 
CR_risk1 = /(0.1+0.9*T)
"""
#P8-L23
#b2=0.9,b1=0.8,c=3,cr*=1.125
#lambda=1.1 - I_lambda:cr<1.2375
#if t is small(<1) and advrsary choose high value(>1)
CR1_high = (0.9+0.1*t+3+0.1*(T-t))/(0.8+0.2*T)
Tadv_high=1
CR1_Tadv_high = (0.9+0.1*t+3+0.1*(Tadv_high-t))/(0.8+0.2*Tadv_high)
#which is 4
#if t is small(<1) and advrsary choose low value(<1)
CR1_low = (0.9+0.1*t+3+0.1*(T-t))/(0.9+0.1*T)
Tadv=t
CR1_Tadv_low = (0.9+0.1*t+3)/(0.9+0.1*t) 
#which is (0.9t+3.1)/(0.1t+0.9) => 4<CR<4.3333 

#if t is high(>1)and adv could choose before t
CR2_low =  (0.9+0.1*T)/(0.8+0.2*T)
Tadv=1
CR2_Tadv_low = 1 #=> 

#if t is high(>1)and adv could choose after t
CR2_high =  (0.9+0.1*t+3+0.1*(T-t))/(0.8+0.2*T)
Tadv=t
CR2_Tadv_high = (0.9+0.1*t+3)/(0.8+0.2*t) # 1.125<CR<4 => t>4.875

#if t  is never then adv choose >1
CR3_high = (0.9+0.1*T)/(0.8+0.2*T) #TADV=Inf,CR=4.5

#I_lambda   0.87<t<20
#now we know that T<1 then 


#P9-L4
u = 1+(c*r)/((a1+b1*r)*exp(r*t)-a1)
v=l*(1+(c*(a1-a2)-c*r*(b2-b1))/(a1*b2-a2*b1))
res=solve(u-v,t)
res=res[0]
res=exp(res*r)
res=res-a1/(a1+b1*r)
res=c*r*(a1*b2-a2*b1)/res
res=res/(a1+b1*r)
res=simplify(res)

#P9-L9
t1=solve(u-v,t)[0]
Ta=-1/r*ln(1-(r*(b2-b1))/(a1-a2))

diff = exp(Ta*r)-exp(t1*r)
diff=simplify(diff*(a1+b1*r)*(a1-a2-r*(b2-b1)))
diff=simplify(diff*(a1*b2*l-a1*b2+a1*c*l-a2*b1*l+a2*b1-a2*c*l+b1*c*l*r-b2*c*l*r))


