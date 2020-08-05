# -*- coding: utf-8 -*-
import math as mt
import numpy as np
import random as rd
import scipy.stats as st
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from pylab import *
#Premier pricer

#Définitipn du coefficient binomial 
def coeff_binom(n,p):
    a = mt.factorial(n)
    b = mt.factorial(p)
    c = mt.factorial(n-p)
    return a/(b*c)

def price1(N,rn,hn,bn,s,f):
    q = (rn-bn)/(hn-bn) # calcul de qn
    const = 1/(1+rn)**N # calcul de la constante devant la somme
    sum = 0
    for i in range(0,N,1):
        coef = coeff_binom(N,i) #calcul du coefficient binomial
        x_i= s*((1+hn)**i)*((1+bn)**(N-i))  #calcul de x_i
        a = q**i 
        b = (1-q)**(N-i)
        sum = sum  + f(x_i)*coef*a*b #somme comme trouvé q2
    return const*sum

#Définition de la fonction max(100-x,0)
def f(x):
    return max(100-x,0)

# Test avec les valeurs du sujet 
    
s = 100
hn = 0.05
bn = -0.05
rn = 0.01
N1 = 30

a = price1(N1,rn,hn,bn,s,f)
print("Le premier pricer vaut:",a)

#Deuxième pricer


def price2(N,rn,hn,bn,s,f):
    
    #initialisation des variables sous la forme d'un tableau
    tableau = np.zeros((N+1,N+1))
    for i in range(0,N+1,1):
        tableau[i][0] = f(s*((1+bn)**i)*((1+hn)**(N-i)))
    
    #calculs des vecteurs successifs
    const = 1/(1+rn)
    q = (rn-bn)/(hn-bn)
    for j in range(1,N+1,1):
        for i in range(j,N+1,1):
            tableau[i][j] = const*(q*tableau[i-1][j-1]+(1-q)*tableau[i][j-1]) 
    res = tableau[N][N]
    return res

    

#Tests pour le pricer 2
N2=3
b = price2(N2,rn,hn,bn,s,f)
print("Le pricer 2 vaut:",b)    
    
#Comparaison avec cette fois ci max(x-100,0)
def f_2(x):
    return max(x-100,0)

a = rd.randint(5,15)

pricer_1 = price1(a,rn,hn,bn,s,f_2)
pricer_2 = price2(a,rn,hn,bn,s,f_2)

print("N=",a)
print("Le pricer vaut:",pricer_1)
print("Le pricer vaut:",pricer_2)
print("L'écart en valeure absolue entre les deux méthode est de",abs(pricer_1-pricer_2))



#Couverture

rn = 0.03

a = 1+hn
b = 1+bn
c=1/(1+rn)
qN = (rn-bn)/(hn-bn)

v1_hn = c*(f_2(s*(a**2))*qN+f_2(s*a*b)*(1-qN))
v1_bn = c*(f_2(s*a*b)*qN+f_2(s*(b**2))*(1-qN))

alpha0 = (v1_bn-v1_hn)/(s*(bn-hn))
beta0  = (v1_hn*(1+bn)-v1_bn*(1+hn))/((1+rn)*(bn-hn))
print("alpha0=",alpha0,"et beta0=",beta0)

#Pour s(1+hn)
alpha1 = (f_2(s*(a**2))-f_2(s*a*b))/(s*a*(hn-bn))
beta1 = (a*f_2(s*b*a)-b*f_2(s*(a**2)))/(((1+rn)**2)*(hn-bn))
print("Pour St1 = 1+hn on a alpha1=",alpha1,"et beta1=",beta1)

#Pour s(1+bn)
alpha1_bis = (f_2(s*a*b)-f_2(s*(b**2)))/(s*b*(hn-bn))
beta1_bis = (a*f_2(s*(b**2))-b*f_2(s*a*b))/(((1+rn)**2)*(hn-bn))
print("Pour St1 = 1+hn on a alpha1=",alpha1_bis,"et beta1=",beta1_bis)


#pricer 3    
def price3(n,s,r,sig,T,f):
    xi = st.norm.rvs(0,1,n)
    c1 = mt.exp(-r*T)
    c2 = mt.exp(T*((r-(sig**2)/2)))
    c3 = sig*mt.sqrt(T)
    somme = 0
    for i in range(0,n,1):
        c4 = mt.exp(c3*xi[i])
        c = s*(c2*c4)
        somme += f(c)
    return (c1*somme)/n

r=0.01
sig=0.1
s=100
T=1
res = np.zeros(10)
n = (10**5)*np.arange(1,11,1)


res = [price3(i,s,r,sig,T,f) for i in n]


plt.plot(n,res)
plt.xlabel('n')
plt.ylabel('Prix')
plt.show()

#Le pricer par formule fermée
    # on utilise st.norm.cdf pour la fonction de répartition d'une loi normale 
def put(s,r,sig,T,K):
    c1 = sig*mt.sqrt(T)
    d = (1/c1)*(mt.log((s/K))+T*(r+(sig**2)/2))
    return -s*st.norm.cdf(-d,0,1)+K*mt.exp(-r*T)*st.norm.cdf(-d+c1,0,1)

r=0.01
sig=0.1
s=100
T=1
K=100


#résultat fonction put

blackScholes = put(s,r,sig,T,K)
print("Par la formule de BS, le put vaut",blackScholes)

#Graphique question 17
#affichage graphique pour comparaison

bs =[blackScholes for i in range(0,10)]

#graphe de comparaison

plt.plot(n,res)
plt.plot(n,bs)
plt.xlabel('n')
plt.ylabel('Prix')
plt.show()

#Graphique question 18

s = [20*i for i in range(1,11,1)]
T = [1/12,1/6,1/4,1/3,1/2,1]


#graphe 3D

fig = figure()
ax = fig.gca(projection='3d')
for i in range(len(s)):
    for j in range(len(T)):
        ax.scatter(s[i],T[j],put(s[i],r,sig,T[j],K))

ax.set_xlabel('s')
ax.set_ylabel('T')
ax.set_zlabel('Prix')
plt.show()

#Question 19

#mouvement brownien
def Bt(T):
    return rd.gauss(0,T)


def S(s,r,sig,t,Bt):
    c1 = (r-(sig**2)/2)*t
    c2 = sig*Bt
    return s*mt.exp(c1+c2)

#fonction max(100-sT,0)    
def sT(t):
    T=1
    s=100
    r=0.02
    sig=0.3
    B=Bt(T)
    return max(100-S(s,r,sig,t,B),0)


N = [10*k for k in range(1,100,1)]
s=100
sigma=0.3
r=0.02
T=1
B = Bt(T)
print(B)

#definition de rN,hN et bN
rN=[(r*T)/N[i] for i in range(len(N))]
hN=[(((1+rN[i])*mt.exp((sigma*mt.sqrt(T))/mt.sqrt(N[i])))-1) for i in range(len(N))]
bN=[((1+rN[i])*mt.exp(-(sigma*mt.sqrt(T))/mt.sqrt(N[i]))-1) for i in range(len(N))]

    
#liste des prix avec price2
price = [price2(N[i],rN[i],hN[i],bN[i],s,sT) for i in range(len(N))]

#liste des prix avec put
p = put(s,r,sigma,T,K)   
pN = [p for i in range(len(N))]


plt.plot(N,price)
plt.plot(N,pN)
plt.xlabel('N')
plt.ylabel('Prix')
plt.show()