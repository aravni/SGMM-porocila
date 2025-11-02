#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Sestavljena gonila v mobilni tehniki
1. poročilo: Dinamična karakteristika vozila
Optimizacija prestavnega razmerja ali
prostih parametrov simulacije (
izkoristek prenosa 'eta' in 
delež maksimalnega pospeška pri speljevanju 'p').
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_simpson
from scipy.optimize import root, root_scalar

# funkcija za branje podatkov
from ReadData import *

# custom funkcija za poravnavo y osi
from alignYaxes import *

# custom nastavitve plottanja
plt.style.use('aljaz')

#####################################################
# podatki
#####################################################

# vozilo
r_st,m_k,f,C_d,A,L,m,eta,i_final,i_ratios,n,M,P_kW = ReadData('./podatki/volvo_480ES.dat')
# r_st,m_k,f,C_d,A,L,m,eta,i_final,i_ratios,n,M,P_kW = ReadData('./podatki/volvo_480TURBO.dat')
# r_st,m_k,f,C_d,A,L,m,eta,i_final,i_ratios,n,M,P_kW = ReadData('./podatki/quattroporte.dat')

r_din = r_st #pribl.!
p = 1 #delež max pospeška pri speljevanju - prosti parameter!

# konstante okolja
alpha=np.deg2rad(0) #naklon, rad
rho=1.225 #gostota zraka, kg/m^3
g=9.81 #gravitacija, m/s^2

#####################################################
# optimizacija prestavnega razmerja
#####################################################

def pospesek(i_ratios):
    '''prestavna razmerja -> željena t do v=100 km/h in v_max'''
    global r_st,m_k,f,C_d,A,m,eta,i_final,n,M,P_kW
    global alpha, rho, g, p, r_din

    ########### 1. hitrosti po prestavah
    v_all = [] #seznam hitr. po prestavah iz meritev obratov
    for i in i_ratios:
        vi = 2*np.pi*n*60 / (i*i_final) *r_din/1e6 #km/h
        v_all.append(vi)
    v_all=np.asarray(v_all)

    ############ 2. upori
    R_f = f*m*g*np.cos(alpha)/1e3 #kotalni, kN
    R_s = m*g*np.sin(alpha)/1e3 #upor strmine, kN
    fR_z = lambda v: 0.5*A/1e6*C_d*rho*(v/3.6)**2/1e3 #zračni upor, kN -> param. v: hitrost [km/h]
    v_range=np.linspace(0,np.max(v_all),len(M)) #to niso meritve!

    # upor vztr. mas
    I_k = 0.5*m_k*(r_st)**2 #MVM kolesa
    k1 = 4*I_k*g*1e-3/(r_st*r_din) #kg kN /s^2 = kN
    k2=0.007 #ocena
    delta_f = lambda i: 1+k1+k2*i**2 #2.način (bolje); i -> prestavno razmerje

    ########### 3. dinamični vozni faktor
    F_K=[]
    D=[] #dinamični vozni faktor
    for i in range(len(i_ratios)):
        FKi= M*i_ratios[i]*i_final*eta/r_din
        F_K.append(FKi)

        Fi = F_K[i] - fR_z(v_all[i]) -R_f-R_s
        D.append(Fi/(m*g*1e-3))
    F_K=np.asarray(F_K)
    D=np.asarray(D)

    ########### 4. pospeški
    a = D*g #pospeški po prestavah
    for i in range(a.shape[0]):
        a[i] /= delta_f(i_ratios[i]) #upoštevanje ocene vztr. rotirajočih mas

    a_real=[] #krivulja največjih pospeškov
    v_real=[]

    a0 = np.max(a[0])*p #speljevanje - konst. pospešek
    v_real.extend(np.linspace(0, v_all[0,np.argmin(np.abs(a[0]-a0))], 3))
    a_real.extend(np.full_like(v_real, a0))

    for i in range(0, v_all.shape[0]-1):
        v_old = v_real[-1] #dosedanja največja hitrost -> gledamo eno prej
        v_new = v_all[i+1,0] #najmanjša hitrost v večji prestavi
        k=0 #index za pospešek v večji prestavi za točke nad v_old
        for j in range(v_all.shape[1]-1):
            if v_all[i,j] > v_all[i+1,k]:
                k+=1 #prištejemo nov indeks da je hitrost usklajena
            if v_all[i,j] > v_old:
                if v_all[i,j] < v_new: #to je edini izračun pospeška pri dani hitrosti
                    v_real.append(v_all[i,j])
                    a_real.append(a[i,j])
                else: #preverba da je to res a_max pri danem v, ki je pozitiven
                    if a[i,j] > a[i+1,k] and a[i,j] >= 0:
                        v_real.append(v_all[i,j])
                        a_real.append(a[i,j])
    #zadnja prestava do konca le, če je pospešek pozitiven
    k=len(a[-1, v_all[-1] <= v_real[-1]])
    for j in range(1,len(a[-1, v_all[-1] > v_real[-1]])):
        if a[-1,k] >= 0:
            v_real.append(v_all[-1,k])
            a_real.append(a[-1,k])
        k+=1
    v_real=np.asarray(v_real)
    a_real=np.asarray(a_real)

    ########## 5. čas in pot pospeševanja
    t = cumulative_simpson(1/a_real, x=v_real/3.6, initial=0) #čas pospeševanja, s

    # t_uporabno = t[:np.argmax(t)] #le naraščajoči del
    # s = cumulative_simpson(v_real[:len(t_uporabno)]/3.6, x=t_uporabno, initial=0)*1e-3 #pot pospeševanja, km
    return t, v_real, a_real
def napaka(i_ratios, t_opt=None, v_opt=None):

    t, v, a = pospesek(i_ratios)
    v_limit = v[np.argmin((v-100)**2)] #najbližje 100 km/h
    n_100 = len(v[v < v_limit]) #indeks ko v==100 km/h
    # print(n_100)
    # print(v_limit, t[n_100])

    out = np.zeros_like(i_ratios)
    if t_opt is not None: 
        out[0] = t[n_100]-t_opt #čas do 100 km/h
        # return np.full_like(i_ratios, t[n_100]-t_opt)
    if v_opt is not None: out[1] = v[-1]-v_opt #v_max
    return out

t_opt = 15 #željen pospešek do 100 km/h, m/s^2 -> prilagodijo se prestavna razmerja
v_opt = 190 #željena max hitrost
sol = root(napaka, i_ratios, args=(t_opt,v_opt),method='lm')
i_ratios_opt = sol.x

t,v,a = pospesek(i_ratios_opt)
v_limit = v[np.argmin((v-100)**2)] #najbližje 100 km/h
n_100 = len(v[v < v_limit]) #indeks ko v==100 km/h

print(sol.success, ':', sol.message, 'Iteracije:',sol.nfev)
# print(f'i_opt:\t {i_ratios_opt}')
print(f"i_opt:\t {', '.join(f'{x:.3f}' for x in i_ratios_opt)}")
print(f'i_0:\t {i_ratios}')
print(f't_100 = {t[n_100]:.5g} s')
print(f'v_max = {v[-1]:.5g} km/h')
if t_opt is not None: print(f't_100 - t_opt = {(t[n_100]-t_opt):.3E} s') #napaka optimizacije - odstopanje časa od t_opt!
if v_opt is not None: print(f'v_max - v_opt = {(v[-1]-v_opt):.3E} km/h') #napaka optimizacije - odstopanje hitrosti od v_opt!
# print(sol)
input('end')

#####################################################
# optimizacija prostih parametrov
#####################################################

def pospesek(eta,p):
    '''prosti parametri (eta,p) -> željena t do v=100 km/h in v_max'''
    global r_st,m_k,f,C_d,A,m,i_final,n,M,P_kW, i_ratios
    global alpha, rho, g, r_din

    ########### 1. hitrosti po prestavah
    v_all = [] #seznam hitr. po prestavah iz meritev obratov
    for i in i_ratios:
        vi = 2*np.pi*n*60 / (i*i_final) *r_din/1e6 #km/h
        v_all.append(vi)
    v_all=np.asarray(v_all)

    ############ 2. upori
    R_f = f*m*g*np.cos(alpha)/1e3 #kotalni, kN
    R_s = m*g*np.sin(alpha)/1e3 #upor strmine, kN
    fR_z = lambda v: 0.5*A/1e6*C_d*rho*(v/3.6)**2/1e3 #zračni upor, kN -> param. v: hitrost [km/h]

    # upor vztr. mas
    I_k = 0.5*m_k*(r_st)**2 #MVM kolesa
    k1 = 4*I_k*g*1e-3/(r_st*r_din) #kg kN /s^2 = kN
    k2=0.007 #ocena
    delta_f = lambda i: 1+k1+k2*i**2 #2.način (bolje); i -> prestavno razmerje

    ########### 3. dinamični vozni faktor
    F_K=[]
    D=[] #dinamični vozni faktor
    for i in range(len(i_ratios)):
        FKi= M*i_ratios[i]*i_final*eta/r_din
        F_K.append(FKi)

        Fi = F_K[i] - fR_z(v_all[i]) -R_f-R_s
        D.append(Fi/(m*g*1e-3))
    F_K=np.asarray(F_K)
    D=np.asarray(D)

    ########### 4. pospeški
    a = D*g #pospeški po prestavah
    for i in range(a.shape[0]):
        a[i] /= delta_f(i_ratios[i]) #upoštevanje ocene vztr. rotirajočih mas

    a_real=[] #krivulja največjih pospeškov
    v_real=[]

    a0 = np.max(a[0])*p #speljevanje - konst. pospešek
    v_real.extend(np.linspace(0, v_all[0,np.argmin(np.abs(a[0]-a0))], 3))
    a_real.extend(np.full_like(v_real, a0))

    for i in range(0, v_all.shape[0]-1):
        v_old = v_real[-1] #dosedanja največja hitrost -> gledamo eno prej
        v_new = v_all[i+1,0] #najmanjša hitrost v večji prestavi
        k=0 #index za pospešek v večji prestavi za točke nad v_old
        for j in range(v_all.shape[1]-1):
            if v_all[i,j] > v_all[i+1,k]:
                k+=1 #prištejemo nov indeks da je hitrost usklajena
            if v_all[i,j] > v_old:
                if v_all[i,j] < v_new: #to je edini izračun pospeška pri dani hitrosti
                    v_real.append(v_all[i,j])
                    a_real.append(a[i,j])
                else: #preverba da je to res a_max pri danem v, ki je pozitiven
                    if a[i,j] > a[i+1,k] and a[i,j] >= 0:
                        v_real.append(v_all[i,j])
                        a_real.append(a[i,j])
    #zadnja prestava do konca le, če je pospešek pozitiven
    k=len(a[-1, v_all[-1] <= v_real[-1]])
    for j in range(1,len(a[-1, v_all[-1] > v_real[-1]])):
        if a[-1,k] >= 0:
            v_real.append(v_all[-1,k])
            a_real.append(a[-1,k])
        k+=1
    v_real=np.asarray(v_real)
    a_real=np.asarray(a_real)

    ########## 5. čas in pot pospeševanja
    t = cumulative_simpson(1/a_real, x=v_real/3.6, initial=0) #čas pospeševanja, s

    # t_uporabno = t[:np.argmax(t)] #le naraščajoči del
    # s = cumulative_simpson(v_real[:len(t_uporabno)]/3.6, x=t_uporabno, initial=0)*1e-3 #pot pospeševanja, km
    return t, v_real, a_real
def napaka(data, t_opt=None, v_opt=None):
    eta=data[0]
    p=data[1]
    t, v, a = pospesek(eta,p)
    v_limit = v[np.argmin((v-100)**2)] #najbližje 100 km/h
    n_100 = len(v[v < v_limit]) #indeks ko v==100 km/h
    # print(n_100)
    # print(v_limit, t[n_100])

    out = np.zeros_like(data)
    if t_opt is not None: 
        # print('rez:',v_limit,t[n_100], t_opt, t[n_100]-t_opt)
        out[0] = t[n_100]-t_opt #čas do 100 km/h
        # return np.full_like(i_ratios, t[n_100]-t_opt)
    if v_opt is not None: out[1] = v[-1]-v_opt #v_max
    return out

alpha=np.deg2rad(0) #naklon, rad

t_opt = 10.5 #željen pospešek do 100 km/h, m/s^2 -> prilagodijo se prosti parametri
v_opt = 190 #željena max hitrost
x0=[eta,p] #zač. pribl, [eta,p]
sol = root(napaka, x0, args=(t_opt,v_opt), method='lm')
eta_opt, p_opt = sol.x

t,v,a = pospesek(eta_opt, p_opt)
v_limit = v[np.argmin((v-100)**2)] #najbližje 100 km/h
n_100 = len(v[v < v_limit]) #indeks ko v==100 km/h

print(sol.success, ':', sol.message, 'Iteracije:',sol.nfev)
print(f'eta_opt:\t {eta_opt}')
print(f'p_opt:\t {p_opt}')
print(f't_100 = {t[n_100]:.5g} s')
print(f'v_max = {v[-1]:.5g} km/h')
if t_opt is not None: print(f't_100 - t_opt = {(t[n_100]-t_opt):.3E} s') #napaka optimizacije - odstopanje časa od t_opt!
if v_opt is not None: print(f'v_max - v_opt = {(v[-1]-v_opt):.3E} km/h') #napaka optimizacije - odstopanje hitrosti od v_opt!
# print(sol)

eta_r=np.linspace(0.8, 1.0, 100)
p_r=np.linspace(0.8, 1.0, 100)

t_opt = 10.5 #željen pospešek do 100 km/h, m/s^2 -> prilagodijo se prosti parametri
v_opt = 190 #željena max hitrost

cenilka=np.zeros((len(eta_r), len(p_r)))

for i,eta_i in enumerate(eta_r):
    for j,p_i in enumerate(p_r):
        err=napaka([eta_i,p_i], t_opt=t_opt, v_opt=v_opt)
        cenilka[i,j] = sum(abs(err))
cenilka/=np.max(cenilka)

plt.imshow(cenilka, cmap='turbo', aspect='equal', origin='lower')
# plt.plot(len(p_r)*p_opt, len(eta_r)*eta_opt, 'kx')

plt.xticks(ticks=np.linspace(0, len(p_r) - 1, 5), labels=np.round(np.linspace(p_r[0], p_r[-1], 5), 2))
plt.yticks(ticks=np.linspace(0, len(eta_r) - 1, 5), labels=np.round(np.linspace(eta_r[0], eta_r[-1], 5), 2))

plt.xlabel("$p$")
plt.ylabel("$\\eta$")

plt.colorbar(label="normalizirana napaka")
plt.show()