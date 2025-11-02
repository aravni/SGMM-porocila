#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Sestavljena gonila v mobilni tehniki
1. poročilo: Dinamična karakteristika vozila
Obravnava zdrsa
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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
r_st,m_k,f,C_d,A,L,m,eta,i_final,i_ratios,n,M,P_kW = ReadData('./podatki/volvo_480TURBO.dat')
r_st,m_k,f,C_d,A,L,m,eta,i_final,i_ratios,n,M,P_kW = ReadData('./podatki/quattroporte.dat')

r_din = r_st #pribl.!
eta = 0.9 #ocena! -> velik vpliv na v_max, manj na t_100 - prosti parameter! NI PODATKA
p = 1 #delež max pospeška pri speljevanju - prosti parameter!
h_t = 750 #višina težišča, mm - prosti parameter! NI PODATKA
m1 = 0.5 #porazdelitev mase na ravnini, delež skupne mase na sprednjem kolesu - prosti parameter! NI PODATKA
mu = 0.6 #max koef trenja (sojemanja)

# konstante okolja
alpha=np.deg2rad(0) #naklon, rad
rho=1.225 #gostota zraka, kg/m^3
g=9.81 #gravitacija, m/s^2

roman=['I','II','III','IV','V','VI','VII','VIII','IX','X']

#####################################################
# zdrs
#####################################################

# sili teže na kolesih
G1= m1*m*g
G2= (1-m1)*m*g

# razdalji koles do težišča
L1 = L*(1-m1)
L2 = L*m1 #=L-L1

# sila trenja = max možna sila na kolesu, kN
F_tr_4x4 = mu*m*g*np.cos(alpha)*1e-3 #za 4x4
F_tr_z = mu*m*g*np.cos(alpha)*((L1-h_t*f)/(L-mu*h_t))*1e-3 #zadnji pogon
F_tr_p = mu*m*g*np.cos(alpha)*((L2+h_t*f)/(L+mu*h_t))*1e-3 #prednji pogon
F_K_max = F_tr_4x4 # izbira glede na pogon!!!
print(F_tr_p,F_tr_z,F_tr_4x4)

v_all = [] #seznam hitr. po prestavah iz meritev obratov
for i in range(len(i_ratios)):
    vi = 2*np.pi*n*60 / (i_ratios[i]*i_final) *r_din/1e6 #km/h
    v_all.append(vi)
v_all=np.asarray(v_all)

def sile(alpha, F_K_max):
    global r_st,m_k,f,C_d,A,L,m,eta,i_final,i_ratios,n,M,P_kW
    global r_din,p,h_t,m1,mu, rho,g
    global L1,L2, v_all
    global F_K1_max, F_K2_max

    # upori
    R_f = f*m*g*np.cos(alpha)/1e3 #kotalni, kN
    R_s = m*g*np.sin(alpha)/1e3 #upor strmine, kN
    fR_z = lambda v: 0.5*A/1e6*C_d*rho*(v/3.6)**2/1e3 #zračni upor, kN -> param. v: hitrost [km/h]

    # upor vztr. mas
    I_k = 0.5*m_k*(r_st)**2 #MVM kolesa
    k1 = 4*I_k*g*1e-3/(r_st*r_din) #kg kN /s^2 = kN
    k2=0.007 #ocena
    delta_f = lambda i: 1+k1+k2*i**2 #2.način (bolje); i: prestavno razmerje
    fR_i_tr = lambda a: m*a*1e-3 #upor translatornih vztr. mas, kN; a: pospešek [m/s^2]
    fR_i_rot = lambda a,i: m*a*(delta_f(i)-1.)*1e-3 #upor rotirajočih vztr. mas, kN; a: pospešek [m/s^2]
    fR_i = lambda a,i: m*a*delta_f(i)/1e3 #upor vztr. mas, kN; a: pospešek [m/s^2]

    # sila na kolesu
    F_K= M*np.transpose([i_ratios])*i_final*eta/r_din
    F_K_dej=np.where(F_K > F_K_max, F_K_max, F_K)
    F_rez = F_K - fR_z(v_all) -R_f-R_s #R_i je upoštevan kasneje pri izračunu pospeškov!
    D=(F_K_dej - fR_z(v_all) -R_f-R_s) / (m*g*1e-3)

    # ravnovesje v smeri vožnje, sili podlage
    Z2 = (m*g*np.cos(alpha)*L1 + h_t*(fR_z(v_all)+F_rez+m*g*np.sin(alpha)))/L
    Z1 = m*g*np.cos(alpha) - Z2

    # realne sile na kolesih, 4x4!
    F_K2 = F_K/(1+Z1/Z2)
    F_K1 = F_K-F_K2

    # max sili na kolesih za 4x4 pogon
    def F_K1_max(alpha):
        return mu*np.cos(alpha)*1e-3 * (m*g*np.cos(alpha)*L1 + h_t*(fR_z(v_all)+F_rez+m*g*np.sin(alpha)))/L
    def F_K2_max(alpha):
        return mu*np.cos(alpha)*1e-3 * (m*g*np.cos(alpha) - (m*g*np.cos(alpha)*L1 + h_t*(fR_z(v_all)+F_rez+m*g*np.sin(alpha)))/L)

    # mejni klanec
    alpha_m = (D+R_f/m/g/1e-3-f*np.sqrt(1-(D+R_f/m/g/1e-3)**2+f**2))/(1+f**2)

    # max hitost je kjer se sekata F_K in R
    a=D*g/np.transpose([delta_f(i_ratios)]) #pospešek
    a = np.where(a < 0, 0, a) #zamenjaj negativne z nič
    R_i = fR_i_rot(a, np.transpose([i_ratios]))
    R = R_s+R_f+fR_z(v_all)+R_i #skupni upori - upošteva se R_i rotacijskih mas!
    v_max = v_all.ravel()[np.abs(F_K-R).argmin()]
    return F_K1, F_K2, v_max, R, alpha_m

#na ravnini!!!
F_K1, F_K2, v_max, R, alpha_m = sile(alpha, F_K_max)
F_K = F_K1+F_K2
F_K1_dej=np.where(F_K1 > F_K1_max(alpha), F_K1_max(alpha), F_K1)
F_K2_dej=np.where(F_K2 > F_K2_max(alpha), F_K2_max(alpha), F_K2)

# mejni klanci
print(np.max(np.tan(alpha_m)*100.))
fig,ax=plt.subplots()
for i in range(len(i_ratios)):
    ax.plot(v_all[i], 100*np.tan(alpha_m[i]), c='k')
    ax.text(v_all[i,-5]+3, 100*np.tan(alpha_m[i][-5])+0.02, roman[i])
# plt.plot(v_all.T, 100*np.tan(alpha_m).T, 'k-')
plt.xlabel('$v$ [km/h]')
plt.ylabel('$\\alpha$ [\%]')
plt.show()

# sile in upori
plt.axhline(F_tr_4x4, c='k', ls='--')
plt.plot(v_all.T, F_K.T, 'k-', label='$F_K$')
plt.plot(v_all.T, F_K1.T, 'b-', label='prednje')
plt.plot(v_all.T, F_K2.T, 'g-', label='zadnje')
plt.plot(v_all.T, F_K1_dej.T, 'b--', label='prednje')
plt.plot(v_all.T, F_K2_dej.T, 'g--', label='zadnje')
plt.plot(v_all.T, R.T, 'r-')

plt.xlabel('$v$ [km/h]')
plt.ylabel('$F_K$ [kN]')
custom_lines = [
    Line2D([0], [0], color='k', ls='--', lw=1, label='$F_{K, max}$'),
    Line2D([0], [0], color='k', lw=1, label='$F_K$'), 
    Line2D([0], [0], color='b', lw=1, label='prednje kolo'), 
    Line2D([0], [0], color='g', lw=1, label='zadnje kolo'),
    Line2D([0], [0], color='r', ls='-', lw=1, label='$\\sum R$')]
plt.legend(handles=custom_lines)
plt.show()

# razmerje pogonskih sil
plt.plot(v_all.T, np.transpose((F_K1/F_K2)), 'k-')
plt.xlabel('$v$ [km/h]')
plt.ylabel('razmerje sile prednje/zadnje kolo')
plt.show()

# odvisnost v_max in F_K od klanca
##############################################################
nakloni_percent=np.linspace(0, 100, 250) #do 100% klanec
nakloni=np.arctan(nakloni_percent/100.) #v rad
F_K1=[]
F_K2=[]
v_max=[]
for alpha in nakloni:
    F_K1i, F_K2i, v_maxi, Ri, alpha_mi = sile(alpha, F_K_max=F_K_max)
    F_K1.append(F_K1i)
    F_K2.append(F_K2i)
    v_max.append(v_maxi)
F_K1=np.asarray(F_K1)
F_K2=np.asarray(F_K2)

# 1. v_max
plt.plot(nakloni_percent, v_max, 'kx-')
plt.xlabel('$\\alpha$ [\%]')
plt.ylabel('$v_{max}$ [km/h]')
plt.show()

# 2. F_K
fig,(ax1,ax2,ax3)=plt.subplots(1, 3, gridspec_kw={'width_ratios': [1, 1, 0.1]}, figsize=(10, 5))
cmap=plt.get_cmap('turbo')

ax1.axhline(F_K_max, c='k')
ax2.axhline(F_K_max, c='k')
for i in range(F_K1.shape[0]):
    ax1.plot(v_all.T, F_K1[i].T, '-', c=cmap(nakloni_percent[i]/100.), alpha=0.5)
    ax2.plot(v_all.T, F_K2[i].T, '-', c=cmap(nakloni_percent[i]/100.), alpha=0.5)

ax1.set_xlabel('$v$ [km/h]')
ax1.set_ylabel('$F_K$ [kN]')
ax2.set_xlabel('$v$ [km/h]')
ax2.set_ylabel('$F_K$ [kN]')
ax1.set_title('prednje kolo')
ax2.set_title('zadnje kolo')
ax1.set_ylim(0,1.1*max([F_K1.max(),F_K2.max(), F_K_max]))
ax2.set_ylim(0,1.1*max([F_K1.max(),F_K2.max(), F_K_max]))

sm = plt.cm.ScalarMappable(cmap=cmap)
cbar = fig.colorbar(sm, cax=ax3, orientation='vertical', label='naklon [\%]')
percent_slopes = np.asarray([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
cbar.set_ticks(percent_slopes/100.)
cbar.set_ticklabels([f"{p}%" for p in percent_slopes])
plt.tight_layout()
plt.show()

# 3. F_K razmerje
fig,ax=plt.subplots(1, 1)
cmap=plt.get_cmap('turbo')
for i in range(F_K1.shape[0]):
    ax.plot(v_all.T, np.transpose((F_K1[i]/F_K2[i])), c=cmap(nakloni_percent[i]/100.), lw=5)
ax.set_xlabel('$v$ [km/h]')
ax.set_ylabel('razmerje sile prednje/zadnje kolo')
sm = plt.cm.ScalarMappable(cmap=cmap)
cbar = fig.colorbar(sm, ax=ax, orientation='vertical', label='naklon [\%]')
cbar.set_ticks(percent_slopes/100.)
cbar.set_ticklabels([f"{p}%" for p in percent_slopes])
plt.tight_layout()
plt.show()