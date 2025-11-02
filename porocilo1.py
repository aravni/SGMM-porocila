#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Sestavljena gonila v mobilni tehniki
1. poročilo: Dinamična karakteristika vozila
Splošni grafi analize dinamične karakteristike vozila
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_simpson

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

# konstante okolja
alpha=np.deg2rad(0) #naklon, rad
rho=1.225 #gostota zraka, kg/m^3
g=9.81 #gravitacija, m/s^2

#####################################################
# grafi
#####################################################

############ 1. karakteristika motorja
fig,ax=plt.subplots()
ax.plot(n,P_kW,'r-x')
ax.set_xlabel('$n$ [1/min]')
ax.set_ylabel('$P$ [kW]', color='r')
ax.minorticks_on()
ax.tick_params(axis='y', colors='r', which='both')

ax1=plt.twinx(ax)
ax1.plot(n,M,'k-x')
ax1.set_ylabel('$M$ [Nm]', color='k')
ax1.minorticks_on()
ax1.tick_params(axis='y', colors='k', which='both')

alignYaxes([ax,ax1], [0,0])
plt.tight_layout()
plt.show()

############ 2. žagasti diagram
n_Pmax = n[np.argmax(P_kW)]
n_Mmax = n[np.argmax(M)]

# linije prestav
v_all = [] #seznam hitr. po prestavah iz meritev obratov
roman=['I','II','III','IV','V','VI','VII','VIII','IX','X']
fig,ax = plt.subplots()
for i_n,i in enumerate(i_ratios):
    vi = 2*np.pi*n*60 / (i*i_final) *r_din/1e6 #km/h
    v_all.append(vi)
    ax.plot(n,vi, c='grey', alpha=0.7)
    ax.text(n[-1]+100, vi[-1], roman[i_n])
v_all=np.array(v_all)
ax.axvline(n_Pmax, c='grey', alpha=0.7, ls='--')
ax.axvline(n_Mmax, c='grey', alpha=0.7, ls='--')
ax.text(n_Pmax-200, np.max(v_all)*0.96, '$n_{P_{max}}$', rotation='vertical')
ax.text(n_Mmax-200, np.max(v_all)*0.95, '$n_{M_{max}}$', rotation='vertical')

# speljevanje
n_data=[] #redkejši podatki od meritev (n in v_all), saj je potek linearen
v_data = []
n_start = np.linspace(0,n_Mmax,2)
v_start = 2*np.pi*n_start*60 / (i_ratios[0]*i_final) *r_din/1e6 #km/h
n_data.append(n_start)
v_data.append(v_start)

# menjave prestav
n_min = n_Mmax
n_max = n_Pmax
for i in range(len(i_ratios)):
    ni = np.linspace(n_min,n_max,2)
    vi = 2*np.pi*ni*60 / (i_ratios[i]*i_final) *r_din/1e6 #km/h
    n_data.append(ni)
    v_data.append(vi)
    #TODO časovni zamik pri menjavi prestav
    if i < len(i_ratios)-1: n_min = n_max * i_ratios[i+1]/i_ratios[i] #neskončno hitra menjava prestave
    if i == len(i_ratios)-2: n_max=max(n) #končna prestava
n_data=np.array(n_data)
v_data=np.array(v_data)

ax.plot(n_data.ravel(), v_data.ravel(), 'k-')
v_max_teo = v_data.max() #brez uporov
ax.axhline(v_max_teo, c='gray', ls='--')
ax.text(100, v_max_teo-12, f'$v_{{max}}={v_max_teo:.3g}$ km/h')
ax.set_xlim(0)
ax.set_ylim(0)
ax.set_xlabel('$n$ [1/min]')
ax.set_ylabel('$v$ [km/h]')
plt.minorticks_on()
plt.tight_layout()
plt.show()

############ 3. upori
R_f = f*m*g*np.cos(alpha)/1e3 #kotalni, kN
R_s = m*g*np.sin(alpha)/1e3 #upor strmine, kN
fR_z = lambda v: 0.5*A/1e6*C_d*rho*(v/3.6)**2/1e3 #zračni upor, kN -> param. v: hitrost [km/h]
v_range=np.linspace(0,np.max(v_all),len(M)) #to niso meritve!
R_z = fR_z(v_range)

# upor vztr. mas
i_m=np.array(i_ratios[:]) #vse prestave v menjalniku -> še ne vemo pospeškov -> rabimo za izračun kasneje

# k=0.04 #osebni avto
# delta = 1.03+k*i_m**2 #1.način (grobo)

I_k = 0.5*m_k*(r_st)**2 #MVM kolesa
k1 = 4*I_k*g*1e-3/(r_st*r_din) #kg km/s^2 = kN
k2=0.007 #ocena
delta_f = lambda i: 1+k1+k2*i**2 #2.način (bolje); i: prestavno razmerje
# k1=0.076 #ocena
# delta = 1+k1+k2*i_m**2 #2.način (bolje)

fR_i_tr = lambda a: m*a*1e-3 #upor translatornih vztr. mas, kN; a: pospešek [m/s^2]
fR_i_rot = lambda a,i: m*a*(delta_f(i)-1.)*1e-3 #upor rotirajočih vztr. mas, kN; a: pospešek [m/s^2]
fR_i = lambda a,i: m*a*delta_f(i)*1e-3 #skupni upor vztr. mas, kN; a: pospešek [m/s^2]
R = R_s+R_f+R_z #skupni upor, brez vztr. mas

# plt.fill_between(v_range, R_s, label='$R_s$')
# plt.fill_between(v_range, R_s, R_f+R_s, label='$R_f$')
# # plt.fill_between(v_range, R_f+R_s, R_f+R_s+R_i, label='$R_i$')
# plt.fill_between(v_range, R_f+R_s, R, label='$R_z$')
# plt.plot(v_range, R, 'k-', lw=3, label='$\sum R$')
# plt.xlabel('$v$ [km/h]')
# plt.ylabel('$R$ [kN]')
# plt.legend(loc='upper left')
# plt.minorticks_on()
# plt.tight_layout()
# plt.show()

############ 4. diagram sile
F_K=[]
n1=0
for i in range(len(i_ratios)):
    FKi= M*i_ratios[i]*i_final*eta/r_din
    F_K.append(FKi)
F_K=np.array(F_K)
F_ideal = max(P_kW)/v_range*3.6 #idealna sila, kN

fig,ax=plt.subplots()
ax.plot([],[],'k-',label='$F_K$')
ax.plot(v_range, F_ideal, '--', c='gray', alpha=0.7, label='$F_{ideal}$')
ax.annotate('$F_{ideal}$',
            xy=(v_range[15], F_ideal[15]), xycoords='data',
            xytext=(v_range[15]+50, F_ideal[15]+1), textcoords='data',
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3"))
ax.plot(v_range, R, 'r-', label='$\sum R$')
ax.text(v_range[30]+20, R[30]-0.3, '$\sum R$', c='red')
for i in range(len(v_all)):
    ax.plot(v_all[i], F_K[i], 'k-')
    ax.text(v_all[i][-10]+5, F_K[i][-10]+0.2, roman[i])

v_max_F = v_all.ravel()[np.abs(F_K-(R_s+R_f+fR_z(v_all))).argmin()] #max hitost je kjer se sekata F_K in R
ax.axvline(v_max_F, c='gray', ls='--')
# ax.axvline(v_max_teo, c='gray', ls='--')
ax.text(v_max_F-10, max(F_K[1]), f'$v_{{max}}={v_max_F:.3g}$ km/h', rotation='vertical')
ax.set_xlim(0)
ax.set_ylim(0,np.max(F_K)*1.1)
ax.set_xlabel('$v$ [km/h]')
ax.set_ylabel('$F$ [kN]')
plt.minorticks_on()
plt.tight_layout()
plt.show()

############ 5. dinamični vozni faktor
F_rez=[]
D=[]
fig,ax=plt.subplots()
for i in range(len((i_ratios))):
    Fi = F_K[i] - fR_z(v_all[i]) -R_f-R_s #R_i je upoštevan kasneje pri izračunu pospeškov!
    F_rez.append(Fi)
    D.append(Fi/(m*g*1e-3))
    ax.plot(v_all[i], D[i], c='k')
    if D[i][-8] < 0:
        ax.text(v_all[i, np.argmax(1/D[i])], 0.02, roman[i])
    else: ax.text(v_all[i,-8]+3, D[i][-8]+0.02, roman[i])
D=np.array(D)
F_rez=np.array(F_rez)
D_ideal=(F_ideal-fR_z(v_range))/(m*g*1e-3)
ax.plot(v_range, D_ideal, '--',c='gray')
ax.set_xlim(0)
ax.set_ylim(0,np.max(D)*1.1)
ax.set_xlabel('$v$ [km/h]')
ax.set_ylabel('$D$ [/]')
plt.minorticks_on()
plt.tight_layout()
plt.show()

############ 6. mejna strmina
alpha_m=[]
fig,ax=plt.subplots()
for i in range(len(i_ratios)):
    alpha_i = (D[i]+R_f/m/g/1e-3-f*np.sqrt(1-(D[i]+R_f/m/g/1e-3)**2+f**2))/(1+f**2)
    alpha_m.append(np.arcsin(alpha_i))
    ax.plot(v_all[i], np.rad2deg(alpha_m[i]), c='k')
    ax.text(v_all[i,-5]+3, np.rad2deg(alpha_m[i][-5])+0.02, roman[i])
    # ax.plot(v_all[i], np.tan(alpha_m[i])*100, c='k')
    # ax.text(v_all[i,-5]+3, np.tan(alpha_m[i][-5])*100+0.02, roman[i])
alpha_m=np.array(alpha_m)
alpha_ideal = (D_ideal-f*np.sqrt(1-D_ideal**2+f**2))/(1+f**2)
ax.plot(v_range, np.rad2deg(alpha_ideal), '--',c='gray')
# ax.plot(v_range, np.tan(alpha_ideal)*100, '--',c='gray')
ax.set_xlabel('$v$ [km/h]')
ax.set_ylabel('$\\alpha$ [°]')
# ax.set_ylabel('$\\alpha$ [\%]')
plt.minorticks_on()
plt.tight_layout()
plt.show()

############ 7. bilanca moči
P=[]
P_rez=[]
fig,ax=plt.subplots()
for i in range(len(i_ratios)):
    Pi = F_K[i]*v_all[i]/3.6
    P.append(Pi)
    Ri=R_f+R_s+fR_z(v_all[i])
    P_rez.append(Pi-Ri*v_all[i]/3.6)
    ax.plot(v_all[i],Pi, c='k')
    ax.text(v_all[i,-5]+3, Pi[-5]+3, roman[i])
P=np.array(P)
P_rez=np.array(P_rez)
ax.axhline(np.max(P_kW), c='gray', ls='--')
ax.text(0, np.max(P_kW)*1.03, f'$P_{{max,motor}}={np.max(P_kW):.3g}$ kW')
ax.axhline(np.max(P), c='gray', ls='--')
ax.text(0, np.max(P)*1.03, f'${np.max(P):.3g}$ kW')
ax.plot(v_range, R*v_range/3.6, 'r')
ax.text(v_range[30]+10, R[30]*v_range[30]/3.6-10, '$\sum P_{R}$', c='red')
# ax.axvline(v_max_F, c='gray', ls='--')
plt.xlabel('$v$ [km/h]')
plt.ylabel('$P$ [kW]')
plt.minorticks_on()
plt.tight_layout()
plt.show()

fig,ax=plt.subplots()
for i in range(len(i_ratios)):
    ax.plot(v_all[i],P_rez[i], c='k')
    ax.text(v_all[i,-5]+3, P_rez[i,-5]+3, roman[i])
ax.plot(v_range, np.max(P_kW)-R*v_range/3.6, '--',c='gray')
ax.plot(v_range, np.max(P)-R*v_range/3.6, '--',c='gray')
# ax.axvline(v_max_F, c='gray', ls='--')
plt.xlabel('$v$ [km/h]')
plt.ylabel('$P_{rez}$ [kW]')
plt.minorticks_on()
plt.tight_layout()
plt.show()

############ 8. pospeški
# a = (D-f)*g #pospeški po prestavah -> en. pdf
a = D*g #R_f je že upoštevan v D!
# a = (F_K - fR_z(v_all))/(m*g) #en. na vajah -> NI OK
for i in range(a.shape[0]):
    a[i] /= delta_f(i_ratios[i]) #upoštevanje ocene vztr. rotirajočih mas

a_real=[] #krivulja največjih pospeškov
v_real=[]

a0 = np.max(a[0])*p #speljevanje - konst. pospešek
v_real.extend(np.linspace(0, v_all[0,np.argmin(np.abs(a[0]-a0))], 3))
a_real.extend(np.full_like(v_real, a0))
v_prestave = [v_real.copy()]
a_prestave = [a_real.copy()]

for i in range(0, v_all.shape[0]-1):
    v_i=[] #po eni prestavi
    a_i=[] #po eni prestavi
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
                v_i.append(v_all[i,j])
                a_i.append(a[i,j])
            else: #preverba da je to res a_max pri danem v, ki je pozitiven
                if a[i,j] > a[i+1,k] and a[i,j] >= 0:
                    v_real.append(v_all[i,j])
                    a_real.append(a[i,j])
                    v_i.append(v_all[i,j])
                    a_i.append(a[i,j])
                    # plt.plot(v_all[i,:], a[i,:],'C0-')
                    # plt.plot(v_all[i+1,:], a[i+1,:],'C1-')
                    # plt.plot(v_all[i,j],a[i,j],'ro')
                    # plt.plot(v_all[i+1,k],a[i+1,k],'yo')
                    # plt.axvline(v_old,c='gray',ls='--')
                    # plt.axvline(v_new,c='k',ls='--')
                    # plt.axvline(v_all[i+1,k],c='y',ls='--')
                    # plt.title(f'{k}')
                    # plt.show()
    if i == 0:
        v_prestave[0].extend(v_i)
        a_prestave[0].extend(a_i)
    else:
        v_prestave.append(v_i)
        a_prestave.append(a_i)
#zadnja prestava do konca le, če je pospešek pozitiven
k=len(a[-1, v_all[-1] <= v_real[-1]])
v_prestave.append(v_all[-1,k:])
a_prestave.append(a[-1,k:])
for j in range(1,len(a[-1, v_all[-1] > v_real[-1]])):
    if a[-1,k] >= 0:
        v_real.append(v_all[-1,k])
        a_real.append(a[-1,k])
    k+=1
v_real=np.array(v_real)
a_real=np.array(a_real)

fig,ax=plt.subplots()
for i in range(len(i_ratios)):
    ax.plot(v_all[i],a[i], 'k-')
    ax.text(v_all[i,-5]+2, a[i,-5]+0.2, roman[i])
ax.plot(v_real, a_real,'r--')
ax.axvline(v_real[-1], c='gray', ls='--')
ax.text(v_real[-1]-10, a_real[-1]+2, f'$v_{{max}} = {v_real[-1]:.3g}$ km/h', rotation='vertical')
plt.xlabel('$v$ [km/h]')
plt.ylabel('$a$ [m/s$^2$]') 
plt.minorticks_on()
plt.tight_layout()
plt.show()

############ 9. čas in pot pospeševanja
t = cumulative_simpson(1/a_real, x=v_real/3.6, initial=0) #čas pospeševanja, s
t_uporabno = t[:np.argmax(t)] #le naraščajoči del
s = cumulative_simpson(v_real[:len(t_uporabno)]/3.6, x=t_uporabno, initial=0)*1e-3 #pot pospeševanja, km

fig,ax=plt.subplots()
ax.plot(v_real, t, 'r-')
ax.set_xlabel('$v$ [km/h]')
ax.set_ylabel('$t$ [s]', color='r')
ax.minorticks_on()
ax.tick_params(axis='y', colors='r', which='both')

ax1=plt.twinx(ax)
ax1.plot(v_real[:len(t_uporabno)], s, 'k-')
ax1.set_ylabel('$s$ [km]', color='k')
ax1.minorticks_on()
ax1.tick_params(axis='y', colors='k', which='both')
ax1.grid(False)

alignYaxes([ax,ax1], [0,0])
plt.tight_layout()

# Zoomed-in plot
v_limit = v_real[np.argmin((v_real-100)**2)] #najbližje 100 km/h
n_plot=len(v_real[v_real < v_limit])
if n_plot>len(t_uporabno): n_plot=t_uporabno
axins = ax.inset_axes([0.1, 0.4, 0.5, 0.5],xlim=(0, v_real[n_plot]),
                      ylim=(0, t[n_plot]), xticklabels=[], yticklabels=[])

axins.plot(v_real[:n_plot], t[:n_plot], 'r-')
axins.minorticks_on()
axins.set_xticks([0,v_real[n_plot]], ['0',f'{v_real[n_plot]:.3g}'])

axins1 = axins.twinx()
axins1.plot(v_real[:n_plot], s[:n_plot], 'k-')
axins1.minorticks_on()

ax.indicate_inset_zoom(axins, edgecolor='k')
alignYaxes([axins,axins1], [0,0])
axins1.set_yticks([0,s[n_plot]], ['0', f'{s[n_plot]:.3g}'])
axins.set_yticks([0,t[n_plot]], ['0',f'{t[n_plot]:.3g}'])
axins.tick_params(axis='y', colors='r', which='both')
plt.show()

# sedaj lahko iz pospeškov izračunam vztrajnostne upore
# R_i = fR_i(a_prestave, np.transpose([i_ratios])) #nehomogena oblika!
a_prestave = [[max(0, x) for x in sublist] for sublist in a_prestave] #zamenjaj negativne z nič
R_i_rot = [np.asarray(fR_i_rot(np.asarray(a), i)) for a, i in zip(a_prestave, i_ratios)]
R_i_tr = [np.asarray(fR_i_tr(np.asarray(a))) for a in a_prestave]
R_z = [np.asarray(fR_z(np.asarray(v))) for v in v_prestave]
R_i_tr = np.asarray([item for sublist in R_i_tr for item in sublist])
R_i_rot = np.asarray([item for sublist in R_i_rot for item in sublist])
R_z = np.abs(np.asarray([item for sublist in R_z for item in sublist]))
R = R_s + R_f + R_i_tr + R_i_rot + R_z
v_prestave_c = np.asarray([item for sublist in v_prestave for item in sublist])
a_prestave_c = np.asarray([item for sublist in a_prestave for item in sublist])

plt.fill_between(v_prestave_c, R_s, color='C0', label='$R_s$')
plt.fill_between(v_prestave_c, R_s, R_f+R_s, color='C1', label='$R_f$')
plt.fill_between(v_prestave_c, R_f+R_s, R_f+R_s+R_i_tr, color='C2', label='$R_{i,tr}$')
plt.fill_between(v_prestave_c, R_f+R_s+R_i_tr, R_f+R_s+R_i_tr+R_i_rot, color='C4', label='$R_{i,rot}$')
plt.fill_between(v_prestave_c, R_f+R_s+R_i_tr+R_i_rot, R, color='C3', label='$R_z$')
plt.plot([ ], [ ], 'C3-', lw=2, label='$\sum R$')
plt.plot(v_all[0], F_K[0], 'k-', lw=2, label='$F_K$')
plt.plot(v_all[1:].T, F_K[1:].T, 'k-', lw=2)
plt.xlabel('$v$ [km/h]')
plt.ylabel('$R$ [kN]')
plt.legend(loc='best')
plt.minorticks_on()
plt.tight_layout()
plt.show()