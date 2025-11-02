#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Sestavljena gonila v mobilni tehniki
1. poročilo: Dinamična karakteristika vozila
Branje podatkov vozila iz datoteke
"""

import numpy as np

def ReadData(filename):
    with open(filename, mode='rt', encoding='ISO-8859-1') as dat:
        a=dat.readline().split(',')
        r_st=float(a[1])

        a=dat.readline().split(',')
        m_k=float(a[1])

        a=dat.readline().split(',')
        f=float(a[1])

        a=dat.readline().split(',')
        C_d=float(a[1])

        a=dat.readline().split(',')
        A=float(a[1])

        a=dat.readline().split(',')
        L=float(a[1])

        a=dat.readline().split(',')
        m=float(a[1])

        a=dat.readline().split(',')
        eta=float(a[1])

        a=dat.readline().split(',')
        i_final=float(a[1])

        a=dat.readline().split(',')
        i_ratios=[]
        for i in range(1,len(a)):
            i_ratios.append(float(a[i]))
        i_ratios = np.array(i_ratios)

        n=[]
        M=[]
        P_kW=[]
        dat.readline()
        for line in dat:
            line=line.split(',')
            n.append(float(line[0]))
            M.append(float(line[1]))
            P_kW.append(float(line[2]))
        n=np.array(n)
        M=np.array(M)
        P_kW=np.array(P_kW)
    return r_st,m_k,f,C_d,A,L,m,eta,i_final,i_ratios,n,M,P_kW