#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  6 13:01:06 2020

@author: boyaolong
"""

import gurobipy as gp
import numpy as np
import cmath
f0=open("ieee33节点配电系统.txt");
#33节点信息
t=33;
l=32;
#前推回代的潮流
#def set_information():
#预置电压和功率
PQ=np.zeros((t,2));
V=np.ones(t,complex)
RX=np.zeros((l,4));
UB=12.66;#电压基准
SB=10;#功率基准
ZB=UB*UB/SB;#阻抗基准
A0=np.zeros((t,t),int);
for line in f0.readlines():
	if line[0]=='#':
		continue;
	if line[:2]=="-1":
		break;
PQ[int(line.split(' ')[0])][0]=float(line.split(' ')[5])/SB/1000;#有功
PQ[int(line.split(' ')[0])][1]=float(line.split(' ')[6])/SB/1000;#无功
RX[int(line.split(' ')[0])-1][0]=int(line.split(' ')[1])-1;
RX[int(line.split(' ')[0])-1][1]=int(line.split(' ')[2])-1;
RX[int(line.split(' ')[0])-1][2]=float(line.split(' ')[3])/ZB;
RX[int(line.split(' ')[0])-1][3]=float(line.split(' ')[4])/ZB;
A0[int(line.split(' ')[1])-1][int(line.split(' ')[2])-1]=1;
A0T=np.transpose(A0);
#print(A0T)
S=np.zeros(t,complex);
ZL=np.zeros(t,complex);
IL=np.zeros(t,complex);
ZL[0]=0;
for i in range(t):
	S[i]=complex(-PQ[i][0],-PQ[i][1]);#复数功率
for i in range(l):
	ZL[i+1]=complex(RX[i][2],RX[i][3]);
#print(ZL)
V[0]=1
IL[t-1]=-np.conjugate(S[t-1]/V[t-1]);#支路电流
max_error=1;#迭代误差
TempV=V;
Vangle=np.zeros((t,2))
#潮流计算
#print(FT)
k=0;
while max_error>0.0001:
	k+=1;
	IN=np.conjugate(S/V)#节点注入电流
	for i in range(l):
		#print(A0[l-i,l-i:])
		#print(IL[l-i+1:])
		#print(A0[l-i-1,l-i:]@IL[l-i:]-IN[l-i-1])
		IL[l-i-1]=A0[l-i-1,l-i:]@IL[l-i:]-IN[l-i-1];#python的矩阵乘法要换成@号
	for j in range(1,t):
		#电压前推过程
		#print(A0T[i,:i]*V[:i]-ZL[i]*IL[i]);
		V[j]=A0T[j,:j]@V[:j]-ZL[j]*IL[j];
	max_error=max(abs(V-TempV));
	TempV = V; # 记忆迭代结果
Vangle[:,0]=abs(V)
for i in range(t):
	Vangle[i,1]=cmath.phase(V[i])/3.1415*180;
print(Vangle)#节点电压和相角


	

