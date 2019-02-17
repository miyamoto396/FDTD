# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 13:38:16 2019

@author: yuta

伝送線路のクラス
コンストラクタの引数
Length   伝送線路の長さ
N        分割数
R        抵抗値
Z        電位係数と誘導係数から計算される特性インピーダンス
mu       誘電率
epsilon  伝導率

"""
import numpy as np
import sympy 
import matplotlib.pyplot as plt
import datetime


sympy.init_printing()

c = 2.99793458e8    #光速

class TransmissionLine:
    def __init__(self,Length,N,R,Z):
        self.Length = Length
        self.N = N
        self.R = R
        self.Z = Z
        self.dx = Length/N
        self.dt = self.dx / c
        self.V = np.array([0]*N,dtype=float)        #電圧のメッシュがN個
        self.I = np.array([0]*(N-1),dtype=float)    #電流のメッシュがN-1個、電圧のメッシュの内側に入るように考える
        self.x_v = np.linspace(0,Length,N)      #vプロット用横軸
        self.x_i = np.linspace(0,Length,N-1)    #iプロット用横軸
        self.t = 0.0
        self.I0 = 0.0    #電源側の境界計算用電流値
        self.In = 0.0

    def Load(self,Rs,Rl):   #Rs:電源側抵抗,Rl:負荷側抵抗
        self.Rs = Rs
        self.Rl = Rl
        self.CoeffA =self.Z/c/self.dt + self.R/2 
        self.CoeffB =self.Z/c/self.dt - self.R/2 

    def calcV(self,E):
        Rs = self.Rs
        Rl = self.Rl
        dx = self.dx
        dt = self.dt
        Z = self.Z
        N = self.N
        
        #電源側での境界条件を解く
        V0sy,I0sy =sympy.symbols('V0sy I0sy')
        ans0 = sympy.solve([ (V0sy - self.V[0])/c/dt + Z*(2*self.I[0]-I0sy-self.I0)/dx, V0sy-E+Rs*I0sy ],[V0sy,I0sy])
        self.V[0] = ans0[V0sy]
        self.I0 = ans0[I0sy]
        
        #負荷側の虚位迂回条件を解く
        Vnsy,Insy =sympy.symbols('Vnsy Insy')
        ansn = sympy.solve([ (Vnsy - self.V[N-1])/c/dt + Z*( Insy+self.In-2*self.I[N-2] )/dx, Vnsy-Rl*Insy ],[Vnsy,Insy])
        self.V[N-1] = ansn[Vnsy]
        self.In = ansn[Insy]
        
        #伝送線路の電圧を計算
        for i in range(1,self.N-1):
            self.V[i]= self.V[i]-Z*(self.I[i]-self.I[i-1])
            
    def calcI(self):
        dx = self.dx        

#        Z = self.Z
#        dt = self.dt
#        R = self.R        
#        ii = sympy.Symbol('ii')

        for i in range(0,self.N-1):
          self.I[i] = self.I[i]*self.CoeffB/self.CoeffA - (self.V[i+1]-self.V[i])/dx/self.CoeffA
          #ans_temp = sympy.solve([ (self.V[i+1]-self.V[i])/dx-Z*(ii-self.I[i])/c/dt-R*(ii+self.I[i])/2],[ii]) 
          #self.I[i] = ans_temp[ii]
            
    def CalcLastTime(self,time):
        LastTimeNumItre = int(self.N/self.Length*c*time)
        print('最終計算時間は:{0}'.format(time))
        print('最終計算時間までの計算回数:{0}'.format(LastTimeNumItre) )

        StartTime = datetime.datetime.now()         #計算開始時刻を取得
        print('計算開始時刻{0}'.format(StartTime))

        self.V_pos = np.array([0.0]*LastTimeNumItre,dtype=float)    #位置の電圧を保存したい変数
        self.I_pos = np.array([0.0]*LastTimeNumItre,dtype=float)    #位置の電圧を保存したい変数        
        self.t_pos = np.array([0.0]*LastTimeNumItre,dtype=float)    #時間を保存しておく
        
        for i in range(0, LastTimeNumItre):
            self.V_pos[i] = self.V[self.N-1]        #ある位置の電位を保存しておく
            self.I_pos[i] = self.I[self.N-2]
            self.t_pos[i] = self.t                  #横軸プロット用時間
            
            self.calcV(30)  #引数が電源電圧
            self.calcI()
            self.t = self.t + self.dt
        
            #計算の実行経過を出力
            if i == int(LastTimeNumItre/10):        #計算が10%終わったときの計算時刻を出力
                Time1 = datetime.datetime.now()
                td = Time1-StartTime
                print('計算が10%終了')
                print('経過時間{0}'.format(td))
            elif i == int(LastTimeNumItre/2):       #計算が50%終わったときの計算時刻を出力
                Time2 = datetime.datetime.now()
                td = Time2-StartTime
                print('計算が50%終了')
                print('経過時間{0}'.format(td))
            elif i == int(LastTimeNumItre/4*3):     #計算が80%終わったときの計算時刻を出力
                Time3 = datetime.datetime.now()
                td = Time3-StartTime
                print('計算が80%終了')
                print('経過時間{0}'.format(td))

    def Printinit(self):
        plt.subplots_adjust(wspace=1.0, hspace=1.5)
        
    def PrintV(self):  
        plt.figure(figsize=(5,10))
        plt.subplot(4, 1, 1)
        plt.ylim([0,40])
        plt.plot(self.x_v,self.V)


        plt.title("V-m at LastTime")
        plt.xlabel("position[m]")
        plt.ylabel("Voltage[V]")
        plt.show()

    def PrintI(self):
        plt.figure(figsize=(5,10))
        plt.subplot(4, 1, 2)
        #plt.ylim(0,)
        plt.plot(self.x_i,self.I)

        plt.title("A-m at LastTime")
        plt.xlabel("position[m]")
        plt.ylabel("Current[A]")
        plt.show()

    def PrintVpos(self):
        plt.figure(figsize=(5,10))
        plt.subplot(4, 1, 3)
        plt.plot(self.t_pos,self.V_pos)

        plt.title("V-s at")
        plt.xlabel("time[s]")
        plt.ylabel("Voltage[V]")
        plt.show()
    
    def PrintIpos(self):
        plt.figure(figsize=(5,10))
        plt.subplot(4, 1, 4)
        plt.plot(self.t_pos,self.I_pos)

        plt.title("A-s at")
        plt.xlabel("time[s]")
        plt.ylabel("Current[A]")
        plt.show()
        
#TransmissionLineクラスここまで

Length1 = 10
N = 20
R = 1.0e-6
Z = 10
Rs1 = 0
Rl1 = 50
LastTime = 1.0e-6

TL1 = TransmissionLine(Length1, N , R , Z)      #伝送線路のインスタンスを作成       
TL1.Load(Rs1, Rl1)                              #負荷の情報をイ入力
TL1.CalcLastTime(LastTime)                       #計算

TL1.Printinit()
TL1.PrintV()
TL1.PrintI()
TL1.PrintVpos()
TL1.PrintIpos()


