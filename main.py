import numpy as np
import math
import matplotlib.pyplot as plt
import random

gTime1=list()
gPhi1=list()
gPhi2=list()
gTime2=list()
def simulation(_Mct,_Rct,_H0,_m,_l,_r,_phi,_rand,mas1,mas2):


   # list_Kg=[80*0.001,160*0.001,400*0.001,4]
    Kg=1
    R0=6371

    H0=_H0
    Mct=_Mct
    Rct=_Rct
    m=_m*(1+random.randint(-_rand,_rand)/100)
    l=_l*(1+random.randint(-_rand,_rand)/100)
    r=_r*(1+random.randint(-_rand,_rand)/100)
    phi=_phi#*(1+random.randint(-_rand,_rand)/100)

    Mu=3.9858*pow(10, 5)


    J0=0.4*Mct*pow(Rct, 2)
    Jm=0.4*m*pow(r, 2)
    Jz=J0+2*m*l*l+2*Jm
    Jy=J0+2*Jm
    Jx=0.9*Jz





    R=R0+H0
    K_theor=Kg/Jz
    C=(3*Mu)/pow(R, 3)
    C=C*(Jx-Jy)/Jz




    def dphi(t, _phi, _omega):
        return _omega

    def domega(t, _phi, _omega):
        res = -C * _phi - K * _omega
        res = res *(1+random.randint(-_rand,_rand)/100)
        return res


    k = np.zeros(4)
    q = np.zeros(4)





    phi=phi*math.pi/180
    omega=0.0
    h=1
    t=0

    T = []
    Phi_list = []
    omega_list =[]

    while t < 3600*24*2:
        Phi_list.append(phi*180/math.pi)
        omega_list.append(omega*180/math.pi)
        mas1.append(phi*180/math.pi)
        mas2.append(t)
        T.append(t)
        K =K_theor * (1 + random.randint(-_rand, _rand) / 100)
        k[0] = h * dphi(t,phi,omega)
        q[0] = h * domega(t,phi,omega)
        k[1] = h * dphi(t+h/2,phi+k[0]/2,omega+q[0]/2)
        q[1] = h * domega(t+h/2,phi+k[0]/2,omega+q[0]/2)
        k[2] = h * dphi(t+h/2,phi+k[1]/2,omega+q[1]/2)
        q[2] = h * domega(t+h/2,phi+k[1]/2,omega+q[1]/2)
        k[3] = h * dphi(t+h,phi+k[2],omega+q[2])
        q[3] = h * domega(t+h,phi+k[2],omega+q[2])
        phi = phi + (k[0]+2*k[1]+2*k[2]+k[3])/6
        omega = omega + (q[0]+2*q[1]+2*q[2]+q[3])/6
        t = t + h

         
    N = len(Phi_list)

    y=np.array(Phi_list)
    yf=np.fft.fft(y)/N
    yf=np.abs(yf)*2
    xf=np.linspace(0,(N-1),N)


    xf=xf[0:int(N/1500)]*pow(10,5)/N
    yf = yf[0:int(N/1500)]
    
    
    
    
    fig, [ax1,ax2,ax3] = plt.subplots(nrows = 3, ncols = 1 )


    ax1.plot(T, Phi_list, label='Возмущения: ' + str(_rand) + '%')
    ax2.plot(xf, yf, '.-',label='АЧХ | Возмущения: ' + str(_rand) + '%')
    ax3.plot(Phi_list, omega_list,label='Фазовый портерт | Возмущения: ' + str(_rand) + '%')
    ax1.grid()
    ax1.legend(loc=2)
    ax1.set_xlabel('время (с)')
    ax1.set_ylabel('Угол (град.)')

    ax2.grid()
    ax2.legend(loc=2)

    ax2.set_xlabel('Частота (герц * 10^(-5))')
    ax2.set_ylabel('Амплитута (град.)')

    ax3.grid()
    ax3.legend(loc=2)
    ax3.set_xlabel('Угол (град.)')
    ax3.set_ylabel('Угловая скорость (град. / c)')
    plt.show()
    

def Asim(_Mct,_Rct,_H0,_m,_l,_r,_phi,):


   # list_Kg=[80*0.001,160*0.001,400*0.001,4]
    Kg=4
    R0=6371

    H0=_H0
    Mct=_Mct
    Rct=_Rct
    m=_m
    l=_l
    r=_r
    phi=_phi

    Mu=3.9858*pow(10, 5)


    J0=0.4*Mct*pow(Rct, 2)
    Jm=0.4*m*pow(r, 2)
    Jz=J0+2*m*l*l+2*Jm
    Jy=J0+2*Jm
    Jx=0.9*Jz





    R=R0+H0
    
    C=(3*Mu)/pow(R, 3)
    C=C*(Jx-Jy)/Jz
    K=Kg/Jz

    phi=phi*math.pi/180

  





    Amp=phi
    
    h=1
    t=0

    T = []
    Phi_list = []

    beta = K/2
    w = math.sqrt(C)
    freq = w/(2*math.pi)
    while t < 3600*72:
        Phi_list.append(phi*180/math.pi)
        T.append(t)
        phi=Amp*math.exp(-beta*t)*math.cos(w*t)
        t = t + h

    fig, ax = plt.subplots()       
    ax.plot(T, Phi_list, label='Аналитическое решение ' )
    ax.grid()
    ax.legend(loc=2)
    ax.set_xlabel('время (с)')
    ax.set_ylabel('Угол (град.)')
    print(freq*(10**5),' Герц 10^-5')
    print(freq, ' Герц')
    plt.show()

     


def main():

    Mct = 400
    Rct=0.6
    H0 = 300
    m=14
    l=20
    r=0.1
    phi=10
    t=0
    simulation(Mct,Rct,H0,m,l,r, phi,50,gPhi1,gTime1)
    simulation(Mct,Rct,H0,m,l,r, phi,0, gPhi2,gTime2)


    plt.plot(gTime1,gPhi1)
    plt.plot(gTime2, gPhi2)
    plt.grid()
    plt.show()

main()
