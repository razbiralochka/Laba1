import numpy as np
import math
import matplotlib.pyplot as plt
import random
import scipy
T = []
Phi_list = []

def simulation(_Mct,_Rct,_H0,_m,_l,_r,_phi,_rand):


   # list_Kg=[80*0.001,160*0.001,400*0.001,4]
    Kg=3
    R0=6371

    H0=_H0
    Mct=_Mct
    Rct=_Rct
    m=_m*(1+random.randint(-_rand,_rand)/100)
    l=_l*(1+random.randint(-_rand,_rand)/100)
    r=_r*(1+random.randint(-_rand,_rand)/100)
    phi=_phi*(1+random.randint(-_rand,_rand)/100)

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









    def f1(t, _phi, _omega):
        return _omega

    def f2(t, _phi, _omega):
        res = -C * _phi - K * _omega
        res = res *(1+random.randint(-_rand,_rand)/100)
        return res


    k = np.zeros(4)
    q = np.zeros(4)





    phi=phi*math.pi/180
    omega=0.0
    h=1
    t=0



    while t < pow(2,17):
        Phi_list.append(phi*180/math.pi)
        T.append(t)
        K =K_theor * (1 + random.randint(-_rand, _rand) / 100)
        k[0] = h * f1(t,phi,omega)
        q[0] = h * f2(t,phi,omega)
        k[1] = h * f1(t+h/2,phi+k[0]/2,omega+q[0]/2)
        q[1] = h * f2(t+h/2,phi+k[0]/2,omega+q[0]/2)
        k[2] = h * f1(t+h/2,phi+k[1]/2,omega+q[1]/2)
        q[2] = h * f2(t+h/2,phi+k[1]/2,omega+q[1]/2)
        k[3] = h * f1(t+h,phi+k[2],omega+q[2])
        q[3] = h * f2(t+h,phi+k[2],omega+q[2])
        phi = phi + (k[0]+2*k[1]+2*k[2]+k[3])/6
        omega = omega + (q[0]+2*q[1]+2*q[2]+q[3])/6
        t = t + h







def main():

    simulation(500,0.7,450,12,35,0.05,10,10) #Данные м варианта + рандом в %
    plt.plot(T, Phi_list, label='случайные отклонения')
    plt.legend(loc=2)
    plt.show()

    N = len(Phi_list)

    y=np.array(Phi_list)
    yf=scipy.fft.fft(y)/(N/2)
    xf = scipy.fft.fftfreq(N, 1/N)/N
    plt.plot(xf,np.abs(yf), label='АЧХ с помехами')
    plt.legend(loc=2)
    plt.xlim([0,0.001])
    plt.show()



    T.clear()
    Phi_list.clear()
    simulation(500, 0.7, 450, 12, 35, 0.05, 10, 0)
    plt.plot(T, Phi_list,label='без помех')
    plt.legend(loc=2)
    plt.show()

    y = np.array(Phi_list)
    yf = scipy.fft.fft(y) / (N / 2)
    xf = scipy.fft.fftfreq(N, 1 / N) / N
    plt.plot(xf, np.abs(yf), label='АЧХ без помех')
    plt.legend(loc=2)
    plt.xlim([0, 0.001])
    plt.show()


if __name__ == "__main__":
	main()
