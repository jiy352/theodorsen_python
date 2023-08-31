import numpy as np
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import scipy.special as scsp

import matplotlib.pyplot as plt

from theodorson import theodorson_function, Lift, generate_kinemtics

# test generate_kinematics
a = -1/2 # quarter chord
# a = -1 # leading edge
b = 1
# V = np.array([10000, 10, 5, 1/10000])[:-1]
V = np.array([10])
omg = 1
alpha_max = np.deg2rad(5)
alpha_0 = np.deg2rad(5)
n_period = 2
t = np.linspace(0, n_period*2*np.pi/omg, 100)

# test Lift
rho = 1.0
k = omg * b / V

H1=scsp.hankel2(1,k)
H0=scsp.hankel2(0,k)

C=H1/(H1+1.j*H0)

L = 2*np.pi*rho*V**2*b*(C ) * alpha_max * (alpha_max * ( np.exp(1j*(np.pi/2-omg*t)) + np.deg2rad(alpha_0)))
L_NC = 2*np.pi*rho*V**2*b*(1j*2*omg*b/(V)) * (alpha_max * ( np.exp(1j*(np.pi/2-omg*t)) + np.deg2rad(alpha_0)))
Cl = L/(rho*V**2*b)

# Cl = (2*np.pi*C)*(alpha_max( np.exp(1j*(np.pi/2-omg*t)) )+ np.deg2rad(alpha_0))
# lgd = ['k=0','k=0.1', 'k=0.2','k=inf']
# plt.rcParams['font.size'] = 12
# plt.rcParams['axes.labelsize'] = 14
# plt.rcParams['axes.titlesize'] = 16
# plt.rcParams['lines.linewidth'] = 2
# plt.figure(figsize=(10, 6))
# for i in range(len(V)):
#     L_NC, L_C, CL = Lift(rho, V[i], a, b, alpha_max, omg, n_period, t)
#     T = 2*np.pi/omg
#     alpha = alpha_max * np.exp(1j*omg*t).real 
#     plt.plot(np.rad2deg(alpha), CL.real, label='Cl')
#     # plt.plot(t/T, CL.real, label='Cl')
#     # plt.plot(t, CL_k1.real/100, label='CL')
#     plt.xlabel('t/T ')
#     plt.xlabel('alpha ')
#     plt.ylabel('Cl')
#     plt.grid(True)
#     plt.legend([lgd[i]])
#     plt.xlim( -1,  1)
#     plt.ylim(-0.1, 0.1)
#     plt.savefig('CL'+str(i)+'.png', dpi=400, transparent=True)
# # plt.ylim(-1.7, 1.7)
# plt.show()
