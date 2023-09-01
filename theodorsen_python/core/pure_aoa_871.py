import numpy as np
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import scipy.special as scsp

import matplotlib.pyplot as plt

from theodorson import theodorson_function, Lift, generate_kinemtics

def lift(k):
    # test generate_kinematics
    a = -1/2 # quarter chord
    c = 2 # chord
    b = c * 0.5 # half chord
    rho = 1.0

    # define kinematics
    # k = 10000
    # k = 1e-10
    # k = 0.2
    omg = 1

    # compute freestream velocity from reduced frequency
    V = omg*b/k

    alpha_max = np.deg2rad(5)
    alpha_0 = np.deg2rad(5)
    n_period = 2
    t = np.linspace(0, n_period*2*np.pi/omg, 100)

    # profile of the angle of attack
    alpha = alpha_max * np.cos(omg*t) + alpha_0
    alpha_dot = -alpha_max * omg * np.sin(omg*t)
    alpha_dot_dot = -alpha_max * omg**2 * np.cos(omg*t)

    alpha = alpha_max*np.exp(1j*omg*t)+ alpha_0
    alpha_dot = -alpha_max * omg * np.sin(omg*t)
    alpha_dot_dot = -alpha_max * omg**2 * np.cos(omg*t)

    # compute theodorsen function
    H1=scsp.hankel2(1,k)
    H0=scsp.hankel2(0,k)
    C=H1/(H1+1.j*H0)
    F = C.real
    G = C.imag

    if k > 1e1:
        F = 0.5
        G = 0
    elif k == 1e-10:
        F = 1
        G = 0
    C = F + 1j*G

    # L_C = (2*np.pi*rho*V**2*b*(C) * alpha_max*np.exp(1j*omg*t)).real
    L_C = (2*np.pi*rho*V**2*b*(C) * alpha_max*np.exp(1j*omg*t)).real
    L_NC = (2*np.pi*rho*V**2*b*(1j*2*omg*b/V) * alpha_max*np.exp(1j*omg*t)).real
    L = L_C + L_NC
    Cl = (L_C/(rho*V**2*b)) + 2*np.pi*alpha_0
    print('V', V)
    print('Cl', L_C/(rho*V**2*b))
    return t, Cl, alpha
# Cl = 2*np.pi*alpha_max*(F*np.cos(omg*t)+G*2*omg*b/V*np.cos(omg*t))

# plottings
fig, ax = plt.subplots(1, 3, figsize=(45, 8))
# ax[0].plot(t, np.rad2deg(alpha), label=r'$\alpha$')
# ax[0].plot(t, np.rad2deg(alpha_dot), label=r'$\dot{\alpha}$')
# ax[0].plot(t, np.rad2deg(alpha_dot_dot), label=r'$\ddot{\alpha}$')
# ax[0].set_xlabel('t')
# ax[0].set_ylabel(r'$\alpha$ [deg]')
# ax[0].legend()
# ax[0].grid()
for i in [1e-8, 0.1, 0.2, 1e8]:
    print(i)
    t, Cl, alpha = lift(k=i)
    ax[0].plot(t, np.rad2deg(alpha.real), label=r'$\alpha$')

    ax[1].plot(t, Cl, label=r'$C_l=$'+str(i))
    ax[2].plot(np.rad2deg(alpha.real), Cl, label=r'$C_l$='+str(i))
# ax[1].plot(lift(1e-10), label=r'$C_l0$')
# ax[1].plot(lift(0.1), label=r'$C_l0.1$')
# ax[1].plot(lift(0.2), label=r'$C_l0.2$')
# ax[1].plot(lift(1e10), label=r'$C_linf$')
ax[1].set_xlabel('t')
ax[1].set_ylabel(r'$C_l$')
ax[1].legend()
ax[2].legend()
ax[1].grid()
ax[2].grid()

# plot cl vs alpha
# ax[2].plot(np.rad2deg(alpha), Cl, label=r'$C_l$')
plt.savefig('pure_aoa.png',transparent=True, bbox_inches='tight', pad_inches=0.1, dpi=400)
plt.show()

