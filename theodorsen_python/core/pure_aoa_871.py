import numpy as np
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import scipy.special as scsp

import matplotlib.pyplot as plt

from theodorson import theodorson_function, Lift, generate_kinemtics

def theo_fun(k):
    '''Returns Theodorsen function at a reduced frequency k'''

    H1=scsp.hankel2(1,k)
    H0=scsp.hankel2(0,k)

    C=H1/(H1+1.j*H0)

    return C

def lift(k, circulatory=True):
    # test generate_kinematics
    a = -1/2 # quarter chord
    c = 2 # chord
    b = c * 0.5 # half chord
    rho = 1.0

    # define kinematics
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

    alpha_complex = alpha_max*np.exp(1j*omg*t)+ alpha_0
    alpha_dot = -alpha_max * omg * np.sin(omg*t)
    alpha_dot_dot = -alpha_max * omg**2 * np.cos(omg*t)

    difference_norm = np.linalg.norm(alpha-alpha_complex.real)
    print('difference_norm', difference_norm)

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
    if circulatory==False:
        L = L_C + L_NC
    else:
        L = L_C
    Cl = (L/(rho*V**2*b)) + 2*np.pi*alpha_0
    print('V', V)
    # print('Cl', L_C/(rho*V**2*b))
    return t, Cl, alpha,alpha_dot,alpha_dot_dot, C


# circulatory = False
circulatory = True
# plottings
fig, ax = plt.subplots(1, 4, figsize=(60, 8))
k_list = [1e-8, 0.1, 0.2, 1e8]
for i in k_list:
# for i in [0.1, 0.4, 0.7, 1]:
    print(i)
    t, Cl, alpha,alpha_dot,alpha_dot_dot, C = lift(k=i,circulatory=circulatory)


    ax[1].plot(t, Cl, label=r'$C_l=$'+str(i))
    ax[2].plot(np.rad2deg(alpha.real), Cl, label=r'$C_l $'+'circulatory = '+str(i))
ax[0].plot(t, np.rad2deg(alpha.real), label=r'$\alpha$')
ax[0].plot(t, np.rad2deg(alpha_dot.real), label=r'$\dot{\alpha}$')
ax[0].plot(t, np.rad2deg(alpha_dot_dot.real), label=r'$\ddot{\alpha}$')

k = np.linspace(0.01, 4, 10000)
C = theo_fun(k)
# subplot the F and G and phase angles
ax[3].plot(k, np.rad2deg(np.angle(C)),label='circulatory angle')
ax[3].plot(k, np.rad2deg(np.angle(1j*k/2)),label='non-circulatory angle')
ax[3].plot(k, np.rad2deg(np.angle(C+1j*k/2)), label='total angle')

ax[1].set_xlabel('t')
ax[1].set_ylabel(r'$C_l$')
ax[1].legend()
ax[0].legend()


ax[2].set_xlabel(r'$\alpha$')
ax[2].set_ylabel(r'$C_l$'+'  circulatory')
ax[2].legend()
ax[3].legend()

ax[1].grid()
ax[2].grid()
ax[3].grid()
# plot cl vs alpha
plt.savefig('pure_aoa.png',transparent=True, bbox_inches='tight', pad_inches=0.1, dpi=400)
plt.show()

