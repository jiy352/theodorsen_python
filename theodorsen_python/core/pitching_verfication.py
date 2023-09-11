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

def lift(k, circulatory=True, a = -1/2):
    # test generate_kinematics
    # quarter chord corresponds to a = -1/2
    c = 2 # chord
    b = c * 0.5 # half chord
    rho = 1.0

    # define kinematics
    omg = 1

    # compute freestream velocity from reduced frequency
    V = omg*b/k

    alpha_max = np.deg2rad(1)
    alpha_0 = np.deg2rad(0)
    n_period = 2
    t = np.linspace(0, n_period*2*np.pi/omg, 100)

    # profile of the angle of attack

    alpha = alpha_max * np.cos(omg*t) + alpha_0
    # alpha_dot = -alpha_max * omg * np.sin(omg*t)
    # alpha_dot_dot = -alpha_max * omg**2 * np.cos(omg*t)

    # alpha_complex = alpha_max*np.exp(1j*omg*t)+ alpha_0
    # alpha_dot = -alpha_max * omg * np.sin(omg*t)
    # alpha_dot_dot = -alpha_max * omg**2 * np.cos(omg*t)

    # alpha_complex = alpha_max*np.exp(1j*omg*t)
    # alpha_dot = alpha_max*np.exp(1j*omg*t) * 1j*omg
    # alpha_dot_dot = (1j*omg)**2 * alpha_max*np.exp(1j*omg*t)

    alpha_complex = alpha_max*np.exp(1j*omg*t) + alpha_0
    alpha_dot = alpha_max*np.exp(1j*omg*t) * 1j*omg
    alpha_dot_dot = (1j*omg)**2 * alpha_max*np.exp(1j*omg*t)

    difference_norm = np.linalg.norm(alpha-alpha_complex.real)
    print('difference_norm', difference_norm)

    # compute theodorsen function
    H1=scsp.hankel2(1,k)
    H0=scsp.hankel2(0,k)
    C=H1/(H1+1.j*H0)
    F = C.real
    G = C.imag


    # L_C = (2*np.pi*rho*V**2*b*(C) * alpha_max*np.exp(1j*omg*t)).real

    Cl_0 = 2*np.pi*alpha_0
    # Cl_C = 2 * np.pi * ((F*(1+1j*k)) + G*(1j - k)) * alpha_complex
    # Cl_NC = np.pi * k * (1j - k/2) * alpha_complex 
    Cl_C = 2 * np.pi * (F+1j*G + 1j*k*(0.5-a)*(F+1j*G)) * alpha_complex
    Cl_NC = np.pi * k * (1j + a * k) * alpha_complex 
    Cl = Cl_C + Cl_NC + Cl_0
    L = Cl * rho * V**2 * b
    print('V', V)
    # print('Cl', L_C/(rho*V**2*b))
    return t, Cl,L, alpha,alpha_dot,alpha_dot_dot, C


# circulatory = False
circulatory = True
# plottings
fig, ax = plt.subplots(1, 4, figsize=(60/3, 8/3))
k_list = [1e-8, 0.1, 0.2, 1e8]
k_list = [0.2, 0.6, 1, 3]
for i in k_list:
# for i in [0.1, 0.4, 0.7, 1]:
    print(i)
    t, Cl,L, alpha,alpha_dot,alpha_dot_dot, C = lift(k=i,circulatory=circulatory)
    np.savetxt('theodorsen/L_1deg'+str(i)+'.txt', L.real)
    np.savetxt('theodorsen/Cl_1deg'+str(i)+'.txt', Cl.real)
    np.savetxt('theodorsen/alpha_1deg'+str(i)+'.txt', alpha.real)

    ax[1].plot(t, Cl, label=r'$C_l=$'+str(i))
    ax[2].plot(np.rad2deg(alpha.real), Cl, label=r'$C_l $'+' = '+str(i))
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
plt.savefig('pure_aoa_qq.png',transparent=True, bbox_inches='tight', pad_inches=0.1, dpi=400)

fig, ax = plt.subplots(1, 4, figsize=(60/3, 8/3))

for i in k_list:
# for i in [0.1, 0.4, 0.7, 1]:
    print(i)
    t, Cl,L, alpha,alpha_dot,alpha_dot_dot, C = lift(k=i,circulatory=circulatory,a=0)
    np.savetxt('theodorsen/L_1deg_middle'+str(i)+'.txt', L.real)
    np.savetxt('theodorsen/Cl_1deg_middle'+str(i)+'.txt', Cl.real)
    np.savetxt('theodorsen/alpha_1deg_middle'+str(i)+'.txt', alpha.real)

    ax[1].plot(t, Cl, label=r'$C_l=$'+str(i))
    ax[2].plot(np.rad2deg(alpha.real), Cl, label=r'$C_l $'+' = '+str(i))
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
plt.savefig('pure_ao_m.png',transparent=True, bbox_inches='tight', pad_inches=0.1, dpi=400)

plt.show()

