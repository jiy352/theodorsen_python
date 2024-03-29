import numpy as np
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import scipy.special as scsp

import matplotlib.pyplot as plt


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

    h_max = 0.1
    alpha_0 = -np.deg2rad(0)
    n_period = 2
    t = np.linspace(0, n_period*2*np.pi/omg, 100)

    # profile of the angle of attack
    alpha = h_max * np.cos(omg*t) + alpha_0
    alpha_dot = -h_max * omg * np.sin(omg*t)
    alpha_dot_dot = -h_max * omg**2 * np.cos(omg*t)

    h_complex = h_max*np.exp(1j*omg*t)+ alpha_0
    alpha_dot = -h_max * omg * np.sin(omg*t)
    alpha_dot_dot = -h_max * omg**2 * np.cos(omg*t)

    h_complex = h_max*np.exp(1j*omg*t)
    alpha_dot = h_max*np.exp(1j*omg*t) * 1j*omg
    alpha_dot_dot = (1j*omg)**2 * h_max*np.exp(1j*omg*t)

    difference_norm = np.linalg.norm(alpha-h_complex.real)
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

    # L_C = (2*np.pi*rho*V**2*b*(C) * h_max*np.exp(1j*omg*t)).real
    L_C = (2*np.pi*rho*V**2*b*(C) * h_max*np.exp(1j*omg*t)).real
    L_NC = (2*np.pi*rho*V**2*b*(1j*2*omg*b/V) * h_max*np.exp(1j*omg*t)).real
    if circulatory==False:
        L = L_C + L_NC
    else:
        L = L_C
    Cl_0 = 2*np.pi*alpha_0
    Cl_C = 2 * np.pi * k * (1j*F-G) / b * h_complex
    Cl_NC = - np.pi * k**2 / b * h_complex 
    Cl = Cl_C + Cl_NC + Cl_0
    L = Cl * rho * V**2 * b
    print('V', V)
    # print('Cl', L_C/(rho*V**2*b))
    np.savetxt('theodorsen/alpha'+str(k)+'.txt', np.rad2deg(alpha))
    np.savetxt('theodorsen/Cl'+str(k)+'.txt', Cl.real)
    np.savetxt('theodorsen/L'+str(k)+'.txt', L.real)
    return t, Cl, h_complex.real,alpha_dot,alpha_dot_dot, C


# circulatory = False
circulatory = True
# plottings
fig, ax = plt.subplots(1, 4, figsize=(60/3, 8/3))
k_list = [1e-8, 0.1, 0.2, 1e8]
k_list = [0.2, 0.6, 1, 3]
k_list = [1e-8, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]
k_list = [0.5]
Cl_list_max = []
for i in k_list:
# for i in [0.1, 0.4, 0.7, 1]:
    print(i)
    t, Cl, alpha,alpha_dot,alpha_dot_dot, C = lift(k=i,circulatory=circulatory)
    Cl_list_max.append(Cl.real.max())


    ax[1].plot(t, Cl.real-np.pi*2*np.deg2rad(5), label=r'$C_l=$'+str(i))
    ax[2].plot(np.rad2deg(alpha.real), Cl.real, label=r'$C_l $'+' = '+str(i))
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

plt.figure()
c = 2 # chord
b = c * 0.5 # half chord
rho = 1.0

# define kinematics
omg = 1

# compute freestream velocity from reduced frequency
V = omg*b/np.array(k_list)
cl = np.array(Cl_list_max)*rho*V**2*b
cl[0] = 0
plt.plot(np.array(k_list), np.array(Cl_list_max),'o-')
plt.show()

