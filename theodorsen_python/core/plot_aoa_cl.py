import numpy as np
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import scipy.special as scsp

import matplotlib.pyplot as plt

from theodorson import theodorson_function, Lift, generate_kinemtics

alpha_file_name = 'theodorsen/alpha_1deg0.2.txt'
alpha_file_name_vast = 'theodorsen/alpha_theodorsen00.2.txt'
# L_file_name = ['theodorsen/L_1deg0.2.txt',
#                 'theodorsen/L_1deg0.6.txt',
#                 'theodorsen/L_1deg1.txt',
#                 'theodorsen/L_1deg3.txt']
Cl_file_name = ['theodorsen/Cl_1deg0.2.txt',
                'theodorsen/Cl_1deg0.6.txt',
                'theodorsen/Cl_1deg1.txt',
                'theodorsen/Cl_1deg3.txt']
Cl_middle_file_name = ['theodorsen/C_L_theodorsen00.2.txt',
                       'theodorsen/C_L_theodorsen00.6.txt',
                       'theodorsen/C_L_theodorsen01.txt',
                       'theodorsen/C_L_theodorsen03.txt']
def plot_cl(alpha_file_name, Cl_file_name, fig, ax, L_file_name=None):


    
    for i in range(len(Cl_file_name)):
        alpha = np.loadtxt(alpha_file_name)
        if alpha.max()<0.2:
            alpha = np.rad2deg(alpha)
        Cl = np.loadtxt(Cl_file_name[i])
        print('i', Cl_file_name[i])

        if L_file_name is not None:
            L = np.loadtxt(L_file_name[i])
        else:
            L = None
        if np.loadtxt(alpha_file_name).max()>0.2:
            ax.plot(alpha, Cl,'--')
        else:
            ax.plot(alpha, Cl)
        # ax[1].plot(alpha, Cl)

if __name__ == '__main__':
    fig, ax = plt.subplots(1, 1, figsize=(60/9, 8/3))
    plot_cl(alpha_file_name, Cl_file_name, fig, ax,L_file_name=None)
    # plt.show()
    # fig, ax = plt.subplots(1, 1, figsize=(60/9, 8/3))

    plot_cl(alpha_file_name_vast, Cl_middle_file_name, fig, ax, L_file_name=None)
    plt.legend(['theodorsen 0.2', 'theodorsen 0.6', 'theodorsen 1', 'theodorsen 3','VAST 0.2', 'VAST 0.6', 'VAST 1', 'VAST 3'])
    plt.xlabel(r'$\alpha$')
    plt.ylabel(r'$C_l$')
    plt.show()


# # circulatory = False
# circulatory = True
# # plottings
# fig, ax = plt.subplots(1, 4, figsize=(60/3, 8/3))
# k_list = [1e-8, 0.1, 0.2, 1e8]
# k_list = [0.2, 0.6, 1, 3]
# for i in k_list:
# # for i in [0.1, 0.4, 0.7, 1]:
#     print(i)
#     t, Cl,L, alpha,alpha_dot,alpha_dot_dot, C = lift(k=i,circulatory=circulatory)
#     np.savetxt('theodorsen/L_1deg'+str(i)+'.txt', L.real)
#     np.savetxt('theodorsen/Cl_1deg'+str(i)+'.txt', Cl.real)
#     np.savetxt('theodorsen/Cl_1deg'+str(i)+'.txt', alpha.real)

#     ax[1].plot(t, Cl, label=r'$C_l=$'+str(i))
#     ax[2].plot(np.rad2deg(alpha.real), Cl, label=r'$C_l $'+' = '+str(i))
# ax[0].plot(t, np.rad2deg(alpha.real), label=r'$\alpha$')
# ax[0].plot(t, np.rad2deg(alpha_dot.real), label=r'$\dot{\alpha}$')
# ax[0].plot(t, np.rad2deg(alpha_dot_dot.real), label=r'$\ddot{\alpha}$')

# k = np.linspace(0.01, 4, 10000)
# C = theo_fun(k)
# # subplot the F and G and phase angles
# ax[3].plot(k, np.rad2deg(np.angle(C)),label='circulatory angle')
# ax[3].plot(k, np.rad2deg(np.angle(1j*k/2)),label='non-circulatory angle')
# ax[3].plot(k, np.rad2deg(np.angle(C+1j*k/2)), label='total angle')

# ax[1].set_xlabel('t')
# ax[1].set_ylabel(r'$C_l$')
# ax[1].legend()
# ax[0].legend()


# ax[2].set_xlabel(r'$\alpha$')
# ax[2].set_ylabel(r'$C_l$'+'  circulatory')
# ax[2].legend()
# ax[3].legend()

# ax[1].grid()
# ax[2].grid()
# ax[3].grid()
# # plot cl vs alpha
# plt.savefig('pure_aoa_qq.png',transparent=True, bbox_inches='tight', pad_inches=0.1, dpi=400)

# fig, ax = plt.subplots(1, 4, figsize=(60/3, 8/3))

# for i in k_list:
# # for i in [0.1, 0.4, 0.7, 1]:
#     print(i)
#     t, Cl,L, alpha,alpha_dot,alpha_dot_dot, C = lift(k=i,circulatory=circulatory,a=0)
#     np.savetxt('theodorsen/L_1deg_middle'+str(i)+'.txt', L.real)
#     np.savetxt('theodorsen/Cl_1deg_middle'+str(i)+'.txt', Cl.real)
#     np.savetxt('theodorsen/Cl_1deg_middle'+str(i)+'.txt', alpha.real)

#     ax[1].plot(t, Cl, label=r'$C_l=$'+str(i))
#     ax[2].plot(np.rad2deg(alpha.real), Cl, label=r'$C_l $'+' = '+str(i))
# ax[0].plot(t, np.rad2deg(alpha.real), label=r'$\alpha$')
# ax[0].plot(t, np.rad2deg(alpha_dot.real), label=r'$\dot{\alpha}$')
# ax[0].plot(t, np.rad2deg(alpha_dot_dot.real), label=r'$\ddot{\alpha}$')

# k = np.linspace(0.01, 4, 10000)
# C = theo_fun(k)
# # subplot the F and G and phase angles
# ax[3].plot(k, np.rad2deg(np.angle(C)),label='circulatory angle')
# ax[3].plot(k, np.rad2deg(np.angle(1j*k/2)),label='non-circulatory angle')
# ax[3].plot(k, np.rad2deg(np.angle(C+1j*k/2)), label='total angle')

# ax[1].set_xlabel('t')
# ax[1].set_ylabel(r'$C_l$')
# ax[1].legend()
# ax[0].legend()


# ax[2].set_xlabel(r'$\alpha$')
# ax[2].set_ylabel(r'$C_l$'+'  circulatory')
# ax[2].legend()
# ax[3].legend()

# ax[1].grid()
# ax[2].grid()
# ax[3].grid()
# # plot cl vs alpha
# plt.savefig('pure_ao_m.png',transparent=True, bbox_inches='tight', pad_inches=0.1, dpi=400)

# plt.show()

