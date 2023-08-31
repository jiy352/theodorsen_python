import numpy as np
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import scipy.special as scsp

import matplotlib.pyplot as plt

def theodorson_function(k):
	'''Returns Theodorsen function at a reduced frequency k'''

	H1=scsp.hankel2(1,k)
	H0=scsp.hankel2(0,k)

	C=H1/(H1+1.j*H0)

	return C.real, C.imag

# def theodorson_function(k):
#     '''
#     The Theodorson function C(K) is defined as:
#         F(K) = Re(C(K))
#         G(K) = Im(C(K))
#     '''
#     F = (0.5005*k**3 + 0.51261*k**2 + 0.21040*k + 0.021573)/ \
#         (k**3 + 1.03538*k**2 + 0.25124*k + 0.02151)

#     G = - (0.00015*k**3 + 0.12240*k**2 + 0.32721*k + 0.001990)/ \
#         (k**3 + 2.48148*k**2 + 0.93453*k + 0.08932)
#     return F, G

def Lift(rho, V, a, b, alpha_max, omg, n_period, t, h_dot=0, h_ddot=0):
    '''
    The lift coefficient is defined as:
        C_L = 2*F(K)*sin(alpha) + 2*G(K)*cos(alpha)
    '''
    k, alpha, alpha_dot, alpha_ddot = generate_kinemtics(a, b, V, alpha_max, omg, n_period,t)
    F, G = theodorson_function(k)
    if abs(G)<1e-4:
        G = 0
    print(F, G)
    L_NC = np.pi*rho*b**2*(h_ddot + V*alpha_dot - b*a*alpha_ddot) 
    L_C = (2*np.pi*rho*V*b*(h_dot + alpha*V + b*alpha_dot*(1/2-a))* (F+G*1j))

    print('coefficients on alpha, alpha_dot, alpha_ddot:', 
            V**2*F*2*np.pi*b,
            np.pi*b**2*V+2*np.pi*b**2*V*F*b*(1/2-a),
            -np.pi*b**2*a*b)
    # CL = (L_NC + L_C)/(1/2*rho*V**2*b*2)
    CL = np.pi*b*(alpha_dot/V + h_ddot/(V**2) - b*a*alpha_ddot/(V**2)) + \
            2*np.pi*(F+G*1j)*(h_dot/V + alpha + b*alpha_dot*(1/2-a)/V)

    # CL = (np.pi*2*(F+G*1j) + 1.j*np.pi*k)*alpha 
    # CL = (CL.real**2 + CL.imag**2)**0.5
    return L_NC, L_C, CL

def generate_kinemtics(a, b, V, alpha_max, omg, n_period, t):
    k = omg * b / V
    print('k:--------------', k)
    print('omg:--------------', omg)
    # alpha = alpha_max * np.exp(1j*omg*t).real
    # alpha_dot = (alpha_max * np.exp(1j*omg*t) * 1j * omg).real
    # alpha_ddot = (alpha_max * np.exp(1j*omg*t) * (1j * omg)**2).real

    # alpha = alpha_max * np.cos(omg*t)+ alpha_max
    # alpha_dot = -omg * alpha_max * np.sin(omg*t)
    # alpha_ddot = -omg**2 * alpha_max * np.cos(omg*t)

    # alpha = alpha_max * np.sin(omg*t)
    # alpha_dot = omg * alpha_max * np.cos(omg*t)
    # alpha_ddot = -omg**2 * alpha_max * np.sin(omg*t)

    # alpha_max = np.deg2rad(5)

    alpha = alpha_max * np.cos(omg*t)
    alpha_dot = -(alpha_max * np.sin(omg*t)* omg)
    alpha_ddot = -(alpha_max * np.cos(omg*t)* omg**2)
    # print('alpha:--------------', alpha)

    return k, alpha, alpha_dot, alpha_ddot

if __name__ == '__main__':

    # test generate_kinematics
    a = -1/2
    b = 1
    V = np.array([1/3, 1/1, 1/0.6, 1/0.2])
    omg = 1
    alpha_max = np.deg2rad(1)
    n_period = 2
    t = np.linspace(0, n_period*2*np.pi/omg, 100)

    # test Lift
    rho = 1.0
    plt.figure()
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['axes.titlesize'] = 16
    plt.rcParams['lines.linewidth'] = 2
    plt.figure(figsize=(10, 6))  # Adjust the figure size
    lgd = ['k=3','k=1', 'k=0.6','k=0.2']

    for i in range(len(V)):
        L_NC, L_C, CL = Lift(rho, V[i], a, b, alpha_max, omg, n_period,t)
        # plt.figure()
        # plt.plot(t, L_C.real, label='L_C_real')
        # plt.plot(t, L_NC.real+L_C.real, label='L_C_real')
        # # plt.plot(t, L_C.imag, label='L_C_imag')
        # plt.xlabel('t')
        # plt.ylabel('L_C')
        T = 2*np.pi/omg
        alpha = alpha_max * np.exp(1j*omg*t).real 
        plt.figure()
        plt.plot(np.rad2deg(alpha), CL.real, label='Cl')
        # plt.plot(t/T, CL.real, label='Cl')
        # plt.plot(t, CL_k1.real/100, label='CL')
        plt.xlabel('t/T ')
        plt.xlabel('alpha ')
        plt.ylabel('Cl')
        plt.grid(True)
        plt.legend([lgd[i]])
        plt.xlim( - 1,  1)
        plt.ylim(-0.1, 0.1)
        plt.savefig('CL'+str(i)+'.png', dpi=400, transparent=True)
    # plt.ylim(-1.7, 1.7)

    # plt.show()
    # exit()

    k, alpha, alpha_dot, alpha_ddot = generate_kinemtics(a, b, V, alpha_max, omg, n_period)
    
    # plot alpha, alpha_dot, alpha_ddot
    plt.figure()
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['axes.titlesize'] = 16
    plt.rcParams['lines.linewidth'] = 2
    plt.figure(figsize=(10, 6))  # Adjust the figure size

    plt.plot(t, np.rad2deg(alpha), label='alpha')
    plt.plot(t, np.rad2deg(alpha_dot.real), label='alpha_dot')
    plt.plot(t, np.rad2deg(alpha_ddot.real), label='alpha_ddot')
    plt.legend(['alpha', 'alpha_dot', 'alpha_ddot'])
    plt.grid(True)
    plt.xlabel('Time')
    plt.ylabel('Degrees')
    plt.title('Angle and Derivatives vs Time')
    plt.legend()
    plt.tight_layout()
    plt.savefig('alpha.png', dpi=400)

    

    # plt.plot(t, np.rad2deg(alpha.imag), label='imag')
    # plt.xlabel('t')
    # plt.ylabel('alpha')
    # plt.show()


    # Plot the Theodorson function
    exit()
    k = np.linspace(0.0, 1000, 100000)
    F, G = theodorson_function(k)

    plt.figure(figsize=(8, 6))

    # Plot the F and G values
    plt.plot(F, G, label='Theodorson Function', color='blue')

    # Annotate specific points
    k_annotate = np.array([0.01, 0.2, 1])
    F_annotate, G_annotate = theodorson_function(k_annotate)

    annotations = ['K = 0.01', 'K = 0.2', 'K = 1']
    for i, annotation in enumerate(annotations):
        plt.scatter(F_annotate[i], G_annotate[i], color='red')
        plt.annotate(annotation, xy=(F_annotate[i], G_annotate[i]), 
                    xytext=(F_annotate[i] + 0.03, G_annotate[i] + 0.03),
                    arrowprops=dict(arrowstyle="->", facecolor='black'),
                    fontsize=10,
                    )

    # Adjust the plot limits to make annotations fit within the plot area
    plt.xlim( - 1,  1)
    plt.ylim(-0.1, 0.1)

    plt.xlabel('F(k)')
    plt.ylabel('G(k)')
    plt.title('Theodorson Function')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()

    plt.savefig('theodorson_plot.png', dpi=300)
    plt.show()
