import numpy as np
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)

import matplotlib.pyplot as plt

def theodorson_function(k):
    '''
    The Theodorson function C(K) is defined as:
        F(K) = Re(C(K))
        G(K) = Im(C(K))
    '''
    F = (0.5005*k**3 + 0.51261*k**2 + 0.21040*k + 0.021573)/ \
        (k**3 + 1.03538*k**2 + 0.25124*k + 0.02151)

    G = - (0.00015*k**3 + 0.12240*k**2 + 0.32721*k + 0.001990)/ \
        (k**3 + 2.48148*k**2 + 0.93453*k + 0.08932)
    return F, G

def Lift(rho, V, a, b, alpha, alpha_dot, alpha_ddot, k, h_dot=0, h_ddot=0):
    '''
    The lift coefficient is defined as:
        C_L = 2*F(K)*sin(alpha) + 2*G(K)*cos(alpha)
    '''
    F, G = theodorson_function(k)
    print(F, G)
    L_NC = np.pi*rho*b**2*(h_ddot + V*alpha_dot - b*a*alpha_ddot)
    L_C = 2*np.pi*rho*V*b*(h_dot + alpha*V + b*alpha_dot*(1/2-a))* (F+G*1j)
    return L_C

def generate_kinemtics(a, b, V, alpha_max, omg, n_period):
    k = omg * b / V
    alpha = alpha_max * np.exp(1j*omg*t).real
    alpha_dot = alpha_max * np.exp(1j*omg*t) * 1j * omg
    alpha_ddot = alpha_max * np.exp(1j*omg*t) * (1j * omg)**2
    return k, alpha, alpha_dot, alpha_ddot

if __name__ == '__main__':
    # test generate_kinematics
    a = -1/2
    b = 1
    V = 10
    omg = 1
    alpha_max = np.deg2rad(5)
    n_period = 2
    t = np.linspace(0, n_period*2*np.pi/omg, 100)
    k, alpha, alpha_dot, alpha_ddot = generate_kinemtics(a, b, V, alpha_max, omg, n_period)

    # plot alpha, alpha_dot, alpha_ddot
    plt.figure()
    plt.plot(t, np.rad2deg(alpha), label='alpha')
    plt.plot(t, np.rad2deg(alpha_dot.real), label='alpha_dot')
    plt.plot(t, np.rad2deg(alpha_ddot.real), label='alpha_ddot')
    plt.legend(['alpha', 'alpha_dot', 'alpha_ddot'])


    

    # plt.plot(t, np.rad2deg(alpha.imag), label='imag')
    plt.xlabel('t')
    plt.ylabel('alpha')
    # plt.show()


    # test Lift
    rho = 1.0
    L_C = Lift(rho, V, a, b, alpha, alpha_dot, alpha_ddot, k)
    plt.figure()
    plt.plot(t, L_C.real, label='L_C_real')
    plt.plot(t, L_C.imag, label='L_C_imag')
    plt.xlabel('t')
    plt.ylabel('L_C')
    plt.show()
    exit()



    # Plot the Theodorson function

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
    plt.xlim(min(F) - 0.1, max(F) + 0.1)
    plt.ylim(min(G) - 0.1, max(G) + 0.1)

    plt.xlabel('F(k)')
    plt.ylabel('G(k)')
    plt.title('Theodorson Function')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()

    plt.savefig('theodorson_plot.png', dpi=300)
    plt.show()
