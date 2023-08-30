
import numpy as np
import scipy.special as scsp
from IPython import embed
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import matplotlib.pyplot as plt
# imaginary variable
j=1.0j


def theo_fun(k):
	'''Returns Theodorsen function at a reduced frequency k'''

	H1=scsp.hankel2(1,k)
	H0=scsp.hankel2(0,k)

	C=H1/(H1+j*H0)

	return C

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

k = np.linspace(0.01, 100, 10000)
C = theo_fun(k)
F, G = theodorson_function(k)
plt.plot(C.real, C.imag)
plt.plot(F, G)
plt.show()