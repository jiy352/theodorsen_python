import numpy as np
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)

import matplotlib.pyplot as plt
omg = 1
t = np.linspace(0, 2*np.pi, 100)
plt.plot(t, np.exp(1j*omg*t).real)