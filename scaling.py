import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
N = []
T = []
with open(file='scaling_results_parallel.txt',mode="r") as file:
    for line in file:
        if(line[0].isnumeric()):
            sp = line.split()
            N.append(int(sp[0]))
            T.append(float(sp[1]))
N = np.array(N)
T = np.array(T)

model = lambda n, a: a * n * np.log(n)

# Fit the model
params, _ = curve_fit(model, N, T)
a = params[0]

print(f"Fitted model: T(N) ≈ {a:.3e} * N log N ")

# Predict and plot
N_fit = np.linspace(min(N), max(N), 200)
T_fit = model(N_fit, a)
fig, ax = plt.subplots(1,1)
ax.scatter(N, T, label="Measured", color="blue")
ax.plot(N_fit, T_fit, label=f"Fit", color="red")
ax.set_xlabel("N")
ax.set_ylabel("Time (s)")
ax.set_xscale("log")
ax.set_yscale("log")

ax.plot(N_fit, a*N_fit*np.log(N_fit), label='NLogN')
#ax.plot(N_fit, c*np.ones(len(N_fit)), label='Constant overhead')
#ax.plot(N_fit, d*N_fit*N_fit, label='N^2')
ax.legend()
plt.show()