import numpy as np
import matplotlib.pyplot as plt



import numpy as np
def precompute_cluster_LRF(poles, weights, p):
    nodes = np.cos((np.arange(p) + 0.5) * np.pi / p)
    g = np.zeros(p)
    for i in range(p):
        for a in range(len(poles)):
            g[i] += weights[a]/(nodes[i]-poles[a])

    return g

def eval_cluster_LRF(x, g, u):
    s = 0
    for i in range(len(u)):
        s+= g[i]*u[i](x)
    return s


def precompute_cheby(p):
    nodes = np.cos((np.arange(p) + 0.5) * np.pi / p)
    denom = np.zeros(p)
    for j in range(p):
        denom[j] = np.prod([nodes[j] - nodes[m] for m in range(p) if m != j])
    u = []
    for i in range(p):
        u.append(lambda x, i=i: np.prod([(x - nodes[m]) for m in range(p) if m != i]) / denom[i])
    return u


def FMM(x, poles, weights, p):
    nodes = np.cos((np.arange(p) + 0.5) * np.pi / p)

    denom = np.zeros(p)
    for j in range(p):
        denom[j] = np.prod([nodes[j] - nodes[m] for m in range(p) if m != j])

    s = 0
    for i in range(p):
        # Compute g(nodes[i]) = sum_a weights[a]/(nodes[i] - poles[a])
        gi = np.sum(weights / (nodes[i] - poles))
        # Compute Lagrange basis u_i(x)
        ui = np.prod([(x - nodes[m]) for m in range(p) if m != i]) / denom[i]
        s += gi * ui
    return s



xs = np.linspace(-1.5,8,1000)
poles = np.array([3.1,4,5,6,7])
weights = np.array([1,2,3,4,5])
tol = 1e-5
p = np.ceil(-np.log(tol)/np.log(5)).astype(int)
g = precompute_cluster_LRF(poles, weights, p)
u = precompute_cheby(p)
interp = [eval_cluster_LRF(x, g, u) for x in xs]
plt.plot(xs, interp, label='FMM Interpolation')
exact = [sum(weights / (x - poles)) for x in xs]
plt.plot(xs, exact, label='Exact Function', linestyle='--')
plt.legend()
plt.xlabel('x')
plt.ylabel('f(x)')
plt.xlim(-1.5, 8)
plt.ylim(-5, 5)
plt.title(f'FMM Interpolation (p = {p}) vs Exact Function.')
plt.grid()
plt.show()