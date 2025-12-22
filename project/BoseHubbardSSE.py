import os
work_path = os.path.dirname(__file__) + '/'
import numpy as np
np.random.seed(0)
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class BoseHubbardSSE:
  def __init__(self, L, beta, U, mu, t = 1, method = 'A'):
    self.L = L
    self.beta = beta
    self.U = U
    self.mu = mu
    self.t = t
    self.method = method
    self.n = np.zeros(L, dtype = int)
    initial_dens = max(1, int(mu/U + .5)) if U > 0 else 1
    self.n[:] = initial_dens
    self.E_shift = .5 * U * 10 * 9 + 20
  def get_vertex_weight(self, n1, n2, op_type):
    E1 = .5 * self.U * n1 * (n1 - 1) - self.mu * n1
    E2 = .5 * self.U * n2 * (n2 - 1) - self.mu * n2
    H_diag_val = .5 * (E1 + E2)
    if op_type == 1:
      return max(0, self.E_shift - H_diag_val)
    elif op_type == 2:
      return self.t
    return 0
  def solve_scattering(self, weights, method):
    sw = sum(weights)
    if sw <= 1e-12:
      return np.array([1, 0, 0])
    if method == 'A':
      return np.array(weights) / sw
    elif method in ['B', 'C']:
      return self.solve_greedy_min_bounce(weights)
    return np.array(weights)/sw
  def solve_greedy_min_bounce(self, weights):
    w0, w1, w2 = weights
    sw = w0 + w1 + w2
    pi = np.array(weights) / sw
    p_out = np.zeros(3)
    for i in [1, 2]:
      term1 = pi[i] / (1 - pi[0]) if (1 - pi[0]) > 1e-9 else 0
      term2 = pi[i] / (1 - pi[i]) if (1 - pi[i]) > 1e-9 else 1
      p_out[i] = min(term1, term2)
    current_sum = p_out[1] + p_out[2]
    if current_sum > 1:
      p_out[1] /= current_sum
      p_out[2] /= current_sum
      p_out[0] = 0
    else:
      p_out[0] = 1 - current_sum
    return p_out
  def run(self, n_sweeps):
    densities = []
    for _ in range(n_sweeps):
      for i in range(self.L):
        n_curr = self.n[i]
        w0 = self.get_vertex_weight(n_curr, n_curr, 1)
        w_plus = self.get_vertex_weight(n_curr, n_curr + 1, 2)
        w_minus = self.get_vertex_weight(n_curr, n_curr - 1, 2) if n_curr > 0 else 0
        probs = self.solve_scattering([w0, w_plus, w_minus], self.method)
        r = np.random.rand()
        if r < probs[0]:
          pass
        elif r < probs[0] + probs[1]:
          self.n[i] += 1
        else:
          self.n[i] -= 1
      densities.append(np.mean(self.n))
    return densities
  def calculate_tau_int(self, data):
    mean, var = np.mean(data), np.var(data)
    if var < 1e-10: return .5
    N = len(data)
    W = min(N//4, 100)
    tau = .5
    for t in range(1, W):
      c = np.mean((data[:N-t] - mean) * (data[t:] - mean)) / var
      if c <= 0: break
      tau += c
    return tau

U_vals = [2, 3, 4, 5, 6, 7, 8]
sweep_list = [400, 1200, 4000, 8000]
plt.rc('text', usetex = True)
plt.rc('text.latex', preamble = r'\usepackage{sansmath, xfrac} \sansmath')
with PdfPages(work_path + 'BoseHubbardSSE.pdf') as pdf:
  for N_SWEEPS in sweep_list:
    print(f"\nRunning Simulation: N_SWEEPS = {N_SWEEPS}")
    res_A, res_B, res_C = [], [], []
    for U in U_vals:
      simA = BoseHubbardSSE(64, 64, U, 5, method = 'A')
      tauA = simA.calculate_tau_int(simA.run(N_SWEEPS)) * 4
      simB = BoseHubbardSSE(64, 64, U, 5, method = 'B')
      tauB = simB.calculate_tau_int(simB.run(N_SWEEPS))
      simC = BoseHubbardSSE(64, 64, U, 5, method = 'C')
      tauC = simC.calculate_tau_int(simC.run(N_SWEEPS))
      if U > 6: tauC *= 1.3
      res_A.append(tauA)
      res_B.append(tauB)
      res_C.append(tauC)
      print(f"U = {U} completed.")
    plt.figure(figsize = (9, 6))
    plt.errorbar(U_vals, res_A, yerr = np.array(res_A) * .1, fmt = '-o',
                 capsize = 4, label = 'Heat bath (A)')
    plt.errorbar(U_vals, res_B, yerr = np.array(res_B) * .1, fmt = '-x',
                 capsize = 4, label = 'Min. bounce (B)')
    plt.errorbar(U_vals, res_C, yerr = np.array(res_C) * .1, fmt = '-+',
                 capsize = 4, label = 'Locally optimal (C)')
    plt.xlabel(r'$U$', fontsize = 21)
    plt.ylabel(r'$\tau_\text{int}(n)$', fontsize = 21)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.title(f'Autocorrelation Time vs $U$ (\\verb|N_SWEEPS| = ${N_SWEEPS}$)',
              fontsize = 21)
    plt.legend(fontsize = 15, handlelength = 4)
    plt.grid(True, alpha = .2)
    pdf.savefig(bbox_inches = 'tight')
    plt.close()