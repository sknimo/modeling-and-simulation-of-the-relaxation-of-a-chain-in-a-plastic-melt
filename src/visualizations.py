import matplotlib.pyplot as plt
import pandas as pd


# plotting the relaxation
# plt.style.use('ggplot')

data_relaxation = pd.read_excel('../data/outputs/stress_relaxation.xlsx')


fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(111)
ax.loglog(data_relaxation.iloc[:, 0],
          data_relaxation.iloc[:, 1], color='chocolate', lw=3, label='HDPE')
ax.set_xlabel(r'$\mathbf{time[s]}$', fontsize=12)
ax.set_ylabel(r'$\mathbf{G[Pa]}$', fontsize=12)
plt.title('Stress relaxation modulus of a plastic melt', fontsize=14)
plt.legend()
# plt.savefig('../reports/figures/relaxation.png', bbox_inches='tight', dpi=800)
plt.show()

# ploting the primitive paths
data_ppf = pd.read_excel('data/ppf_297831.xlsx')
plt.style.use('seaborn-v0_8-bright')
fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111)
ax.plot(data_ppf.time, data_ppf.early_ppf, 'b', lw=3, label=r'$\tau_{early}$')
ax.plot(data_ppf.time, data_ppf.late_ppf, 'g', lw=3, label=r'$\tau_{late}$')
ax.plot(data_ppf.time, data_ppf.ppf, 'k:', label='ppf')
ax.set_xlabel('$s$', fontsize=14)
ax.set_ylabel(r'$\tau_{s}$', fontsize=14)
plt.semilogy()
plt.legend()
# plt.savefig('../reports/figures/ppf.png', bbox_inches='tight', dpi=800)
plt.show()
