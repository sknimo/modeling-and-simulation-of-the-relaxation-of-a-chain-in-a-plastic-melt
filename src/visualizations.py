import matplotlib.pyplot as plt
import pandas as pd


plt.style.use('ggplot')

data_relaxation = pd.read_excel('../data/outputs/stress_relaxation.xlsx')


fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(111)
ax.loglog(data_relaxation.iloc[:, 0], data_relaxation.iloc[:, 1], label='HDPE')
ax.set_xlabel(r'$\mathbf{time[s]}$', fontsize=12)
ax.set_ylabel(r'$\mathbf{G[Pa]}$', fontsize=12)
plt.title('Stress relaxation modulus of a plastic melt', fontsize=14)
plt.legend()
plt.savefig('../reports/figures/relaxation.png', bbox_inches='tight', dpi=800)
plt.show()
