import numpy as np
import matplotlib.pyplot as plt

D = np.loadtxt("conductances.txt", delimiter=',')
Dl = np.loadtxt("conductances_largeInit.txt", delimiter=',')

fs = 16

cmap = plt.get_cmap('ocean')
colors = cmap(np.linspace(0, 1, D.shape[1]))

fig,ax = plt.subplots(1, 2, figsize=(20,6), sharey=True, gridspec_kw={'wspace' : 0.02})

for i in range(D.shape[1]):
    # plt.plot(D[:,i], color=colors[i])
    ax[0].plot(D[:,i])
    ax[1].plot(Dl[:,i])

# plt.yscale('log')
for j in range(2):
    ax[j].set_xlim(right=700)
    ax[j].set_xlabel('Iteration', fontsize=fs)

ax[0].set_title(r'$\langle D_0\rangle = 0.1$', fontsize=fs)
ax[1].set_title(r'$\langle D_0\rangle = 3.1$', fontsize=fs)

ax[0].set_ylabel('Edge conductance', fontsize=fs)
# plt.savefig('./figs/D_SI.pdf', dpi=500, bbox_inches='tight')
plt.show()
