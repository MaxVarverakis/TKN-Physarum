import numpy as np
import matplotlib.pyplot as plt

data = np.concatenate(tuple(np.loadtxt(f'./kdata/{i}.txt') for i in range(3)))

print(f'Mean : {np.mean(data):.2e} \t Sigma : {np.std(data):.2e}')

plt.figure()
plt.hist(data, bins=60)
plt.xlabel('$\Delta \mathcal{F}$', fontsize=16)
# plt.xlabel('Approximate $\lambda$')
plt.ylabel('Count', fontsize=16)
plt.xscale('log')
plt.show()
