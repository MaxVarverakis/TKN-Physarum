import numpy as np
import matplotlib.pyplot as plt

resolution = 15
data = np.concatenate(tuple(np.loadtxt(f'./parallel_data_{resolution}/{i}.txt') for i in range(1,8)))

plt.figure(1)
h,b,_ = plt.hist(data, bins=40)
b = b[ b >= 0 ]
logbins = np.logspace(np.log10(b[0]),np.log10(b[-1]),len(b)-1)
plt.close()

# print(f'Mean : {np.mean(data):.2e} \t Sigma : {np.std(data):.2e}')

plt.figure(2)
plt.hist(data, bins=logbins)
plt.xlabel('$\Delta \mathcal{F}$', fontsize=16)
# plt.xlabel('Approximate $\lambda$')
plt.ylabel('Count', fontsize=16)
plt.xscale('log')
plt.show()
