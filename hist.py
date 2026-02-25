import numpy as np
import matplotlib.pyplot as plt

resolution = 20
data = np.concatenate(tuple(np.loadtxt(f'./parallel_data_1e-4_{resolution}/{i}_{j}.txt') for i in range(1,8) for j in range(10)))
# data = np.concatenate(tuple(np.loadtxt(f'./parallel_data_1e-4_{resolution}_edge/{i}_{j}.txt') for i in range(1,8) for j in range(9)))
# data = np.loadtxt('./data/data5.txt')

plt.plot(np.sort(data)[::-1], ls='', marker='.', c='b', alpha=0.5)
plt.xlabel('Rank', fontsize=16)
plt.ylabel(r'$\frac{v^\top H v}{v^\top v} \quad \left[\left(\mathcal{F}^*\right)^{-1}\right]$', fontsize=16, labelpad=1)
# plt.ylabel('$\Delta \mathcal{F}$', fontsize=16)
plt.yscale('log')
plt.show()

# print(np.sort(data)[-10:])

h,b,_ = plt.hist(data, bins=60)
b = b[ b >= 0 ]
logbins = np.logspace(np.log10(b[0]),np.log10(b[-1]),len(b)-1)
plt.close()

# print(f'Mean : {np.mean(data):.2e} \t Sigma : {np.std(data):.2e}')
print(f'Max : {max(data):.2e} \t Min : {min(data):.2e}')
# print(len(data))
# n=100
# print(data[np.argsort(data)[-n:]])

# plt.hist(data[abs(data)<1e2], bins=100)
# plt.hist(data, bins=60)
plt.hist(data, bins=logbins)
plt.xlabel(r'$\frac{v^\top H v}{v^\top v} \quad \left[\left(\mathcal{F}^*\right)^{-1}\right]$', fontsize=16, labelpad=1)
# plt.xlabel('Approximate $\lambda/\mathcal{F}^*$', fontsize=16)
# plt.xlabel('$\Delta \mathcal{F}$', fontsize=16)
plt.ylabel('Count', fontsize=16)
plt.xscale('log')
plt.yscale('log')
# plt.title('Edge removal fitness change')

plt.show()
