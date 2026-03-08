import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import wasserstein_distance

fs = 16
resolution = 20
# data = np.loadtxt(f'./parallel_data_1e-4_many_graphs_20_clamp/1.txt')
data = np.concatenate(tuple(np.loadtxt(f'./parallel_data_1e-4_many_graphs_{resolution}_clamp_5/{8 * j + i}.txt') for i in range(8) for j in range(112))) # cost nonzero
# data = np.concatenate(tuple(np.loadtxt(f'./parallel_data_1e-4_many_graphs_{resolution}_clamp_4/{8 * j + i}.txt') for i in range(8) for j in range(254))) # cost nonzero
# data = np.concatenate(tuple(np.loadtxt(f'./parallel_data_1e-4_many_graphs_{resolution}_clamp_3/{8 * j + i}.txt') for i in range(8) for j in range(330))) # cost zero
# data = np.concatenate(tuple(np.loadtxt(f'./parallel_data_1e-4_many_graphs_{resolution}_clamp_2/{8 * j + i}.txt') for i in range(8) for j in range(11)))
# data = np.concatenate(tuple(np.loadtxt(f'./parallel_data_1e-4_many_graphs_{resolution}_clamp/{8 * j + i}.txt') for i in range(8) for j in range(360)))

data = data[data > 0]

# plt.semilogy(np.sort(data)[::-1], ls='', marker='.', c='b', alpha=0.5)
# plt.xlabel('Rank', fontsize=fs)
# plt.ylabel(r'$\frac{v^\top H v}{v^\top v} \quad \left[\left(\mathcal{F}^*\right)^{-1}\right]$', fontsize=fs, labelpad=1)
# # plt.ylabel('$\Delta \mathcal{F}$', fontsize=fs)
# # plt.yscale('log')
# plt.show()

# print(f'Mean : {np.mean(data):.2e} \t Sigma : {np.std(data):.2e}')
# print(f'Max : {max(data):.2e} \t Min : {min(data):.2e}')
# print(len(data))
# n=100
# print(data[np.argsort(data)[-n:]])

sorted_data = np.sort(data)
cumulative = np.ones(len(sorted_data))
cumulative_data = np.cumsum(cumulative) / len(sorted_data)

c1 = 'steelblue'
c2 = 'tab:red'
bn = np.logspace(-7,2,120)

counts, edges, _ = plt.hist(data, bins=bn, log=True, weights = np.ones_like(data) / len(data), color=c1, alpha=0.85)
# plt.hist(data, bins=np.logspace(-6,2,100), weights = np.ones_like(data) / len(data), cumulative=True, histtype='step', color='k')


ax = plt.gca()
ax.set_xlabel('Rayleigh quotient ' + r'$\left[\left(\mathcal{F}^*\right)^{-1}\right]$', fontsize=fs, labelpad=1)
ax.set_ylabel('Density', fontsize=fs) #, color=c1)
# ax.spines["right"].set_color(c1)
# ax.tick_params(axis='y', colors=c1)
ax.set_xscale('log')
ax.xaxis.set_minor_locator(ticker.LogLocator(subs='all', numticks=10))
# ax.set_yscale('log')

# cnt, _ = np.histogram(data, bins=bn)
# centers = 0.5 * (edges[1:] + edges[:-1])
# widths = np.diff(edges)
# sigma = [ np.sqrt(c) / (width * len(data)) if (c > 0) else 0 for c,width in zip(cnt, widths) ]
# # sigma = np.sqrt(counts) / len(data)
# ax.errorbar(centers, counts, yerr=sigma, fmt='.')
# print(1/np.sqrt(max(cnt)))

ax2 = ax.twinx()
ax2.plot(sorted_data, cumulative_data, c=c2)
ax2.set_ylim(-0.01,1.01)
# ax2.set_xlim(min(sorted_data), max(sorted_data))
ax2.spines["right"].set_color(c2)
ax2.tick_params(axis='y', colors=c2)
# ax2.set_yscale('log')
ax2.set_ylabel('Cumulative probability', fontsize=fs, color=c2)

# plt.savefig('./figs/HDensity_5.svg', dpi=500, bbox_inches='tight')
# plt.title('Edge removal fitness change')

plt.show()


# plt.plot(sorted_data, cumulative_data)
# plt.hist(data, bins=np.logspace(-6,3,100), weights = np.ones_like(data) / len(data), cumulative=True, histtype='step')
# plt.xlabel('Rayleigh quotient ' + r'$\left[\left(\mathcal{F}^*\right)^{-1}\right]$', fontsize=fs, labelpad=1)
# plt.ylabel('Cumulative probability', fontsize=fs)
# plt.xscale('log')
# plt.show()


# plt.hist(np.log10(data[data>0]), log=True, density=True, bins=100, color=c1, alpha=0.85)
# # plt.xscale('log')
# plt.show()

recompute_convergence = False

n_values = 24
n_trials = 50
n_samples = np.unique(np.logspace(2, 5, n_values).astype(int))
if recompute_convergence:
    mean_D = np.zeros(n_values)
    std_D = np.zeros(n_values)
    log_data = np.log10(data)
    for i,n in enumerate(n_samples):
        distances = np.zeros(n_trials)
        for trial in range(n_trials):
            idx = np.random.choice(len(log_data), size=n, replace=False)
            sample = log_data[idx]

            distances[trial] = wasserstein_distance(sample, log_data)

        mean_D[i] = np.mean(distances)
        std_D[i] = np.std(distances)

    with open('convergence.txt', 'w') as file:
        for mean,std in zip(mean_D, std_D):
            file.write(f'{mean},{std}\n')
else:
    D = np.loadtxt('convergence.txt', delimiter=',')
    mean_D = D[:, 0]
    std_D  = D[:, 1]

plt.plot(n_samples, mean_D, marker='o', color=c1)
plt.fill_between(
    n_samples,
    mean_D - std_D,
    mean_D + std_D,
    alpha=0.3,
)
ref = mean_D[0] * (n_samples / n_samples[0])**(-0.5)
plt.plot(n_samples, ref, ls='--', color=c2, label=r'$n^{-1/2}$')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'Wasserstein distance in $\log_{10}(x)$', fontsize=fs)
# plt.ylabel('Log-space Wasserstein distance', fontsize=fs)
plt.xlabel('Sample size', fontsize=fs)
# plt.savefig('./figs/convergence_w_line_5.svg', dpi=500, bbox_inches='tight')

plt.show()


quantiles = [0.5, 0.9, 0.99]

q_mean = np.zeros((len(quantiles), len(n_samples)))
q_std  = np.zeros((len(quantiles), len(n_samples)))

recompute_quantiles = False

if recompute_quantiles:
    for i,n in enumerate(n_samples):
        q_trials = np.zeros((n_trials, len(quantiles)))

        for trial in range(n_trials):

            idx = np.random.choice(len(data), size=n, replace=False)
            sample = data[idx]

            q_trials[trial] = np.quantile(sample, quantiles)

        q_mean[:,i] = np.mean(q_trials, axis=0)
        q_std[:,i]  = np.std(q_trials, axis=0)

    with open('quantiles_mean.txt', 'w') as file:
        for i in range(len(quantiles)):
            for j in range(len(n_samples)-1):
                file.write(f'{q_mean[i,j]},')
            file.write(f'{q_mean[i,j+1]}\n')
    
    with open('quantiles_std.txt', 'w') as file:
        for i in range(len(quantiles)):
            for j in range(len(n_samples)-1):
                file.write(f'{q_std[i,j]},')
            file.write(f'{q_std[i,j+1]}\n')
else:
    q_mean = np.loadtxt('quantiles_mean.txt', delimiter=',')
    q_std = np.loadtxt('quantiles_std.txt', delimiter=',')

for i in range(len(quantiles)):
    plt.errorbar(n_samples, q_mean[i,:], yerr=q_std[i,:], marker='.', label=f'{quantiles[i]*100:.0f}%', capsize=4)

plt.xlabel('Sample size', fontsize=fs)
plt.ylabel('Quantile value', fontsize=fs)
plt.xscale('log')
plt.yscale('log')

plt.legend(title='Quantile', title_fontproperties={'weight' : 'bold'})
# plt.savefig('./figs/quantiles_5.svg', dpi=500)

plt.show()
