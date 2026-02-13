import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('data2.txt')

print(f'Mean : {np.mean(data)} \t Sigma : {np.std(data)}')

plt.figure()
plt.hist(data, bins=60)
plt.xlabel('Approximate $\lambda$')
plt.ylabel('Count')
plt.show()
