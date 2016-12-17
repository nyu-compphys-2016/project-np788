import math, pylab, random
import numpy as np
import matplotlib.pyplot as plt

J = 1.0  
L = 100
nsweep = 1000

N = L**2
kT = np.arange(3.269, 3.79)
e = np.zeros(len(kT))
m = np.zeros(len(kT))
c = np.zeros(len(kT))

# initial data
s = np.ones((L, L))
E = 0.0
M = 0.0
em = []
for i in range(L):
    for j in range(L):
        E -= J*s[i,j]*(s[(i+1)%L,j]+s[i,(j+1)%L])
        M += s[i,j]

for t in range(len(kT)):

    for sweep in range(nsweep):
        print("Now on sweep ",t, sweep) 


        # sweep over all particles in lattice
        for i in range(L):
            for j in range(L):

                # compute energy required to flip spin
                dE = s[(i+1)%L,j]+s[(i-1)%L,j]+s[i,(j+1)%L]+s[i,(j-1)%L]
                dE *= 2.0*J*s[i,j]

                # Metropolis algorithm to see if we should accept trial
                if dE <= 0.0 or random.random() <= math.exp(-dE/kT[t]):
                    s[i,j] *= -1
                    M += 2.0*s[i,j]
                    E += dE

        deltae = E-e[t]
        deltam = M-m[t]
        e[t] += deltae/(sweep+1)
        m[t] += deltam/(sweep+1)
        c[t] += deltae*(E-e[t])
        em.append(m[-1])
    e[t] /= N
    m[t] /= N
    c[t] /= nsweep*N*kT[t]**2
plt.figure()
plt.plot(np.arange(0,nsweep), em)
plt.xlabel('Number of sweeps')
plt.ylabel('Magnetization ($J/\mu$)')
plt.title('Magnetization as a function of sweeps at T = 3.269')
plt.show()

