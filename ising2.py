import math, pylab, random
import numpy as np
import matplotlib.pyplot as plt

J = 1.0  
L = 100
nsweep = 1000
#seed = input('random number seed -> ')
#random.seed(seed)
N = L**2
kT = np.arange(0.1, 5.0, 0.05)
e = np.zeros(len(kT))
m = np.zeros(len(kT))
c = np.zeros(len(kT))

# initial data
s = np.ones((L, L))
E = 0.0
M = 0.0
for i in range(L):
    for j in range(L):
        E -= J*s[i,j]*(s[(i+1)%L,j]+s[i,(j+1)%L])
        M += s[i,j]

### prepare animated plot
##pylab.ion()
##image = plt.imshow(s, vmax=1, vmin=-1)

# slowly warm up
for t in range(len(kT)):
    pylab.ion()
    image = plt.pcolormesh(s, cmap='bwr', vmax=1, vmin=-1)
    # average nsweep sweeps
    for sweep in range(nsweep):
        print("Now on sweep ",t, sweep) 
        # update animated plot
        #image.set_data(s)
        plt.title('Temperature = %g'%kT[t])
        plt.draw()
        plt.show()

        # sweep over all particles in lattice
        for i in range(L):
            for j in range(L):

                # compute energy required to flip spin
                dE = s[(i+1)%L,j]+s[(i-1)%L,j]+s[i,(j+1)%L]+s[i,(j-1)%L]
                dE *= 2.0*J*s[i,j]

                # Metropolis algorithm to see if we should accept trial
                if dE <= 0.0 or random.random() <= math.exp(-dE/kT[t]):
                    # accept trial: reverse spin; return dE and dM
                    s[i,j] *= -1
                    M += 2.0*s[i,j]
                    E += dE

        # update running means and variances
        deltae = E-e[t]
        deltam = M-m[t]
        e[t] += deltae/(sweep+1)
        m[t] += deltam/(sweep+1)
        c[t] += deltae*(E-e[t])
    e[t] /= N
    m[t] /= N
    c[t] /= nsweep*N*kT[t]**2
    pylab.ioff()
    plt.savefig("frame" + str(t) + ".png")
# produce plots
#pylab.ioff()
#plt.show()
plt.figure()
plt.plot((2.269, 2.269), (-2.0, -0.4), label = '$T_C = 2.269$')
plt.plot(kT, e, 'o', label = 'Energy')
plt.title('Energy as temp increases')
plt.xlabel('Temperature ($J/k_B$)')
plt.ylabel('Energy per atom ($J$)')
plt.legend()
plt.grid()
plt.figure()
plt.plot((2.269, 2.269), (-0.2, 1.0), label = '$T_C = 2.269$')
plt.plot(kT, m, 'o', label = 'Magnetization')
plt.title('Magnetization as temp increases')
plt.xlabel('Temperature ($J/k_B$)')
plt.ylabel('Magnetization per atom ($J/\mu$)')
plt.legend()
plt.grid()
plt.figure()
plt.plot((2.269, 2.269), (0.0, 4.0), label = '$T_C = 2.269$')
plt.plot(kT, c, 'o', label = 'Specific heat')
plt.xlabel('Temperature ($J/k_B$)')
plt.ylabel('Specific heat per atom ($J/k_B^2$)')
plt.legend()
plt.grid()
plt.show()
