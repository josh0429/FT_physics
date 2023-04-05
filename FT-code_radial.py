# Quick import cell
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'notebook'])
import scipy as sp
from scipy.fft import fftfreq, fftshift, ifftshift
from scipy.fft import fft, ifft, fft2, ifft2

N = 100 # Number of radial points
Qrange = 10 # Range in Q in units of GeV

# Define radial function here in units of GeV
def func(Q):
    A = 1.70
    B = 3.30
    a = [1, 2.33]
    b = [1, 14.72, 24.20, 84.1]
    tau = (Q**2)/(4*0.938**2)
    GEn = A*tau / (1 + B*tau) / (1 + Q**2 / 0.71)**2
    GMn = -1.91*(a[0] + a[1]*tau) / (b[0] + b[1]*tau + b[2]*tau**2 + b[3]*tau**3)
    return (GEn + tau*GMn) / (1 + tau)

def FTradial(N,Qrange):
    x = np.concatenate((np.linspace(0,Qrange,N,endpoint=False), np.linspace(-Qrange,0,N,endpoint=False)))
    y = np.concatenate((np.linspace(0,Qrange,N,endpoint=False), np.linspace(-Qrange,0,N,endpoint=False)))
    xvals, yvals = np.meshgrid(x, y)
    rho_Q = func(np.sqrt((xvals)**2 + (yvals)**2))
    rho_b = np.real(fftshift(np.diff(x)[0]*np.diff(y)[0]*fft2(rho_Q) / ((2*np.pi)**2)))
    b = fftshift(fftfreq(len(x), np.diff(x)[0] / (2*np.pi)))
    b = b / 5.068 # Converts to fm; comment out if you do not want unit conversion
    rho_b = rho_b * ((5.068)**2) # Converts to fm; comment out if you do not want unit conversion
    return np.transpose([b[N:2*N], rho_b[N,N:2*N]])

output = FTradial(N,Qrange)
np.savetxt('output.csv', output, delimiter=",")

plt.plot(output[:,0], output[:,1])
plt.xlabel(r'$b$ [fm]')
#plt.xlabel(r'$b$ [GeV]')
plt.ylabel(r'$ \rho(b)$ [fm$^{-2}$]')
#plt.ylabel(r'$ \rho(b)$ [GeV$^{-2}$]')
plt.show()