import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib as mpl

qtd = 0
global tk_prev, vn_prev, J_prod, V_es
def findpeak(est):
    c = 410
    q = 2.8
    a = 0.585
    b = -0.092

    if (0 <= est < 1):
        peak = (c * (1 - ((1 - est) ** (q - 2.4))) ** (1 / (q - 0.07)))
        print("0-1")

    elif (1 <= est <= 10):
        peak = (c * (1 / (a * (np.exp(np.log(est) * b)))))
        print("1-10")

    else:
        print("wrong number")

    return (peak)


est = float(input('Please insert the does of estradiol?'))
p = findpeak(est)

# function that returns dy/dt
def model(x, t):
	C = x[0]
	h = x[1]
	I = x[2]

	# Parameters
	param = {'sigma0': 0.05, 'ko': 0.5, 'kf': 0.5, 'vm2': 15, 'k2': 0.1, 'n': 2.02, 'kc1': 0.15,
	'kc2': 0.15, 'ki': 0.1, 'vm3': 40, 'm': 2.2, 'Odelta': 0.15, 'Kdelta': 0.5, 'alpha': 0.8,
	'O3K': 0.01, 'K3K': 1.0, 'Omega_5P': 0.01, 'vp': 0.05, 'kp': 0.3, 'kd': 0.08}

	# Terms
	sigma1 = 4 * param['vm3'] * ( (param['kc1']**param['n'] * C**param['n']) / ( (C**param['n'] + param['kc1']**param['n']) * (C**param['n'] + param['kc2']**param['n']) ) ) * (I**param['m'] / (param['ki']**param['m'] + I**param['m']) ) * (h - C)
	sigma2 = param['vm2'] * (C**2 / (param['k2']**2 + C**2) )
	sigma3 = param['vp'] * (C**2 / (param['kp']**2 + C**2))
	#J_delta = param['Odelta'] * ((C / (C + param['Kdelta'])) + ((1 - param['alpha']) * (param['Kdelta'] / (param['Kdelta'] + C) ) ) )
	#J_3K = param['O3K'] / param['K3K'] * I
	#J_5P = param['Omega_5P'] * I
	
	# est
	if (t < 30):
		V_es = ((p / 20) * np.exp(t / 10)) / 1000
	if (30 <= t < 90):
		V_es = (4.5 * p * np.exp(-(t / 20))) / 1000
	if (t >= 90):
		V_es = 0
	# C - ODE Citosolic Calcium Concentration
	dCdt = param['sigma0'] - (param['ko'] * C) + sigma1 - sigma2 + (param['kf'] * (h - C)) + V_es

	# h - ODE further one for Ca 2+ -mediated deinactivation
	dhdt = sigma2 - sigma1 - (param['kf'] * (h - C))

	# I - ODE intracellular (cytosolic) IP3
	dIdt = sigma3 - (param['kd'] * I) #J_delta - J_3K - J_5P

	return [dCdt, dhdt, dIdt]

# initial condition
C0 = 0.1 # ca
I0 = 0.1 # IP3
h0 = 1.5 # endoplasmic reticulum

x0 = [C0, h0, I0] # vector ca in sytosol endoplasmic reticulum IP3

# time points (time slap)
t = np.arange(0,1200,0.01)

x = odeint(model, x0, t)

#obtaining the results from x
C = x[:,0]
h = x[:,1]
I = x[:,2]

fig = plt.figure()

plot_C, = plt.plot(t, C, 'r', label="plot_C")
#plot_h, = plt.plot(t, h, 'g', label="plot_h")
plot_I, = plt.plot(t, I, 'b', label="plot_I")
plt.legend([plot_C, plot_I], ['Cytosol', 'IP3'])
#plt.legend([plot_C], ['Cytosol'])
#plt.axis([0, 600, 0, np.max(C) * 1.2])
plt.title('Healthy Model')
plt.xlabel('Time (s)')
plt.ylabel('Calcium Concentration (uM)')




plt.show()



from scipy import signal
freqs, times, spectrogram = signal.spectrogram(C)
plt.pcolormesh( times/100,freqs, spectrogram, shading='gouraud')
plt.colorbar()
plt.title('Spectrogram for Healthy Cell')
plt.ylabel('Frequency band')
plt.xlabel('Time window')
plt.tight_layout()
plt.show()

plt.figure(figsize=(5, 2))
plt.imshow(spectrogram, aspect='auto', cmap='hot_r', origin='lower')
plt.colorbar()
plt.title('Spectrogram')
plt.ylabel('Frequency band')
plt.xlabel('Time window')
plt.tight_layout()
plt.show()
print('Cytosol:', C[-1])
print('RE:', h[-1])
