import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from random import uniform

global tk_prev, vn_prev, J_prod, V_es
tk_prev = 0
vn_prev = 0


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
frq = float(input('Please insert the frequency time?'))

# function that returns dy/dt
def model(x, t):
    C = x[0]
    E = x[1]
    IP3R = x[2]
    IP3 = x[3]

    global tk_prev, vn_prev, J_prod

    # Parameters
    param = {'k0': 0.03, 'k1': 0.0004, 'k2': 0.2, 'k3': 0.5, 'k5': 0.25, 'k6': 4, 'k7': 0.08,
             'k9': 0.08, 'v7': 0.02, 'Kip3': 0.3, 'Ka': 0.2, 'Ki': 0.2, 'Kca': 0.3, 'beta': 35,
             'Hcce': 10, 'kcce': 0.01, 'kir': 0.08, 'Hir': 0.9, 'kmr': 0.5, 'Kd': 10, 'ATP': 0.01,
             'O_delta': 0.15, 'tk': 0.3, 'vn': 30}

    # Terms
    V_LM = param['k0']
    V_CCE = (param['kcce'] * param['Hcce'] ** 2) / (param['Hcce'] ** 2 + E ** 2)
    V_IR = param['kir'] * (param['ATP'] ** 1.4 / (param['ATP'] ** 1.4 + param['Hir']))
    V_OUT = param['k5'] * C
    V_ER_leak = param['k1'] * (E - C)
    V_ER_rel = ((param['k2'] * IP3R * C ** 2 * IP3 ** 2) / (
                (param['Ka'] ** 2 + C ** 2) * (param['Kip3'] ** 2 + IP3 ** 2))) * (E - C)
    V_SERCA = param['k3'] * C
    V_IP3R_Rec = (param['k6'] * param['Ki'] ** 2) / (param['Ki'] ** 2 + C ** 2)
    V_IP3R_Inact = param['k6'] * IP3R
    V_PLC_beta = param['kmr'] * (param['ATP'] / (param['Kd'] + param['ATP']))
    V_PLC_delta = (param['v7'] * C ** 2) / (param['Kca'] ** 2 + C ** 2)
    V_IP3_Deg = param['k9'] * IP3
    # est


    if (125 <= t % frq <150):
        a = (t % frq) -125
        V_es = ((p / 20) * np.exp((a) / 10))/1000

    elif (150 <= t % frq <= 215):
        a = (t % frq) -125
        V_es = (4.5 * p * np.exp(-((a) / 20)))/1000

    else :
        V_es = 0

    #print (t, V_es)


    # ODE IP3 Terms
    if (t == 0):
        J_prod = (param['O_delta'] * uniform(0, 1)) + (param['O_delta'] * uniform(0, 1))
    if (t - tk_prev >= param['tk']):
        tk_prev = t
        J_prod = (param['O_delta'] * uniform(0, 1))
    if (t - vn_prev >= param['vn']):
        vn_prev = t
        J_prod = (param['O_delta'] * uniform(0, 1))

    # C - ODE Citosolic Calcium Concentration
    dCdt = V_LM + V_CCE + V_IR - V_OUT + V_ER_leak + V_ER_rel - V_SERCA + V_es

    # h - ODE further one for Ca 2+ -mediated deinactivation
    dEdt = param['beta'] * (V_SERCA - V_ER_leak - V_ER_rel)

    # I - ODE intracellular (cytosolic) IP3
    dIP3Rdt = V_IP3R_Rec - V_IP3R_Inact
    dIP3dt = V_PLC_beta + V_PLC_delta - V_IP3_Deg
    # dIP3dt = J_prod - V_IP3_Deg

    return [dCdt, dEdt, dIP3Rdt, dIP3dt]


# initial condition
C0 = 0.1  # 0.06
IP30 = 0.1  # 0.0096
IP3R0 = 0.1  # 0.9174e-6
E0 = 1.5  # 72

x0 = [C0, E0, IP3R0, IP30]

# time points
t = np.arange(0, 1200, 0.01)

x = odeint(model, x0, t)

C = x[:, 0]
print(type(C))

E = x[:, 1]
IP3R = x[:, 2]
IP3 = x[:, 3]

fig = plt.figure()

plot_C, = plt.plot(t, C, 'r', label="plot_C")
#plot_h, = plt.plot(t, E, 'g', label="plot_h")
plot_IP3, = plt.plot(t, IP3, 'b', label="plot_I")
plt.legend([plot_C, plot_IP3], ['Cytosol', 'IP3'])
plt.legend([plot_C , plot_IP3], ['Cytosol', 'IP3'])
plt.axis([0, np.max(t), 0, np.max(C) * 1.2])
plt.title('Alzheimer Model + Estradiol (0.2345 uM)')
plt.xlabel('Time (s)')
plt.ylabel('Calcium Concentration (uM)')
plt.show()


from scipy import signal
freqs, times, spectrogram = signal.spectrogram(C)
plt.pcolormesh( times/100,freqs, spectrogram, shading='gouraud')
plt.colorbar()
plt.title('Spectrogram for AD cell + Estradiol (0.2345 uM)')
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

