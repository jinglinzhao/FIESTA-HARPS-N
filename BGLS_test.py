""" Make a dummy test and simple a figure for visualization """
import pylab as plt
from bgls import * 

# make signal
period = 60.
nperiods = 3
npoints = 1000
nsamples = 100

# lb pars
ofac = 200
plow = 1
phigh = 200

SNR = 2
t = np.linspace(0, nperiods * period, npoints)
y = np.sin(2 * np.pi / float(period) * t)

# take samples
s = np.sort(np.random.randint(0, len(t), nsamples))
ts = t[s]
ys = y[s]
yerr = ys * np.random.normal(0, 1. / SNR, len(ys))
ysnoise = ys + yerr

plt.figure()
periods, probs = bgls(ts, ysnoise, ys / SNR, plow=plow, phigh=phigh, ofac=ofac)

ax = plt.subplot(211)
# ax.plot(t, y, 'k-')
ax.errorbar(ts, ysnoise, yerr=ys / SNR, color='b', linestyle='None')
ax.plot(t, np.sin(2 * pi / float(periods[probs.argmax()]) * t), 'r-')
ax.set_xlabel('Time')
ax.set_ylabel('Signal')

ind = probs > 0.01
v = np.random.choice(periods[ind], 10, p=probs[ind] / sum(probs[ind]))
ax.plot(t, np.sin(2 * pi / v[:, None] * t).T, 'g-', alpha=0.3)

ax = plt.subplot(212)
ax.plot(periods, np.log(probs), 'k-')
# ax.vlines([period], 0, 1, color='b')
ax.set_xlabel('Periods')
ax.set_ylabel('Probabilities')
plt.show()