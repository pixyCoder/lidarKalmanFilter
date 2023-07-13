#!/usr/bin/python

import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
import simulate_measurements
import ekfilter


f_x = []
f_y = []
f_x_sig = []
f_y_sig =[]
f_xv = []
f_yv = []
f_xv_sig = []
f_yv_sig =[]

z = simulate_measurements.simulate_measurements()
for iii in range(0, len(z[0])):
  f = ekfilter.ekfilter(z, iii)
  f_x.append(f[0])
  f_y.append(f[1])
  f_xv.append(f[3])
  f_yv.append(f[4])
  f_x_sig.append(np.sqrt(f[2][0][0]))
  f_y_sig.append(np.sqrt(f[2][1][1]))


plot1 = plt.figure(figsize=(8, 6), dpi=80)
plt.grid(True)
plt.plot(z[5], z[3], linestyle='solid')
plt.plot(z[5], z[0], linestyle='dashed')
plt.title('True range and Measured range vs Update number', fontsize=14)
plt.ylabel('Range (m)', fontsize=18)
plt.xlabel('Time (seconds)', fontsize=18)
plt.savefig('TrueAndMeasuredRange.png')


plot2 = plt.figure(figsize=(8, 6), dpi=80)
plt.grid(True)
plt.plot(z[5], z[4], linestyle='solid')
plt.plot(z[5], z[1], linestyle='dashed')
plt.title('True azimuth and Measured azimuth vs Update number', fontsize=14)
plt.ylabel('Azimuth (degrees)', fontsize=18)
plt.xlabel('Time (seconds)', fontsize=18)
plt.savefig('TrueAndMeasuredAzimuth.png')


plot3 = plt.figure(figsize=(8, 6), dpi=80)
plt.grid(True)
plt.plot(z[5], f_xv, color="red", linestyle="dashed", label="Vx")
plt.plot(z[5], f_yv, color="blue", linestyle="dashed", label="Vy")
plt.title('Velocity estimate on each measurement update', fontsize=14, fontweight="bold")
plt.ylabel('Velocity (m/s)', fontsize=18)
plt.xlabel('Time (seconds)', fontsize=18)
plt.legend(fontsize=18)
plt.savefig('KalmanVelocityEstimates.png')






# Compute Range Error
e_x_err = []
e_x_3sig = []
e_x_3sig_neg = []
e_y_err = []
e_y_3sig = []
e_y_3sig_neg = []
for m in range(0, len(z[0])):
    e_x_err.append(f_x[m]-z[6][m])
    e_x_3sig.append(3*f_x_sig[m])
    e_x_3sig_neg.append(-3*f_x_sig[m])
    e_y_err.append(f_y[m]-z[7][m])
    e_y_3sig.append(3*f_y_sig[m])
    e_y_3sig_neg.append(-3*f_y_sig[m])
    print('ex, ey : ', e_x_err[m], ', ', e_y_err[m])   


plot4 = plt.figure(4), plt.grid(True)
line1 = plt.scatter(z[5], e_x_err)
line2, = plt.plot(z[5], e_x_3sig, color='green')
plt.plot(z[5], e_x_3sig_neg, color='green')
plt.ylabel('Position Error (meters)')
plt.xlabel('Time (seconds)')
plt.title('X Position localization errors \n', fontweight="bold")
plt.legend([line1, line2,], ['X Position Error', '3 Sigma Error Bound'])
plt.savefig('XpositionErrors.png')


plot5 = plt.figure(5), plt.grid(True)
yline1 = plt.scatter(z[5], e_y_err)
yline2, = plt.plot(z[5], e_y_3sig, color='green')
plt.plot(z[5], e_y_3sig_neg, color='green')
plt.ylabel('Position Error (meters)')
plt.xlabel('Time (seconds)')
plt.title('Y Position localization errors \n', fontweight="bold")
plt.legend([yline1, yline2,], ['Y Position Error', '3 Sigma Error Bound'])
plt.savefig('YpositionErrors.png')

plt.show()
