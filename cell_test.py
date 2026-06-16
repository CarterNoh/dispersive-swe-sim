import matplotlib.pyplot as plt
import numpy as np

cellsize = np.array([0.1, 0.15, 0.25, 0.5, 0.75, 1])
d0 = np.array([1.696, 3.965, 11.252, 45.60, 102.6, 182.5]) # max stable depth for grid = 512 with disturbance of 0.1, timestep 1/60
d1 = np.array([1.712, 3.972, 11.260, 45.61, 102.8, 182.9]) # max stable depth for grid = 256 with disturbance of 0.1, timestep 1/60
d_avg = np.mean((d0, d1), axis=0)
A = np.vstack((cellsize**2, cellsize, np.ones(cellsize.shape))).T
[b2, b1, b0] = np.linalg.pinv(A) @ d_avg
# print(b2, b1, b0)
d_est = b2*cellsize**2 + b1*cellsize + b0
# plt.plot(cellsize, d0)
# plt.plot(cellsize, d1)
# plt.plot(cellsize, d_est)
# plt.legend(("512", "256", "Estimate"))
# plt.show()
# Result is highly quadratic, a very good fit!

def solve_quadratic(d):
    a = b2 # 182.80027907467993 
    b = b1 # 0.045464332332812774
    c = b0 #-0.14717654147795045
    return (-b + np.sqrt(b**2 - 4*a*(c-d))) / (2 * a)

print(solve_quadratic(12))

# Note: stability also seems to be somewhat dependant on the timestep, but I didn't have time to do more testing. 

