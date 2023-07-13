import numpy as np
from numpy.linalg import inv

def ekfilter(z, updateNumber): # z = [r, b]
    dt = 1.0
    j = updateNumber
    # Initialize State
    if updateNumber == 0: # First Update
        # compute position values from measurements
        temp_x = z[0][j]*np.sin(z[1][j]*np.pi/180) # x = r*sin(b)
        temp_y = z[0][j]*np.cos(z[1][j]*np.pi/180) # y = r*cos(b)
        # State vector - initialize position values
        ekfilter.x = np.array([[temp_x],
                            [temp_y],
                            [0],
                            [0]])
        # State covariance matrix - initialized to zero for first update
        ekfilter.P = np.array([[0, 0, 0, 0],
                                 [0, 0, 0, 0],
                                 [0, 0, 0, 0],
                                 [0, 0, 0, 0]])
        # State transistion matrix - linear extrapolation assuming constant velocity
        ekfilter.A = np.array([[1, 0, dt, 0],
                             [0, 1, 0, dt],
                             [0, 0, 1, 0],
                             [0, 0, 0, 1]])
        # Measurement covariance matrix
        ekfilter.R = z[2][j]
        # System error matrix - initialized to zero matrix for first update
        ekfilter.Q = np.array([[0, 0, 0, 0],
                                 [0, 0, 0, 0],
                                 [0, 0, 0, 0],
                                 [0, 0, 0, 0]])
        # Residual and kalman gain
        # not computed for first update but initialized so it could be output
        residual = np.array([[0, 0],
                      [0, 0]])
        K = np.array([[0, 0],
                      [0, 0],
                      [0, 0],
                      [0, 0]])

    # Reinitialize State
    if updateNumber == 1: # Second Update
        # Get previous state vector values
        prev_x = ekfilter.x[0][0]
        prev_y = ekfilter.x[1][0]
        temp_x = z[0][j]*np.sin(z[1][j]*np.pi/180) # x = r*sin(b)
        temp_y = z[0][j]*np.cos(z[1][j]*np.pi/180) # y = r*cos(b)
        #  Compute velocity - vel = (pos2 - pos1)/deltaTime
        temp_xv = (temp_x - prev_x)/dt
        temp_yv = (temp_y - prev_y)/dt
        # State vector - reinitialized with new position and computed velocity
        ekfilter.x = np.array([[temp_x],
                            [temp_y],
                            [temp_xv],
                            [temp_yv]])
        # state covariance matrix - initialized to large values
        ekfilter.P = np.array([[1, 0, 0, 0],
                                 [0, 1, 0, 0],
                                 [0, 0, 1, 0],
                                 [0, 0, 0, 1]])
        # State transistion matrix - linear extrapolation assuming constant velocity
        ekfilter.A = np.array([[1, 0, dt, 0],
                             [0, 1, 0, dt],
                             [0, 0, 1, 0],
                             [0, 0, 0, 1]])
        # Measurement covariance matrix - provided by the measurment source
        ekfilter.R = z[2][j]
        # System error matrix
        # adds 4 m std dev in position and 20 m/s std dev in velocity
        ekfilter.Q = np.array([[0.1, 0, 0, 0],
                                 [0, 0.1, 0, 0],
                                 [0, 0, 0.025, 0],
                                 [0, 0, 0, 0.025]])
        # Residual and kalman gain-  initialized so it could be output
        residual = np.array([[0, 0],
                      [0, 0]])
        K = np.array([[0, 0],
                      [0, 0],
                      [0, 0],
                      [0, 0]])

    if updateNumber > 1: # Third Update and Subsequent Updates
      # Predict state and state covariance forward
      x_prime = ekfilter.A.dot(ekfilter.x)
      P_prime = ekfilter.A.dot(ekfilter.P).dot(ekfilter.A.T) + ekfilter.Q
      # Form state to measurement transition matrix
      x1 = x_prime[0][0]
      y1 = x_prime[1][0]
      x_sq = x1*x1
      y_sq = y1*y1
      den = x_sq+y_sq
      den1 = np.sqrt(den)
      ekfilter.H = np.array([[  x1/den1,    y1/den1, 0, 0],
                           [y1/den, -x1/den, 0, 0]])
      ekfilter.HT = np.array([[x1/den1, y1/den],
                              [y1/den1, -x1/den],
                              [0, 0],
                              [0, 0]])
      # Measurement covariance matrix
      ekfilter.R = z[2][j]
      # Compute Kalman Gain
      S = ekfilter.H.dot(P_prime).dot(ekfilter.HT) + ekfilter.R
      K = P_prime.dot(ekfilter.HT).dot(np.linalg.inv(S))
      # Estimate State
      temp_z = np.array([[z[0][j]],
                         [z[1][j]]])
      # Convert the predicted cartesian state to polar range and azimuth
      pred_x = x_prime[0][0]
      pred_y = x_prime[1][0]
      sumSquares = pred_x*pred_x + pred_y*pred_y
      pred_r = np.sqrt(sumSquares)
      pred_b = np.arctan2(pred_x, pred_y) * 180/np.pi
      h_small = np.array([[pred_r],
                       [pred_b]])
      # Compute residual difference between state and measurement for data time
      residual = temp_z - h_small
      # Compute new estimate for state vector using the Kalman Gain
      ekfilter.x = x_prime + K.dot(residual)
      # Compute new estimate for state covariance using the Kalman Gain
      ekfilter.P = P_prime - K.dot(ekfilter.H).dot(P_prime)
    return [ekfilter.x[0], ekfilter.x[1], ekfilter.P, ekfilter.x[2], ekfilter.x[3], K, residual];
#               0             1              2           3               4          5     6
