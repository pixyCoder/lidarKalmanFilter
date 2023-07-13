import numpy as np
from numpy.linalg import inv

def simulate_measurements():
    # Measurements are taken 1 time a second
    t = np.linspace(0, 30, num=30)
    numOfMeasurements = len(t)
    # Define x and y position and velocity components
    x = 100.0
    y = 100.0
    velx = 2.5
    vely = -2.5

    # Create storage arrays for true position data
    t_time = []
    t_x = []
    t_y = []
    t_r = []
    t_b = []
    # Compute the Real Position Data - This simulation use the cartesian
    # coordinate frame for computing new positional points and then converts
    # them to the polar coordinates.
    for i in range(0,numOfMeasurements):
        # Set the delta time to 1.0 second
        dT = 1.0
        # Store off the update time for this position update
        t_time.append(t[i])
        # Compute the new x and y position data with the assumption all of the
        # velocity motion is in the x direction i.e. linear motion
        x = x+dT*velx
        y = y+dT*vely
        # Store off the computed x and y data
        t_x.append(x)
        t_y.append(y)
        # Compute the sum of the squares of x and y as an intermediate step
        # before computing the range and storing it off
        temp = x*x + y*y
        r = np.sqrt(temp)
        t_r.append(r)
        # Compute the azimuth (or bearing) with the arctan2 function and convert
        # it to degrees. Then store this azimuth data
        b = np.arctan2(x, y) * 180/np.pi
        t_b.append(b)

    # Create storage containers for polar measurement data
    m_r = []
    m_b = []
    m_cov = []
    # Bearing standard deviation is in milliradians
    # Range standard deviation is in meters
    sig_b = 0.004*180/np.pi
    sig_r = 0.1 
    # Storage containers for cartesian measurements - for analysis purposes
    m_x = []
    m_y = []
    for ii in range(0, len(t_time)):
        # Compute the error for each measurement
        # By taking the max between .25 of the defined standard deviation and
        # the randomly generated normal error, it guarantees an error
        temp_sig_b = np.maximum(sig_b * np.random.randn(), 1.0*sig_b)
        temp_sig_r = np.maximum(sig_r * np.random.randn(), 1.0*sig_r)
        

        # Save the measurement values for bearing and range as a Function
        # of the true value + the error generated above
        temp_b = t_b[ii] + temp_sig_b
        temp_r = t_r[ii] + temp_sig_r
        # Save off the measurement data
        m_b.append(temp_b)
        m_r.append(temp_r)
        m_cov.append(np.array([[temp_sig_r*temp_sig_r, 0],
                      [0, temp_sig_b*temp_sig_b]]))
        m_x.append(temp_r*np.sin(temp_b*np.pi/180))
        m_y.append(temp_r*np.cos(temp_b*np.pi/180))

    return [m_r, m_b, m_cov, t_r, t_b, t_time, t_x, t_y, m_x, m_y]
    #        0    1     2      3    4    5      6   7    8    9
