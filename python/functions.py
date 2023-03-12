import numpy as np
import matplotlib.pyplot as plt


def dmaclin(p, m, w, xi, dt):
    """
    Calculation of the response of a simple system by the linear acceleration 
    method.

    traduced by:
        Michael Heredia Pérez
        Universidad Nacional de Colombia

    written in matlab by:
        Jorge E. Hurtado G.
        Universidad Nacional de Colombia
    
    Input:
        p  : columns vector of external load.
        m  : system mass.
        w  : natural frequency of the system.
        xi : viscous damping fraction.
        dt : time step

    Output:
        t : time vector
        d : response displacement
        v : response velocity
        a : response acceleration

    Tested data:
        load tokachikenoki.txt
        dmaclin(-tokachikenoki,1,2*pi,0.05,0.02)

        load elcentro.txt
        dmaclin(-elcentro,1,2*pi,0.05,0.02)

        load mexico.txt
        dmaclin(-mexico,1,2*pi,0.05,0.02)

        load vinadelmar.txt
        dmaclin(-vinadelmar,1,2*pi,0.05,0.005)
    """

    n    = p.shape[0]               # Number of points
    tmax = dt*n                     # max time
    t    = np.linspace(0, tmax, n)  # time domain
    d0, v0, a0 = 0, 0, 0            # initial conditions for x, v=dot(x), a=dot(v)

    k = m*w**2                      # Stiffnes
    c = 2*m*w*xi                    # viscosity coeff.
    kbar = 3 + 3*c/dt + 6*m/(dt**2) 
    ikbar = 1/kbar  

    # Memory reserved
    d, v, a = np.zeros(n), np.zeros(n), np.zeros(n)

    for i in range(n):
        p1   = p[i,:]
        dp   = m*(6*d0/dt**2 + 6*v0/dt + 2*a0)
        dp  += c*(3*d0/dt + 2*v0 + dt*a0/2)
        pbar = p1 + dp
        
        d1 = ikbar*pbar
        v1 = 3*(d1-d0)/dt-2*v0-dt*a0/2
        a1 = 6*(d1-d0)/dt^2 - 6*v0/dt - 2*a0
        
        d[i] = d1; d0 = d1
        v[i] = v1; v0 = v1
        a[i] = a1; a0 = a1
        

    # Plotable data
    data = [-p/m, d, v, a]
    info = ["Soil acceleration", "Displacement", "Velocity", "Acceleration"]
    
    # Plotting
    fig, axs = plt.subplots(2, 2, layout = "constrained")
    #gridspec = axs[0, 0].get_subplotspec().get_gridspec()
    
    count = 0

    for r in range(2):  
        for c in range(2): # For each ax in axs matrix:

            ax = axs[r, c]
            ax.plot(t, data[count])
            ax.grid(b=True, which='major', linestyle='-')
            ax.set_xlabel("Time")
            ax.set_ylabel(info[count])
            
            count += 1
    
    plt.show()





def desplin(acc, Tmin, Tmax, DT, vxi, dt):
    """
    Calculation of the response spectra of a simple linear system by the linear 
    acceleration method.

    traduced by: 
        Michael Heredia Pérez
        Universidad Nacinal de Colombia

    written in Matlab by: 
        Jorge E. Hurtado G.
        Universidad Nacional de Colombia
    
    Input:
        acc  : column vector (array) soil acceleration 
        Tmin : minimum calculation period
        Tmax : maximun calculation period
        DT   : period increment
        vxi  : vector containing the viscous damping fractions for which the 
                spectra are to be calculated
        dt   : time step of the accelerogram

    Output:
        Sd: displacement spectrum
        Sv: velocity spectrum
        Sa: acceleration spectrum

    Tested data:
        load vinadelmar.txt
        desplin(vinadelmar, 0.05, 4, 0.05, [0.02 0.05], 0.005)
        
        load elcentro.txt
        desplin(elcentro, 0.05, 4, 0.05, [0.02 0.05], 0.02)
        
        load mexico.txt
        desplin(mexico, 0.05, 4, 0.05, [0.02 0.05], 0.02)
    """

    # Number of points in the data
    nacc = acc.shape[0]                    # Accelerations
    nvxi = len(vxi)                        # Vicous damping
    
    # Number of instants to be calculated.
    m = (Tmax - Tmin)/DT + 1

    # range of periods from the given limits
    T = np.linspace(Tmin, Tmax, m)
    
    # Frecuencies all over the periods.
    W = 2*np.pi / T

    # Memory reserved
    Sd = np.zeros((nvxi, m))
    Sv = np.zeros((nvxi, m))
    Sa = np.zeros((nvxi, m))

    for i in range(nvxi):
        # For each viscous damping given...
        xi = vxi[i]
        for j in range(m):
            # For each instant ...
            w = W[j]
            d, v, a = dmaclin(-acc, 1, w, xi, dt)
            
            Sd[i, j] += np.max(abs(d))
            Sd[i, j] += np.max(abs(v))
            Sd[i, j] += np.max(abs(acc + a))      
    
    t = np.linspace(1, dt*nacc, nacc)

    # Plotable data
    data = [acc, Sd, Sv, Sa]
    info = ["Acceleration", "Displ. spectrum", "Vel. spectrum", "Acc. spectrum"]
    labl = [r"$\xi = 0.02$", r"$\xi = 0.05$"]

    # Plotting
    fig, axs = plt.subplots(2, 2, layout = "constrained")
    #gridspec = axs[0, 0].get_subplotspec().get_gridspec()
    
    count = 0

    for r in range(2):  
        for c in range(2): # For each ax in axs matrix:

            ax = axs[r, c]
            for vx in range(nvxi):
                ax.plot(t, data[count][vx], label=labl[count])
            ax.label(loc=1)
            ax.grid(b=True, which='major', linestyle='-')
            if count == 0:
                ax.set_xlabel("Time")
            else:
                ax.set_xlabel("Period")
            ax.set_ylabel(info[count])
            
            count += 1
    
    plt.show()