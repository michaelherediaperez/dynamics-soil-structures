import numpy as np
import matplotlib.pyplot as plt


def dmaclin():
    pass


def desplin(acc, Tmin, Tmax, DT, vxi, dt):
    """
    Calculation of the response spectra of a simple linear system by the linear 
    acceleration method.

    traduced by: 
        Michael Heredia Pérez
        Universidad Nacinal de Colombia sede Manizales

    written in Matlab: 
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

    Recommended data:
        load vinadelmar.txt
        desplin(vinadelmar, 0.05, 4, 0.05, [0.02 0.05], 0.005)
        
        load elcentro.txt
        desplin(elcentro, 0.05, 4, 0.05, [0.02 0.05], 0.02);
        
        load mexico.txt
        desplin(mexico, 0.05, 4, 0.05, [0.02 0.05], 0.02);
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

    # Memory reservoir.
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

    # For each ax in axs matrix:
    for r in range(2):
        for c in range(2):

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

            # Increment the counter
            count += 1
    
    plt.show()