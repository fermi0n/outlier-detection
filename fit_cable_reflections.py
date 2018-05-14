import numpy as np
from numpy import polyfit,poly1d, pi, arange, zeros,exp, pi, array, std, mean, random,real,imag, shape
import matplotlib.pyplot as plt
from scipy import optimize

# pads missing fine channels with zeros/gaussian noise for FFTs

def pad_edge_channels(bp_data,zero_padded=False):

    sigma = std(bp_data)
    avg = mean(bp_data)

    flagged_channels = [0,1,16,30,31]

    chan_count = 0
    data_count = 0
    padded_bp = zeros(24 * 32)
    for i in range(24):
        for j in range(32):
            if(j in flagged_channels):
                if(zero_padded):
                    padded_bp[data_count] = avg 
                else:
                    padded_bp[data_count] = avg + sigma * random.rand()
            else:
                padded_bp[data_count] = bp_data[chan_count]
                chan_count +=1
            data_count += 1

    return padded_bp

def clip_edge_channels(bp_data,clip_width=1):

    chan_count = 0
    data_count = 0
    clipped_bp = zeros(24 * (27-2*clip_width),dtype=complex)
    clipped_x = zeros(24 * (27-2*clip_width))
    
    channels = range(32)
    channels.pop(16)
    channels = channels[2+clip_width:-(2+clip_width)]

    for i in range(24):
        for j in range(32):
            if(j in channels):
                clipped_bp[data_count] = bp_data[chan_count]
                clipped_x[data_count] = chan_count
                data_count +=1
            chan_count += 1

    return clipped_bp, clipped_x

def bp_xvals():
    # returns list of unflagged channel numbers
    flagged_channels = [0,1,16,30,31]
    chan_count = 0
    xvals = []
    for i in range(24):
        for j in range(32):
            if(j not in flagged_channels):
                xvals.append(chan_count)
            chan_count +=1 
        
    return array(xvals)

def fit_reflection(data):

    sigma = np.std(data)
    avg = np.mean(data)
    phase = pi/2.0

    padded = pad_edge_channels(data)

    freq = 1.0 / (np.fft.fftfreq(len(padded))[np.argmax(abs(np.fft.fft(padded))[1:len(padded)/2])+1])

#    print freq
#    freq = 20.0
    #freq = 31.0
    freq = 24.0

    clipped, clipped_x = clip_edge_channels(padded,clip_width=4)

    sigma = np.std(clipped)
    avg = np.mean(clipped)

    #tX = data
    #Tx = arange(len(data))
    tX = clipped
    Tx = clipped_x

    p0 = [sigma,freq,phase,avg]
#    fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/p[1]*x+p[2]) + p[3] # Target function
    fitfunc = lambda p, x: (p[0]*sigma)*np.cos(2*np.pi/(p[1]*freq)*x+(p[2]*phase)) + (p[3]*avg)
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    
#    print p0

    initial_vals = [1.0,1.0,0.0,1.0]

    p1, success = optimize.leastsq(errfunc, initial_vals, args=(Tx, tX))

#    print p1

    model_vals = [i * p for i,p in zip(p0,p1)]

#    print model_vals

    time = np.linspace(Tx.min(), Tx.max(), 1000)
#    plt.plot(Tx, tX, "ro", time, fitfunc(p1, time), "r-")
#    plt.plot(Tx, tX, "ro", time, fitfunc(p1, time), "r-")
#    plt.plot(time, fitfunc(p0, time), "b-")

    model_points = fitfunc(p1,bp_xvals())
    model_points = fitfunc(p0,bp_xvals())

    return model_points, model_vals

def fit_bandpass_and_reflection(bp_data):

    if(len(shape(bp_data)) == 1):

        xvals = arange(len(bp_data))
        
        clipped, clipped_x = clip_edge_channels(bp_data,clip_width=0)
    
        z_real = polyfit(clipped_x,real(clipped),3)
        z_imag = polyfit(clipped_x,imag(clipped),3)

        p_real = poly1d(z_real)
        p_imag = poly1d(z_imag)

        return p_real(xvals), p_imag(xvals)
    
    else:
        if(len(shape(bp_data)) == 3 and shape(bp_data)[0] == 2 and shape(bp_data)[1] == 2):
            xvals = arange(shape(bp_data)[2])
            models_out = [[[],[]],[[],[]]]
            for x in range(2):
                for y in range(2):
                    clipped, clipped_x = clip_edge_channels(bp_data[x][y],clip_width=0)  
                    z_real = polyfit(clipped_x,real(clipped),3)
                    z_imag = polyfit(clipped_x,imag(clipped),3)

                    p_real = poly1d(z_real)
                    p_imag = poly1d(z_imag)   
                    
                    models_out[x][y] = p_real(xvals) + 1.0j*p_imag(xvals) 

            return models_out

        else:
            print 'Bandpass data of unsupported dimensions: %s' % shape(bp_data)
            exit(1)
    
