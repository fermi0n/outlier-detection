from numpy import polyfit,poly1d,arange,ravel,array

def fit_fsum(fsum):

    n_bands = 24
    freqs = arange(27)
    freqs[14:] += 1
    fit_sum = []
    
    for i in range(24):
        y = fsum[i*27+2:(i+1)*27-2]
        z = polyfit(freqs[2:-2],y,1)
        p = poly1d(z)
        fit_sum.append(p(freqs))

    return ravel(fit_sum)

def fit_bandpass(bp_data):

    freqs = arange(len(bp_data))
    y = array(bp_data)
    z = polyfit(freqs,bp_data,2)
    p = poly1d(z)

    return p(freqs)

def fit_sigma(bp_data,model):

    sigma = sum((bp_data - model)**2) / (float) (len(bp_data))
    return sigma




