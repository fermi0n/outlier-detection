import sys,os
import numpy as np
from matplotlib import pyplot as plt
from glob import glob

class cal_source():
	def __init__(self):
		self.name = ''
		self.flux1 = []
                self.flux2 = []
                self.iono_amps = []
                self.iono_offsets = []
                self.max_flux = 0.0

def plot_source(cal_list,index,n_bands):
        c = cal_list[index]
        plt.clf()
        plt.title('Source: %s' % c.name)
        plt.subplot2grid((2, 2), (0, 0))
        plt.title('Iono Amps')
        for b in range(n_bands):
                plt.plot(c.iono_amps[int(b)])
        plt.subplot2grid((2, 2), (0, 1))
        plt.title('Iono Offsets')
        for b in range(n_bands):
                plt.plot(c.iono_offsets[int(b)])
        plt.subplot2grid((2, 2), (1, 0))
        plt.title('FLUX Values')
        plt.xlabel('Cadence Number')
        for b in range(n_bands):
                plt.plot(c.flux1[b])
		plt.show()
        #plt.savefig('rts_logs_%s.png' % c.name)


print len(sys.argv)
#logfiles = glob(sys.argv[1])
logfiles = sys.argv[1:]

if(len(logfiles) != 24):
        print 'Only Found %s Matching Log Files' % len(logfiles)

cal_names = []
cal_list = []
cal_index = 0
cal_dict = {}

found_nan = 0

for l in logfiles:
        log_file = open(l)
        band = int(l.split('node')[1][:3]) - 1
        for line in log_file:
                if(line.find('IONO') > 0):
                        entries = line.split()
                        cal_name = entries[3]
                        if(line.find('fit looks bad') > 0):
                                offset = float(entries[-5][:-1])
                                amp = float(entries[-1][:-1])
                                if(cal_name == 'J032750+023316'):
                                        print cal_name,offset,amp,band
                        else:
                                if(line.find('scaling amps') > 0):
                                        offset = float(entries[-6][1:])
                                        amp = float(entries[-1])
                                else:
                                        offset = float(entries[-2][1:])
                                        amp = 1.0

                        if(cal_name not in cal_names):
                                cal_names.append(cal_name)
                                cal = cal_source()
                                cal.index = cal_index
                                cal.name = cal_name
                                cal_dict[cal_name] = cal_index
                                cal_index += 1
                                cal.iono_offsets = [[] for ll in logfiles]
                                cal.iono_amps = [[] for ll in logfiles]
                                cal.flux1 = [[] for ll in logfiles]
                                cal.flux2 = [[] for ll in logfiles]
                                cal.iono_amps_prod = np.ones(len(logfiles))
                                cal.iono_offsets[band].append(offset)
                                cal.iono_amps[band].append(amp)
                                if(amp > 0.0):
                                        cal.iono_amps_prod[band] = amp
                                cal_list.append(cal)

                        else:
                                cal_list[cal_dict[cal_name]].iono_offsets[band].append(offset)
                                cal_list[cal_dict[cal_name]].iono_amps[band].append(amp)
                                if(amp > 0.0):
                                        cal_list[cal_dict[cal_name]].iono_amps_prod[band] *= amp

                        if(np.isnan(offset) or np.isnan(amp)):
                                if(found_nan == 0):
                                        print 'First Nan %s %s %s' % (cal_name, cal_dict[cal_name],band+1)
                                        found_nan = 1


                if(line.find('FLUX') > 0):
                        cal_name = (line.split())[3]
                        f2 = line.split('{')[1]
                        flux1 = (f2.split('}'))[0].split()[0]
                        flux2 = (f2.split('}'))[0].split()[1]
                        cal_list[cal_dict[cal_name]].flux1[band].append(float(flux1))
                        cal_list[cal_dict[cal_name]].flux2[band].append(float(flux2))

# prints out some eye-catching values

for c in cal_list:
        c.max_flux = max(np.ravel(c.flux1))
        c.max_amp = max(np.ravel(c.iono_amps))
        c.std_amp = np.std(np.ravel(c.iono_amps))
        c.med_flux = np.median(np.ravel(c.flux1))
        c.med_amp = np.median(np.ravel(c.iono_amps))
        c.max_amp_prod = max(c.iono_amps_prod)
#        for b in range(len(logfiles)):
#                if(c.iono_amps[b][0] < 0):
#                        print c.name,c.index,b,'iono amps 0 ',c.iono_amps[b][0]
        #        plt.plot(np.array(c.iono_amps[0]))
#        if(float(c.max_flux) > 10):
#                print c.name,c.index,'max FLUX',c.max_flux
#        if(c.max_amp > 50):
#                print c.name,c.index,'max_amp',c.max_amp

#        if(c.max_amp_prod > 10):
#                print c.name,c.index,'max amp prod', c.max_amp_prod
#        if(c.med_amp < 0.8):
#                print c.name,c.index,'med amp ', c.med_amp

# make a plot of all of the values for the 10th ranked source
# note that the first five sources have

n_bands = len(logfiles)
#c = cal_list[20]
#plt.clf()
#plt.title('Source: %s' % c.name)
#plt.subplot2grid((2, 2), (0, 0))
#plt.title('Iono Amps')
#for b in range(n_bands):
#        plt.plot(c.iono_amps[b])
#plt.subplot2grid((2, 2), (0, 1))
#plt.title('Iono Offsets')
#for b in range(n_bands):
#        plt.plot(c.iono_offsets[b])
#plt.subplot2grid((2, 2), (1, 0))
#plt.title('FLUX Values')
#plt.xlabel('Cadence Number')
#for b in range(n_bands):
#        plt.plot(c.flux1[b])

#plt.savefig('rts_logs_%s.png' % c.name)

#plot_source(cal_list,20,n_bands)

max_fluxes = [c.max_flux for c in cal_list]
max_amps = [c.max_amp for c in cal_list]
med_fluxes = [c.med_flux for c in cal_list]
med_amps = [c.med_amp for c in cal_list]
std_amps = [c.std_amp for c in cal_list]
