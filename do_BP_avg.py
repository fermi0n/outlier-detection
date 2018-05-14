#from pylab import clf, plot, subplot, text, figure, title, close,savefig, matplotlib
#matplotlib.use('Agg')
import RTS_cals
from astropy.io import fits
from optparse import OptionParser,OptionGroup
from matplotlib import pyplot as plt

from numpy import conj, sqrt,zeros, fft, array,std, mean, polyfit, poly1d, arange, real, imag, median

import glob
from collections import Counter
import os,sys
from fit_cable_reflections import fit_bandpass_and_reflection, bp_xvals

usage = 'Usage: do_BP_fits.py [text file of obsIDs] [metafits]'

parser = OptionParser(usage=usage)

parser.add_option('--goodlist',dest='goodlist',type='string',default='',help='List of good obsids to be processed')

parser.add_option('--tag',dest='tag',type='string',default=None,help='Label to tag output plots')

(options, args) = parser.parse_args()

#obs_list = glob.glob('%s/1*' % args[0])
obs_list = open(args[0]).readlines()
obs_list = [l.split()[0] for l in obs_list]
meta_file = fits.open(args[1])

mwa_dir = os.getenv('MWA_DIR','/astro/mwaeor/MWA/')

cent_chan = meta_file[0].header['CENTCHAN']
rts_inputs = meta_file[1].data['Input']
cables = meta_file[1].data['Flavors']
tilenames = meta_file[1].data['Tile']
cable_types = set(cables)

cable_counts = dict(Counter(cables[::2])) 
rts_ants = rts_inputs / 2
rts2cables = dict(zip(rts_ants[::2],cables[::2]))
rts2tiles = dict(zip(rts_ants[::2],tilenames[::2]))

if(options.tag is None):
    tag = 'avg_cals'
else:
    tag = options.tag

if('EDA' in cable_types):
    cable_types.remove('EDA')

band_x = arange(32)

all_cals = []
med_cal = RTS_cals.rts_cal()
avg_cal = RTS_cals.rts_cal()
fit_cal = RTS_cals.rts_cal()
all_xx = [None] * 128
all_yy = [None] * 128

if(len(obs_list) < 1):
    print 'ERROR: No observations found'
    sys.exit(1)

if(options.goodlist):
    good_obs = []
    for line in open(options.goodlist):
        good_obs.append(line.split()[0])

for obs in obs_list:
    
    print obs
    if(len(glob.glob('%s/data/%s/%s/Band*' % (mwa_dir,obs,options.tag))) != 24):
        print 'No Bandpass files %s \n' % obs
        continue
    if(options.goodlist):
        if(obs.split('/')[-1] not in good_obs):
            print 'Bad Observation %s \n' % obs
            continue
    raw_cal = RTS_cals.rts_cal()
    raw_cal.load_all_BP_jones(path='%s/data/%s/%s/' % (mwa_dir,obs,options.tag), raw=True)
    raw_cal.load_all_DI_jones(path='%s/data/%s/%s/' % (mwa_dir,obs,options.tag))
    raw_cal.form_single_jones()

    for i,a in enumerate(raw_cal.antennas):
        if(raw_cal.antennas[i].BP_jones[0] is not None):
#            xx_abs,yy_abs = raw_cal.all_single_jones(antenna=i,cent_chan=cent_chan)
#            xx_abs = (xx_abs * conj(xx_abs))
#            yy_abs = (yy_abs * conj(yy_abs))
            single_jones = raw_cal.all_single_jones(antenna=i,cent_chan=cent_chan)
            xx_abs = single_jones[0][0]
            yy_abs = single_jones[1][1]
            if(all_xx[i] is None):
                all_xx[i] = [array(xx_abs)]
                all_yy[i] = [array(yy_abs)]
            else:
                all_xx[i].append(xx_abs)
                all_yy[i].append(yy_abs)


#    fit_cal = RTS_cals.rts_cal()
#    fit_cal.load_all_BP_jones(path='%s/' % obs, raw=False)
#    fit_cal.load_all_DI_jones(path='%s/' % obs)
#    fit_cal.form_single_jones()

    all_cals.append(raw_cal)

n_bands = 24
n_antennas = 128
xvals = bp_xvals()

for b in range(n_bands):
    for n in range(n_antennas):
        di = [a.antennas[n].DI_jones[b] for a in all_cals if a.antennas[n].DI_jones[b] is not None]
        # the di_ref is the same for all antennas but this way we handle flags
        di_ref = [a.antennas[n].DI_jones_ref[b] for a in all_cals if a.antennas[n].DI_jones_ref[b] is not None]
        # real / imag??
        if(all(d is None for d in di)):
#        if(not any(di)):
            med_cal.antennas[n].DI_jones[b] = None
            avg_cal.antennas[n].DI_jones[b] = None
            med_cal.antennas[n].DI_jones_ref[b] = None
            avg_cal.antennas[n].DI_jones_ref[b] = None
            
        else:
            med_cal.antennas[n].DI_jones[b] = median(di,axis=0)
            avg_cal.antennas[n].DI_jones[b] = mean(di,axis=0)
            med_cal.antennas[n].DI_jones_ref[b] = median(di_ref,axis=0)
            avg_cal.antennas[n].DI_jones_ref[b] = mean(di_ref,axis=0)
            
        bp = [a.antennas[n].BP_jones[b] for a in all_cals if a.antennas[n].BP_jones[b] is not None]
        if(all(ch is None for ch in bp)):
#        if(not any(di)):
            med_cal.antennas[n].BP_jones[b] = None
            avg_cal.antennas[n].BP_jones[b] = None
        else:
            med_cal.antennas[n].BP_jones[b] = median(bp,axis=0)
            avg_cal.antennas[n].BP_jones[b] = median(bp,axis=0)

med_cal.form_single_jones()
avg_cal.form_single_jones()



for i,a in enumerate(med_cal.antennas):
    if(med_cal.antennas[i].BP_jones[0] is not None):
#            xx_abs,yy_abs = med_cal.all_single_jones(antenna=i,cent_chan=cent_chan)
#            xx_abs = (xx_abs * conj(xx_abs))
#            yy_abs = (yy_abs * conj(yy_abs))
            single_jones = raw_cal.all_single_jones(antenna=i,cent_chan=cent_chan)
            xx_abs = single_jones[0][0]
            yy_abs = single_jones[1][1]
        
            if(all_xx[i] is None):
                all_xx[i] = [array(xx_abs)]
                all_yy[i] = [array(yy_abs)]
            else:
                all_xx[i].append(xx_abs)
                all_yy[i].append(yy_abs)

xx_fit = [None] * 128
yy_fit = [None] * 128

for i,a in enumerate(avg_cal.antennas):
    if(avg_cal.antennas[i].BP_jones[0] is not None):
#            xx_abs,yy_abs = avg_cal.all_single_jones(antenna=i,cent_chan=cent_chan)
            single_jones = avg_cal.all_single_jones(antenna=i,cent_chan=cent_chan)
            xx_abs = single_jones[0][0]
            yy_abs = single_jones[1][1]

            single_jones_fit = fit_bandpass_and_reflection(single_jones)

            xx_fit[i] = fit_bandpass_and_reflection(xx_abs)
            yy_fit[i] = fit_bandpass_and_reflection(yy_abs)

            # load single jones values into fit cal BP and load DI terms with 
            # identity

            fit_cal.load_cals_from_single_jones(single_jones_fit,antenna=i,cent_chan=cent_chan)

#            xx_abs = (xx_abs * conj(xx_abs))
#            yy_abs = (yy_abs * conj(yy_abs))
            if(all_xx[i] is None):
                all_xx[i] = [array(xx_abs)]
                all_yy[i] = [array(yy_abs)]
            else:
                all_xx[i].append(xx_abs)
                all_yy[i].append(yy_abs)


RTS_cals.write_BP_files(avg_cal,avg_cal,filename='%s_avg' % options.tag)
RTS_cals.write_DI_files(avg_cal,filename='%s_avg' % options.tag)

# plot antenna 0
for a in all_xx[0]:
    plt.plot(a,color='gray',alpha=0.3)

plt.plot(all_xx[0][-1],color='black')    
    
#plt.title('%s' % c)
plt.savefig('ant0_%s.png' % (tag))
plt.close()
    
    

#for c in cable_types:
#    plt.clf()
#    plt.figure(figsize=(10,2*cable_counts[c]))
#    plot_number = 1
#    for i,b in enumerate(all_xx):
#        if b is not None:
#            if(rts2cables[i] == c):
#                plt.subplot(cable_counts[c],1,plot_number)
#                plt.plot(xvals,(real(all_xx[i][-1]) - (xx_fit[i][0]))[xvals])
#                plt.plot(xvals,(imag(all_xx[i][-1]) - (xx_fit[i][1]))[xvals])
#                plt.text(700,1,'%s %d ' % (rts2tiles[i],i))
#                plot_number += 1
#    plt.title('%s' % c)
#    plt.savefig('%s_%s.png' % (tag,c))
#    plt.close()
