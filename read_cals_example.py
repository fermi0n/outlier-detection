import RTS_cals
from astropy.io import fits
from optparse import OptionParser,OptionGroup
from matplotlib import pyplot as plt

from numpy import conj, sqrt,zeros, fft, array,std, mean, polyfit, poly1d, arange, real, imag, median, shape, ravel

import glob
from collections import Counter
import os,sys
from fit_cable_reflections import fit_bandpass_and_reflection, bp_xvals

usage = 'Usage: read_cals_example.py [text file of obsIDs] [metafits]'

parser = OptionParser(usage=usage)

parser.add_option('--goodlist',dest='goodlist',type='string',default='',help='List of good obsids to be processed')

parser.add_option('--subdir',dest='subdir',type='string',default='',help='Label to subdir output plots')

parser.add_option('--mwa_dir',dest='mwa_dir',type='string',default=None,help='Base path to data directory (if different from $MWA_DIR)')

parser.add_option('--tile', dest='tile', type='int',default=5,help='Tile number (antenna number) to analyse')

(options, args) = parser.parse_args()

obs_list = open(args[0]).readlines()
obs_list = [l.split()[0] for l in obs_list]
meta_file = fits.open(args[1])
cal = 0

if(options.mwa_dir is None):
    mwa_dir = os.getenv('MWA_DIR','/astro/mwaeor/MWA/')
else:
    mwa_dir = options.mwa_dir

# the script currently uses a representative single metafits file to read things like the tilenames and cable types

cent_chan = meta_file[0].header['CENTCHAN']
rts_inputs = meta_file[1].data['Input']
cables = meta_file[1].data['Flavors']
tilenames = meta_file[1].data['Tile']
cable_types = set(cables)

cable_counts = dict(Counter(cables[::2]))
rts_ants = rts_inputs / 2
rts2cables = dict(zip(rts_ants[::2],cables[::2]))
rts2tiles = dict(zip(rts_ants[::2],tilenames[::2]))

if(options.subdir is ''):
    subdir = 'example'
else:
    subdir = options.subdir

all_cals = []
all_xx = [None] * 128
all_yy = [None] * 128

tile = options.tile

for obs in obs_list:
    raw_cal = RTS_cals.rts_cal()
    # Actually reading of calibration files. The BP files contain both the 'raw'
    # BP calibrations and the functional fits across each coarse channel. In
    # this case we read the raw values which are probably more interesting
    # for QA

    raw_cal.load_all_BP_jones(path='%s/data/%s/%s/' % (mwa_dir,obs,options.subdir), raw=True)
    raw_cal.load_all_DI_jones(path='%s/data/%s/%s/' % (mwa_dir,obs,options.subdir))
    # Forms a product of the BP and DI Jones terms
    raw_cal.form_single_jones()

    # The single jones matrices are still carried in the normal RTS per
    # coarse channel representation. For ease of inspection, we call a
    # helper function which loads all of the gains into a single list per
    # antenna

    #bp_nonzeros = bp_xvals()
    #if options.tile not in raw_cal.antennas:
#        print("FATAL ERROR: Tile not found")
#        print(raw_cal.antennas)
#        print(options.tile)
    for i,a in enumerate(raw_cal.antennas):
#    i = options.tile

        if(raw_cal.antennas[i].BP_jones[0] is not None):
            single_jones = raw_cal.all_single_jones(antenna=i,cent_chan=cent_chan)
        # Futher we separately gather just the XX and YY values into lists
            xx_abs = single_jones[0][0]
            yy_abs = single_jones[1][1]

        # if (all_xx is None):
        #     all_xx = [[sqrt(a.real**2 + a.imag**2) for a in xx_abs if a != 0.0]]
        #     all_yy = [[sqrt(a.real**2 + a.imag**2) for a in yy_abs if a != 0.0]]
        #         #all_yy[i] = [[sqrt(yy_abs[i].real**2 + yy_abs[i].imag**2) for i in bp_nonzeros]]
        # else:
        #     all_xx.append([sqrt(a.real**2 + a.imag**2) for a in xx_abs if a != 0.0])
        #     all_yy.append([sqrt(a.real**2 + a.imag**2) for a in yy_abs if a != 0.0])

            if (all_xx[i] is None):
                all_xx[i] = [[sqrt(a.real**2 + a.imag**2) for a in xx_abs if a != 0.0]]
                all_yy[i] = [[sqrt(a.real**2 + a.imag**2) for a in yy_abs if a != 0.0]]
            else:
                all_xx[i].append([sqrt(a.real**2 + a.imag**2) for a in xx_abs if a != 0.0])
                all_yy[i].append([sqrt(a.real**2 + a.imag**2) for a in yy_abs if a != 0.0])
#
#
    # in case we want them later, hold all of the calibrations in a single list
    # as well

    all_cals.append(raw_cal)

    filename = "observations-" + str(tile) + ".txt"

    f = open(filename, "a")
    f.write("%s,%s,X," %(obs, tile))
    for value in all_xx[tile][cal]:
        f.write("%f," %value)
    f.write("\n")
    f.write("%s,%s,Y," %(obs, tile))
    for value in all_yy[tile][cal]:
        f.write("%f," %value)
    f.write("\n")
    f.flush()
    f.close()
    cal += 1
    #except Exception, err:
    #    print("Error writing observation %s, Error is: %s" %(obs,err))

# all_xx is now a list of lists; [n_antennas] with each antenna containing [n_cals]. Maybe there is a better format for this?
print(all_xx[0][0])
# plot abs(xx) and abs(yy) for antenna 0 for cal 0 and 1
plt.plot(abs(ravel(all_xx[0][0])),color='red',alpha=0.3)
plt.plot(abs(ravel(all_yy[0][0])),color='blue',alpha=0.3)
plt.plot(abs(ravel(all_xx[0][1])),color='green',alpha=0.3)
plt.plot(abs(ravel(all_yy[0][1])),color='yellow',alpha=0.3)

plt.show()

#plt.savefig('example.png')
