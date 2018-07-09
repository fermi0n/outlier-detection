import RTS_cals
from astropy.io import fits
from optparse import OptionParser,OptionGroup
from matplotlib import pyplot as plt

from numpy import conj, sqrt,zeros, fft, array,std, mean, polyfit, poly1d, arange, real, imag, median, shape, ravel, where

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

parser.add_option('--split', action="store_true")

parser.add_option('--metadata', action="store_true", help="If flag set, store the metadata in a separate file")

parser.add_option('--nogains', action="store_true", help="If flag set, don't store gains")

(options, args) = parser.parse_args()

obs_list = open(args[0]).readlines()
obs_list = [[l.split()[0], l.split()[1]] for l in obs_list]

cal = 0

if(options.mwa_dir is None):
    mwa_dir = os.getenv('MWA_DIR','/astro/mwaeor/MWA/')
else:
    mwa_dir = options.mwa_dir

if options.subdir is '':
    subdir = 'example'
else:
    subdir = options.subdir

all_cals = []
all_xx = [None] * 128
all_yy = [None] * 128

tile = options.tile

for val in obs_list:
    if options.split:
        subdir = val[1].split('/')[-1]
    else:
        subdir = options.subdir

    obs = val[0]
    files = glob.glob('%s/data/%s/%s/*.fits' %(mwa_dir, obs, subdir))
    meta_file = fits.open(files[0])

    cent_chan = meta_file[0].header['CENTCHAN']
    gridnum = meta_file[0].header['GRIDNUM']
    rts_inputs = meta_file[1].data['Input']
    cables = meta_file[1].data['Flavors']
    tilenames = meta_file[1].data['Tile']
    dipole_delays = meta_file[1].data['Delays']
    cable_types = set(cables)

    cable_counts = dict(Counter(cables[::2]))
    rts_ants = rts_inputs / 2
    rts2cables = dict(zip(rts_ants[::2], cables[::2]))
    rts2tiles = dict(zip(rts_ants[::2], tilenames[::2]))

    raw_cal = RTS_cals.rts_cal()

    ##Order of delays is XX1,YY1,XX2,YY2 etc so select just XX   (JUST DO XX FOR NOW, TODO: ADD YY)
    XX = arange(0,256,2).astype(int)

    ##If a delay is set to 32, the dipole is flagged
    tile_inds,dipole_flags = where(dipole_delays[XX,:] == 32)

    tile_flags = tilenames[tile_inds]

    ##Make flags direction
    flags_dict = {}
    ##Set up lists to contain flags for tiles
    #for tile in set(tile_flags):
        #flags_dict['Tile%d' %int(tile)] = []

    ##Populate flags. Allows for more than one flag per tile
    for tilenum,dipole in zip(tile_flags,dipole_flags):
        flags_dict['Tile%d' %int(tilenum)].append(int(dipole))

    print(flags_dict)
    #If flagged to save the metadata
    if (options.metadata):
        f = open("metadata-" + str(tile) + ".txt", "a")
        f.write("%s,%s,%s\n" %(obs, cent_chan, gridnum))
	if 'Tile%s'%tile in flags_dict:
            f.write("Tile %d," %tile)
       	    for dipole in flags_dict.get('Tile%s'%tile):
                f.write("%d," %dipole)
            f.write("\n")
        f.flush()
        f.close()
    #If not flagged to skip the gains
    if (options.nogains):
        continue
    # Actually reading of calibration files. The BP files contain both the 'raw'
    # BP calibrations and the functional fits across each coarse channel. In
    # this case we read the raw values which are probably more interesting
    # for QA

    #Adding some error handling. If this fails, then just skip this Observation

    try:
        raw_cal.load_all_BP_jones(path='%s/data/%s/%s/' % (mwa_dir,obs,subdir), raw=True)
        raw_cal.load_all_DI_jones(path='%s/data/%s/%s/' % (mwa_dir,obs,subdir))
        # Forms a product of the BP and DI Jones terms
        raw_cal.form_single_jones()
    except IOError as e:
        print("Couldn't load this one")
        print e
        continue
    except:
        continue


    # The single jones matrices are still carried in the normal RTS per
    # coarse channel representation. For ease of inspection, we call a
    # helper function which loads all of the gains into a single list per
    # antenna

    for i,a in enumerate(raw_cal.antennas):
#    i = options.tile
        #Do things in a hacked way just to try to not modify Bart's script as much as possible
        if (i != tile):
            continue

        if raw_cal.antennas[i].BP_jones[0] is not None:
            single_jones = raw_cal.all_single_jones(antenna=i,cent_chan=cent_chan)
        # Futher we separately gather just the XX and YY values into lists
            xx_abs = single_jones[0][0]
            yy_abs = single_jones[1][1]

            xx_modified_vals = [sqrt(a.real**2 + a.imag**2) for a in xx_abs if a != 0.0]
            yy_modified_vals = [sqrt(a.real**2 + a.imag**2) for a in yy_abs if a != 0.0]

            filename = "observations-" + str(tile) + ".txt"

            f = open(filename, "a")
            f.write("%s,%s,X," %(obs, tile))
            for value in xx_modified_vals:
                f.write("%f," %value)
            f.write("\n")
            f.write("%s,%s,Y," %(obs, tile))
            for value in yy_modified_vals:
                f.write("%f," %value)
            f.write("\n")
            f.flush()
            f.close()
