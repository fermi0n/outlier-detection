import glob
from pylab import *
matplotlib.use('Agg')
from astropy.io import fits
from fit_coarse_channels import fit_bandpass, fit_sigma

import RTS_cals

class tile_loader:
    def __init__(self):
        self.obs_list = None   #not sure if I can use sys.argv here
        self.all_xx = [None] * 128
        self.all_yy = [None] * 128

        self.allx_obs_dict = {}
        self.ally_obs_dict = {}
        self.missing_list = []

    def load_observations(self, basepath = './', obsids = None):

        if obsids == None:
            self.obs_list = glob.glob('%s/1*' % basepath)
        else:
            self.obs_list = ['%s/%s' % (basepath, obsid) for obsid in obsids]

        print self.obs_list

        for obs in self.obs_list:

            metafile = (glob.glob('%s/*_metafits*' % obs))[0]
            meta_file = fits.open(metafile)
            flags = meta_file[1].data['Flag']
            flags = flags[::2]
            tiles = meta_file[1].data['Tilename']
            tiles = tiles[::2]
            di_files = glob.glob('%s/DI*.dat' % obs)
            bp_files = glob.glob('%s/Band*.dat' % obs)
            if(len(di_files) != 24 or len(bp_files) !=24):
                print 'Missing cal files: ', obs, len(di_files), len(bp_files)
                self.missing_list.append(obs)
                continue

            rts_cal = RTS_cals.rts_cal()
            rts_cal.load_all_DI_jones(path='%s/' % obs)
            rts_cal.load_all_BP_jones(path='%s/' % obs, raw=True,add_flagged=False)
            rts_cal.load_metadata(metafile)
            rts_cal.apply_metafits_gains()
            rts_cal.form_single_jones()
            low_gain = False
            high_gain = False
            high_tiles = []
            low_tiles = []

            for i,a in enumerate(rts_cal.antennas):
                if(rts_cal.antennas[i].BP_jones[0] is not None):
                    xx_abs,yy_abs = rts_cal.all_single_jones(antenna=i,reverse_bands=True,correct_gain_jump=False,conjugate_jones=False)
                    xx_abs = (xx_abs * conj(xx_abs))
                    yy_abs = (yy_abs * conj(yy_abs))
                    if(self.all_xx[i] == None):
                        self.all_xx[i] = [xx_abs]
                        self.all_yy[i] = [yy_abs]
                    else:
                        self.all_xx[i].append(xx_abs)
                        self.all_yy[i].append(yy_abs)

            self.allx_obs_dict[obs[-10:]] = self.all_xx
            self.ally_obs_dict[obs[-10:]] = self.all_yy

        print self.allx_obs_dict

tileloader = tile_loader()
tileloader.load_observations(basepath='/lustre/projects/p048_astro/MWA/data', obsids = ['1061306296', ])

#obs_list = glob.glob('%s/1*' % sys.argv[1])
#obs_list = ['%s/1061311664'% sys.argv[1], '%s/1061311784'% sys.argv[1], '%s/1061312032'% sys.argv[1], '%s/1061312152'% sys.argv[1]]
# obs_list = ['%s/1061311664'% sys.argv[1]]
# all_xx = [None] * 128
# all_yy = [None] * 128
# all_std_xx = [None] * 128
# all_std_yy = [None] * 128
#
# allx_obs_dict = {}
# ally_obs_dict = {}
#
# missing_list = []
# low_gain_list = []
# high_gain_list = []
# good_list = []
#
# for obs in obs_list:
#
#     metafile = (glob.glob('%s/*_metafits*' % obs))[0]
#     meta_file = fits.open(metafile)
#     flags = meta_file[1].data['Flag']
#     flags = flags[::2]
#     tiles = meta_file[1].data['Tilename']
#     tiles = tiles[::2]
#     di_files = glob.glob('%s/DI*.dat' % obs)
#     bp_files = glob.glob('%s/Band*.dat' % obs)
#     if(len(di_files) != 24 or len(bp_files) !=24):
#         print 'Missing cal files: ', obs, len(di_files), len(bp_files)
#         missing_list.append(obs)
#         continue
#
#     rts_cal = RTS_cals.rts_cal()
#     rts_cal.load_all_DI_jones(path='%s/' % obs)
#     rts_cal.load_all_BP_jones(path='%s/' % obs, raw=True,add_flagged=False)
#     rts_cal.load_metadata(metafile)
#     rts_cal.apply_metafits_gains()
#     rts_cal.form_single_jones()
#     low_gain = False
#     high_gain = False
#     high_tiles = []
#     low_tiles = []
#     for i,a in enumerate(rts_cal.antennas):
#         if(rts_cal.antennas[i].BP_jones[0] is not None):
#             xx_abs,yy_abs = rts_cal.all_single_jones(antenna=i,reverse_bands=True,correct_gain_jump=False,conjugate_jones=False)
#             xx_abs = (xx_abs * conj(xx_abs))
#             yy_abs = (yy_abs * conj(yy_abs))
#             if(all_xx[i] == None):
#                 all_xx[i] = [xx_abs]
#                 all_yy[i] = [yy_abs]
#             else:
#                 all_xx[i].append(xx_abs)
#                 all_yy[i].append(yy_abs)
#             if(flags[i] == 0):
#                 mean_xx = mean(xx_abs)
#                 mean_yy = mean(yy_abs)
#
#                 # I have read in the bandpass amplitudes as xx_abs and yy_abs
#                 # Now I apply two simple test to see if any of the tiles
#                 # have either very high or very low amplitudes
#
#                 # Fit 2nd order polynomial to find 'unsmooth' BPs
#
#                 model_xx = fit_bandpass(xx_abs)
#                 resid_xx = xx_abs - model_xx
#                 norm_xx = xx_abs / model_xx
#                 std_xx = std(norm_xx)
#                 model_yy = fit_bandpass(yy_abs)
#                 resid_yy = yy_abs - model_yy
#                 norm_yy = yy_abs / model_yy
#                 std_yy = std(norm_yy)
#
#                 # This variable is for collecting data over all observations if I want
#                 # to later plot distribution in interactively in ipython
#                 if(all_std_xx[i] == None):
#                     all_std_xx[i] = [max(abs(1.0 - norm_xx))]
#                     all_std_yy[i] = [max(abs(1.0 - norm_yy))]
#
#                 else:
#                     all_std_xx[i].append(max(abs(1.0 - norm_xx)))
#                     all_std_yy[i].append(max(abs(1.0 - norm_yy)))
#
#                 # Standard deviations are taken from value over all obs to
#                 # avoid self-contamination of bad obsids
#                 std_xx = 0.1
#                 std_yy = 0.1
#
#                 # if(mean_xx < 0.1 or mean_yy < 0.1):
#                 #     print 'Low Gain', i,mean_xx,mean_yy,tiles[i]
#                 #     low_gain = True
#                 #     low_tiles.append(i)
#                 #     print i, max(xx_abs), max(yy_abs)
#                 # else:
#                 #     if(max(abs(1.0 - norm_xx)) > 5.0 * std_xx or max((1.0 - norm_yy)) > 5.0 * std_yy):
#                 #         print 'High Gain', i,max(norm_xx),std_xx,max(norm_yy),std_yy,tiles[i]
#                 #         clf()
#                 #         plot(model_xx, label='Model xx')
#                 #         plot(xx_abs, label='XX abs')
#                 #         plot(norm_xx, label='Norm xx')
#                 #         plot(model_yy, label='Model y')
#                 #         plot(yy_abs, label='Y Abs')
#                 #         plot(norm_yy, label='Norm y')
# 			    #         legend(loc="upper center")
#                 #         title('High Gain %s %s' % (obs,tiles[i]))
#                 #         savefig('%s_high_gain.png' % (obs.split('/')[-1]))
#                 #         high_gain = True
#                 #         high_tiles.append(i)
#
#         # This is to store all the observations
#     allx_obs_dict[obs[-10:]] = all_xx
#     ally_obs_dict[obs[-10:]] = all_yy
#     print 'printing all_x'
#     print all_xx
#     print(type(all_xx))
#     print "Printing all_xx i (antenna) = 3, channel = 4"
#     print all_xx[3][0][4]
#     # if(low_gain):
#     #     low_gain_list.append((obs,low_tiles,high_tiles))
#     # else:
#     #     if(high_gain):
#     #         high_gain_list.append((obs,low_tiles,high_tiles))
#     #     else:
#     #         good_list.append(obs)
# #
# # tiles2rts = dict(zip(tiles,range(len(tiles))))
# #
# # low_counts = dict()
# # for l in low_gain_list:
# #     if(len(l[1]) < 3):
# #         for t in l[1]:
# #             if(t in low_counts):
# #                 low_counts[t] += 1
# #             else:
# #                 low_counts[t] = 1
#
# #for key in sorted(mydict.iterkeys()):
#
# obs_listing = [key for key in sorted(allx_obs_dict.iterkeys())]
# example = [value[3][0][4] for (key, value) in sorted(allx_obs_dict.items())]  #antenna = 3, channel = 4
#
# #print len(example[])
#
# print "printing obs example"
# print obs_listing
# print example
# obs_listing.append('1061311784')
# example.append(1.0242324632995+0j)
# #Graph it for antenna = 3, channel = 4
# clf()
# plot(obs_listing, example, label='xx_abs')
# legend(loc="upper center")
# title('Antenna = 3, channel = 4')
# show()
#
# #print len(example[0])
#
# #In [227]: for a in all_std_xx:
# #    if(a is not None):
# #        plot(a)
