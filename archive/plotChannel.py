import glob
import sys
from pylab import *
matplotlib.use('Agg')
from astropy.io import fits
from fit_coarse_channels import fit_bandpass, fit_sigma
import numpy as np
import matplotlib.pyplot as plt
import math

import CalsLoader

class Observations:
    def __init__(self):
        self.obs_list = None   #not sure if I can use sys.argv here
        self.all_xx = [None] * 128
        self.all_yy = [None] * 128

        self.allx_obs_dict = {}
        self.ally_obs_dict = {}
        self.missing_list = []
        self.flaggedlist = []
        self.steepnesses = []

        self.scaled_x_dict = {}
        self.scaled_y_dict = {}

        self.minval = 0.0
        self.maxval = 100.0

        self.discretisedAmps = {}

    def load_observations(self, basepath = './', obsids = None):

        for obs in obsids:
            cals_loader = CalsLoader.CalsLoader()
            cals_loader.obtainAmplitudeforObservation(basepath + '/'+ obs)

            self.allx_obs_dict[obs] = cals_loader.JPX
            self.ally_obs_dict[obs] = cals_loader.JQY

            self.discretisedAmps[obs] = [None]*128

        self.obs_list = obsids

    def saveObservations(self, filename):
        #Print all the observations to csv files in a simple format:
        #ObsId, antenna, channel, x-amp, Y-amp
        f= open(filename,"a")

        for obs in self.obs_list:
            for antenna in range(128):
                if self.allx_obs_dict[obs][antenna] is not None:
                    f.write("%s,%s,X, " %(obs, antenna))
                    for channel, value in enumerate(self.allx_obs_dict[obs][antenna]):
                        f.write("%f, "%self.allx_obs_dict[obs][antenna][channel])
                    f.write("\n")
                    f.write("%s,%s,Y, " %(obs, antenna))
                    for channel, value in enumerate(self.ally_obs_dict[obs][antenna]):
                        f.write("%f, "%self.ally_obs_dict[obs][antenna][channel])
                    f.write("\n")

    def plotChannel(self, channel):

        colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

        sp = 0
        ppl = 16
        maxv = 1.5

        fig = plt.figure(figsize=(18.0, 10.0))

        #Create the xtick labels
        timelist = []
        for obs in self.obs_list:
            date = datetime.datetime.fromtimestamp(float(obs) + 315964782)
            timelist.append(str(date.day) + '/' + str(date.month) + ':' + str(date.hour) + ':' + str(date.minute) + '.' + str(date.second))

        for ii in range(128):
            xlist = []
            ylist = []
            for obs in self.obs_list:
                xlist.append(self.allx_obs_dict[obs][ii][channel])
                ylist.append(self.ally_obs_dict[obs][ii][channel])

            ax = fig.add_subplot(8,16,sp+1,)
            ax.plot(self.obs_list, xlist, color=colours[1])
            ax.plot(self.obs_list, ylist, color=colours[2])

            ax.set_xticklabels(timelist)
            #ax.tick_params('x', labelrotation=90)
            plt.title('Tile %d'%ii)

            if ppl != 16:
                plt.setp(ax.get_xticklabels(), visible=False) # plot setup
                plt.setp(ax.get_yticklabels(), visible=False)

            if ppl == 16:
                ppl = 0
                plt.setp(ax.get_xticklabels(), visible=False)
            ppl += 1

            plt.ylim([-0.1,maxv])

            if sp == 15:
                XX_amp, = ax.plot([],[], colours[1],label='XX',linewidth=3.0)
                YY_amp, = ax.plot([],[], colours[2],label='YY',linewidth=3.0)
                ax.legend((XX_amp,YY_amp),('XX','YY'), bbox_to_anchor=(0, 2, .12, .12),prop={'size':14})

            sp += 1

        plt.tight_layout()
        fig.subplots_adjust(top=0.9)
        plt.suptitle('Amps | Channel %d' %(channel),fontsize=14)
        print(timelist)
        plt.show()


    def plotObservation(self, obs):

        colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

        sp = 0
        ppl = 16
        maxv = 1.5

        fig = plt.figure(figsize=(18.0, 10.0))

        for ii in range(128):

#ax.plot(freq[freq_idx], self.JPX[ii], colours[1])
#ax.plot(freq[freq_idx], self.JQY[ii], colours[2])
#ax.set_xlim(min(freq[freq_idx]), max(freq[freq_idx]))

            ax = fig.add_subplot(8,16,sp+1,)
            ax.plot(self.allx_obs_dict[obs][ii], colours[1])
            ax.plot(self.ally_obs_dict[obs][ii], colours[2])
            plt.title('Tile %d'%ii)

            if ppl != 16:
                plt.setp(ax.get_xticklabels(), visible=False) # plot setup
                plt.setp(ax.get_yticklabels(), visible=False)

            if ppl == 16:
                ppl = 0
                plt.setp(ax.get_xticklabels(), visible=False)
            ppl += 1

            plt.ylim([-0.1,maxv])

            if sp == 15:
                XX_amp, = ax.plot([],[], colours[1],label='XX',linewidth=3.0)
                YY_amp, = ax.plot([],[], colours[2],label='YY',linewidth=3.0)
                ax.legend((XX_amp,YY_amp),('XX','YY'), bbox_to_anchor=(0, 2, .12, .12),prop={'size':14})

            sp += 1

        plt.tight_layout()
        fig.subplots_adjust(top=0.9)
        plt.suptitle('Amps | Observation %s' %(obs),fontsize=18)

        if sys.argv[2] == 'save':
            savefig('Amps_%s.png'%obs, bbox_inches='tight')
        else:
            plt.show()

tileloader = Observations()
tileloader.load_observations(basepath='/lustre/projects/p048_astro/MWA/data', obsids = [sys.argv[1]])
tileloader.saveObservations("observations.txt")
#tileloader.plotObservation(obs=sys.argv[1])
