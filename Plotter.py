#import glob
import sys
from pylab import *
import argparse
#matplotlib.use('Agg')
#from astropy.io import fits
#from fit_coarse_channels import fit_bandpass, fit_sigma
import numpy as np
import matplotlib.pyplot as plt
#import math
import TileData as td
import time as time

class Plotter:
    def __init__(self):
        random = 0

    def plotObservation(self, tiledata, obs, output='display'):

        colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

        sp = 0
        ppl = 16
        maxv = 1.5

        fig = plt.figure(figsize=(18.0, 10.0))

        for ii in range(128):

            ax = fig.add_subplot(8,16,sp+1,)
            ax.plot(tiledata.allx_obs_dict[obs][ii], colours[1])
            ax.plot(tiledata.ally_obs_dict[obs][ii], colours[2])
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

        if output == 'save':
            savefig('%sAmps_Obs_%s.png'%(time.strftime("%d%b_%H:%M:%S", time.localtime()), obs), bbox_inches='tight')
        else:
            plt.show()

    def plotChannel(self, tiledata, channel, output='display'):

        colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

        sp = 0
        ppl = 16
        maxv = 1.5

        fig = plt.figure(figsize=(18.0, 10.0))

        #Converts GPS time to readable time
        timelist = []
        for obs in tiledata.obs_list:
            date = datetime.datetime.fromtimestamp(float(obs) + 315964782)
            timelist.append(str(date.day) + '/' + str(date.month) + ':' + str(date.hour) + ':' + str(date.minute) + '.' + str(date.second))

        for ii in range(128):
            xlist = []
            ylist = []
            for obs in tiledata.obs_list:
                xlist.append(tiledata.allx_obs_dict[obs][ii][channel])
                ylist.append(tiledata.ally_obs_dict[obs][ii][channel])

            ax = fig.add_subplot(8,16,sp+1,)
            ax.plot(tiledata.obs_list, xlist, color=colours[1])
            ax.plot(tiledata.obs_list, ylist, color=colours[2])

            ax.set_xticklabels(timelist)
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

        if output == 'save':
            savefig('%s_Amps_Channel_%s.png'%(time.strftime("%d%b_%H:%M:%S", time.localtime()),channel), bbox_inches='tight')
        else:
            plt.show()

if __name__ == "__main__":

    tileloader = td.TileData()
    plotter = Plotter()

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help="Obtain amplitude data from file", nargs='?')
    parser.add_argument('-s', '--save', action='store_true', help="Save plot to file rather than display")
    parser.add_argument("data", help="Either channel number or observation ID")
    args = parser.parse_args()
    #Check if its a channel
    if (int(args.data) in range(800)):
        if args.file == None:
            print("Usage: For speed reasons, want to load this data from file")
        else:
            tileloader.loadObservationsFromFile(args.file)
            #Plot channel
            #First check that it is in fact a channel
            if (args.save):
                plotter.plotChannel(tileloader, int(args.data), 'save')
            else:
                plotter.plotChannel(tileloader, int(args.data))

    else:
        if args.file == None:
            tileloader.load_observations(basepath='/lustre/projects/p048_astro/MWA/data/', obsids=[args.data])
        else:
	    print("Loading observations from file:")
	    print(args.file)
            tileloader.loadObservationsFromFile(args.file)
        if (args.data in tileloader.obs_list):
            #Plot observation
            if (args.save):
                plotter.plotObservation(tileloader, args.data, 'save')
            else:
                plotter.plotObservation(tileloader, args.data)
        else:
            print("Usage: You must provide an observation ID or a channel number")


    #tileloader.loadObservationsFromFile('observations.txt')
    #tileloader.load_observations(basepath='/lustre/projects/p048_astro/MWA/data', obsids = [sys.argv[1]])
    #tileloader.saveObservations("observations.txt")
