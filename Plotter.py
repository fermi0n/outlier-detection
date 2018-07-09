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

    def plot_observation_tile(self, tiledata, obs, tile, output='display'):

        colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

        maxv = max (max(tiledata.allx_obs_dict[obs][tile], tiledata.ally_obs_dict[obs][tile]))

        plt.plot(tiledata.allx_obs_dict[obs][tile], colours[1])
        plt.plot(tiledata.ally_obs_dict[obs][tile], colours[2])
        plt.title('Tile %d' %tile)

        if (output == 'save'):
            savefig('Amps_Obs_%s_Tile_%s.png' %(obs, tile), bbox_inches='tight')
        else:
            plt.show()

    def plot_channel_tile(self, tiledata, channel, tile, output='display', metafits=None):

        colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

        timelist = []
        sortedlist = sorted([int(obs) for obs in tiledata.obs_list])
        for obs in sortedlist:
            #date = datetime.datetime.fromtimestamp(float(obs) + 315964782)
            #timelist.append(str(date.day) + '/' + str(date.month) + ':' + str(date.hour) + ':' + str(date.minute) + '.' + str(date.second))
            timelist.append((obs - sortedlist[0]) / 248.0)

        xlist = []
        ylist = []
        for obs in tiledata.obs_list:
            xlist.append(tiledata.allx_obs_dict[str(obs)][tile][channel])
            ylist.append(tiledata.ally_obs_dict[str(obs)][tile][channel])

        fig, ax1 = plt.subplots()
        ax1.plot(timelist, xlist, color='green', marker='o')
        ax1.plot(timelist, ylist, color='red', marker='o')
        ax1.set_title('Tile %d Channel %d' %(tile, channel))
        ax1.set_xlabel("Observation")
        ax1.set_ylabel("Amps")

        if (metafits is not None):
            pointings=[tiledata.metadata_dict[x][1] for x in sortedlist]
            centchannels = [tiledata.metadata_dict[x][2] / 121.0 for x in sortedlist]

            ax2 = ax1.twinx()
            ax2.plot(timelist, pointings, color='yellow', marker='o')
            ax2.plot(timelist, centchannels, color='blue', marker='o')
        if (output == 'save'):
            savefig('Amps_Channel_%s_Tile_%s.png' %(channel, tile), bbox_inches='tight')
        else:
            plt.show()

    def plotObservation(self, tiledata, obs, output='display'):

        colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

        sp = 0
        ppl = 16
        #maxv = 1000.0
        maxv = max(max([max(tiledata.allx_obs_dict[obs][ii]) for ii in range(128)]), max([max(tiledata.ally_obs_dict[obs][ii]) for ii in range(128)]))
        print (maxv)

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
        plt.suptitle('Amps | Channel %d' %channel, fontsize=14)
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
    parser.add_argument('-m', '--metafits', help="Obtain metafits information from file", nargs='?')
    parser.add_argument('-p', '--path', help="Base path where the observations are kept", default='/lustre/projects/p048_astro/MWA/data')
    parser.add_argument('-s', '--save', action='store_true', help="Save plot to file rather than display")
    parser.add_argument('-t', '--tile', type=int, help="Tile number to plot")
    parser.add_argument("data", help="Either channel number or observation ID")
    args = parser.parse_args()
    #Check if its a channel
    if args.save:
        param = 'save'
    else:
        param = 'display'
    if int(args.data) in range(800):
        if args.file is None:
            print("Usage: For speed reasons, if plotting a channel, must load data from file")
        else:
            tileloader.load_observations_from_file(args.file)
            if (args.metafits is not None):
                print("Loading metafits")
                tileloader.load_metafits_from_file(args.metafits)
            #Plot channel
            if args.tile is not None:
                plotter.plot_channel_tile(tileloader, int(args.data), args.tile, param, metafits=args.metafits)
            else:
                plotter.plotChannel(tileloader, int(args.data), param)

    else:  # Is an observation ID
        if args.file is None:  #i.e. loading from path
            tileloader.load_observations(basepath=args.path, obsids=[args.data])
        else:  #Loading from file
	    #    print("Loading observations from file:")
	     #   print(args.file)
            tileloader.load_observations_from_file(args.file)

        if args.data in tileloader.obs_list:  #All is good - the observation ID was in the file or was found in the filesystem
            #Plot observation
            # Check whether is a single tile or all
            if args.tile is not None:
                plotter.plot_observation_tile(tileloader, args.data, args.tile, param)
            else:
                plotter.plotObservation(tileloader, args.data, param)
        else:
            if args.file is None:
                print("Usage: You did not provide a valid observation ID that we could find in path selected:")
                print(args.path)
            else:
                print("Usage: You must provide an observation ID which is present in the file you have selected")
                print(args.file)
