import sys
from pylab import *
import argparse
import numpy as np
import matplotlib.pyplot as plt
import TileData as td
import time as time


    def plotObservation(self, tiledata, obs, output='display'):

        colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

        sp = 0
        ppl = 16
        #maxv = 1.5
        maxv = max(max([max(tiledata.allx_obs_dict[obs][ii]) for ii in range(128)]), max([max(tiledata.ally_obs_dict[obs][ii]) for ii in range(128)]))
        print maxv

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

def plotTile(self, distancedata, channel, output='display'):

    colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

    sp = 0
    ppl = 16

    maxv = 1.2
#    maxv = max(max([max(tiledata.allx_obs_dict[obs][ii]) for ii in range(128)]), max([max(tiledata.ally_obs_dict[obs][ii]) for ii in range(128)]))

    fig = plt.figure(figsize=(18.0, 10.0))

    #Converts GPS time to readable time
    timelist = []
    for obs in tiledata.obs_list:
        date = datetime.datetime.fromtimestamp(float(obs) + 315964782)
        timelist.append(str(date.day) + '/' + str(date.month) + ':' + str(date.hour) + ':' + str(date.minute) + '.' + str(date.second))


    for obsid, tile, distance in distancedata:

        #Get the observational data




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

if __name__ == "__main__":

    tileloader = td.TileData()
    plotter = Plotter()

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help="Obtain amplitude data from file", nargs='?')
    parser.add_argument('-p', '--path', help="Base path where the observations are kept", default='/lustre/projects/p048_astro/MWA/data')
    parser.add_argument('-s', '--save', action='store_true', help="Save plot to file rather than display")
    parser.add_argument("data", help="Either channel number or observation ID")
    args = parser.parse_args()
    #Check if its a channel
    if (int(args.data) in range(800)):
        if args.file == None:
            print("Usage: For speed reasons, if plotting a channel, must load data from file")
        else:
            tileloader.load_observations_from_file(args.file)
            #Plot channel
            if (args.save):
                plotter.plotChannel(tileloader, int(args.data), 'save')
            else:
                plotter.plotChannel(tileloader, int(args.data))

    else:  #Is an observation ID
        if args.file == None:  #i.e. loading from path
            tileloader.load_observations(basepath=args.path, obsids=[args.data])
        else:  #Loading from file
	    #    print("Loading observations from file:")
	     #   print(args.file)
            tileloader.load_observations_from_file(args.file)

        if (args.data in tileloader.obs_list):  #All is good - the observation ID was in the file or was found in the filesystem
            #Plot observation
            if (args.save):
                plotter.plotObservation(tileloader, args.data, 'save')
            else:
                plotter.plotObservation(tileloader, args.data)
        else:
            if (args.file==None):
                print("Usage: You did not provide a valid observation ID that we could find in path selected:")
                print(args.path)
            else:
                print("Usage: You must provide an observation ID which is present in the file you have selected")
                print(args.file)
