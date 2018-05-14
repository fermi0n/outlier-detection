import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss
import scipy.spatial.distance as ssd
from core import kshape, zscore, _sbd
from operator import itemgetter
import DataAnalysis as da
import CalsLoader
from TileDataAnalysis import DataAnalyserByTile
#Uncomment below line to use DTW
#import mlpy as mlpy
import argparse
import pickle
import datetime


if __name__ == '__main__':

    #Unpickle this data
    tile = 20

    with open ('/lustre/projects/p048_astro/dmendonca/%s_data.pickle'%(tile), 'rb') as f:
        tiledata = pickle.load(f)

    #distancedata = (hardcoded distances)
    #distancedata = tileloader.distances

    colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

    sp = 0
    ppl = 16
    maxv = 12

    fig = plt.figure(figsize=(18.0, 10.0))

    #Converts GPS time to readable time
    timelist = []
    modifiedobslist = []
    for obs in tiledata.obs_list:
        #if int(tiledata.metadata_dict[obs][1]) < 127 and max(tiledata.allx_obs_dict[obs]) > 8 and int(obs) > 1129399688:
        if int(tiledata.metadata_dict[obs][1]) < 127:
            date = datetime.datetime.fromtimestamp(float(obs) + 315964782)
            timelist.append(str(date.day) + '/' + str(date.month) + ':' + str(date.hour) + ':' + str(date.minute) + '.' + str(date.second))
            modifiedobslist.append(obs)
    print("Left out the following observations:")
    print([x for x in tiledata.obs_list if x not in modifiedobslist])

    for channel in range(127):
        xlist = []
        ylist = []
        for obs in modifiedobslist:
            xlist.append(tiledata.allx_obs_dict[obs][channel])
            ylist.append(tiledata.ally_obs_dict[obs][channel])
                #print ("Amp greater than 7 for %s"%obs)

        ax = fig.add_subplot(8,16,sp+1,)
        #ax.plot(modifiedobslist, xlist, color=colours[1], marker='o', markersize=2)
        #ax.plot(modifiedobslist, ylist, color=colours[2], marker='o', markersize=2)
        if (channel < 20):
            ax.plot(modifiedobslist, xlist, color=colours[1], marker='o', markersize=2)
            ax.plot(modifiedobslist, ylist, color=colours[2], marker='o', markersize=2)
        else:
            ax.plot(modifiedobslist, xlist, color=colours[1])
            ax.plot(modifiedobslist, ylist, color=colours[2])

        ax.set_xticklabels(timelist)
        plt.title('Channel %d'%channel, fontsize='8')

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
    plt.suptitle('Amps | Tile %d' %(tile),fontsize=14)

    plt.show()
