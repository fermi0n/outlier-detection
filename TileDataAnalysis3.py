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


if __name__ == '__main__':

    #Get plotting ready
    colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

    sp = 0
    ppl = 16

    maxv = 16.0
    #    maxv = max(max([max(tiledata.allx_obs_dict[obs][ii]) for ii in range(128)]), max([max(tiledata.ally_obs_dict[obs][ii]) for ii in range(128)]))

    fig = plt.figure(figsize=(18.0, 10.0))


    #Unpickle this data
    tile = '21'

    with open ('/lustre/projects/p048_astro/dmendonca/%s_data.pickle'%(tile), 'rb') as f:
        tileloader = pickle.load(f)

    #distancedata = (hardcoded distances)
    distancedata = tileloader.distances

    for obsid, tile, distance in distancedata[:127]:

        ax = fig.add_subplot(8,16,sp+1,)
        ax.plot(tileloader.allx_obs_dict[obsid], colours[1])
        ax.plot(tileloader.ally_obs_dict[obsid], colours[2])
        plt.title('Obs %s\n, distance %.1f'%(obsid, distance), fontsize=6)

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
        plt.suptitle('Amps | Tile %s' %(tile),fontsize=18)

    plt.show()
