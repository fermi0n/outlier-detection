import matplotlib
#matplotlib.use('Agg')
#Uncomment below line to use DTW
#import mlpy as mlpy
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss
import scipy.spatial.distance as ssd
from core import kshape, zscore, _sbd
from operator import itemgetter
import DataAnalysis as da
#from TileDataAnalysis import DataAnalyserByTile
#import CalsLoader
import argparse
import pickle
import itertools
import time
import TileDataAnalysiswithClustering

if __name__ == '__main__':
#    dm = da.DistanceMeasures()
    tileloader = TileDataAnalysiswithClustering.DataAnalyserByTile()

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--numclusters', help="Number of clusters to plot", type=int, default=2)
    parser.add_argument("files", nargs='*', help="Files to load - assume that the filename contains the tile number")
    parser.add_argument('-m', '--masterlist', help='where the master list of observation IDs is')

    args = parser.parse_args()

    #Going to represent this as a dict of dicts. first dict key is the tilenumber, second is the observation ID
    Matrix_dict = {}
    obs_list = []

    #Load observations IDs from master list.
    # f = open (args.masterlist, "r")
    # for line in f:
    #     values = line.split(',')
    #     if (values[0] not in obs_list):
    #         self.obs_list.append(values[0])
    # f.close()

    for f in args.files:
        print ("Doing file %s"%f)
        tilenumber = f[-6:-4]  #Need to update this to handle > 2 digits
        tileloader = TileDataAnalysiswithClustering.DataAnalyserByTile()
        tileloader.load_observationsfromFile(f)

        listofsilhouettes, listofclusters = tileloader.TileKShapeClustering(tilenumber)
        Matrix_dict[tilenumber] = {}

        for key, cluster in enumerate(listofclusters[args.numclusters]):
            print cluster
            for observation in cluster:
                Matrix_dict[tilenumber][observation] = key

    print Matrix_dict

    #Save this dict so that I don't need to go through this again.
    pickle.dump(Matrix_dict, open( "matrix.p", "wb" ) )

    data = np.array(Matrix_dict)

    fig, ax = plt.subplots(1, 1)
    cax = ax.imshow(data, cmap='viridis', interpolation='nearest', origin='lower')
    ax.set_title('Heat map of observation clusters')
    ax.set_xlabel("tile")
    ax.set_ylabel("observations")
    plt.show()
