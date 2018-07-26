import matplotlib
#matplotlib.use('Agg')
#Uncomment below line to use DTW
#import mlpy as mlpy
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
import TileDataAnalysiswithClustering
#from TileDataAnalysis import DataAnalyserByTile
#import CalsLoader
import argparse
import pickle
import itertools
import time
import pandas as pd
import glob
from collections import Counter


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--epsilon', help="Value of epsilon, the density range", default=0.005, type=float)
    parser.add_argument('-d', '--dir', help="Load amplitude data from this directory. Assume filenames are observations-*.txt",
        default='/home/student.unimelb.edu.au/dxm/observations-chrisjordan/')


    args = parser.parse_args()

    outliers = []
    tiledict = {}

    for filename in glob.glob('%s/observations-*.txt' %args.dir):

        tileloader = TileDataAnalysiswithClustering.DataAnalyserByTile()
        tileloader.load_observationsfromFile(filename)

        clusters = tileloader.TileDBScanClustering(float(args.epsilon))
        for index, cluster in enumerate(clusters):
            if len(cluster) < 4:
                outliers.extend(cluster)
                print(cluster)
                for c in cluster:
                    if c in tiledict:
                        tiledict[c].append(filename[len('%s/observations'%args.dir):-4])
                    else:
                        tiledict[c] = [filename[len('%s/observations'%args.dir):-4]]

                name= "Tile:" + filename[len('%s/observations'%args.dir):-4] + '-' + str(index)
                tileloader.plotCluster(cluster, None, name, True)
        print(tiledict)

    print(tiledict)

    counts = dict(Counter(outliers))

    names = list(counts.keys())
    values = list(counts.values())

    plt.bar(range(len(counts)), values, align='center')
    plt.xticks(range(len(counts)), names, rotation='vertical')
    plt.show()
