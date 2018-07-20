import matplotlib
#matplotlib.use('Agg')
#Uncomment below line to use DTW
#import mlpy as mlpy
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss
from sklearn.cluster import DBSCAN
import scipy.spatial.distance as ssd
from core import kshape, zscore, _sbd
from operator import itemgetter
import DataAnalysis as da
import TileDataAnalysiswithClustering
#from TileDataAnalysis import DataAnalyserByTile
#import CalsLoader
import argparse
import pickle
import itertools
import time
import pandas as pd


if __name__ == '__main__':
    dm = da.DistanceMeasures()
    tileloader = TileDataAnalysiswithClustering.DataAnalyserByTile()

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help="Load amplitude data from this file")
    #parser.add_argument('-n', '--numclusters', help="Number of clusters to plot", type=int, default=2)
    parser.add_argument('-t', '--tile', help="Tile to perform analysis on", type=int)
    parser.add_argument('-e', '--epsilon', help="Value of epsilon, the density range", default=0.005, type=float)

    args = parser.parse_args()

    tileloader.load_observationsfromFile(args.file)

    clusters = tileloader.TileDBScanClustering(float(args.epsilon))

    name= "Tile:" + str(args.tile) + "Cluster:"

    for index, cluster in enumerate(clusters):
        tileloader.plotCluster(cluster, None, name + str(index), True)
        for obs in cluster:
            print("%s,%s"%(obs, index))
