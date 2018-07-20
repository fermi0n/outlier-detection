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


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--epsilon', help="Value of epsilon, the density range", default=0.005, type=float)

    args = parser.parse_args()

    for i in range(100, 127):
        outliers = []
        tileloader = TileDataAnalysiswithClustering.DataAnalyserByTile()
        tileloader.load_observationsfromFile('/home/student.unimelb.edu.au/dxm/observations-chrisjordan/observations-%d.txt'%i)

        clusters = tileloader.TileDBScanClustering(float(args.epsilon))
        for cluster in clusters:
            if len(cluster) < 4:
                outliers.extend(cluster)
        #print("%d,"%i, end='')
        #print (outliers)
        for outlier in outliers:
            print('%s,'%outlier, end='')
        print('')
