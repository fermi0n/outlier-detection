import matplotlib
#matplotlib.use('Agg')
#Uncomment below line to use DTW
#import mlpy as mlpy
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
#from TileDataAnalysis import DataAnalyserByTile
#import CalsLoader
import argparse
import pickle
import itertools
import time
import pandas as pd
from collections import Counter


if __name__ == '__main__':


    outliers = []
    f = open('obs.csv')

    for line in f:
        values = line.strip(',').split(',')
        outliers.extend(values)

    counts = dict(Counter(outliers))
    counts.pop('')
    counts.pop('\n')
    print(counts)
    names = list(counts.keys())
    values = list(counts.values())

    plt.bar(range(len(counts)), values, align='center')
    plt.xticks(range(len(counts)), names, rotation='vertical')
    plt.show()
