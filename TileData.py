import glob
from pylab import *
#matplotlib.use('Agg')

import numpy as np
import DataAnalysis
#import matplotlib.pyplot as plt
import math
#import CalsLoader
import datetime
import time
#import vatcode as vat
#import seaborn as sns
#from scipy.spatial.distance import pdist, squareform

obsIDs = [
'1061311664',
'1061312392',
'1061312520',
'1061312880',
'1061313008',
'1061313248',
'1061312392',
'1061312520',
'1061312760',
'1061312880',
'1061313008',
'1061313248',
'1061313496',
'1061313736',
'1061313856',
'1061313984',
'1061314104',
'1061314224',
'1061314344',
'1061314472',
'1061314592',
'1061314712',
'1061314832',
'1061314960',
'1061316296']


class TileData:
    def __init__(self):
        self.obs_list = []
        self.all_xx = [None] * 128
        self.all_yy = [None] * 128

        self.allx_obs_dict = {}
        self.ally_obs_dict = {}

        self.scaled_x_dict = {}
        self.scaled_y_dict = {}

        self.minval = 0.0
        self.maxval = 100.0

        self.discretisedAmps = {}

    def load_observations(self, basepath = './', obsids = None):

        for obs in obsids:
            cals_loader = CalsLoader.CalsLoader()
            cals_loader.obtainAmplitudeforObservation(basepath + '/'+ obs)

            self.allx_obs_dict[obs] = cals_loader.JPX
            self.ally_obs_dict[obs] = cals_loader.JQY

            self.discretisedAmps[obs] = [None]*128

        self.obs_list = obsids

    def scaleObservations(self):

        #First calculate the maximum of the x and y values, in order to do minMax scaling
        maxval = 0.0
        minval = 100.0
        for key, value in self.allx_obs_dict.iteritems():
            for x in value:
                observationmax = max(x)
                observationmin = min([100] + [a for a in x if a != 0])
                if observationmax > maxval:
                    maxval = observationmax
                if observationmin < minval:
                    minval = observationmin


        for key, value in self.ally_obs_dict.iteritems():
            for y in value:
                observationmax = max(y)
                observationmin = min([100] + [a for a in y if a != 0])
                if observationmax > maxval:
                    maxval = observationmax
                if observationmin < minval:
                    minval = observationmin

        self.minval = minval
        self.maxval = maxval

        for key, value in self.allx_obs_dict.iteritems():
            self.scaled_x_dict[key] = [None]*128
            for index, x in enumerate(value):
                self.scaled_x_dict[key][index] = [(a - minval) / (maxval - minval) for a in x]

        for key, value in self.ally_obs_dict.iteritems():
            self.scaled_y_dict[key] = [None]*128
            for index, y in enumerate(value):
                self.scaled_y_dict[key][index] = [(a - minval) / (maxval - minval) for a in y]

    def saveObservations(self):
        #Print all the observations to csv files in a simple format:
        #ObsId, antenna, channel, x-amp, Y-amp
        f= open("observations.txt","w+")

        for obs in self.obs_list:
            for antenna in range(128):
                if self.allx_obs_dict[obs][antenna] is not None:
                    f.write("%s,%s,X, " %(obs, antenna))
                    for channel, value in enumerate(self.allx_obs_dict[obs][antenna]):
                        f.write("%f, "%self.allx_obs_dict[obs][antenna][channel])
                    f.write("\n")
                    f.write("%s,%s,Y, " %(obs, antenna))
                    for channel, value in enumerate(self.ally_obs_dict[obs][antenna]):
                        f.write("%f, "%self.ally_obs_dict[obs][antenna][channel])
                    f.write("\n")

    #load observations from file written in the format given in the above method
    def loadObservationsFromFile(self, filename):

        f = open(filename, "r")

        for line in f:
            values = line.split(',')
            if (values[0] not in self.obs_list):
                self.obs_list.append(values[0])
                self.allx_obs_dict[values[0]] = [None]*128
                self.ally_obs_dict[values[0]] = [None]*128
            if values[2] == "X":
                self.allx_obs_dict[values[0]][int(values[1])] = [float(i) for i in values[3:-1]]
            else:
                self.ally_obs_dict[values[0]][int(values[1])] = [float(i) for i in values[3:-1]]
        self.obs_list.sort()


if __name__ == "__main__":

    tileloader = TileData()
    tileloader.loadObservationsFromFile("observationsmodified.txt")

    #print(tileloader.allx_obs_dict['1061313248'][126])
    knn = DataAnalysis.KNNAnalysis()
    histogramdata, distances = knn.KShapeHistogram(tileloader)

    print(distances[:30])

    plt.hist(histogramdata, bins=100) #, bins=list(np.arange(0.0, 0.1, 0.001)))
    plt.show()

#    ax=sns.heatmap(R,vmax=0.1, xticklabels=False,yticklabels=False)
#    ax.set(xlabel='Objects', ylabel='Objects')
#    plt.show(ax)

#tileloader.load_observations(basepath='/lustre/projects/p048_astro/MWA/data', obsids = obsIDs)
#tileloader.saveObservations()
#tileloader.scaleObservations()
#tileloader.discretiseAmplitudes()
#print tileloader.CalculateKNearestNeighbourScore(10, obsIDs[0], 20)
#
# knn = knearestneighbour.knearestneighbour()
#
# scores = []
# for obs in obsIDs:
#     for i in range (128):
# #        score = knn.CalculateKNNScoreEuclidean(10, obs, i, tileloader.scaled_x_dict, obsIDs)
#         score = knn.CalculateKNNScoreDistanceCorrelation(10, obs, i, tileloader.scaled_x_dict, obsIDs) + knn.CalculateKNNScoreDistanceCorrelation(10, obs, i, tileloader.scaled_y_dict, obsIDs)
#         scores.append([obs, i, score])
#
# scores.sort(key=takethird, reverse=True)
# print "k-NN Scores are:"
# print scores
#
# #tileloader.calculate_steepness('1061311664', 20)
# #tileloader.identify_features()
#    tileloader.plotObservation(obs=obsIDs[0])
# tileloader.plotFlaggedTiles(scores[:10])
