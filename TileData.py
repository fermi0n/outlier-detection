import glob
#from pylab import *
import argparse
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

class TileData:
    def __init__(self):
        self.obs_list = []
#        self.all_xx = [None] * 128
#        self.all_yy = [None] * 128

        self.allx_obs_dict = {}
        self.ally_obs_dict = {}

        self.scaled_x_dict = {}
        self.scaled_y_dict = {}

        self.minval = 0.0
        self.maxval = 100.0

        self.problem_obs_list = []

        self.discretisedAmps = {}
        self.metadata_dict = {}

    def load_observations(self, basepath = '.', obsids = None):

        for obs in obsids:
            try:
                cals_loader = CalsLoader.CalsLoader()
                cals_loader.obtainAmplitudeforObservation(basepath + '/'+ obs)
                print("Loaded from cals_loader")
                self.allx_obs_dict[obs] = cals_loader.JPX
                self.ally_obs_dict[obs] = cals_loader.JQY
                self.obs_list.append(obs)
                self.discretisedAmps[obs] = [None]*128
                self.metadata_dict[obs] = cals_loader.metadata

            except Exception, err:
                print("Error loading observation %s, Error is: %s"%(obs, err))
                self.problem_obs_list.append(obs)

        print("Problem Obs IDs are ", self.problem_obs_list)
        print("Good Obs IDs are ", self.obs_list)

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

    def saveObservations(self, filename):
        #Print all the observations to csv files in a simple format:
        #ObsId, antenna, channel, x-amp, Y-amp
        f= open(filename,"w+")

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
        f.close()
        f = open((filename.split('.')[0])+".meta", "w+")
        for obs in self.obs_list:
            f.write("%s,%s,%s" %(obs,self.metadata_dict[obs][0], self.metadata_dict[obs][1]))
            f.write("\n")
        f.close()

    #load observations from file written in the format given in the above method
    def load_observations_from_file(self, filename):

        f = open(filename, "r")

        for line in f:
            values = line.split(',')
            if values[0] not in self.obs_list:
                self.obs_list.append(values[0])
                self.allx_obs_dict[values[0]] = [None]*128
                self.ally_obs_dict[values[0]] = [None]*128
            if values[2] == "X":
                self.allx_obs_dict[values[0]][int(values[1])] = [float(i) for i in values[3:-1]]
            else:
                self.ally_obs_dict[values[0]][int(values[1])] = [float(i) for i in values[3:-1]]
        self.obs_list.sort()
        f.close()

    def load_metafits_from_file(self, filename):

        f = open(filename, "r")
        #Don't do the flagged dipoles just yet
        for line in f:
            values = line.split(',')
            self.metadata_dict[values[0]] = [values[1], values[2]]
        f.close()

if __name__ == "__main__":

    tileloader = TileData()

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help="Save amplitude data in this file")
    parser.add_argument('-p', '--path', help="Base path where the observations are kept", default='/lustre/projects/p048_astro/MWA/data')
    #parser.add_argument('-s', '--save', action='store_true', help="Save plot to file rather than display")
    parser.add_argument("observations", nargs='*', help="Observation IDs")
    args = parser.parse_args()
    print("Loading all these observations into memory:")
    print(args.observations)
    if (len(args.observations) > 20):
        print ("Usage: Don't provide more than 20 observations. Later need to think about allowing this and splitting the load function into multiples")
        sys.exit(0)
    for obs in args.observations:
	tileloader.load_observations(basepath=args.path, obsids=[obs])
    	print("Saving observations to file")
    	print(args.file)
    	tileloader.saveObservations(args.file)
#    tileloader.load_observations(basepath=args.path, obsids=args.observations)
    #print("Saving observations to file")
    #print(args.file)
    #tileloader.saveObservations(args.file)

    #print(tileloader.allx_obs_dict['1061313248'][126])
    #knn = DataAnalysis.KNNAnalysis()
    #histogramdata, distances = knn.KShapeHistogram(tileloader)

    #print(distances[:30])

    #plt.hist(histogramdata, bins=100) #, bins=list(np.arange(0.0, 0.1, 0.001)))
    #plt.show()




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
