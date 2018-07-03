import numpy as np
from operator import itemgetter
import DataAnalysis as da
import RTS_calibration_data as rcd
import argparse
import pickle
import itertools
import time

class Flagger:

    time_stretch = 1.0
    freq_stretch = 1.0
    RADIUS = 10.0

    def __init__(self):
        self.obs_list = []
        self.problem_obs_list = []
        self.NaN_list = []
        self.distances = []
        self.distancematrix_dict = {}


    def flag_measurements(self, cals_data, single_tile=True):
        #Note: this function will perform the flagging of measurements according to outlier detection methods
        return 1

    def DistanceMatrixCalculator(self, cals_data):

        amplitudes = []
        print cals_data.obs_list
        self.obs_list = sorted(cals_data.obs_list)

        for obs in cals_data.obs_list:
            if len([x for x in cals_data.allx_obs_dict[obs] if np.isnan(x)]) > 0:
                print("Obs %s has NaNs"%obs)
                self.NaN_list.append(obs)
            if len([x for x in cals_data.ally_obs_dict[obs] if np.isnan(x)]) > 0:
                if (obs not in self.NaN_list):
                    print("Obs %s has NaNs"%obs)
                    self.NaN_list.append(obs)

        for obs in self.NaN_list:
            self.obs_list.remove(obs)

        min_obs = min(self.obs_list)
        for obs in self.obs_list:
            for channel in range(128):
                amplitudes.append([(obs - min_obs)*time_stretch] + [channel*freq_stretch] + [self.allx_obs_dict[obs][channel] + self.ally_obs_dict[obs][channel]])

        print ("NaN list is")
        print(self.NaN_list)
        print("Obs list is")
        print(self.obs_list)
        print("List of tiles is")
        print([obs for obs, value in listoftiles])

        #Then create a matrix of amp->amp distances
        #Only look at amp->amp distances where the channel and time are within a radius of RADIUS - i.e. outside that, set distance to zero
        #TODO: Up to here
        # R = {}
        # for (obs, tilevalue) in listoftiles:
        #     R[obs] = {}
        #
        # for counter, (obs, tilevalue) in enumerate(listoftiles):
        #     #            listofscores = []
        #     #Do this cleverly so we don't need to calculate kshapedistance twice for each pair, but we still retain a symmetric matrix
        #     #listofscores = [R[i][counter] for i in range(counter)]
        #     for obs2,tilevalue2 in listoftiles[counter:]:
        #         distance = da.DistanceMeasures.KShapeDistance(tilevalue, tilevalue2)
        #         if (np.isnan(distance)):
        #             print("NaN distance discovered for observations %s, %s"%(obs, obs2))
        #             exit(0)
        #         R[obs][obs2] = R[obs2][obs] = distance
        #
        # #print (R)
        # self.distancematrix_dict = R

    def distances(self, vector_one, vector_two):
        #Assume each vector is four numbers - time (obs ID), channel (frequency), x-amplitude and y-amplitude, with the time and channel already stretched
        #Return the distance between the x's and distance between the y's
        #TODO: Can speed this up slightly
        time_dist = vector_one[0] - vector_two[0]
        freq_dist = vector_one[1] - vector_two[1]
        if ((time_dist**2 + freq_dist**2) >= RADIUS):
            return 0.0, 0.0
        else:
            x-distance = np.sqrt(time_dist**2 + freq_dist**2 + (vector_one[2] - vector_two[2])**2))
            y-distance = np.sqrt(time_dist**2 + freq_dist**2 + (vector_one[3] - vector_two[3])**2))
            return x-distance, y-distance

if __name__ == '__main__':

    calgains = rcd.CalGains_data()

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help="Load amplitude data from this file")

    args = parser.parse_args()

    calgains.load_observationsfromFile(args.file) #assume single_tile is true

    flagger = Flagger()
    flagger.DistanceMatrixCalculator(calgains)
