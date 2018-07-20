import numpy as np
from operator import itemgetter
import DataAnalysis as da
import RTS_calibration_data as rcd
import argparse
import pickle
import itertools
import time
#from sklearn.neighbors import LocalOutlierFactor
import pandas as pd

class Flagger:

    time_stretch = 1.0/240.0   #Each observation is 0.5 apart
    freq_stretch = 0.5 #Each channel is 0.5 apart
    RADIUS = 10.0
    EXPECTEDVALS = 2
    GRAD = 0.4

    def __init__(self):
        self.obs_list = []
        self.problem_obs_list = []
        self.NaN_list = []
        self.distances = []
        self.distancematrix_dict = {}


    def calc_density_matrix(self, R):
        #R is the matrix of distances between points.
        #Need to calculate, for each point, the sum of distances to neighbouring points. Since I've already set the ones I don't care about to zero, can just sum up each row!!
        densities = {}
        for key, value in R:
            densities[key] = sum(value.values())
        return densities

    def flag_measurements(self, Rx, Ry, single_tile=True):
        #Note: this function will perform the flagging of measurements according to outlier detection methods
        #clf = LocalOutlierFactor(n_neighbours = 10, metric=self.calc_distances)
        #Calculate a density matrix for Rx and Ry
        densityx = self.calc_density_matrix(Rx)
        densityy = self.calc_density_matrix(Ry)
        #Then find any outlier densities
        average = np.mean(densityx.values())
        stdev = np.std(densityx.values())
        xoutliers = [x for x in densityx if (abs(x-average) >= 2.0*stdev)]
        print xoutliers
        average = np.mean(densityy.values())
        stdev = np.std(densityy.values())
        youtliers = [y for y in densityy if (abs(y-average) >= 2.0*stdev)]
        print youtliers
        return xoutliers, youtliers

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
                amplitudes.append([(int(obs) - int(min_obs))*self.time_stretch] + [channel*self.freq_stretch] + [cals_data.allx_obs_dict[obs][channel]] + [cals_data.ally_obs_dict[obs][channel]])

        print ("NaN list is")
        print(self.NaN_list)
        print("Obs list is")
        print(self.obs_list)
        print("First ten amplitudes are:")
        print(amplitudes[:10])
        print("Length of amplitudes is %d" %len(amplitudes))

        #amplitudes = amplitudes[:200]
        #Then create a matrix of amp->amp distances
        Rx = {}
        Ry = {}
        for (obs, channel, _, _) in amplitudes:
            Rx[str(obs)+':'+str(channel)] = {}
            Ry[str(obs)+':'+str(channel)] = {}
        for counter, (obs, channel, x_amp, y_amp) in enumerate(amplitudes):
            #print ("Up to %d" %counter)
            for obs2, channel2, x_amp2, y_amp2 in amplitudes[counter:]:
                distancex, distancey = self.calc_distances([obs, channel, x_amp, y_amp], [obs2, channel2, x_amp2, y_amp2])
                Rx[str(obs)+':'+str(channel)][str(obs2)+':'+ str(channel2)] = Rx[str(obs2)+':'+ str(channel2)][str(obs)+':'+str(channel)] = distancex
                Ry[str(obs)+':'+str(channel)][str(obs2)+':'+ str(channel2)] = Ry[str(obs2)+':'+ str(channel2)][str(obs)+':'+str(channel)] = distancey

        return Rx, Ry

    def calc_distances(self, vector_one, vector_two):
        #Assume each vector is four numbers - time (obs ID), channel (frequency), x-amplitude and y-amplitude, with the time and channel already stretched
        #Return the distance between the x's and distance between the y's
        #TODO: Can speed this up slightly
        time_dist = vector_one[0] - vector_two[0]
        freq_dist = vector_one[1] - vector_two[1]
        if ((time_dist**2 + freq_dist**2) >= self.RADIUS):
            return 0.0, 0.0
        else:
            x_distance = np.sqrt(time_dist**2 + freq_dist**2 + (vector_one[2] - vector_two[2])**2)
            y_distance = np.sqrt(time_dist**2 + freq_dist**2 + (vector_one[3] - vector_two[3])**2)
            return x_distance, y_distance


    def find_outliers(self, calgains_data):

        #Convert things to a numpy 2D array rather than the dict of lists
        #listoflists = map(list, calgains.allx_obs_dict.items())
        #listoflists = [[int(key), value] for key, value in calgains.allx_obs_dict.items()]
        #print(listoflists)
        #df = pd.DataFrame([calgains.allx_obs_dict])
        #print (df)
        sortedkeys = sorted(calgains.allx_obs_dict.iterkeys())
        flatlist = [calgains.allx_obs_dict[key] for key in sortedkeys]
        #print(flatlist)
        arraylist = np.array(flatlist)
        #print(arraylist)
        updatedsortedkeys = [(int(key) - int(sortedkeys[0])) / 248.0 for key in sortedkeys]
        #print(updatedsortedkeys)

        outliers = []

        #SIMPLIFIED MODEL!!! JUST FIND AREAS WHERE THE GRADIENT IN THE TIME OR THE FREQ DIRECTION IS > SOME values
        #First check if all lengths are the same
        listoflengths = [len(rows) for rows in flatlist]
        #print(listoflengths)
        for i in range(1, len(updatedsortedkeys)-1):
            for j in range(1, len(arraylist[i])-1):
                #print ("%d, %d"%(i,j))
                if abs(arraylist[i][j] - arraylist[i][j-1]) > self.GRAD:
                    outliers.append([i, j, "Freq", sortedkeys[i], updatedsortedkeys[i], arraylist[i][j], abs(arraylist[i][j] - arraylist[i][j-1])])
                if (j < len(arraylist[i-1])-1) and (abs((arraylist[i][j] - arraylist[i-1][j]) / (updatedsortedkeys[i] - updatedsortedkeys[i-1])) > self.GRAD):
                    outliers.append([i, j, "Time", sortedkeys[i], updatedsortedkeys[i], arraylist[i][j], abs(arraylist[i][j] - arraylist[i-1][j])])

        new_outliers = sorted(outliers, key=itemgetter(6), reverse=True)
        print(new_outliers)
        #Only try nearest 5 items
        # #Stop once we have identified two n_neighbours
        # for i in range(2, len(sortedkeys)-2):
        #     for j in range(2, len(arraylist[i])-2):
        #         sum = 0
        #         for counterj in [1,2]:
        #             if np.sqrt((arraylist[i][j-counterj] - arraylist[i][j])**2 + counterj**2) < self.RADIUS:
        #                 sum += 1
        #             if np.sqrt((arraylist[i][j+counterj] - arraylist[i][j])**2 + counterj**2) < self.RADIUS:
        #                 sum += 1
        #         for counteri in [1,2]:
        #             if (np.sqrt((arraylist[i-counteri][j] - arraylist[i][j])**2 + (updatedsortedkeys[i] - updatedsortedkeys[i-counteri])**2) < self.RADIUS:)
        #                 sum += 1
        #             if (np.sqrt((arraylist[i+counteri][j] - arraylist[i][j])**2 + (updatedsortedkeys[i] - updatedsortedkeys[i+counteri])**2) < self.RADIUS:)
        #                 sum += 1
        #         if (sum < 2):
        #
        #             outliers.append([i, j])
        #
        # print(outliers)

        #listoflists = sorted([[int(key), value] for key, value in calgains.allx_obs_dict.items()])
        #print(listoflists)
        #xarray = [np.ravel(x) for x in listoflists]
        #print(xarray)
        #print(xarray)
        #np.array(list(calgains.allx_obs_dict)))

if __name__ == '__main__':

    calgains = rcd.CalGains_data()

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help="Load amplitude data from this file")

    args = parser.parse_args()

    calgains.load_observationsfromFile(args.file) #assume single_tile is true

    flagger = Flagger()
    flagger.find_outliers(calgains)
    #Rx, Ry = flagger.DistanceMatrixCalculator(calgains)
    #xoutlier, youtlier = flagger.flag_measurements(Rx, Ry)
    #print("Xoutliers are:")
    #print(xoutlier)
    #print("Youtlier are:")
    #print(youtlier)
