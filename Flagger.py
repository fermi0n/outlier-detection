import numpy as np
from operator import itemgetter
import DataAnalysis as da
import RTS_calibration_data as rcd
import argparse
import pickle
import itertools
import time
import matplotlib.pyplot as plt
#from sklearn.neighbors import LocalOutlierFactor
import pandas as pd
import glob

class Flagger:

    time_stretch = 1.0/240.0   #Each observation is 0.5 apart
    freq_stretch = 0.5 #Each channel is 0.5 apart
    RADIUS = 10.0
    EXPECTEDVALS = 2
    GRAD = 0.4

    def __init__(self, tile='NS'):
        self.obs_list = []
        self.problem_obs_list = []
        self.NaN_list = []
        self.distances = []
        self.distancematrix_dict = {}
        self.tile = tile


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
        #print xoutliers
        average = np.mean(densityy.values())
        stdev = np.std(densityy.values())
        youtliers = [y for y in densityy if (abs(y-average) >= 2.0*stdev)]
        #print youtliers
        return xoutliers, youtliers

    def DistanceMatrixCalculator(self, cals_data):

        amplitudes = []
        print (cals_data.obs_list)
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


    def find_outliers(self, calgains_data, direction='x'):

        if (direction=='x'):
            obs_dict = calgains.allx_obs_dict
        else:
            obs_dict = calgains.ally_obs_dict

        sortedkeys = sorted(obs_dict.keys())
        flatlist = [obs_dict[key] for key in sortedkeys]
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
        for i in range(1, len(updatedsortedkeys)-1):  #i is the index of the obsIDs
            for j in range(1, len(arraylist[i])-1):  #j is the channel
                #print ("%d, %d"%(i,j))
                if abs(arraylist[i][j] - arraylist[i][j-1]) > self.GRAD:
                    outliers.append([i, j, "Freq", sortedkeys[i], updatedsortedkeys[i], arraylist[i][j], abs(arraylist[i][j] - arraylist[i][j-1])])
                if (j < len(arraylist[i-1])-1) and (abs((arraylist[i][j] - arraylist[i-1][j]) / (updatedsortedkeys[i] - updatedsortedkeys[i-1])) > self.GRAD):
                    outliers.append([i, j, "Time", sortedkeys[i], updatedsortedkeys[i], arraylist[i][j], abs(arraylist[i][j] - arraylist[i-1][j])])

        new_outliers = sorted(outliers, key=itemgetter(6), reverse=True)
        #print(new_outliers)
        #outliers is given as: "1) Index of ObSID; 2) CHannel, 3) "Time" or "Freq" 4) OBSID, 5) Updated OBSID - difference in time between now and first obs, divided by 248
        #6) Gain at this point, 7) Gradient between here and previous point
        return new_outliers, sortedkeys


    def plotOutliers(self, outliers, observations, tiledata):

        for index, (obsindex, channel, typeofissue, obsid, updatedobsid, gain, gradient) in enumerate(outliers[0:10]):

            #First set up the plot using subplot
            # fig = plt.figure()
            # ax11 = fig.add_subplot(221)
            # ax21 = fig.add_subplot(223)
            # ax12 = fig.add_subplot(222, yticklabels=[], sharey=ax11)
            # #ax12.set_ylabel()
            # ax22 = fig.add_subplot(224, sharey=ax21)
            fig, axarray = plt.subplots(2, 2, sharey='row')

            ax11 = axarray[0,0]
            ax12 = axarray[0,1]
            ax21 = axarray[1,0]
            ax22 = axarray[1,1]
            ax11.set_xticks([])
            ax12.set_xticks([])
            #ax, ax2, ax3, ax4) = plt.subplots(2,2)

            plt.subplots_adjust(wspace=0.05, hspace=0.1)
            colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

            timelist = []
            xlist = []
            ylist = []
            Annotate = True
            for index, obs in enumerate(observations):
                if (len(tiledata.allx_obs_dict[str(obs)]) >= channel):
                    try:
                        xlist.append(tiledata.allx_obs_dict[str(obs)][channel])
                        ylist.append(tiledata.ally_obs_dict[str(obs)][channel])
                        timelist.append((float(obs) - float(observations[0])) / 2480000.0)
                    except IndexError as ie:
                        print(ie)
                        print('channel is %d, obs is %s'%(channel, obs))
                        print(tiledata.allx_obs_dict[str(obs)])
                        print ('length is %d' %len(tiledata.allx_obs_dict[str(obs)]))
                elif obsindex >= index - 5:
                    Annotate = False

            ax11.plot(timelist, xlist, color=colours[1], marker='o',label='xx')
            ax11.plot(timelist, ylist, color=colours[2], marker='o', label='yy')
            try:
                ax11.annotate('Outlier', xy=(timelist[obsindex], xlist[obsindex]), xycoords='data', xytext=(0.77, 0.85), textcoords='axes fraction',
                    arrowprops=dict(facecolor='grey', shrink=0.05, width=2, headwidth=5), fontsize=7)
            except IndexError as ie:
                print(ie)
                print("obsindex is %d" %obsindex)
            #ax11.set_title('Tile %s Channel %d' %(self.tile, channel))
            #ax11.set_xlabel("Observation ID (time) in 248-sec chunks from %s" %observations[0])
            #ax11.set_ylabel("Gains")
            #ax1.legend(loc=0)
            #insetax = plt.axes([0.65, .6, .2, .2])
            ax12.plot(timelist[obsindex-2:obsindex+2], xlist[obsindex-2:obsindex+2], color=colours[1],marker='o')
            for i in range(int(obsindex-2), int(obsindex+2)):
                try:
                    ax12.annotate(str(observations[i]), xy=(timelist[i], xlist[i]), fontsize=7)
                except IndexError as ie:
                    print(ie)
                    print("index is %d" %i)

            obs = observations[obsindex]
            for (obs, axis) in [(str(observations[obsindex-1]), ax21), (str(observations[obsindex]), ax22)]:
                maxv = max (max(tiledata.allx_obs_dict[obs], tiledata.ally_obs_dict[obs]))

                axis.plot(tiledata.allx_obs_dict[obs], colours[1], label='XX')
                axis.plot(tiledata.ally_obs_dict[obs], colours[2], label='YY')
                #axis.title('Observation %s' %(str(obs)))
                #axis.set_xlabel('Channel')
                #axis.set_ylabel('Gains Amplitude')
            plt.suptitle('%s outlier at Observation %s, channel %d for tile %s' %(typeofissue, str(obs), channel, self.tile))
            plt.savefig('Outlier-%s-%d.png' %(self.tile, index), bbox_inches='tight')
            plt.close()

if __name__ == '__main__':

    calgains = rcd.CalGains_data()

    parser = argparse.ArgumentParser()
    #parser.add_argument('-f', '--file', help="Load amplitude data from this file")
    #parser.add_argument('-t', '--tile', help="This says what the time number is", default="NS")
    parser.add_argument('-d', '--dir', help="Load amplitude data from this directory. Assume filenames are observations-*.txt",
        default='/home/student.unimelb.edu.au/dxm/observations-chrisjordan/')

    args = parser.parse_args()

    xoutliersdict = {}
    youtliersdict = {}

    outputfile = open('outlierresults.txt', 'a')

    for filename in glob.glob('%s/observations-*.txt' %args.dir):  #i.e. load each tile separately.
        print(filename)
        tile = filename[len(args.dir) + len('/observations'):-4]
        calgains.load_observationsfromFile(filename) #assume single_tile is true

        flagger = Flagger(tile)
    #flagger.tile = args.tile
        outliers, observations = flagger.find_outliers(calgains, 'x')
        xoutliersdict[tile] = (outliers, observations)
        outliers, observations = flagger.find_outliers(calgains, 'y')
        youtliersdict[tile] = (outliers, observations)

        outputfile.write('')


        #flagger.plotOutliers(outliers, observations, calgains)
    #NEED TO sort outliers dict now, and create the heatmap

    #First get the top N outliers (start with all)
    #Then for that,create a 2D numpy array with channel vs time, and the value being the sum (absolute) of the gradients across all tiles
    numpyarray = np.zeros((len(calgains.obs_list), 700))
    #print(outliersdict)
    for outlier_dict in [xoutliersdict, youtliersdict]:  #Note: We aren't printing whether its X or Y - since it doesn't matter, it is just flagged
        for (outliers, observations) in outlier_dict.values():

            for (obsindex, channel, typeofissue, obsid, updatedobsid, gain, gradient) in outliers:
                try:
                    numpyarray[obsindex][channel] += abs(gradient)
                    print ("Added %f" %abs(gradient))
                    outputfile.write("%s, %d\n"%(obsid, channel))

                except IndexError as ie:
                    print(ie)
                    print (obsindex, channel, typeofissue, obsid, updatedobsid, gain, gradient)

        print(numpyarray)

    fig, ax = plt.subplots(1, 1)
    cax = ax.imshow(numpyarray, cmap='viridis', interpolation='nearest', origin='lower')
    ax.set_xlabel('channel')
    ax.set_ylabel('observation index')
    fig.colorbar(cax, orientation='vertical')
    plt.show()
