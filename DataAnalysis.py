import numpy as np
import scipy.stats as ss
import scipy.spatial.distance as ssd
from core import kshape, zscore, _sbd
from operator import itemgetter
import TileData as tl
#Uncomment below line to use DTW 
#import mlpy as mlpy
import matplotlib.pyplot as plt

class DistanceMeasures:

    def __init__(self):
        #insert some variables here
        sampel = 0

    @staticmethod
    def HammingDistance(sequence1, sequence2):
        #Return the Hamming distance between equal-length sequences"""
        if len(sequence1) != len(sequence2):
            raise ValueError("Undefined for sequences of unequal length")

        return sum(el1 != el2 for el1, el2 in zip(sequence1, sequence2))

    @staticmethod
    def EuclideanDistance(list1, list2):

        #print "In Euclidean distance"
        mean1 = np.mean(list1)
        mean2 = np.mean(list2)
        #print mean1
        #print mean2

        difflist = [abs(l1 - mean1 - (l2 - mean2)) for l1, l2 in zip(list1, list2)]
        #print sum(difflist)

        return sum(difflist)

    @staticmethod
    def DTWDistance(list1, list2):

        #Uncomment below line to run DTW (can't use on gstar)
        #result = mlpy.dtw_std(list1, list2, dist_only=True)
        result = 5
    #    print (result)
        return result

    @staticmethod
    def KShapeDistance(list1, list2):

        if (len(list1) != len(list2)):
            print("Unequal lengths")
        result = _sbd(list1, list2)
        ##print(result)
        return result[0]

class KNNAnalysis:

    def CalculateKNearestNeighbourScore(self, k, refobs, reftile, discretisedAmps, obs_list):
        #Just do brute force, fuck it
        #First calculate all distances to all other sequences
        distances = []
        for obs in obs_list:
            for i in range(128):
                distances.append(DistanceMeasures.HammingDistance(discretisedAmps[refobs][reftile], discretisedAmps[obs][i]))

        #Now order those
        distances.sort()
#        #print distances
        return distances[k]

    def KShapeHistogram(self, tiledata):


        #First create an array of tile amplitudes
        listoftiles = []
        print tileloader.obs_list[-3:]
        for obs in tileloader.obs_list[-3:]:
            listoftiles.extend([[obs, i] + [tileloader.allx_obs_dict[obs][i] + tileloader.ally_obs_dict[obs][i]] for i in range(128)])

        #Then create a matrix of tile->tile distances
        R = []
        for obs, tileindex, tilevalue in listoftiles:
            listofscores = []
            for obs2, tileindex2, tilevalue2 in listoftiles:
                listofscores.append([obs, tileindex, obs2, tileindex2] + [DistanceMeasures.KShapeDistance(tilevalue, tilevalue2)])
            #Sort this list
            listofscores.sort(key=itemgetter(4))
            R.append(listofscores)

        #Then sum the distances of each tile to other tiles (Note: this will give global anomalies, not local - i.e. won't take clustering into account)
        distances = []
        for line in R:
            summation = 0
            for obs, tileindex, _, _, distance in line:
                summation = summation + distance
            distances.append([obs, tileindex, summation])

        distances.sort(key=itemgetter(2), reverse=True)
        #print(distances)

        #Just get the summation values for purposes of histogram plotting
        histogramlist = [summation for obs, index, summation in distances]  #or could do histogramlist = list(map(itemgetter(3), distances)) I think
        return histogramlist, distances


    def CalculateKNNScoreEuclidean(self, k, refobs, reftile, Amps, obs_list):

        #print Amps[refobs][reftile]
        ##print Amps[obs][i]
        distances = []
        for obs in obs_list:
            for i in range(128):
                distances.append(DistanceMeasures.EuclideanDistance(Amps[refobs][reftile], Amps[obs][i]))

        #Now order these distances
        distances.sort()
        return distances[k]


    def CalculateKNNScoreSpearman(self, k, refobs, reftile, Amps, obs_list):

        distances = []
        for obs in obs_list:
            for i in range(128):
                distances.append(ss.spearmanr(Amps[refobs][reftile], Amps[obs][i]))

        #Now order these distances
        distances.sort()
        return sum(distances[:k])

    def CalculateKNNScoreDistanceCorrelation(self, k, refobs, reftile, Amps, obs_list):

        distances = []
        for obs in obs_list:
            for i in range(128):
                distances.append(ssd.correlation(Amps[refobs][reftile], Amps[obs][i]))

        #Now order these distances
        distances.sort()
        return sum(distances[:k])

#    def CalculateDissimilarityMatrix(self, Amps):

# knn = knearestneighbour()
# difference = knn.EuclideanDistance([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
# #print difference
if __name__ == '__main__':
    dm = DistanceMeasures()
    tileloader = tl.TileData()
    knn = KNNAnalysis()
    tileloader.loadObservationsFromFile("observationsmodified.txt")

    histogramdata, distances = knn.KShapeHistogram(tileloader)

    print("Printing largest 30 distances")
    print(distances[:30])
    print("Printing smallest 30 distances")
    print(distances[-10:])

    plt.hist(histogramdata, bins=1000) #, bins=list(np.arange(0.0, 0.1, 0.001)))
    plt.show()
