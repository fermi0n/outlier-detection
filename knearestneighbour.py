import numpy as np
import scipy.stats as ss
import scipy.spatial.distance as ssd
from core import kshape, zscore, _sbd
from operator import itemgetter
import load_tiles as tl
import mlpy as mlpy

class knearestneighbour:

    def __init__(self):
        #insert some variables here
        sampel = 0

    def HammingDistance(self, sequence1, sequence2):
        #Return the Hamming distance between equal-length sequences"""
        if len(sequence1) != len(sequence2):
            raise ValueError("Undefined for sequences of unequal length")

        return sum(el1 != el2 for el1, el2 in zip(sequence1, sequence2))

    def CalculateKNearestNeighbourScore(self, k, refobs, reftile, discretisedAmps, obs_list):
        #Just do brute force, fuck it
        #First calculate all distances to all other sequences
        distances = []
        for obs in obs_list:
            for i in range(128):
                distances.append(self.HammingDistance(discretisedAmps[refobs][reftile], discretisedAmps[obs][i]))

        #Now order those
        distances.sort()
#        #print distances
        return distances[k]

    def EuclideanDistance(self, list1, list2):

        #print "In Euclidean distance"
        mean1 = np.mean(list1)
        mean2 = np.mean(list2)
        #print mean1
        #print mean2

        difflist = [abs(l1 - mean1 - (l2 - mean2)) for l1, l2 in zip(list1, list2)]
        #print sum(difflist)

        return sum(difflist)

    def DTWDistance(self, list1, list2):

        result = mlpy.dtw_std(list1, list2, dist_only=True)
        print (result)
        return result

    def KShapeDistance(self, list1, list2):

        if (len(list1) != len(list2)):
            print("Unequal lengths")
        result = _sbd(list1, list2)
        ##print(result)
        return result[0]

    def KShapeHistogram(self, tileloader):

        #First create an array of tile amplitudes
        listoftiles = []
        for obs in tileloader.obs_list[:3]:
            listoftiles.extend([[obs, i] + [tileloader.allx_obs_dict[obs][i] + tileloader.ally_obs_dict[obs][i]] for i in range(127)])

        #Then create a matrix of tile->tile distances
        R = []
        for obs, tileindex, tilevalue in listoftiles:
            listofscores = []
            for obs2, tileindex2, tilevalue2 in listoftiles:
                listofscores.append([obs, tileindex, obs2, tileindex2] + [self.DTWDistance(tilevalue, tilevalue2)])
            #Sort this list
            listofscores.sort(key=itemgetter(4))
            R.append(listofscores)

        #print (R)
        #Then sum the first N (say 10 to start with) distances of each tile to other tiles
        distances = []
        for line in R:
            summation = 0
            for obs, tileindex, _, _, distance in line[:9]:
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
                distances.append(self.EuclideanDistance(Amps[refobs][reftile], Amps[obs][i]))

        #Now order these distances
        distances.sort()
        ##print distances
        #print distances
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
        #

# knn = knearestneighbour()
# difference = knn.EuclideanDistance([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
# #print difference
