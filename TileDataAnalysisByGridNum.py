import numpy as np
import scipy.stats as ss
import scipy.spatial.distance as ssd
from core import kshape, zscore, _sbd
from operator import itemgetter
import DataAnalysis as da
import CalsLoader
from TileDataAnalysis import DataAnalyserByTile
#Uncomment below line to use DTW
#import mlpy as mlpy
import matplotlib.pyplot as plt
import argparse
import pickle


if __name__ == '__main__':

    fig = plt.figure(figsize=(18.0, 10.0))

    #Unpickle this data
    tile = '4'

    with open ('/lustre/projects/p048_astro/dmendonca/%s_data.pickle'%(tile), 'rb') as f:
        tileloader = pickle.load(f)

    #For each gridnum, plot a box-and-whiskers plot of each observation distances to all OTHER observations
    #List of unique gridnums
    uniquegridnums = []

    #First create an array of tile amplitudes
    listoftiles = []
    #print self.obs_list
    for obs in tileloader.obs_list:
        listoftiles.append([obs, tile, tileloader.metadata_dict[obs][0]] + [tileloader.allx_obs_dict[obs] + tileloader.ally_obs_dict[obs]])
        if tileloader.metadata_dict[obs][0] not in uniquegridnums:
            uniquegridnums.append(tileloader.metadata_dict[obs][0])
    #Then create a matrix of tile->tile distances
    R = []
    for counter, (obs, tileindex, gridnum, tilevalue) in enumerate(listoftiles):

        #Do this cleverly so we don't need to calculate kshapedistance twice for each pair, but we still retain a symmetric matrix
        listofscores = [R[i][counter] for i in range(counter)]
        for obs2, tileindex2, gridnum2, tilevalue2 in listoftiles[counter:]:
            listofscores.append([obs, tileindex, gridnum, obs2, tileindex2, gridnum2] + [da.DistanceMeasures.KShapeDistance(tilevalue, tilevalue2)])
        R.append(listofscores)

    print (R)

    #OK, now grab the tile->tile distances that have the same gridnum
    distances = []
    for line in R:
        intrasums = 0
        extrasums = 0
        for obs, tileindex, gridnum, _, _, gridnum2, distance in line:
            if (gridnum == gridnum2):
                intrasums = intrasums + distance
            else:
                extrasums = extrasums + distance
        distances.append([obs, tileindex, intrasums, extrasums])

    #distances.sort(key=itemgetter(2), reverse=True)
    #histogramlist = [summation for obs, index, summation in distances]  #or could do histogramlist = list(map(itemgetter(3), distances)) I think
    #fig, ax1 = plt.subplots(figsize=(2, 1))
    #fig.canvas.set_window_title('Box Plot of Tile %s with intrasums and extrasums' %tile)
    #plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

    #bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
#    plt.setp(bp['boxes'], color='black')
#plt.setp(bp['whiskers'], color='black')
#plt.setp(bp['fliers'], color='red', marker='+')

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
#ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
#               alpha=0.5)

# Hide these grid behind plot objects
#ax1.set_axisbelow(True)
#ax1.set_title('Comparison of IID Bootstrap Resampling Across Five Distributions')
#ax1.set_xlabel('Distribution')
#ax1.set_ylabel('Value')



    # plot box plots
    plt.boxplot([[intrasums for _, _, intrasums, _ in distances],[extrasums for _, _, _, extrasums in distances]], vert=True)
    #plt.boxplot([extrasums for _, _, _, extrasums in distances])
    #plt.set_title('box plot')
    plt.show()

    print(distances)
    print(uniquegridnums)
