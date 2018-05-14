import matplotlib
matplotlib.use('Agg')
#Uncomment below line to use DTW
#import mlpy as mlpy
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss
import scipy.spatial.distance as ssd
from core import kshape, zscore, _sbd
from operator import itemgetter
import DataAnalysis as da
import CalsLoader
import argparse
import pickle


class DataAnalyserByTile:

    def __init__(self):
        #Note: In this class, this dict is not a double array but a single array (as it only loads for single tile, not for all tiles)
        self.allx_obs_dict = {}
        self.ally_obs_dict = {}
        self.obs_list = []
        self.problem_obs_list = []
        self.metadata_dict = {}
        self.distances = []

    def TileKShapeHistogram(self, tile):

        #First create an array of tile amplitudes
        listoftiles = []
        print self.obs_list
        for obs in self.obs_list:
            listoftiles.append([obs, tile] + [self.allx_obs_dict[obs] + self.ally_obs_dict[obs]])

        #Then create a matrix of tile->tile distances
        R = []
        for counter, (obs, tileindex, tilevalue) in enumerate(listoftiles):
            #            listofscores = []
            #Do this cleverly so we don't need to calculate kshapedistance twice for each pair, but we still retain a symmetric matrix
            listofscores = [R[i][counter] for i in range(counter)]
            for obs2, tileindex2, tilevalue2 in listoftiles[counter:]:
                listofscores.append([obs, tileindex, obs2, tileindex2] + [da.DistanceMeasures.KShapeDistance(tilevalue, tilevalue2)])
            R.append(listofscores)

        print (R)

        #Then sum the distances of each tile to other tiles (Note: this will give global anomalies, not local - i.e. won't take clustering into account)
        distances = []
        for line in R:
            summation = 0
            for obs, tileindex, _, _, distance in line:
                summation = summation + distance
            distances.append([obs, tileindex, summation])

        distances.sort(key=itemgetter(2), reverse=True)
        #print(distances)

        self.distances = distances

        #Just get the summation values for purposes of histogram plotting
        histogramlist = [summation for obs, index, summation in distances]  #or could do histogramlist = list(map(itemgetter(3), distances)) I think
        return histogramlist, distances

    def load_observations(self, tile, obsids, basepath):

        for obs in obsids:
            try:
                cals_loader = CalsLoader.CalsLoader()
                cals_loader.obtainAmplitudeforObservation(basepath + '/'+ obs)
                print("Obtained amps")
                self.allx_obs_dict[obs] = cals_loader.JPX[tile]
                self.ally_obs_dict[obs] = cals_loader.JQY[tile]
                print("Set amps to TileDataAnalysis class")
                self.obs_list.append(obs)
                self.metadata_dict[obs] = cals_loader.metadata
                print("Set metadata dict")

            except Exception, err:
                print("Error loading observation ", obs, "Error is:", err)
                self.problem_obs_list.append(obs)

        print("Problem Obs IDs are ", self.problem_obs_list)
        print("Good Obs IDs are ", self.obs_list)


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


    def loadandsaveObservations(self, tile, obsids, basepath, filename):

    	for obs in obsids:
            try:
                f = open(filename, "a")
                cals_loader = CalsLoader.CalsLoader()
                cals_loader.obtainAmplitudeforObservation(basepath + '/' + obs)
                if cals_loader.JPX[tile] is not None:
            #        print cals_loader.JPX[tile]
                    f.write("%s,%s,X," %(obs, tile))
                    for value in cals_loader.JPX[tile]:
                        f.write("%f," %value)
                    f.write("\n")
                if cals_loader.JQY[tile] is not None:
                    f.write("%s,%s,Y," %(obs, tile))
                    for value in cals_loader.JQY[tile]:
                        f.write("%f," %value)
                    f.write("\n")
                f.flush()
                f.close()
            except Exception, err:
                print("Error loading observation ", obs, "Error is:", err)
                self.problem_obs_list.append(obs)



if __name__ == '__main__':
    dm = da.DistanceMeasures()
    tileloader = DataAnalyserByTile()

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help="Load amplitude data from this file")
    parser.add_argument('-p', '--path', help="Base path where the observations are kept", default='/lustre/projects/p048_astro/MWA/data')
    parser.add_argument('-t', '--tile', help="Tile to perform analysis on", type=int)
    parser.add_argument('-o', '--output', help="Basepath and filename where to save the observations", default='/lustre/projects/p048_astro/dmendonca/observations.txt')
    parser.add_argument("observations", nargs='*', help="Observation IDs to load")

    args = parser.parse_args()

    if (args.file == None):
        if (args.observations == None):
            print("Usage: Must provide list of observations to assess on")
            sys.exit(0)
        else:
            tileloader.loadandsaveObservations(int(args.tile), args.observations, args.path, args.output)

    else:
        #tileloader.loadObservationsFromFile(args.file)  Need to update this to not load all tile data, just the tile we are interested in.
        print("Usage: Not available yet - TODO")


#    histogramdata, distances = tileloader.TileKShapeHistogram(args.tile)

#    print("Printing all distances")
#    print(distances)

    #Pickle this data
    #with open ('%s/%s_data.pickle'%(args.output, args.tile), 'wb') as f:
    #    pickle.dump(tileloader, f, pickle.HIGHEST_PROTOCOL)

    # plt.hist(histogramdata, bins=500) #, bins=list(np.arange(0.0, 0.1, 0.001)))
    # plt.savefig('/lustre/projects/p048_astro/dmendonca/%s_Tile_%s.png'%(time.strftime("%d%b_%H:%M:%S", time.localtime()),args.tile))
    #
    # for obsid, tile, distance in distances[:127]:
    #
    #     ax = fig.add_subplot(8,16,sp+1,)
    #     ax.plot(tileloader.allx_obs_dict[obsid], colours[1])
    #     ax.plot(tileloader.ally_obs_dict[obsid], colours[2])
    #     plt.title('Obs %s\n, distance %.1f'%(obsid, distance), fontsize=6)
    #
    #     if ppl != 16:
    #         plt.setp(ax.get_xticklabels(), visible=False) # plot setup
    #         plt.setp(ax.get_yticklabels(), visible=False)
    #
    #     if ppl == 16:
    #         ppl = 0
    #         plt.setp(ax.get_xticklabels(), visible=False)
    #
    #     ppl += 1
    #
    #     plt.ylim([-0.1,maxv])
    #
    #     if sp == 15:
    #         XX_amp, = ax.plot([],[], colours[1],label='XX',linewidth=3.0)
    #         YY_amp, = ax.plot([],[], colours[2],label='YY',linewidth=3.0)
    #         ax.legend((XX_amp,YY_amp),('XX','YY'), bbox_to_anchor=(0, 2, .12, .12),prop={'size':14})
    #
    #     sp += 1
    #
    #     plt.tight_layout()
    #     fig.subplots_adjust(top=0.9)
    #     plt.suptitle('Amps | Tile %s' %(tile),fontsize=18)
    #
    # plt.savefig('/lustre/projects/p048_astro/dmendonca/Tile_%s_TopDistances.png'%(args.tile))
