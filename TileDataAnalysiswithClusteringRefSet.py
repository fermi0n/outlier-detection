import matplotlib
#matplotlib.use('Agg')
#Uncomment below line to use DTW
#import mlpy as mlpy
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss
import scipy.spatial.distance as ssd
from core import kshape, zscore, _sbd
from operator import itemgetter
import DataAnalysis as da
from TileDataAnalysis import DataAnalyserByTile
import CalsLoader
import argparse
import pickle
import itertools
import time
import TileData as tl

##Reference set obsIDs: 1061316296 (high-band) and 1127246440 (low-band)

class DataAnalyserByTile:

    def __init__(self):
        #Note: In this class, this dict is not a double array but a single array (as it only loads for single tile, not for all tiles)
        self.allx_obs_dict = {}
        self.ally_obs_dict = {}
        self.obs_list = []
        self.problem_obs_list = []
        self.NaN_list = []
        self.metadata_dict = {}
        self.distances = []
        self.distancematrix_dict = {}
        self.listofreferences = []
        self.save = False

    def Linkage(self, setA, setB, distancematrix_dict):
        #Set A is one set of obsIDs, set B is the next set
        #distancematrix_dict is a dict of dicts
        #Use the average linkage algorithm, described in https://en.wikipedia.org/wiki/Hierarchical_clustering
        sum = 0
        for A in setA:
            for B in setB:
                sum = sum + distancematrix_dict[A][B]
        average = sum / float(len(A) * len(B))
        return average

    def AverageDistanceObsToCluster(self, obsID, setObs, R):
        #SetObs is the set against which we are calculating this quantity.
        #Use algorithm described in P.J. RoFseeuw: Graphical aid to cluster analysis.
        #This function will be used to calculate a(i) and d(i)
        sum = 0
        #Another possible Error
        if (len(setObs) < 1):
            print("Serious issue - set obs has less than one element")
            print(setObs)
        for obs in setObs:
            sum = sum + R[obsID][obs]
            if (np.isnan(sum)):
                print ("sum is nan, distancematrix is%.2f, observation pair is %s, %s"%(R[obsID][obs], obsID, obs))
        average = sum / float(len(setObs))
        return average

    def CalculateSilhouettes (self, clusterlist, R):

        silhouettes = []
        #For each cluster:
        for cluster in clusterlist:
            #Get the rest of the clusters
            #otherclusters = [clust in clusterlist where clust != cluster]
            #For each obs in that cluster:
            for obs in cluster:
                #Calculate distances to each other cluster, find min
                #This is b(i)
                neighbourdistance = min([self.AverageDistanceObsToCluster(obs, clust, R) for clust in clusterlist if clust != cluster])
                #Calculate a(i)
                WithinClusterAverageDistance = self.AverageDistanceObsToCluster(obs, cluster, R)
                #Calculate s(i):
                s = (neighbourdistance - WithinClusterAverageDistance) / float(max(neighbourdistance, WithinClusterAverageDistance))
                if np.isnan(s):
                    print("neighbourdistance = %.2f, withinclusterdistance = %.2f"%(neighbourdistance, WithinClusterAverageDistance))
                silhouettes.append((obs, cluster, s))
        #Now we have the list of all s(i)'s!
        return silhouettes


    def InitialiseTileList(self):

        #First create an array of tile amplitudes
        listoftiles = []
        print self.obs_list
        for obs in self.obs_list:
            if len([x for x in self.allx_obs_dict[obs] if np.isnan(x)]) > 0:
                print("Obs %s has NaNs"%obs)
                self.NaN_list.append(obs)
            if len([x for x in self.ally_obs_dict[obs] if np.isnan(x)]) > 0:
                if (obs not in self.NaN_list):
                    print("Obs %s has NaNs"%obs)
                    self.NaN_list.append(obs)

        for obs in self.NaN_list:
            self.obs_list.remove(obs)

        for obs in self.obs_list:
            listoftiles.append([obs] +  [self.allx_obs_dict[obs] + self.ally_obs_dict[obs]])
        print ("NaN list is")
        print(self.NaN_list)
        print("Obs list is")
        print(self.obs_list)
        print("List of tiles is")
        print([obs for obs, value in listoftiles])

        return listoftiles


    def CreateDistanceMatrix(self, listoftiles):
        #Create a matrix of tile->tile distances
        R = {}
        for (obs, tilevalue) in listoftiles:
            R[obs] = {}
        for counter, (obs, tilevalue) in enumerate(listoftiles):
            #            listofscores = []
            #Do this cleverly so we don't need to calculate kshapedistance twice for each pair, but we still retain a symmetric matrix
            #listofscores = [R[i][counter] for i in range(counter)]
            for obs2,tilevalue2 in listoftiles[counter:]:
                distance = da.DistanceMeasures.KShapeDistance(tilevalue, tilevalue2)
                if (np.isnan(distance)):
                    print("NaN distance discovered for observations %s, %s"%(obs, obs2))
                    exit(0)
                R[obs][obs2] = R[obs2][obs] = distance

#        print (R)
        return R
        #self.distancematrix_dict = R

    def TileKShapeClusteringA(self, listoftiles):

        R = self.CreateDistanceMatrix(listoftiles)

        clusterlist = [[obs] for obs, _ in listoftiles]
        print ("Initial clusterlist is")
        print(clusterlist)
        listofsilhouettes = {}
        listofclusters = {}
        while (len(clusterlist) >= 2):

            combinations = [[obs1, obs2] for counter, obs1 in enumerate(clusterlist) for obs2 in clusterlist[counter+1:]]
            #combinations = [[obs1, obs2] ]
            #print("Combinations are")
            #print(combinations)
            distances = [self.Linkage(c[0], c[1], R) for c in combinations]
            mincombination = combinations[distances.index(min(distances))]

            clusterlist.remove(mincombination[0])
            clusterlist.remove(mincombination[1])
            #Add mincombination to clusterlist
            clusterlist.append(mincombination[0] + mincombination[1])

            #Calculate the silhouettes, only for 2 to 10
            if (1 < len(clusterlist) < 10): #Note that we have already done a clusterlist reduction
                #Calculate silhouettes
                silhouettes = self.CalculateSilhouettes(clusterlist, R)
                listofsilhouettes[len(clusterlist)] = silhouettes
                listofclusters[len(clusterlist)] = list(clusterlist)

        #Calculate averages
        # silhouettestats = []
        # for i in range(2, 10):
        #     average_s = np.mean([s for (_, _, s) in listofsilhouettes[i]])
        #     silhouettestats.append([i, average_s])
        #
        # #Print average_s
        # print("Printing silhouette statistics")
        # print(silhouettestats)

        return listofsilhouettes, listofclusters




    def load_observations(self, tile, obsids, basepath):

        for obs in obsids:
            try:
                cals_loader = CalsLoader.CalsLoader()
                cals_loader.obtainAmplitudeforObservation(basepath + '/'+ obs)

                self.allx_obs_dict[obs] = cals_loader.JPX[tile]
                self.ally_obs_dict[obs] = cals_loader.JQY[tile]
                self.obs_list.append(obs)
                self.metadata_dict[obs] = cals_loader.metadata

            except Exception, err:
                print("Error loading observation ", obs, "Error is:", err)
                self.problem_obs_list.append(obs)

        print("Problem Obs IDs are ", self.problem_obs_list)
        print("Good Obs IDs are ", self.obs_list)


    def plotCluster(self, cluster, distancesdict, clustername, saveflag=False):

        colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

        sp = 0
        ppl = 16

    #    maxv = 10.0
        maxv = max(max([max(tiledata.allx_obs_dict[obs][ii]) for ii in range(128)]), max([max(tiledata.ally_obs_dict[obs][ii]) for ii in range(128)]))

        fig = plt.figure(figsize=(18.0, 10.0))

        observations = []
        if distancesdict != None:
            observations = [(obsid, distancesdict[obsid][0], distancesdict[obsid][1]) for obsid in cluster]
            observations.sort(key=itemgetter(1))
        else:  #To handle the NaN case
            observations = [(obs, 0, None) for obs in cluster]

        for obsid, s, refob in observations[:128]:

            ax = fig.add_subplot(8,16,sp+1,)
            ax.plot(self.allx_obs_dict[obsid], colours[1])
            ax.plot(self.ally_obs_dict[obsid], colours[2])
            plt.title('Obs %s\n, d = %.2f, %s'%(obsid, s*100, refob), fontsize=6)

            if ppl != 16:
                plt.setp(ax.get_xticklabels(), visible=False) # plot setup
                plt.setp(ax.get_yticklabels(), visible=False)

            if ppl == 16:
                ppl = 0
                plt.setp(ax.get_xticklabels(), visible=False)

            ppl += 1

            plt.ylim([-0.1,maxv])

            if sp == 15:
                XX_amp, = ax.plot([],[], colours[1],label='XX',linewidth=3.0)
                YY_amp, = ax.plot([],[], colours[2],label='YY',linewidth=3.0)
                ax.legend((XX_amp,YY_amp),('XX','YY'), bbox_to_anchor=(0, 2, .12, .12),prop={'size':14})

            sp += 1

            plt.tight_layout()
            fig.subplots_adjust(top=0.9)
            plt.suptitle('Amps | %s' %(clustername),fontsize=18)

        if (saveflag):
            plt.savefig('%s_%s.png'%(clustername, time.strftime("%d%b_%H:%M:%S", time.localtime())))
        else:
            plt.show()

    def plotDistancesHeatMap(self, listofclusters, listoftiles, name=''):

        #Create 2x2 matrix (instead of the current dict)
        Matrix = []
        labels = []
        xaxislabels = ['']
        yaxislabels = ['']
        xticklocations = [0,]
        for index, cluster in enumerate(listofclusters):
            labels.extend(cluster)
            #xaxislabels = ([''] * (len(cluster)-1)) + ['|']
            #yaxislabels = ([''] * (len(cluster)-1)) + ['-']
            prevloc = xticklocations[-1]
            xticklocations.append(prevloc+len(cluster)/2)
            xticklocations.append(prevloc+len(cluster))
            xaxislabels.append('Cluster %d'%index)
            yaxislabels.append('Cluster %d'%index)
            xaxislabels.append('')
            yaxislabels.append('')
        #xaxislabels=['Cluster %d'%i, '|' for i in range(len(cluster))]
        #xaxislabels=['|']*len(xticklocations)
        #xaxislabels=[''] + ['Cluster %d |'%i for i in range(len(xticklocations)-1)]
        #yaxislabels = ['-'] * len(xticklocations)
        #print("Labels are")
        #print(labels)
        R = self.CreateDistanceMatrix(listoftiles)

    #    labels = labels[len(cluster[0])+ len(cluster[1]):]
        print("Clusters")
        print(listofclusters)
        print("Labels")
        print (labels)

        for obs in labels:
            distances = []
            for obs2 in labels:
                distances.append(R[obs][obs2])
                #distances.insert(0, R[obs][obs2])
            Matrix.append(distances)
        print("Matrix")
        print(Matrix)
        data = np.array(Matrix)

        #print(data)
        #print(data.ndim)
        fig, ax = plt.subplots(1, 1)
        cax = ax.imshow(np.array(Matrix), cmap='viridis', interpolation='nearest', origin='lower')
        ax.set_xticks(xticklocations)
        ax.set_yticks(xticklocations)
        ax.set_xticklabels(xaxislabels, rotation=30)
        ax.set_yticklabels(yaxislabels)
        ax.tick_params(axis='both', which='both', direction='inout')
        ax.set_title(name)
        fig.colorbar(cax, orientation='vertical')
        if (self.save):
            plt.savefig('%s_%s.png'%(name, time.strftime("%d%b_%H:%M:%S", time.localtime())))
        else:
            plt.show()


if __name__ == '__main__':
    dm = da.DistanceMeasures()
    tileloader = DataAnalyserByTile()

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help="Load amplitude data from this file")
    parser.add_argument('-p', '--path', help="Base path where the observations are kept", default='/lustre/projects/p048_astro/MWA/data')
    parser.add_argument('-n', '--numgoodclusters', help="Number of good clusters to plot", type=int, default=2)
    parser.add_argument('-o', '--outlierclusters', help="Number of outlier clusters to plot", type=int, default=4)
    parser.add_argument('-t', '--tile', help="Tile to perform analysis on", type=int)
    parser.add_argument('-s', '--save', help="If this flag is on, will save to file", action='store_true')
    parser.add_argument("observations", nargs='*', help="Observation IDs to load")

    args = parser.parse_args()

    tileloader.save = args.save

    if (args.file == None):
        if (args.observations == None):
            print("Usage: Must provide list of observations to assess on")
            sys.exit(0)
        else:
            tileloader.load_observations(args.tile, args.observations, args.path)

    else:

        with open (args.file, 'rb') as f:
            tileload = pickle.load(f)

        tileloader.allx_obs_dict = tileload.allx_obs_dict
        tileloader.ally_obs_dict = tileload.ally_obs_dict
        tileloader.obs_list = tileload.obs_list

    print("Pickle file loaded, now loading reference data")
    #Now, load the reference setObs into self.listofeferences
    goldenobsids = ['1061316296', '1127246440']
    referencetileloader = tl.TileData()
    referencetileloader.load_observations(basepath='/lustre/projects/p048_astro/MWA/data', obsids=[goldenobsids[0]])
    referencetileloader.load_observations(basepath='/lustre/projects/p006_astro/MWA/data', obsids=[goldenobsids[1]])

    for i in range(128):
        for gobid in goldenobsids:
            tileloader.listofreferences.append([gobid + '_' + str(i), referencetileloader.allx_obs_dict[gobid][i] + referencetileloader.ally_obs_dict[gobid][i]])

    print("Reference tiles loaded")
    #Now, split the outliers and the good observations
    #Look for observations which are > threshold distance from any of the reference set.

    #Initialise tile list
    listoftiles = tileloader.InitialiseTileList()

    listofoutliers = []
    listofgoodobs = []
    THRESHOLD = 0.005
    dictofdistances = {}

    #Evaluate min distance to the reference set.
    for (obs, tilevalue) in listoftiles:
        distance = 1
        tile = None
        for tileindex2, tilevalue2 in tileloader.listofreferences:
            newdistance = da.DistanceMeasures.KShapeDistance(tilevalue, tilevalue2)
            if (newdistance < distance):
                distance = newdistance
                tile = tileindex2[-5:]
        dictofdistances[obs] = [distance, tile]
        if distance > THRESHOLD:
            print("Outlier %s"%obs)
            listofoutliers.append([obs, tilevalue])
        else:
            listofgoodobs.append([obs, tilevalue])


    plt.hist([distance for distance, _ in dictofdistances.values()], 500)
    plt.title('Distribution of minimum distances to reference set')

    if (args.save):
        plt.savefig('Hist_Tile:%d_%s.png'%(args.tile, time.strftime("%d%b_%H:%M:%S", time.localtime())))
    else:
        plt.show()

    #Now, cluster the good observations
    listofgoodsilhouettes, listofgoodclusters = tileloader.TileKShapeClusteringA(listofgoodobs)

    #Now, cluster the outliers.
    listofoutliersilhouettes, listofoutlierclusters = tileloader.TileKShapeClusteringA(listofoutliers)

    #now plot these
    if (args.tile==None):
        name="Good Cluster:"
    else:
        name= "Tile:" + str(args.tile) + ",Good Cluster:"

    for index, cluster in enumerate(listofgoodclusters[args.numgoodclusters]):
        #tileloader.plotCluster(cluster, listofgoodsilhouettes[args.numgoodclusters], name + str(index), False)
        tileloader.plotCluster(cluster, dictofdistances, name + str(index), args.save)

    if (args.tile==None):
        name="Outlier Cluster:"
    else:
        name= "Tile:" + str(args.tile) + ",Outlier Cluster:"

    for index, cluster in enumerate(listofoutlierclusters[args.outlierclusters]):
        tileloader.plotCluster(cluster, dictofdistances, name + str(index), args.save)

    tileloader.plotCluster(tileloader.NaN_list, None, "Tile:" + str(args.tile) + ",Cluster:NaNs", args.save)

    #tileloader.plotDistancesHeatMap(listofgoodclusters[args.numgoodclusters] + listofoutlierclusters[args.outlierclusters], listoftiles)
    tileloader.plotDistancesHeatMap(listofoutlierclusters[args.outlierclusters], listoftiles, name='%d-Outliers_v_outliers'%args.tile)
    tileloader.plotDistancesHeatMap(listofgoodclusters[args.numgoodclusters], listoftiles, name='%d-Good_obs_v_good_obs'%args.tile)
