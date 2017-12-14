import glob
from pylab import *
matplotlib.use('Agg')
from astropy.io import fits
from fit_coarse_channels import fit_bandpass, fit_sigma
import numpy as np
import matplotlib.pyplot as plt
import math

import RTS_cals
import CalsLoader

obsIDs = ['1061311664','1061314592']

class tile_loader:
    def __init__(self):
        self.obs_list = None   #not sure if I can use sys.argv here
        self.all_xx = [None] * 128
        self.all_yy = [None] * 128

        self.allx_obs_dict = {}
        self.ally_obs_dict = {}
        self.missing_list = []
        self.flaggedlist = []
        self.steepnesses = []

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

            #Just to check I'm getting the right data from the CalsLoader file.
            print "testing - printing allx_obs_dict[" + obs + "][12]"
            print self.allx_obs_dict[obs][12]

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

        print maxval
        print minval
        print "doing ys"

        for key, value in self.ally_obs_dict.iteritems():
            for y in value:
                observationmax = max(y)
                observationmin = min([100] + [a for a in y if a != 0])
                if observationmax > maxval:
                    maxval = observationmax
                if observationmin < minval:
                    minval = observationmin
        print maxval
        print minval
        self.minval = minval
        self.maxval = maxval

        for key, value in self.allx_obs_dict.iteritems():
            for x in value:
                self.scaled_x_dict[key] = [(a - minval) / (maxval - minval) for a in x]

        for key, value in self.ally_obs_dict.iteritems():
            for y in value:
                self.scaled_y_dict[key] = [(a - minval) / (maxval - minval) for a in y]

        print "printing unscaled vs scaled values"
        print "Allx_obs_dict['1061311664'][12]"
        print self.allx_obs_dict['1061311664'][12]
        print "Scaled version"
        print self.scaled_x_dict['1061311664'][12]

    def convertAmptoSequence(self, number):
        if number == 0: return 'z'
        interval = (self.maxval - self.minval) / 20.0
        sequence = 'abcdefghijklmnopqrst'
        index = int(math.floor((number - self.minval ) / interval))
        if (index < 0 or index > 19):
            print "yikes - we have a problem. index is:" + str(index)
            return 't'
        return sequence[index]

    def discretiseAmplitudes(self):

        #Calculate table of cutoffs
        for obs in self.obs_list:
            for i in range(128):
                xsequence = ''.join(map(self.convertAmptoSequence, self.allx_obs_dict[obs][i]))
                ysequence = ''.join(map(self.convertAmptoSequence, self.ally_obs_dict[obs][i]))
                print xsequence + ysequence
                self.discretisedAmps[obs][i] = xsequence + ysequence


    def HammingDistance(self, sequence1, sequence2):
        #Return the Hamming distance between equal-length sequences"""
        if len(sequence1) != len(sequence2):
            raise ValueError("Undefined for sequences of unequal length")

        return sum(el1 != el2 for el1, el2 in zip(sequence1, sequence2))

    def CalculateKNearestNeighbourScore(self, k, refobs, reftile):
        #Just do brute force, fuck it
        #First calculate all distances to all other sequences
        distances = []
        for obs in self.obs_list:
            for i in range(128):
                distances.append(self.HammingDistance(self.discretisedAmps[refobs][reftile], self.discretisedAmps[obs][i]))

        #Now order those
        distances.sort()
#        print distances
        return distances[k]

    #Calculate how steep the bandpass is
    def calculate_steepness(self, obsid, antenna):

        #Fit a one-dim polynomial (i.e.a line)
        #Calculate how steep it is (i.e. the gradient)
        print "Calculating for obsid = " + str(obsid) + "antenna = " + str(antenna)
        print "Where length of obsid for x is " +  str(len(self.allx_obs_dict[obsid]))
        print "and length of obsid for y is " + str(len(self.ally_obs_dict[obsid]))
        if self.allx_obs_dict[obsid][antenna] is None or self.ally_obs_dict[obsid][antenna] is None:
            print "Can't continue - antenna turns out to be None"
            return (0,0)
        print 'length of dict is %d' %len(self.allx_obs_dict[obsid][antenna][0])
        print 'x dict is'
        print self.allx_obs_dict[obsid][antenna][0]
        zx = np.polyfit(arange(763), self.allx_obs_dict[obsid][antenna][0], 1)
        print zx
        zy = np.polyfit(arange(763), self.ally_obs_dict[obsid][antenna][0], 1)
        print zy
        return (np.real(zx[0]), np.real(zy[0]))

    #Copied in case I screw this up - adding y's to this calculation
    # def calculate_steepness(self, obsid, antenna):
    #
    #     #Fit a one-dim polynomial (i.e.a line)
    #     #Calculate how steep it is (i.e. the gradient)
    #     print "Calculating for obsid = " + str(obsid) + "antenna = " + str(antenna)
    #     print "Where length of obsid for x is " +  str(len(self.allx_obs_dict[obsid]))
    #     print "and length of obsid for y is " + str(len(self.ally_obs_dict[obsid]))
    #     if self.allx_obs_dict[obsid][antenna] is None:
    #         print "Can't continue - antenna turns out to be None"
    #         return 0
    #     print 'length of dict is %d' %len(self.allx_obs_dict[obsid][antenna][0])
    #     print 'dict is'
    #     print self.allx_obs_dict[obsid][antenna][0]
    #     z = np.polyfit(arange(763), self.allx_obs_dict[obsid][antenna][0], 1)
    #     print z
    #     return z[0]

    #This method will run over the observations collected, and calculate a "feature set" for each observation. I can then use that
    #to identify observations with outlier features
    def identify_features(self):

        self.steepnesses = [[self.calculate_steepness(obs, antenna) for antenna in range(128)] for obs in obsIDs]
        steepnesses = self.steepnesses
        print "steepnesses are:"
        print steepnesses
        xs = [x for i in steepnesses for (x,y) in i ]
        ys = [y for i in steepnesses for (x, y) in i ]
        print xs
        print ys

        xstddev = np.std(xs)
        xavg = np.mean(xs)
        ystddev = np.std(ys)
        yavg = np.mean(ys)


        for obsindex, obsvalue in enumerate(obsIDs):
        #    steepnesses = [self.calculate_steepness(obs, antenna) for antenna in range(128)]
            #print "steepnesses are:"
            #print steepnesses

            #Find anything which is more than 2 std devs from average
        #    xs = [x for (x,y) in steepnesses]
        #    ys = [y for (x, y) in steepnesses]

        #    stddev = np.std(xs)
        #    avg = np.mean(xs)
            print "X Stddev is " + str(xstddev) + "and average is " + str(xavg)
            print "Y Stddev is " + str(ystddev) + "and average is " + str(yavg)
            for antennaindex, antennasteepness in enumerate(steepnesses[obsindex]):
                #xsteepness = antennasteepness[0]
                #if steepness is None:
                #    continue
                if np.abs(antennasteepness[0]) >= (np.abs(xavg) + 3.0 * xstddev):
                    print "FLAGGED x steepness is" + str(antennasteepness[0]) + "for obsID " + obsvalue
                    plot(tileloader.allx_obs_dict[obsvalue][antennaindex][0], label='x amps')
                    plot(tileloader.ally_obs_dict[obsvalue][antennaindex][0], label='y amps')
                    title('X amp flagged: Antenna %d for OBSID %s' %(antennaindex, obsvalue))
                    savefig('X_antenna_%s_obs_%s_too_steep.png' % (antennaindex, obsvalue))

                if np.abs(antennasteepness[1]) >= (np.abs(yavg) + 3.0 * ystddev):
                    print "FLAGGED y steepness is" + str(antennasteepness[1]) + "for obsID " + obsvalue
                    plot(tileloader.allx_obs_dict[obsvalue][antennaindex][0], label='x amps')
                    plot(tileloader.ally_obs_dict[obsvalue][antennaindex][0], label='y amps')
                    title('Y amp flagged: Antenna %d for OBSID %s' %(antennaindex, obsvalue))
                    savefig('Y_antenna_%s_obs_%s_too_steep.png' % (antennaindex, obsvalue))
            # stddev = np.std(ys)
            # avg = np.mean(ys)
            # print "Y Stddev is " + str(stddev) + "and average is " + str(avg)
            # for i, steepness in enumerate(ys):
            #     #if steepness is None:
            #     #    continue
            #     if np.abs(np.real(steepness)) >= (np.abs(avg) + 3.0 * stddev):
            #         print "FLAGGED Steepness is" + str(steepness)
            #         plot(tileloader.allx_obs_dict['1061311664'][i][0], label='x amps')
            #         plot(tileloader.ally_obs_dict['1061311664'][i][0], label='y amps')
            #         title('Y amp flagged: Antenna %d for OBSID 1061311664' %(i))
            #         savefig('Y_Antenna_%s_too_steep.png' % (i))




    def plotObservation(self, obs):

        colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

        sp = 0
        ppl = 16
        maxv = 1.5

        fig = plt.figure(figsize=(18.0, 10.0))

        for ii in range(128):

#ax.plot(freq[freq_idx], self.JPX[ii], colours[1])
#ax.plot(freq[freq_idx], self.JQY[ii], colours[2])
#ax.set_xlim(min(freq[freq_idx]), max(freq[freq_idx]))

            ax = fig.add_subplot(8,16,sp+1,)
            ax.plot(self.allx_obs_dict[obs][ii], colours[1])
            ax.plot(self.ally_obs_dict[obs][ii], colours[2])
            plt.title('Tile %d'%ii)

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
        plt.suptitle('Amps | Observation %s' %(obs),fontsize=18)
        plt.show()
		#savefig('Amps_%s.png'%obs, bbox_inches='tight')

def takethird(elem):
    return elem[2]

tileloader = tile_loader()
tileloader.load_observations(basepath='/lustre/projects/p048_astro/MWA/data', obsids = obsIDs)
tileloader.scaleObservations()
tileloader.discretiseAmplitudes()
print tileloader.CalculateKNearestNeighbourScore(10, obsIDs[0], 20)

scores = []
for obs in obsIDs:
    for i in range (128):
        score = tileloader.CalculateKNearestNeighbourScore(10, obs, i)
        scores.append([obs, i, score])

scores.sort(key=takethird, reverse=True)
print "k-NN Scores are:"
print scores

#tileloader.calculate_steepness('1061311664', 20)
#tileloader.identify_features()
tileloader.plotObservation(obs=obsIDs[1])


#Plot a 'good' one:
#plot(tileloader.allx_obs_dict['1061311664'][32][0], label='x amps')
#plot(tileloader.ally_obs_dict['1061311664'][32][0], label='y amps')
#title('Good tile: Antenna 32 for OBSID 1061311664')
#show()
#Let's plot all the ones that are a bit off and have a look-see, shall we
#plot(tileloader.allx_obs_dict['1061311664'][1][0])
#title("Antenna 2 for OBSID 1061311664")
#show()
