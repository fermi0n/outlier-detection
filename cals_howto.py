from pylab import *
import os,sys
import RTS_cals
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def euclideandistance(instance1, instance2, length):
    distance = 0
    for x in range(length):
        distance += (instance1[x] - instance2[x])**2
    return np.sqrt(distance)

rts_cal = RTS_cals.rts_cal()

data_dir = sys.argv[1]
meta_file = fits.open(sys.argv[2])

cent_chan = meta_file[0].header['CENTCHAN']

rts_cal.load_all_BP_jones(path=data_dir, raw=True)
rts_cal.load_all_DI_jones(path=data_dir)
rts_cal.form_single_jones()

all_xx = [None] * 128
all_yy = [None] * 128

Data = []
xx_all = []
yy_all = []

FlaggedData = []
K = 5
RADIUS = 0.01

def flagdatapoint(index):
    FlaggedData.append(Data[index])
    print 'Data flagged'
    print Data[index]

for i,a in enumerate(rts_cal.antennas):
    if(rts_cal.antennas[i].BP_jones[0] is not None):
        xx_abs,yy_abs = rts_cal.all_single_jones(antenna=i,reverse_bands=True,correct_gain_jump=False,conjugate_jones=False)
        xx_abs = real((xx_abs * conj(xx_abs)))
        yy_abs = real(yy_abs * conj(yy_abs))
        all_xx[i] = xx_abs
        #print ("Single value of xx_abs is", xx_abs)
        all_yy[i] = yy_abs
        for j, a in enumerate(xx_abs):
            if not xx_abs[j] == 0 and not yy_abs[j] == 0:
                Data.append([xx_abs[j], yy_abs[j], i, j])
                xx_all.append(xx_abs[j])
                yy_all.append(yy_abs[j])

plt.scatter(xx_all, yy_all)
plt.show()

print 'starting the flagging'
for index, value in enumerate(Data):
    print "Up to value ", index
    neighbours = 0
    for j, v in enumerate(Data):
        if not index == j:
            distance = euclideandistance(value, v, 2)
            if distance <= RADIUS:
                neighbours += 1
            if neighbours >= K:
                #print 'too many neighbours'
                break
    if neighbours < K:
        flagdatapoint(index)
    if index > 20000:
        break

print 'flagging done'

#print "printing all_xx[0]"
#print all_xx[0]

print FlaggedData
flagged_xs = [i[0] for i in FlaggedData]
flagged_ys = [i[1] for i in FlaggedData]
plt.scatter(xx_all, yy_all,color='blue')
plt.scatter(flagged_xs, flagged_ys, color='red')
plt.show()
plt.savefig('FlaggedData.png')

# xs = all_xx[(FlaggedData[0])[2]]
# ys = all_yy[(FlaggedData[0])[2]]
# for index, value in enumerate(xs):
#     if value == 0.0 and index > 0:
#         xs[index] = xs[index-1]
# for index, value in enumerate(ys):
#     if value == 0.0 and index > 0:
#         ys[index] = ys[index-1]


#fig = plt.subplots(nrows=np.ceil(np.sqrt(len(FlaggedData))), ncols=np.ceil(np.sqrt(len(FlaggedData)))

def plotSingleFlaggedTile(axes, r, c, currentindex):
    xs = all_xx[(FlaggedData[currentindex])[2]]
    ys = all_yy[(FlaggedData[currentindex])[2]]
    for index, value in enumerate(xs):
        if value == 0.0 and index > 0:
            xs[index] = xs[index-1]
    for index, value in enumerate(ys):
        if value == 0.0 and index > 0:
            ys[index] = ys[index-1]
    axes[r, c].plot(xs, label='XX abs')
    axes[r, c].plot(ys, label='YY abs')
    axes[r, c].legend(loc='upper right')
    #axes.set_title('Flagged Tile %d', %(r*c))
    axes[r, c].annotate('Flagged point', xy=((FlaggedData[r*c])[3],FlaggedData[r*c][1]-0.02), xytext=((FlaggedData[r*c])[3]+1, 0.1), arrowprops=dict(facecolor='black', shrink=0.05),)

rows = 2
cols = 2
currentantenna = FlaggedData[0][2]
currentindex = 0

fig, axes = plt.subplots(nrows=rows, ncols=cols)  #Sharey=true (consider using)
for r in range(0, rows):
    for c in range(0, cols):
        print 'plotting using value %d' %currentindex
        xs = all_xx[(FlaggedData[currentindex])[2]]
        ys = all_yy[(FlaggedData[currentindex])[2]]
        for index, value in enumerate(xs):
            if value == 0.0 and index > 0:
                xs[index] = xs[index-1]
        for index, value in enumerate(ys):
            if value == 0.0 and index > 0:
                ys[index] = ys[index-1]
        axes[r, c].plot(xs, label='XX abs')
        axes[r, c].plot(ys, label='YY abs')
        axes[r, c].set_title('Flagged Antenna %d' % currentantenna)
        #axes[r, c].legend(loc='upper right')

        for i, v in enumerate(FlaggedData[currentindex:]):
            if (v[2] == currentantenna):
                axes[r, c].annotate('Flag', xy=((FlaggedData[currentindex+i])[3],min(FlaggedData[currentindex+i][0], FlaggedData[currentindex+i][1])-0.05), xytext=((FlaggedData[currentindex+i])[3], 0.1), arrowprops=dict(facecolor='black', width=0.1, headwidth=0.1, shrink=0.01),)
                print 'annotated using value:'
                print v
            else:
                currentantenna = v[2]
                currentindex += i
                print 'breaking out of enumeration loop with currentindex = ', currentindex, 'and currentantenna = ', currentantenna
                print v
                break

fig.tight_layout()
plt.show()

#plt.plot(xs, label='XX abs')
#plt.plot(ys, label='YY abs')
#plt.legend(loc="upper center")
#plt.title('Flagged tile')
#plt.annotate('Flagged point', xy=((FlaggedData[0])[3],FlaggedData[0][1]-0.02), xytext=((FlaggedData[0])[3]+1, 0.1), arrowprops=dict(facecolor='black', shrink=0.05),)
#plt.show()

#plt.savefig('Plot.png')
#print all_xx
