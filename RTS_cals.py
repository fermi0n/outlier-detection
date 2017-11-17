from numpy import cos, sin,array, matrix, inner, conj, real, imag,zeros,arange,shape, sqrt, ones
from astropy.io import fits

class rts_cal():
    def __init__(self):
        self.n_ants = 128
        self.n_bands = 24
        self.antennas=[]
        self.flagged_antennas = None
        for i in range(self.n_ants):
            self.antennas.append(rts_antenna(i))
        self.freqs = None
        self.metafile = None
        self.BP_weights = [None] * self.n_bands


    # Each BP file contains all antennas, so load must occur at top level

    def load_BP_jones(self,band,path='./',raw=False,add_flagged=True,chan_bw=0.04):
        calsfile = path + 'BandpassCalibration_node%03d.dat' % band
        cals_file = open(calsfile)

        present_freqs = ((cals_file.readline()).split(','))
        present_freqs = [int(round(float(p)/chan_bw)) for p in present_freqs]
        flagged_channels = array([i for i in range(32) if i not in present_freqs])
        all_gains = []

        all_antennas = []
        for line in cals_file:
            gains = line.split(',')
            antenna = int(gains[0])
            all_antennas.append(antenna-1)
            amp_gains = array(map(float,gains[1::2]))
            phase_gains = array(map(float,gains[2::2]))
            # These gains are in (surely) amp-phase, so let's convert to real,imag
#            r_gains = [float(amp_gains[i]) * cos(float(phase_gains[i])) for i in range(len(amp_gains))]
#            i_gains = [float(amp_gains[i]) * sin(float(phase_gains[i])) for i in range(len(amp_gains))]
#            c_gains = [r + 1.0j * i for r,i in zip(r_gains,i_gains)]
            r_gains = amp_gains * cos(phase_gains)
            i_gains = amp_gains * sin(phase_gains)
            c_gains = r_gains + 1.0j * i_gains
            all_gains.append(c_gains)

        if(self.flagged_antennas is None):
            present_antennas = set(all_antennas)
            self.flagged_antennas = [a for a in range(128) if a not in present_antennas]

        # gains are in order (data, fits) for XX, XY, YX, YY

        if(raw==True):
            xx_gains = all_gains[::8]
            xy_gains = all_gains[2::8]
            yx_gains = all_gains[4::8]
            yy_gains = all_gains[6::8]
        else:
            xx_gains = [all_gains[i] for i in range(len(all_gains)) if i%8 == 1]
            xy_gains = [all_gains[i] for i in range(len(all_gains)) if i%8 == 3]
            yx_gains = [all_gains[i] for i in range(len(all_gains)) if i%8 == 5]
            yy_gains = [all_gains[i] for i in range(len(all_gains)) if i%8 == 7]

        antenna_number = [all_antennas[i] for i in range(len(all_antennas)) if i%8 == 0]

        # Need to treat separate frequency channels

        all_jones = []
        flagged_jones = matrix([[0,0],[0,0]])

        for i in range(len(xx_gains)):
            jones = [matrix([[xx,xy],[yx,yy]]) for xx,xy,yx,yy in zip(xx_gains[i],xy_gains[i],yx_gains[i],yy_gains[i])]
            self.antennas[antenna_number[i]].BP_jones[band-1] = jones
            if(add_flagged):
                # Add in flagged channels
                for ch in range(32):
                    if(ch in flagged_channels):
                        self.antennas[antenna_number[i]].BP_jones[band-1].insert(ch,flagged_jones)    
                self.BP_weights[band-1] = ones(32)
                self.BP_weights[band-1][flagged_channels] = 0.0
            else:
                self.BP_weights[band-1] = ones(len(present_channels))

#                self.antennas[antenna_number[i]].BP_jones[band-1].insert(0,flagged_jones)
#                self.antennas[antenna_number[i]].BP_jones[band-1].insert(0,flagged_jones)
#                self.antennas[antenna_number[i]].BP_jones[band-1].append(flagged_jones)
#                self.antennas[antenna_number[i]].BP_jones[band-1].append(flagged_jones)
                # center channel
#                self.antennas[antenna_number[i]].BP_jones[band-1].insert(16,flagged_jones)

            
#            all_jones.append(jones)

        cals_file.close()

    def load_DI_jones(self,band,path='./'):
        calsfile = path + 'DI_JonesMatrices_node%03d.dat' % band
        cals_file = open(calsfile)
        flux_scale = float(cals_file.readline())
        all_jones = []

        for line in cals_file:
            xxr, xxi, xyr, xyi, yxr, yxi, yyr, yyi = line.split(',')
            jones = matrix([[float(xxr) + 1.0j * float(xxi),float(xyr) + 1.0j * float(xyi)],[float(yxr) + 1.0j * float(yxi),float(yyr) + 1.0j * float(yyi)]]) 
            all_jones.append(jones)        
        cals_file.close()

        inv_ref0 = all_jones[0].I

        DI_jones = [j * inv_ref0 for j in all_jones[1:]]        

        for a,dj in zip(self.antennas,DI_jones):
            a.DI_jones[band-1] = dj

    def average_BP_jones(self):
        for idx,a in enumerate(self.antennas):
            if(a!=None):
                self.antennas[idx].average_BP_jones()
        freqs = [(self.freqs[2*i] + self.freqs[2*i+1])/2.0 for i in range(len(self.freqs)/2)]
        self.freqs = freqs

    def form_single_jones(self):
        for idx,a in enumerate(self.antennas):
            if(a!=None):
                self.antennas[idx].form_single_jones()

    def load_all_BP_jones(self,path='./',raw=False,add_flagged=True):
        n_bands = 24
        for i in range(1,n_bands+1):
            self.load_BP_jones(i,path=path,raw=raw,add_flagged=add_flagged)

    def load_all_DI_jones(self,path='./'):
        n_bands = 24
        for i in range(1,n_bands+1):
            self.load_DI_jones(i,path=path)

    def load_metadata(self,metafile):
        from astropy.time import Time
        self.metafile = metafile
        meta_file = fits.open(metafile)
        self.meta_gains = (meta_file[1].data)['Gains']
        t = Time(meta_file[0].header['DATE-OBS'])
        self.mjd = 2.4e6 + t.mjd
        # set frequencies across the whole band
        nchans = meta_file[0].header['nchans']
        fine_chan = meta_file[0].header['FINECHAN']
        freqcent = meta_file[0].header['FREQCENT']

        freqs = (fine_chan * 1e3) * arange(nchans)
        self.freqs = freqs + freqcent * 1e6 - (freqs[nchans/2] + freqs[nchans/2 -1]) / 2.0 - (fine_chan * 1e3) / 2.0 
        meta_file.close()
 

    def apply_metafits_gains(self):
        for i in range(self.n_ants):
            g0 = float(self.meta_gains[2*i,0])
            norm_gains = [float(g)/g0 for g in self.meta_gains[2*i]]
            d = [di / g for di,g in zip(self.antennas[i].DI_jones,reversed(norm_gains))]
            self.antennas[i].DI_jones = d
        


    def all_BP_jones(self,antenna=0,reverse_bands=False,cent_chan=12):
        all_BP_xx = []
        all_BP_yy = []
        a = self.antennas[antenna]
        if(a!=None):
            for b,B in enumerate(a.BP_jones):
                band_number = b + cent_chan - 12
                if(band_number < 129):
                   if(B is not None):
                        for chan in B:
                            all_BP_xx.append(chan[0,0]) 
                else:
                    if(a.BP_jones[self.n_bands-(band_number-128)] is not None):
                        for chan in (a.BP_jones[self.n_bands-(band_number-128)]):
                            all_BP_xx.append(chan[0,0])
                

#            if(reverse_bands):
#                for B in reversed(a.BP_jones):
#                    if(B is not None):
#                        for chan in B:
#                            all_BP_xx.append(chan[0,0])
#            else:
#                for B in a.BP_jones:
#                    if(B is not None):
#                        for chan in B:
#                            all_BP_xx.append(chan[0,0])
        return all_BP_xx

    def BP_jones_amps(self,antenna=0):
        bp_jones_amps = []
        a = self.antennas[antenna]
        for B in a.BP_jones:
            bp_band_amps = []
            if(B!=None):
                for chan in B:
                    bp_band_amps.append(abs(chan))
            bp_jones_amps.append(bp_band_amps)
            
        return bp_jones_amps

    def all_single_jones(self,antenna=0,reverse_bands=False,correct_gain_jump=True,conjugate_jones=True,pol='xx'):
        all_single_xx = []
        all_single_yy = []
        a = self.antennas[antenna]
        if(a!=None):
            if(reverse_bands):
                for idx,B in enumerate(reversed(a.Single_jones)):
                    if(B!=None):
                        for chan in B:
                            if(conjugate_jones):
                                chan_j = inner(chan,conj(chan))
                            else:
                                chan_j = chan
                            if(correct_gain_jump):
                                if(idx >=16):
                                    all_single_xx.append(0.25 * chan_j[0,0])
                                    all_single_yy.append(0.25 * chan_j[1,1])
                                else:
                                    all_single_xx.append(chan_j[0,0])
                                    all_single_yy.append(chan_j[1,1])
                            else:
                               all_single_xx.append(chan_j[0,0]) 
                               all_single_yy.append(chan_j[1,1])
            else:
                for B in a.Single_jones:
                    if(B!=None):
                        for chan in B:
                            all_single_xx.append(chan[0,0])
                            all_single_yy.append(chan[1,1])
        
        return all_single_xx, all_single_yy

    def write_UCLA_cal(self,filename=None,reverse_bands=True):

        if(self.metafile is None):
            print "Error: Use load_metadata to define metafits file"
        else:    
            meta_file = fits.open(self.metafile)
        antenna_n = meta_file[1].data['Antenna']
        tile_n = meta_file[1].data['Tile']
        obsid = meta_file[0].header['GPSTIME']

        if(filename is None):
            out_file = open('RTS_%s_ucla.txt' % obsid,'w+')
        else:
            out_file = open(filename,'w+')

        # convert Cotter/FHD antenna number into rts antenna number
        ant2rts = dict(zip(antenna_n[::2],range(self.n_ants)))
        ant2tile = dict(zip(antenna_n[::2],tile_n[::2]))

        out_file.write("#Program of origin: RTS\n")
        out_file.write("# Convention: Divide uncalibrated data by these gains to obtain calibrated data.\n")

        pols = [['EE','EN'],['NE','NN']]

    #    for ai,a in enumerate(self.antennas):
        for p1 in [0,1]:
            for p2 in [0,1]:
                n_bands = 24
                n_chan = 16
                if(reverse_bands is True):
                    bands = reversed(range(n_bands))
                else:
                    bands = range(n_bands)
                for i_n,i in enumerate(bands):
                    for j in range(n_chan):
                        for k in range(self.n_ants):
                           if(ant2rts[k] not in self.flagged_antennas):
                               flag = 0
                               single_jones = self.antennas[ant2rts[k]].Single_jones[i][j]
                           else:
                               flag = 1
                               single_jones = zeros((2,2))
                               
                           out_file.write("%s %s %6.3f %s %.5f %f %f %d\n" % (ant2tile[k],k,self.freqs[i_n*n_chan+ j]/1.0e6,pols[p1][p2],self.mjd,real(single_jones[p1,p2]),imag(single_jones[p1,p2]),flag))

#                for i in range(self.n_ants):
#                    if(ant2rts[i] not in self.flagged_antennas):
#                        flag = 0 
#                        single_jones = self.antennas[ant2rts[i]].Single_jones
#                    else:
#                        flag = 1
#                        single_jones = zeros((24,16,2,2))
#                    for B in single_jones: 
#                        if(B!=None):
#                            for chan in B:
#                                out_file.write("%s %s %6.3f %s %.5f %f %f %d\n" % (i,tile_n[2*i],freq,pols[p1][p2],self.mjd,real(chan[p1,p2]),imag(chan[p1,p2]),flag))
        out_file.close()

class rts_antenna():
    def __init__(self,ant_i,n_bands=24):
        self.ant_n = ant_i
        self.n_bands = n_bands
        self.BP_jones = [None] * n_bands
        self.DI_jones = [None] * n_bands
        self.Single_jones = [None] * n_bands

    def average_BP_jones(self,fscrunch=2):
        for idx,B in enumerate(self.BP_jones):
            if(B!=None):
                avg_jones = [(B[2*i] + B[2*i+1])/2.0 for i in range(len(B)/2)]
                avg_jones[len(avg_jones)/2] *= 2.0
                self.BP_jones[idx] = avg_jones

    def form_single_jones(self):
        for idx,[B,DI] in enumerate(zip(self.BP_jones,self.DI_jones)):
            if(B!=None):
                self.Single_jones[idx] = [DI * bp for bp in B]

                

def write_BP_files(raw_cal,fit_cal,filename='test'):
    from cmath import phase
    n_bands = 24
    n_tiles = 128
    flagged_channels = [0,1,16,30,31]
    bp_freqs = "0.080000, 0.120000, 0.160000, 0.200000, 0.240000, 0.280000, 0.320000, 0.360000, 0.400000, 0.440000, 0.480000, 0.520000, 0.560000, 0.600000, 0.680000, 0.720000, 0.760000, 0.800000, 0.840000, 0.880000, 0.920000, 0.960000, 1.000000, 1.040000, 1.080000, 1.120000, 1.160000\n"
    for n in range(n_bands):
        band_file = filename + '%03d.dat' % (n+1)
        fp = open(band_file,'w+')
        fp.write(bp_freqs)
        for i in range(n_tiles):
            if(raw_cal.antennas[i] is not None):
                if(raw_cal.antennas[i].BP_jones[n] is not None):
                    for pol1 in [0,1]:
                        for pol2 in [0,1]:
                            fp.write('%d' % (i+1))
                            for ch_n,ch in enumerate(raw_cal.antennas[i].BP_jones[n]):
                                if(ch_n not in flagged_channels):
                                    fp.write(', %f, %f' % (abs(ch[pol1,pol2]), phase(ch[pol1,pol2])))
                            fp.write('\n')
                            fp.write('%d' % (i+1))
                            for ch_n,ch in enumerate(fit_cal.antennas[i].BP_jones[n]):
                                if(ch_n not in flagged_channels):
                                    fp.write(', %f, %f' % (abs(ch[pol1,pol2]), phase(ch[pol1,pol2])))
                            fp.write('\n')
                    
        fp.close()

                
            
        
        
            
        
        


            
