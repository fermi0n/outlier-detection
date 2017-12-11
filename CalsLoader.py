#!/usr/bin/python
#
# script that reads in a TLE file, estimates the power output by the beamformer, and plots the result.

import matplotlib
matplotlib.use('Agg')
import sys, getopt, string, re, os
from astropy.io import fits as pyfits
from pylab import *
import glob
import numpy as np
import matplotlib.pyplot as plt
import math
import optparse
from matplotlib.legend_handler import HandlerLine2D


class CalsLoader:
	def __init__(self):

		self.sel_offset = 1
		self.plot_chan = 0
		self.plot_raw  = 0
		self.obs = '1061314592'
		self.JPX = []
		self.JQY = []
		self.chan_sel = []

	def DI_factors(self, DI_file):
		sel_offset = 1
		plot_chan = 0
		plot_raw  = 0

		try:
			 fid = open(DI_file, 'r')
		except IOError:
			 print 'Can\'t open file \'' + DI_file + '\' for reading.'
			 sys.exit(0)

		# read the file into an array
		lines = fid.readlines()
		# close the file
		fid.close()
		# check that something was read in
		if len(lines)==0:
			print "Erorr reading DI jones matrix file: no lines read from the file."
			sys.exit(0)

		N_ant = len(lines)-2

		J0 = zeros([2,2], "complex")
		G0 = zeros([2,2], "complex")
		J  = zeros([N_ant,2,2], "complex")
		G  = zeros([N_ant,2,2], "complex")
		JJ  = zeros([N_ant,2,2], "complex")

		StokesI = float(lines[0]);

		tmp1 = string.split( lines[1], ',' )
		if len(lines)==0:
			print "Erorr reading DI jones matrix file: line 1 does not appear to have 8 columns."
			sys.exit(0)
		J0[0,0] = float(tmp1[0]) + float(tmp1[1])*1j
		J0[0,1] = float(tmp1[2]) + float(tmp1[3])*1j
		J0[1,0] = float(tmp1[4]) + float(tmp1[5])*1j
		J0[1,1] = float(tmp1[6]) + float(tmp1[7])*1j

		inv_J0 = inv(J0)

		G0 = inner( J0, conj(J0) )

		good_ant = 0

		for ant in range(0, N_ant):
			lineIndex = ant + 2
			tmp1 = string.split( lines[lineIndex], ',' )
			if len(lines)==0:
				print "Erorr reading DI jones matrix file: line %d does not appear to have 8 columns." % lineIndex
				sys.exit(0)

			# this check is not foolproof and can lead to bad things
			# when there is an antenna which should be flagged but is not
	#		if float(tmp1[1])==0 and float(tmp1[3])==0 and float(tmp1[5])==0 and float(tmp1[7])==0:
				#print "flag for antenna %d?" % ant
	#			girolamo = 1
	#		else:
			if(1):

				jones = matrix([[float(tmp1[0]) + float(tmp1[1])*1j,float(tmp1[2]) + float(tmp1[3])*1j],[float(tmp1[4]) + float(tmp1[5])*1j,float(tmp1[6]) + float(tmp1[7])*1j]])
				#jones = jones * inv_J0

				J[good_ant,0,0] = jones[0,0]
				J[good_ant,0,1] = jones[0,1]
				J[good_ant,1,0] = jones[1,0]
				J[good_ant,1,1] = jones[1,1]
				#J[good_ant,0,0] = float(tmp1[0]) + float(tmp1[1])*1j
				#J[good_ant,0,1] = float(tmp1[2]) + float(tmp1[3])*1j
				#J[good_ant,1,0] = float(tmp1[4]) + float(tmp1[5])*1j
				#J[good_ant,1,1] = float(tmp1[6]) + float(tmp1[7])*1j

				G[good_ant] = inner( J[good_ant], conj(J[good_ant]) )
				good_ant = good_ant + 1
		N_ant = good_ant

		factor1 = abs(G0[0,0])/mean(abs(G[:,0,0]))
		factor2 = abs(G0[0,1])/mean(abs(G[:,0,1]))
		factor3 = abs(G0[1,0])/mean(abs(G[:,1,0]))
		factor4 = abs(G0[1,1])/mean(abs(G[:,1,1]))

		fattori = J
		return fattori

	def phase_or_amp(self, sel_offset, data):
		argo,lines,bins,dgains,ggo,flg_t,tiles,calid,JD = data
		#-------------------------------------------------------------
		# initialise the body list
		PX_lsq = []
		PY_lsq = []
		QX_lsq = []
		QY_lsq = []
		PX_fit = []
		PY_fit = []
		QX_fit = []
		QY_fit = []

		tmp1 = string.split( lines[0], ',' )
		tmp2 = []
		for val in tmp1:
			tmp2.append( float(val) )
		N_ch = len(tmp2)

		freq = zeros(N_ch)
		for k in range(0,N_ch):
			freq[k] = tmp2[k]

		if self.plot_raw:
			ch = range(0,N_ch)
		else:
			ch = zeros(N_ch)
			for k in range(0,N_ch):
				ch[k] = freq[k]/0.04

		freq_idx = argsort(freq)

		chan_sel = sel_offset - 1 + 2*array(freq_idx)

		for lineIndex in range(1, len(lines), 8):  # get 0, 8, 16, ...
			tmp = string.split( lines[lineIndex+0], ',' ); PX_lsq.append([]); PX_lsq[-1]=zeros(N_ch*2+1)
			for k in range(0,len(tmp)): PX_lsq[-1][k] = float(tmp[k])
			tmp = string.split( lines[lineIndex+1], ',' ); PX_fit.append([]); PX_fit[-1]=zeros(N_ch*2+1)
			for k in range(0,len(tmp)): PX_fit[-1][k] = float(tmp[k])
			tmp = string.split( lines[lineIndex+2], ',' ); PY_lsq.append([]); PY_lsq[-1]=zeros(N_ch*2+1)
			for k in range(0,len(tmp)): PY_lsq[-1][k] = float(tmp[k])
			tmp = string.split( lines[lineIndex+3], ',' ); PY_fit.append([]); PY_fit[-1]=zeros(N_ch*2+1)
			for k in range(0,len(tmp)): PY_fit[-1][k] = float(tmp[k])
			tmp = string.split( lines[lineIndex+4], ',' ); QX_lsq.append([]); QX_lsq[-1]=zeros(N_ch*2+1)
			for k in range(0,len(tmp)): QX_lsq[-1][k] = float(tmp[k])
			tmp = string.split( lines[lineIndex+5], ',' ); QX_fit.append([]); QX_fit[-1]=zeros(N_ch*2+1)
			for k in range(0,len(tmp)): QX_fit[-1][k] = float(tmp[k])
			tmp = string.split( lines[lineIndex+6], ',' ); QY_lsq.append([]); QY_lsq[-1]=zeros(N_ch*2+1)
			for k in range(0,len(tmp)): QY_lsq[-1][k] = float(tmp[k])
			tmp = string.split( lines[lineIndex+7], ',' ); QY_fit.append([]); QY_fit[-1]=zeros(N_ch*2+1)
			for k in range(0,len(tmp)): QY_fit[-1][k] = float(tmp[k])

		N_ant = len(PX_lsq)

		PXX_lsq = zeros([N_ant,len(PX_lsq[1])-1])
		PYY_lsq = zeros([N_ant,len(PY_lsq[1])-1])
		QXX_lsq = zeros([N_ant,len(QX_lsq[1])-1])
		QYY_lsq = zeros([N_ant,len(QY_lsq[1])-1])

		# rearrange format for the DI_JM

		gg = zeros([len(bins),N_ant,4,2], "float")

		# BP
		# gg and ggo have different number of antennas, leading to some
		# badness from flagged tiles, lets remove flagged tiles from ggo

	        #print 'ggo',ggo[1][126][0,0], ggo[1][127][0,0]


		ggo = np.delete(ggo,[int(f) for f in flg_t if f is not ''],axis=1)

		for sb in range(0,len(bins)):
			for nn in range(0,int(N_ant)):
				gg[sb][nn][0,0] = sqrt(ggo[sb][nn][0,0].real**2 + ggo[sb][nn][0,0].imag**2)
				gg[sb][nn][0,1] = math.atan2(ggo[sb][nn][0,0].imag,ggo[sb][nn][0,0].real)
				gg[sb][nn][1,0] = sqrt(ggo[sb][nn][0,1].real**2 + ggo[sb][nn][0,1].imag**2)
				gg[sb][nn][1,1] = math.atan2(ggo[sb][nn][0,1].imag,ggo[sb][nn][0,1].real)
				gg[sb][nn][2,0] = sqrt(ggo[sb][nn][1,0].real**2 + ggo[sb][nn][1,0].imag**2)
				gg[sb][nn][2,1] = math.atan2(ggo[sb][nn][1,0].imag,ggo[sb][nn][1,0].real)
				gg[sb][nn][3,0] = sqrt(ggo[sb][nn][1,1].real**2 + ggo[sb][nn][1,1].imag**2)
				gg[sb][nn][3,1] = math.atan2(ggo[sb][nn][1,1].imag,ggo[sb][nn][1,1].real)


		# multipling by the DI_JM
		for nn in range(0,int(N_ant)):
			oi = 0
			for rr in range(0,len(bins)):
				for el in range(bins[rr]/2):
					PXX_lsq[nn][oi] = PX_lsq[nn][oi+1]*gg[rr][nn][0,0]
					PXX_lsq[nn][oi+1] = PX_lsq[nn][oi+2]+gg[rr][nn][0,1]
					PYY_lsq[nn][oi] = PY_lsq[nn][oi+1]*gg[rr][nn][1,0]
					PYY_lsq[nn][oi+1] = PY_lsq[nn][oi+2]+gg[rr][nn][1,1]
					QXX_lsq[nn][oi] = QX_lsq[nn][oi+1]*gg[rr][nn][2,0]
					QXX_lsq[nn][oi+1] = QX_lsq[nn][oi+2]+gg[rr][nn][2,1]
					QYY_lsq[nn][oi] = QY_lsq[nn][oi+1]*gg[rr][nn][3,0]
					QYY_lsq[nn][oi+1] = QY_lsq[nn][oi+2]+gg[rr][nn][3,1]
					oi += 2

		PX_lsq = []
		PY_lsq = []
		QX_lsq = []
		QY_lsq = []

		# rearrange the array and update for flagged tiles
		k = 0
		ann = 0

		for nn in range(0,128):
			PX_lsq.append([]); PX_lsq[-1]=zeros(N_ch*2)
			PY_lsq.append([]); PY_lsq[-1]=zeros(N_ch*2)
			QX_lsq.append([]); QX_lsq[-1]=zeros(N_ch*2)
			QY_lsq.append([]); QY_lsq[-1]=zeros(N_ch*2)
			if str(nn) in flg_t:
				finto = 12
			else:
				PX_lsq[-1][:] = PXX_lsq[ann][:]
				PY_lsq[-1][:] = PYY_lsq[ann][:]
				QX_lsq[-1][:] = QXX_lsq[ann][:]
				QY_lsq[-1][:] = QYY_lsq[ann][:]
				ann += 1
		# ===========================
		# test for transpose conjugate of J
		tra = 1
		if tra == 1:
			qtt = len(freq[freq_idx])
			Jp = zeros([128,qtt,2,2], "complex")
			JpH = zeros([128,qtt,2,2], "complex")
			JJ = zeros([128,qtt,2,2], "complex")
			JPX = []
			JQY = []
			for at in range(0,128):
				qq = 0
				gel = 0
				fc = 0
				for qt in range(0,qtt):
					if qt == bins[gel]/2+fc:
						gel += 1
						fc = qt
					Jp[at,qt,0,0] = (PX_lsq[at][qq]*(math.cos(PX_lsq[at][qq+1]) + math.sin(PX_lsq[at][qq+1])*1j))/dgains[gel]
					Jp[at,qt,0,1] = (PY_lsq[at][qq]*(math.cos(PY_lsq[at][qq+1]) + math.sin(PY_lsq[at][qq+1])*1j))/dgains[gel]
					Jp[at,qt,1,0] = (QX_lsq[at][qq]*(math.cos(QX_lsq[at][qq+1]) + math.sin(QX_lsq[at][qq+1])*1j))/dgains[gel]
					Jp[at,qt,1,1] = (QY_lsq[at][qq]*(math.cos(QY_lsq[at][qq+1]) + math.sin(QY_lsq[at][qq+1])*1j))/dgains[gel]
					qq += 2


			for at in range(0,128):
				for qt in range(0,qtt):
					JpH[at,qt] = np.conjugate(np.transpose(Jp[at,qt]))

			for at in range(0,128):
				for qt in range(0,qtt):
					JJ[at,qt] = Jp[at,qt]*JpH[at,qt]

			for at in range(0,128):
				JPX.append([]); JPX[-1]=zeros(N_ch*2)
				JQY.append([]); JQY[-1]=zeros(N_ch*2)
				tq = 0
				for qt in range(0,qtt):

						#print qt, gel, bins[gel], dgains[gel]
					JPX[-1][tq] = sqrt(JJ[at,qt,0,0].real**2 + JJ[at,qt,0,0].imag**2)
					JPX[-1][tq+1] = math.atan2(JJ[at,qt,0,0].imag,JJ[at,qt,0,0].real)
					JQY[-1][tq] = sqrt(JJ[at,qt,1,1].real**2 + JJ[at,qt,1,1].imag**2)
					JQY[-1][tq+1] = math.atan2(JJ[at,qt,1,1].imag,JJ[at,qt,0,0].real)
					tq += 2

		# ===========================

		print "plotting Jones amps for %d frequency channels and %d antennas" % (N_ch, N_ant)
		print "Printing argo"
		print argo
		self.JPX = [JPX[ii][chan_sel] for ii in argo]
		self.JQY = [JQY[ii][chan_sel] for ii in argo]
		self.chan_sel = chan_sel
		print "Printing JPX[2][chan_sel]:"
		print JPX[2][chan_sel]
		# --------------------------------------------------------------------------------- #

		if self.plot_chan==1: freq=ch

		band_start = freq[0]
		quart_bw = ( freq[-1] - band_start ) / 4.0

		rc('font', size=8)

		# --------------------------------------------------------------------------------- #
		# Additional Infos file

		outtile = ('Add_Info_%s.txt' %self.obs)
		out_tile = open(outtile,"w+")
		out_tile.write('#Flagged tiles before the RTS run:\n')
		for ii in flg_t[:-1]:
			if ii != '': out_tile.write('%s\n' %tiles[int(ii)]) #out_tile.write('%s\n' %tiles[argo[int(ii)]])
		out_tile.write('\n')
		#
		# =================================================
		# plotting stdev of Amps - Ants grouped by rx
		# if sel_offset == 1:
		# 	fig = matplotlib.pyplot.figure(figsize=(22.0, 14.0))
		#
		# 	fake=range(1,25)
		#
		# 	ant = 0
		# 	cc = 0
		# 	sp = 0
		# 	sspp = 0
		# 	ppll = 8
		# 	maxo = 0.2
		# 	ccl = ['b','g','r','c','m','y','k','0.75']
		# 	for ii in argo:
		# 		strt = 0
		# 		if ppll ==8:
		# 			ax = fig.add_subplot(4,4,sspp+1);
		# 			if ppll == 8:
		#         			ppll = 0
		#         			sspp +=1
		# 		ppll += 1
		# 		if str(ii) in flg_t: # plot flagged boxes
		# 			ax.plot(fake,'-k')
		# 		else: # plot data
		# 			for l in range(0,len(bins)):
		# 				ax.plot(l+1,np.std(JPX[ii][chan_sel][strt:strt+bins[l]/2]), 's', markeredgecolor=ccl[ppll-1], markersize=6,markerfacecolor='none')
		# 				ax.plot(l+1,np.std(JQY[ii][chan_sel][strt:strt+bins[l]/2]), 'o', c=ccl[ppll-1], markersize=5,alpha=.5)
		# 				if np.std(JPX[ii][chan_sel][strt:strt+bins[l]/2]) > maxo: maxo = np.std(JPX[ii][chan_sel][strt:strt+bins[l]/2])+0.02
		# 				if np.std(JQY[ii][chan_sel][strt:strt+bins[l]/2]) > maxo: maxo = np.std(JQY[ii][chan_sel][strt:strt+bins[l]/2])+0.02
		# 				strt += bins[l]/2
		# 		fig.subplots_adjust(hspace=.2)
		# 		fig.subplots_adjust(wspace=.2)
		# 		plt.title('Rx %i'%sspp)
		# 		sp += 1
		# 		plt.ylim([-0.01,maxo]);plt.xlim([-1,24])
		# 		plt.suptitle('Stdev | color/antenna order = b-g-r-c-m-y-k-grey | Squares = X pol',fontsize=14)
		# 	savefig('StDev_plots_%s.png'%self.obs, bbox_inches='tight', pad_inches=1) #,bbox_inches='tight'
		# 	plt.close(fig)
		# 	#
			# ================================================
			# Check for spikes
			#
			# out_tile.write('#Check for outliers \n')
			# for i in argo:
			# 	m = 0
			# 	strt = 0
			# 	for l in range(0,len(bins)):
			# 		if np.std(JPX[i][chan_sel][strt:strt+bins[l]/2]) > 0.2: out_tile.write('%s - Chan %i - Pol XX\n' %(tiles[i],l+1))
			# 		strt += bins[l]/2
			# 		m += 1
			# 	strt = 0
			# 	m = 0
			# 	for l in range(0,len(bins)):
			# 		if np.std(JQY[i][chan_sel][strt:strt+bins[l]/2]) > 0.2: out_tile.write('%s - Chan %i - Pol YY\n' %(tiles[i],l+1))
			# 		strt += bins[l]/2
			# 		m += 1
		#
		# out_tile.close()
		# -----------------------
		# Plotting

		fig = plt.figure(figsize=(18.0, 10.0))

		if sel_offset == 1: thresh = 0.5

		ant = 0
		cc = 0
		sp = 0
		ppl = 16
		refa = 2 # reference antenna
		# set upper ylim if plotting amps
		maxv = 1.5
		if sel_offset == 1:
			for ii in argo:
				if max(JPX[ii][chan_sel]) > maxv or max(JQY[ii][chan_sel]) > maxv: maxv = max(max(JPX[ii][chan_sel]),max(JQY[ii][chan_sel])) + 0.1
		#
		colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']
		for ii in argo:
			anto = ant + cc
			if str(ii) in flg_t: # plot flagged boxes
				ax = fig.add_subplot(8,16,sp+1,alpha=0.5,axisbg='gray')
				ax.patch.set_alpha(0.5)
				plt.title('Flagged %s'% tiles[ii])
			else: # plot data
				ax = fig.add_subplot(8,16,sp+1,)
				if sel_offset == 1:
					if tra == 1:
						print "Have reached where I think the plotting is happening:"
						print "chan_sel is:"
						print chan_sel
						print "ii is" + str(ii)
						print "frequency(freq_idx) is"
						print freq[freq_idx]
						print JPX[ii][chan_sel]
						print "OK, printing just JPX[ii]"
						print JPX[ii]
						print "Now trying to do it slowly"
						print JPX[ii][0]
						print JPX[ii][1]
						print JPX[ii][2]
						print JPX[ii][3]
						print JPX[ii][4]
						print JPX[ii][5]
						ax.plot(freq[freq_idx],JPX[ii][chan_sel], colours[1])
						ax.plot(freq[freq_idx],JQY[ii][chan_sel], colours[2])
					else:
						ax.plot(freq[freq_idx],PX_lsq[ii][chan_sel], colours[1])
						ax.plot(freq[freq_idx],QY_lsq[ii][chan_sel], colours[2])
					ax.set_xlim(min(freq[freq_idx]), max(freq[freq_idx]))
				plt.title('%s'%tiles[ii])
				ant += 1

			if ppl != 16:
				plt.setp(ax.get_xticklabels(), visible=False) # plot setup
				plt.setp(ax.get_yticklabels(), visible=False)
			if ppl == 16:
				ppl = 0
				plt.setp(ax.get_xticklabels(), visible=False)
			ppl += 1
			if sel_offset == 1: plt.ylim([-0.1,maxv])

			if sp == 15:
				if sel_offset == 1:
					XX_amp, = ax.plot([],[], colours[1],label='XX',linewidth=3.0)
					YY_amp, = ax.plot([],[], colours[2],label='YY',linewidth=3.0)
					ax.legend((XX_amp,YY_amp),('XX','YY'), bbox_to_anchor=(0, 2, .12, .12),prop={'size':14})
				else:
					lgnd = ax.legend(bbox_to_anchor=(0, 2, .12, .12),handler_map={XX_phase: HandlerLine2D(numpoints=4),YY_phase: HandlerLine2D(numpoints=4)},
						prop={'size':14},markerscale=5)
			sp += 1
		plt.tight_layout()
		fig.subplots_adjust(top=0.9)

		if sel_offset == 1:
			plt.suptitle('Amps | Calibrator %s | JD %s' %(calid,JD),fontsize=18)
			plt.show()
			savefig('Amps_%s.png'%self.obs, bbox_inches='tight')

		plt.close(fig)

	def make_the_plot(self, cal_loc):
		if cal_loc[-1] == '/': pass
		else: cal_loc = cal_loc + '/'

		bpp_file = glob.glob("%sBandpassCalibration_*" %cal_loc)
		dii_file = glob.glob("%sDI_JonesMatrices_*" %cal_loc)

		if(1==2):
			flg_t = []
		else:
			try:
				text=open('%sflagged_tiles.txt' %cal_loc,"r").read()
				flg_t=text.split('\n')
			except:
				flg_t = []
				print "flagged_tiles.txt not found, reading flags from *instr_config* file and the logs"

		luu=len(bpp_file)

		# --- grep ChanGains and JD from metafits file
		try:
			hdulist = pyfits.open(cal_loc+self.obs+'.metafits')
		except:
			hdulist = pyfits.open(cal_loc+self.obs+'_metafits_ppds.fits')
		gains = hdulist[1].data['Gains']
		JD = hdulist[0].header['MJD']
		#if(options.apply_digital_gains):
		if (1==2):
			dg = list(gains[0,:]/64.0)
		else:
			#dg = np.ones(len(gains[0,:]))
			dg = [1] * len(gains[0,:])

		##Look for the master log, either with a speific search path or not
		#if search:
			#logs001 = glob.glob("%s*%s*_master.log" %(cal_loc,search))
		#else:
		logs001 = glob.glob("%s*_master.log" %cal_loc)
		logs001.sort(key=lambda x: os.path.getmtime(x))
		logg=logs001[-1]
		text=open(logg,"r").readlines()
		calid = 'NA'
		dgains = []
		for lina in text:
			if 'Freq' in lina:
				friq = float(lina[lina.find("Freq")+5:lina.find("\n")])/10**6
				basech,scrt=divmod((friq+0.005)/1.28+0.5, 1)
				if basech > 128:
					bp_file = sorted(bpp_file, reverse=True)
					di_file = sorted(dii_file, reverse=True)
					#dgains = dg[::-1]
					dgains = dg
				elif basech < 105:
					bp_file = sorted(bpp_file)
					di_file = sorted(dii_file)
					#dgains = dg
					dgains = dg[::-1]
				else:
					diff = int(128-basech)
					dritto = sorted(bpp_file)
					rove = sorted(bpp_file, reverse=True)
					bp_file = dritto[0:(diff+1)]+rove[0:(23-diff)]
					dritto2 = sorted(dii_file)
					rove2 = sorted(dii_file, reverse=True)
					di_file = dritto2[0:(diff+1)]+rove2[0:(23-diff)]
					dgrev = dg[::-1]
					dgains = dgrev[0:23-diff]+dg[0:diff+1]
					#dgains = dg[0:diff+1]+dgrev[0:23-diff]
			elif 'Calibrator   1' in lina:
				if '<' in lina:
					calid = lina[lina.find("<")+1:lina.find(">")]

		#------ loop on files
		lines = []
		counter = 1
		ggo = []
		bins = []
		for bbb in range(0,int(luu)):
			try:
				fid = open(bp_file[bbb], 'r')
			except IOError:
				print 'Can\'t open file \'' + bp_file[bbb] + '\' for reading.'
				sys.exit(0)
			line = fid.readlines()
			lott=len(line)
			tempo = string.split(line[1],',')
			if counter == 1:
				bins.append(len(tempo)-1)
				lines = line
				for el in range(0,lott):
					lines[el] = lines[el].strip()
			elif counter > 1: #+1.04
				bins.append(len(tempo)-1)
				line[0] = line[0].strip()
				a = string.split(line[0],',')
				llu = len(a)
				for ell in range(0,llu):
					a[ell] = str(float(a[ell]) + 1.04*(counter-1))
					b = ','.join(a)
				lines[0] = lines[0] +','+ b
				for el in range(1,lott):
					linea = line[el].split(',')
					lineaa = ','.join(linea[1:])
					lines[el] = lines[el] +','+ lineaa.strip()
			if len(line)==0:  # check that something was read in
				print "Error reading bandpass file: no lines read from the file."
				sys.exit(0)
			elif len(line)%8!=1:
				print "Error reading bandpass file: number of lines should be 1 plus a multiple of 8."
				sys.exit(0)
			fid.close()  # close the file
			counter += 1
		# loop on DI-files (same number of files of BP files)
			ggo.append(self.DI_factors(di_file[bbb]))

		#-------------------------------------------------------------
		#read in the instr_config
		tiles=[]
		try:
			##Find the first thing that resembles a config file
			config_files = glob.glob('%s*instr_config*.txt' %cal_loc)
			if len(config_files) == 0:
				print os.getcwd()
				config_files = glob.glob('%s/*instr_config*.txt' %os.getcwd())
			config_file = config_files[0]

			print 'Using %s as config file' %config_file
			text=open(config_file,"r").readlines()
			yes = 1
			for line in text:
				if line.startswith('#'): continue
				if yes == 2:
					entr = line.split('\t')
					tiles.append(entr[-1].rstrip('\n').split(' ')[-1])
					flagg = entr[-1].rstrip('\n').split(' ')[0]
					if int(flagg) == 1:
						if str(int(entr[0])/2) not in flg_t: flg_t.append(str(int(entr[0])/2))
					yes=0
				yes += 1
			print flg_t, tiles
		except:
			print 'Instr_config file not found - reading Tile flags from metafits file'
			flags = hdulist[1].data['FLAG'][::2]
			tilesn = hdulist[1].data['Tile'][::2]
			for i in range(0,len(tilesn)):
				tiles.append('Tile'+str(tilesn[i]))
			for l in range(0,len(flags)):
				if flags[l] == 1: flg_t.append(str(l))

		##Read in a node log file and check for flagged tiles - may be more flags due to Bart's
		##super flagging

		#if search:
			#node_log = glob.glob('%s*%s*node*.log' %(cal_loc,search))[0]
		#else:
		node_log = glob.glob('%s*node*.log' %cal_loc)[0]

		node_lines = open(node_log,'r').read().split('\n')
		for line in node_lines:
			if 'TileFlags' in line:
				if len(line.split()) == 6:
					flgs = line.split()[2]
					for flg in flgs.split(','):
	                                       # if(options.ignore_flagged_tiles_text == False):
										    if flg not in flg_t: flg_t.append(flg)

		supr = []
		for e in tiles:
			subb = e.strip('Tile')
			supr.append(int(subb))
		argo = sorted(range(len(supr)),key=lambda x:supr[x])

		#----------------------------------

		data = [argo,lines,bins,dgains,ggo,flg_t,tiles,calid,JD]

		##Create a plot for either/both phase or amplitude
		if self.sel_offset == 0 or self.sel_offset == 1: self.phase_or_amp(1,data)
		#if sel_offset == 0 or sel_offset == 2: phase_or_amp(2,data)


	def plotAmplitudes(self, obsdir):
		print "Plotting amplitudes"
		self.make_the_plot(obsdir)

# parser = optparse.OptionParser()
#
# parser.add_option('-b','--base_dir',
# 	help='MANDATORY. Either where the cal solutions live, or where the observation dirs live')
# #parser.add_option('-f','--file', help='Optional. Single bandpass file to plot.')
# parser.add_option('-p','--phases', default=False,action='store_true',help='Optional. Just plot the phases')
# parser.add_option('-a','--amps', default=False,action='store_true',help='Optional. Just plot the amplitudes')
# parser.add_option('-o','--obs', default=False,help='Optional. Text file containing a list of the obs IDs located in the base_dir')
# parser.add_option('-s','--search', default=False,
# 	help='Optional. Extra string to search for particular rts logs, i.e. putting --search=string will search for rts*string*log rather than rts*log - useful if you have multiplelogs hanging around')
# parser.add_option('-g','--apply_digital_gains', default=False,action='store_true',help='Apply digital gains read from metafits. Default is to assume they have already been applied by the RTS')
# parser.add_option('-i','--ignore_flagged_tiles_text', default=False,action='store_true',help='Ignore the flagged_tiles.txt files in situations where is does not match the flagging of the calibration solutions')
# parser.add_option('-n','--name', default=False,help='Optional. Assign obsid on command line. Useful for non-standard directory structures')
# options,args = parser.parse_args()

# sel_offset = 0
# plot_chan = 0
# plot_raw  = 0
# obs = 'NoObsID'
#----------------------------------------------------------------------------------------------------------------------#
# read command-line arguments

# base_dir = options.base_dir
# if base_dir[-1] == '/':
# 	pass
# elif base_dir[-1] != '/' :
# 	if str(options.obs) == 'False':
# 		if(options.name):
# 			obs = options.name
# 		else:
# 			obs = base_dir.split('/')[-1]
# 	else:
# 		base_dir = base_dir + '/'

#if options.file: bp_file = options.file
# if options.phases: sel_offset = 2
# if options.amps: sel_offset = 1
# if options.obs: obs = options.obs
# search = options.search

#cwd = os.getcwd()


##If more than one obs ID to process, iterate through the directories, otherwise
##just plot a single plot
# if options.obs:
# 	print '----------------------------------------------------------------------------'
# 	obsIDs = open(obs,'r').read().split('\n')
# 	del obsIDs[-1]
# 	for obsID in obsIDs:
# 		obs = obsID
# 		print "Plotting %s" %obs
# 		make_the_plot(base_dir+obsID)
# 		print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
# 	print '------------------------------- Finito -------------------------------'
# else:
# 	make_the_plot(base_dir)
