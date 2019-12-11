import ROOT as ro
import numpy as np
from optparse import OptionParser
import os, ipdb, time
import progressbar

fit_func_options = ['erf', 'tanh', 'weibull']
fit_functions = {'erf': '[3]*(TMath::Erf((x-[0])/[1])+[2])',
                 'tanh': '[3]*(TMath::TanH((x-[0])/[1])+[2])',
                 'weibull': '[3]*([2]-exp(-(x/[0])**[1]))'}
                 # 'weibull': '[3]*([2]-exp(-(x**[1])/[0]))'}

fit_parms = {'erf': np.array((250, 320, 0.75, 125), 'f8'),
             'tanh': np.array((250, 320, 0.75, 125), 'f8'),
             'weibull': np.array((250, 1.1, 10, 25), 'f8')}
             # 'weibull': np.array((250 ** 1.5, 1.5, 10, 25), 'f8')}

fit_params_low = {'erf': {0: 0, 1: 1, 2: -10, 3: 10},
                  'tanh': {0: 0, 1: 1, 2: -10, 3: 10},
                  'weibull': {0: 0, 1: 1.001, 2: 0, 3: 1}}

fit_params_high = {'erf': {0: 10000, 1: 10000, 2: 10, 3: 10000},
                  'tanh': {0: 10000, 1: 10000, 2: 10, 3: 10000},
                  'weibull': {0: 10000, 1: 2, 2: 128, 3: 1000}}
                  # 'weibull': {0: 1000000, 1: 2, 2: 1000, 3: 1000}}

fit_method = ('Minuit2', 'Migrad',)  # found to be the best for this case
# fit_method = ('Minuit2', 'Scan',)
# fit_method = ('Minuit2', 'Simplex',)
# fit_method = ('Fumili2',)
# fit_method = ('GSLMultiFit', 'GSLLM',)
# fit_method = ('GSLMultiFit', 'BFGS2')
# fit_method = ('GSLSimAn', 'SimAn',)


def Correct_Path(path, times=2):
	"""
	This method is used to get the correct path to the file regardless to how it is passed
	:param path: path to file. can be absolute or relative
	:param times: number of times to iterate the method. sometimes 2 iterations are required
	:return: returns the full absolute path of the one given
	"""
	abs_path = ''
	if path[0] == '~':
		abs_path += os.path.expanduser('~')
		abs_path += path[1:]
	elif os.path.isabs(path):
		abs_path += path
	else:
		abs_path += os.path.abspath(path)
	if times != 1:
		return Correct_Path(abs_path, 1)
	return abs_path


def RoundInt(n, nptype='int32'):
	"""
	Method used to round to the nearest integer. It can be used on arrays or single numbers
	:param n: The number to be rounded
	:param nptype: The type used before returning. It is used for truncation. It follows numpy dtype
	:return: a number or a list of the rounded number in the defined type (int or float)
	"""
	val = np.floor(np.add(n, 0.5, dtype='f8'), dtype='f8').astype(nptype)
	if 'i' in nptype.lower():
		if val.size == 1:
			return int(val)
		return list(val)
	elif 'f' in nptype.lower():
		if val.size == 1:
			return float(val)
		return list(val)
	return val


def CreateProgressBarUtils(maxVal=1):
	widgets = [
		'Processed: ', progressbar.Counter(),
		' out of {mv} '.format(mv=maxVal), progressbar.Percentage(),
		' ', progressbar.Bar(marker='>'),
		' ', progressbar.Timer(),
		' ', progressbar.ETA()
		# ' ', progressbar.AdaptativeETA(),
		#  ' ', progressbar.AdaptativeTransferSpeed()
	]
	bar = progressbar.ProgressBar(widgets=widgets, maxval=maxVal)
	return bar


class phCalibrationFits:
	def __init__(self, ph_file, i2c, trim_val=35, ntriggs=20, after2019=False, fit_func='erf'):
		self.ph_file = Correct_Path(ph_file)
		self.i2c = i2c
		self.n_triggs = RoundInt(ntriggs / 2.0)
		self.trim_val = trim_val
		self.ph_file_dir = '/'.join(self.ph_file.split('/')[:-1])
		self.ph_file_name = self.ph_file.split('/')[-1]
		self.after2019 = after2019
		self.fit_func = fit_func if fit_func in fit_func_options else 'erf'
		self.vcals_low = np.zeros(1, 'uint16')
		self.vcals_high = np.zeros(1, 'uint16')
		self.vcals_all = np.zeros(1, 'uint16')
		self.high_scale = 7
		self.is_data_stacked = False
		self.pix_data = {}  # col, row dic of pix data
		self.pix_graph = {}
		self.pix_fit = {}
		self.pix_fit_probab = {}
		self.pix_fit_ndf = {}
		self.pix_masked = {}
		ro.Math.MinimizerOptions.SetDefaultMinimizer(*fit_method)
		ro.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(1000000)
		ro.Math.MinimizerOptions.SetDefaultTolerance(0.00001)
		ro.gStyle.SetOptFit(1111)
		self.ReadPhFile()
		self.SaveFitsFile()
		self.PrintBadPixelsCalibration()
		self.canvas = {}

	def ReadPhFile(self):
		"""
		This method reads and fits the phCalibration file given to the class. It uses the information on number of triggers used, function to fit, and trimming value.
		:return: Nada :D
		"""
		if os.path.isfile(self.ph_file):
			with open(self.ph_file) as ff0:
				lines = ff0.readlines()
				bar = CreateProgressBarUtils(len(lines))
				bar.start()
				for it, line in enumerate(lines):
					if 'Low range' in line:
						vcalsl = line.replace('Low range:  ', '').replace(' \n', '').split()
						self.vcals_low = np.array([int(v) for v in vcalsl], 'uint16')
						self.is_data_stacked = self.vcals_low.size > 100 and not self.after2019
					if 'High range' in line:
						vcalsh = line.replace('High range:  ', '').replace(' \n', '').split()
						self.vcals_high = np.array([int(v) * self.high_scale for v in vcalsh], 'uint16')
					if 'Pix' in line:
						templine = line.replace('\n', '')
						(col, row) = [int(rc) for rc in templine.split('Pix')[-1].split()]
						if col not in self.pix_data.keys():
							self.pix_data[col] = {}
						self.pix_data[col][row] = np.array([int(v) for v in templine.split('Pix')[0].split()], 'uint16')
						self.EstimateLowHighIfStacked(col, row)
						self.SetGraph(col, row)
						self.FitPixel(col, row)
					bar.update(it + 1)
				bar.finish()

	def EstimateLowHighIfStacked(self, col, row):
		"""
		If the data is stacked, that is before the December 2019 fix after desy test beam and there are more than 100 vcals in the lower range, then this method estimates the mixed information
		The corrected data is updated in the dictionary pix_data[col][row]
		:param col: column of the pixel
		:param row: row of the pixel
		:return: Nada :D
		"""
		if self.is_data_stacked:
			estimated_missing_low = np.zeros(self.vcals_high.size)
			if self.vcals_low.size > 100 + self.vcals_high.size:
				# there are still vcals for the low range after the mix with the stacking with the high range. approximate with a line between known points
				estimated_missing_low = RoundInt([(self.pix_data[col][row][100 + self.vcals_high.size].astype('int32') - self.pix_data[col][row][99].astype('int32')) * (vcal - self.vcals_low[99].astype('int32')) / float(self.vcals_low[100 + self.vcals_high.size].astype('int32') - self.vcals_low[99].astype('int32')) + self.pix_data[col][row][99].astype('int32') for vcal in self.vcals_low[100: 100 + self.vcals_high.size]])
			elif self.vcals_low.size > 100:
				# there are no vcals for the low range after the mix with the stacking. assume last known value
				estimated_missing_low = [self.pix_data[col][row][99] for vcal in self.vcals_low[100: 100 + self.vcals_high.size]]
			estimated_compensated_high = np.array(RoundInt(np.subtract(self.pix_data[col][row][-self.vcals_high.size:], estimated_missing_low))).astype('uint16')
			estimated_compensated_high = np.where(estimated_compensated_high > 255, 255, estimated_compensated_high)  # the roc does not return more than 255
			for it in xrange(self.vcals_high.size):
				self.pix_data[col][row][100 + it] = estimated_missing_low[it]
				self.pix_data[col][row][it - self.vcals_high.size] = estimated_compensated_high[it]

	def SetGraph(self, col, row):
		"""
		Sets up the graph for fitting the raw data. The object is stored in a dictionary pix_graph[col][row]
		:param col: column of the pixel
		:param row: row of the pixel
		:return: Nada :D
		"""
		if col not in self.pix_graph.keys():
			self.pix_graph[col] = {}
		if self.vcals_all.size == 1:
			self.vcals_all = np.concatenate([self.vcals_low, self.vcals_high]).astype('f8')
		yerrsl = np.sqrt(np.divide(self.pix_data[col][row], self.n_triggs, dtype='f8'), dtype='f8')
		# yerrsh = np.sqrt(np.divide(self.pix_data[col][row], self.n_triggs, dtype='f8'), dtype='f8')
		yerrsh = np.where(self.pix_data[col][row] < 255, np.sqrt(np.divide(self.pix_data[col][row], self.n_triggs, dtype='f8'), dtype='f8'), 200).astype('f8')
		self.pix_graph[col][row] = ro.TGraphAsymmErrors(self.vcals_all.size, self.vcals_all, self.pix_data[col][row].astype('f8'), np.zeros(self.vcals_all.size, 'f8'), np.zeros(self.vcals_all.size, 'f8'), yerrsl, yerrsh)
		self.pix_graph[col][row].GetYaxis().SetRangeUser(0, 255)
		self.pix_graph[col][row].GetXaxis().SetTitle('Vcal [ADC]')
		self.pix_graph[col][row].GetYaxis().SetTitle('PH [ADC]')
		nameTitle = 'Pix_{c}_{r}_{f}_fit'.format(c=col, r=row, f=self.fit_func.capitalize())
		self.pix_graph[col][row].SetNameTitle(nameTitle, nameTitle)
		self.pix_graph[col][row].SetMarkerColor(ro.kBlue)
		self.pix_graph[col][row].SetLineColor(ro.kBlue)
		self.pix_graph[col][row].SetMarkerStyle(ro.TAttMarker.kFullDotMedium)

	def FitPixel(self, col, row):
		"""
		This method does the fitting for a given pixel. the fitted data is stored in the dictionary pix_fit_probab[col][row]. Dictionary pix_masked[col][row] is initialized
		:param col: column of the pixel
		:param row: row of the pixel
		:return: Nada :D
		"""
		if col not in self.pix_fit.keys():
			self.pix_fit[col] = {}
		self.pix_fit[col][row] = ro.TF1('{f}_{c}_{r}'.format(f=self.fit_func, c=col, r=row), fit_functions[self.fit_func], 0, 255 * self.high_scale)
		self.pix_fit[col][row].SetNpx(int(255 * self.high_scale))
		self.pix_fit[col][row].SetParameters(fit_parms[self.fit_func])
		for param in xrange(len(fit_parms[self.fit_func])):
			self.pix_fit[col][row].SetParLimits(param, fit_params_low[self.fit_func][param], fit_params_high[self.fit_func][param])
		xfitmax = self.vcals_all[np.abs(self.pix_data[col][row] - 255).argmin()] if (self.pix_data[col][row] > 254).any() else self.vcals_all.max()
		self.pix_graph[col][row].Fit('{f}_{c}_{r}'.format(f=self.fit_func, c=col, r=row), 'Q0', '', self.trim_val, xfitmax)
		if self.pix_fit[col][row].GetProb() < 0.9:
			self.pix_graph[col][row].Fit('{f}_{c}_{r}'.format(f=self.fit_func, c=col, r=row), 'Q0', '', self.trim_val, xfitmax)
		if self.pix_fit[col][row].GetProb() < 0.9:
			self.pix_graph[col][row].Fit('{f}_{c}_{r}'.format(f=self.fit_func, c=col, r=row), 'Q0', '', self.trim_val, xfitmax)
		if col not in self.pix_fit_probab.keys():
			self.pix_fit_probab[col] = {}
		self.pix_fit_probab[col][row] = self.pix_fit[col][row].GetProb()
		if col not in self.pix_fit_ndf.keys():
			self.pix_fit_ndf[col] = {}
		self.pix_fit_ndf[col][row] = self.pix_fit[col][row].GetNDF()
		if col not in self.pix_masked.keys():
			self.pix_masked[col] = {}
		self.pix_masked[col][row] = False

	def SaveFitsFile(self):
		"""
		This method saves the new fitted parameters stored in pix_fit dictionary. If the file already exists, a time stamp is appended to avoid overwriting. The user is then free to decide which file to keep
		:return: Nada :D
		"""
		fit_file_name = self.ph_file_name.split('Calibration')[0] + 'CalibrationFit' + self.fit_func.capitalize() + self.ph_file_name.split('Calibration')[1]
		if os.path.isfile('{d}/{f}'.format(d=self.ph_file_dir, f=fit_file_name)):
			date = time.localtime()
			fit_file_name = fit_file_name.replace('.dat', '')
			fit_file_name += '_{y}{m}{d}_{h}_{min}'.format(y=date.tm_year, m=date.tm_mon, d=date.tm_mday, h=date.tm_hour, min=date.tm_min)
			fit_file_name += '.dat'
		with open('{d}/{f}'.format(d=self.ph_file_dir, f=fit_file_name), 'w') as f1:
			f1.write('Parameters of the vcal vs. pulse height fits\n')
			f1.write('{f}\n\n'.format(f=fit_functions[self.fit_func]))
			for col, rowfit in self.pix_fit.iteritems():
				for row, fit in rowfit.iteritems():
					self.pix_fit_probab[col][row] = fit.GetProb()
					self.pix_fit_ndf[col][row] = fit.GetNDF()
					for pari in range(fit.GetNpar()):
						f1.write('{p} '.format(p=fit.GetParameter(pari)))
					f1.write('    Pix {c} {r}\n'.format(c=col, r=row))

	def PrintBadPixelsCalibration(self, th=0.9, minNDF=9):
		"""
		This method is used to identify bad pixels using how good the fitting was. The pixels are printed and a description of why they are flagged is shown. The inefficiency is shown to guide the user in the coarse selection of the parameters
		A more detailed study using PlotPixelFit is suggested to define the final parameters
		:param th: the threshold used by the fitter to estimate how good was the fit taking into accound the chi2 and the ndf. It is highly suceptible to the errors (n_triggers). It should be changed accordingly
		:param minNDF: minimum number of degrees of freedon for the fit. It is up to the user to define what this minimum should be. Usually less than 8 is bad for any of the fitted functions as they have 4 parameters
		:return: Nada :D
		"""
		ineff = 0
		tot = 0
		for col, rowpf in self.pix_fit_probab.iteritems():
			for row, pf in rowpf.iteritems():
				tot += 1
				self.pix_masked[col][row] = False
				if self.pix_fit_ndf[col][row] < minNDF:
					self.pix_masked[col][row] = True
					ineff += 1
					print 'pixel: ', col, row, 'has a NDF of', self.pix_fit_ndf[col][row], '. Fit Not reliable. Its probab is', pf
				elif pf < th:
					self.pix_masked[col][row] = True
					ineff += 1
					print 'pixel: ', col, row, 'has a fit probab of', pf, 'and has a NDF of', self.pix_fit_ndf[col][row]
		print 'Total calibration inefficiency: {v}% ({v2} pixels)'.format(v=ineff * 100. / tot, v2=ineff)

	def PlotPixelFit(self, col, row):
		"""
		This method is used to plot the raw data with the fit and show information of the fit. It is used to help decide whether the particular pixel should survive the masking or not
		This method should be used in tandem with PrintBadPixelsCalibration to make correct decisions for the mask file
		:param col: Column of the pixel
		:param row: Row of the pixel
		:return: Nada :D
		"""
		if col not in self.canvas.keys():
			self.canvas[col] = {}
		if row not in self.canvas[col].keys():
			self.canvas[col][row] = ro.TCanvas('c_pix_{c}_{r}'.format(c=col, r=row), 'c_pix_{c}_{r}'.format(c=col, r=row), 1)
		self.canvas[col][row].cd()
		self.canvas[col][row].SetGridx()
		self.canvas[col][row].SetGridy()
		self.canvas[col][row].SetTicky()
		ro.gPad.Update()
		self.pix_graph[col][row].Draw('ap')
		self.pix_fit[col][row].Draw('same')
		ro.gPad.Update()
		if self.pix_graph[col][row].FindObject('stats'):
			self.pix_graph[col][row].FindObject('stats').SetX1NDC(0.54)
			self.pix_graph[col][row].FindObject('stats').SetX2NDC(0.9)
			self.pix_graph[col][row].FindObject('stats').SetY1NDC(0.1)
			self.pix_graph[col][row].FindObject('stats').SetY2NDC(0.38)
			ro.gPad.Update()
		self.canvas[col][row].BuildLegend(0.54, 0.38, 0.9, 0.48)
		ro.gPad.Update()

	def SaveMaskFile(self, th=0.9, minNDF=8, maskBorders=True):
		"""
		This method creates a mask file excluding the pixels that are cutted by the parameters "th" and "minNDF". If the file already exists, a time stamp is appended to let the user decide which file to keep
		:param th: minimum allowed probab of the fitted parameters given the data.
		:param minNDF: minimum allowed number of degrees of freedom used for the fitting
		:param maskBorders: boolean variable. If true, the mask file will also contain the borders
		:return: Nada :D
		"""
		self.PrintBadPixelsCalibration(th, minNDF)
		mask_file_name = 'Mask_' + self.ph_file_name.split('Calibration')[0] + 'CalibrationFit' + self.fit_func.capitalize() + self.ph_file_name.split('Calibration')[1].replace('.dat', '') + '_th_' + str(th) + '_minDNF_' + str(minNDF) + '.dat'
		if os.path.isfile('{d}/{f}'.format(d=self.ph_file_dir, f=mask_file_name)):
			date = time.localtime()
			mask_file_name = mask_file_name.replace('.dat', '')
			mask_file_name += '_{y}{m}{d}_{h}_{min}'.format(y=date.tm_year, m=date.tm_mon, d=date.tm_mday, h=date.tm_hour, min=date.tm_min)
			mask_file_name += '.dat'
		with open('{d}/{f}'.format(d=self.ph_file_dir, f=mask_file_name), 'w') as f1:
			f1.write('# mask file generated by low fit probability and not enough points for a reliable fit\n')
			for col, rowmasked in self.pix_masked.iteritems():
				for row, ismasked in rowmasked.iteritems():
					if ismasked:
						f1.write('pix {i} {c} {r}\n'.format(i=i2c, c=col, r=row))
					elif maskBorders and (col == 0 or col == 51 or row == 0 or row == 79):
						f1.write('pix {i} {c} {r}\n'.format(i=i2c, c=col, r=row))
		print 'Mask file', mask_file_name, 'saved inside directory', self.ph_file_dir

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-f', '--file', dest='file', type='str', help='file containing the vcal adc for each pixel')
	parser.add_option('-i', '--i2c', dest='i2c', type='int', help='i2c of the device. Used for saving the mask file if desired')
	parser.add_option('-n', '--ntriggs', dest='ntriggs', type='int', default=20, help='number of triggers sent per vcal. Used to estimate the uncertainties on the raw data')
	parser.add_option('-t', '--trimmval', dest='trimmval', type='int', default=35, help='trimm value used. Used to give priority to the fit for vcals starting from trimmval')
	parser.add_option('--after2019', dest='after2019', default=False, action='store_true', help='If the data was taken after 2019, then no overlap exists in the vcals histogram. Do not use it if data was taken before December 2019 test beam')
	parser.add_option('--fitfunction', default='erf', type='str', help='fit function to use. The options implemented are: "erf" (default), "tanh", "weibull"')

	(options, args) = parser.parse_args()
	ffile = str(options.file)
	i2c = int(options.i2c)
	after2019 = bool(options.after2019)
	fitFunc = str(options.fitfunction)
	ntriggs = int(options.ntriggs)
	trimmv = int(options.trimmval)

	z = phCalibrationFits(ffile, i2c, trimmv, ntriggs, after2019, fitFunc)
