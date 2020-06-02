#!usr/bin/env python3

import argparse
import sys
import os
import csv
import pysam
from datetime import datetime
import math
import pandas as pd

extended_help = """
"""

# Input command line arguments
parser = argparse.ArgumentParser(
	description='CNV caller.',
	formatter_class=argparse.RawDescriptionHelpFormatter,
	epilog=extended_help)
parser.add_argument(
	'--samplesAll', 
	required=False, 
	type=str, 
	default= None, 
	metavar='<path>', 
	help='path to a directory of all bam files')

parser.add_argument(
	'--samplesControl', 
	required=False, 
	type=str, 
	default= None,
	metavar='<path>', 
	help='path to a directory of bam control files')
parser.add_argument(
	'--samplesExperimental', 
	required=False, 
	type=str, 
	default= None,
	metavar='<path>', 
	help='path to a directory of bam experimental files')
parser.add_argument(
	'--window', 
	required=False, 
	type=int, 
	default= 1000,
	metavar='<int>', 
	help='size of the window')
parser.add_argument(
	'--sampleInfo',
	required=False,
	type=str,
	metavar='<path>',
	help='path to sample info csv file')
parser.add_argument(
	'--controlName',
	required=False,
	type=str,
	default = 'Control',
	metavar='<str>',
	help = 'name of control class')

arg = parser.parse_args()

# Write print statement outputs to file
sys.stdout = open(datetime.now().strftime('%I:%M%p_%b%d') + 'CNV_Me.print', 'w')

# Chromosomes of interest
chromosomes = []
for i in range (1, 23): # 1..22
        chromosomes.append('chr' + str(i))
chromosomes.extend(["chrX", "chrY"])

def getBamfiles(pathAll, pathControl, pathExp, sampleInfo):
        bamfilePaths = []
	numControl = 0
	numControlM = 0
	numControlF = 0
	numExp = 0
        if pathAll == None:
                print("samplesAll argument is empty. Now using samplesControl and samplesExperimental arguments.")
                if pathControl == None or pathExp == None:
                        print("samplesControl and samplesExperimental arguments cannot be empty if samplesAll argument is not provided.")
                else:
                        for bamfile in os.listdir(pathControl):
                                if bamfile.endswith('.bam'):
                                        bamfilePaths.append(pathControl + '/' + bamfile)
					numControl += 1
                        for bamfile in os.listdir(pathExp):
                                if bamfile.endswith('.bam'):
                                        bamfilePaths.append(pathExp + '/' + bamfile)
					numExp += 1
        else:
		for bamfile in os.listdir(pathAll):
                        if bamfile.endswith('.bam'):
                        	bamfilePaths.append(pathExp + '/' + bamfile)
				samplename = str(bamfile.split("_")[0])
				if sampleInfo.loc[samplename, 'Diagnosis'] == arg.controlName:
					numControl += 1
					if sampleInfo.loc[samplename, 'Sex'] == 'M':
						numControlM += 1
					elif sampleInfo.loc[samplename, 'Sex'] == 'F':
						numControlF += 1 
				else:
					numExp += 1
				
        return (bamfilePaths, numControl, numControlM, numControlF, numExp)


def populateCnvReads(bamfilePaths, sampleInfo):
	"""
	Function: 
	Read BAM files and store coverage for each chrm:bin:sample.
	
	Parameters:
	bamfilePaths (list): list containing all BAM file paths.
	
	Returns:
	cnv (dict): dictionary containing cnv info of each chrm:bin:sample.
	bamfiles (list): list containing BAM files. 
	stats (dict): dictionary containing the total reads and total read lengths for each sample.
	"""
	cnv = {chrm: {} for chrm in chromosomes}
	stats = {}
	print(chromosomes)
	bamfiles = []
	samfileStats = []
	
	for bamfile in bamfilePaths:	
		bamfiles.append(bamfile)
		print('bamfile = %s' % bamfile)
		tstart = datetime.now()
		samfile = pysam.AlignmentFile(bamfile, 'rb') 
					
		stats[bamfile] = {key: 0 for key in [*chromosomes, "totalReads", "totalLen", "sex"]}	
		
		# File name is last component of full path	
		filename = bamfile.split("/")[-1]
		print(filename)
		# Sample name includes characters before 1st underscore
		samplename = str(filename.split("_")[0])

		recordedSex = sampleInfo.loc[samplename, 'Sex']	
	
		stats[bamfile]["sex"] = recordedSex
	
		# stats[bamfile]["sex"] = "male" or "female"
		for read in samfile:
			# Skip PCR and optical duplicate reads
			if read.is_duplicate == False:
				chrm = read.reference_name
				readLength = read.reference_length
				pos = read.reference_start
				if chrm in cnv:
					# Increment stats components
					stats[bamfile][chrm] += 1
					stats[bamfile]['totalReads'] += 1				
					stats[bamfile]['totalLen'] += readLength		
					
					# Set binkey	
					binkey = int(pos/arg.window)
					if binkey in cnv[chrm]:
						if bamfile in cnv[chrm][binkey]:
						# binkey and sample exist
							cnv[chrm][binkey][bamfile]['sumReadCounts'] += 1
							cnv[chrm][binkey][bamfile]['sumReadLengths'] += (1*readLength)
						else:
						# binkey exists but add new sample
							cnv[chrm][binkey].update({bamfile: {\
								'sumReadCounts': 1,\
								'sumReadLengths': 1*readLength,\
								'copynum': 0
							}})
					else: 
						# add new binkey and new sample
						cnv[chrm].update({binkey: {bamfile: {\
							'sumReadCounts': 1,\
							'sumReadLengths': 1*readLength,\
							'copynum': 0
						}}})
				
		print('time for bamfile = %s' % (datetime.now() - tstart))
	return (cnv, bamfiles, stats)

def normalizeCnv(cnv, chrm, binkey, stats, numMaleSamples, numFemaleSamples, numControlM, numControlF):
	"""
	Function:
	Find normalization factors for autosomal and sex chromosomes for each chrm:bin:sample

	Parameters:
	cnv (dict): current cnv dictionary
	chrm (str): chromosome key for cnv dictionary
	binkey (int): binkey key for cnv:chrm dictionary
	numMaleSamples (int): number of male samples
	numFemaleSamples (int): number of female samples
	numControlM (int): number of male control samples
	numControlF (int): number of female control samples

	Returns:
	normfactorAut (int): normalization factor for all autosomal chromosomes
	normfactorMaleX (int): normalization factor for all male X chromosomes
	normfactorMaleY (int): normalization factor for all male Y chromosomes
	normfactorFemaleX (int): normalization factor for all female X chromosomes
	normfactorFemaleY (int): normalization factor for all female Y chromosomes
	"""
	normfactorAut = 0
	normfactorMaleX = 0
	normfactorFemaleX = 0
	normfactorMaleY = 0
	normfactorFemaleY = 0

	for bamfile in cnv[chrm][binkey]:
		if "control" in bamfile: # TODO: change to if bamfile diagnosis == "control"
			if bamfile in cnv[chrm][binkey].keys():
				# Coverage      = sum of (reads * read length) in each bin / bin size
				# To normalize: divide by total coverage of sample = sum of (reads * read length) in sample
				binCoverage = cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.window
				totalReadsAndLengths = stats[bamfile]['totalLen']
				normalizedBinCoverage = binCoverage / totalReadsAndLengths
				
				# if bamfile is from female, stays the same, diploid X, no Y
				# if bamfile is from male, use normfactorSexChrm and normalize to 1x each

				sampleSex = stats[bamfile]["sex"]
				if chrm == "chrX":
					if sampleSex == "M":
						normfactorMaleX += normalizedBinCoverage
					if sampleSex == "F":
						normfactorFemaleX += normalizedBinCoverage
				elif chrm == "chrY":
					if sampleSex == "M":
						normfactorMaleY += normalizedBinCoverage
					if sampleSex == "F":
						normfactorFemaleY += normalizedBinCoverage
				else:
					normfactorAut += normalizedBinCoverage

	# Normalize coverage values to a copy number of ~2 for autosomal chromosomes (for normal diploid)       
	normfactorAut = (normfactorAut) / (numControlM + numControlF)
	# Normalize coverage value to a copy number of ~1 for male sex chromosomes
	if numControlM != 0:
		normfactorMaleX = normfactorMaleX / numControlM
		normfactorMaleY = normfactorMaleY / numControlM
	if numControlF != 0:
		normfactorFemaleX = normfactorFemaleX / numControlF
		normfactorFemaleY = normfactorFemaleY / numControlF

	return (normfactorAut, normfactorMaleX, normfactorMaleY, normfactorFemaleX, normfactorFemaleY)

	
def populateCnvCopynum(cnv, stats, numMaleSamples, numFemaleSamples, numControlM, numControlF):
	"""
	Function: 
	Update "cnv" dictionary with normalized copy number estimations.
	
	Parameters:
	cnv (dict): dictionary containing raw reads and lengths in each chrm:bin:sample.
	numControl (int): number of control BAM files (samples).
	
	Returns:
	cnv (dict): updated dictionary containing copy number estimations in each chrm:bin:sample.
	"""
	for chrm in cnv:
		for binkey in cnv[chrm]:
			normfactorAut, normfactorMaleX, normfactorMaleY, normfactorFemaleX, normfactorFemaleY = \
				normalizeCnv(cnv, chrm, binkey, stats, numMaleSamples, numFemaleSamples, numControlM, numControlF)		
			
			for bamfile in cnv[chrm][binkey]:
				if bamfile in cnv[chrm][binkey].keys():
					sampleSex = stats[bamfile]["sex"]
					if chrm == "chrX":
						if sampleSex == "M" and normfactorMaleX != 0:
							cnv[chrm][binkey][bamfile]['copynum'] = \
								(cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.window) / \
								stats[bamfile]['totalLen'] / normfactorMaleX
						if sampleSex == "F" and normfactorFemaleX != 0:
							cnv[chrm][binkey][bamfile]['copynum'] = \
								(cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.window) / \
								stats[bamfile]['totalLen'] / (0.5 * normfactorFemaleX)
					elif chrm == "chrY":
						if sampleSex == "M" and normfactorMaleY != 0:
							cnv[chrm][binkey][bamfile]['copynum'] = \
								(cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.window) / \
								stats[bamfile]['totalLen'] / normfactorMaleY
						if sampleSex == "F" and normfactorFemaleY != 0:
							cnv[chrm][binkey][bamfile]['copynum'] = \
								(cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.window) / \
								stats[bamfile]['totalLen'] / normfactorFemaleY
					else:
						if normfactorAut != 0:
							cnv[chrm][binkey][bamfile]['copynum'] = \
						 		(cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.window) / \
								stats[bamfile]['totalLen'] / (0.5 * normfactorAut)
	return cnv


def outputCnv(cnv, samples):
	"""
	Function: 
	Output CNV info to file.
	
	Parameters:
	cnv (dict): dictionary containing CNV estimations
	samples (list): list contanining BAM files paths

	Returns:
	None
	"""
	with open(datetime.now().strftime('%I:%M%p_%b%d') + 'CNV_Me_output.txt', 'w', newline='') as outfile:
		outfile = csv.writer(outfile, delimiter='\t')
		samples.sort()
		outfile.writerow(['chr', 'start', 'end', *samples])
        
		for chrm in cnv:
			for binkey in cnv[chrm]:
				readsList = []
				copynumList = []
				for sample in sorted (cnv[chrm][binkey]):
					copynumList.append(cnv[chrm][binkey][sample]['copynum'])
				copynumList = [round(elem, 3) for elem in copynumList]

				if len(cnv[chrm][binkey]) == len(samples):
					outfile.writerow([
						chrm,
						binkey*arg.window,
						(binkey+1)*arg.window,
						*copynumList
					])
	return None


if __name__ == '__main__':
	timestart = datetime.now()	
	print('timestart = %s' % timestart)
	print(arg.controls)
	print(arg.experimentals)
	print(arg.window)
	print(arg.sampleInfo)

	sampleInfo = pd.read_csv(arg.sampleInfo, converters={i: str for i in range(len(pd.read_csv(arg.sampleInfo)))})
	sampleInfo.index = list(sampleInfo['Name'])
	
	# Convert sex label to M and F
	for i in range(0, len(sampleInfo)):
		if ('M' or 'm' or 'Male' or 'male') in sampleInfo['Sex'][i]:
			sampleInfo.loc[i, 'Sex'] = 'M'
		if ('F' or 'f' or 'Female' or 'female') in sampleInfo['Sex'][i]:
			sampleInfo.loc[i, 'Sex'] = 'F'
	
	numMaleSamples = len(sampleInfo[sampleInfo['Sex'] == 'M'])
	numFemaleSamples = len(sampleInfo[sampleInfo['Sex'] == 'F']) 	
	
	bamfilePaths, numControl, numControlM, numControlF, numExp = getBamfiles(arg.samplesAll, arg.samplesControl, arg.samplesExperimental, sampleInfo)
	print(*bamfilePaths, sep='\n')
	print(numControl)
	print(numExp)
	end_getBamfiles = datetime.now()
	print('getBamfiles() %s' % (end_getBamfiles - timestart))

	cnv, bamfiles, stats = populateCnvReads(bamfilePaths, sampleInfo)
	end_populateCnvReads = datetime.now()
	print('populateCnvReads() %s' % (end_populateCnvReads - end_getBamfiles))

	cnv = populateCnvCopynum(cnv, stats, numMaleSamples, numFemaleSamples, numControlM, numControlF)
	end_populateCnvCopynum = datetime.now()
	print('populateCnvCopynum() %s' % (end_populateCnvCopynum - end_populateCnvReads))
		
	outputCnv(cnv, bamfiles)
	end_outputCnv = datetime.now()
	print('outputCnv() %s' % (end_outputCnv - end_populateCnvCopynum))

	print('timeend = %s' % (datetime.now() - timestart))
