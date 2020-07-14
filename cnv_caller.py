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
Copy number variation (CNV) caller.

Usage example:
python3 cnv_caller.py --samplesDir ./DS_TD_bam --bin 5000 --sampleInfo ./DS_TD_sample_info.csv --controlName Control --testCovariate Genotype

"""

# Input command line arguments
parser = argparse.ArgumentParser(
    description = 'CNV caller.',
    formatter_class = argparse.RawDescriptionHelpFormatter,
    epilog = extended_help)
parser.add_argument(
    '--samplesDir', 
    required = True, 
    type = str, 
    default = None, 
    metavar = '<path>', 
    help = 'Path to a directory containing all BAM file samples or subdirectories of BAM files.')
parser.add_argument(
    '--binSize', 
    required = False, 
    type = int, 
    default = 1000,
    metavar = '<int>', 
    help = 'Size of bin.')
parser.add_argument(
    '--sampleInfo',
    required = True,
    type = str,
    metavar = '<path>',
    help = 'Path to sample info .csv file.')
parser.add_argument(
    '--controlName',
    required = False,
    type = str,
    default = 'Control',
    metavar = '<str>',
    help = 'Name of control class.')
parser.add_argument(
    '--testCovariate',
    required = False,
    type = str,
    default = 'Diagnosis',
    metavar = '<str>',
    help = 'Name of test covariate class (e.g. Diagnosis, Genotype).')

arg = parser.parse_args()

# TODO: If sampleInfo has no Sex column, then call python script to run SexChecker and read in predicted sex column as Sex
# OR normalize with autosomal chromosomes, don't do male and female normalization factor

# TODO: List control samples first in output file, currently done alphabetically

# TODO: Option to read in excel sample info file

# TODO: Messages to console and delete .print file?

# Write print statement outputs to file
sys.stdout = open('CNV_Me_' + datetime.now().strftime('%I:%M%p_%b%d') + '.print', 'w')

# Chromosomes of interest for hg38
chromosomes = []
for i in range (1, 23): # 1..22
        chromosomes.append('chr' + str(i))
chromosomes.extend(["chrX", "chrY"])


def getBamfiles(path, bamfilePaths):
    """
    Function:
    Traverse though path recursively to find all BAM files.

    Parameters:
    path (str): path to directory of BAM files or subdirectories of BAM files.
    bamfilePaths (list): list containing currently found BAM file paths.

    Returns:
    bamfilePaths (list): list containing all BAM file paths.
    """
    for item in os.scandir(path):
        if item.is_file() and item.path.endswith('.bam'):
            bamfilePaths.append(item.path)
        elif item.is_dir():
            bamfilePaths = getBamfiles(item.path, bamfilePaths)
    return (bamfilePaths)

def getSampleNumbers(bamfilePaths):
    """
    Function:
    Find sample names and number of control and experimental samples.

    Parameters:
    bamfilePaths (list): list containing all BAM file paths.
    
    Returns:
    samplenames (list): list of sample names.
    numControl (int): number of control samples.
    numControlM (int): number of male control samples.
    numControlF (int): number of female control samples.
    numExp (int): number of experimental samples.

    """
    numControl = 0
    numControlM = 0
    numControlF = 0
    numExp = 0
    samplenames = []
    for bamfile in bamfilePaths:
        filename = bamfile.split("/")[-1]
        samplename = str(filename.split("_")[0])
        samplenames.append(samplename)
        if sampleInfo.loc[samplename, arg.testCovariate] == arg.controlName:
            numControl += 1 
            if sampleInfo.loc[samplename, 'Sex'] == 'M':
                numControlM += 1
            elif sampleInfo.loc[samplename, 'Sex'] == 'F':
                numControlF += 1
        else:
            numExp += 1
    return (samplenames, numControl, numControlM, numControlF, numExp)


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
        # stats[bamfile]["diagnosis"] = sampleInfo.loc[samplename, 'Diagnosis']   
        stats[bamfile][arg.testCovariate] = sampleInfo.loc[samplename, arg.testCovariate]   

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
                    binkey = int(pos/arg.binSize)
                    if binkey in cnv[chrm]:
                        if bamfile in cnv[chrm][binkey]:
                        # Binkey and sample exist
                            cnv[chrm][binkey][bamfile]['sumReadCounts'] += 1
                            cnv[chrm][binkey][bamfile]['sumReadLengths'] += (1*readLength)
                        else:
                        # Binkey exists but add new sample
                            cnv[chrm][binkey].update({bamfile: {\
                                'sumReadCounts': 1,\
                                'sumReadLengths': 1*readLength,\
                                'copynum': 0
                            }})
                    else: 
                        # Add new binkey and new sample
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
        # Control samples
        if stats[bamfile][arg.testCovariate] == arg.controlName:
            if bamfile in cnv[chrm][binkey].keys():
                # Coverage      = sum of (reads * read length) in each bin / bin size
                # To normalize: divide by total coverage of sample = sum of (reads * read length) in sample
                binCoverage = cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.binSize
                totalReadsAndLengths = stats[bamfile]['totalLen']
                normalizedBinCoverage = binCoverage / totalReadsAndLengths
                
                # If bamfile is from female, stays the same, diploid X, no Y
                # If bamfile is from male, use normfactorSexChrm and normalize to 1x each

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
                                (cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.binSize) / \
                                stats[bamfile]['totalLen'] / normfactorMaleX
                        if sampleSex == "F" and normfactorFemaleX != 0:
                            cnv[chrm][binkey][bamfile]['copynum'] = \
                                (cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.binSize) / \
                                stats[bamfile]['totalLen'] / (0.5 * normfactorFemaleX)
                    elif chrm == "chrY":
                        if sampleSex == "M" and normfactorMaleY != 0:
                            cnv[chrm][binkey][bamfile]['copynum'] = \
                                (cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.binSize) / \
                                stats[bamfile]['totalLen'] / normfactorMaleY
                        if sampleSex == "F" and normfactorFemaleY != 0:
                            cnv[chrm][binkey][bamfile]['copynum'] = \
                                (cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.binSize) / \
                                stats[bamfile]['totalLen'] / normfactorFemaleY
                    else:
                        if normfactorAut != 0:
                            cnv[chrm][binkey][bamfile]['copynum'] = \
                                (cnv[chrm][binkey][bamfile]['sumReadLengths'] / arg.binSize) / \
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
    with open('CNV_Me_output_' + datetime.now().strftime('%I:%M%p_%b%d') + '.txt', 'w', newline='') as outfile:
        outfile = csv.writer(outfile, delimiter='\t')
        samples.sort()
        # TODO: done alphabetically, list control samples first
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
                        binkey*arg.binSize,
                        (binkey+1)*arg.binSize,
                        *copynumList
                    ])
    return None


if __name__ == '__main__':
    timestart = datetime.now()  
    print('timestart = %s' % timestart)
    print(arg.samplesDir)
    print(arg.binSize)
    print(arg.sampleInfo)

    sampleInfo = pd.read_csv(arg.sampleInfo, converters={i: str for i in range(len(pd.read_csv(arg.sampleInfo)))})
    sampleInfo.index = list(sampleInfo['Name'])

    # Convert sex label to M and F
    for i in range(0, len(sampleInfo)):
        if ('M' or 'm' or 'Male' or 'male') in sampleInfo.iloc[i]['Sex']:
            sampleInfo.iloc[i]['Sex'] = 'M'
        if ('F' or 'f' or 'Female' or 'female') in sampleInfo.iloc[i]['Sex']:
            sampleInfo.iloc[i]['Sex'] = 'F'

    numMaleSamples = len(sampleInfo[sampleInfo['Sex'] == 'M'])
    numFemaleSamples = len(sampleInfo[sampleInfo['Sex'] == 'F'])    

    bamfilePaths = []
    bamfilePaths = getBamfiles(arg.samplesDir, bamfilePaths)

    samplenames, numControl, numControlM, numControlF, numExp = getSampleNumbers(bamfilePaths)
    print(*bamfilePaths, sep='\n')
    print(numControl)
    print(numExp)
    print(samplenames)
    print(numControlM)
    print(numControlF)

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
