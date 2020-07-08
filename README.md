# CNV_Me
Copy number variation (CNV) caller to detect large-scale structural variation using BAM (.bam) files.

## Command line arguments
`--samplesAll` <path> Path to all BAM samples. <br>
`--samplesControl` <path> Path to a directory of BAM files for control samples. <br>
`--samplesExperimental` <path> Path to a directory of BAM files for experimental samples. <br>
`--window` <int> Size of window. Default is 1000. <br>
`--sampleInfo` <path> Path to a sample info csv file. <br>
`--controlName` <str> Name of the control class. Default is "Control".

### Workflow
![bin cov](https://latex.codecogs.com/gif.latex?Coverage%5C%3Aper%5C%3Abin%20%3D%20%5Cfrac%7B%5Csum_%7Bi%20%3D%20bin_%7Bfirst%7D%7D%5E%7Bbin_%7Blast%7D%7D%5B%28Read%5C%3Acount%29%20%5Ctimes%20%28Read%5C%3Alength%29%5D_i%7D%7BBin%5C%3Asize%7D)

![normalized bin cov](https://latex.codecogs.com/gif.latex?Normalized%5C%3Abin%5C%3Acoverage%20%3D%20%5Cfrac%7BCoverage%5C%3Aper%5C%3Abin%7D%7B%5Csum_%7Bi%20%3D%20Sample_%7Bfirst%7D%7D%5E%7BSample_%7Blast%7D%7D%5B%28Read%5C%3Acount%29%20%5Ctimes%20%28Read%5C%3Alength%29%5D_i%7D)

![normalization factor autosome](https://latex.codecogs.com/gif.latex?Normalization%5C%3Afactor_%7BAutosome%7D%20%3D%20%5Cfrac%7B%5Csum_%7Bi%20%3D%20ControlSample_%7Bfirst%7D%7D%5E%7BControlSample_%7Blast%7D%7D%5B%280.5%29%20%5Ctimes%20%28Normalized%5C%3Abin%5C%3Acoverage%29%5D_i%7D%7BNumber%5C%3Acontrol%5C%3Asamples%7D)

![normalization factor M ChrX](https://latex.codecogs.com/gif.latex?Normalization%5C%3Afactor_%7BMaleChrX%7D%20%3D%20%5Cfrac%7B%5Csum_%7Bi%20%3D%20MaleControlSample_%7Bfirst%7D%7D%5E%7BMaleControlSample_%7Blast%7D%7D%5BNormalized%5C%3Abin%5C%3Acoverage%5D_i%7D%7BNumber%5C%3Amale%5C%3Acontrol%5C%3Asamples%7D)

![normalization factor M ChrY](https://latex.codecogs.com/gif.latex?Normalization%5C%3Afactor_%7BMaleChrY%7D%20%3D%20%5Cfrac%7B%5Csum_%7Bi%20%3D%20MaleControlSample_%7Bfirst%7D%7D%5E%7BMaleControlSample_%7Blast%7D%7D%5BNormalized%5C%3Abin%5C%3Acoverage%5D_i%7D%7BNumber%5C%3Amale%5C%3Acontrol%5C%3Asamples%7D)

![normalization factor F ChrX](https://latex.codecogs.com/gif.latex?Normalization%5C%3Afactor_%7BFemaleChrX%7D%20%3D%20%5Cfrac%7B%5Csum_%7Bi%20%3D%20FemaleControlSample_%7Bfirst%7D%7D%5E%7BFemaleControlSample_%7Blast%7D%7D%5B%280.5%29%20%5Ctimes%20%28Normalized%5C%3Abin%5C%3Acoverage%29%5D_i%7D%7BNumber%5C%3Afemale%5C%3Acontrol%5C%3Asamples%7D)

![normalized copy number](https://latex.codecogs.com/gif.latex?Normalized%5C%3Acopy%5C%3Anumber%20%3D%20%5Cfrac%7BNormalized%5C%3Abin%5C%3Acoverage%7D%7BNormalization%5C%3Afactor%7D)

![Approach flowchart](cnv_flowchart.svg)

## Output
`CNV_Me.print` All print statments are recorded in this file. <br>
`CNV_Me_output.txt` Output tab-delimited file of CNV values for each window of each chromosome of each sample. The columns are 

### Example of output file
| chr | start | end | control_sample_1 | control_sample_2 | experimental_sample_1 | experimental_sample_2 |
| --------------- | --------------- | --------------- | --------------- | --------------- | --------------- | --------------- |
| chr1 | 10000 | 15000 | 2.217 | 2.283 | 1.627 | 1.504 |
| ... | ... | ... | ... | ... | ... | ... |
| chrY | 2785000 | 2790000 | 1.11 |	1.038 | 1.054	| 0.814|
