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
