## 1.2.0 - 2024-03-26
- Workflow requires an array of fastq files with read-groups as input. Single fastq-pairs (or single fastq files) must also be inputted as an array. 
- runDragen task outputs a merged bam file
- Added new task makeCSV and headerFormat, and removed readGroupFormat task
- Replaced mode parameter with isRNA parameter (false by default)
## 1.1.0 - 2024-03-05
- Changes the zippedOut output into a zipped directory. This ensures that extraction creates a new directory instead of tarbombing the working directory.
## 1.0.0 - 2024-03-04
- Completes lane level alignments using Dragen
- Supports whole transcriptome alignment