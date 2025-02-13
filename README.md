# dragenAlign

This workflow will align sequence data (WG or WT) provided as fastq files to the reference sequence using Illumina Dragen. Adapter trimming is optional. The bam file will be sorted and indexed.

## Overview

## Dependencies

* [dragen](https://developer.illumina.com/dragen)
* [gsi modules : dragen-scripts 0.3](https://gitlab.oicr.on.ca/ResearchIT/modulator)

## Usage

### Cromwell
```
java -jar cromwell.jar run dragenAlign.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputGroups`|Array[InputGroup]|Array of fastq files to align using Dragen. Read-group information is required for fastq files, with the following fields being non-optional: RGID, RGSM, RGLB, RGPU. Each FASTQ file can only be referenced once.
`outputFileNamePrefix`|String|Prefix for output files
`reference`|String|The genome reference build. For example: hg19, hg38, mm10


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`adapterTrim`|Boolean|true|Should adapters be trimmed, [true, trimmed]
`isRNA`|Boolean|false|Specifies whether to complete transcriptomic analysis, [false, genomic]


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`extractInfoLine.parsingScript`|String|"$DRAGEN_SCRIPTS_ROOT/bin/composeList.py"|Script for parsing inputs into a line
`extractInfoLine.timeout`|Int|4|Timeout for the job
`extractInfoLine.jobMemory`|Int|4|Job allocated RAM
`extractInfoLine.modules`|String|"dragen-scripts/0.1"|dependency modules
`composeList.listWritingScript`|String|"$DRAGEN_SCRIPTS_ROOT/bin/writeFile.py"|Script for writing out list of inputs
`composeList.jobMemory`|Int|4|Job allocated RAM
`composeList.timeout`|Int|4|Timeout for the job
`composeList.modules`|String|"dragen-scripts/0.1"|dependency modules
`runDragen.adapter1File`|String|"/staging/data/resources/ADAPTER1"|Adapters to be trimmed from read 1
`runDragen.adapter2File`|String|"/staging/data/resources/ADAPTER2"|Adapters to be trimmed from read 2
`runDragen.jobMemory`|Int|500|Memory allocated for this job
`runDragen.timeout`|Int|96|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`bam`|File|BAM file with alignments|vidarr_label: bam
`bamIndex`|File|index of BAM file with alignments|vidarr_label: bamIndex
`zippedOut`|File|Zipped .csv and .tab files (additional outputs)|vidarr_label: zippedOut
`outputChimeric`|File?|Optional output file with chimeric junctions|vidarr_label: outputChimeric


## Commands
This section lists command(s) run by dragenAlign workflow
 
* Running dragenAlign
 
### Ensures the read-group information is valid, and outputs a header for the input CSV.
 
```
     python3 ~{parsingScript} -i ~{write_json(fastqInput)}
```
 
### Compose a list of inputs for dragen:wq

```
    python3 ~{listWritingScript} -o ~{outputFileName} -l "~{sep=';' inputLines}"
```
 
```
     set -euo pipefail
 
     dragen -f \
     -r ~{dragenRef} \
     --fastq-list ~{csv} \
     --fastq-list-all-samples true \
     --enable-map-align true \
     --enable-map-align-output true \
     --output-directory ./ \
     --output-file-prefix ~{prefix} \
     ~{if (adapterTrim) then "--read-trimmers adapter" +
                             " --trim-adapter-read1 ~{adapter1File}" +
                             " --trim-adapter-read2 ~{adapter2File}" else ""} \
     --trim-min-length 1 \
     --enable-bam-indexing true \
     --enable-sort true \
     --enable-duplicate-marking false \
     ~{if (isRNA) then "--enable-rna true" else ""}
     
     mkdir ~{zipFileName}
     cp -t ~{zipFileName} $(ls | grep '~{prefix}.*.csv\|~{prefix}.*.tab' | tr '\n' ' ')
     zip -r ~{zipFileName}.zip ~{zipFileName}
```
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
