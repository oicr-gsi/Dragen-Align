# dragenAlign

This workflow will align sequence data (WG or WT) provided as fastq files to the reference sequence using Illumina Dragen. Adapter trimming is optional. The bam file will be sorted and indexed.

## Overview

## Dependencies

* [dragen](https://developer.illumina.com/dragen)


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
`headerFormat.jobMemory`|Int|1|Memory allocated for this job
`headerFormat.timeout`|Int|5|Hours before task timeout
`makeCSV.jobMemory`|Int|1|Memory allocated for this job
`makeCSV.timeout`|Int|5|Hours before task timeout
`runDragen.adapter1File`|String|"/staging/data/resources/ADAPTER1"|Adapters to be trimmed from read 1
`runDragen.adapter2File`|String|"/staging/data/resources/ADAPTER2"|Adapters to be trimmed from read 2
`runDragen.jobMemory`|Int|500|Memory allocated for this job
`runDragen.timeout`|Int|96|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`bam`|File|Output bam aligned to genome|
`bamIndex`|File|Index for the aligned bam|
`zippedOut`|File|Zip file containing the supporting .csv and .tab outputs from Dragen|
`outputChimeric`|File?|Output chimeric junctions file, if available|


## Commands
This section lists command(s) run by dragenAlign workflow
 
* Running dragenAlign
 
### Ensures the read-group information is valid, and outputs a header for the input CSV.
 
``` 
     set -euo pipefail 
 
     headerString="Read1File,Read2File"
     
     # Split the string into an array of key-value pairs
     IFS=, read -ra rgArray <<< ~{readGroupString}
 
     # Adds valid keys (for Dragen) to headerString
     for field in "${rgArray[@]}"; do
       tag=${field:0:5}
       if [ "$tag" == "RGID=" ] || [ "$tag" == "RGLB=" ] || [ "$tag" == "RGPL=" ] || \
          [ "$tag" == "RGPU=" ] || [ "$tag" == "RGSM=" ] || [ "$tag" == "RGCN=" ] || \
          [ "$tag" == "RGDS=" ] || [ "$tag" == "RGDT=" ] || [ "$tag" == "RGPI=" ]
       then
         headerString+=",${field:0:4}"
       else
         # Redirect error message to stderr
         echo "Invalid tag: '$tag'" >&2  
         exit 1
       fi
     done
 
     # Ensures the required header information is present
     if [ "$(echo "$headerString" | grep -c "RGID")" != 1 ] || \
        [ "$(echo "$headerString" | grep -c "RGSM")" != 1 ] || \
        [ "$(echo "$headerString" | grep -c "RGLB")" != 1 ] || \
        [ "$(echo "$headerString" | grep -c "RGPU")" != 1 ]; then
       echo "Missing required read-group information from header" >&2  
       exit 1
     fi
 
     echo "$headerString"
```
 
### Format input CSV file for Dragen.
 
``` 
     set -euo pipefail 
     
     echo ~{csvHeader} > ~{csvResult}
 
     # Load arrays into bash variables
     arrRead1s=(~{sep=" " read1s})
     if ~{isPaired}; then arrRead2s=(~{sep=" " read2s}); fi
     arrReadGroups=(~{sep=" " readGroups})
     
     # Iterate over the arrays concurrently
     for (( i = 0; i < ~{arrayLength}; i++ ))
     do
       read1="${arrRead1s[i]}"
       if ~{isPaired}; then read2="${arrRead2s[i]}"; else read2=""; fi
       readGroup=$(echo "${arrReadGroups[i]}" | sed 's/RG..=//g')
       echo "$read1,$read2,$readGroup" >> ~{csvResult}
     done
```
 
### Align to reference using Dragen.
 
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
