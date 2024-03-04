# dragenAlign

This workflow will align sequence data (WG or WT) provided as fastq files to the reference sequence using Illumina Dragen. Adapter trimming is optional. The bam file will be sorted and indexed

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
`fastqR1`|File|Read 1, gzipped
`outputFileNamePrefix`|String|Prefix for output files
`reference`|String|The genome reference build. For example: hg19, hg38, mm10
`mode`|String|Specifies whether to complete genomic or transcriptomic analysis. Possible options are 'genome' or 'transcriptome'


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`fastqR2`|File?|None|Read 2 for paired-end reads, gzipped
`adapterTrim`|Boolean|true|Should adapters be trimmed, [true, trimmed]
`rgInfo`|String|"ID=1"|Comma separated list of key value pairs representing the read-group information, possible keys are (ID, SM, LB, PU, PL, CN, DS, DT, PI)


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`readGroupFormat.jobMemory`|Int|1|Memory allocated for this job
`readGroupFormat.timeout`|Int|5|Hours before task timeout
`runDragen.adapter1File`|String|"/staging/data/resources/ADAPTER1"|Adapters to be trimmed from read 1
`runDragen.adapter2File`|String|"/staging/data/resources/ADAPTER2"|Adapters to be trimmed from read 2
`runDragen.jobMemory`|Int|500|Memory allocated for this job
`runDragen.timeout`|Int|96|Hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`bam`|File|Output bam aligned to genome
`bamIndex`|File|Index for the aligned bam
`zippedOut`|File|Zip file containing the supporting .csv and .tab outputs from Dragen
`outputChimeric`|File?|Output chimeric junctions file, if available


## Commands
 This section lists command(s) run by dragenAlign workflow
 
 * Running dragenAlign
 
 === Ensures the read-group information is valid, and outputs it in the correct format ===.
 
 ``` 
     set -euo pipefail 
 
     fieldNames=("ID=" "LB=" "PL=" "PU=" "SM=" "CN=" "DS=" "DT=" "PI=") 
 
     # Split the string into an array of fields (key-value pairs)
     IFS=, read -ra rgArray <<< ~{rgInfo}
 
     # Declares an associative array and appends the values of the valid keys in rgArray
       # If duplicate fields names are present, the right-most value will be used.
     declare -A fieldsArray
 
     for field in "${rgArray[@]}"; do
       tag=${field:0:3}
       validTag=false
       for name in "${fieldNames[@]}"; do
         if [ "$tag" == "$name" ]; then
           validTag=true
           fieldsArray[${field:0:2}]=$(echo "$field" | cut -d '=' -f2)
         fi
       done
       if ! $validTag; then
         # Redirect error message to stderr
         echo "Invalid tag: '$tag'" >&2  
         exit 1
       fi
     done
 
     # Outputs the read-group information in the proper format for Dragen
     readGroupString=""
     for key in "${!fieldsArray[@]}"; do
       readGroupString+=" --RG${key} ${fieldsArray[$key]}"
     done
 
     echo "$readGroupString"
 ```
 
 === Align to reference using Dragen ===.
 
 ```
     set -euo pipefail
 
     dragen -f \
     -r ~{dragenRef} \
     -1 ~{read1} \
     ~{if (defined(read2)) then "-2 ~{read2}" else ""} \
     ~{rgInfoString} \
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
     
     zip -m ~{zipFileName} $(ls | grep '~{prefix}.*.csv\|~{prefix}.*.tab' | tr '\n' ' ')
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
