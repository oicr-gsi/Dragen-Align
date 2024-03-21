version 1.0

struct InputGroup {
  File fastqR1
  File fastqR2
  String readGroup
}

workflow dragenAlign {

  input {
    Array[InputGroup] inputGroups
    String outputFileNamePrefix
    String reference
    Boolean adapterTrim = true
    String isRNA
  }

  scatter (ig in inputGroups) {
    File read1s = ig.fastqR1
    File read2s = ig.fastqR2
    String readGroups = ig.readGroup
  }
  
  Map[String,String] dragenRef_by_genome = { 
    "hg38": "/staging/data/references/hg38fa.v9"
  }

  String dragenRef = dragenRef_by_genome[reference]
  
  parameter_meta {
    inputGroups: "Array of fastq files to align using Dragen. Read-group information is required for fastq files, with the following fields being non-optional: RGID, RGSM, RGLB, RGPU. Each FASTQ file can only be referenced once."
    outputFileNamePrefix: "Prefix for output files"
    reference: "The genome reference build. For example: hg19, hg38, mm10"
    adapterTrim: "Should adapters be trimmed, [true, trimmed]"
    isRNA: "Specifies whether to complete transcriptomic analysis. Possible options are 'true' or 'false'"
  }

  call makeCSV {  
    input: 
    read1s = read1s,
    read2s = read2s,
    prefix = outputFileNamePrefix,
    readGroups = readGroups
  }

  call runDragen  { 
    input: 
    csv = makeCSV.outCSV,
    dragenRef = dragenRef,
    adapterTrim = adapterTrim,
    prefix = outputFileNamePrefix,
    isRNA = isRNA
  }

  meta {
    author: "Lawrence Heisler and Muna Mohamed"
    email: "lheisler@oicr.on.ca and mmohamed@oicr.on.ca"
    description: "This workflow will align sequence data (WG or WT) provided as fastq files to the reference sequence using Illumina Dragen. Adapter trimming is optional. The bam file will be sorted and indexed."
    dependencies: [
      {
        name: "dragen",
        url: "https://developer.illumina.com/dragen"
      }
    ]
  }

  output {
    File bam = runDragen.bam
    File bamIndex = runDragen.bamIndex
    File zippedOut = runDragen.zippedOut
    File? outputChimeric = runDragen.outputChimeric
  }

}

task makeCSV { 
  input { 
    Array[File] read1s
    Array[File]? read2s #Check if it should be optional?
    Array[File] readGroups
    String prefix
    Int jobMemory = 1
    Int timeout = 5
  }

  parameter_meta { 
    read1s: "Array of read 1 fastq files" 
    read2s: "Array of read 2 fastq files. May be empty if single-read" 
    readGroups: "Array of read-group information to be added into the bam file header" 
    prefix: "Prefix for output files"
    jobMemory: "Memory allocated for this job" 
    timeout: "Hours before task timeout" 
  } 
  
  String csvResult = "{prefix}_dragenInput.csv"
  #Required in CSV: 
  #RGID,RGSM,RGLB,Lane (RGPU),Read1File,Read2File 

  command <<< 
    set -euo pipefail 

    # Allowed RG fields in Dragen
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
  >>> 

  runtime { 
    memory: "~{jobMemory} GB" 
    timeout: "~{timeout}" 
  } 

  output { 
    File outCSV = ~{csvResult}
  }

  meta { 
    output_meta: { 
      outCSV: "Formatted csv input for Dragen, containing fastq files and read-group information" 
    } 
  } 
} 

task runDragen {
  input {
    File csv
    String dragenRef
    String prefix
    Boolean isRNA
    Boolean adapterTrim
    String adapter1File = "/staging/data/resources/ADAPTER1"
    String adapter2File = "/staging/data/resources/ADAPTER2"
    Int jobMemory = 500
    Int timeout = 96
  }

  parameter_meta {
    csv: "Formatted csv input for Dragen, containing fastq files and read-group information"
    dragenRef: "The reference genome to align the sample with by Dragen"
    prefix: "Prefix for output files"
    isRNA: "True/False, whether to complete transcriptomic analysis"
    adapterTrim: "True/False for adapter trimming"
    adapter1File: "Adapters to be trimmed from read 1"
    adapter2File: "Adapters to be trimmed from read 2"
    jobMemory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
  }
  
  String zipFileName = "~{prefix}_additional_outputs"

  command <<<
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
  >>>

  runtime {
    timeout: "~{timeout}"
    backend: "DRAGEN"
  }
  
  output {
    File bam = "~{prefix}.bam"
    File bamIndex = "~{prefix}.bam.bai"
    File zippedOut = "~{zipFileName}.zip"
    File? outputChimeric = "~{prefix}.Chimeric.out.junction"
  }

  meta {
    output_meta: {
      bam: "Output bam aligned to genome",
      bamIndex: "Index for the aligned bam",
      zippedOut: "Zip file containing the supporting .csv and .tab outputs from Dragen",
      outputChimeric: "Output chimeric junctions file, if available"
    }
  }
}