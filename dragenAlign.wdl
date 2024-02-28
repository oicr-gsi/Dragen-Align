version 1.0

workflow dragenAlign {

  input {
    File fastqR1
    File? fastqR2
    String outputFileNamePrefix
    String reference
    Boolean adapterTrim = true
    String rgInfo = "ID=1"
    String mode
  }

  parameter_meta {
    fastqR1: "Read 1 of the fastq pair, gzipped"
    fastqR2: "Read 2 of the fastq pair, gzipped"
    outputFileNamePrefix: "Prefix for output files"
    reference: "The genome reference build. For example: hg19, hg38, mm10"
    adapterTrim: "Should adapters be trimmed, [true, trimmed]"
    rgInfo: "Comma separated list of key value pairs representing the read-group information, possible keys are (ID=, SM=, LB=, PU=, PL=, CN=, DS=, DT=, PI=)"
    mode: "Specifies whether to complete genomic or transcriptomic analysis. Possible options are 'genome' or 'transcriptome'"
  }

  Map[String,String] dragenRef_by_genome = { 
    "hg38": "/staging/data/references/hg38fa.v9"
  }

  String dragenRef = dragenRef_by_genome [ reference ]
  
  # Validating the read-group information
  call readGroupFormat {  
    input: 
    rgInfo = rgInfo
  }

  call runDragen  { 
    input: 
    read1 = fastqR1,
    read2 = fastqR2,
    dragenRef = dragenRef,
    adapterTrim = adapterTrim,
    prefix = outputFileNamePrefix,
    mode = mode,
    rgInfoString = readGroupFormat.rgInfoValid
  }

  meta {
    author: "Lawrence Heisler and Muna Mohamed"
    email: "lheisler@oicr.on.ca and mmohamed@oicr.on.ca"
    description: "This workflow will align sequence data provided as fastq files to the reference sequence using Illumina Dragen.  Adapter trimming is optional.  The bam file will be sorted and indexed"
    dependencies: [
      {
        name: "dragen",
        url: "https://developer.illumina.com/dragen"
      }
    ]
  }

  output {
    File? bam = runDragen.bam
    File? bamIndex = runDragen.bamIndex
    File? metrics = runDragen.metrics
  }

}

task readGroupFormat { 
  input { 
    String rgInfo
    Int jobMemory = 1
    Int timeout = 5
  }

  parameter_meta { 
    rgInfo: "The read-group information to be added into the bam file header" 
    jobMemory: "Memory allocated for this job" 
    timeout: "Hours before task timeout" 
  } 
  
  command <<< 
    set -euo pipefail 

    fieldNames=("ID=" "LB=" "PL=" "PU=" "SM=" "CN=" "DS=" "DT=" "PI=") 

    # Split the string into an array of fields (key-value pairs)
    IFS=, read -ra rgArray <<< ~{rgInfo}

    # Declares an associative array to append the values of the valid keys in rgArray
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
    String rgInfoValid = read_string(stdout())
  }

  meta { 
    output_meta: { 
      rgInfoValid: "Formatted string containing the validated read-group information" 
    } 
  } 
} 

task runDragen {
  input {
    File read1
    File? read2
    String dragenRef
    String prefix
    String mode
    Boolean adapterTrim
    String adapter1File = "/staging/data/resources/ADAPTER1"
    String adapter2File = "/staging/data/resources/ADAPTER2"
    String rgInfoString
    Int jobMemory = 500
    Int timeout = 96
  }

  parameter_meta {
    read1: "Fastq file for read 1"
    read2: "Fastq file for read 2"
    dragenRef: "The reference genome to align the sample with by Dragen"
    prefix: "Prefix for output files"
    mode: "Specifies whether to complete genomic or transcriptomic analysis. Possible options are 'genome' or 'transcriptome'"
    adapterTrim: "True/False for adapter trimming"
    adapter1File: "Adapters to be trimmed from read 1"
    adapter2File: "Adapters to be trimmed from read 2"
    rgInfoString: "Formatted string containing the validated read-group information"
    jobMemory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
  }
  
  # Boolean indicating whether to enable transcriptomic analysis
  Boolean enableRNA = if mode == "transcriptome" then true else false

  command <<<
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
    --enable-rna ~{enableRNA}
  >>>

  runtime {
    timeout: "~{timeout}"
    backend: "DRAGEN"
  }
  
  output {
    File bam = "~{prefix}.bam"
    File bamIndex = "~{prefix}.bam.bai"
    File metrics = "~{prefix}.mapping_metrics.csv"
  }

  meta {
    output_meta: {
      bam: "Output bam aligned to genome",
      bamIndex: "Index for the aligned bam",
      metrics: "Mapping metrics"
    }
  }
}