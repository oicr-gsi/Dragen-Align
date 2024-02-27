version 1.0

workflow dragenAlign {

  input {
    File fastqR1
    File? fastqR2
    String outputFileNamePrefix
    String reference
    Boolean adapterTrim = true
    String rgInfoString = "ID=1"
    String mode
  }

  parameter_meta {
    fastqR1: "Read 1 of the fastq pair, gzipped"
    fastqR2: "Read 2 of the fastq pair, gzipped"
    outputFileNamePrefix: "Prefix for output files"
    reference: "The genome reference build. For example: hg19, hg38, mm10"
    adapterTrim: "Should adapters be trimmed, [true, trimmed]"
    rgInfoString: "Comma separated list of key value pairs representing the read-group information, possible keys are (ID=, SM=, LB=, PU=, PL=, CN=, DS=, DT=, PI=)"
    mode: "Specifies whether to complete genomic or transcriptomic analysis. Possible options are 'genome' or 'transcriptome'"
  }

  Map[String,String] dragenRef_by_genome = { 
    "hg38": "/staging/data/references/hg38fa.v9"
  }

  String dragenRef = dragenRef_by_genome [ reference ]
  
  # Validating the read-group information
  call readGroupCheck {  
      input: 
      rgInfoString = rgInfoString 
  }

  if (readGroupCheck.validReadGroups) {
    call runDragen  { 
      input: 
      read1 = fastqR1,
      read2 = fastqR2,
      dragenRef = dragenRef,
      adapterTrim = adapterTrim,
      prefix = outputFileNamePrefix,
      rgInfoString = rgInfoString
    }
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

task readGroupCheck { 

  input { 
    String rgInfoString 
    Int jobMemory = 1
    Int timeout = 1 
  } 

  parameter_meta { 
    rgInfoString: "The read-group information to be added into the bam file header" 
    jobMemory: "Memory allocated for this job" 
    timeout: "Hours before task timeout" 
  } 
  
  command <<< 
    set -euo pipefail 

    fieldNames=("ID=" "LB=" "PL=" "PU=" "SM=" "CN=" "DS=" "DT=" "PI=") 
    
    # Split the string into an array 
    IFS=, read -ra readFields <<< ~{rgInfoString}

    for field in "${readFields[@]}"; do 
        tag=${field:0:3}
        validTag=false
        for name in "${fieldNames[@]}"; do
            if [ "$tag" == "$name" ]; then
                validTag=true
            fi
        done 
        if ! $validTag; then
            # Redirect error message to stderr
            echo "Invalid tag: '$tag'" >&2  
            exit 1
        fi
    done 
  >>> 

  runtime { 
      memory: "~{jobMemory} GB" 
      timeout: "~{timeout}" 
  } 

  output { 
      Boolean validReadGroups = true
  } 

  meta { 
      output_meta: { 
          validReadGroups: "Boolean specifying if the read-group information is valid" 
      } 
  }  

} 

task runDragen {
    input {
      File read1
      File? read2
      String dragenRef
      String prefix
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
      adapterTrim: "True/False for adapter trimming"
      adapter1File: "Adapters to be trimmed from Read1"
      adapter2File: "Adapters to be trimmed from Read2"
      jobMemory: "Memory allocated for this job"
      timeout: "Hours before task timeout"
    }

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
                              "--trim-adapter-read1 ~{adapter1File}" +
                              "--trim-adapter-read2 ~{adapter2File}" else ""} \
      --trim-min-length 1 \
      --enable-bam-indexing true \
      --enable-sort true \
      --enable-duplicate-marking false \
      ~{if (mode == "transcriptome") then "--enable-rna true" else ""}
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