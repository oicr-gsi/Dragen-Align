version 1.0

workflow dragenAlign {

  input {
    File fastqR1
    File fastqR2
    String outputFileNamePrefix
    String reference
    Boolean adapterTrim = true
    String rgInfoString = "--RGID 1"
  }

  parameter_meta {
    fastqR1: "R1 of the fastq paired, gipped"
    fastqR2: "R2 of the fastq pair, gzipped"
    outputFileNamePrefix: "Prefix for output files"
    reference: "The genome reference build. For example: hg19, hg38, mm10"
    adapterTrim: "Should adapters be trimmed"
    rgInfoString: "space separated list of space separated key value pairs, possible keys are --RGID,--RGSM,--RGLB,--RGPU"
  }

  Map[String,String] dragenRef_by_genome = { 
    "hg38": "/staging/data/references/hg38fa.v9"
  }

  String dragenRef = dragenRef_by_genome [ reference ]

  call runDragen  { 
    input: 
    read1 = fastqR1,
    read2 = fastqR2,
    dragenRef = dragenRef,
    adapterTrim = adapterTrim,
    prefix = outputFileNamePrefix,
    rgInfoString = rgInfoString
  }

  meta {
    author: "Lawrence Heisler"
    email: "lheisler@oicr.on.ca"
    description: "This workflow will align a fastq pair to the reference seqeunce using Illumina Dragen.  Adapter trimming is optional.  The bam file will be sorted and indexed"
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
    File metrics = runDragen.metrics
  }
}

task runDragen {
    input {
      File read1
      File read2
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
      -1 $read1 -2 $read2 \
      ~{rgInfoString} \
      --enable-map-align true \
      --enable-map-align-output true \
      --output-directory /staging/data/scratch \
      --output-file-prefix ~{prefix} \
      --read-trimmers adapter \
      --trim-adapter-read1 ~{adapter1File} \
      --trim-adapter-read2 ~{adapter2File} \
      --trim-min-length 1 \
      --enable-bam-indexing true \
      --enable-sort true \
      --enable-duplicate-marking false
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
        bam: "output bam aligned to genome",
        bamIndex: "index for the aligned bam",
        metrics: "mapping metrics"
      }
    }
  }
  



