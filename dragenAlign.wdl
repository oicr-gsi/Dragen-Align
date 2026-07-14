version 1.0

struct InputGroup {
  File fastqR1
  File? fastqR2
  String readGroup
}

struct GenomeResources {
    String referenceDirectory
    String dragenVersion
}

workflow dragenAlign {

  input {
    Array[InputGroup] inputGroups
    String outputFileNamePrefix
    String reference
    Boolean adapterTrim = true
    Boolean isRNA = false
  }

  Boolean isPaired = if (defined(inputGroups[0].fastqR2)) then true else false

  scatter (ig in inputGroups) {
    File read1s = ig.fastqR1
    String readGroups = ig.readGroup
  }

  if(isPaired) {
    scatter (ig in inputGroups) {
      # Workaround for converting File? type to File
      File read2s = select_all([ig.fastqR2])[0]
    }
  }
  
  Map[String,GenomeResources] dragenRef_by_genome = { 
    "hg38": {
      "referenceDirectory": "/.mounts/labs/gsiprojects/gsi/Dragen/reference/hg38fa.p12/",  # /staging/data/references/hg38-p12.v9
      "dragenVersion": "4.2.4"
    },
    "hg38_noAlt": {
      "referenceDirectory": "/.mounts/labs/gsiprojects/gsi/Dragen/reference/hg38_noAlt-p12",
      "dragenVersion": "4.2.4"
    }
  }

  String dragenRef = dragenRef_by_genome[reference].referenceDirectory
  String dragen_version = dragenRef_by_genome[ reference ].dragenVersion
  
  parameter_meta {
    inputGroups: "Array of fastq files to align using Dragen. Read-group information is required for fastq files, with the following fields being non-optional: RGID, RGSM, RGLB, RGPU. Each FASTQ file can only be referenced once."
    outputFileNamePrefix: "Prefix for output files"
    reference: "The genome reference build. For example: hg19, hg38, mm10"
    adapterTrim: "Should adapters be trimmed, [true, trimmed]"
    isRNA: "Specifies whether to complete transcriptomic analysis, [false, genomic]"
  }

  scatter(t in inputGroups) {
    call extractInfoLine {
      input:
      fastqInput = object{fastqR1: t.fastqR1, fastqR2: t.fastqR2, readGroup: t.readGroup}
    }
  }

  call composeList {
    input:
      inputLines = extractInfoLine.outputLine,
      outputFileName = "dragen_inputs.csv"
  }

  call runDragen  {
    input: 
    csv = composeList.inputList,
    dragenRef = dragenRef,
    dragenVersion = dragen_version,
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
      },
      {
        name: "gsi modules : dragen-scripts/0.3",
        url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
      }
    ]
    output_meta: {
        bam: {
            description: "BAM file with alignments",
            vidarr_label: "bam"
        },
        bamIndex: {
            description: "index of BAM file with alignments",
            vidarr_label: "bamIndex"
        },
        zippedOut: {
            description: "Zipped .csv and .tab files (additional outputs)",
            vidarr_label: "zippedOut"
        },
        outputChimeric: {
            description: "Optional output file with chimeric junctions",
            vidarr_label: "outputChimeric"
        }
   }
  }

  output {
    File bam = runDragen.bam
    File bamIndex = runDragen.bamIndex
    File zippedOut = runDragen.zippedOut
    File? outputChimeric = runDragen.outputChimeric
  }

}

# =====================================================================
# A scripted extraction of info from RG line to dragen-compliant string
# =====================================================================
task extractInfoLine {
   input {
       InputGroup fastqInput
       String parsingScript = "$DRAGEN_SCRIPTS_ROOT/bin/composeList.py"
       Int timeout = 4
       Int jobMemory = 4
       String modules = "dragen-scripts/0.3"
   }

   parameter_meta {
     fastqInput: "InputGroup struct entry with fastq files"
     parsingScript: "Script for parsing inputs into a line"
     timeout: "Timeout for the job"
     jobMemory: "Job allocated RAM"
     modules: "dependency modules"
   }

   command <<<
    python3 ~{parsingScript} -i ~{write_json(fastqInput)}
   >>>

   runtime {
     timeout: "~{timeout}"
     modules: "~{modules}"
     memory:  "~{jobMemory} GB"
   }

   output {
     String outputLine = read_string(stdout())
   }

   meta {
     output_meta: {
       outputLine: "Output line to use in a list of fastq files in dragen-compliant format"
     }
   }
}

# =====================================================================
#  Compose a dragen-compliant list of inputs to use with snv caller
# =====================================================================
task composeList {
   input  {
      Array[String] inputLines
      String listWritingScript = "$DRAGEN_SCRIPTS_ROOT/bin/writeFile.py"
      String outputFileName
      Int jobMemory = 4
      Int timeout = 4
      String modules = "dragen-scripts/0.3"
   }

   parameter_meta {
     inputLines: "Array of input lines to print"
     listWritingScript: "Script for writing out list of inputs"
     outputFileName: "Name of an output file, list of inputs"
     jobMemory: "Job allocated RAM"
     timeout: "Timeout for the job"
     modules: "dependency modules"
   }

   command<<<
   python3 ~{listWritingScript} -o ~{outputFileName} -l "~{sep=';' inputLines}"
   >>>


   runtime {
      timeout: "~{timeout}"
      modules: "~{modules}"
      memory:  "~{jobMemory} GB"
   }

   output {
     File inputList = "~{outputFileName}"
   }

   meta {
     output_meta: {
       inputList: "Output file to use with dragen SNV caller"
     }
   }
}


# ================================================================
# Main task for generating SNV calls in somatic mode (DRAGEN mode)
#
# we need CSV files with a header and data lines organized as:
#
# RGID Read Group
# RGSM Sample ID
# RGLB Library
# Lane Flow cell lane
# Read1File - Full path to a valid FASTQ input file
# Read2File - Full path to a valid FASTQ input file. Required for paired-end input. If not using paired-end input, leave empty.
# Each FASTQ file can only be referenced once in the CSV list.
# All values in the Read2File column must be reference valid files or must all be empty.
# ================================================================
task runDragen {
  input {
    File csv
    String dragenRef
    String dragenVersion
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
    dragenVersion: "Expected version of dragen software on the DRAGEN node"
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
    dragen_version: "~{dragenVersion}"
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
