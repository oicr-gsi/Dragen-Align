version 1.0

workflow dragenSomaticVariantCaller {
    input {
        File fastqR1
        File fastqR2
        String outputFileNamePrefix
        String reference
        Boolean adapterTrim = true,
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
    "hg38": "/.mounts/labs/gsiprojects/gsi/Dragen/reference/hg38fa.v9"
    }


    String bwaMem_ref = dragenRef_by_genome [ reference ]
    call runDragen  { 
                input: 
                read1 = fastqR1,
                read2 = fastqR2,
                bwaRef = bwaMem_ref,
                adapterTrim = adapterTrim,
                prefix = outputFileNamePrefix,
                rgInfoString = rgInfoString
        }    
    }

    meta {
        author: "Lawrence Heisler"
        email: "lheisler@oicr.on.ca"
        description: "This workflow will align a fastq pair to the reference seqeunce using Illumina Dragen.  Adapter trimming is optional.  The bam file will be sorted and indexed"
        dependencies: [
        {
            name: "dragen/",
            url: ""
        }
      ]
    }

    output {
        File bam = runDragen.bam
        File bamIndex = runDragen.bamIndex
    }
}

task runDragen {
    input {
        File read1
        File read2
        String dragenRef
        String prefix
        Boolean adapterTrim
        String Adapter1File=/.mounts/labs/gsiprojects/gsi/Dragen/resources/ADAPTER1
        String Adapter2File=/.mounts/labs/gsiprojects/gsi/Dragen/resources/ADAPTER2
        String rgInfoString
        Int jobMemory = 500
        Int timeout = 96
    }

    parameter_meta {
        read1: "Fastq file for read 1"
        read2: "Fastq file for read 2"
        dragenRef: "The reference genome to align the sample with by Dragen"
        adapterTrim: "True/False for adapter trimming"
        Adapter1File: "Adapters to be trimmed from Read1"
        Adapter2File: "Adapters to be trimmed from Read2"
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
       --output-directory $OUT_DIR \
       --output-file-prefix ~{prefix} \
       --read-trimmers adapter \
       --trim-adapter-read1 /.mounts/labs/gsiprojects/gsi/Dragen/ADAPTER1 \
       --trim-adapter-read2 /.mounts/labs/gsiprojects/gsi/Dragen/ADAPTER2 \
       --trim-min-length 1 \
       --enable-bam-indexing true \
       --enable-sort true \
       --enable-duplicate-marking false
    >>>

    runtime {
        memory:  "~{jobMemory} GB"
        timeout: "~{timeout}"
    }  
    
    output {
        File outputBam = "~{resultBam}"
        File outputBamIndex
		WHAT OTHER FILES SHOULD BE PROVISIONED OUT
    }

    meta {
        output_meta: {
            outputBam: "output bam aligned to genome"
            outputBamIndex: "index for the aligned bam"
        }
    }

}



}
