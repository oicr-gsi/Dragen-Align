#!/bin/bash
set -o nounset
set -o pipefail

cd $1

module load samtools/1.9 2>/dev/null

ls | sort

find -name *.zip -xtype f -exec unzip -l {} \; | awk '{print $4}' | sed -e '/^\(Name\|----\)$/d' -e '/^[[:space:]]*$/d' | sort

find -name *.bam -xtype f -exec samtools view -H {} \; | grep '^@RG' | sort

find -name *.bam -xtype f -exec samtools flagstat {} \; | sort

find -name *.bam -xtype f -exec /bin/bash -c "samtools view {} | md5sum" \; | sort

