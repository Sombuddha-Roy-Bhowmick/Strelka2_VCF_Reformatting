# Strelka2 VCF Reformatting
Reformatting of the somatic and germline VCFs produced by Strelka2 variant caller.
Since the allele frequency (AF) is not directly available in the VCF output. So, to extract such values from the Strelka2 VCF, the python scripts are available here to do so.
For full instructions on how to interpret results, please visit the official github page of Strelka2. (https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md)

For Somatic SNVs:
refCounts = Value of FORMAT column $REF + “U” (e.g. if REF="A" then use the value in FOMRAT/AU)
altCounts = Value of FORMAT column $ALT + “U” (e.g. if ALT="T" then use the value in FOMRAT/TU)
tier1RefCounts = First comma-delimited value from $refCounts
tier1AltCounts = First comma-delimited value from $altCounts
Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)

For Somatic indels:
tier1RefCounts = First comma-delimited value from FORMAT/TAR
tier1AltCounts = First comma-delimited value from FORMAT/TIR
Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)

In addition, the reformatting script for Strelka2 Germline calling is also provided here, along with an additional script for calculating Ti/Tv ratio, using multithreading option is also available here.

