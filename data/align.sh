#!/bin/bash

# The static build of 2.4 segfaults on various >4 linux kernels (debian / ubuntu at least). Local compile from source works though:

mkdir -p '2.4/'
/data/disco-refs/2.4/STAR/source/STAR \
	--genomeDir '/data/disco-refs/2.4/index' \
	--runThreadN 7 \
	--outSAMtype BAM SortedByCoordinate \
	--alignEndsType Local \
	--alignIntronMax 200000 \
	--alignMatesGapMax 200000 \
	--alignSJDBoverhangMin 10 \
	--chimJunctionOverhangMin 12 \
	--chimSegmentMin 12 \
	--outFilterIntronMotifs None \
	--outSAMstrandField intronMotif \
	--outFileNamePrefix '2.4/' \
	--readFilesIn 'R1.fq' 'R2.fq'



# Expensive:
#	--twopass1readsN -1 \
#	--twopassMode Basic \


# From EGFRvIII manuscript:
#--outSAMtype BAM Unsorted
#--outSAMstrandField intronMotif
#--outFilterIntronMotifs None
#--alignIntronMin 20
#--alignIntronMax 200000 \
#--alignMatesGapMax 200000 \
#--alignSJoverhangMin 10
#--alignSJDBoverhangMin 1
#--chimSegmentMin 12
#--chimJunctionOverhangMin 12
#--chimOutType WithinBAM SeparateSAMold



mkdir -p '2.7/'
/data/disco-refs/2.7/STAR/bin/Linux_x86_64_static/STAR \
	--genomeDir '/data/disco-refs/2.7/index' \
	--runThreadN 7 \
	--outSAMtype BAM Unsorted \
	--outSAMstrandField intronMotif \
	--outFilterIntronMotifs None \
	--alignIntronMin 20 \
	--alignIntronMax 200000 \
	--alignMatesGapMax 200000 \
	--alignSJoverhangMin 10 \
	--alignSJDBoverhangMin 1 \
	--chimSegmentMin 12 \
	--chimJunctionOverhangMin 12 \
	--chimOutType WithinBAM SeparateSAMold \
	--outFileNamePrefix '2.7/' \
	--readFilesIn 'R1.fq' 'R2.fq'


	
	#--outSAMtype BAM SortedByCoordinate \
	#--alignEndsType Local \
	#--alignIntronMax 200000 \
	#--alignMatesGapMax 200000 \
	#--alignSJDBoverhangMin 10 \
	#--chimJunctionOverhangMin 12 \
	#--chimSegmentMin 12 \
	#--outFilterIntronMotifs None \
	#--outSAMstrandField intronMotif \
	#--twopass1readsN -1 \
	#--twopassMode Basic \
	#--outFileNamePrefix '2.7/' \
	#--readFilesIn 'protocol/fastq/C6U5LANXX_7046-004-132_TCCGCGAA-GTACTGAC_L004_R1.fastq.gz,protocol/fastq/C6VVJANXX_7046-004-132_TCCGCGAA-GTACTGAC_L007_R1.fastq.gz,protocol/fastq/C7192ANXX_7046-004-132_TCCGCGAA-GTACTGAC_L005_R1.fastq.gz,protocol/fastq/C7192ANXX_7046-004-132_TCCGCGAA-GTACTGAC_L007_R1.fastq.gz,protocol/fastq/C75EJANXX_7046-004-132_TCCGCGAA-GTACTGAC_L005_R1.fastq.gz,protocol/fastq/C75EJANXX_7046-004-132_TCCGCGAA-GTACTGAC_L006_R1.fastq.gz,protocol/fastq/C75JHANXX_7046-004-132_TCCGCGAA-GTACTGAC_L001_R1.fastq.gz,protocol/fastq/C75JHANXX_7046-004-132_TCCGCGAA-GTACTGAC_L002_R1.fastq.gz,protocol/fastq/C75JHANXX_7046-004-132_TCCGCGAA-GTACTGAC_L003_R1.fastq.gz,protocol/fastq/C75JHANXX_7046-004-132_TCCGCGAA-GTACTGAC_L004_R1.fastq.gz,protocol/fastq/C75JHANXX_7046-004-132_TCCGCGAA-GTACTGAC_L005_R1.fastq.gz,protocol/fastq/C77NWANXX_7046-004-132_TCCGCGAA-GTACTGAC_L005_R1.fastq.gz,protocol/fastq/C77W2ANXX_7046-004-132_TCCGCGAA-GTACTGAC_L005_R1.fastq.gz,protocol/fastq/C77W2ANXX_7046-004-132_TCCGCGAA-GTACTGAC_L006_R1.fastq.gz,protocol/fastq/C77W2ANXX_7046-004-132_TCCGCGAA-GTACTGAC_L007_R1.fastq.gz,protocol/fastq/C77W2ANXX_7046-004-132_TCCGCGAA-GTACTGAC_L008_R1.fastq.gz' 'protocol/fastq/C6U5LANXX_7046-004-132_TCCGCGAA-GTACTGAC_L004_R2.fastq.gz,protocol/fastq/C6VVJANXX_7046-004-132_TCCGCGAA-GTACTGAC_L007_R2.fastq.gz,protocol/fastq/C7192ANXX_7046-004-132_TCCGCGAA-GTACTGAC_L005_R2.fastq.gz,protocol/fastq/C7192ANXX_7046-004-132_TCCGCGAA-GTACTGAC_L007_R2.fastq.gz,protocol/fastq/C75EJANXX_7046-004-132_TCCGCGAA-GTACTGAC_L005_R2.fastq.gz,protocol/fastq/C75EJANXX_7046-004-132_TCCGCGAA-GTACTGAC_L006_R2.fastq.gz,protocol/fastq/C75JHANXX_7046-004-132_TCCGCGAA-GTACTGAC_L001_R2.fastq.gz,protocol/fastq/C75JHANXX_7046-004-132_TCCGCGAA-GTACTGAC_L002_R2.fastq.gz,protocol/fastq/C75JHANXX_7046-004-132_TCCGCGAA-GTACTGAC_L003_R2.fastq.gz,protocol/fastq/C75JHANXX_7046-004-132_TCCGCGAA-GTACTGAC_L004_R2.fastq.gz,protocol/fastq/C75JHANXX_7046-004-132_TCCGCGAA-GTACTGAC_L005_R2.fastq.gz,protocol/fastq/C77NWANXX_7046-004-132_TCCGCGAA-GTACTGAC_L005_R2.fastq.gz,protocol/fastq/C77W2ANXX_7046-004-132_TCCGCGAA-GTACTGAC_L005_R2.fastq.gz,protocol/fastq/C77W2ANXX_7046-004-132_TCCGCGAA-GTACTGAC_L006_R2.fastq.gz,protocol/fastq/C77W2ANXX_7046-004-132_TCCGCGAA-GTACTGAC_L007_R2.fastq.gz,protocol/fastq/C77W2ANXX_7046-004-132_TCCGCGAA-GTACTGAC_L008_R2.fastq.gz'

##	--readFilesCommand zcat \


