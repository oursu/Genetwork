enhfile=$1
genefile=$2
genename=$3
out=$4
bedpe=$5
w=$6

enhfile=/srv/scratch/oursu/cas9-3d/results/2015-09-28/K562.SE_removeBlacklist.bed
genefile=/srv/scratch/oursu/cas9-3d/results/2015-09-28/gencode.v19.annotation.PC.lincRNA.TSS.bed.gz
genename=TP53
out=/srv/scratch/oursu/cas9-3d/results/2015-09-28/test.gz
bedpe=/srv/scratch/oursu/3Dgenome/results/processed_data/HiC/counts/intra/GM12878_combined/5kb/CHROMO/CHROMO_5kb.RAWobserved.norm_SQRTVC.obsOverExp_SQRTVCexpected.bedpe.gz

#1. Make bed file with TSS +-500bp
genecoords=${out}.${genename}.bed
zcat -f ${genefile} | grep -w ${genename} | awk '{s=$2-500}{e=$3+500}{print "chr"$1"\t"s"\t"e"\t"$4}' | head -n1 > ${genecoords}
zcat -f ${genefile} | grep -w ${genename} | awk -v w=${w}'{s=$2-w}{e=$3+w}{print "chr"$1"\t"s"\t"e"\t"$4}' | head -n1 > ${genecoords}_window.bed


#2. Get interactions within a reasonable window around the gene TSS
chromo=$(zcat -f ${genecoords} | cut -f1)
bedpe_chr=$(echo ${bedpe} | sed 's/CHROMO/'${chromo}'/g')
bedtools pairtobed -type both -a ${bedpe_chr} -b ${genecoords}_window.bed | \
cut -1-8 | gzip > ${out}