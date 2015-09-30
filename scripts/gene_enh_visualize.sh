enhfile=$1
genefile=$2
genename=$3
out=$4
bedpe=$5
w=$6
geneexpression=$7

#1. Make bed file with TSS +-500bp
genecoords=${out}.${genename}.bed
out=${out}.${genename}
s=${out}.sh
zcat -f ${genefile} | grep -w ${genename} | awk '{s=$2-500}{e=$3+500}{print "chr"$1"\t"s"\t"e"\t"$4}' | head -n1 > ${genecoords}
zcat -f ${genefile} | grep -w ${genename} | awk -v w=${w} '{s=$2-w}{e=$3+w}{print "chr"$1"\t"s"\t"e"\t"$4}' | head -n1 > ${genecoords}_window.bed

#also plot gene expression profile
ensembl=$(zcat -f ${genecoords} | sed 's/;/\t/g' | sed 's/[.]/\t/g' | cut -f5)
echo ${ensembl}
echo "zcat -f ${geneexpression} | head -n1 | cut --complement -f1 > ${out}.expression" >> ${s}
echo "zcat -f ${geneexpression} | grep -w ${ensembl} | cut --complement -f1 >> ${out}.expression" >> ${s}
echo "Rscript /srv/scratch/oursu/code/Genetwork/scripts/plot_2col_barchart.R ${out}.expression ${out}.expression.pdf" >> ${s}

#2. Get interactions within a reasonable window around the gene TSS
chromo=$(zcat -f ${genecoords} | cut -f1)
bedpe_chr=$(echo ${bedpe} | sed 's/CHROMO/'${chromo}'/g')
mini=$(zcat -f ${genecoords}_window.bed | cut -f2)
maxi=$(zcat -f ${genecoords}_window.bed | cut -f3)
#prepare xbeds,ybeds
echo "zcat -f ${enhfile} | awk '{print \$1\"\t\"\$2\"\t\"\$3\"\tenhancers\"}' | gzip > ${out}.enhancers.gz " >> ${s}
echo "bedtools pairtobed -type both -a ${bedpe_chr} -b ${genecoords}_window.bed | cut -f1-8 | gzip > ${out}" >> ${s}

#3. Plot!
RCODE=/srv/scratch/oursu/code/Genetwork/Genetwork/R/
echo "Rscript ${RCODE}/plot_bedpe_pdf.R ${out} ${chromo} ${mini} ${maxi} ${out}_plot ${out}.enhancers.gz,${genecoords} ${out}.enhancers.gz,${genecoords}" >> ${s}
echo "rm ${out}.enhancers.gz ${genecoords} ${genecoords}_window.bed ${out}.expression" >> ${s}
chmod 755 ${s}
qsub -o ${s}.o -e ${s}.e ${s}






