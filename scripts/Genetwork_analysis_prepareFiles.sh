

#For TADs and large-scale structures
w=20000000
for res in 25kb;
do
 for celltype in GM12878_combined;
 do 
  for chromosome in {1..22};
  do
   chromo=chr${chromosome}
   DI=/srv/scratch/oursu/3Dgenome/results/processed_data/HiC/counts/
   INFILE=${DI}/intra/${celltype}/${res}/${chromo}/${chromo}_${res}.RAWobserved.norm_SQRTVC.obsOverExp_SQRTVCexpected.bedpe.gz
   OUT=/srv/scratch/oursu/3Dgenome/results/Genetwork_analysis/HiC/counts/intra/${res}/${chromo}/${chromo}_${res}.w${w}
   mkdir -p $(dirname ${OUT})
   s=${OUT}_script.sh
   echo "source /srv/scratch/oursu/code/genome_utils/3Dutils/bashrc_3D" > ${s}
   echo "divide_bedpe_in_windows.sh ${INFILE} ${OUT} ${w} \${chrSizes}" >> ${s}
   cat ${s}
   qsub -l h_vmem=20G -o ${s}.o -e ${s}.e ${s}
  done
 done
done

