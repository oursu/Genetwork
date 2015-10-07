
require(Genetwork)
chromosome="chr1"
start=0
end=80000000 
wsize=20000000
outdir="/srv/scratch/oursu/3Dgenome/results/Genetwork_analysis/clicks/HiC/counts/intra/25kb/"
outpref="TADs_random"
res="25kb"
bedpefile=paste('/srv/scratch/oursu/3Dgenome/results/Genetwork_analysis/HiC/counts/intra/',res,'/',chromosome,'/',chromosome,'_',res,'.wWINDOW.wMINI.gz',sep='')
'_',res,'.RAWobserved.norm_SQRTVC.obsOverExp_SQRTVCexpected.wMINI.gz',sep='')


dev.off()
plot.new()
bedpe_to_rectangles_byChromosome(bedpefile,chromosome,start,end,wsize,outdir,outpref)