get_click_coords=function(){
  coords <- locator(type="l")
  return(coords)
}

#' Convert click coordinates to rectangles.
#' This function assumes that click coordinates are given in pairs
#' The first click is the lower left, and the second click is the top right
#' 
#' @param coords The coordinates object received from pairs of clicks identifying rectangles.
#' @param chromo Current chromosome for the rectangles (for making a bedpe file of rectangles).
#' @param mini Minimum coordinate in the figure (important for the conversion from pixels to genomic coordinates)
#' @param maxi Maximum coordinate in the figure
#' @return A bedpe data frame with rectangles (first bedpe entry is the y axis, second bedpe entry is the x axis)
get_rectangularized_coords=function(coords,chromo,mini,maxi){
  print('get_rectangularized_coords')
  #need to have some coordinates
  stopifnot(length(coords[['x']])!=0)
  stopifnot(length(coords[['y']])!=0)
  #should have an even number of entries for rectangularizations
  stopifnot(length(coords[['x']])%%2==0)
  stopifnot(length(coords[['y']])%%2==0)
  #maximum should be larger than minimum
  stopifnot(maxi>mini)
  
  rectangles=data.frame(chr1=character(),
                        start1=numeric(),
                        end1=numeric(),
                        chr2=character(),
                        start2=numeric(),
                        end2=numeric(),
                        stringsAsFactors=FALSE)
  
  gendist=maxi-mini #genomic distance covered
  print('gendist')
  print(gendist)
  xst=coords[['x']][1]
  xend=coords[['x']][2]
  yst=coords[['y']][1]
  yend=coords[['y']][2]
  xlen=xend-xst
  ylen=yend-yst
  n=length(coords[['x']])/2
  for (r in c(2:n)){
    r_idx=2*(r-1)+1
    rectangle=data.frame(chr1=chromo,
                         start1=coords[['y']][r_idx],
                         end1=coords[['y']][r_idx+1],
                         chr2=chromo,
                         start2=coords[['x']][r_idx],
                         end2=coords[['x']][r_idx+1])
    rectangles=rbind(rectangles,rectangle)
  }
  rectangles[,'start2']=floor(((rectangles[,'start2']-xst)/xlen)*gendist+mini)
  rectangles[,'end2']=ceiling(((rectangles[,'end2']-xst)/xlen)*gendist+mini)
  rectangles[,'start1']=floor(((rectangles[,'start1']-yst)/ylen)*gendist+mini)
  rectangles[,'end1']=ceiling(((rectangles[,'end1']-yst)/ylen)*gendist+mini)
  return(rectangles)
}

#Expects the bedpe file to have 8 columns, that is, to have a score
read_bedpe=function(bedpefile){
  print('read_bedpe')
  BEDPE_ENTRIES_N=8
  bedpe=read.table(bedpefile)
  stopifnot(dim(bedpe)[2]>=BEDPE_ENTRIES_N)
  bedpe=data.frame(chr1=bedpe[,1],
                   start1=bedpe[,2],
                   end1=bedpe[,3],
                   chr2=bedpe[,4],
                   start2=bedpe[,5],
                   end2=bedpe[,6],
                   name=bedpe[,7],
                   value=bedpe[,8])
  return(bedpe)
}

read_bed=function(bedfile,featurename){
  bed=read.table(bedfile)
  bed=data.frame(chr=bed[,1],start=bed[,2],end=bed[,3],featurename)
  return(bed)
}

annotate_bedpe_with_widths=function(bedpe){
  print('annotate widths')
  bedpe=cbind(bedpe,w1=abs(as.numeric(bedpe$start1)-as.numeric(bedpe$end1)),
              w2=abs(as.numeric(bedpe$start2)-as.numeric(bedpe$end2)),
              mid1=(as.numeric(bedpe$start1)+as.numeric(bedpe$end1))/2,
              mid2=(as.numeric(bedpe$start2)+as.numeric(bedpe$end2))/2)
  return(bedpe)
}

plot_bedpe_scores=function(bedpe,midpoint,chromo,mini,maxi,
                           xbeds=data.frame(chr=character(),start=numeric(),
                                            end=numeric(),name=character()),
                           ybeds=data.frame(chr=character(),start=numeric(),
                                            end=numeric(),name=character())){
  require("ggplot2",lib.loc="/srv/scratch/oursu/code/")
  require("gridExtra",lib.loc="/srv/scratch/oursu/code/")
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
  }
  
  #mini=min(c(bedpe$start1,bedpe$start2))
  #maxi=max(c(bedpe$end1,bedpe$end2))
  MIN=-3
  MAX=3  
  
  if (dim(xbeds)[1]==0 && (dim(ybeds)[1]==0)){
    print(ggplot(bedpe, aes(x=mid1, y=mid2, fill = log(value,base=2))) + 
            geom_tile(aes(width = w2, height=w1))+
            xlab('X coordinate')+ 
            ylab('Y coordinate')+
            geom_vline(xintercept = mini)+
            geom_vline(xintercept = maxi)+
            geom_hline(yintercept = mini)+
            geom_hline(yintercept = maxi)+
            xlim(mini,maxi)+
            ylim(mini,maxi)+
            theme_bw()+ 
            scale_fill_gradient2(low="blue", midpoint=midpoint,high="red",limits=c(MIN,MAX)))
  }
  if ((dim(xbeds)[1]>0) && (dim(ybeds)[1]>0)){
    print(chromo)
    print(as.character(unique(xbeds$chr)))
    chromox=which(as.character(xbeds[,1])==as.character(chromo))
    chromoy=which(as.character(ybeds[,1])==as.character(chromo))
    stopifnot(length(chromox)>0)
    stopifnot(length(chromoy)>0)
    xbeds=xbeds[chromox,]
    ybeds=ybeds[chromoy,]
    print('1')
    #I'll assume that in this case the user has inputted a good bed file
    bottom_start=ggplot(xbeds,aes(name,start,ymin=start,ymax=end,colour=name))+geom_linerange(size=10)+theme_bw()+xlab("")+ylab("")+coord_flip()
    bottom=bottom_start+guides(colour=FALSE)+theme(axis.ticks=element_blank(), axis.text.y=element_blank())+ylim(c(mini,maxi))
    #side plot
    print('2')
    theheatmap=ggplot(bedpe, aes(x=mid1, y=mid2, fill = log(value,base=2))) + 
      geom_tile(aes(width = w2, height=w1))+
      xlab('X coordinate')+ 
      ylab('Y coordinate')+
      geom_vline(xintercept = mini)+
      geom_vline(xintercept = maxi)+
      geom_hline(yintercept = mini)+
      geom_hline(yintercept = maxi)+
      xlim(mini,maxi)+
      ylim(mini,maxi)+
      guides(fill=FALSE)+
      theme_bw()+
      theme(axis.ticks= element_blank(),axis.text.y = element_blank())+
      scale_fill_gradient2(low="blue", midpoint=midpoint,high="red",limits=c(MIN,MAX))
    print('3')
    side=ggplot(ybeds,aes(name,start,ymin=start,ymax=end,colour=name))+geom_linerange(size=10)+theme_bw()+xlab("")+ylab("")+guides(colour=FALSE)+ylim(c(mini,maxi))
    print('4')
    grid.arrange(side,theheatmap, g_legend(bottom_start),bottom,ncol=2, nrow=2,
                 heights=c(5,1), widths=c(1,5))    
    print('5')
  }
}

bedpe_to_rectangles=function(bedpefile,chromo,mini,maxi,centervalue){
  bedpe=read_bedpe(bedpefile)
  bedpe=annotate_bedpe_with_widths(bedpe)
  plot_bedpe_scores(bedpe,centervalue,chromo,mini,maxi)
  rectangles=get_rectangularized_coords(get_click_coords(),chromo,mini,maxi)
  return(rectangles)
}

plot_bedpe_to_pdf=function(bedpefile,chromo,mini,maxi,centervalue,out,
                           xbeds=data.frame(chr=character(),start=numeric(),
                                            end=numeric(),name=character()),
                           ybeds=data.frame(chr=character(),start=numeric(),
                                            end=numeric(),name=character())){
  print(dim(xbeds))
  bedpe=read_bedpe(bedpefile)
  bedpe=annotate_bedpe_with_widths(bedpe)
  pdf(paste(out,'.pdf',sep=''),height=20,width=23)
  plot_bedpe_scores(bedpe,centervalue,chromo,mini,maxi,xbeds,ybeds)
  dev.off()
}

#' Function to show you the whole chromosome, one window at a time.
#' It assumes that you started from a bedpe format file, 
#' ran divide_bedpe_in_windows.sh to divide the file into genomic windows 
#' window-level visualization, and are now ready to annotate each window.
#' 
#' @param bedpefile A template for the bedpe files that will be read. 
#' Substitute WINDOW for the window size, and MINI for the beginning of the current window.
#' As example is: /srv/scratch/oursu/3Dgenome/results/processed_data/HiC/counts/intra/GM12878_combined/5kb/chr21/windowData/wWINDOW/chr21_5kb.RAWobserved.norm_SQRTVC.obsOverExp_SQRTVCexpected.wMINI.gz
#' @param chromo Name of the chromosome. This is for making the output file.
#' @param mini Starting position of the genomic windows to show. 
#' @param maxi Ending position of the genomic windows to show
#' @param w Size of windows that will be shown to you.
#' @param outdir The name of the output directory where the output data will be stored
#' @param pref The prefix of the files that will be created. Useful for when the same method 
#' will be run for finding different types of elements, e.g. TADs, loops, etc., then 
#' you can just set the prefix to "TADs", "Loops", etc.
#' @examples bedpe_to_rectangles_byChromosome('/srv/scratch/oursu/3Dgenome/results/processed_data/HiC/counts/intra/GM12878_combined/5kb/chr21/windowData/wWINDOW/chr21_5kb.RAWobserved.norm_SQRTVC.obsOverExp_SQRTVCexpected.wMINI.gz','chr21',30000000,39000000,1000000,outdir,outpref)
bedpe_to_rectangles_byChromosome=function(bedpefile,chromo,mini,maxi,w,outdir,pref){
  for (cur_mini in seq(from=mini,to=maxi,by=w)){
    print(cur_mini)
    cur_bedpefile=gsub('WINDOW',format(w,scientific=FALSE),gsub('MINI',format(cur_mini,scientific=FALSE),bedpefile))
    print(cur_bedpefile)
    rectangles=bedpe_to_rectangles(cur_bedpefile,chromo,cur_mini,cur_mini+w,0)
    outfile=paste(outdir,pref,basename(cur_bedpefile),'.bedpe',sep='')
    write.table(rectangles,file=outfile,
                sep='\t',quote=F,row.names=F,col.names=F)
    system(paste('zcat -f ',outfile,' | gzip > ',outfile,'.gz',sep=''))
    system(paste('rm ',outfile))
  }
}




