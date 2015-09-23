get_click_coords=function(){
  coords <- locator(type="l")
  return(coords)
}

get_coords_heatmap=function(df){
  require(pheatmap)
  pheatmap(as.matrix(df),cluster_cols=F,cluster_rows=F,scale='none')
  return(get_click_coords())
}

#' Convert click coordinates to rectangles.
#' This function assumes that click coordinates are given in pairs
#' The first click is the top left, and the second click is the bottom right
#' 
#' @param coords The coordinates object received from pairs of clicks identifying rectangles.
#' @param chromo Current chromosome for the rectangles (for making a bedpe file of rectangles).
#' @return A bedpe data frame with rectangles (first bedpe entry is the y axis, second bedpe entry is the x axis)
#' @examples
#' a=data.frame(a=c(1,2,3),b=c(2,3,4),d=c(3,4,5))
#' coords_out=get_coords_heatmap(a)
#' get_rectangularized_coords(coords_out,'chr1')
get_rectangularized_coords=function(coords,chromo,gendist,mini){
  #need to have some coordinates
  stopifnot(length(coords[['x']])!=0)
  stopifnot(length(coords[['y']])!=0)
  #should have an even number of entries for rectangularizations
  stopifnot(length(coords[['x']])%%2==0)
  stopifnot(length(coords[['y']])%%2==0)
  
  rectangles=data.frame(chr1=character(),
                        start1=numeric(),
                        end1=numeric(),
                        chr2=character(),
                        start2=numeric(),
                        end2=numeric(),
                        stringsAsFactors=FALSE)
  
  xst=coords[['x']][1]
  xend=coords[['x']][2]
  yst=coords[['y']][2]
  yend=coords[['y']][1]
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

read_bedpe=function(bedpefile){
  bedpe=read.table(bedpefile)
  bedpe=data.frame(chr1=bedpe[,1],
                   start1=bedpe[,2],
                   end1=bedpe[,3],
                   chr2=bedpe[,4],
                   start2=bedpe[,5],
                   end2=bedpe[,6],
                   v=bedpe[,7])
  return(bedpe)
}

annotate_bedpe_with_widths=function(bedpe){
  bedpe=cbind(bedpe,w1=abs(as.numeric(bedpe$start1)-as.numeric(bedpe$end1)),
              w2=abs(as.numeric(bedpe$start2)-as.numeric(bedpe$end2)),
              mid1=(as.numeric(bedpe$start1)+as.numeric(bedpe$end1))/2,
              mid2=(as.numeric(bedpe$start2)+as.numeric(bedpe$end2))/2)
  return(bedpe)
}

divide_bedpe_by_exp=function(bedpe,w,divisionfile,colToDivide,newColName,pseudocount){
  divisiondata=read.table(divisionfile)
  bedpe_dist=abs(bedpe[,'start1']-bedpe[,'start2'])/w+1
  dvals=divisiondata[,1][bedpe_dist]
  bedpe=cbind(bedpe,newcol=(as.numeric(bedpe[,colToDivide]+pseudocount)/(dvals+pseudocount)))
  colnames(bedpe)[ncol(bedpe)]=newColName
  return(bedpe)
}

plot_bedpe_by_col=function(bedpe,colName,midpoint,chromo,mini,maxi){
  require(ggplot2)
  #mini=min(c(bedpe$start1,bedpe$start2))
  #maxi=max(c(bedpe$end1,bedpe$end2))
  print(ggplot(bedpe, aes(x=mid1, y=mid2, fill = log(obsOverExp,base=2))) + 
          geom_tile(aes(width = w2, height=w1))+xlab('X coordinate')+ 
          ylab('Y coordinate')+geom_vline(xintercept = mini)+geom_vline(xintercept = maxi)+
          geom_hline(yintercept = mini)+geom_hline(yintercept = maxi)+xlim(mini,maxi)+theme_bw()+ ylim(mini,maxi)+
          scale_fill_gradient2(low="yellow", midpoint=midpoint,high="purple",limits=c(-10,10)))
  rectangles=get_rectangularized_coords(get_click_coords(),chromo,abs(maxi-mini),mini)
  return(rectangles)
}

bedpe_to_rectangles=function(infile,chromo,expectedFile,mini,maxi,w){
  bedpe=read_bedpe(infile)
  bedpe=annotate_bedpe_with_widths(bedpe)
  bedpe=divide_bedpe_by_exp(bedpe,w,expectedFile,
                            'v','obsOverExp',0.1)
  plot_bedpe_by_col(bedpe,'obsOverExp',0,chromo,mini,maxi)
}

analyze_files_by_chr=function(indir,chromo,expectedFile,largeW,outdir,smallW){
  fs=Sys.glob(paste(indir,'/',chromo,'/*window*gz',sep=''))
  for (infile in fs){
    wstart=as.numeric(gsub('.gz','',gsub('window.','',unlist(regmatches(infile,gregexpr("window.*.gz",infile))))))
    rectangles=bedpe_to_rectangles(infile,chromo,expectedFile,as.numeric(wstart),as.numeric(wstart+largeW),smallW)
    write.table(rectangles,file=paste(outdir,'/SelectedRectangles_',chromo,'_',round(smallW/1000),'kb.RAWobserved.bedpe.window.',
                                      wstart,'.bedpe.txt',sep=''),
                sep='\t',quote=F,col.names=F,row.names=F) 
  }
}

quick_start=function(){
  analyze_files_by_chr(indir='/Users/oursu/',
                       chromo='chr22',
                       expectedFile='/Users/oursu/chr22_5kb.RAWexpected',
                       1000000,
                       '/Users/oursu/',
                       5000)
}






