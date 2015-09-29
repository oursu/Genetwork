#this is for using LA files.
#Not used at the moment
divide_bedpe_by_exp=function(bedpe,w,divisionfile,colToDivide,newColName,pseudocount){
  divisiondata=read.table(divisionfile)
  bedpe_dist=abs(bedpe[,'start1']-bedpe[,'start2'])/w+1
  dvals=divisiondata[,1][bedpe_dist]
  bedpe=cbind(bedpe,newcol=(as.numeric(bedpe[,colToDivide]+pseudocount)/(dvals+pseudocount)))
  colnames(bedpe)[ncol(bedpe)]=newColName
  return(bedpe)
}
bedpe_to_rectangles_old=function(bedpefile,chromo,expectedFile,mini,maxi,w){
  bedpe=read_bedpe(bedpefile)
  bedpe=annotate_bedpe_with_widths(bedpe)
  bedpe=divide_bedpe_by_exp(bedpe,w,expectedFile,
                            'v','obsOverExp',0.1)
  plot_bedpe_by_col(bedpe,'obsOverExp',0,chromo,mini,maxi)
}

#testing bedpe_to_rectangles('/Users/oursu/testbedpe.gz','chr21',9830000,10205000,0)


gb2way=function(m,xbeds,ybeds,xbedName,ybedName,titleName,tickScale,xlimits1,xlimits2,ylimits1,ylimits2){
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
  }
  xbeds$start=xbeds$start/tickScale
  xbeds$end=xbeds$end/tickScale
  ybeds$start=ybeds$start/tickScale
  ybeds$end=ybeds$end/tickScale
  require(reshape2)
  require(ggplot2)
  require(gridExtra)
  #bottom plot
  bottom_start=ggplot(xbeds,aes(name,start,ymin=start,ymax=end,colour=name))+geom_linerange(size=10)+theme_bw()+xlab("")+ylab("")+coord_flip()
  bottom=bottom_start+guides(colour=FALSE)+theme(axis.ticks=element_blank(), axis.text.y=element_blank())+ylim(c(xlimits1,xlimits2))
  #side plot
  side=ggplot(ybeds,aes(name,start,ymin=start,ymax=end,colour=name))+geom_linerange(size=10)+theme_bw()+xlab("")+ylab("")+guides(colour=FALSE)+ylim(c(ylimits1,ylimits2))
  legend_bed=g_legend(bottom_start)
  #setup sizes for the matrix representation
  xbed=xbeds[which(as.character(xbeds$name)==xbedName),]
  xbedsizes=abs(xbed$start-xbed$end)
  print(xbedsizes)
  ybed=ybeds[which(as.character(ybeds$name)==ybedName),]
  ybedsizes=abs(ybed$start-ybed$end)
  if (nrow(ybed)!=nrow(m)){
    print('ybedName does not have the same number of items as m for y axis')
  }
  if (nrow(xbed)!=ncol(m)){
    print('xbedName does not have the same number of items as m for x axis')
  }
  rownames(m)=as.character((ybed$start+ybed$end)/2)
  colnames(m)=as.character((xbed$start+xbed$end)/2)
  mmelt=melt(m)
  #matrix picture
  w=rep(xbedsizes,each=length(ybedsizes))
  h=rep(ybedsizes,times=length(xbedsizes))
  mmelt_heatmap=ggplot(mmelt,aes(Var2,Var1,fill=value))+theme_bw()+geom_tile(width=w,height=h)+guides(fill=FALSE)+ylab("")+xlab("")+theme(axis.ticks= element_blank(),axis.text.y = element_blank())+xlim(c(xlimits1,xlimits2))+ylim(c(ylimits1,ylimits2))
  grid.arrange(side,mmelt_heatmap, legend_bed,bottom,ncol=2, nrow=2,main =paste(titleName,'\nScale=',as.character(tickScale/1000),' kb',sep=''))#,heights=c(5,1), widths=c(5,1))
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

bedpe_to_rectangles_with_beds=function(bedpefile,chromo,mini,maxi,centervalue,xbeds,ybeds){
  bedpe=read_bedpe(bedpefile)
  bedpe=annotate_bedpe_with_widths(bedpe)
  plot_bedpe_score_with_beds(bedpe,centervalue,chromo,mini,maxi,xbeds,ybeds)
}

