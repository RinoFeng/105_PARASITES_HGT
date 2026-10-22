args<-commandArgs(TRUE);o<-normalizePath(args[1]);library(grid)
d<-read.delim(file.path(o,'data/FigureS2_species_sources.tsv'),check.names=FALSE)
dom<-c('Bacteria','Eukaryota','Archaea','Viruses');cols<-c('#4DBBD5','#00A087','#E64B35','#8266B0','#8C8C8C');names(cols)<-dom
life<-c(Intracellular='#DC0000B2',Extracellular='#3C5488B2','Free-living'='#00A087B2',Unknown='#D9D9D9',Conflicting='#555555');grp<-c('Microbial eukaryotes'='#2E86AB','Metazoan parasites'='#F18F01',Unknown='#D9D9D9')
stopifnot(nrow(d)==105,sum(as.matrix(d[,dom]))==8989,sum(d$Viruses)==54)
textg<-function(t,x,y,size=8,...)grid.text(t,x,y,gp=gpar(fontsize=size,fontfamily='Arial'),...)
panel<-function(page){
 z<-d[d$page==page,];n<-nrow(z);x0<-.445;cw<-.218;px<-.750;pw<-.158;nt<-.710;mx<-1500;top<-.89;bot<-.025;rh<-(top-bot)/53
 textg(if(page==1)'A' else 'B',.004,.985,10,just=c('left','top'));textg('Candidate count',x0+cw/2,.975,8);textg('N',nt,.975,8);textg('Source proportion',px+pw/2,.975,8);textg('L',.945,.975,7.8);textg('G',.980,.975,7.8)
 for(v in c(0,750,1500)){x<-x0+cw*v/mx;grid.lines(c(x,x),c(bot,top),gp=gpar(col='#dddddd',lwd=.4));textg(v,x,.932,7.6)}
 for(v in c(0,50,100))textg(paste0(v,'%'),px+pw*v/100,.932,7.6)
 for(i in seq_len(n)){
  y<-top-(i-.5)*rh;v<-as.numeric(z[i,dom]);N<-sum(v)
  if(i%%2==0)grid.rect(.5,y,1,rh,gp=gpar(fill='#F4F4F4',col=NA))
  grid.text(z$species[i],.432,y,just='right',gp=gpar(fontsize=7.8,fontfamily='Arial',fontface='italic'))
  cur<-x0;cp<-px
  for(j in seq_along(v)){ww<-cw*v[j]/mx;if(ww>0)grid.rect(cur+ww/2,y,ww,rh*.70,gp=gpar(fill=cols[j],col=NA));cur<-cur+ww;ww<-if(N>0)pw*v[j]/N else 0;if(ww>0)grid.rect(cp+ww/2,y,ww,rh*.70,gp=gpar(fill=cols[j],col=NA));cp<-cp+ww}
  textg(N,nt+.022,y,7.8,just='right');if(N==0)textg('—',px+pw/2,y,7.8)
  lc<-life[z$lifestyle[i]];if(is.na(lc))lc<-life['Unknown'];gc<-grp[z$group_label[i]];if(is.na(gc))gc<-grp['Unknown']
  grid.rect(.944,y,.014,rh*.70,gp=gpar(fill=lc,col=NA));grid.rect(.980,y,.014,rh*.70,gp=gpar(fill=gc,col=NA))
 }
}
plotone<-function(){
 grid.newpage();grid.rect(gp=gpar(fill='white',col=NA));textg('Source profiles across 105 sampled genomes',.002,.988,10,just=c('left','top'));textg('8,989 candidates; 54 current viral-source annotations',.995,.988,8.5,just=c('right','top'))
 pushViewport(viewport(.242,.513,.482,.895));panel(1);popViewport();pushViewport(viewport(.751,.513,.482,.895));panel(2);popViewport()
 grid.lines(c(.497,.497),c(.065,.967),gp=gpar(col='#BBBBBB',lwd=.5))
 xx<-c(.012,.120,.235,.337,.420)
 for(j in seq_along(dom)){grid.rect(xx[j],.045,.007,.013,gp=gpar(fill=cols[j],col=NA));textg(dom[j],xx[j]+.007,.045,7.8,just='left')}
 textg('L: Lifestyle:',.585,.045,7.8,just='left');xx<-c(.664,.754,.844,.934)
 for(j in 1:4){grid.rect(xx[j],.045,.007,.013,gp=gpar(fill=life[j],col=NA));textg(names(life)[j],xx[j]+.007,.045,7.8,just='left')}
 textg('G: Group:',.012,.015,7.8,just='left');for(j in 1:2){x<-c(.073,.245)[j];grid.rect(x,.015,.007,.013,gp=gpar(fill=grp[j],col=NA));textg(names(grp)[j],x+.007,.015,7.8,just='left')}
 textg('A and B continue the same species order; identical count scales',.995,.015,7.8,just='right')
}
f<-file.path(o,'figures/FigureS2_single_page');cairo_pdf(paste0(f,'.pdf'),width=270/25.4,height=150/25.4,family='Arial');plotone();dev.off();png(paste0(f,'.png'),width=270,height=150,units='mm',res=450,type='cairo');plotone();dev.off();svglite::svglite(paste0(f,'.svg'),width=270/25.4,height=150/25.4,fix_text_size=FALSE);plotone();dev.off()
