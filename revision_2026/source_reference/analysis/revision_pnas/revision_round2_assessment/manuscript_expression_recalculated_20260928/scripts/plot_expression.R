args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else '.')
library(grid)
read_ppm <- function(file) {
  con <- file(file, 'rb'); on.exit(close(con))
  stopifnot(readLines(con, n=1)=='P6')
  dims <- as.integer(strsplit(readLines(con,n=1),' ')[[1]])
  stopifnot(readLines(con,n=1)=='255')
  bytes <- as.integer(readBin(con,'raw',n=prod(dims)*3))
  colors <- rgb(bytes[seq(1,length(bytes),3)],bytes[seq(2,length(bytes),3)],bytes[seq(3,length(bytes),3)],maxColorValue=255)
  as.raster(matrix(colors,nrow=dims[2],ncol=dims[1],byrow=TRUE))
}
d <- read.delim(file.path(root,'data/Ld_signal_primary_reverse_gene.tsv'),check.names=FALSE,row.names=1)
meta <- read.delim(file.path(root,'data/Ld_sample_metadata.tsv'))
order <- c('SRR32312282','SRR32312281','SRR32312280','SRR32312279','SRR32312278','SRR32312277')
v <- log10(as.matrix(d[,order])); stopifnot(all(is.finite(v)))
palette <- colorRampPalette(c('#F7FBFF','#6BAED6','#08306B'))(256)
lo <- floor(min(v)*10)/10; hi <- ceiling(max(v)*10)/10
color_for <- function(x) palette[pmax(1,pmin(256,round((x-lo)/(hi-lo)*255)+1))]
panel <- function(standalone=FALSE) {
  grid.rect(gp=gpar(fill='white',col=NA))
  if(!standalone) grid.text('D',x=.018,y=.985,just=c('left','top'),gp=gpar(fontsize=13,fontface='bold'))
  grid.text('Nascent transcription (PRO-seq)',x=.51,y=.965,gp=gpar(fontsize=8,fontface='bold'))
  grid.text('Promastigote',x=.18,y=.912,gp=gpar(fontsize=7.8))
  grid.text('Axenic amastigote-like',x=.52,y=.912,gp=gpar(fontsize=7.8))
  x0 <- .035; heatwidth <- .66; y0 <- .11; heatheight <- .745
  cw <- heatwidth/6; rh <- heatheight/nrow(v)
  for(j in 1:6) {
    grid.rect(x=x0+(j-.5)*cw,y=.878,width=cw,height=.018,gp=gpar(fill=if(j<=3)'#1B9E77' else '#D95F02',col=NA))
    for(i in 1:nrow(v)) grid.rect(x=x0+(j-.5)*cw,y=y0+heatheight-(i-.5)*rh,width=cw,height=rh,gp=gpar(fill=color_for(v[i,j]),col='white',lwd=.35))
    grid.text(c('P1','P2','P3','A1','A2','A3')[j],x=x0+(j-.5)*cw,y=.085,gp=gpar(fontsize=7.8))
  }
  for(i in 1:nrow(v)) grid.text(rownames(v)[i],x=.71,y=y0+heatheight-(i-.5)*rh,just='left',gp=gpar(fontsize=7.8,fontface=if(rownames(v)[i]=='LDBPK_171590')'bold' else 'plain'))
  # Compact numerical legend; all genes exceed signal 10 in every library.
  yy <- .035
  for(k in 1:100)grid.rect(x=.10+(k-.5)*.40/100,y=yy,width=.40/100,height=.017,gp=gpar(fill=color_for(lo+(k-.5)*(hi-lo)/100),col=NA))
  grid.text(sprintf('%.1f',lo),x=.07,y=yy,gp=gpar(fontsize=7.8))
  grid.text(sprintf('%.1f',hi),x=.535,y=yy,gp=gpar(fontsize=7.8))
  grid.text('log10(signal)',x=.59,y=yy,just='left',gp=gpar(fontsize=7.8))
}
bg <- read_ppm(file.path(root,'data/SI_FigureS8_original.ppm'))
whole <- function() {
  grid.newpage();grid.raster(bg,width=1,height=1,interpolate=FALSE)
  pushViewport(viewport(x=.710,y=.515,width=.575,height=.310,gp=gpar(fontfamily='Arial')))
  panel();popViewport()
}
width <- 180/25.4; height <- width*nrow(bg)/ncol(bg)
png(file.path(root,'figures/SI_FigureS8.png'),width=width,height=height,units='in',res=400,type='cairo',bg='white');whole();dev.off()
cairo_pdf(file.path(root,'figures/SI_FigureS8.pdf'),width=width,height=height,family='Arial');whole();dev.off()
cairo_pdf(file.path(root,'figures/FigureS8D_PROseq.pdf'),width=110/25.4,height=100/25.4,family='Arial')
grid.newpage();pushViewport(viewport(gp=gpar(fontfamily='Arial')));panel(TRUE);popViewport();dev.off()
svg(file.path(root,'figures/FigureS8D_PROseq.svg'),width=110/25.4,height=100/25.4,family='Arial')
grid.newpage();pushViewport(viewport(gp=gpar(fontfamily='Arial')));panel(TRUE);popViewport();dev.off()
cat('Generated corrected panel D and full Figure S8; fixed alphabetical rows; no significance selection.\n')
