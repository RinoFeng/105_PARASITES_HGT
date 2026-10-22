args<-commandArgs(TRUE);o<-normalizePath(args[1]);d<-read.delim(file.path(o,'data/Figure3A_ring_counts.tsv'));domains<-c('Bacteria','Eukaryota','Archaea','Viruses');cols<-c(Bacteria='#E74C3C',Eukaryota='#2ECC71',Archaea='#F39C12',Viruses='#8266B0');d<-d[order(match(d$domain,domains),-d$n),]
arc<-function(a,b,r0,r1,col){t<-seq(a,b,length.out=max(5,ceiling(500*abs(a-b)/(2*pi))));polygon(c(r1*cos(t),rev(r0*cos(t))),.19+c(r1*sin(t),rev(r0*sin(t))),col=col,border='white',lwd=.7)}
draw<-function(){
 par(mar=c(0,0,0,0),family='Arial',xaxs='i',yaxs='i');plot.new();plot.window(xlim=c(-1.50,1.50),ylim=c(-1.27,1.23),asp=1);text(-1.43,1.14,'A',font=2,cex=1.35,adj=0)
 start<-pi/2
 for(dom in domains){n<-sum(d$n[d$domain==dom]);end<-start-2*pi*n/8989;arc(start,end,.30,.67,cols[dom]);mid<-(start+end)/2
  if(n/8989>.02)text(.49*cos(mid),.19+.49*sin(mid),sprintf('%s\n(%.1f%%)',dom,100*n/8989),cex=.92,font=2)
  start<-end
 }
 start<-pi/2
 for(i in seq_len(nrow(d))){end<-start-2*pi*d$n[i]/8989;c<-col2rgb(cols[d$domain[i]])/255;c<-.7*c+.3;arc(start,end,.677,1.035,rgb(c[1],c[2],c[3]));mid<-(start+end)/2
  # Text is drawn after all sectors, so later polygons cannot cover labels.
  start<-end
 }
 start<-pi/2
 for(i in seq_len(nrow(d))){end<-start-2*pi*d$n[i]/8989;mid<-(start+end)/2;if(d$n[i]/8989>=.05)text(.865*cos(mid),.19+.865*sin(mid),paste0(d$lineage_group[i],'\n(',format(d$n[i],big.mark=','),')'),cex=.79);start<-end}
 text(0,.19,'Total\n(8,989)',font=2,cex=1.15)
 yy<-c(-.95,-1.08);xx<-c(-1.12,.28)
 for(i in 1:4){j<-(i-1)%%2+1;k<-(i-1)%/%2+1;n<-sum(d$n[d$domain==domains[i]]);rect(xx[j]-.06,yy[k]-.025,xx[j]-.015,yy[k]+.025,col=cols[i],border=NA);text(xx[j],yy[k],sprintf('%s: %s (%.1f%%)',domains[i],format(n,big.mark=','),100*n/8989),adj=0,cex=.80)}
}
f<-file.path(o,'figures/Figure3A_current_sources');png(paste0(f,'.png'),width=1320,height=1100,res=200,type='cairo',pointsize=14);draw();dev.off();cairo_pdf(paste0(f,'.pdf'),width=6.6,height=5.5,pointsize=14,family='Arial');draw();dev.off();svglite::svglite(paste0(f,'.svg'),width=6.6,height=5.5,pointsize=14,fix_text_size=FALSE);draw();dev.off()
