# based on paintCytobands from the quantsmooth package by David Duffy
# Colby Chiang (cc2qe@virginia.edu)

drawIdeogram<-function(chrom, pos=c(0,0), units="bases", genome=c('hg18', 'hg19', 'mm9', 'rn4'), width=0.4, length.out, bands="major", orientation=c("h","v"), legend=TRUE, cex.leg=0.7, bleach=0, ...) {
  bleacher<-function(x) { (x * (1-bleach)) + bleach}
  chrom.bands<-NULL;rm(chrom.bands) # trick to satisfy R check
  
  # set the cytoband annotation file
  #data(chrom.bands,package="quantsmooth",envir=environment())
  
  cytoFile <- paste("~/code/seque/drawIdeogram/annot/", genome, ".cytoR.txt.gz", sep='')
  
  chrom.bands <- read.table(cytoFile, sep="\t", col.names=c("chr", "arm", "band", "ISCN.top", "ISCN.bot", "bases.top", "bases.bot", "stain", "cM.top", "cM.bot", "n.markers", "p.markers"))
  
  chrom<-switch(as.character(chrom),
         "98"="X",
         "99"="Y",
         as.character(chrom))
  units<-match.arg(units)
  orientation<-match.arg(orientation)
  # original function only required ypos
  if (length(pos)==1) pos<-c(0,pos)
  chromdata<-subset(chrom.bands, chrom.bands$chr==chrom)
  if (nrow(chromdata)>0){
    lc<-nchar(as.character(chromdata$band))
    sel<-!(substr(chromdata$band,lc,lc) %in% letters)
    if (bands!="major") sel<-!sel
    chromdata<-chromdata[sel,]
    rm(lc,sel)
    bandpos<-switch(units,
           cM =chromdata[,c("cM.top","cM.bot")],
           bases = chromdata[,c("bases.top","bases.bot")],
           ISCN =  chromdata[,c("ISCN.top","ISCN.bot")])

    type.b<-match(chromdata$stain,c("acen","gneg", "gpos", "gvar", "stalk"))
    bandcol<-gray(bleacher(c(0.5,1,0.2,0.6,0.75)))[type.b]
    banddens<-c(30,-1,-1,-1,10)[type.b]
    bandbord<-gray(bleacher(c(0,0,0,0,1)))[type.b]
    if (!missing(length.out)) {
      bandpos<-(bandpos/max(bandpos))*length.out
    }
    n<-nrow(chromdata)
    centromere<-which(chromdata$arm[-n]!=chromdata$arm[-1])
    if (length(centromere) == 0) { centromere <- 0 }
    idx<-c(2:(max(centromere-1,2)), (centromere+2):(n-1))
    if (orientation=="h") {
      rect(pos[1]+bandpos[idx,1],pos[2],pos[1]+bandpos[idx,2],pos[2]-width, col=bandcol[idx], density=banddens[idx], border=bandbord[idx])
      
      # centromeric semicircles only if it's not acrocentric
      if (2 != 1) {
      	qs.semicircle(pos[1]+bandpos[1,2], pos[2]-width, width,
                 bandpos[1,2]-bandpos[1,1], 2, col=bandcol[1], density=banddens[1], border=bandbord[1],...)

      	qs.semicircle(pos[1]+bandpos[centromere,1], pos[2]-width, width,
                 bandpos[centromere,2]-bandpos[centromere,1],
                 4, col=bandcol[centromere], density=banddens[centromere], border=bandbord[centromere],...)
      # q-ter semicircle
      qs.semicircle(pos[1]+bandpos[n,1], pos[2]-width, width,
                 bandpos[n,2]-bandpos[n,1], 4, col=bandcol[n], density=banddens[n], border=bandbord[n],...)
      # p-ter semicircle
      qs.semicircle(pos[1]+bandpos[centromere+1,2], pos[2]-width, width,
                 bandpos[centromere+1,2]-bandpos[centromere+1,1],
                 2, col=bandcol[centromere+1], density=banddens[centromere+1], border=bandbord[centromere+1],...)
      }
      else {
      # p-ter semicircle
      qs.semicircle(pos[1]+bandpos[1,2], pos[2]-width, width,
                 bandpos[1,2]-bandpos[1,1],
                 2, col=bandcol[1], density=banddens[1], border=bandbord[1],...)
      # q-ter semicircle
      qs.semicircle(pos[1]+bandpos[n,1], pos[2]-width, width,
                 bandpos[n,2]-bandpos[n,1], 4, col=bandcol[n], density=banddens[n], border=bandbord[n],...)
      	
      }

      centromere.size=0.6*0.5*width/yinch(1)
      symbols(pos[1]+bandpos[centromere+1,1], pos[2]-0.5*width,circles=1,inches=centromere.size, add=TRUE,fg=gray(bleacher(0)),bg="white",...)
      if (legend) text(pos[1]+(bandpos[,1]+bandpos[,2])/2,pos[2]+0.5*width,paste(chromdata[,"arm"],chromdata[,"band"],sep=""),adj=c(0,0.5),srt=90,cex=cex.leg,...)
    } else {
      rect(pos[1],pos[2]-bandpos[idx,1],pos[1]-width,pos[2]-bandpos[idx,2], col=bandcol[idx], density=banddens[idx], border=bandbord[idx],...)
      qs.semicircle(pos[1]-width, pos[2]-bandpos[1,2], width,
                 bandpos[1,2]-bandpos[1,1], 3, col=bandcol[1], density=banddens[1], border=bandbord[1],...)
      qs.semicircle(pos[1]-width, pos[2]-bandpos[n,1], width,
                 bandpos[n,2]-bandpos[n,1], 1, col=bandcol[n], density=banddens[n], border=bandbord[n],...)
      qs.semicircle(pos[1]-width, pos[2]-bandpos[centromere,1], width,
                 bandpos[centromere,2]-bandpos[centromere,1],
                 1, col=bandcol[centromere], density=banddens[centromere], border=bandbord[centromere],...)
      qs.semicircle(pos[1]-width, pos[2]-bandpos[centromere+1,2], width,
                 bandpos[centromere+1,2]-bandpos[centromere+1,1],
                 3, col=bandcol[centromere+1], density=banddens[centromere+1], border=bandbord[centromere+1],...)
      centromere.size=0.6*0.5*width/xinch(1)
      symbols(pos[1]-0.5*width, pos[2]-bandpos[centromere,2],circles=1,inches=centromere.size, add=TRUE,fg=gray(bleacher(0)),bg="white",...)
      if (legend) text(pos[1]+0.5*width,pos[2]-(bandpos[,1]+bandpos[,2])/2,paste(chromdata[,"arm"],chromdata[,"band"],sep=""),adj=c(0,0.5),srt=0,cex=cex.leg,...)
    }
  } else {
    warning(paste("Chromosome",chrom,"is not plotted because cytoband data is not available"))
  }
}

qs.semicircle <- function(base.x, base.y, base.length, height=base.length, side=1, orientation=NULL,plottype="poly",...) {
  # based on lodplot package
  # David Duffy <David.Duffy@qimr.edu.au>
  # URL: http://www.qimr.edu.au/davidD
  # - col is now propagated through ..., other plotting parameters can now also be given
  # - different types poly/line
  radius<-base.length/2
  x<-radius*seq(-1,1,length=40)
  y<-height/radius*sqrt(radius^2-x^2)
  if (is.null(orientation)) {
    co<-as.integer(cos(pi*(3-side)/2))
    so<-as.integer(sin(pi*(3-side)/2))
  }else{
    co<-cos(orientation)
    so<-sin(orientation)
  }
  tx<-co*x - so*y
  ty<-so*x + co*y
  if (is.null(orientation)) {
    if (side==1 || side==3) {
      base.x<-base.x+radius
    }else if (side==2 || side==4) {
      base.y<-base.y+radius
    }
  }
  x<-base.x+tx
  y<-base.y+ty
  switch(plottype,
    poly=polygon(x,y,...),
    line=lines(x,y,...)
  )
}

getChromLength <- function(chrom, genome=c('hg18', 'hg19', 'mm9', 'rn4')) {
  # set the cytoband annotation file
  #data(chrom.bands,package="quantsmooth",envir=environment())
  if (genome == 'hg19') {
  	cytoFile <- "~/code/seque/drawIdeogram/annot/hg19.cytoR.txt.gz"
  }
  else if (genome == 'hg18') {
  	cytoFile <- "~/code/seque/drawIdeogram/annot/hg18.cytoR.txt.gz"
  }
  else if (genome == 'mm9') {
  	cytoFile <- "~/code/seque/drawIdeogram/annot/mm9.cytoR.txt.gz"
  }
  else if (genome == 'rn4') {
  	cytoFile <- "~/code/seque/drawIdeogram/annot/rn4.cytoR.txt.gz"
  }

  chrom.bands <- read.table(cytoFile, sep="\t", col.names=c("chr", "arm", "band", "ISCN.top", "ISCN.bot", "bases.top", "bases.bot", "stain", "cM.top", "cM.bot", "n.markers", "p.markers"))
  
  chromdata<-subset(chrom.bands, chrom.bands$chr==chrom)
  
  return(chromdata$bases.bot[nrow(chromdata)])
}