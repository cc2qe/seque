paintCytobands<-function(chrom, pos=c(0,0), units=c("cM","bases","ISCN"), width=0.4, length.out, bands="major", orientation=c("h","v"), legend = TRUE, cex.leg=0.7, bleach = 0,...) {
  # Based on paintCytobands from quantsmooth package
  # added:
  #  -bleach
  #  -length.out
  #  -using all of cM,bases,ISCN
  #  -using hatches for stalk, acen
  #  -legend + cex.leg
  #  -orientation
  #  extracted semicircle for general use
  bleacher<-function(x) { (x * (1-bleach)) + bleach}
  chrom.bands<-NULL;rm(chrom.bands) # trick to satisfy R check
  
  # this is important!!!!
  data(chrom.bands,package="quantsmooth",envir=environment())
  
  
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
    lc<-nchar(chromdata$band)
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
    idx<-c(2:(centromere-1), (centromere+2):(n-1))
    if (orientation=="h") {
      rect(pos[1]+bandpos[idx,1],pos[2],pos[1]+bandpos[idx,2],pos[2]-width, col=bandcol[idx], density=banddens[idx], border=bandbord[idx])
      qs.semicircle(pos[1]+bandpos[1,2], pos[2]-width, width,
                 bandpos[1,2]-bandpos[1,1], 2, col=bandcol[1], density=banddens[1], border=bandbord[1],...)
      qs.semicircle(pos[1]+bandpos[n,1], pos[2]-width, width,
                 bandpos[n,2]-bandpos[n,1], 4, col=bandcol[n], density=banddens[n], border=bandbord[n],...)
      qs.semicircle(pos[1]+bandpos[centromere,1], pos[2]-width, width,
                 bandpos[centromere,2]-bandpos[centromere,1],
                 4, col=bandcol[centromere], density=banddens[centromere], border=bandbord[centromere],...)
      qs.semicircle(pos[1]+bandpos[centromere+1,2], pos[2]-width, width,
                 bandpos[centromere+1,2]-bandpos[centromere+1,1],
                 2, col=bandcol[centromere+1], density=banddens[centromere+1], border=bandbord[centromere+1],...)

      centromere.size=0.6*0.5*width/yinch(1)
      symbols(pos[1]+bandpos[centromere,2], pos[2]-0.5*width,circles=1,inches=centromere.size, add=TRUE,fg=gray(bleacher(0)),bg="white",...)
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