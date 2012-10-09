# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Sep. 2012
# Version 1.0
# Licence GPL v3

.level <- function(levels) {
  if (length(levels) == 1) {
    if (levels == 0) levels <- c(-6,-3,0,3,6) # default
    else levels <- c(-abs(levels),0,abs(levels))
  } else {
    if (length(levels[levels > 0]) > 0 & length(levels[levels < 0]) > 0) {
      if(length(levels[levels > 0]) == length(levels[levels < 0])) levels <- sort(c(levels[levels < 0],0,levels[levels > 0]))
      else if(length(levels[levels > 0]) > length(levels[levels < 0])) levels <- sort(c(-levels[levels > 0],0,levels[levels > 0]))
      else levels <- sort(c(levels[levels < 0],0,-levels[levels < 0]))
    }
    else if (length(levels[levels > 0]) > 0) levels <- sort(c(-levels[levels > 0],0,levels[levels > 0]))
    else levels <- sort(c(levels[levels < 0],0,-levels[levels < 0]))
  }
  levels
}

plot.specisLISA <- function(x, y, cex=2,levels, xyLegend, xlab="X Coordinates",ylab="Y Coordinates", main, ...) {
  op <- par(mar = par()$mar)
  
  par(mar=par()$mar + c(0,2,0,2))
  
  orig.cex <- cex
  if (!missing(y)) {
    dx <- ((as.vector(bbox(y)[1,2]) - as.vector(bbox(y)[1,1])) * 0.04) /2
    limx <-c(as.vector(bbox(y)[1,1]) - dx,as.vector(bbox(y)[1,2])+dx) 
    dy <- ((as.vector(bbox(y)[2,2]) - as.vector(bbox(y)[2,1])) * 0.04) /2
    limy <-c(as.vector(bbox(y)[2,1]) - dy, as.vector(bbox(y)[2,2])+dy)
    
    plot(0,0,xlim=limx,ylim=limy,xlab=xlab,ylab=ylab,main=main)
  }
  else {
    dx <- ((as.vector(bbox(x@species)[1,2]) - as.vector(bbox(x@species)[1,1])) * 0.1) /2
    limx <-c(as.vector(bbox(x@species)[1,1]) - dx,as.vector(bbox(x@species)[1,2])+dx) 
    dy <- ((as.vector(bbox(x@species)[2,2]) - as.vector(bbox(x@species)[2,1])) * 0.1) /2
    limy <-c(as.vector(bbox(x@species)[2,1]) - dy,as.vector(bbox(x@species)[2,2])+dy)
    
    plot(0,0,xlim=limx,ylim=limy,xlab="X Coordinate",ylab="Y Coordinate",main=main)
  }
  cx <- rep(NA,length(x@LISA))
  cx <- ifelse(x@LISA <= levels[1] | x@LISA >= levels[length(levels)],cex,cx)
  cx.d <- (cex - 0.5) / trunc(length(levels)/2)
  for(i in 2:(trunc(length(levels)/2)+1)) {
    cex <-  cex - cx.d
    cx <- ifelse((x@LISA > levels[i-1] & x@LISA <= levels[i]) | (x@LISA >= levels[length(levels)+1-i] & x@LISA < levels[length(levels)+2-i]),cex,cx)
  }
  pch <- ifelse(x@LISA >= 0,16,1)
  xy <- coordinates(x@species)
  points(xy[,1],xy[,2],cex=cx,pch=pch)
  if (!missing(y)) plot(y,add=T)
  
  txt <- paste("< ",levels[1],sep="")
  for (i in 2:length(levels)) txt <- c(txt,paste(levels[i-1]," : ",levels[i],sep=''))
  txt <- c(txt,paste("> ",levels[length(levels)],sep=''))
  
  cex <- orig.cex + cx.d
  cx <- c()
  while (cex != 0.5) {
    cex <- cex - cx.d
    cx <- c(cx,cex)
  }
  while (cex != (orig.cex+cx.d)) {
    cx <- c(cx,cex)
    cex <- cex + cx.d
  }
  pch <- c(rep(1,length(txt)/2),rep(16,length(txt)/2))
  legend(xyLegend[1],xyLegend[2],legend=txt,pt.cex=cx,pch=pch,title='LISA')
  par(op)
}



if (!isGeneric("plot")) {
  setGeneric("plot", function(x,y,...)
    standardGeneric("plot"))
}	


setMethod("plot", signature(x='speciesLISA',y="SpatialPolygons"), 
          function(x,y,cex=2,levels=c(0,3,6), xyLegend, xlab="X Coordinates",ylab="Y Coordinates", main, ...) {
            if (missing(xyLegend)) xyLegend <- c(bbox(y)[1,2] - (bbox(y)[1,2]-bbox(y)[1,1]) * 0.16,bbox(y)[2,1] + (bbox(y)[2,2]-bbox(y)[2,1]) * 0.25)
            else if(length(xyLegend) != 2 | class(xyLegend) != 'numeric') xyLegend <- c(bbox(y)[1,2] - (bbox(y)[1,2]-bbox(y)[1,1]) * 0.16,bbox(y)[2,1] + (bbox(y)[2,2]-bbox(y)[2,1]) * 0.25)
            
            if (missing(main)) main <- "Impact of positional uncertainty based on LISA"
            levels <- .level(levels)            
            
            plot.specisLISA(x=x,y=y,levels=levels,xyLegend=xyLegend,xlab=xlab,ylab=ylab, main=main, ...)
          }
)

setMethod("plot", signature(x='speciesLISA',y="SpatialPolygonsDataFrame"), 
          function(x,y,cex=2,levels=c(0,3,6), xyLegend, xlab="X Coordinates",ylab="Y Coordinates", main, ...) {
            if (missing(xyLegend)) xyLegend <- c(bbox(y)[1,2] - (bbox(y)[1,2]-bbox(y)[1,1]) * 0.16,bbox(y)[2,1] + (bbox(y)[2,2]-bbox(y)[2,1]) * 0.25)
            else if(length(xyLegend) != 2 | class(xyLegend) != 'numeric') xyLegend <- c(bbox(y)[1,2] - (bbox(y)[1,2]-bbox(y)[1,1]) * 0.16,bbox(y)[2,1] + (bbox(y)[2,2]-bbox(y)[2,1]) * 0.25)
            
            if (missing(main)) main <- "Impact of positional uncertainty based on LISA"
            levels <- .level(levels)            
            
            plot.specisLISA(x=x,y=y,levels=levels,xyLegend=xyLegend,xlab=xlab,ylab=ylab, main=main, ...)
          }
)

setMethod("plot", signature(x='speciesLISA',y="missing"), 
          function(x,y,cex=2,levels=c(0,3,6), xyLegend, xlab="X Coordinates",ylab="Y Coordinates", main, ...) {
            if (missing(xyLegend)) xyLegend <- c(bbox(x@species)[1,2] - (bbox(x@species)[1,2]-bbox(x@species)[1,1]) * 0.16,bbox(x@species)[2,1] + (bbox(x@species)[2,2]-bbox(x@species)[2,1]) * 0.25)
            else if(length(xyLegend) != 2 | class(xyLegend) != 'numeric') xyLegend <- c(bbox(x@species)[1,2] - (bbox(x@species)[1,2]-bbox(x@species)[1,1]) * 0.16,bbox(x@species)[2,1] + (bbox(x@species)[2,2]-bbox(x@species)[2,1]) * 0.25)
            
            if (missing(main)) main <- "Impact of positional uncertainty based on LISA"
            levels <- .level(levels)            
            
            plot.specisLISA(x=x,levels=levels,xyLegend=xyLegend,xlab=xlab,ylab=ylab, main=main, ...)
          }
)