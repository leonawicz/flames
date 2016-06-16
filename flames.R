setwd("C:/github/flames")
library(raster)
library(jpeg)

jpg1 <- readJPEG("source_image.jpg")
m <- 0.21*jpg1[,,1] + 0.71*jpg1[,,2] + 0.07*jpg1[,,3] # RGB to grayscale (luminosity method)
r <- raster(m, xmn=0, xmx=ncol(m), ymn=0, ymx=nrow(m))

base.width <- dim(r)[2]
r.agg <- aggregate(r, rep(base.width/40, 2))
r.agg2 <- aggregate(r, rep(base.width/24, 2))
m.agg <- as.matrix(r.agg)
m.agg2 <- as.matrix(r.agg2)

# Support functions for main plotting function
.gridlines <- function(v,h,...){
  abline(v=v, col="#FFFFFF20")
  abline(h=h, col="#FFFFFF20")
}

.plot_pic <- function(x, m, return.seq=T, type="rgb", ...){
  colvec <- switch(type, grey=grey(x), rgb=rgb(x[,,1], x[,,2], x[,,3]))
  colors <- unique(colvec)
  colmat <- array(match(colvec, colors), dim=dim(x)[1:2])

  seq.x <- seq(0-0.5/ncol(m),ncol(m)+0.5/ncol(m),length=ncol(colmat))
  seq.y <- seq(0-0.5/nrow(colmat),1+0.5/nrow(colmat),length=nrow(colmat))
  par(mar=c(0,0,0,0))
  image(x=seq.x, y=seq.y,
        z=t(colmat[nrow(colmat):1,]), col=colors, xlab="", ylab="", axes=FALSE, xlim=range(seq.x), ylim=range(seq.y), ...)
  if(return.seq) return(list(x=seq.x, y=seq.y)) else return(NULL)
}

.plot_heatmap <- function(m, xy, showgrid=T, return.seq=T, alpha){
  seq.x <- seq(xy$x[1], tail(xy$x,1), length=ncol(m))
  seq.y <- seq(xy$y[1], tail(xy$y,1), length=nrow(m))
  image(x=seq.x, y=seq.y,
        z=t(m[nrow(m):1,]), col=paste0(colorRampPalette(c("black","red","yellow","white"))(100),alpha), xlim=range(seq.x), ylim=range(seq.y), add=T)
  if(showgrid) .gridlines(v=seq.x, h=seq.y, col="#FFFFFF20")
  if(return.seq) return(list(x=seq.x, y=seq.y)) else return(NULL)
}

.sample_under_curve <- function(m, n, f){
  smp <- c()
  for(i in 1:(nrow(m))) smp <- c(smp, sample(x=1:ncol(m), size=n, prob=f(m[i,]^(i/nrow(m))), rep=T))
  table(smp)
}

.plot_streaks_over_segments <- function(bp, xy, tb, exp.alpha=1.5){
  bp2_vals <- 1-tb/max(tb)
  bp2.y0 <- -(bp2_vals) + diff(c(0,max(xy$y)))
  pts.m <- c()
  for(i in 1:length(bp2.y0)){
    pts.x <- runif(50, bp[i,1]-0.5, bp[i,1]+0.5)
    pts.y <- abs(rnorm(50, 0, 0.1)) + bp2.y0[i]
    pts.m <- rbind(pts.m, cbind(pts.x,pts.y))
  }
  pts.alpha <- (101 - cut(pts.m[,2], 100, labels=F))^exp.alpha
  pts.alpha <- paste0(0, round(99*( (pts.alpha-min(pts.alpha))/max(pts.alpha-min(pts.alpha)) )))
  pts.alpha <- substr(pts.alpha, nchar(pts.alpha)-1, nchar(pts.alpha))
  points(pts.m[,1], pts.m[,2], col=paste0("#00FFFF",pts.alpha), pch="|", cex=runif(50, 0.1,0.6))
  segments(x0=bp[,1]-0.5, x1=bp[,1]+0.5, y0 = bp2.y0, col="#00FFFF")
}

# Main plotting function
plot_pic_mat_hist <- function(x, mat, mat.agg, alpha=60, n=1000, f.list=list(function(x) x^1.2, function(x) x^0.5), col.hist=c("#ADFF2F30","#00FFFF05"), seed=745, mirror=T, ...){
  set.seed(seed)
  xy <- .plot_pic(x=x, m=mat)
  xy <- .plot_heatmap(m=mat.agg, xy=xy, showgrid=T, alpha=alpha)
  tb <- .sample_under_curve(m=mat, n=n, f=f.list[[1]])
  bp <- barplot(tb/max(tb), space=0, offset=diff(c(0,xy$y[1])), col=col.hist[1], add=T)
  interp <- approx(bp[,1], tb/max(tb), seq(bp[1,1],tail(bp[,1],1),length=1000))
  spl <- smooth.spline(interp$x, interp$y, df=25)
  lines(spl$x, spl$y, lwd=2, col="#FF4500")
  if(mirror){
    tb <- .sample_under_curve(m=mat, n=n, f=f.list[[2]])
    .plot_streaks_over_segments(bp=bp, xy=xy, tb=tb)
    interp <- approx(bp[,1], tb/max(tb), seq(bp[1,1],tail(bp[,1],1),length=1000))
    spl2 <- smooth.spline(interp$x, interp$y, df=25)
    lines(spl2$x, spl2$y, lwd=2, col="#00FFFF")
    lines(spl2$x, (spl$y + spl2$y)/2, lwd=2, col="white")
  }
}

jpeg("plots/flames_plotzzz.jpg", width=ncol(jpg1)/2, height=nrow(jpg1)/2, res=300)
plot_pic_mat_hist(x=jpg1, mat=m.agg, mat.agg=m.agg2)
dev.off()
