static char *smu_plot =
"#!/usr/bin/env Rscript \n\
 \n\
require(\"argparse\") \n\
require(\"ggplot2\") \n\
require(\"scales\") \n\
require(\"viridis\") \n\
require(\"cowplot\") \n\
 \n\
parser <- ArgumentParser(description = \"Make spectra-cn plots. Line, filled, and stacked spectra-cn plots will be generated.\") \n\
parser$add_argument(\"-f\", \"--file\", type=\"character\", help=\".spectra-cn.hist file (required)\", default=NULL) \n\
parser$add_argument(\"-a\", \"--assign\", type=\"character\", help=\".spectra-cn.hist file (required)\", default=NULL) \n\
parser$add_argument(\"-o\", \"--output\", type=\"character\", help=\"output prefix (required)\") \n\
parser$add_argument(\"-x\", \"--xdim\", type=\"double\", default=5, help=\"width of plot [default %(default)s]\") \n\
parser$add_argument(\"-y\", \"--ydim\", type=\"double\", default=5, help=\"height of plot [default %(default)s]\") \n\
parser$add_argument(\"-t\", \"--type\", type=\"character\", default=\"all\", help=\"available types: line, fill, stack, or all. [default %(default)s]\") \n\
parser$add_argument(\"-p\", \"--pdf\", dest='pdf', default=FALSE, action='store_true', help=\"get output in .pdf. [default .png]\") \n\
parser$add_argument(\"-s\", \"--source\", type=\"character\", help=\"source .ktab file (required)\", default=NULL) \n\
args <- parser$parse_args() \n\
 \n\
fancy_scientific <- function(d) { \n\
  if (d[2] > 1000000) { \n\
    for (i in 1:length(d)) { \n\
      if (is.na(d[i])) { \n\
        next \n\
      } \n\
      d[i] <- paste( as.character(as.integer(d[i])/1000000), \"M\", sep=\"\") \n\
    } \n\
  } else if (d[2] > 1000) { \n\
    for (i in 1:length(d)) { \n\
      if (is.na(d[i])) { \n\
        next \n\
      } \n\
      d[i] <- paste( as.character(as.integer(d[i])/1000), \"K\", sep=\"\") \n\
    } \n\
  } else { \n\
    for (i in 1:length(d)) { \n\
      if (is.na(d[i])) { \n\
        next \n\
      } \n\
      d[i] <- as.character(as.integer(d[i])) \n\
    } \n\
  } \n\
  d \n\
} \n\
 \n\
format_theme <- function() { \n\
    theme(legend.text = element_text(size=8), \n\
          legend.title = element_text(hjust=.5, size=8), \n\
          legend.margin = margin(t=3, b=3, l=3, r=3, unit='pt'), \n\
          plot.title = element_text(face=\"bold\",size=14,hjust=.5,vjust=.5, \n\
                                    margin=margin(t=12,b=18,unit=\"pt\")), \n\
          plot.margin =unit(c(0,0,6,6),\"pt\"), \n\
          axis.title=element_text(size=12), \n\
          axis.text=element_text(size=11), \n\
          plot.background = element_blank(), \n\
          panel.background = element_blank(), \n\
          panel.grid.major = element_blank(), \n\
          panel.grid.minor = element_blank(), \n\
	  panel.border = element_rect(colour=\"black\", fill=NA, size=2)) \n\
} \n\
 \n\
layer_on_contour <- function(with) { \n\
  if (with) \n\
    geom_contour(aes(z=Count), colour=\"white\", bins=6, show.legend=FALSE) \n\
} \n\
 \n\
plot_heat <- function(dat, h, with, ybrk, ylab) { \n\
 \n\
  y_max <- max(dat[,1]) + .5; \n\
  x_max <- max(dat[,2]) + .5; \n\
 \n\
  M <- ggplot(data=dat, aes(x=KF2,y=KF1,z=Count)) + \n\
    geom_tile(aes(fill=Count)) + \n\
    layer_on_contour(with) + \n\
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max), expand=FALSE) + \n\
    scale_fill_viridis(name=\"kmer pairs\", labels=fancy_scientific) + \n\
    scale_x_continuous(breaks=c(12,17,25,33,49), labels=c(\"1/8\",\"1/6\",\"1/4\",\"1/3\",\"1/2\")) + \n\
    scale_y_continuous(breaks=ybrk, labels=ylab) + \n\
    guides(fill = guide_colourbar(title.position = \"top\", frame.colour=\"black\", \n\
                                  ticks.colour=\"black\", barheight = unit(.14*h,\"in\"))) + \n\
    xlab(paste(\"Normalized minor kmer coverage: B/(A+B)\", sep=\"\")) + \n\
    ylab(paste(\"Total coverage of the kmer pair: A+B\", sep=\"\")) + \n\
    format_theme(); \n\
 \n\
  M \n\
} \n\
 \n\
plot_contour <- function(dat, h, ybrk, ylab) { \n\
 \n\
  y_max <- max(dat[,1]) + .5; \n\
  x_max <- max(dat[,2]) + .5; \n\
 \n\
  M <- ggplot(data=dat, aes(x=KF2,y=KF1,z=Count)) + \n\
    geom_contour(aes(z=Count, colour=..level..), bins=12, show.legend=TRUE) + \n\
    coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max), expand=FALSE) + \n\
    scale_color_viridis(option=\"turbo\", name=\"kmer pairs\", \n\
                        labels=fancy_scientific, breaks = 6) + \n\
    scale_x_continuous(breaks=c(12,17,25,33,49), labels=c(\"1/8\",\"1/6\",\"1/4\",\"1/3\",\"1/2\")) + \n\
    scale_y_continuous(breaks=ybrk, labels=ylab) + \n\
    guides(color = guide_legend(title.position = \"top\", frame.colour=\"black\", \n\
                                ticks.colour=\"black\", reverse=TRUE, \n\
                                keyheight = unit(1.5*h,\"pt\"))) + \n\
    xlab(paste(\"Normalized minor kmer coverage: B/(A+B)\", sep=\"\")) + \n\
    ylab(paste(\"Total coverage of the kmer pair: A+B\", sep=\"\")) + \n\
    format_theme() \n\
 \n\
  M \n\
} \n\
 \n\
save_plot <- function(name, type, outformat, h, w) { \n\
  ggsave(file = paste(name, type, outformat, sep = \".\"), height = h, width = w) \n\
} \n\
 \n\
my_plot <- function(hist, asn, name, h=5, w=5, x_max, type, pdf=FALSE, s) { \n\
 \n\
  dat=read.table(hist, header=TRUE) \n\
  sms=read.table(asn, header=TRUE) \n\
 \n\
  outformat=\"png\" \n\
  if (pdf) { \n\
    outformat=\"pdf\" \n\
  } \n\
 \n\
  cover  <- sms[1,5]; \n\
  ploidy <- sms[1,6]+1; \n\
  Ploid  <- c( \"Diploid\", \"Triploid\", \"Tetraploid\", \"Hexaploid\", \"Octaploid\" ); \n\
  Nhaps  <- c( 2, 3, 4, 6, 8 ); \n\
 \n\
  icov <- sms[1,3]; \n\
  if (icov/2 <= 10 || ploidy == 1) \n\
    { unit <- icov/2; \n\
      fact <- Nhaps[ploidy]/2; \n\
    } \n\
  else if (ploidy == 2 || (ploidy == 4 && icov/3 <= 10)) \n\
    { unit <- icov/3; \n\
      fact <- Nhaps[ploidy]/3; \n\
    } \n\
  else if (ploidy == 3 || (ploidy == 5 && icov/4 <= 10)) \n\
    { unit <- icov/4; \n\
      fact <- Nhaps[ploidy]/4; \n\
    } \n\
  else if (ploidy == 4) \n\
    { unit <- icov/6; \n\
      fact <- Nhaps[ploidy]/6; \n\
    } \n\
  else \n\
    { unit <- icov/8; \n\
      fact <- Nhaps[ploidy]/8; \n\
    } \n\
 \n\
  ybrk <- c() \n\
  ylab <- c() \n\
  for (i in 1:(50/unit)) \n\
    { ybrk <- append(ybrk, unit*i); \n\
      ylab <- append(ylab, paste(i*fact,\"n\",sep=\"\")); \n\
    } \n\
 \n\
  if (type == \"combo\") \n\
    M <- plot_heat(dat, h, TRUE, ybrk, ylab) \n\
  else if (type == \"heat\") \n\
    M <- plot_heat(dat, h, FALSE, ybrk, ylab) \n\
  else # type == \"contour\" \n\
    { M <- plot_contour(dat, h, ybrk, ylab) \n\
 \n\
      bord <- rep(50,50); \n\
      for (i in 1:nrow(dat)) \n\
        { x <- trunc(dat[i,1])+1; \n\
          if (dat[i,2] < bord[x]) \n\
            bord[x] <- dat[i,2]; \n\
        } \n\
      for (x in 1:50) \n\
        bord[x] <- trunc(bord[x]); \n\
 \n\
      bfam <- as.data.frame(bord); \n\
 \n\
      B <- ggplot(bfam, aes(x=(1:50)-.5,y=bord)) +  \n\
             geom_col(fill=\"light grey\", color=\"light grey\", width=1.0) + \n\
             coord_flip(xlim=c(0,50), expand=FALSE) + \n\
             theme(axis.title=element_blank(), axis.text=element_blank(), \n\
                                               axis.ticks=element_blank()) + \n\
             theme(plot.margin=unit(c(0,0,0,0),\"pt\")) + \n\
             theme(panel.grid=element_blank(), panel.background=element_blank()); \n\
    } \n\
 \n\
  hc = rgb(0.8352, 0.2431, 0.3098); \n\
 \n\
  E <- get_legend(M); \n\
 \n\
  M <- M + theme(legend.position='none'); \n\
 \n\
  cvr <- aggregate(x=dat$Count,by=list(dat$KF1),FUN=sum); \n\
  C <- ggplot(cvr, aes(x=Group.1,y=x)) + geom_col(fill=hc, color=\"black\", width=1.0) + \n\
         coord_flip(xlim=c(0,50), expand=FALSE) + \n\
         theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) + \n\
         theme(plot.margin=unit(c(0,3,0,3),\"pt\")) + \n\
         theme(panel.grid=element_blank(), panel.background=element_blank()); \n\
 \n\
  hpr <- aggregate(x=dat$Count,by=list(dat$KF2),FUN=max) \n\
  H <- ggplot(hpr, aes(x=Group.1,y=x)) + geom_col(fill=hc, color=\"black\", width=1.0) + \n\
         coord_cartesian(xlim=c(0,50), expand=FALSE) + \n\
         theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) + \n\
         theme(plot.margin=unit(c(3,0,3,0),\"pt\")) + \n\
         theme(panel.grid=element_blank(), panel.background=element_blank()); \n\
 \n\
 \n\
  if (trunc(50/unit)*fact < 10) \n\
    { lmargin <- (5*(.082+.0166))/w; } \n\
  else \n\
    { lmargin <- (5*(.10+.0166))/w; } \n\
  bmargin <- (5*(.077+.0166))/h; \n\
 \n\
  deco <- c() \n\
  xpos <- c() \n\
  ypos <- c() \n\
  for (i in 1:nrow(sms)) \n\
    { deco <- append(deco, paste(sms[i,1],sprintf(\"%2.0f%%\",sms[i,4]))) \n\
      xpos <- append(xpos, .83); \n\
      ypos <- append(ypos, .76-(.03*(i-1))*(5/h)); \n\
    } \n\
 \n\
  smud <- c() \n\
  slab <- c() \n\
  xsmu <- c() \n\
  ysmu <- c() \n\
  for (i in 1:nrow(sms)) \n\
    { smud <- append(smud, \"+\"); \n\
      if (sms[i,2] >= 47 && ploidy >= 2) \n\
        { slab <- append(slab, paste(sms[i,1], \" \", sep=\"\")); } \n\
      else \n\
        { slab <- append(slab, sms[i,1]); } \n\
      xsmu <- append(xsmu, lmargin + (0.8-lmargin)*(sms[i,2]+.5)/50.); \n\
      ysmu <- append(ysmu, bmargin + (0.8-bmargin)*(sms[i,3]+.5)/50.); \n\
    } \n\
 \n\
  if (type == \"contour\") \n\
     { ggdraw() + \n\
       draw_plot(B,x=lmargin,y=bmargin,width=0.80-lmargin,height=0.80-bmargin) + \n\
       draw_plot(M,x=0.0,y=0.0,width=0.8,height=0.8) + \n\
       draw_plot(H,x=lmargin,y=0.8,width=0.80-lmargin,height=0.20) + \n\
       draw_plot(C,x=0.80,y=bmargin,width=0.20,height=0.80-bmargin) + \n\
       draw_plot(E,x=0.80,y=0.80,width=0.20,height=0.20) + \n\
       draw_label(paste(s), x=0.15/w, y=.95, hjust=0, fontface=\"bold.italic\", size=40/.pt) + \n\
       draw_text(paste(Ploid[ploidy]), x=0.25/w, y=.90, hjust=0, size=32/.pt) + \n\
       draw_text(paste(\"1n =\",sprintf(\"%.1f\",cover/Nhaps[ploidy])), x=lmargin+(.1/w), \n\
                                             y=bmargin+(.1/h), hjust=0, vjust=0, size=32/.pt) + \n\
       # draw_text(\"+\", x=lmargin, y=bmargin, hjust=.5, vjust=.5, size=32/.pt) + \n\
       draw_text( deco, xpos, ypos, hjust=0, vjust=0, size=28/.pt, family=\"mono\") + \n\
       draw_text( smud, xsmu, ysmu, hjust=.5, vjust=.5, size=28/.pt, family=\"mono\") + \n\
       draw_text( slab, xsmu, ysmu, hjust=.5, vjust=-1., size=(20/.pt)*(w/5), family=\"mono\"); \n\
    } \n\
  else \n\
     { ggdraw() + \n\
       draw_plot(M,x=0.0,y=0.0,width=0.8,height=0.8) + \n\
       draw_plot(H,x=lmargin,y=0.8,width=0.80-lmargin,height=0.20) + \n\
       draw_plot(C,x=0.80,y=bmargin,width=0.20,height=0.80-bmargin) + \n\
       draw_plot(E,x=0.80,y=0.80,width=0.20,height=0.20) + \n\
       draw_label(paste(s), x=0.15/w, y=.95, hjust=0, fontface=\"bold.italic\", size=40/.pt) + \n\
       draw_text(paste(Ploid[ploidy]), x=0.25/w, y=.90, hjust=0, size=32/.pt) + \n\
       draw_text(paste(\"1n =\",sprintf(\"%.1f\",cover/Nhaps[ploidy])), x=lmargin+(.1/w), \n\
                                             y=bmargin+(.1/h), hjust=0, vjust=0, size=32/.pt) + \n\
       # draw_text(\"+\", x=lmargin, y=bmargin, hjust=.5, vjust=.5, size=32/.pt) + \n\
       draw_text( deco, xpos, ypos, hjust=0, vjust=0, size=28/.pt, family=\"mono\") + \n\
       draw_text( smud, xsmu, ysmu, hjust=.5, vjust=.5, size=28/.pt, family=\"mono\") + \n\
       draw_text( slab, xsmu, ysmu, hjust=.5, vjust=-1., size=(20/.pt)*(w/5), family=\"mono\"); \n\
    } \n\
 \n\
  if (type == \"combo\") \n\
    save_plot(name=name, type=\"st\", outformat, h=h, w=w) \n\
  else if (type == \"heat\") \n\
    save_plot(name=name, type=\"fi\", outformat, h=h, w=w) \n\
  else  # type == \"contour\" \n\
    save_plot(name=name, type=\"ln\", outformat, h=h, w=w) \n\
} \n\
 \n\
my_plot(hist = args$file, asn = args$assign, name = args$output, h = args$ydim, w = args$xdim, type = args$type, pdf = args$pdf, s = args$source) \n\
";
