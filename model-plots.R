# visualisations from simulation and plots of geometric models

library(ggplot2)
library(ggforce)
library(gridExtra)

############## simulation snapshot visualisations

# read in arrangements of stromule segments
df = read.csv("geom-snaps.csv")
sub = df[df$sample==0,]

# loop through experiments (stromule structures) and inter-plastid distances, retrieving the simulation output for each and making a plot
# units are 2um / pixel
plot.ir = plot.pa = list()
counter = 0
for(expt in c(0,1,2,3)) {
  for(dist in c(0, 10, 20, 30, 40, 50)) {
    # plastid positions
    plastids = data.frame(x0 = c(0, dist), y0 = c(0, 0), r = c(4, 4))
    counter = counter+1
    # read in corresponding file
  fname = paste(c("cellmodelmany-", expt, "-", dist, ".txt"), collapse="")
  pdf = read.csv(fname, sep=" ", header=F)
  colnames(pdf) = c("x", "y", "IR", "PA")
  if(counter %in% c(1, 7, 13, 19)) {
    x1 = -20; x2 = 40; y1 = -20; y2 = 20
  } else {
    x1 = -20; x2 = 60; y1 = -30; y2 = 30
  }
  # interaction range plot
  plot.ir[[counter]] = ggplot() + 
    geom_raster(data=pdf, aes(x=x-100,y=y-100,fill=IR)) + scale_fill_gradientn(limits = c(0,38), colours = terrain.colors(10))  + 
    geom_segment(data = sub[sub$expt==expt & sub$dist==dist,], aes(x=x1,y=y1,xend=x2,yend=y2)) +
    geom_circle(data= plastids, aes(x0=x0, y0=y0, r=r), fill="#448844") + 
    theme_void() + xlim(x1,x2) + ylim(y1,y2) +
    theme(legend.position="none")
  # plastid access plot
  plot.pa[[counter]] = ggplot() + 
     geom_raster(data=pdf, aes(x=x-100,y=y-100,fill=PA)) + scale_fill_gradientn(limits = c(0,0.77), colours = terrain.colors(10))  + 
    geom_segment(data = sub[sub$expt==expt & sub$dist==dist,], aes(x=x1,y=y1,xend=x2,yend=y2)) +
    geom_circle(data= plastids, aes(x0=x0, y0=y0, r=r), fill="#448844") + 
    theme_void() + xlim(x1,x2) + ylim(y1,y2) +
    theme(legend.position="none")
  }
}

# example plots for all conditions
sf = 2
png("si-fig-ir.png", width=1000*sf, height=500*sf, res=72*sf)
grid.arrange(grobs = plot.ir, ncol=6)
dev.off()
png("si-fig-pa.png", width=1000*sf, height=500*sf, res=72*sf)
grid.arrange(grobs = plot.pa, ncol=6)
dev.off()

# a couple of specific example plots
png("pair-plastids.png", width=500*sf, height=150*sf, res=72*sf)
g.p.1 = plot.ir[[14]] + xlim(-30, 60) + ylim(-40, 40)
g.p.2 = plot.ir[[17]] + xlim(-30, 60) + ylim(-40, 40)
grid.arrange(g.p.1, g.p.2, nrow=1)
dev.off()

# single-plastid cases
png("single-plastid.png", width=600*sf, height=200*sf, res=72*sf)
grid.arrange(plot.ir[[1]], plot.ir[[7]], plot.ir[[13]], plot.ir[[19]], 
             plot.pa[[1]], plot.pa[[7]], plot.pa[[13]], plot.pa[[19]], nrow=2)
dev.off()

############## simulation statistics

## read in simulation histograms
df = read.csv("hists-many.csv")
# label experiments
df$Stromule = ""
df$Stromule[df$expt==0] = "None"
df$Stromule[df$expt==1] = "Single"
df$Stromule[df$expt==2] = "Branched"
df$Stromule[df$expt==3] = "Double"
df$Stromule = factor(df$Stromule, levels = c("None", "Single", "Double", "Branched"))
# set length scale to um
df$dist = df$dist/2

# plots as a function of distance, faceted by inter-plastid distance
g1 = ggplot(df[df$hist==1,], aes(x=val*50, y=n, color=Stromule)) + 
  geom_line() + xlim(0,20) + theme_classic() +
  xlab("Distance from plastid") + ylab("Area") + facet_wrap(~dist, nrow=1)
g2 = ggplot(df[df$hist==2 & df$val != 0,], aes(x=val, y=n, color=Stromule)) + 
  geom_line() + theme_classic() +
  xlab("Angle fraction to plastid") + ylab("Area") + facet_wrap(~dist, nrow=1)

png("si-fig-hists.png", width=1000*sf, height=300*sf, res=72*sf)
grid.arrange(g1, g2, nrow=2)
dev.off()

# plots for fixed distance
g1.a = ggplot(df[df$hist==1 & df$val == 0,], aes(x=dist, y=n, fill=Stromule)) + 
  geom_col(position="dodge") +
  theme_classic() + xlab("Distance between plastids / µm") + ylab("Interaction area / AU")
g2.a = ggplot(df[df$hist==2 & df$val > 0.5,], aes(x=dist, y=n, fill=Stromule)) + 
  geom_col(position="dodge") + 
  theme_classic() + xlab("Distance between plastids / µm") + ylab("Plastid access / AU")
#grid.arrange(g1.a, g2.a, nrow=2)

# plots as a function of distance, just for single plastid
g1.b = ggplot(df[df$hist==1 & df$dist==0,], aes(x=val*50, y=n, color=Stromule)) + 
  geom_line() + xlim(0,20) + theme_classic() +
  xlab("Distance from plastid") + ylab("Area") 
g2.b = ggplot(df[df$hist==2 & df$dist==0 & df$val != 0,], aes(x=val, y=n, color=Stromule)) + 
  geom_line() + theme_classic() +
  xlab("Angle fraction to plastid") + ylab("Area") 

png("fig-hists.png", width=600*sf, height=260*sf, res=72*sf)
grid.arrange(g1.b, g1.a, g2.b, g2.a, nrow=2)
dev.off()

############## geometric models

# set physical quantities
r0 = 2
l = 10
h = 0.5
rho = 0.25
r = sqrt(r0**2 - rho*l/2)

# area as a function of distance d and stromule length l, for conserved (i) area (ii) volume
a.i = function(l, d) {
  return(pi*d*(d + sqrt(4*r0**2 - 2*l*rho)) + pi/2 * (d+rho)**2 + 2*d*l - 2*d**2)
}
a.ii = function(l, d) {
  return(pi*d*(d + (8*r0**3 - 6*l*rho**2)**(1/3)) + pi/2 * (d+rho)**2 + 2*d*l - 2*d**2)
}

eqs = data.frame()
eqs = rbind(eqs, data.frame(l=0:10, expt="Linear area", val=a.i(0:10, 1)))
eqs = rbind(eqs, data.frame(l=0:10, expt="Linear vol", val=a.ii(0:10, 1)))
eqs = rbind(eqs, data.frame(l=0:10, expt="Linear area ratio", val=a.i(0:10, 1)/a.i(0, 1)))
eqs = rbind(eqs, data.frame(l=0:10, expt="Linear vol ratio", val=a.ii(0:10, 1)/a.ii(0, 1)))

# plot
eq.g.1 = ggplot(eqs[eqs$expt %in% c("Linear area ratio", "Linear vol ratio"),], aes(x=l, y=val, color=expt)) + 
  geom_line() + theme_classic()

# simplified model: region with/without stromule in 2D and 3D
a.2d.no = function(d) {
  return(pi*(r0+d)**2 - pi*r0**2)
}
a.2d.yes= function(d) {
  return(pi*(r+d)**2 - pi*r**2 - 2*d**2 + 2*l*d + (1/2)*pi*d**2)
}
v.3d.no = function(d) {
  return( pi*(r0+d)**2*h - pi*r0**2*h +2*pi*(r0+d)**2*d )
  #return( (4/3)*pi*(r0+d)**3 - (4/3)*pi*r0**3 )
}
v.3d.yes = function(d) {
  return( pi*(r+d)**2*h - pi*r**2*h  +2*pi*(r+d)**2*d + pi*d**2*(l-d) +(2/3)*pi*d**3 )
  #return( (4/3)*pi*(r+d)**3 - (4/3)*pi*r**3 - pi*d**3 + pi*d**2*l + (2/3)*pi*d**3 )
}

eqs = data.frame()
eqs = rbind(eqs, data.frame(d=0:100, expt="2D no", val=a.2d.no(0:100)))
eqs = rbind(eqs, data.frame(d=0:100, expt="2D yes", val=a.2d.yes(0:100)))
eqs = rbind(eqs, data.frame(d=0:100, expt="2D ratio", val=a.2d.yes(0:100)/a.2d.no(0:100)))

eqs = rbind(eqs, data.frame(d=0:100, expt="3D no", val=v.3d.no(0:100)))
eqs = rbind(eqs, data.frame(d=0:100, expt="3D yes", val=v.3d.yes(0:100)))
eqs = rbind(eqs, data.frame(d=0:100, expt="3D ratio", val=v.3d.yes(0:100)/v.3d.no(0:100)))

# plot
eq.g.2 = ggplot(eqs[eqs$expt %in% c("2D ratio", "3D ratio"),], aes(x=d, y=val, color=expt)) +
  geom_line() + theme_classic()

# produce figure
png("si-geom-plots.png", width=600, height=300)
grid.arrange(eq.g.1, eq.g.2, nrow=1)
dev.off()

## overlap calculations
# "beak" -- the overlapping beak-shaped region between shells of width d surrounding plastids impinging on a point with angle alpha
beak = function(d, alpha) {
  return( 0.5*d**2 * (sin(alpha) + 0.5*(sin(alpha)**2)/(sin(alpha/2)**2) * tan((pi-alpha)/2) - pi + alpha) )
}

# transformation functions for use in optimisation (enforcing a given parameter range)
transfun = function(x) {
  return(1/(1+exp(-x)))
}
invtransfun = function(x) {
  return(log(x / (1-x)))
}

# total overlapping region if we have three segments impinging on a point
# most of the code is bookkeeping to reflect the relative angles in a consistent format
total.beak = function(theta, r=1, verbose=FALSE) {
  alpha.1 = 2*pi*transfun(theta[1])
  alpha.2 = 2*pi*transfun(theta[2])
  
  if(alpha.1 > pi) { 
    d.1 = 2*pi - alpha.1 
  } else {
    d.1 = alpha.1
  }
  if(alpha.2 > pi) { 
    d.2 = 2*pi - alpha.2
  } else {
    d.2 = alpha.2
  }
  if(d.1 > d.2) {  # the angle to line 1 is greater than the angle to line 2
    bigger = 1
    if((alpha.1 > pi & alpha.2 > alpha.1) | (alpha.1 < pi & alpha.2 < alpha.1)) { # line 2 falls inside the smaller angle to line 1
      inside = TRUE
    } else {
      inside = FALSE
    }
  } else { # the angle to line 2 is greater than the angle to line 1
    bigger = 2
    if((alpha.2 > pi & alpha.1 > alpha.2) | (alpha.2 < pi & alpha.1 < alpha.2)) { # line 1 falls inside the smaller angle to line 2
      inside = TRUE
    } else {
      inside = FALSE
    }
  }
  if(alpha.1 > alpha.2) {
    d.12 = alpha.1 - alpha.2
  } else {
    d.12 = alpha.2 - alpha.1
  }
  if(d.12 > pi) {
    d.12 = d.12 - pi
  }
  
  if(verbose == TRUE) {
    print(paste(c("Arguments ", alpha.1, " ", alpha.2, " Angles ", d.1, " and ", d.2, " between ", d.12, " bigger is ", bigger, " inside is ", inside), collapse=""))
  }
  if(inside == TRUE) {
    if(bigger == 1) {
      return( beak(r, d.2) + beak(r, d.12) - beak(r, d.1) )
    } else {
      return( beak(r, d.1) + beak(r, d.12) - beak(r, d.2) )
    }
  } else {
    return( beak(r, d.1) + beak(r, d.12) + beak(r, d.2) )
  }
}

# optimise (minimise) the total overlap area as a function of the two governing angles
best = optim(c(1,1.1), total.beak)
alpha.1 = 2*pi*transfun(best$par[1])
alpha.2 = 2*pi*transfun(best$par[2])
total.beak(best$par, verbose=TRUE)

# build dataframe storing the total overlap as a function of the two angles
beaks = data.frame()
for(a1 in 0:62) {
  for(a2 in 0:62) {
    beaks = rbind(beaks, data.frame(alpha1=a1/10, alpha2=a2/10, A=total.beak(c(invtransfun((a1/10)/(2*pi)),invtransfun((a2/10)/(2*pi))))))
  }
}

# plot this overlap 
png("si-overlap-plot.png", width=800*sf, height=600*sf, res=72*sf)
ggplot(beaks[beaks$alpha1 < pi & beaks$alpha2 < beaks$alpha1+pi,], aes(x=alpha1,y=alpha2,fill=log(A))) + 
  geom_tile() + scale_fill_gradientn(limits = c(-2,4), colours = terrain.colors(10, rev=TRUE) ) +
  theme_classic()
dev.off()

