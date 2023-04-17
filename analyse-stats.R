# data analysis for physical measurements of plastids
library(ggplot2)
library(readxl)
library(gridExtra)

# various likelihood functions for model fits
lnorm.neg.llik = function(theta, data) {
  return(-sum(dlnorm(data, mean=theta[1], sd=abs(theta[2]), log=TRUE)))
}
gamma.neg.llik = function(theta, data) {
  return(-sum(dgamma(data, shape=abs(theta[1]), rate=abs(theta[2]), log=TRUE)))
}
poisson.neg.llik = function(theta, data) {
  return(-sum(dpois(data, lambda=theta, log=TRUE)))
}
poisson.neg.llik.counts = function(theta, data) {
  llik = 0
  for(k in 0:(length(data)-1)) {
    llik = llik + data[k+1]* (k*log(theta) - theta - log(factorial(k))) 
  }
  return(-llik)
}

################### stromule length
#### stressed case
din = read_excel("Data/stromule length.xlsx")
df = data.frame(din)
lens = df$total.length.in.µm
lens = lens[!is.na(lens)]

# bin the observations
lens.bins = c()
for(i in 0:35) {
  lens.bins = c(lens.bins, length(which(lens >= i & lens < i+1))/length(lens))
}
# fit a log-normal distribution and produce a dataframe reflecting the fitted distribution
lnormfit = nlm(lnorm.neg.llik, c(1,1), data=lens)
lens.est = dlnorm(0:35+0.5, mean=lnormfit$estimate[1], sd=abs(lnormfit$estimate[2]))
lens.df = data.frame(l = rep(0:35, 2), val=c(rep("Observed", 36), rep("Estimate", 36)), p = c(lens.bins, lens.est))

# gather various statistics
mn.est = mean(log(lens))
sn.est = sd(log(lens))
m.est = exp(mn.est)
hi.est = exp(mn.est+1.96*sn.est)
lo.est = exp(mn.est-1.96*sn.est)

# plot observed and fitted distributions
g.1 = ggplot(lens.df) + 
  geom_col(aes(x=l, y=p, fill=val), position="dodge") + 
  theme_classic() + xlab("Stromule length / µm") + ylab("Frequency") 

#### unstressed case
din2 = read_excel("Data/stromule length-NI.xlsx")
df = data.frame(din2)
lens.ni = df$total.length.in.µm

# bin the observations
lens.bins.ni = c()
for(i in 0:35) {
  lens.bins.ni = c(lens.bins.ni, length(which(lens.ni >= i & lens.ni < i+1))/length(lens.ni))
}
# fit a log-normal distribution and produce a dataframe reflecting the fitted distribution
lnormfit.ni = nlm(lnorm.neg.llik, c(1,1), data=lens.ni)
lens.est.ni = dlnorm(0:35+0.5, mean=lnormfit.ni$estimate[1], sd=abs(lnormfit.ni$estimate[2]))
lens.df.ni = data.frame(l = rep(0:35, 2), val=c(rep("Observed", 36), rep("Estimate", 36)), p = c(lens.bins.ni, lens.est.ni))

# plot observed and fitted distributions
g.1a = ggplot(lens.df.ni) + 
  geom_col(aes(x=l, y=p, fill=val), position="dodge") + 
  theme_classic() + xlab("Stromule length / µm") + ylab("Frequency") 

## combined stressed/unstressed
# build and plot dataframe combining the (binned) stressed and unstressed observations
lens.both = data.frame(l = c(0:35, 0:35), expt=c(rep("Stress", 36), rep("Control", 36)), p = c(lens.bins, lens.bins.ni))
g.1b = ggplot(lens.both) + 
  geom_col(aes(x=l, y=p, fill=expt), position="dodge") + 
  theme_classic() + xlab("Stromule length / µm") + ylab("Frequency")

# build and plot dataframe combining the (individual) stressed and unstressed observations
lens.ind.df = data.frame(label = c(rep("Unstressed", length(lens.ni)), rep("Stressed", length(lens))), val=c(lens.ni, lens))
#ggplot(lens.ind.df, aes(x=label, y=val)) + geom_jitter() + geom_boxplot()
#ggplot(lens.ind.df, aes(x=label, y=log(val))) + geom_jitter() + geom_boxplot()
g.1c = ggplot(lens.ind.df, aes(x=label, y=val)) + 
  geom_jitter(width=0.1, size = 0.5, alpha = 0.5) + 
  geom_boxplot(width=0.25, outlier.shape=NA, position=position_nudge(x=-0.25)) +
  scale_y_continuous(trans = "log10", breaks=c(0.5, 1, 2, 4, 8, 16, 32)) +
  theme_classic() + xlab("Experiment") + ylab("Length / µm")

################### cytoplasm thickness
din = read_excel("Data/cytoplasm thickness.xlsx", sheet="Tabelle1")
df = data.frame(din)
# extract and plot statistics
thickness.At = df$cytoplasm.thickness.in.A..thaliana.in.nm
thickness.Nb = df$cytoplasm.thickness.in.N..benthamiana.in.nm
thick.df = data.frame(species = c(rep("Arabidopsis", length(thickness.At)), rep("Nicotiana", length(thickness.Nb))), val = c(thickness.At, thickness.Nb))
g.2 = ggplot(thick.df, aes(x=species, y=val)) + 
  geom_jitter(width=0.1, size = 0.5, alpha = 0.1) + 
  geom_boxplot(width=0.25, outlier.shape=NA, position=position_nudge(x=-0.25)) + 
  theme_classic() + xlab("Species") + ylab("Cytoplasm\nthickness / nm")

################### plastid size
din = read_excel("Data/plastid size.xlsx")
df = data.frame(din)
# extract and plot statistics
lengths = df$length.in.µm
widths = df$width.in.µm
heights = df$hight.in.µm
size.df = data.frame(measurement = c(rep("length", length(lengths)), rep("width", length(widths)), rep("height", length(heights))), val = c(lengths, widths, heights))
g.3 = ggplot(size.df, aes(x=measurement, y=val)) + 
  geom_jitter(width=0.1, size = 0.5, alpha = 0.1) + 
  geom_boxplot(width=0.25, outlier.shape=NA, position=position_nudge(x=-0.25)) +
  theme_classic() + xlab("Plastid dimension") + ylab("Size / µm")

################### plastid density
din = read_excel("Data/plastid density-final sizes.xlsx", skip=4)
df = data.frame(din)
xs = seq(200, 15000, by=200)

# linear model in log space
my.lm = lm(log(df$plastid.density...pl.per.µm.2) ~ log(df$cell.area.in.µm.2))
my.lm.stats = summary(my.lm)
R2 = my.lm.stats$r.squared
# produce plot summarising fit in log space
coef.a = exp(my.lm.stats$coefficients[1,1])
coef.b = my.lm.stats$coefficients[2,1]
plot.label = paste(c("Density = ", round(coef.a, digits=2), 
                     " Area ^ ", round(coef.b, digits=2), 
                     "\nR² = ", round(R2, digits=3)), collapse="")
g.4a = ggplot(df, aes(x=log(cell.area.in.µm.2), y=log(plastid.density...pl.per.µm.2))) + 
  geom_point() + geom_smooth(method="lm", color="#000000") +
  theme_classic() + xlab("log( Cell area / µm² )") + ylab("log( Plastid density / µm⁻² )") +
  geom_text(data=data.frame(x=7,y=-5.5,label=plot.label), aes(x=x, y=y, label=label))

# overlay fitted line on plot of data
g.4 = ggplot(df, aes(x=cell.area.in.µm.2, y=plastid.density...pl.per.µm.2)) + 
  geom_point() + geom_line(data=data.frame(x=xs, y=coef.a*xs**coef.b), aes(x=x, y=y)) +
  theme_classic() + xlab("Cell area / µm²") + ylab("Plastid density / µm⁻²")


################### stromule counts by plastid
din = read_excel("Data/stromules per plastid.xlsx", skip=1)
df = data.frame(din)
df = df[!is.na(df$all.plastids) & !is.na(df$...1),]

# convert entries to numbers and sum
mysum = function(l) {
  return(sum(as.numeric(l)))
}
# build dataframe combining all observations of stromule count
sum.df = data.frame()
sum.df = data.frame(label="Stressed", k = 0:4, val=c(mysum(df[1:5,3]), mysum(df[1:5,4]), mysum(df[1:5,5]), mysum(df[1:5,6]), mysum(df[1:5,7])))
sum.df = rbind(sum.df, data.frame(label="Unstressed", k = 0:4, val=c(mysum(df[6:10,3]), mysum(df[6:10,4]), mysum(df[6:10,5]), mysum(df[6:10,6]), mysum(df[6:10,7]))))
sum.df$val[sum.df$label == "Stressed"] = sum.df$val[sum.df$label == "Stressed"]/sum(sum.df$val[sum.df$label == "Stressed"])
sum.df$val[sum.df$label == "Unstressed"] = sum.df$val[sum.df$label == "Unstressed"]/sum(sum.df$val[sum.df$label == "Unstressed"])
# plot
g.5a = ggplot(sum.df, aes(x=k, y=val, fill=label)) + 
  geom_col(position="dodge", width=0.5) +
  theme_classic() + xlab("Number of stromules") + ylab("Probability")

# loop through different observations
g.5 = list()
for(i in 1:nrow(df)) {
  # extract counts
  counts.by.stromules = as.numeric(df[i,3:7])
  # fit poisson distribution
  pois.fit = nlm(poisson.neg.llik.counts, 1, data=counts.by.stromules)
  print(pois.fit$estimate)
  # produce dataframe from fitted distribution
  pois.est = dpois(0:4, pois.fit$estimate)
  pois.df = data.frame(k = rep(0:4, 2), val = c(rep("Observed", 5), rep("Estimate", 5)), p = c(pois.est, counts.by.stromules/sum(counts.by.stromules)))
  # plot estimated and observed distributions
  g.5[[i]] = ggplot(pois.df) + 
    geom_col(aes(x=k, y=p, fill=val), position="dodge") + 
    theme_classic() + xlab("Stromules per plastid") + ylab("Probability")
}

sf = 2
png("main-plot-stats.png", width=1000*sf, height=700*sf, res=72*sf)
grid.arrange(grid.arrange(g.2, g.3, g.4, nrow=1), grid.arrange(g.1c, g.5a, nrow=1), nrow=2)
dev.off()

png("si-plot-stats.png", width=1000*sf, height=700*sf, res=72*sf)
grid.arrange(g.4a, g.1, g.1a, g.5[[1]], g.5[[6]], nrow=2)
dev.off()
