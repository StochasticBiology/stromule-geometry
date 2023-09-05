library(ggplot2)

# read in output from larger-scale simulation
df = read.csv("geom-snaps-vmany.csv")
df.grid = read.csv("vmany-grid.csv")

# various visualisation options. 
# 1 -- just distance to plastid
# 2 -- montage of distance and plastid access
# 3 -- just plastid access
vis = 1
if(vis == 1) {
  df.grid$plottable =-df.grid$ir/max(df.grid$ir)
  midpoint = mean(df.grid$plottable)-0.05
} else if(vis == 2) {
  df.grid$plottable = ifelse(df.grid$y > 0, 0.3-df.grid$ir/max(df.grid$ir) , df.grid$pa/max(df.grid$pa))
  midpoint = 0.3
} else if(vis == 3) {
  df.grid$plottable = df.grid$pa/max(df.grid$pa)
  midpoint = mean(df.grid$plottable)
} 
  
# plot size and colour palette
dim = 2*max(df.grid$y)
mycol0 = "#AAAAAA"
mycol = "#000000"
mycol1 = "#00CCCC"
mycol2 = "#000000"
mycol3 = "#CEB300"
# "shadow" offset
mdx=-0.5
mdy=0.5
g1 = ggplot() + geom_tile(data=df.grid, aes(x=x, y=y, fill=plottable)) +                          # colour map
  geom_segment(data=df, aes(x=x1,y=y1,xend=x2,yend=y2), color=mycol, linewidth=3, alpha=0.5) +    # shadow of stromule
  geom_point(data=df, aes(x=x1, y=y1), size=4, color=mycol, alpha=0.5, stroke=NA) +               # shadow of plastid body
  geom_segment(data=df, aes(x=x1+mdx,y=y1+mdy,xend=x2+mdx,yend=y2+mdy), color=mycol0) +           # stromule
  geom_point(data=df, aes(x=x1+mdx, y=y1+mdy), size=3, color=mycol0) +                            # plastid body
  coord_cartesian(xlim=c(-dim/2,dim/2),ylim=c(-dim/2,dim/2)) + 
  theme_light() + scale_fill_gradient2(low=mycol1, mid=mycol2, high=mycol3, midpoint=midpoint) + 
  theme_void() + theme(legend.position="none") 

# output to screen and file
g1
sf = 3
png(paste(c("view-", vis, ".png"), collapse=""), width=800*sf, height=800*sf, res=72*sf)
g1
dev.off()
