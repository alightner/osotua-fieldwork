library(igraph)
library(tidyverse)
library(data.table)
library(Rcpp)
library(microbenchmark)
library(ggsci)
library(gridExtra)

# --- Plots comparing strategic defection vs. standard NBT -----
## pure strategy dyads (psd)
df <- read_csv('abm-datasets/dyad_df.csv')
dk <- read_csv('abm-datasets/dyad_cdf.csv')
dk <- subset(dk, condition != 'esile')
dk$condition <- factor(dk$condition, levels=c('standard', 
  'stingy (1)', 'no exchange', 'greedy (1)', 'greedy (2)'
))
## pairwise plots, part 1 analysis (defection strategies vs. standard NBT)
plscope <- 1000
sub_plscope <- 100

# short-term benefit for defection in NBT
psd1 <- ggplot(dk, aes(x=time, y=cdf, linetype=condition)) + geom_line() +
  xlim(c(0,sub_plscope)) + theme_classic() + labs(y='proportion herds surviving')

# long-term cost for defection NBT
psd2 <- ggplot(dk, aes(x=time, y=cdf, linetype=condition)) + geom_line() +
  xlim(c(0,plscope)) + theme_classic() + labs(y='proportion herds surviving') +
  geom_vline(xintercept = sub_plscope, alpha=0.2) + 
  theme(legend.position = 'none')

# long-term cost for defection NBT (clearer view, log-scaled y-axis)
psd3 <- ggplot(dk, aes(x=time, y=cdf, linetype=condition)) + geom_line() +
  xlim(c(0,plscope)) + scale_y_log10() + theme_classic() + labs(y='proportion herds surviving') +
  geom_vline(xintercept = sub_plscope, alpha=0.2)

# long-term cost for defection NBT (zoomed in by lim. y-axis y=0.5)
# supplementary material plots? ('smplots')
smplot1 <- ggplot(dk, aes(x=time, y=cdf, linetype=condition)) + geom_line() +
  xlim(c(0,plscope)) + theme_classic() + labs(y='proportion herds surviving') +
  geom_rect(xmin=0, xmax=plscope, ymin=0, ymax=0.5, fill=NA, colour='black', linetype=3)
smplot2 <- ggplot(dk, aes(x=time, y=cdf, linetype=condition)) + geom_line() +
  xlim(c(0,plscope)) + ylim(c(0,0.5)) + theme_classic() + labs(y='proportion herds surviving') +
  geom_vline(xintercept = sub_plscope, alpha=0.2)

# 2d density plots 
kde2plot <- function(condx){
  gpl <- ggplot(subset(df, condition==condx), aes(x=t_fin1,y=t_fin2)) + 
    stat_density2d(geom="polygon", aes(fill=..level..)) + 
                #   color="black",linetype=1) + 
    # stat_contour(binwidth=1) +
    scale_x_log10() + scale_y_log10() + theme_classic() + 
    ggtitle(condx) + labs(x='duration, P1', y='duration, P2') +
    guides(fill=FALSE)
  return(gpl)
}
cxlist <- levels(droplevels(data.frame(table(df$condition))$Var1))
cxlist <- cxlist[cxlist!="esile"]
cxlist <- factor(cxlist, 
                 levels=c('standard', 'no exchange', 'stingy (1)', 
                          'greedy (1)', 'greedy (2)'))
plist <- as.list(sort(cxlist))
plist <- lapply(plist, kde2plot)
#nCol <- floor(sqrt(length(plist)))
nCol = 3
topPlot1 <- do.call("grid.arrange", c(plist, ncol=nCol))

# summary stats
summarydf1 <- df %>% subset(condition != 'esile') %>% 
  group_by(condition) %>% 
  summarise(mean(st1), mean(st2),
            median(t_fin1), median(t_fin2))

summarycdf1 <- dk %>% group_by(condition) %>% 
  summarise(t100 = mean(cdf[time==100]),
            t250 = mean(cdf[time==250]),
            t500=mean(cdf[time <= 510 & time >= 490]),
            t1000=mean(cdf[time ==1000]))


