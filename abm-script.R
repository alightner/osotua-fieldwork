library(igraph)
library(tidyverse)
library(data.table)
library(Rcpp)
library(microbenchmark)
library(ggsci)

## model script (cpp)
sourceCpp('osotua-model-script.cpp')

## function for data analysis/ datasets (CDF)
pr_vec <- function(x){    # function for CDF plotting
  x1 <- table(x)/(sum(table(x)))
  y1 <- rep(NA, length(x1))
  incl <- as.numeric(levels(droplevels(as.data.frame(x1)$x)))
  for(i in incl){
    y1[which(incl==i)] <- sum(x1[which(incl==i):length(incl)])
  }
  return(list(time=incl, survival=y1))
}

## global parameters set
t = 10000
duration = 1000
init_st = 70    # min = 64, coded into cpp script
vtrate = 0.1; vt_mean = 0.3; vt_sd = 0.1
gr_mean=0.034; gr_sd=0.0253

# ----- DYAD/PARTNERSHIP MODEL -----
## network (instance) parameters set
g = graph(c(1,2), n=1, directed=FALSE)
#g = graph(c(1,2,2,3,3,1), n=1, directed=FALSE)
n = vcount(g)
adj = lapply(adjacent_vertices(g, V(g)), function(x){x-1})
OsotuaStandard(t, duration, n, init_st, adj, vtrate, vt_mean, vt_sd,
               gr_mean, gr_sd, newfile = TRUE)

# 'esile' was coded out of curiosity, optional to run
# and **only works for dyads right now**
creditor <- unlist(adj)
pr_repay <- c(0,1)
tolerated <- 5
credit_limit <- 500
EsileStandard(t, duration, n, init_st, adj, creditor, vtrate,
              vt_mean, vt_sd, gr_mean, gr_sd, pr_repay, 
              tolerated, credit_limit, newfile=TRUE)

# "stingy" defect strategy
## defect parameters added (NOTE: 1. wdefect index starts at 0, and
## 2. "defect" is stingy/withheld benefits)
wdefect = 0
prdefect = 0.5
OsotuaStingy(t, duration, n, init_st, adj, vtrate, vt_mean, vt_sd,
             gr_mean, gr_sd, wdefect, prdefect, newfile = TRUE)
## pairwise defect (stingy) == control, no exchange occurs
wdefect = c(0,1)
prdefect = 0.5
OsotuaDefect2(t, duration, n, init_st, adj, vtrate, vt_mean, vt_sd,
             gr_mean, gr_sd, wdefect, prdefect, newfile = TRUE)

# "stingy" defect strategy
wdefect = 0
prdefect = 0.5
OsotuaGreedy(t, duration, n, init_st, adj, vtrate, vt_mean, vt_sd,
             gr_mean, gr_sd, wdefect, prdefect, newfile = TRUE)
wdefect = c(0,1)
prdefect = 0.5
OsotuaGreedy2(t, duration, n, init_st, adj, vtrate, vt_mean, vt_sd,
              gr_mean, gr_sd, wdefect, prdefect, newfile = TRUE)

## creating datasets
files <- c('OsotuaStandard.csv',
           'OsotuaStingy.csv',
           'OsotuaDefect2.csv',
           'OsotuaGreedy.csv',
           'OsotuaGreedy2.csv'
           #'EsileStandard.csv'
           )
df <- data.table(rbindlist(lapply(files, read_csv, col_names=FALSE)))
colnames(df) <- c('st1','st2','t_fin1','t_fin2')
df$condition <- c(rep('standard', t),
                  rep('stingy (1)', t),
                  rep('no exchange', t),
                  rep('greedy (1)', t),
                  rep('greedy (2)', t)
                  #rep('esile', t)
                  )
dk <- data.table(
  cdf = c(
    pr_vec(c(df$t_fin1[df$condition=='standard'],df$t_fin2[df$condition=='standard']))$survival,
    pr_vec(c(df$t_fin1[df$condition=='no exchange'],df$t_fin2[df$condition=='no exchange']))$survival,
    pr_vec(df$t_fin1[df$condition=='stingy (1)'])$survival,
    pr_vec(df$t_fin1[df$condition=='greedy (1)'])$survival,
    pr_vec(c(df$t_fin1[df$condition=='greedy (2)'],df$t_fin2[df$condition=='greedy (2)']))$survival
#    pr_vec(c(df$t_fin1[df$condition=='esile'],df$t_fin2[df$condition=='esile']))$survival
  ), time = c(
    pr_vec(c(df$t_fin1[df$condition=='standard'],df$t_fin2[df$condition=='standard']))$time,
    pr_vec(c(df$t_fin1[df$condition=='no exchange'],df$t_fin2[df$condition=='no exchange']))$time,
    pr_vec(df$t_fin1[df$condition=='stingy (1)'])$time,
    pr_vec(df$t_fin1[df$condition=='greedy (1)'])$time,
    pr_vec(c(df$t_fin1[df$condition=='greedy (2)'],df$t_fin2[df$condition=='greedy (2)']))$time
#    pr_vec(c(df$t_fin1[df$condition=='esile'],df$t_fin2[df$condition=='esile']))$time
  ), condition = c(
    rep('standard', length(pr_vec(c(df$t_fin1[df$condition=='standard'],df$t_fin2[df$condition=='standard']))$time)),
    rep('no exchange', length(pr_vec(c(df$t_fin1[df$condition=='no exchange'],df$t_fin2[df$condition=='no exchange']))$time)),
    rep('stingy (1)', length(pr_vec(df$t_fin1[df$condition=='stingy (1)'])$time)),
    rep('greedy (1)', length(pr_vec(df$t_fin1[df$condition=='greedy (1)'])$time)),
    rep('greedy (2)', length(pr_vec(c(df$t_fin1[df$condition=='greedy (2)'],df$t_fin2[df$condition=='greedy (2)']))$time))
#    rep('esile', length(pr_vec(c(df$t_fin1[df$condition=='esile'],df$t_fin2[df$condition=='esile']))$time))
  ))
dk$condition <- factor(dk$condition, levels=c(
  'standard', 
  #'esile', 
  'stingy (1)', 'no exchange', 'greedy (1)', 'greedy (2)'
))
write.table(df, file='dyad_df_pr50.csv', sep=',', row.names = FALSE)
write.table(dk, file='dyad_cdf_pr50.csv', sep=',', row.names = FALSE)

# ----- TRIAD PARTNERSHIP MODEL -------
## global parameters set
t = 10000
duration = 1000
init_st = 70    # min = 64, coded into cpp script
vtrate = 0.1; vt_mean = 0.3; vt_sd = 0.1
gr_mean=0.034; gr_sd=0.0253

# probability of defection
prdefect = 0.5

## network (instance) parameters set
g = graph(c(1,2,2,3,3,1), n=1, directed=FALSE)
n = vcount(g)
adj = lapply(adjacent_vertices(g, V(g)), function(x){x-1})
OsotuaStandard(t, duration, n, init_st, adj, vtrate, vt_mean, vt_sd,
               gr_mean, gr_sd, newfile = TRUE)

# "stingy" defect strategy
## defect parameters added (NOTE: 1. wdefect index starts at 0, and
## 2. "defect" is stingy/withheld benefits)
wdefect = 0
#prdefect = 0.5
OsotuaStingy(t, duration, n, init_st, adj, vtrate, vt_mean, vt_sd,
             gr_mean, gr_sd, wdefect, prdefect, newfile = TRUE)
wdefect = c(0,1)
#prdefect = 0.5
OsotuaStingy(t, duration, n, init_st, adj, vtrate, vt_mean, vt_sd,
             gr_mean, gr_sd, wdefect, prdefect, newfile = FALSE)
## pairwise defect (stingy) == control, no exchange occurs
wdefect = c(0,1,2)
#prdefect = 0.5
OsotuaDefect2(t, duration, n, init_st, adj, vtrate, vt_mean, vt_sd,
              gr_mean, gr_sd, wdefect, prdefect, newfile = TRUE)

# "greedy" defect strategy
wdefect = 0
#prdefect = 0.5
OsotuaGreedy(t, duration, n, init_st, adj, vtrate, vt_mean, vt_sd,
             gr_mean, gr_sd, wdefect, prdefect, newfile = TRUE)
wdefect = c(0,1)
#prdefect = 0.5
OsotuaGreedy(t, duration, n, init_st, adj, vtrate, vt_mean, vt_sd,
              gr_mean, gr_sd, wdefect, prdefect, newfile = FALSE)
wdefect = c(0,1,2)
#prdefect = 0.5
OsotuaGreedy(t, duration, n, init_st, adj, vtrate, vt_mean, vt_sd,
             gr_mean, gr_sd, wdefect, prdefect, newfile = FALSE)

## creating datasets
files <- c('OsotuaStandard.csv',
           'OsotuaStingy.csv',
           'OsotuaDefect2.csv',
           'OsotuaGreedy.csv'
           #'OsotuaGreedy2.csv'
           #'EsileStandard.csv'
           )
df <- data.table(rbindlist(lapply(files, read_csv, col_names=FALSE)))
colnames(df) <- c('st1','st2', 'st3','t_fin1','t_fin2', 't_fin3')
df$condition <- c(rep('standard', t),
                  rep('stingy (1)', t),
                  rep('stingy (2)', t),
                  rep('no exchange', t),
                  rep('greedy (1)', t),
                  rep('greedy (2)', t),
                  rep('greedy (3)', t)
                  #rep('esile', t)
                  )
dk <- data.table(
  cdf = c(
    pr_vec(c(df$t_fin1[df$condition=='standard'],
             df$t_fin2[df$condition=='standard'],
             df$t_fin3[df$condition=='standard']))$survival,
    pr_vec(c(df$t_fin1[df$condition=='no exchange'],
             df$t_fin2[df$condition=='no exchange'],
             df$t_fin3[df$condition=='no exchange']))$survival,
    pr_vec(df$t_fin1[df$condition=='stingy (1)'])$survival,
    pr_vec(c(df$t_fin1[df$condition=='stingy (2)'],
             df$t_fin2[df$condition=='stingy (2)']))$survival,
    pr_vec(df$t_fin1[df$condition=='greedy (1)'])$survival,
    pr_vec(c(df$t_fin1[df$condition=='greedy (2)'],
             df$t_fin2[df$condition=='greedy (2)']))$survival,
    pr_vec(c(df$t_fin1[df$condition=='greedy (3)'],
             df$t_fin2[df$condition=='greedy (3)'],
             df$t_fin3[df$condition=='greedy (3)']))$survival
  #  pr_vec(c(df$t_fin1[df$condition=='esile'],df$t_fin2[df$condition=='esile']))$survival
  ), time = c(
    pr_vec(c(df$t_fin1[df$condition=='standard'],
             df$t_fin2[df$condition=='standard'],
             df$t_fin3[df$condition=='standard']))$time,
    pr_vec(c(df$t_fin1[df$condition=='no exchange'],
             df$t_fin2[df$condition=='no exchange'],
             df$t_fin3[df$condition=='no exchange']))$time,
    pr_vec(df$t_fin1[df$condition=='stingy (1)'])$time,
    pr_vec(c(df$t_fin1[df$condition=='stingy (2)'],
             df$t_fin2[df$condition=='stingy (2)']))$time,
    pr_vec(df$t_fin1[df$condition=='greedy (1)'])$time,
    pr_vec(c(df$t_fin1[df$condition=='greedy (2)'],
             df$t_fin2[df$condition=='greedy (2)']))$time,
    pr_vec(c(df$t_fin1[df$condition=='greedy (3)'],
             df$t_fin2[df$condition=='greedy (3)'],
             df$t_fin3[df$condition=='greedy (3)']))$time
  ), condition = c(
    rep('standard', length(
      pr_vec(c(df$t_fin1[df$condition=='standard'],
               df$t_fin2[df$condition=='standard'],
               df$t_fin3[df$condition=='standard']))$time
    )),
    rep('no exchange', length(
      pr_vec(c(df$t_fin1[df$condition=='no exchange'],
               df$t_fin2[df$condition=='no exchange'],
               df$t_fin3[df$condition=='no exchange']))$time
    )),
    rep('stingy (1)', length(
      pr_vec(df$t_fin1[df$condition=='stingy (1)'])$time
    )),
    rep('stingy (2)', length(
      pr_vec(c(df$t_fin1[df$condition=='stingy (2)'],
               df$t_fin2[df$condition=='stingy (2)']))$time
    )),
    rep('greedy (1)', length(
      pr_vec(df$t_fin1[df$condition=='greedy (1)'])$time
    )),
    rep('greedy (2)', length(
      pr_vec(c(df$t_fin1[df$condition=='greedy (2)'],
               df$t_fin2[df$condition=='greedy (2)']))$time
    )),
    rep('greedy (3)', length(
      pr_vec(c(df$t_fin1[df$condition=='greedy (3)'],
               df$t_fin2[df$condition=='greedy (3)'],
               df$t_fin3[df$condition=='greedy (3)']))$time
    ))
  ))
dk$condition <- factor(dk$condition, levels=c(
  'standard', 'stingy (1)', 'stingy (2)',  'no exchange', 'greedy (1)', 'greedy (2)', 'greedy (3)'
))

write.table(df, file='triad_df_pr50.csv', sep=',', row.names = FALSE)
write.table(dk, file='triad_cdf_pr50.csv', sep=',', row.names = FALSE)


