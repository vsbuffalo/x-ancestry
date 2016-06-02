library(ggplot2)
library(scales)
library(dplyr)
library(purrr)
library(reshape2)
library(wesanderson)
library(tidyr)
library(readr)
library(grid)
library(RColorBrewer)

source("dist-functions.R")

# EDA/draft plots were done in ggplot2, the final publication plots were done
# base graphics. I've kept the original ggplot2 code here for record.

LINESIZE=1.2
GRAPH_DIR = "paper/images/"

graphdir <- function(file, dev='eps') 
  paste(file.path(GRAPH_DIR, file), dev, sep='.')

gsave <- function(p, file, dev='eps') {
  # wrapper to save an image; dev is device and extension.
  # don't include extension.
  ggsave(filename=graphdir(file, dev), plot=p, device='eps')
}

simple_theme <- function(size=20) {
    theme(legend.background = element_blank(), legend.key = element_blank(),
          panel.background = element_blank(), panel.border = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(size=0.1, colour="black"),
          # font sizes
          strip.text.x = element_text(size=size),
          axis.text=element_text(size=size),
          legend.text=element_text(size=size),
          axis.title=element_text(size=size + 2),
          # spacing
          axis.text.x = element_text(vjust=-3),
          axis.title.x = element_text(vjust=-3))
}

## ====== Autosomes Results
# these are primarily for testing the simulation routine, not for the paper

# we compare the prototype Python sim routine (which has a bug: incorrectly
# enforces dependencies between independely assorting chromosomes) to the newer C version.
# The affect of these incorrect dependencies can be seen in these plots.
dfap <- counts2prob(read_tsv("data/auto-shared-ancestor-blockcounts-py.txt"))[, c(3, 1, 2, 4)]
dfa <- counts2prob(read_tsv("data/auto-shared-ancestor-blockcounts.txt"))

# The Python routine has a defaultdict(list) of sims; the
# get_blocks_across_gen() function grabs the gen^th - 1 entry for generation
# gen. This is to to follow the convetion mentioned in the paper.  We take the
# generation in the tab-delimited file (which is one *ahead* in the sims due to
# this offsetting) and subtract it to re-match up the sims and their generation.
dfap$gen <- dfap$gen - 1 
dfa$type = 'c'
dfap$type = 'py'
dfa <- rbind(dfa, dfap)
# now, put the data *back* into the convention (add one generation), since this is what the 
# analytic results in the paper describe.
dfa$gen <- dfa$gen + 1
p <- ggplot(dfa) + geom_point(aes(x=nblocks,y=prob,color=type)) + facet_wrap(~gen)
dfa_analytic <- fit_analytic(0:40, 1:10, auto_blocknum_pmf) 
p <- p + geom_line(dfa_analytic, mapping=aes(nblocks, prob), color='blue', size=0.8)
p <- p + facet_wrap(~gen, scales="free_y", ncol=2)
p

# autosome block length distributions
dfal <- read_tsv("data/auto-shared-ancestor-blocklens.txt")
dfal_analytic <- fit_analytic(seq(0, 2.5, 0.0001), 1:10, auto_blocklen_pdf, xname='lens')
p <- ggplot(dfal) + geom_histogram(aes(x=lens, y=..density..), fill="gray67") 
p <- p + geom_line(data=dfal_analytic, mapping=aes(x=lens, y=prob), color='blue')
p <- p + facet_wrap(~gen, scales="free")
p


## ====== X Chromsome Results
# NOTE:

# A guide to interpretting the various - 1 and + 1 around this code.  The
# original Python prototype simulation procedure followed the convention of
# saying that blocks were counted in the *offspring* not the parent.  Thus,
# when the blocks for the kth generation were counted, the were counting the
# actual sim results from the (k-1)th generation.
#
# To prevent confusion, the new simulation procedure outputs blocks as they are
# in the individuals in the kth generation. So to follow this same convention
# we need to add 1 to the simulation results. Note that the analytic results
# are always in terms of counting the blocks in the offspring, not the parents
# (see the k-1 in the math)
#
# To align the two sims together, we take gen' = gen - 1, such that the label
# in the Python sims is matched up as the blocks in the individual, e.g. the
# original Python sims take the 5th gen individuals and label them as fourth;
# here the block counts labelled 4th generation (from the 3rd generation sim
# individuals) are relabelled as 4-1=3rd generation. 
#
# Then, all simulations are converted *back* to gen' = gen + 1, such that
# blocks *are* counted in the offspring. These then match the analytic.
# results for the kth generation.
#
# All of this is done so the offsetting is done in the analysis side, not the
# simulation side.

## Below are the results comparing C/Python routines.
dfx <- counts2prob(read_tsv("data/x-shared-ancestor-blockcounts.txt"))
dfxp <- counts2prob(read_tsv("data/x-shared-ancestor-blockcounts-py.txt"))[, c(3, 1, 2, 4)]
dfxp$gen <- dfxp$gen - 1  # this puts the python sims in generation notation without convention
dfxc <- dfx
dfxc$type <- 'c'
dfxp$type <- 'py'
dfxc <- rbind(dfxc, dfxp)
dfxc$gen <- dfxc$gen + 1
dfxc_analytic <- fit_analytic(0:6, 1:11, x_blocknum_pmf) 

## Diagnostic plot for C/Python routine comparisons
p <- ggplot(dfxc) + geom_point(aes(x=nblocks, y=prob, color=type), size=0.8, alpha=0.7)
p <- p + facet_wrap(~gen)
# p <- p + geom_line(dfxc_analytic, mapping=aes(nblocks, prob), color='blue'), alpha=0.4
p

## X chromosome segment number figure for paper
dfx <- counts2prob(read_tsv("data/x-shared-ancestor-blockcounts.txt"))
dfx_analytic <- fit_analytic(0:7, 1:7, x_blocknum_pmf) 
dfx_analytic2 <- fit_analytic(0:7, 1:7, x_blocknum_thinned_pmf) 
dfx$gen <- dfx$gen + 1 # count blocks in offspring
p <- ggplot(dfx %>% filter(gen < 10)) + geom_point(aes(x=nblocks, y=prob), color='gray48')
p <- p + facet_wrap(~gen, nrow=2)
p <- p + geom_line(dfx_analytic %>% filter(gen > 1, gen < 10), mapping=aes(nblocks, prob), color='blue', size=0.8)
p <- p + geom_line(dfx_analytic2 %>% filter(gen > 1, gen < 10), mapping=aes(nblocks, prob), color='green', size=0.8)
p <- p + xlab("number of blocks") + ylab("probability") 
p <- p + simple_theme(11)
p <- p + xlim(c(0, 6))
p
# gsave(p, 'x-ancestor-blockcounts')

plot_sims(dfx, dfx_analytic, gens=2:7, mar=0.2, analytic2=dfx_analytic2,
          legnames=c('poisson-binomial', 'poisson thinning'))

plot_sims(dfx, dfx_analytic, gens=2:7, cex=12, mar=0.2,
          analytic2=dfx_analytic2,
          lgx=3.6,
          legnames=c('poisson-binomial', 'poisson thinning'),
          filename=graphdir("x-ancestor-blockcounts"))

## Take a dataframe of lengths and bin them in nbins, returning a dataframe of
## frequencies.
lengths2dataframe  <- function(d, nbins=100) {
  group_by(d, gen, grp=cut(length, seq(0, XGENLEN, length.out=nbins))) %>% 
    summarize(count=n()) %>% group_by(gen) %>% mutate(freq=count/sum(count))  %>%
      mutate(length=gsub("\\(|\\]", "", as.character(grp))) %>%
      separate(length, into=c('start', 'end'), sep=',', convert=TRUE) %>%
      mutate(midpoint=(start+end)/2) %>% select(gen, freq, grp, midpoint) %>%
      mutate(length=midpoint)
}

## X Chromsome segment length figure for paper 
#
# NOTE: Here we are counting blocks in the PARENT, not offspring. No
# generational offset.
dfxl <- read_tsv("data/x-shared-ancestor-blocklens.txt")
# note that we chop off last bin, which is biased due to too few samples
dfxl_analytic <- fit_analytic(seq(0, 1.95, by=0.01), 2:7, x_blocklen_pdf, 'mid')
p <- ggplot(dfxl %>% filter(gen > 1, gen < 8, lens < XGENLEN))
p <- p + geom_point(aes(lens, y=..density..), stat="bin", breaks=seq(0, XGENLEN, length.out=30), color='gray48')
p <- p + facet_wrap(~gen)
p <- p + geom_line(data=dfxl_analytic, mapping=aes(x=mid, y=prob), size=0.8,
                   color='blue')
p <- p + simple_theme(11)
p <- p + xlab("block length (Morgans)") + ylab("frequency") 
p <- p + xlim(c(0, 1.95))
# gsave(p, 'x-ancestor-blocklens')

dfxl_binned <-  binlens(dfxl)
# plot_sims(dfxl_binned, dfxl_analytic, gens=2:7, xcol='mid', ycol='prob')
plot_sims(dfxl_binned, dfxl_analytic, gens=2:7, cex=12, xcol='mid', ycol='prob',
          xlab="length of IBD segments (Morgans)", ylab="density",
          filename=graphdir("x-ancestor-blocklens"))



## ====== Figure 1 Plots
## Figure 1 (a): number of autosomal, X, and genetic X, genetic autosomal autosomes
## Data munging
FIG1_SIZE = 16
gen <- 1:15
nanc_auto <- function(k) 2^k
nanc_x <- function(k) fib(k+2)
# expected number of ancestors  
ngenanc_auto <- Vectorize(function(k, c=NAUTO, v=AUTOLEN)  (1 - exp(-(c + v*(k - 1))/2^(k-1)))*nanc_auto(k), 'k')
ngenanc_x_auto <- Vectorize(function(k, c=1, v=XGENLEN)  (1 - exp(-(c + v*(k - 1))/2^(k-1)))*nanc_auto(k), 'k')
ngenanc_x <- Vectorize(function(k) (1 - x_blocknum_pmf(0, k))*fib(k+2), 'k')

nanc_d <- data.frame(gen=gen, nanc_auto=nanc_auto(gen),
                     nanc_x=nanc_x(gen),
                     ngenanc_auto=ngenanc_auto(gen),
                     ngenanc_auto_alt=ngenanc_x_auto(gen),
                     ngenanc_x=ngenanc_x(gen))

nanc_dm <- melt(nanc_d, id.vars='gen')
labels <- list(`autosome ancestor` = "nanc_auto",
               `X ancestor` = "nanc_x",
               `genetic autosomes` = "ngenanc_auto",
               `genetic autosome length of X` = "ngenanc_auto_alt",
               `genetic X` = "ngenanc_x")
levels(nanc_dm$variable) <- labels

## Figure 1 (b): proportion of X ancestors
## Data munging
# three lines:
# 1. prob of genetic X ancestor | X genalogical ancestor
# 2. prob of a genetic X ancestor | autosome genealogical ancestor

# load in sim results
da <- read_tsv("data/genetic-auto-all-ancs-blockcounts.txt") %>% 
  group_by(gen) %>% mutate(prob=count[gen_anc == 1]/sum(count))
dx_a <- read_tsv("data/genetic-x-all-ancs-blockcounts.txt") %>%
  filter(gen > 1) %>%
  group_by(gen) %>% mutate(prob=count[gen_anc == 1]/sum(count))
dx_x <- read_tsv("data/genetic-x-x-ancs-blockcounts.txt") %>%
  filter(gen > 1) %>%
  group_by(gen) %>% mutate(prob=count[gen_anc == 1]/sum(count))

# probability of autosomal ancestor given genealogical ancestor 
pgenanc_auto <- Vectorize(function(k) 1-auto_blocknum_pmf(0, k))
# probability of genetic X given X genealogical ancestor
pgenanc_x <- Vectorize(function(k) (1 - x_blocknum_pmf(0, k)))
# probability of genetic X given autosomal genealogical ancestor
pgenanc_x_auto <- Vectorize(function(k) (fib(k+2)/2^(k))*(pgenanc_x(k)), 'k')
px_auto <- Vectorize(function(k) fib(k+2)/2^(k), 'k')

panc_d <- data.frame(gen=gen, 
                     pgenanc_auto=pgenanc_auto(gen),
                     pgenanc_x=pgenanc_x(gen),
                     pgenanc_x_auto=pgenanc_x_auto(gen), 
                     px_auto=px_auto(gen))
panc_dm <- melt(panc_d, id.vars='gen')
levels(panc_dm$variable) <- list(
                        `P(genetic autosome)` = "pgenanc_auto",
                        `P(genetic X | X ancestor)` = "pgenanc_x",
                        `P(genetic X | ancestor)` = "pgenanc_x_auto",
                        `P(is X ancestor)` = "px_auto")

plabels <- c(expression(P(N[auto] > 0)), 
             expression(paste("P(", N[x], " > 0 | X ancestor)")),
             expression(paste("P(", N[x], " > 0 | ancestor)")),
             expression(paste("P(X ancestor | ancestor)")))
             # expression(P(N[x] > 0 | X ancestor)), expression(P(X ancestor)))


## Figure 1 Panel Plot
setEPS()
HEIGHT <- 24
postscript(graphdir("num-ancestors"), width=phi^2*HEIGHT, height=HEIGHT)

opar <- par(no.readonly=TRUE)
CEX <- 6
lwd <- 2*CEX
layout(matrix(c(1, 2),nrow=1), widths=c(0.94, 1))
par(mar=c(2, 2, 1, 1),
    # mfrow=c(1, 2), 
    oma = c(1, 1, 0.5, 0.5), mgp = c(2, 0.4, 0), cex=CEX)
yax <- c(1, 10, 100, 1000, 10000, 100000) 
axs_col <- "grey38"
title_col <- "grey10"
subtitle_col <- "grey34"  # bold font, lighter color looks better
cols <- wes_palette("Darjeeling")
# Figure (A)
plot(value ~ gen, data=nanc_dm, log="y", type='n', axes=FALSE, 
     xlab="", ylab="", bty='n', xlim=c(1, 14), ylim=range(yax))
i <- 1
for (type in unique(nanc_dm$variable)) {
  lines(value ~ gen, data=nanc_dm[nanc_dm$variable == type & nanc_dm$gen <= 14, ], 
        col=cols[i],
        lwd=lwd)
  i <- i + 1
}
axis(1, col = title_col, col.axis = axs_col, tck=-0.018,
     at=seq(2, 14, 2), tck=0.018, lwd=CEX)
loglabel <- function(x) {
  parse(text=sprintf("10^%d", log10(x)))
}
axis(2, col = title_col, col.axis = axs_col, tck=-0.018, at=yax,
     labels=loglabel(yax),
     las=1, tck=0.018, lwd=CEX)
mtext("generation", side=1, line = 1.6, col=title_col, cex=CEX*1.2)
mtext("expected number of ancestors", side=2, line = 2, col=title_col, cex=CEX*1.2)
mtext("A", side=2, line=-1.2, las=1, at=10^(log10(10^5)*1.08), 
      cex=1.6*CEX, col=subtitle_col, font=2)
legend(1, 90000, names(labels), bty='n', fill=cols, border=0, cex=0.82,
       text.col=title_col)

# Figure (B)
cols <- wes_palette("Moonrise3")[-4]
plot(value ~ gen, data=panc_dm, type='n', axes=FALSE, 
     xlab="", ylab="", bty='n', xlim=c(0, 14), ylim=c(0, 1))
i <- 1
dfs <- list(da, dx_x, dx_a, NULL)
for (type in unique(panc_dm$variable)) {
  lines(value ~ gen, data=panc_dm[panc_dm$variable == type & panc_dm$gen <= 14, ],
        col=cols[i], lwd=lwd)
  sim <- dfs[[i]]
  if (!is.null(sim)) {
      points(prob ~ gen, data=sim, cex=0.6, pch=19, col=cols[i])
  }
  i <- i + 1
}
axis(1, col = title_col, col.axis = axs_col, tck=-0.018, 
     at=seq(2, 14, 2), tck=0.018, lwd=CEX)
axis(2, col = title_col, col.axis = axs_col, tck=-0.018, line=-1.6,
     las=1, tck=0.018, lwd=CEX)
mtext("generation", side=1, line = 1.6, col=title_col, cex=CEX*1.2)
mtext("probability", side=2, line=0.4, col=title_col, cex=CEX*1.2)
mtext("B", side=2, line=-2.8, las=1, at=1.08, cex=1.6*CEX, col=subtitle_col, font=2)
legend(8.6, 1, plabels, bty='n', fill=cols, border=0, cex=0.82, 
       title.col=title_col)
par(opar)
dev.off()

p <- ggplot(nanc_dm)
p <- p + geom_line(aes(x=gen, y=value, color=variable), size=LINESIZE)
p <- p + scale_y_log10()
p <- p + xlab("generations") + ylab("number of ancestors")
p <- p + simple_theme(FIG1_SIZE) + scale_color_discrete(guide=guide_legend(title=NULL))
p <- p + theme(legend.position = c(0.01, 1), legend.justification = c(0, 1),
               legend.background = element_rect(colour="white"))
# gsave(p, 'num-ancestors')

cols <- hue_pal()(4)
p <- ggplot(panc_dm) 
p <- p + theme(legend.position = c(0.51, 1), legend.justification = c(0, 1),
               legend.background = element_rect(colour="white"))
p <- p + geom_line(aes(x=gen, y=value, color=variable),
                                size=LINESIZE)
p <- p + geom_point(data=da, mapping=aes(x=gen, y=prob), color=cols[1])
p <- p + geom_point(data=dx_a, mapping=aes(x=gen, y=prob), color=cols[2])
p <- p + geom_point(data=dx_x, mapping=aes(x=gen, y=prob), color=cols[3])
p <- p + simple_theme(FIG1_SIZE) + scale_color_discrete(guide=guide_legend(title=NULL))
p <- p + xlab("generation") + ylab("probability") 
p <-  p + scale_x_continuous(breaks=seq(0, 14, 2))
# gsave(p, 'prob-ancestors')

## ===== Shared Ancestor is Female
# some cute math shows this converges to (1/10)* (sqrt(5) + 5)
shared_x_anc_is_female <- function(x) (fib(k+1)/fib(k+2))^2/((fib(k)/fib(k+2))^2 + (fib(k+1)/fib(k+2))^2)
plot(k <- seq(1, 15), shared_x_anc_is_female(k))
abline(a=(1/10)*(sqrt(5) + 5), b=0, col='red')


## ===== Half-Cousins
dxhc <- counts2prob(read_tsv("data/x-halfcousin-blockcounts.txt"))
dxhc_analytic <- fit_analytic(0:10, 1:10, half_cousins_x_blocknum_pmf)
dxhc_analytic2 <- fit_analytic(0:10, 1:10, half_cousins_x_blocknum_thinned_pmf)

p <- ggplot(dxhc %>% filter(gen < 8)) + geom_point(aes(x=nblocks, y=prob), color='gray48')
p <- p + facet_wrap(~gen)
p <- p + geom_line(data=dxhc_analytic %>% filter(gen < 8), aes(x=nblocks, y=prob), color='cornflowerblue', size=0.8)
p <- p + geom_line(data=dxhc_analytic2 %>% filter(gen < 8), aes(x=nblocks, y=prob), color='green', size=0.8)
p <- p + simple_theme(FIG1_SIZE) 
p <- p + xlab("number of blocks") + ylab("probability") 
p <- p + simple_theme(11)
p <- p + scale_x_continuous(breaks=seq(0, 10, 2))
# p
# gsave(p, 'x-halfcousins-blockcounts')

plot_sims(dxhc, dxhc_analytic, gens=2:7)

plot_sims(dxhc, dxhc_analytic, gens=2:7, cex=12, 
          analytic2=dxhc_analytic2,
          filename=graphdir("x-halfcousins-blockcounts"))


## ===== Half-Cousins Lengths

dfxhcl <- read_tsv("data/x-halfcousin-lens.txt")
dfxhcl_analytic <- fit_analytic(seq(0, 1.95, by=0.01), 2:7,
                                x_halfcousins_blocklen_pdf, 'mid')
# plot_sims(binlens(dfxhcl), dfxhcl_analytic, gen=2:7, xcol='mid', ycol='prob')

# plot_sims(binlens(dfxhcl), dfxhcl_analytic, sims2=dfxl_binned,
#           analytic2=dfxl_analytic,
#           gen=2:7, xcol='mid', ycol='prob', 
#           xlab="length of IBD segments", ylab="density")

plot_sims(binlens(dfxhcl, breaks=12), dfxhcl_analytic, 
          gen=2:7, xcol='mid', ycol='prob', xrng=c(0, 2),  # smin=TRUE,
          xlab="length of IBD segments (Morgans)", ylab="density")

plot_sims(binlens(dfxhcl, breaks=12), dfxhcl_analytic, 
          gen=2:7, xcol='mid', ycol='prob', cex=12, xrng=c(0, 2),
          smin=0.1,
          xlab="length of IBD segments (Morgans)", ylab="density",
          filename=graphdir("x-halfcousins-blocklens"))

## Full Cousins
dxfc <- counts2prob(read_tsv("data/x-fullcousin-blockcounts.txt"))
dxfc_analytic <- fit_analytic(0:10, 1:10, full_cousins_x_blocknum_pmf)
p <- ggplot(dxfc %>% filter(gen < 8)) + geom_point(aes(x=nblocks, y=prob), color='gray48')
p <- p + facet_wrap(~gen)
p <- p + geom_line(data=dxfc_analytic %>% filter(gen < 8), aes(x=nblocks, y=prob), color='blue', size=0.8)
p <- p + simple_theme(10) 
p <- p + xlab("number of blocks") + ylab("probability") 
p <- p + scale_x_continuous(breaks=seq(0, 10, 2))
p
# gsave(p, 'x-fullcousins-blockcounts')

plot_sims(dxfc, dxfc_analytic, gens=2:7)

plot_sims(dxfc, dxfc_analytic, gens=2:7, cex=12, 
          filename="paper/images/x-fullcousins-blockcounts.eps")


## ====== Full and Half Cousin distributions in one plot

# poisson-binomial
plot_sims(dxhc, dxhc_analytic, sims2=dxfc, analytic2=dxfc_analytic, gens=2:7)

# figure for paper
plot_sims(dxhc, dxhc_analytic, sims2=dxfc, analytic2=dxfc_analytic, gens=2:7,
          pointcols=c('#AAC3F0', '#CFC1EE'), cex=12,
          legnames=c('half-cousin', 'full-cousin'),
          lgx=4, lgy=0.78,
          filename=graphdir("x-full-half-blockcounts"))

## ======  Inference
dhsi <- read_tsv("data/x-halfcousin-blockcounts-rms.txt")

ggplot(dhsi) + geom_bar(aes(x=rm, y=count), stat='identity') + facet_wrap(~gen)

ggplot(dhsi) + geom_bar(aes(x=nblocks, y=count), stat='identity') + facet_wrap(~gen)

ggplot(dhsi) + geom_line(aes(x=rm, y=count, color=as.factor(gen), group=as.factor(gen))) + facet_wrap(~nblocks, scales="free")

pp <- dhsi %>% group_by(gen) %>% mutate(prob=count/sum(count))

p <- ggplot(pp %>% filter(nblocks > 0))
p <- p + geom_line(aes(x=rm, y=prob, color=as.factor(nblocks), group=nblocks))
p <- p + facet_wrap(~gen, scales="free")
p

# create df of likelihood and posterior using some sweet purrr goodness.
dd <- data.frame(k=2:12) %>% split(.$k) %>% 
  map(~ data.frame(k=.$k, r=rms_sharedanc(.$k))) %>%
  map(function(.) {
    d <- expand.grid(k=.$k, r=.$r, nblocks=0:20)
    d$ll <- with(d, half_cousins_blocknum_rm_ll(nblocks, r, k))
    d$post <- with(d, half_cousins_blocknum_rm_posterior(nblocks, r, k))
    d
  }) %>%
  rbind_all() %>% group_by(k, nblocks) %>% mutate(ll=ll/max(ll))

# posterior plot for figure 7. This does *not* condition on N > 0
HEIGHT <- 48
setEPS()
postscript(graphdir("rm-posterior"), width=phi*HEIGHT, height=HEIGHT)
ll_plot(dd, 2:7, ncol=3, cex=6, coln='post', title='number of IBD\nsegments')
dev.off()

## ====== Validation on RMs
# if we have two types of balls and arrange them randomly (equal probability)
# what are the odds we arrange them in a way such that no two white balls are
# adjacent? This is the proportion of X ancestors over all ancestors,
# F(k+2)/2^k. Simulation:
simLineage <- function(k) paste0(sample(c('m', 'f'), k, replace=TRUE),collapse="")
isX <- function(x) regexpr('mm', x) == -1
k <- 1:20
N <- 1000
sims <- do.call(rbind, lapply(k, function(k) {
                data.frame(k=k, p=sum(replicate(N, isX(simLineage(k))))/N)
}))
p <- ggplot(sims) + geom_line(aes(k, p))
p <- p + geom_point(data=data.frame(k=k, prob=fib(k+2)/2^k), aes(k, prob))
p  # hey, cute, it works.


## ====== Prob full cousins given female is X

k <- 1:10
prob_full_cousins <- function(k) ((fib(k)^2)/(fib(k)^2 + fib(k+1)^2))^2
plot(k, prob_full_cousins(k), type='l')
plot(k, prob_full_cousins(k)/(0.5*(fib(k+1)^2/2^(k-1)) + 0.5*(fib(k)/2^(k-1))),
     type='l')
abline(a=(1/10)*(3 - sqrt(5)), b=0, col='cornflowerblue')


## ====== Pascal's Triangle Figure

## ====== Prob full cousins given female is X

plotbifib(9)
plotbifib(9, cex=4, filename=graphdir('bifib'))
dev.off()

## ====== Alternate MC Process

ibdmc <- function(k, nu=1.96) {
  u <- rexp(100, nu*k*2^-k)
  up <- rexp(100, nu*k)
  flip <- rbinom(1, 1, 2^-k)
  is_ibd <- seq_along(u) %% 2 == 0
  if (flip) {
    pp <- ifelse(!is_ibd, u, up) 
  } else {
    pp <- ifelse(is_ibd, up, u) 
  }
  i <- cumsum(pp) <= 1.96
  if (!any(i)) return(list(segs=nu, state=flip))
  segs <- pp[1:max(which(i))]
  state <- is_ibd[1:(1+max(which(i)))]
  tot <- sum(segs)
  segs <- c(segs, nu - tot)
  return(list(segs=segs, states=state))
}

ibdMC <- function(n, k, nu=1.96) {
  replicate(n, {ibdmc(k, nu)}, simplify=FALSE)
}

nsegs <- function(res) {
  data.frame(n=sapply(res, function(x) sum(x$states)))
}

nseg_table <- function(nsegs) {
  pt <- prop.table(table(nsegs))
  setNames(as.data.frame(pt), c('nsegs', 'prob'))
}

ddd <- do.call(rbind, lapply(2:9, function(k) {
  d <- nseg_table(nsegs(ibdMC(1000, k, AUTOLEN)))
  d$gen <- k
  d$nsegs <- as.numeric(as.character(d$nsegs))
  d
}))

p <- ggplot(dfa) + geom_point(aes(x=nblocks,y=prob,color=type)) + facet_wrap(~gen)
dfa_analytic <- fit_analytic(0:40, 1:10, auto_blocknum_pmf) 
p <- p + geom_line(dfa_analytic, mapping=aes(nblocks, prob), color='blue', size=0.8)
p <- p + facet_wrap(~gen, scales="free_y", ncol=2)
p + geom_point(data=ddd, aes(nsegs, prob), color='purple')


## ====== Pedigree Collapse
# both of these model pedigree collapse in a single individual's pedigree
pedcollapse <- Vectorize(function(k, n=1e5) prod((1-((1:2^k)/n))), 'k')
xpedcollapse <- Vectorize(function(k, n=1e5) prod((1-(2*(1:(fib(k)))/n)))*prod((1-(2*(1:(fib(k+1)))/n))), 'k')
# xpedcollapse2 <- Vectorize(function(k, n=1e5) prod((1-((1:(fib(k)))/n))), 'k') # wrong; for comparison only, not accounting for sex
plot(k <- 1:20, pedcollapse(k), type='l')
lines(k, xpedcollapse(k), col='red')
# lines(k, xpedcollapse2(k), col='green')


axs_col <- "grey38"
title_col <- "grey10"
cex <- 14
axis_cexmul <- 1.6
HEIGHT <- 80
setEPS()
postscript(graphdir("ped-collapse"), width=phi*HEIGHT, height=HEIGHT)
opar <- par(no.readonly=TRUE)
cols <- wes_palette("Darjeeling")
par(mar=c(5, 5, 0, 0), cex=cex)
plot.new()
maxk <- 20 
k <- 1:maxk
plot.window(xlim=c(0, maxk), ylim=c(0, 1))
lines(k, pedcollapse(k), col=cols[3], lwd=2.1*cex)
lines(k, xpedcollapse(k), col=cols[5], lwd=2.1*cex)
axis(1, col = axs_col, col.axis = axs_col, cex=2.1*cex, tck=0.018, lwd=cex,
     at=seq(0, maxk, 5), labels=seq(0, maxk, 5))
axis(2, col = axs_col, col.axis = axs_col, cex=2.1*cex, tck=0.018, lwd=cex,
     at=seq(0, 1, 0.2), las=1)
mtext("generation", side=1, line=3, col=title_col, cex=cex*axis_cexmul)
mtext("probability all ancestors distinct", side=2, line=3,
      col=title_col, cex=cex*axis_cexmul)
legend(14, 1, c("genealogical ancestors", "X ancestors"), bty='n',
       cex=1.2,
       fill=cols[c(3, 5)], border=0, text.col=title_col)
par(opar)
dev.off()

## ==== Prob(X anc | N = 0)
setEPS()
postscript(graphdir("prob-xanc-n0"), width=phi*HEIGHT, height=HEIGHT)
opar <- par(no.readonly=TRUE)
cols <- wes_palette("Darjeeling")
maxk <- 15
maxp <- 1
cex <- 14
k <- 1:maxk
par(mar=c(5, 5, 0, 0), cex=cex)
plot.new()
plot.window(xlim=c(0, maxk), ylim=c(0, maxp))
lines(k, prob_xanc_N0(k), col=cols[3], lwd=2.1*cex, type='b', pch=20)
lines(k, prob_x_anc(k), col=cols[3], lwd=2.1*cex, lty=2, type='b', pch=20)
lines(k, prob_x_sharedanc_N0(k), col=cols[5], lwd=2.1*cex, type='b', pch=20)
lines(k, prob_x_shared_anc(k), col=cols[5], lwd=2.1*cex, lty=2, type='b', pch=20)
axis(1, col = axs_col, col.axis = axs_col, cex=2.1*cex, tck=0.018, lwd=cex,
     at=c(1, seq(5, maxk, 5)), labels=c(1, seq(5, maxk, 5)))
axis(2, col = axs_col, col.axis = axs_col, cex=2.1*cex, tck=0.018, lwd=cex,
     at=seq(0, maxp, 0.2), las=1)
mtext("generation", side=1, line=3, col=title_col, cex=cex*axis_cexmul)
mtext(expression(paste("P(X ancestor | ", N[x], " = 0)")),
      side=2, line=3, col=title_col, cex=cex*axis_cexmul)
legend(11, 1, c("X ancestor", "X ancestor prior",
                "half-cousins", "half-cousins prior"),
       bty='n', lty=c(1, 2, 1, 2), lwd=axis_cexmul*cex,
       col=cols[c(3, 3, 5, 5)], border=0, text.col=title_col)
par(opar)
dev.off()


# ==== Thinning total variation distance
maxk <- 10
maxp <- 1
k <- 1:maxk
k_auto <- 1:maxk
x <- 0:100
tvd <- sapply(k, function(k)
       0.5*sum(abs(x_blocknum_pmf(x, k) - x_blocknum_thinned_pmf(x, k))))
tvd_auto <- sapply(k_auto, function(k) {
       0.5*sum(abs(auto_blocknum_pmf(x, k) - 
                   auto_blocknum_thinned_pmf(x, k)))
       })

setEPS()
postscript(graphdir("total-var-dist"), width=phi*HEIGHT, height=HEIGHT)
opar <- par(no.readonly=TRUE)
cols <- wes_palette("Darjeeling")
cex <- 14
par(mar=c(5, 5, 0, 0), cex=cex)
plot.new()
plot.window(xlim=c(0, maxk), ylim=c(0, maxp))
lines(k, tvd, col=cols[3], lwd=2.1*cex, type='b', pch=20)
lines(k_auto, tvd_auto, col=cols[5], lwd=2.1*cex, type='b', pch=20)
axis(1, col = axs_col, col.axis = axs_col, cex=2.1*cex, tck=0.018, lwd=cex,
     at=c(1, seq(5, maxk, 5)), labels=c(1, seq(5, maxk, 5)))
axis(2, col = axs_col, col.axis = axs_col, cex=2.1*cex, tck=0.018, lwd=cex,
     at=seq(0, maxp, 0.2), las=1)
mtext("generation", side=1, line=3, col=title_col, cex=cex*axis_cexmul)
mtext("total variation distance",
      side=2, line=3, col=title_col, cex=cex*axis_cexmul)
legend(7.8, 0.94, c("X chromosome", "autosomes"), 
       bty='n', fill=cols[c(3, 5)], border=0, text.col=title_col,
       cex=1.2)
par(opar)
dev.off()
