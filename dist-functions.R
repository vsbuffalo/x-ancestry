# X ancestry distributions
library(colorspace)

## Constants
XGENLEN = 1.96
MAXBPOINTS = 1e3
NAUTO = 22
AUTOLEN = 35.24

phi <- (1 + sqrt(5))/2
psi <- (1 - sqrt(5))/2
fib <- function(k) (phi^{k} - psi^{k})/sqrt(5)

is_pmf <- function(probs, EPS=1e-6) {
  stopifnot(all(probs >= 0))
  stopifnot(all(probs <= 1))
  stopifnot(abs(sum(probs) - 1) < EPS)
}


tuple2midpoint <- function(x) {
  # extract midpoint from stringified tuple (from Pandas).
  sapply(strsplit(gsub("\\(([^,]+), ([^\\)]+)\\)", "\\1;;\\2", x), ";;"), 
         function(y) mean(as.numeric(y)))
}

counts2prob <- function(d, col='nblocks') {
  col_i <- which(colnames(d) == col)
  colnames(d)[col_i] <- "YYY"
  out <- d %>% group_by(gen) %>% mutate(prob=count/sum(count))
  col_i <- which(colnames(out) == "YYY")
  colnames(out)[col_i] <- col
  out
}

fit_analytic <- function(x, gens, fun, xname='nblocks') {
  out <- expand.grid(x, gen=gens) %>% group_by(gen) %>% mutate(prob=fun(x, gen))
  colnames(out)[1] <- xname
  return(out)
}

binlens <- function(d, breaks='sturges') {
  map_df(d %>% split(.$gen), function(x) {
    h <- hist(x$lens, breaks)
    data.frame(gen=x$gen[1], mid=h$mids, prob=h$density)
  })
}



hcols <- function (n, h = c(-98, 147), c. = c(100, 87), l = c(66, 70),
            power = c(0.428571428571429, 2.72448979591837), fixup = TRUE,
            gamma = NULL, alpha = 1, ...) {
    if (!is.null(gamma))
      warning("'gamma' is deprecated and has no effect")
    if (n < 1L)
      return(character(0L))
    h <- rep(h, length.out = 2L)
    c <- rep(c., length.out = 2L)
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, 0, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L],
                         C = c[2L] - diff(c) * rval^power[1L], H = h[2L] - diff(h) *
                           rval), fixup = fixup, ...)
                         if (!missing(alpha)) {
                           alpha <- pmax(pmin(alpha, 1), 0)
                           alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)),
                                           width = 2L, upper.case = TRUE)
                           rval <- paste(rval, alpha, sep = "")
                         }
    return(rval)
}

simulateAlpha <- Vectorize(function(col, alpha) {
  cols <- rev(colorRampPalette(c(col, '#FFFFFF'))(100))
  cols[floor(alpha*100)]
}, 'col')


addAlpha <- Vectorize(function(col, alpha) {
  args <- as.list(col2rgb(col)[,1]/255)
  args$alpha <- alpha
  do.call(rgb, args)
}, c('col', 'alpha'))


plot_sims <- function(sims, analytic, sims2=NULL, analytic2=NULL,
                      gens=NULL, cex=1, margin=0.62,
                      # linecols=c('cornflowerblue', 'mediumpurple2'),
                      # linecols=c(rgb(12, 161, 137, maxColorValue=255), rgb(104, 188, 215, maxColorValue=255)),
                      # linecols=wes_palette("Darjeeling")[c(5, 1)],
                      linecols=wes_palette("Darjeeling")[c(5, 3)],
                      pointcols=c('gray66'),
                      pointalpha=0.5,
                      filename=NULL, xcol="nblocks", ycol="prob",
                      xrng=NULL, smin = FALSE,  # sims bound analytic
                      lg_gen=4, legnames=NULL, lgx=2.4, lgy=0.68,
                      ylab="probability", xlab="number of IBD segments") {
  opar <- par(no.readonly = TRUE)
  HEIGHT <- 80
  if (!is.null(filename)) {
    setEPS()
    postscript(filename, width=phi*HEIGHT, height=HEIGHT)
  }
  axs_col <- "grey38"
  title_col <- "grey10"
  sims <- as.data.frame(sims)
  analytic <- as.data.frame(analytic)
  if (!is.null(sims2))
    sims2 <- as.data.frame(sims2)
  if (!is.null(analytic2))
    analytic2 <- as.data.frame(analytic2)
  if (!is.null(sims2)) {
    # adjust point color
    pointcols <- simulateAlpha(linecols, pointalpha)
  }
  if (is.null(gens))
    gens <- unique(sims$gen)
  ncols <- floor(length(gens)/2) + (length(gens) %% 2 > 0)
  par(mfrow=c(2, ncols), mar=rep(margin, 4),
      oma = c(4, 4, 0.5, 0.5), mgp = c(2, 0.4, 0), cex=cex)
  yrng <- range(sims[, ycol])
  if (is.null(xrng)) {
    xrng <- range(sims[, xcol])
  } 
  i <- 0
  for (k in gens) {
    d <- sims[sims$gen == k, ]
    da <- analytic[analytic$gen == k, ] 
    if (smin) {
      smax <- max(d[, xcol]) + smin
      da <- da[da[, xcol] <= smax, ]
    }
    plot(d[, xcol], d[, ycol], ylim=yrng, xlim=xrng, axes=FALSE, pch=19, cex=0.8,
         col=pointcols[1])
    if (!is.null(sims2)) {
      d2 <- sims2[sims2$gen == k, ]
      points(d2[, xcol], d2[, ycol], ylim=yrng, xlim=xrng, pch=19, cex=0.8,
             col=pointcols[2])
    }
    lines(da[da[, xcol] <= xrng[2], xcol], 
          da[da[, xcol] <= xrng[2], ycol],
          ylim=yrng, xlim=xrng, type='l', col=linecols[1], lwd=cex*2)
    if (!is.null(analytic2)) {
      da2 <- analytic2[analytic2$gen == k, ] 
      lines(da2[da2[, xcol] <= xrng[2], xcol], da2[da2[, xcol] <= xrng[2], ycol],
            ylim=yrng, xlim=xrng, type='l', col=linecols[2], lwd=cex*2)
    }
    if (!is.null(legnames) && k == lg_gen) {
      legend(lgx, lgy, legnames, bty='n', fill=linecols,
             border=0, cex=0.86,
             text.col=title_col)
    }


    if (i >= 2*ncols - ncols) {
      # x-axis
      axis(1, col = axs_col, col.axis = axs_col, cex=cex, tck=0.018, lwd=cex)
    } else {
      axis(1, col = axs_col, col.axis = axs_col, labels=FALSE, tick=TRUE, 
           tck = .01, cex=cex, lwd=cex)
    }
    if (i %% ncols == 0) {
      # y-axis
      axis(2, col = axs_col, col.axis = axs_col, cex=cex, lwd=cex, tck=0.018,
           las=1)
      # if (i == 2*ncols - ncols)
        # box(bty="L")
    }
    mtext(sprintf("k = %d", k), side=3, line=-2.4, col=axs_col, cex=cex*1.1)
    i <- i + 1
  }
  mtext(xlab, side=1, outer=TRUE, line = 2.1, col=title_col, cex=cex*1.4)
  mtext(ylab, side=2, outer=TRUE, line = 2.1, col=title_col, cex=cex*1.4)
  if (!is.null(filename)) 
    dev.off()
  on.exit(par(opar))
}

ll_plot <- function(x, gens, ncol=4, cex=1, margin=1.2,
                    coln='post', lg_gen=4,
                    cols=NULL,
                    title=NULL,
                    xlab="number of recombinational meioses",
                    ylab="posterior") {
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  axs_col <- "grey38"
  title_col <- "grey10"
  nrow <- floor(length(gens)/ncol) + (length(gens) %% ncol > 0)
  par(mfrow=c(nrow, ncol), oma=c(4, 4, 0.5, 0.5), mar=rep(margin, 4),
      mgp = c(2, 0.4, 0), cex=cex)
  nblocks <- 0:8
  if (is.null(cols)) {
    cols <- brewer.pal(length(nblocks)+2L, "Spectral")[-c(1, 6)]
    # cols <- hcols(length(nblocks)+2L)[-4][-4]
  } else if(is.function(cols)) {
    cols <- cols(length(nblocks)) # for palette
  }
  for (gen in gens) {
    d <- x %>% filter(k == gen)
    plot.new()
    rms <- rms_sharedanc(gen)
    plot.window(xlim=range(rms), ylim=c(0, 1))
    axis(1, col = axs_col, col.axis = axs_col, cex=cex, tck=0.018, lwd=cex, at=rms)
    axis(2, col = axs_col, col.axis = axs_col, cex=cex, tck=0.018, lwd=cex, las=1)
    sapply(seq_along(unique(d$nblocks)), function(i) {
      n <- unique(d$nblocks)[i]
      dn <- d %>% filter(nblocks == n) %>% as.data.frame
      lines(dn$r, dn[, coln], col=cols[i], lwd=1.4*cex, type='b', pch=19, cex=0.1*cex)
    })
  if (gen == lg_gen)
    legend(rms[length(rms)-1]-0.1, 0.94, nblocks, bty='n', fill=cols, border=0, cex=0.86,
           text.col=title_col, title=title)

    # RMS line
    lines(rms, recomb_meioses_pmf(2*gen-1, rms), lwd=1.4*cex, pch=19, cex=0.1*cex, 
          lty=2, type='b', col='gray28')
    mtext(sprintf("k = %d", gen), side=3, line=-1.8, col=axs_col, cex=cex*1.1)
  }
  mtext(xlab, side=1, outer=TRUE, line = 2.1, col=title_col, cex=cex*1.4)
  mtext(ylab, side=2, outer=TRUE, line = 2.1, col=title_col, cex=cex*1.4)
}



plotbifib <- function(h, cex=1, filename=NULL) {
  # warning, all who dare to enter; this code is a land of weird constants
  # to make the figure look right.
  CEX <- cex
  if (!is.null(filename)) {
    setEPS()
    postscript(filename, height=70, width=120)
  }
  opar <- par(no.readonly=TRUE)
  plot.new()
  par(mar=c(1, 1, 1, 1), cex=CEX) 
  bifib <- function(h, buf=0.5, xbuf=0, ybuf=0, ymul=0.5, xmul=1, font=1,
                    cex=CEX, col='black') {
    lx <- numeric(h)
    ly <- numeric(h)
    for(k in 0:(h-2)) {
      r <- 0:k
      rms <- choose(r+1, k-r)
      x <- xbuf - (length(rms)+1.5)/2 + seq_along(rms) - 0.25
      y <- ymul*k + ybuf
      lx[k+1] <- x[1]
      ly[k+1] <- y[1] - 0.03
      text(x=xmul*x, y=y, label=rms, cex=cex, col=col, font=font)
    }
    return(list(lx, ly))
  }
  ymul <- 0.3
  plot.window(ylim=c(h/2 - 2.3, 0), xlim=c(-(h-1)/2 - 0.5, (h-1)/2))
  # returns positions of left most text
  pos <- bifib(h, buf=0, ymul=ymul, font=2, col='gray10')
  red <- rgb(222, 44, 38, maxColorValue=255)
  blue <- rgb(67, 162, 202, maxColorValue=255)
  # red numbers are numbers shift down and right
  bifib(h-1, xbuf=0.3, ybuf=ymul - ymul/3, buf=0, ymul=ymul, col=red,
        cex=0.8*CEX)
  # red numbers are numbers shift down twice and right
  bifib(h-2, xbuf=0.18, ybuf=2*ymul - ymul/3, buf=0, ymul=ymul, col=blue, 
        cex=0.8*CEX)
  # below adds a one to top left
  bifib(2, xbuf=0.18, ybuf=ymul - ymul/3, buf=0, ymul=ymul, xmul=2.5, col=blue, 
        cex=0.8*CEX)
  for (k in 2:(h-1)) text((h-1)/2, ymul*(k-1), sprintf("k = %d", k-1), col='gray38', cex=0.8*CEX)
  # points(pos[[1]]-0.2, pos[[2]]-0.5)
  offset <- -0.22; xoffset <- -1.5
  pos <- list(pos[[1]][-h], pos[[2]][-h])
  pos <- list(pos[[1]][-1], pos[[2]][-1])
  segments(-(pos[[1]]+offset - xoffset),
           pos[[2]] + offset,
          -(seq((h-1)/2, -(h-1)/2)[2:(h-1)] - 0.76 - xoffset),
          h/2 - 2.36,
           col='gray52', lwd=c(1, 10)[!is.null(filename) + 1L])
  text((h-1)/2 - 1:(h-1) - 0.12, h/2-2.25, sprintf("r = %d", k:1 - 1), col='gray38',
       cex=0.8*CEX)
  if (!is.null(filename)) {
    dev.off()
  }
  par(opar)
}



#' Return the possible females over a genealogy back k generations, from a present-day female.
#'
#'
#' @param k Number of generations.
rms <- function(k) seq(floor(k/2), k)


#' Return the possible females between two female cousins over a genealogy back k generations.
#'
#'
#' @param k Number of generations.
rms_sharedanc <- function(k, fullcousins=FALSE, shared_ancestor=NULL,
                          two_daughters=FALSE) {
  if (fullcousins)
    return(2 + seq(2*floor((k-2)/2), 2*k-2))
  if (is.null(shared_ancestor))
    return(seq(floor((2*k-1)/2), 2*k-1))
  if (shared_ancestor == 'f' && !two_daughters)
    return(1 + seq(2*floor((k-1)/2), 2*(k-1)))
  if (shared_ancestor == 'm' || (shared_ancestor == 'f' && two_daughters)) {
    offset <- 2 + (shared_ancestor == 'f')
    return(offset + seq(2*floor((k-2)/2), 2*(k-2)))
  }
  stop('unimplemented parameter choice')
}


#' PMF for number of recombinational meioses (females)
#' 
#' This PMF returns the probability of r recombinational meioses occuring in k
#' generations (e.g. to an ancestor of unknown sex in the kth generation). The
#' present-day individual is a female, which does not constrain the sex of the
#' next individual further back in the genealogy.
#'
#' @param k Number of generations.
#' @param r Number of recombinational meioses (females).
recomb_meioses_pmf <- function(k, r) choose(r + 1, k - r)/fib(k+2)


#' PMF for number of IBD segments shared between an ancestor in the kth generation
#' and present-day female.
#'
#' This PMF models the number of blocks between an ancestor in the kth
#' generation and a present-day female (which does not constrain the sex of the
#' next individual futher back in the genealogy). This models recombination as
#' a Poisson process and treats each block as surviving segregation
#' independently.
#'
#' @param nblocks Number of blocks
#' @param k Number of generations to ancestor.
#' @param genlen Genetic length (in Morgans) of the sex chromosome.
x_blocknum_pmf <- Vectorize(function(nblocks, k, genlen=XGENLEN, nchrom=1L) {
  # TODO change name to xblocknum
  if (nblocks < 0)
    return(NA)
  recombs <- seq(floor((k-1)/2), k-1)
  breakpoints <- seq(0, MAXBPOINTS)
  prob <- sum(sapply(recombs, function(r) dpois(breakpoints, r*genlen)*
              dbinom(nblocks, breakpoints+nchrom, prob=1/2^r)*
              recomb_meioses_pmf(k-1, r)))
  return(sum(prob))
}, c('nblocks', 'k'))

x_blocknum_thinned_pmf<- Vectorize(function(nblocks, k, genlen=XGENLEN, nchrom=1L) {
  # TODO change name to xblocknum
  if (nblocks < 0)
    return(NA)
  recombs <- seq(floor((k-1)/2), k-1)
  prob <- sum(sapply(recombs, function(r) dpois(nblocks, (nchrom + r*genlen)/2^r)*recomb_meioses_pmf(k-1, r)))
  return(sum(prob))
}, c('nblocks', 'k'))




#' PDF for IDB segmenth length between an ancestor in the kth generation
#' and a present-day female.
#'
#'
#' @param u Segment length (in Morgans)
#' @param k Number of generations to ancestor
x_blocklen_pdf <- Vectorize(function(u, k) {
  if (u < 0) return(NA)
  recombs <- seq(floor((k)/2), k)
  return(sum(dexp(u, recombs)*recomb_meioses_pmf(k, recombs)))
}, c('u', 'k'))

#' PDF for IDB segmenth length between two present-day female half-cousins
#'
#'
#' @param u Segment length (in Morgans)
#' @param k Number of generations to ancestor
x_halfcousins_blocklen_pdf <- Vectorize(function(u, k) {
  if (u < 0) return(NA)
  recombs <- rms_sharedanc(k)
  return(sum(dexp(u, recombs)*recomb_meioses_pmf(2*k-1, recombs)))
}, c('u', 'k'))




auto_blocknum_pmf <- Vectorize(function(nblocks, k) {
  breakpoints <- seq(0, MAXBPOINTS)
  sum(dbinom(nblocks, breakpoints+NAUTO, prob=1/2^(k-1))*
      dpois(breakpoints, AUTOLEN*(k-1)))
}, c('nblocks', 'k'))

auto_blocknum_thinned_pmf <- Vectorize(function(nblocks, k) {
  dpois(nblocks, (22 + AUTOLEN*(k-1))/2^{k-1})
}, c('nblocks', 'k'))


auto_blocklen_pdf <- Vectorize(function(lens, k) {
   k*exp(-k*lens)
}, c('lens', 'k'))


#' The number of recombinational meioses from a female in the kth generation to
#' a present-day female
#'
#' @param k Number of generations.
#' @param r Number of recombinational meioses.
lineage_recomb_meioses_pmf <- function(k, r) choose(r, k - r)/fib(k+1)

#' PMF for number of blocks shared between two half-cousins with
#' ancestor k generations ago, with known sex
#'
#' This PMF marginalizes over the number of recombinational meioses down each
#' lineage from a female or male shared ancestor separetely (these are r1 and
#' r2). Both present-day cousins are assumed to be female.
#'
#' @param nblocks Number of blocks shared between two half-cousins k
#' generations ago.
#' @param k Number of generations until shared ancestor.
#' @param genlen The genetic length (in Morgans) of the sex chromosome.
#' @param shared_ancestor either 'F' or 'M' for female or male.
half_cousins_blocknum_r1r2_pmf <- function(nblocks, k, genlen=XGENLEN, 
                                           shared_ancestor="F") {
  args = c('nblocks', 'k')
  breakpoints <- seq(0, MAXBPOINTS)
  fun <- list(f = Vectorize(function(nblocks, k, genlen=XGENLEN) {
                if (nblocks < 0) return(NA)
                recombs <- (floor((gen-1)/2) + 1):gen
                x <- 0
                for (r1 in recombs) {
                  for (r2 in recombs) {
                    x <- x + sum(dpois(breakpoints, (r1 + r2)*genlen)*
                        dbinom(nblocks, breakpoints+1L, p=1/2^(r1 + r2 - 1))*
                        lineage_recomb_meioses_pmf(gen, r1)*
                        lineage_recomb_meioses_pmf(gen, r2))
                  }
                } 
                return(x)
              }, args),
              m = function() stop("not implemented"))
  fun[[tolower(shared_ancestor)]](nblocks, k, genlen)
}

# TODO: these need to be checked
prob_x_female_given_r_v1 <- function(r, k) {
  s <- seq(floor((k-1)/2), k-1)
  b <- sum(choose(s+1, k - s - 1)*choose(r - s, k - r + s))
  (fib(2*k-1)/(choose(r+1, 2*k-r-1)*(fib(k)^2 + fib(k+1)^2))) * b
}

prob_r_female_convolution <- Vectorize(function(r, k) {
  # the number of ways to distribute r females over a genealogy
  # with an ancestor in the kth generation.
  if (!(r %in% rms_sharedanc(k, shared_ancestor='f')))
    return(0)
  s <- seq(floor((k-1)/2), k-1)
  sum(choose(s+1, k-s-1)*choose(r-s, k-r+s))
}, 'r')

prob_r_female_convolution_2daughters <- Vectorize(
function(r, k) {
  if (!(r %in% rms_sharedanc(k, shared_ancestor='f', two_daughters=TRUE)))
    return()
  s <- seq(floor((k-2)/2), k-2)
  sum(choose(s+1, k-s-2)*choose(r-s-2, k-r+s+1))
}, 'r')

prob_r_male_convolution <- Vectorize(function(r, k) {
  if (!(r %in% rms_sharedanc(k, shared_ancestor='m')))
    return(0)
  s <- seq(floor((k-2)/2), k-2)
  sum(choose(s+1, k-s-2)*choose(r-s-1, k-r+s))
}, 'r')

prob_r_given_female <- Vectorize(function(r, k) {
  prob_r_female_convolution(r, k)/fib(k+1)^2
}, 'r')

prob_r_given_male <- Vectorize(function(r, k) {
  prob_r_male_convolution(r, k)/fib(k)^2
}, 'r')

prob_r_given_female_2daughters <- function(r, k) {
  prob_r_female_convolution_2daughters(r, k)/fib(k)^2
}

prob_x_female_given_r_v2 <- function(r, k) {
  pf <- fib(k+1)^2/(fib(k)^2 + fib(k+1)^2)
  return((pf*prob_r_given_female(r, k))/recomb_meioses_pmf(2*k-1, r))
}

prob_x_female_given_r3 <- function(r, k) {
  (1/recomb_meioses_pmf(2*k-1, r))*prob_r_female_convolution(r, k)/(fib(k+1)^2 + fib(k)^2)
}

#' PMF for IBD segments shared between two half-cousins with ancestor k
#' generations ago
#'
#' This PMF marginalizes over the number of recombinational meioses between two
#' present-day females, using a convolution approach to finding the probability
#' of the ancestor's sex given the number of recombinational meioses.
#'
#' @param nblocks Number of blocks shared between two half-cousins k
#' generations ago.
#' @param k Number of generations until shared ancestor.
#' @param genlen The genetic length (in Morgans) of the X chromosome.
#'
half_cousins_x_blocknum_pmf <- Vectorize(function(nblocks, k, genlen=XGENLEN) {
  breakpoints <- seq(0, MAXBPOINTS)
  recombs <- seq(floor((2*k-1)/2), 2*k-1)
  sum(sapply(recombs, function(r) {
           # prob_x_female_given_r_v2 and prob_x_female_given_r3 are identical
           (prob_x_female_given_r_v2(r, k)*sum(dpois(breakpoints, (r+1)*genlen)*
                                           dbinom(nblocks, breakpoints+1, 1/2^r)) + 
         (1-prob_x_female_given_r_v2(r, k))*sum(dpois(breakpoints, r*genlen)*
                                             dbinom(nblocks, breakpoints+1, 1/2^r)))*
             recomb_meioses_pmf(2*k-1, r)
  }))
}, c('nblocks', 'k'))


#' Thinned PMF for IBD segments shared between two half-cousins with ancestor k
#' generations ago
#'
#' This PMF marginalizes over the number of recombinational meioses between two
#' present-day females, using a convolution approach to finding the probability
#' of the ancestor's sex given the number of recombinational meioses.
#'
#' @param nblocks Number of blocks shared between two half-cousins k
#' generations ago.
#' @param k Number of generations until shared ancestor.
#' @param genlen The genetic length (in Morgans) of the X chromosome.
#'
half_cousins_x_blocknum_thinned_pmf <- Vectorize(function(nblocks, k, genlen=XGENLEN) {
  recombs <- seq(floor((2*k-1)/2), 2*k-1)
  sum(sapply(recombs, function(r) {
           # prob_x_female_given_r_v2 and prob_x_female_given_r3 are identical
           (prob_x_female_given_r_v2(r, k)*dpois(nblocks, (1 + (r+1)*genlen)/2^r) + 
         (1-prob_x_female_given_r_v2(r, k))*dpois(nblocks, (1 + r*genlen)/2^r))*recomb_meioses_pmf(2*k-1, r)
  }))
}, c('nblocks', 'k'))


prob_r_fullcousins <- Vectorize(function(r, k) {
  # how to arrange r females in the two lineages from the two daughters of 
  # and X full sib family k generations back
  if (k == 2) {
    if (r == 2) return(1)
    return(0)
  }
  if (!(r %in% rms_sharedanc(k, fullcousins=TRUE)))
    return(0)
  s <- seq(floor((k-2)/2), k-2)
  sum(choose(s+1, k-s-2)*choose(r-s-1, k-r+s))/fib(k)^2
}, 'r')

full_cousins_x_male <- Vectorize(function(nblocks, k, r, genlen=XGENLEN) {
  # Half cousin blocknum PMF given male shared ancestor and known r
  breakpoints <- seq(0, MAXBPOINTS)
  # this range is across all females; prob_r_given_male() returns 0 when not possible.
  sum(dpois(breakpoints, r*genlen)*dbinom(nblocks, breakpoints+1, 1/2^r))
}, c('nblocks', 'k'))


full_cousins_x_female <- Vectorize(
function(nblocks, k, r, genlen=XGENLEN) {
  # female as shared ancestor, fixed on her having two daughters. This is male case + 4 RMs
  breakpoints <- seq(0, MAXBPOINTS)
  sum(dpois(breakpoints, (r+2)*genlen)*dbinom(nblocks, breakpoints+1, 1/2^(r+1)))
}, c('nblocks', 'k'))


#' PMF for IBD segments shared between two full-cousins with ancestor k
#' generations ago
#'
#'
#' @param nblocks Number of blocks shared between two full-cousins k
#' generations ago.
#' @param k Number of generations until shared ancestor
#' @param genlen The genetic length (in Morgans) of the X chromosome.
full_cousins_x_blocknum_pmf <- Vectorize(function(nblocks, k, genlen=XGENLEN) {
  s <- seq(0, nblocks)
  recombs <- rms_sharedanc(k, fullcousins=TRUE)
  sum(sapply(recombs, function(r) {
     sum(full_cousins_x_male(s, k, r)*
         full_cousins_x_female(nblocks-s, k, r))*prob_r_fullcousins(r, k)
  }))
}, c('nblocks', 'k'))


half_cousins_x_blocknum_r <- Vectorize(function(nblocks, k, r, genlen=XGENLEN) {
  # P(N = n | R = r, k)
  breakpoints <- seq(0, MAXBPOINTS)
  recombs <- rms_sharedanc(k)
  stopifnot(r %in% recombs)
  (prob_x_female_given_r_v2(r, k)*sum(dpois(breakpoints, (r+1)*genlen)*
                                    dbinom(nblocks, breakpoints+1, 1/2^r)) + 
   (1-prob_x_female_given_r_v2(r, k))*sum(dpois(breakpoints, r*genlen)*
                                        dbinom(nblocks, breakpoints+1, 1/2^r)))
}, c('nblocks', 'k', 'r'))

half_cousins_blocknum_rm_ll <- Vectorize(function(nblocks, r, k) {
  # just a wrapper for half_cousins_x_blocknum_r  
  half_cousins_x_blocknum_r(nblocks, k, r)
}, c('nblocks', 'r', 'k'))


half_cousins_blocknum_rm_posterior <- Vectorize(function(nblocks, r, k) {
  # posterior P(R = r | N = n)P(R=r)/P(N=n)
  half_cousins_x_blocknum_r(nblocks, k, r)*recomb_meioses_pmf(2*k-1, r)/half_cousins_x_blocknum_pmf(nblocks, k)
}, c('nblocks', 'r', 'k'))

# quick sum-to-one sanity test.
local({
  k <- 5
  is_pmf(half_cousins_blocknum_rm_posterior(4, rms_sharedanc(k), k))
  is_pmf(half_cousins_blocknum_rm_posterior(5, rms_sharedanc(k), k))

  is_pmf(prob_r_given_female(0:30, k))
})

## probability of shared X ancestor
prob_x_shared_anc <- function(k) 0.5*(fib(k+1)/2^{k-1})^2 + 0.5*(fib(k)/2^{k-1})^2

## probability of x ancestor
prob_x_anc <- function(k) fib(k+2)/2^k

## probability of X shared anc | N = 0
prob_x_sharedanc_N0 <- Vectorize(function(k) {
  num <- prob_x_shared_anc(k)*half_cousins_x_blocknum_pmf(0, k)
  denom <- prob_x_shared_anc(k)*half_cousins_x_blocknum_pmf(0, k) + (1-prob_x_shared_anc(k))
  num/denom
}, 'k')

## probability of X anc | N = 0
prob_xanc_N0 <- Vectorize(function(k) {
  if (k == 1) return(NA)  # all present-day female's ancestors are x anc and N>1
  num <- prob_x_anc(k)*x_blocknum_pmf(0, k)
  denom <- prob_x_anc(k)*x_blocknum_pmf(0, k) + (1-prob_x_anc(k))
  num/denom
}, 'k')

