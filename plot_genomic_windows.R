#!/usr/bin/env Rscript

# (c) 2023 Angel G. Rivera-Colon
# Plot the distribution of coverage and variant sites along a genome

# --- User-defined paths --- #
# TODO: Argparse

# Working directory
work_dir <- '/path/to/workdir'
# Chromosome info table
chr_info_f <- '/path/to/chromosomes.tsv'
# Variant sites:
# File
var_sites_f <- '/path/to/variant_site_tally.tsv'
# Title for variant sites plot
title_var <- 'Variant Sites'
# Output directory
var_out <- '/path/to/output'

# Depth of coverage
# File
coverage_f <- '/path/to/window_average_coverage.tsv'
# Tilte for coverage plot
title_cov <- 'Depth of Coverage'
# Output directory
cov_out <- '/path/to/output'

# -------------------------- #

# Libraries
library(ggplot2)
library(RColorBrewer)
library(scales)

# Constants
SCALE <- 1e6 # 1 Mbp
WIN.SIZE <- 250e3 # 250 Kbp

# Load variant sites
load.var.sites <- function(varSites.f, scale=SCALE){
  varSites.df <- read.table(varSites.f, header=TRUE)
  varSites.df$Mb <- varSites.df$Pos/scale
  return(varSites.df)
}

# Load coverage file
load.coverage.f <- function(coverage.f, scale=SCALE){
  coverage.df <- read.table(coverage.f, header=TRUE)
  coverage.df$Mb <- coverage.df$Pos/scale
  return(coverage.df)
}

# Load the chromosome info
load.chr.info <- function(chrInfo.f, scale=SCALE){
  chrInfo.df <- read.delim(chrInfo.f)
  # Convert to Mbp to simplify X axis
  chrInfo.df$Mb <- chrInfo.df$length/SCALE
  # Create an index on the Y axis based on the row of the given chr
  chrInfo.df$Y <- (nrow(chrInfo.df)+1) - row(chrInfo.df)[,1]
  return(chrInfo.df)
}

# Function to add the corresponsing Y axis value to df
add.chr.y <- function(df, chrInfo.df){
  # Create a named list
  vals <- chrInfo.df$Y
  scaf <- chrInfo.df$scaffold
  names(vals) <- scaf
  # Add those values to the DF
  df$Y <- vals[df$Chrom]
  return(df)
}

# Find the step of the x-axis
find.step <- function(max.mb){
  step <- 10
  if (max.mb < 10){
    step <- 1
  } else if (max.mb < 20) {
    step <- 2
  } else if (max.mb < 50) {
    step <- 5
  } else if (max.mb < 100) {
    step <- 10
  } else if (max.mb < 250) {
    step <- 20
  } else {
    step <- 50
  }
  return(step)
}

# Find the x-axis max boundary
find.x.boundary <- function(mb.vector, pad=0.01){
  max.mb <- max(mb.vector)
  boundary <- max.mb + (max.mb*pad)
  return(boundary)
}

# Plot a chromosome boundary box
plot.chr.boundary <- function(fig, max.mb, y.val, box.height=0.6,
                              border.col='grey25', fill.col=NA){
  # Prepare input data
  bh <- box.height/2
  x <- c(0, max.mb, max.mb, 0)
  y <- c((y.val+bh), (y.val+bh), (y.val-bh), (y.val-bh))
  df <- data.frame(x=x, y=y)
  # Draw the rectangle
  fig <- fig +
    geom_rect(data=NULL,
              aes(xmin=0, xmax=max.mb, ymin=(y.val-bh), ymax=(y.val+bh)),
              fill=fill.col,
              color=border.col)
  return(fig)
}

# Draw the chromosome background squares
chrom.squares <- function(fig, chrInfo.df, box.height=0.65, border.col=NA, fill.col='grey99'){
  for (r in 1:nrow(chrInfo)){
    row <- chrInfo.df[r,]
    fig <- plot.chr.boundary(fig, row$Mb, row$Y, box.height=box.height,
                             border.col=border.col, fill.col=fill.col)
  }
  return(fig)
}

# Prepare the base plot
make.base.plot <- function(chrInfo.df, name, box.height=0.65, pad.x=0.01, pad.y=0.5){
  # Calculate coordinate boundaries
  max.x <- find.x.boundary(chrInfo.df$Mb, pad.x)
  min.x <- 0#-(max.x*pad.x)
  step.x <- find.step(max.x)
  max.y <- max(chrInfo.df$Y)
  min.y <- min(chrInfo.df$Y)
  labels = rev(chrInfo.df$label)
  # Base plot
  fig <- ggplot(chrInfo, aes(x=Mb, y=Y)) +
    labs(x='Position (Mbp)',
         y='Scaffold',
         title=name) +
    scale_x_continuous(limits = c(min.x, max.x),
                       breaks=seq(0,(max.x+step.x),step.x),
                       expand = c(0.025, 0)) +
    scale_y_continuous(limits = c((min.y-pad.y), (max.y+pad.y)),
                       breaks=seq(min.y, max.y, 1),
                       labels=labels,
                       expand=c(0.01,0.01)) +
    # theme_linedraw() +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5))
  return(fig)
}


# Make variant sites plot
variant.sites.plot <- function(chrInfo.df, varSites.df, title, outdir='.'){
  # Make base plot
  fig <- make.base.plot(chrInfo.df, title)
  # Add the chromosome boundaries
  fig <- chrom.squares(fig, chrInfo.df, fill.col='grey99')
  # Plot the variant sites
  fig <- plot.var.sites(fig, varSites.df, chrInfo.df, 10)
  # Save the plot
  plot_f <- paste(outdir, '/variant_sites.pdf', sep='')
  ggsave(plot_f, plot=fig, width=8, height=6)
  return(fig)  
}


# Plot the variant sites
plot.var.sites <- function(figure, varSites.df, chrInfo.df,
                           nBins=15, bh=0.325){
   # Add the Y to the variant sites
  varSites.df <- add.chr.y(varSites.df, chrInfo.df)
  # Sort by the value count
  varSites.df <- varSites.df[order(varSites.df$NumSitesWindow),]
  # Calculate the distance boundaries
  min.legend <- 0
  # max.legend <- max(varSites.df$NumSitesWindow)
  max.legend <- as.numeric(quantile(varSites.df$NumSitesWindow, 0.9995))
  # Add to plot
  figure <- figure +
    geom_segment(data=varSites.df, aes(x=Mb, xend=Mb,
                                    y=(Y-bh), yend=(Y+bh),
                                    # color=FreqSitesWindow),
                                    color=NumSitesWindow),
                 alpha=1.0,
                 size=0.1) +
    # Use the binned gradient
    scale_color_binned(type='viridis',
                       name='Variant sites\nper window',
                       n.breaks=nBins,
                       # breaks=bins,
                       limits=c(min.legend, max.legend),
                       alpha=0.9) +
    theme(legend.box.background=element_rect(),
          legend.text = element_text(size=8),
          legend.key.height= unit(1, 'cm'))
  figure <- chrom.squares(figure, chrInfo.df, border.col='grey25', fill.col=NA)
  return(figure)
}


# Plot the site coverage
plot.site.coverage <- function(figure, coverage.df, chrInfo.df,
                               nBins=15, bh=0.325){
  # Add the Y to the coverage df
  coverage.df <- add.chr.y(coverage.df, chrInfo.df)
  # Sort by the value count
  coverage.df <- coverage.df[order(coverage.df$MeanCov),]
  # Calculate the distance boundaries
  min.legend <- 0
  # max.legend <- max(coverage.df$MeanCov)
  max.legend <- as.numeric(quantile(coverage.df$MeanCov, 0.9995))
  # Add to plot
  figure <- figure +
    geom_segment(data=coverage.df, aes(x=Mb, xend=Mb,
                                       y=(Y-bh), yend=(Y+bh),
                                       color=MeanCov),
                 alpha=1.0,
                 size=0.1) +
    # Use the binned gradient
    scale_color_binned(type='viridis',
                       name='Mean depth\nof coverage',
                       n.breaks=nBins,
                       # breaks=bins,
                       limits=c(min.legend, max.legend),
                       alpha=0.9) +
    theme(legend.box.background=element_rect(),
          legend.text = element_text(size=8),
          legend.key.height= unit(1, 'cm'))
  figure <- chrom.squares(figure, chrInfo.df, border.col='grey25', fill.col=NA)
  return(figure)
}

# Make coverage plot
coverage.plot <- function(chrInfo.df, coverage.df, title, outdir='.'){
  # Make base plot
  fig <- make.base.plot(chrInfo.df, title)
  # Add the chromosome boundaries
  fig <- chrom.squares(fig, chrInfo.df, fill.col='grey99')
  # Plot the coverage
  fig <- plot.site.coverage(fig, coverage.df, chrInfo.df)
  # Save the plot
  plot_f <- paste(outdir, './depth_of_coverage.pdf', sep='')
  ggsave(plot_f, plot=fig, width=8, height=6)
  return(fig)  
}

#
# Run:
#

# Prepare environment and load data
setwd(work_dir)
chrInfo  <- load.chr.info(chr_info_f)
varSites <- load.var.sites(var_sites_f)
coverage <- load.coverage.f(coverage_f) 

# Plot the variant sites
varSites.fig <- variant.sites.plot(chrInfo, varSites,
                                   title=title_var, outdir = var_out)

# Plot the coverage
coverage.fig <- coverage.plot(chrInfo, coverage, 
                              title=title_cov, outdir = cov_out)



