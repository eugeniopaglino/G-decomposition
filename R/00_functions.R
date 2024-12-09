# Here we have the key functions implementing the groupwise decomposition
# There are more functions that those described in the paper. Before
# realizing that computing the full decomposition was not too computationally
# intensive if coded efficiently, I explored ways of obtaining an approximate
# decomposition based on averaging the two external contributions (see paper)
# and then performing some additional adjustments.
# It turns out this is not straightforward and though you can identify good
# adjustments for the case with three groups things become more complicated
# with more groups. I still leave the functions here in case I find time to 
# do some more work.

# Load packages
library(tidyverse)
library(DemoDecomp)

# Compute difference in life-expectancy at birth between two populations.
# Takes as input a matrix of exposures and one of death counts. The set
# argument is the other groups we are adding, the baseline is the baseline
# group and the target is the group whose contribution we want to compute.
# The function automatically aggregates deaths and exposures and compute
# the correct mortality rates before feeding them into the function that
# compute life expectancy.

# NOTE: We can replace LTabr (from DemoDecomp) with any other function taking
# as input age-specific mortality rates and giving back some indicator.
# This could also be extended to fertility applications with some additional
# work to figure out how to adapt this function.
compute_contrib <- function(nEx_mat,nDx_mat,set,baseline,target_group) {
  set_with <- c(baseline,set)
  set_without <- set_with[set_with != target_group]
  
  nEx_with <- apply(nEx_mat[,set_with],M=1,sum)
  nDx_with <- apply(nDx_mat[,set_with],M=1,sum)
  nmx_with <- nDx_with/nEx_with
  e0_with <- LTabr(nmx_with)
  
  nEx_without <- apply(nEx_mat[,set_without],M=1,sum)
  nDx_without <- apply(nDx_mat[,set_without],M=1,sum)
  nmx_without <- nDx_without/nEx_without
  e0_without <- LTabr(nmx_without)
  
  contrib <- e0_with - e0_without
  contrib
}

# Compute the contribution of one group, return just the contribution
compute_group_contrib_simple <- function(nEx_mat,nDx_mat,baseline,target_group,G,type='full') {
  group_contrib <- 0
  
  # Here we need to compute the first "external" contribution (see paper)
  set_with <- c(baseline,target_group)
  set_without <- set_with[set_with != target_group]
  
  nEx_with <- apply(nEx_mat[,set_with],M=1,sum)
  nDx_with <- apply(nDx_mat[,set_with],M=1,sum)
  nmx_with <- nDx_with/nEx_with
  e0_with <- LTabr(nmx_with)
  
  nEx_without <- nEx_mat[,set_without] # It's not nice but these are a bit different
  nDx_without <- nDx_mat[,set_without]
  nmx_without <- nDx_without/nEx_without
  e0_without <- LTabr(nmx_without)
  
  # We set the weights and compute all other contributions except the last
  w <- if_else(type=='full',(1/G),1/2)
  contrib <- e0_with - e0_without
  group_contrib <- group_contrib + contrib*w
  
  if (type=='full') {
    for (n in 2:(G-1)) {
      sets <- combn(contrib_groups,n)
      for (s in 1:ncol(sets)) {
        if (target_group %in% sets[,s]) {
          set <- sets[,s]
          w <- (1/choose(G-1,n-1))*(1/G)
          contrib <- compute_contrib(nEx_mat,nDx_mat,set,baseline,target_group)
          group_contrib <- group_contrib + contrib*w
        }
      }
    }
  }
  
  # Here we need to compute the last "external" contribution (see paper)
  set <- contrib_groups
  w <- if_else(type=='full',(1/G),1/2)
  contrib <- compute_contrib(nEx_mat,nDx_mat,set,baseline,target_group)
  group_contrib <- group_contrib + contrib*w
  
  group_contrib
}

# Compute the contribution of one group, return all intermediate terms (useful for checks).
# This function is almost identical to the one above but returns more information.
# I kept them separate to avoid making a very long unique function hard to read.
compute_group_contrib_full <- function(nEx_mat,nDx_mat,baseline,target_group,G,type='full') {
  group_contribs <- c()
  weights <- c()
  ms <- c()
  
  # Here we need to compute the first "external" contribution (see paper)
  set_with <- c(baseline,target_group)
  set_without <- set_with[set_with != target_group]
  
  nEx_with <- apply(nEx_mat[,set_with],M=1,sum)
  nDx_with <- apply(nDx_mat[,set_with],M=1,sum)
  nmx_with <- nDx_with/nEx_with
  e0_with <- LTabr(nmx_with)
  
  nEx_without <- nEx_mat[,set_without] # It's not nice but these are a bit different
  nDx_without <- nDx_mat[,set_without]
  nmx_without <- nDx_without/nEx_without
  e0_without <- LTabr(nmx_without)
  
  # We set the weights and compute all other contributions except the last
  w <- if_else(type=='full',(1/G),1/2)
  contrib <- e0_with - e0_without
  group_contribs <- c(group_contribs,contrib)
  weights <- c(weights,w)
  ms <- c(ms,1)
  
  if (type=='full') {
    for (n in 2:(G-1)) {
      sets <- combn(contrib_groups,n)
      for (s in 1:ncol(sets)) {
        if (target_group %in% sets[,s]) {
          set <- sets[,s]
          w <- (1/choose(G-1,n-1))*(1/G)
          contrib <- compute_contrib(nEx_mat,nDx_mat,set,baseline,target_group)
          group_contribs <- c(group_contribs,contrib)
          weights <- c(weights,w)
          ms <- c(ms,n)
        }
      }
    }
  }
  
  # Here we need to compute the last "external" contribution (see paper)
  set <- contrib_groups
  w <- if_else(type=='full',(1/G),1/2)
  contrib <- compute_contrib(nEx_mat,nDx_mat,set,baseline,target_group)
  group_contribs <- c(group_contribs,contrib)
  weights <- c(weights,w)
  ms <- c(ms,G)
  
  data.frame(
    C = group_contribs,
    w = weights,
    m = ms
  )
}

# Combines contributions for all groups, performin the decomposition
groupwise_decomp <- function(nEx_mat,nDx_mat,baseline,contrib_groups,type='full',return_type='simple') {
  
  G <- length(c(baseline,contrib_groups)) - 1
  
  if (type=='full' & return_type=='full') {
    decomp_results <- map_dfr(
      contrib_groups,
      ~ compute_group_contrib_full(nEx_mat,nDx_mat,baseline,.x,G,type) %>%
        mutate(group=.x)
    )
  } else {
    group_contribs <- c()
    
    for (target_group in contrib_groups) {
      group_contrib <- compute_group_contrib_simple(nEx_mat,nDx_mat,baseline,target_group,G,type)
      group_contribs <- c(group_contribs,group_contrib)
    }
    
    decomp_results <- data.frame(
      group = contrib_groups,
      C = group_contribs
    )
  }
  
  if (type=='adjusted') {
    nEx_B <- nEx_mat[,baseline]
    nDx_B <- nDx_mat[,baseline]
    nmx_B <- nDx_B/nEx_B
    e0_B <- LTabr(nmx_B)
    
    nEx_T <- apply(nEx_mat,M=1,sum)
    nDx_T <- apply(nDx_mat,M=1,sum)
    nmx_T <- nDx_T/nEx_T
    e0_T <- LTabr(nmx_T)
    
    actual_diff <- e0_T - e0_B
    decomp_diff <- decomp_results %>% pull(C) %>% sum()
    
    decomp_results <- decomp_results %>%
      mutate(C = C + (actual_diff-decomp_diff)/G)
  }
  
  decomp_results
}

total_e0 <- function(pars,groups,ages) {
  
  G <- length(unique(groups))
  A <- length(unique(ages))
  N <- G*A
  
  split_pars <- split(pars,rep(c('nEx','nDx'),each=N))
  nEx <- split_pars$nEx
  nDx <- split_pars$nDx
  
  lt_data <- tibble(
    nEx=nEx,
    nDx=nDx,
    group=groups,
    age=ages
  )
  
  nmx <- lt_data %>%
    group_by(age) %>%
    summarise(
      nEx=sum(nEx),
      nDx=sum(nDx)
    ) %>%
    ungroup() %>%
    mutate(nmx=nDx/nEx) %>%
    pull(nmx)
  
  e0 <- LTabr(nmx)
  
  e0
}

# DEPRECATED
compute_group_contrib <- function(nEx_mat,nDx_mat,baseline,target_group,G,type='full') {
  group_contribs <- c()
  weights <- c()
  
  set_with <- c(baseline,target_group)
  set_without <- set_with[set_with != target_group]
  
  nEx_with <- apply(nEx_mat[,set_with],M=1,sum)
  nDx_with <- apply(nDx_mat[,set_with],M=1,sum)
  nmx_with <- nDx_with/nEx_with
  e0_with <- LTabr(nmx_with)
  
  nEx_without <- nEx_mat[,set_without] # It's not nice but these are a bit different
  nDx_without <- nDx_mat[,set_without]
  nmx_without <- nDx_without/nEx_without
  e0_without <- LTabr(nmx_without)
  
  w <- if_else(type=='full',(1/G),1/2)
  contrib <- e0_with - e0_without
  group_contribs <- c(group_contribs,contrib)
  weights <- c(weights,w)
  
  if (type=='full') {
    for (n in 2:(G-1)) {
      sets <- combn(contrib_groups,n)
      for (s in 1:ncol(sets)) {
        if (target_group %in% sets[,s]) {
          set <- sets[,s]
          w <- (1/choose(G-1,n-1))*(1/G)
          contrib <- compute_contrib(nEx_mat,nDx_mat,set,baseline,target_group)
          group_contribs <- c(group_contribs,contrib)
          weights <- c(weights,w)
        }
      }
    }
  }
  
  set <- contrib_groups
  w <- if_else(type=='full',(1/G),1/2)
  contrib <- compute_contrib(nEx_mat,nDx_mat,set,baseline,target_group)
  group_contribs <- c(group_contribs,contrib)
  weights <- c(weights,w)
  
  data.frame(
    C = group_contribs,
    w = weights
  )
}
