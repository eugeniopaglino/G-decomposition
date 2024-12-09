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

# The functions below replicate those in the 00_functions.R file but extend
# them by nesting generalized decomposition within the groupwise decomposition
# to obtain group-specific elements by an additional variable (here age).

# Load packages
library(tidyverse)
library(DemoDecomp)

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

compute_contrib_by_age <- function(nEx_mat,nDx_mat,set,baseline,target_group) {
  set_with <- c(baseline,set)
  set_without <- set_with[set_with != target_group]
  
  nEx_with <- apply(nEx_mat[,set_with],M=1,sum)
  nDx_with <- apply(nDx_mat[,set_with],M=1,sum)
  nmx_with <- nDx_with/nEx_with
  
  nEx_without <- apply(nEx_mat[,set_without],M=1,sum)
  nDx_without <- apply(nDx_mat[,set_without],M=1,sum)
  nmx_without <- nDx_without/nEx_without
  
  # Decomposing the group contribution with the DemoDecomp methods
  contrib <- stepwise_replacement(
    func = LTabr,
    direction = 'both',
    pars1 = nmx_without,
    pars2 = nmx_with
  )
}

compute_group_contrib_simple <- function(nEx_mat,nDx_mat,baseline,target_group,G,type='full') {
  group_contrib <- 0
  
  set_with <- c(baseline,target_group)
  set_without <- set_with[set_with != target_group]
  
  nEx_with <- apply(nEx_mat[,set_with],M=1,sum)
  nDx_with <- apply(nDx_mat[,set_with],M=1,sum)
  nmx_with <- nDx_with/nEx_with
  
  nEx_without <- nEx_mat[,set_without] # It's not nice but these are a bit different
  nDx_without <- nDx_mat[,set_without]
  nmx_without <- nDx_without/nEx_without
  
  w <- if_else(type=='full',(1/G),1/2)
  
  contrib <- stepwise_replacement(
    func = LTabr,
    direction = 'both',
    pars1 = nmx_without,
    pars2 = nmx_with
  )
  
  group_contrib <- group_contrib + contrib*w
  
  if (type=='full') {
    for (n in 2:(G-1)) {
      sets <- combn(contrib_groups,n)
      for (s in 1:ncol(sets)) {
        if (target_group %in% sets[,s]) {
          set <- sets[,s]
          w <- (1/choose(G-1,n-1))*(1/G)
          contrib <- compute_contrib_by_age(nEx_mat,nDx_mat,set,baseline,target_group)
          group_contrib <- group_contrib + contrib*w
        }
      }
    }
  }
  
  set <- contrib_groups
  w <- if_else(type=='full',(1/G),1/2)
  contrib <- compute_contrib_by_age(nEx_mat,nDx_mat,set,baseline,target_group)
  group_contrib <- group_contrib + contrib*w
  
  group_contrib
}

compute_groupwise_decomp <- function(nEx_mat,nDx_mat,baseline,contrib_groups,type='full') {
  
  G <- length(c(baseline,contrib_groups)) - 1
  
  group_contribs <- c()
  
  for (target_group in contrib_groups) {
    group_contrib <- compute_group_contrib_simple(nEx_mat,nDx_mat,baseline,target_group,G,type)
    group_contribs <- c(group_contribs,group_contrib)
  }
  
  groupwise_decomp <- tibble(
    group = rep(contrib_groups,each=length(unique(ages))),
    age = rep(unique(ages),length(contrib_groups)),
    C = group_contribs
  )
  
  groupwise_decomp
}
