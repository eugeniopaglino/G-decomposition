# Loading necessary packages
library(here)
library(colorspace)
library(ggrepel)
library(tidyverse)

# Do not rely on this to completely clean your environment
# Better to do a full restart of R before running
rm(list=ls())

i_am('R/12_groupwise_decomp_extended.R')

in_dir <- here('data')
out_dir <- here('output')
fig_dir <- here('figures')

source(here('R','00_functions_extended.R'))
mort_data <- read_csv(here(out_dir,'mort_data_clean.csv'))

mort_data_nested <- mort_data %>%
  rename(group=urbanization) %>%
  group_by(year,sex) %>%
  nest()

mort_data_nested_data <- mort_data_nested$data
baseline <- 'Large Central Metro'
groups <- mort_data_nested_data[[1]] %>% pull(group)
ages <- mort_data_nested_data[[1]] %>% pull(x)
contrib_groups <- unique(groups)[unique(groups)!=baseline]

# Set up matrices of age-specific exposures by group
nEx_mats <- map(
  mort_data_nested_data, 
  ~ .x %>%
    select(group,x,nEx) %>%
    pivot_wider(values_from = nEx, names_from = group) %>%
    select(-x) %>%
    as.matrix()
)

# Set up matrices of age-specific deaths by group
nDx_mats <- map(
  mort_data_nested_data, 
  ~ .x %>%
    select(group,x,nDx) %>%
    pivot_wider(values_from = nDx, names_from = group) %>%
    select(-x) %>%
    as.matrix()
)

groupwise_decomp_full <- pmap_dfr(
  list(nEx_mat=nEx_mats,
       nDx_mat=nDx_mats,
       sex=mort_data_nested$sex,
       year=mort_data_nested$year),
  function(nEx_mat,nDx_mat,year,sex) {
    compute_groupwise_decomp(nEx_mat,nDx_mat,baseline,contrib_groups,type='full') %>%
      mutate(year=year,sex=sex)
  }
)

contributions_by_age <- groupwise_decomp_full %>%
  mutate(
    group=fct_relevel(
      group,
      c("Large Fringe Metro","Medium Metro","Small Metro",
        "Micropolitan (Nonmetro)","NonCore (Nonmetro)")
    ),
    ageStr = glue('{age}-{age+4}'),
    ageStr = case_when(
      age==0 ~ '<1',
      age==1 ~ '1-4',
      age==85 ~ '85+',
      T ~ ageStr
    ),
    ageStr = fct_reorder(ageStr,age)
  ) %>%
  filter(year==2019) %>%
  ggplot() +
  geom_col(aes(y=ageStr, x=C, fill=group),width=0.8) +
  scale_fill_discrete_sequential(palette='SunsetDark') +
  facet_wrap(~sex,ncol=2) +
  labs(y='',
       x='Contribution to the Life Expectancy Gap between the US and Large Central Metros (Years)',
       fill='') +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  guides(
    fill=guide_legend(nrow=2)
  )

contributions_by_age

svg(here('figures','12_contributions_by_age.svg'),width=7,height=6)
contributions_by_age
dev.off()

nExs1 <- map(
  mort_data_nested_data,
  ~ .x %>%
    mutate(nEx=if_else(group==baseline,nEx,0)) %>%
    pull(nEx)
)

nDxs1 <- map(
  mort_data_nested_data,
  ~ .x %>%
    mutate(nDx=if_else(group==baseline,nDx,0)) %>%
    pull(nDx)
)

nExs2 <- map(
  mort_data_nested_data,
  ~ .x %>%
    pull(nEx)
)

nDxs2 <- map(
  mort_data_nested_data,
  ~ .x %>%
    pull(nDx)
)

# Decomposing this difference with the two methods
horiuchi_decomp_by_age <- pmap_dfr(
  list(
    nEx1 = nExs1,
    nEx2 = nExs2,
    nDx1 = nDxs1,
    nDx2 = nDxs2,
    sex = mort_data_nested$sex,
    year = mort_data_nested$year
  ),
  function(nEx1,nEx2,nDx1,nDx2,sex,year) {
    tibble(
      sex = sex,
      year = year,
      age = rep(ages,2),
      group = rep(groups,2),
      pars1 = c(nEx1,nDx1),
      pars2 = c(nEx2,nDx2),
      C = horiuchi(
        func = total_e0,
        pars1 = pars1,
        pars2 = pars2,
        groups = groups,
        ages = ages,
        N=5
      )
    )
  }
)

# Replicate with Horiuchi
horiuchi_contribs_by_age <- horiuchi_decomp_by_age %>%
  filter(group!=baseline) %>%
  group_by(year,sex,group,age) %>%
  summarise(C=sum(C)) %>%
  ungroup() %>%
  mutate(
    group=fct_relevel(
      group,
      c("Large Fringe Metro","Medium Metro","Small Metro",
        "Micropolitan (Nonmetro)","NonCore (Nonmetro)")
    ),
    ageStr = glue('{age}-{age+4}'),
    ageStr = case_when(
      age==0 ~ '<1',
      age==1 ~ '1-4',
      age==85 ~ '85+',
      T ~ ageStr
    ),
    ageStr = fct_reorder(ageStr,age)
  ) %>%
  filter(year==2019) %>%
  ggplot() +
  geom_col(aes(y=ageStr, x=C, fill=group),width=0.8) +
  scale_fill_discrete_sequential(palette='SunsetDark') +
  facet_wrap(~sex,ncol=2) +
  labs(y='',
       x='Contribution to the Life Expectancy Gap between the US and Large Central Metros (Years)',
       fill='') +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  guides(
    fill=guide_legend(nrow=2)
  )

horiuchi_contribs_by_age

svg(here('figures','12_horiuchi_contribs_by_age.svg'),width=7,height=6)
horiuchi_contribs_by_age
dev.off()

# Run execution benchmark agains Horiuchi
res <- microbenchmark::microbenchmark(
  groupwise_decomp_full <- pmap_dfr(
    list(nEx_mat=nEx_mats,
         nDx_mat=nDx_mats,
         sex=mort_data_nested$sex,
         year=mort_data_nested$year),
    function(nEx_mat,nDx_mat,year,sex) {
      compute_groupwise_decomp(nEx_mat,nDx_mat,baseline,contrib_groups,type='full') %>%
        mutate(year=year,sex=sex)
    }
  ),
  horiuchi_decomp_by_age <- pmap_dfr(
    list(
      nEx1 = nExs1,
      nEx2 = nExs2,
      nDx1 = nDxs1,
      nDx2 = nDxs2,
      sex = mort_data_nested$sex,
      year = mort_data_nested$year
    ),
    function(nEx1,nEx2,nDx1,nDx2,sex,year) {
      tibble(
        sex = sex,
        year = year,
        age = rep(ages,2),
        group = rep(groups,2),
        pars1 = c(nEx1,nDx1),
        pars2 = c(nEx2,nDx2),
        C = horiuchi(
          func = total_e0,
          pars1 = pars1,
          pars2 = pars2,
          groups = groups,
          ages = ages,
          N=5
        )
      )
    }
  ),
  times = 5
)

save(res,file=here(out_dir,'12_age_specific_res.RData'))
