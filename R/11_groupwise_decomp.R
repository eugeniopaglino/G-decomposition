# Loading necessary packages
library(here)
library(colorspace)
library(ggrepel)
library(tidyverse)

# Do not rely on this to completely clean your environment
# Better to do a full restart of R before running
rm(list=ls())

i_am('R/11_groupwise_decomp.R')

in_dir <- here('data')
out_dir <- here('output')
fig_dir <- here('figures')

source(here('R','00_functions.R'))
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

# Run decomposition (with simple return)
groupwise_full <- pmap_dfr(
  list(nEx_mat=nEx_mats,
       nDx_mat=nDx_mats,
       sex=mort_data_nested$sex,
       year=mort_data_nested$year),
  function(nEx_mat,nDx_mat,year,sex) {
    groupwise_decomp(nEx_mat,nDx_mat,baseline,contrib_groups,type='full') %>%
      mutate(year=year,sex=sex)
  }
)

# Run decomposition (returning all sub-contributions)
groupwise_full_full <- pmap_dfr(
  list(nEx_mat=nEx_mats,
       nDx_mat=nDx_mats,
       sex=mort_data_nested$sex,
       year=mort_data_nested$year),
  function(nEx_mat,nDx_mat,year,sex) {
    groupwise_decomp(nEx_mat,nDx_mat,baseline,contrib_groups,type='full',return_type = 'full') %>%
      mutate(year=year,sex=sex)
  }
)

# Run approximate decomposition with just a simple average of external contributions
groupwise_avg <- pmap_dfr(
  list(nEx_mat=nEx_mats,
       nDx_mat=nDx_mats,
       sex=mort_data_nested$sex,
       year=mort_data_nested$year),
  function(nEx_mat,nDx_mat,year,sex) {
    groupwise_decomp(nEx_mat,nDx_mat,baseline,contrib_groups,type='average') %>%
      mutate(year=year,sex=sex)
  }
)

# Run approximate decomposition with just a simple average of external contributions
# but with the additional adjustment forcing the contributions to sum to the total
# difference
groupwise_adj <- pmap_dfr(
  list(nEx_mat=nEx_mats,
       nDx_mat=nDx_mats,
       sex=mort_data_nested$sex,
       year=mort_data_nested$year),
  function(nEx_mat,nDx_mat,year,sex) {
    groupwise_decomp(nEx_mat,nDx_mat,baseline,contrib_groups,type='adjusted') %>%
      mutate(year=year,sex=sex)
  }
)

# Combine decomposition results with different methds for comparisons
decomp_comparisons <- groupwise_full_full %>%
  group_by(year,sex,group,m) %>%
  summarise(Cm = sum(C*w)/(sum(w))) %>%
  ungroup() %>%
  left_join(groupwise_full %>% rename(Cfull=C),by=c('year','sex','group')) %>%
  left_join(groupwise_avg %>% rename(Cavg=C),by=c('year','sex','group')) %>%
  left_join(groupwise_adj %>% rename(Cadj=C),by=c('year','sex','group'))

# Visualise the results
comparison_plot <- decomp_comparisons %>%
  mutate(
    group=fct_relevel(
      group,
      c("Large Fringe Metro","Medium Metro","Small Metro",
        "Micropolitan (Nonmetro)","NonCore (Nonmetro)")
    )
  ) %>%
  ggplot() +
  geom_point(mapping=aes(y=group,x=Cm,shape='Ordered Contributions')) +
  geom_text_repel(mapping=aes(y=group,x=Cm,label=m),nudge_y = 0.3,min.segment.length = 0) +
  geom_point(mapping=aes(y=group,x=Cfull,shape='Exact')) +
  geom_point(mapping=aes(y=group,x=Cavg,shape='External Average')) +
  geom_point(mapping=aes(y=group,x=Cadj,shape='External Average + Adjustment')) +
  scale_shape_manual(values=15:18,breaks = c('Exact','External Average','External Average + Adjustment','Ordered Contributions')) +
  scale_y_discrete(expand = expansion(mult = 0.25)) +
  facet_grid(year~sex,switch = 'y') +
  labs(shape='',
       x='Contribution',
       y='') +
  theme_bw() +
  theme(legend.position = 'bottom',
        strip.text.y.left = element_text(angle=0),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

comparison_plot

svg(here(fig_dir,'11_comparison_plot.svg'),width=8,height=8)
comparison_plot
dev.off()

groupwise_full <- groupwise_full %>%
  mutate(
    group=fct_relevel(
      group,
      c("Large Fringe Metro","Medium Metro","Small Metro",
        "Micropolitan (Nonmetro)","NonCore (Nonmetro)")
    )
  ) %>%
  arrange(sex,group)

# Create table with results from the full decomposition
groupwise_full %>%
  pivot_wider(names_from = year,values_from = C) %>%
  group_by(sex) %>%
  gt() %>%
  cols_label(
    group='Metro Category',
    sex=''
  ) %>%
  fmt_number(
    `2017`:`2019`,
    decimals = 3
  ) %>%
  tab_spanner(
    `2017`:`2019`,
    label='Contribution to National Life Expectancy'
  ) %>%
  gtsave(.,filename = here(fig_dir,'11_contributions_table.docx'))

# Visualize contributions
contributions_plots <- groupwise_full %>%
  ggplot() +
  geom_hline(yintercept=0) +
  geom_col(mapping=aes(x=year,y=C,fill=group)) +
  scale_fill_discrete_sequential(palette = 'SunsetDark') +
  facet_nested(~sex) +
  labs(x='',
       fill='',
       y='Contribution to The Difference in Life Expectancy\nbetween Large Central Metros and the Whole US') +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle=90))

contributions_plots

svg(here(fig_dir,'11_contributions_plots.svg'),width=4,height=4)
contributions_plots
dev.off()

# Set up data to replicate the decomposition results using the general
# decomposition methods as implemented in DemoDecomp

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

# Performing the decomposition with the line-integral decomposition method.
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

# As implemented, there are separate components for deaths and exposures
horiuchi_decomp <- horiuchi_decomp_by_age %>%
  filter(group!=baseline) %>%
  group_by(year,sex,group) %>%
  summarise(C=sum(C)) %>%
  ungroup()

# Performing the decomposition with the life table response experiment 
# decomposition method.
ltre_decomp_by_age <- pmap_dfr(
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
      C = ltre(
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

ltre_decomp <- ltre_decomp_by_age %>%
  filter(group!=baseline) %>%
  group_by(year,sex,group) %>%
  summarise(C=sum(C)) %>%
  ungroup()

# Performing the decomposition with the stepwise decomposition method.
stepwise_decomp_by_age <- pmap_dfr(
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
      C = stepwise_replacement(
        func = total_e0,
        direction = 'both',
        pars1 = pars1,
        pars2 = pars2,
        groups = groups,
        ages = ages
      )
    )
  }
)

stepwise_decomp <- stepwise_decomp_by_age %>%
  filter(group!=baseline) %>%
  group_by(year,sex,group) %>%
  summarise(C=sum(C)) %>%
  ungroup()

# Combine results for comparison
decomp_methods_comparisons <- groupwise_full %>%
  left_join(horiuchi_decomp %>% rename(Choriuchi=C),by=c('year','sex','group')) %>%
  left_join(ltre_decomp %>% rename(Cltre=C),by=c('year','sex','group')) %>%
  left_join(stepwise_decomp %>% rename(Cstepwise=C),by=c('year','sex','group'))

# Plot the results
decomp_methods_comparison_plot <- decomp_methods_comparisons %>%
  mutate(
    group=fct_relevel(
      group,
      c("Large Fringe Metro","Medium Metro","Small Metro",
        "Micropolitan (Nonmetro)","NonCore (Nonmetro)")
    )
  ) %>%
  ggplot() +
  geom_point(mapping=aes(y=group,x=C,shape='Groupwise')) +
  geom_point(mapping=aes(y=group,x=Cstepwise,shape='Stepwise')) +
  geom_point(mapping=aes(y=group,x=Choriuchi,shape='Horiuchi')) +
  geom_point(mapping=aes(y=group,x=Cltre,shape='LTRE')) +
  facet_grid(year~sex,switch = 'y') +
  scale_shape_manual(values=21:24) +
  labs(shape='',
       x='Contribution',
       y='') +
  theme_bw() +
  theme(legend.position = 'bottom',
        strip.text.y.left = element_text(angle=0),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

decomp_methods_comparison_plot

svg(here(fig_dir,'11_decomp_methods_comparison_plot.svg'),width=6,height=6)
decomp_methods_comparison_plot
dev.off()

# Repeat plot removing stepwise which does not seem to work
decomp_methods_comparison_no_stepwise_plot <- decomp_methods_comparisons %>%
  mutate(
    group=fct_relevel(
      group,
      c("Large Fringe Metro","Medium Metro","Small Metro",
        "Micropolitan (Nonmetro)","NonCore (Nonmetro)")
    )
  ) %>%
  ggplot() +
  geom_point(mapping=aes(y=group,x=C,shape='Groupwise')) +
  geom_point(mapping=aes(y=group,x=Choriuchi,shape='Horiuchi')) +
  geom_point(mapping=aes(y=group,x=Cltre,shape='LTRE')) +
  facet_grid(year~sex,switch = 'y') +
  scale_shape_manual(values=21:23) +
  labs(shape='',
       x='Contribution',
       y='') +
  theme_bw() +
  theme(legend.position = 'bottom',
        strip.text.y.left = element_text(angle=0),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

decomp_methods_comparison_no_stepwise_plot

svg(here(fig_dir,'11_decomp_methods_comparison_no_stepwise_plot.svg'),width=6,height=6)
decomp_methods_comparison_no_stepwise_plot
dev.off()

# Run execution benchmarking 
res <- microbenchmark::microbenchmark(
  
  groupwise_full <- pmap_dfr(
    list(
      nEx_mat=nEx_mats,
      nDx_mat=nDx_mats,
      sex=mort_data_nested$sex,
      year=mort_data_nested$year
      ),
    function(nEx_mat,nDx_mat,year,sex) {
      groupwise_decomp(nEx_mat,nDx_mat,baseline,contrib_groups,type='full') %>%
        mutate(year=year,sex=sex)
    }
  )
  ,
  
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
  
  ,
  
  ltre_decomp_by_age <- pmap_dfr(
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
        C = ltre(
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
  
  ,
  times = 5
)

save(res,file=here(out_dir,'11_group_specific_res.RData'))
