# Retrospective analysis of model dynamics


# Load libraries ----------------------------------------------------------

# data wrangling
library(tidyverse)
library(lubridate)

# modeling with pomp
library(pomp)

# plotting
library(ggthemes)
library(cowplot)
library(geofacet)

# set plotting theme
ggplot2::theme_set(ggthemes::theme_few(base_size = 12))

# source helper functions for derived quantities
source("../code/result-exploration/getReff.R")



# Extract latent trends and compute Re ------------------------------------

all_files <- list.files("../output")
all_mif_files <- all_files[grep("mif_res_2.RDS", all_files)]
all_filter_files <- all_files[grep("filtered_states.RDS", all_files)]
all_pomp_files <- all_files[grep("pomp_model.RDS", all_files)]

all_states_filters <- tibble()
all_states_data <- tibble()
 
for(i in 1:length(all_mif_files)) {
  this_mif_file <- all_mif_files[i]
  this_state <- str_split(this_mif_file, "-", simplify = TRUE)[1]
  this_filter_file <- all_filter_files[grep(this_state, all_filter_files)]
  this_pomp_file <- all_pomp_files[grep(this_state, all_pomp_files)]
  
  this_fit <- readRDS(paste0("../output/", this_mif_file))
  this_filter <- readRDS(paste0("../output/", this_filter_file))
  this_pomp <- readRDS(paste0("../output/", this_pomp_file))
  
  this_data <- this_pomp$pomp_data %>%
    filter(time > 0)
  
  this_population <- this_filter %>%
    filter(time == 1, pf_rep == 1) %>%
    pull(S)
  
  this_mobility <- this_pomp$pomp_covar@table %>%
    t() %>%
    as_tibble() %>%
    dplyr::select(rel_beta_change) %>%
    mutate(time = this_pomp$pomp_covar@times) %>%
    rename("mobility" = rel_beta_change)
  
  this_params <- this_fit$partable_natural %>%
    rownames_to_column(var = "parameter") %>%
    pivot_longer(cols = X1:X100) %>%
    dplyr::select(-is_fitted) %>%
    pivot_wider(names_from = parameter, values_from = value) %>%
    filter(LogLik == max(LogLik))
  
  # calculate latent trend from basis functions and bspline coefs
  beta <- this_params[names(this_params)[startsWith(names(this_params), "b")]]
  beta <- beta %>% dplyr::select(-beta_s, -base_detect_frac)
  beta <- as.numeric(beta)
  B <- pomp::bspline.basis(
    x = this_pomp$pomp_data$time,
    degree = 3,
    nbasis = length(beta)
  )
  mu <- B %*% beta  # expectation on estimation scale
  trend <- exp(mu) / (1 + exp(mu))

  # make it a dataframe for joining
  trend_df <- tibble(
    date = unique(this_data$date),
    latent_trend = trend[2:length(trend)]  # ignore time 0
  ) 
  
  this_filter <- this_filter %>%
    left_join(this_data %>% dplyr::select(date, time), by = "time") %>%
    left_join(trend_df, by = "date") %>%
    # mutate(latent_trend = exp(trendO) / (1+exp(trendO))) %>%
    left_join(this_mobility, by = "time") %>%
    mutate(combined_trend = mobility*latent_trend) %>%
    mutate(N = this_population) %>%
    mutate(cumulative_infections = N - S) %>%
    group_by(pf_rep) %>%
    arrange(pf_rep, time) %>%
    group_by(pf_rep) %>%
    mutate(cumulative_cases = cumsum(C_new),
           cumulative_deaths = cumsum(D_new)) %>%
    ungroup() %>%
    mutate(omega = combined_trend * this_params$beta_s) %>%
    group_by(pf_rep) %>%
    mutate(q = get_q(t = time,
                     params = this_params)) %>%
    mutate(s = get_s(t = time, 
                     params = this_params)) %>%
    mutate(R_e = getReff(S = S, 
                         omega = omega,
                         q = q,
                         s = s,
                         params = this_params)) %>%
    ungroup() %>%
    mutate(location = this_state) %>%
    filter(!is.na(date))
  
  # resample trajectories proportional to log likelihood
  samples <- this_filter %>%
    dplyr::select(pf_rep, loglik) %>%
    distinct() %>%
    mutate(weights = exp(loglik - mean(loglik))) %>%
    mutate(weights = weights / sum(weights))
  
  set.seed(20220711)
  sampled_filters <- this_filter %>%
    left_join(samples, by = c("pf_rep", "loglik")) %>%
    slice_sample(n = nrow(samples), replace = TRUE, weight_by = weights) %>%
    dplyr::select(-weights)

  all_states_filters <- bind_rows(all_states_filters, this_filter)
  all_states_data <- bind_rows(all_states_data, this_data)
}



# Plot GA fits and MASE ---------------------------------------------------

example_state <- "Maryland"

all_states_data %>%
  filter(location == example_state) %>%
  dplyr::select(date, cases, deaths, time) %>%
  pivot_longer(cols = c(cases, deaths)) %>%
  filter(time > 0) -> data_for_plot

all_states_filters %>%
  filter(location == example_state) %>%
  dplyr::select(time, pf_rep, C_new, D_new) %>%
  pivot_longer(cols = c(C_new, D_new)) %>%
  group_by(time, name) %>%
  summarise(median_value = median(value),
            lower_95_value = quantile(value, 0.025),
            upper_95_value = quantile(value, 0.975),
            .groups = "drop") %>%
  mutate(name = case_when(
    name == "C_new" ~ "cases",
    TRUE ~ "deaths"
  )) %>%
  left_join(data_for_plot, by = c("time", "name")) -> data_states_for_plot

my_color <- "black"
my_alpha <- 0.2
# cases
data_states_for_plot %>%
  filter(name == "cases") %>%
  filter(value < 20000) %>%
  ggplot(., aes(x = date)) +
  geom_point(aes(y = value), size = 1) +
  geom_ribbon(aes(ymin = lower_95_value, ymax = upper_95_value), 
              fill = my_color, alpha = my_alpha) +
  geom_line(aes(y = median_value), size = 0.75, color = my_color) +
  labs(x = NULL, y = "Daily cases", tag = "A") +
  scale_x_date(date_breaks = "1 months", date_labels = "%b") -> cases_plot

# deaths
data_states_for_plot %>%
  filter(name == "deaths") %>%
  filter(value < 400) %>%
  ggplot(., aes(x = date)) +
  geom_point(aes(y = value), size = 1) +
  geom_ribbon(aes(ymin = lower_95_value, ymax = upper_95_value), 
              fill = my_color, alpha = my_alpha) +
  geom_line(aes(y = median_value), size = 0.75, color = my_color) +
  labs(x = NULL, y = "Daily deaths", tag = "B") +
  scale_x_date(date_breaks = "1 months", date_labels = "%b") -> deaths_plot

# trends
all_states_filters %>%
  filter(location == example_state) %>%
  dplyr::select(date, pf_rep, mobility, latent_trend) %>%
  pivot_longer(cols = c(mobility, latent_trend)) %>%
  group_by(date, name) %>%
  summarise(median_value = mean(value, na.rm = TRUE),
            lower_95_value = quantile(value, 0.025, na.rm = TRUE),
            upper_95_value = quantile(value, 0.975, na.rm = TRUE), 
            .groups = "drop") %>%
  ggplot(., aes(x = date, color = name, fill = name)) +
  geom_ribbon(aes(ymin = lower_95_value, ymax = upper_95_value), 
              alpha = my_alpha, color = NA) +
  geom_line(aes(y = median_value)) +
  labs(x = NULL, y = "Relative transmission reduction", tag = "C") +
  scale_x_date(date_breaks = "1 months", date_labels = "%b") +
  scale_color_brewer(palette = "Set1", name = NULL) +
  scale_fill_brewer(palette = "Set1", name = NULL) +
  theme(legend.position = c(0.15,0.2)) -> trend_plot

# effective R
all_states_filters %>%
  filter(location == example_state) %>%
  dplyr::select(date, pf_rep, R_e) %>%
  group_by(date) %>%
  summarise(median_value = median(R_e, na.rm = TRUE),
            lower_95_value = quantile(R_e, 0.025, na.rm = TRUE),
            upper_95_value = quantile(R_e, 0.975, na.rm = TRUE), 
            .groups = "drop") %>%
  ggplot(., aes(x = date)) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_ribbon(aes(ymin = lower_95_value, ymax = upper_95_value), 
              fill = my_color, alpha = my_alpha, color = NA) +
  geom_line(aes(y = median_value), color = my_color) +
  labs(x = NULL, y = "Reproduction number", tag = "D") +
  scale_x_date(date_breaks = "1 months", date_labels = "%b") -> re_plot

cowplot::plot_grid(cases_plot,
                   deaths_plot,
                   trend_plot,
                   re_plot,
                   ncol = 2,
                   align = "v") -> state_example
ggsave(
  "../figures/example-fits.png",
  state_example,
  height = 4,
  width = 8.5,
  units = "in",
  dpi = 360,
  scale = 1.5
)



# Plot mandates and mobility ----------------------------------------------

mobility <- tibble()
mob_files <- list.files("../data/ucmobility/state_breakdown/")
for(i in 1:length(mob_files)) {
  this_file <- mob_files[i]
  this_loc_abb <- str_split(this_file, "[_]", simplify = TRUE)[2]
  this_loc_abb <- gsub(".rds", "", this_loc_abb)
  if(this_loc_abb == "DC") {
    this_loc <- "DistrictOfColumbia"
  } else {
    this_loc <- state.name[state.abb == this_loc_abb]
  }
  this_mob <- readRDS(paste0("../data/ucmobility/state_breakdown/", this_file))
  this_mob <- this_mob %>%
    mutate(location = this_loc,
           date = Date) %>%
    dplyr::select(location, date, rel_beta_change)
  mobility <- bind_rows(mobility, this_mob)
}

mandates <- read_csv("../data/ihme-mandates-data.csv") %>%
  # why oh why IHME...get your dates right
  separate(Date, into = c("month", "day", "year")) %>%
  mutate(month = str_pad(month, 2, side = "left", pad = "0"),
         day = str_pad(day, 2, side = "left", pad = "0")) %>%
  mutate(Date = paste(year, month, day, sep = "-")) %>%
  mutate(Date = as_date(Date)) %>%
  dplyr::select(Date, Number_of_Mandates_On, State_Name) %>%
  # well, actually "State_Name" is state abbreviation %>%
  rename("State_Abbv" = State_Name) %>%
  # now add in State_Name
  left_join(tibble(
    State_Name = state.name,
    State_Abbv = state.abb,
  ), by = "State_Abbv")

# fill in dates
mandates <- mandates %>%
  group_by(State_Abbv, State_Name) %>%
  arrange(Date) %>%
  complete(Date = seq.Date(min(Date), max(Date), by="day")) %>%
  ungroup()

# pivot for merge in with the other trends
mandates <- mandates %>%
  dplyr::select(State_Name, Date, Number_of_Mandates_On) %>%
  rename("location" = State_Name,
         "date" = Date,
         "num_mandates" = Number_of_Mandates_On)

mobility_with_mandates <- mobility %>%
  left_join(mandates, by = c("location", "date")) %>%
  group_by(location) %>%
  mutate(relative_mandates = num_mandates/max(num_mandates, na.rm = TRUE)) %>%
  dplyr::select(-num_mandates) %>%
  pivot_longer(cols = c(rel_beta_change, relative_mandates)) %>%
  drop_na() %>%
  mutate(value = case_when(
    value > 1 ~ 1,
    TRUE ~ value
  ))

means <- mobility_with_mandates %>%
  group_by(date, name) %>%
  summarise(mean_value = mean(value),
            .groups = "drop")

ggplot(data = mobility_with_mandates, aes(x = date)) +
  geom_line(aes(y = value,
                color = name,
                linetype = location),
            alpha = 0.1) +
  geom_line(data = means, aes(y = mean_value, color = name), size = 1) +
  scale_linetype_manual(values = rep(1, 51)) +
  guides(linetype = "none") +
  geom_vline(aes(xintercept = as.Date("2020-04-24")),
             linetype = 2,
             color = "grey35") +
  annotate(
    geom = "text",
    x = as.Date("2020-05-15"),
    y = 0.1,
    label = "2020-04-24",
    color = "grey35",
    size = 3
  ) +
  scale_color_brewer(
    palette = "Set1",
    name = NULL,
    labels = c("mobility", "mandates")
  ) +
  labs(x = NULL, y = "Relative value", tag = "A") +
  scale_x_date(date_breaks = "1 months", date_labels = "%b") +
  theme(legend.position = c(0.85,0.2)) -> mob_mands_plot

all_states_filters %>%
  dplyr::select(location, pf_rep, date, R_e) %>%
  group_by(location, date) %>%
  summarise(med_re = median(R_e),
            .groups = "drop") %>%
  ggplot(., aes(x = date, y = med_re, group = location)) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_line(alpha = 0.3) +
  geom_vline(aes(xintercept = as.Date("2020-04-24")),
             linetype = 2,
             color = "grey35") +
  annotate(
    geom = "text",
    x = as.Date("2020-05-15"),
    y = 0.1,
    label = "2020-04-24",
    color = "grey35",
    size = 3
  ) +
  labs(x = NULL, y = "Reproduction number", tag = "B") +
  scale_x_date(date_breaks = "1 months", date_labels = "%b") -> re_plot_all

mob_mands_re <- cowplot::plot_grid(mob_mands_plot, re_plot_all)

ggsave(
  "../figures/mobility-mandates-re.png",
  mob_mands_re,
  width = 10,
  height = 3,
  units = "in",
  dpi = 360,
  scale = 1.5
)



# Correlation analysis ----------------------------------------------------

## PLACEHOLDER: READ IN ARCHIVED RESULTS
all_files <- list.files(path = "../../COVID-stochastic-fitting/output/current/", pattern = ".csv")
param_files <- list.files(path = "../../COVID-stochastic-fitting/output/current/", pattern = "params-natural.rds")
state_summaries <- tibble()
state_parameters <- tibble()
state_logliks <- tibble()
statedf <-readRDS("../../COVID-stochastic-fitting/output/current/statedf.rds")
statevec <- gsub(".csv","",all_files)
allstates_pop <- statedf %>% filter(state_full %in% statevec) %>% pull(total_pop) %>% sum()
for(i in 1:length(all_files)) {
  do_file <- all_files[i]
  location <- sub(".csv", "", do_file)
  state_metadata <- statedf %>% filter(state_full == sub(".csv", "", do_file))
  state_pop <- state_metadata %>% pull(total_pop)
  state_initR0 <- state_metadata %>% pull(initR0)
  state_beta_s <- (state_initR0*.1)
  
  tmpparamfile <- paste0("../../COVID-stochastic-fitting/output/current/", param_files[i])
  # tmp state params
  statepars <- readRDS(tmpparamfile) 
  statepars_fixed <- statepars %>% 
    rownames_to_column(var = "param") %>% 
    dplyr::filter(is_fitted == "no") %>% 
    select(param,X1) %>% 
    pivot_wider(values_from = X1, names_from = param) %>% 
    select(-c(MIF_ID, LogLik, LogLik_SE))
  statepars_fitted <- statepars %>% 
    rownames_to_column(var = "param") %>% 
    dplyr::filter(is_fitted == "yes") %>% 
    select(param,X1) %>% 
    pivot_wider(values_from = X1, names_from = param)
  
  statepars_allmle <- bind_cols(statepars_fixed,statepars_fitted)
  
  # results
  tmpfile <- paste0("../../COVID-stochastic-fitting/output/current/", do_file)
  tmp <- read.csv(tmpfile) %>% mutate(date = as.Date(date))
  firstcasedate <- tmp$date %>% min()
  tmp <- tmp %>%
    filter(sim_type == "status_quo" | is.na(sim_type),
           variable %in% c("daily_cases", "daily_deaths", "daily_all_infections", 
                           "actual_daily_cases", "actual_daily_deaths",
                           "mobility_trend", "latent_trend", "combined_trend",
                           "cumulative_all_infections", "cumulative_deaths")) %>%
    dplyr::select(location, sim_type, period, date, variable, mean_value) %>% 
    pivot_wider(names_from = variable, values_from = mean_value) %>%
    # calculate prevalence
    mutate(prevalence = daily_all_infections / (state_pop-cumulative_deaths)) %>% 
    # calculate omega
    mutate(omega = combined_trend * statepars_allmle$beta_s) %>% 
    # calculate mean S
    mutate(S = state_pop - cumulative_all_infections) %>% 
    mutate(susceptible_fraction = S / (state_pop - cumulative_deaths)) %>% 
    # calculate q
    mutate(q = get_q(t = as.numeric(date-firstcasedate), params = statepars_allmle)) %>%
    # calculate s
    mutate(s = get_s(t = as.numeric(date-firstcasedate), params = statepars_allmle)) %>%
    # calculate mean R_e
    mutate(R_e = getReff(S = S, 
                         omega = omega,
                         q = q,
                         s = s,
                         params = statepars_allmle)) %>%
    pivot_longer(cols = !c(location, sim_type, period, date), 
                 names_to = "variable", 
                 values_to = "mean_value",
                 values_drop_na = TRUE)
  state_summaries <- bind_rows(state_summaries, tmp)
}

# reorder matrix
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
  return(cormat)
}
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
# heatmap function
make_heatmap <- function(cm){
  ggheatmap <- cm %>% 
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    scale_y_discrete(position = "right") +
    theme_minimal() + # minimal theme
    coord_fixed() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, size = 5, hjust = 1),
          axis.text.y = element_text(angle = 0, vjust = .5, size = 5, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          legend.justification = c(1, 0),
          legend.position = c(0.6, 0.7),
          legend.direction = "horizontal") +
    guides(fill = guide_colorbar(
      barwidth = 7, barheight = 1,title.position = "top", title.hjust = 0.5))
  # geom_text(aes(Var2, Var1, label = value), color = "black", size = 4)
  return(ggheatmap)
}
GetCorr <- function(df) {
  corrmatrix <- cor(df, use = "complete.obs")
  # Reorder correlation matrix
  corrmatrix.reordered <- reorder_cormat(corrmatrix)
  
  # get triangles
  corrmatrix.lower <- get_lower_tri(corrmatrix)
  corrmatrix.upper <- get_upper_tri(corrmatrix)
  corrmatrix.reordered.lower <- get_lower_tri(corrmatrix.reordered)
  corrmatrix.reordered.upper <- get_upper_tri(corrmatrix.reordered)
  
  # melt
  library(reshape2)
  corrmatrix.upper.melt <- reshape2::melt(corrmatrix.upper, na.rm = TRUE)
  corrmatrix.reordered.upper.melt <- reshape2::melt(corrmatrix.reordered.upper, na.rm = TRUE)
  
  return(make_heatmap(corrmatrix.reordered.upper.melt))
}
# START HERE
myvars <- c("mobility_trend", "latent_trend", "R_e")
all_trends <- state_summaries %>%
  filter(variable %in% myvars) %>%
  dplyr::select(location, date, variable, mean_value) %>%
  filter(date <= as_date("2020-12-31"))

trends_for_cor <- all_trends %>%
  filter(!location %in% c("California", "Connecticut", "Illinois", "Massachusetts"))
all_re_wide <- trends_for_cor %>%
  filter(variable == "R_e") %>%
  pivot_wider(names_from = variable, values_from = mean_value) %>%
  select(date, location, R_e) %>% 
  pivot_wider(names_from = location, values_from = R_e) 
re_p1_df <- all_re_wide %>%
  filter(date < as.Date("2020-04-24")) %>%
  dplyr::select(-date)
re_p2_df <- all_re_wide %>%
  filter(date >= as.Date("2020-04-24")) %>%
  dplyr::select(-date)
all_mo_wide <- trends_for_cor %>%
  filter(variable == "mobility_trend") %>%
  pivot_wider(names_from = variable, values_from = mean_value) %>%
  select(date, location, mobility_trend) %>% 
  pivot_wider(names_from = location, values_from = mobility_trend) 
mo_p1_df <- all_mo_wide %>%
  filter(date < as.Date("2020-04-24")) %>%
  dplyr::select(-date)
mo_p2_df <- all_mo_wide %>%
  filter(date >= as.Date("2020-04-24")) %>%
  dplyr::select(-date)
p1 <- GetCorr(re_p1_df) + ggtitle("Re, Phase 1")
p2 <- GetCorr(re_p2_df) + ggtitle("Re, Phase 2")
p3 <- GetCorr(mo_p1_df) + ggtitle("Mobility, Phase 1")
p4 <- GetCorr(mo_p2_df) + ggtitle("Mobility, Phase 2")
cowplot::plot_grid(p3, p4, p1, p2, ncol = 2, labels = c("A","C","B","D"))
ggsave("../figures/cross-corrs.png", height = 8, width = 10, units = "in", dpi = 360)



# Analysis of variance ----------------------------------------------------

trends <- all_trends %>%
  pivot_wider(names_from = variable, values_from = mean_value) %>%
  mutate(TimePeriod = case_when(
    date < "2020-04-24" ~ "a",
    TRUE ~ "b"))
results <- tibble()
for (state in unique(trends$location)) {
  tmp <- trends %>%
    filter(location == state)
  for (tp in unique(tmp$TimePeriod)) {
    df <- tmp %>%
      filter(TimePeriod == tp)
    m1 <- broom::tidy(anova(lm(R_e ~ 1, data = df))) %>%
      mutate(term = "null")
    m2 <- broom::tidy(anova(lm(R_e ~ mobility_trend, data = df)))
    m3 <- broom::tidy(anova(lm(R_e ~ mobility_trend + latent_trend, data = df)))
    res <- bind_rows(m1, m2, m3) %>%
      filter(term != "Residuals") %>%
      mutate(TimePeriod = tp) %>%
      mutate(location = state)
    results <- bind_rows(results, res)
  }
}
output <- results %>%
  group_by(location, TimePeriod) %>%
  mutate(PropVar = sumsq/max(sumsq)) %>%
  dplyr::select(location, TimePeriod, term, PropVar) %>% 
  distinct() %>%
  filter(term != "null") %>%
  mutate(Label = case_when(
    TimePeriod == "a" ~ "Before May 1, 2020",
    TRUE ~ "After May 1, 2020"
  )) %>%
  mutate(Label = factor(Label, levels = c("Before May 1, 2020", "After May 1, 2020"), labels = c("Ph. 1", "Ph. 2"))) 

output <- output %>%
  left_join(data.frame(
    location = state.name,
    state = state.abb
  ))
ggplot(output, aes(x = Label, y = PropVar)) +
  # geom_col(aes(fill = term), position = position_dodge(), width = 0.5) +
  geom_point(aes(color = term), position = position_dodge(width = 0.5), size= 2) + 
  geom_linerange(aes(color = term, x=Label, ymin=0, ymax=PropVar),
               position = position_dodge(width = 0.5)) +
  facet_geo(~ state, grid = "us_state_grid2") +
  labs(x = NULL, y = NULL) +
  scale_color_brewer(palette = "Set1", name = NULL, labels = c("Latent trend", "Relative mobility"))  +
  scale_y_continuous(breaks = c(0, 1)) +
  ggthemes::theme_few() +
  theme(legend.position = c(0.5, 0.95))
ggsave("../figures/anova_re_states.png", width = 10, height = 7, units = "in", dpi = 360)



# Re summaries ------------------------------------------------------------

re_summaries <- all_states_filters %>%
  dplyr::select(location, date, pf_rep, R_e) %>%
  group_by(location, date) %>%
  summarise(mean_r = mean(R_e),
            median_r = median(R_e),
            lwr95_r = quantile(R_e, 0.025),
            upr95_r = quantile(R_e, 0.975),
            min_r = min(R_e),
            max_r = max(R_e),
            .groups = "drop") %>%
  mutate(phase = case_when(
    date < "2020-04-24" ~ "a",
    TRUE ~ "b"))

# highest in each state
re_summaries %>%
  group_by(location, phase) %>%
  filter(mean_r == max(mean_r)) %>%
  arrange(phase, mean_r)

# lowest in each state
re_summaries %>%
  group_by(location, phase) %>%
  filter(mean_r == min(mean_r)) %>%
  arrange(phase, mean_r)

# average by phase across all states
all_states_filters %>%
  dplyr::select(location, date, pf_rep, R_e) %>%
  mutate(phase = case_when(date < "2020-04-24" ~ "a",
                           TRUE ~ "b")) %>%
  group_by(phase) %>%
  summarise(
    mean_r = mean(R_e),
    median_r = median(R_e),
    lwr95_r = quantile(R_e, 0.025),
    upr95_r = quantile(R_e, 0.975),
    min_r = min(R_e),
    max_r = max(R_e),
    .groups = "drop"
  )
