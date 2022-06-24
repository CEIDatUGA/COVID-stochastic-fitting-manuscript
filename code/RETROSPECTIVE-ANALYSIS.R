# Retrospective analysis of model dynamics


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(lubridate)
library(pomp)
library(ggthemes)
library(cowplot)

ggplot2::theme_set(ggthemes::theme_few(base_size = 12))
source("../code/result-exploration/getReff.R")



# Extract latent trends and compute Re ------------------------------------

all_files <- list.files("../output")
all_mif_files <- all_files[grep("mif_res_2.RDS", all_files)]
all_filter_files <- all_files[grep("filtered_states.RDS", all_files)]
all_pomp_files <- all_files[grep("pomp_model.RDS", all_files)]

all_states_filters <- tibble()
 
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
    pivot_longer(cols = X1:X10) %>%
    dplyr::select(-is_fitted) %>%
    pivot_wider(names_from = parameter, values_from = value) %>%
    filter(LogLik == min(LogLik))
  
  this_filter <- this_filter %>%
    left_join(this_data %>% dplyr::select(date, time), by = "time") %>%
    mutate(latent_trend = exp(trendO) / (1+exp(trendO))) %>%
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
  
  all_states_filters <- bind_rows(all_states_filters, this_filter)
}


# Plot GA fits and MASE ---------------------------------------------------

ga_pomp_list <- readRDS("../output/Georgia-pomp_model.RDS")
ga_data <- ga_pomp_list$pomp_data
ga_states <- readRDS("../output/Georgia-filtered_states.RDS")

ga_data %>%
  dplyr::select(date, cases, deaths, time) %>%
  pivot_longer(cols = c(cases, deaths)) %>%
  filter(time > 0) -> data_for_plot

ga_states %>%
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
  filter(location == "Georgia") %>%
  dplyr::select(date, pf_rep, mobility, latent_trend) %>%
  pivot_longer(cols = c(mobility, latent_trend)) %>%
  group_by(date, name) %>%
  summarise(median_value = median(value),
            lower_95_value = quantile(value, 0.025),
            upper_95_value = quantile(value, 0.975), 
            .groups = "drop") %>%
  ggplot(., aes(x = date, color = name, fill = name)) +
  geom_ribbon(aes(ymin = lower_95_value, ymax = upper_95_value), 
              alpha = my_alpha) +
  geom_line(aes(y = median_value), size = 0.75) +
  labs(x = NULL, y = "Relative transmission reduction", tag = "C") +
  scale_x_date(date_breaks = "1 months", date_labels = "%b") +
  scale_color_brewer(palette = "Set1", name = NULL) +
  scale_fill_brewer(palette = "Set1", name = NULL) +
  theme(legend.position = c(0.15,0.2)) -> trend_plot

# effective R
all_states_filters %>%
  filter(location == "Georgia") %>%
  dplyr::select(date, pf_rep, R_e) %>%
  group_by(date) %>%
  summarise(median_value = median(R_e),
            lower_95_value = quantile(R_e, 0.025),
            upper_95_value = quantile(R_e, 0.975), 
            .groups = "drop") %>%
  ggplot(., aes(x = date)) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_ribbon(aes(ymin = lower_95_value, ymax = upper_95_value), 
              fill = my_color, alpha = my_alpha) +
  geom_line(aes(y = median_value), size = 0.75, color = my_color) +
  labs(x = NULL, y = "Reproduction number", tag = "D") +
  scale_x_date(date_breaks = "1 months", date_labels = "%b") -> re_plot

cowplot::plot_grid(cases_plot,
                   deaths_plot,
                   trend_plot,
                   re_plot,
                   ncol = 2,
                   align = "v") -> ga_example
ggsave(
  "../figures/ga-example-fits.png",
  ga_example,
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
