# Script to calculate forecast coverage for each scenario over time.
# This uses the archived results of forecasts "as produced" at the time.
# Thus, the model changes over time. But this is a quantification of
# our real time forecast coverage.
#
# Forecast "coverage" is defined as:
#   percentage of observations that fall within the 95% prediction interval
#   across all states at different forecast horizons.


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(rvest)

potential_dates <- seq.Date(from = as_date("2020-06-09"),
                            to = as_date("2020-12-31"),
                            by = "1 day")
potential_folders <- paste0(
  "https://github.com/CEIDatUGA/COVID-stochastic-fitting/tree/master/output/",
  as.character(potential_dates)
)

potential_folders[1] %>%
  read_html() %>%
  html_nodes(xpath = '//*[@role="rowheader"]') %>%
  html_nodes('span a') %>%
  html_attr('href') %>%
  sub('blob/', '', .) %>%
  paste0('https://raw.githubusercontent.com', .) %>%
  purrr::map_df(read.csv) ->  combined_data

forecast_date <- combined_data %>% 
  filter(period == "Past") %>%
  summarise(max_date = max(date)) %>%
  pull(max_date)

combined_data$forecast_date <- forecast_date

horizon_dates <- c(
  as_date(forecast_date) + days(1),
  as_date(forecast_date) + days(7),
  as_date(forecast_date) + days(14),
  as_date(forecast_date) + days(28)
) 

mark_data <- combined_data %>%
  filter(variable %in% c("daily_cases", "daily_deaths")) %>%
  filter(date %in% as.character(horizon_dates)) %>%
  dplyr::select(location, forecast_date, date, sim_type, 
                period, date, variable, lower_95, upper_95)
