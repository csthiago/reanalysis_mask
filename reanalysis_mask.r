pacman::p_load(tidytable, lubridate, ggplot2, lme4, fixest, scales)
masks20 <- arrow::read_parquet("masks20.parquet")
masks21 <- arrow::read_parquet("masks21.parquet")
# masks22 <- fread("data_download_file_best_masks_2022.csv")
owd <- arrow::read_parquet("owd.parquet") # all data owd
bs_dt <- data.table::fread("fulldataset.csv")
# get only locations from the paper
masks20 <- masks20 |> filter(location_name %in% unique(bs_dt$location))
masks21 <- masks21 |> filter(location_name %in% unique(bs_dt$location))
# check all 24 appears
length(unique(masks20$location_name))
length(unique(masks21$location_name))

# select only necessary variables
masks <- bind_rows(masks20, masks21) |>
  select(
    location_id,
    date,
    location_name,
    mandates_mean,
    mask_use_obs,
    mask_use_mean
  )
# select variables from owd data
owd1 <- owd |>
  filter(location %in% unique(bs_dt$location)) |>
  select(
    location, date, people_fully_vaccinated_per_hundred, stringency_index,
    new_cases_per_million,
    new_deaths_per_million,
    total_deaths_per_million,
    total_cases_per_million,
    excess_mortality_cumulative_absolute:excess_mortality_cumulative_per_million
  ) |>
  group_by(location) |>
  # carry over people fully vaccinated from the last week
  fill(people_fully_vaccinated_per_hundred, .direction = "down")
# Keep only weeks with excess mortality estimates
dt <- owd1 |>
  filter(!is.na(excess_mortality_cumulative)) |>
  left_join(masks,
    by = c(
      "location" = "location_name",
      "date"
    )
  )

dt <- dt |>
  group_by(location) |>
  mutate(
    cum_excess_lead_1 = lead(excess_mortality_cumulative, n = 1),
    cum_excess_lead_2 = lead(excess_mortality_cumulative, n = 2),
    new_deaths_covid_lag1 = lag(new_deaths_per_million, n = 1),
    new_deaths_covid_1 = lead(new_deaths_per_million, n = 1),
    new_deaths_covid_lag2 = lag(new_deaths_per_million, n = 2),
    total_deaths_covid_2 = lead(total_deaths_per_million, n = 2),
    new_cases_lag2 = lag(new_cases_per_million, n = 2),
    stringency_index_lag2 = lag(stringency_index, n = 2),
    mask_lag = lag(mask_use_mean, n = 1),
    mask_lag_2 = lag(mask_use_mean, n = 2),
    excess_lead_1 = lead(excess_mortality),
    excess_lead_2 = lead(excess_mortality, n = 2),
    people_fully_vaccinated_per_hundred = replace_na(people_fully_vaccinated_per_hundred, 0) # 2020 didnt had any fully vaccinated
  )


dt21 <- dt |>
  filter(date < "2022-01-06") |>
  filter(!is.na(excess_mortality_cumulative_per_million)) |>
  mutate(
    mask_use_mean_sq = mask_use_mean * mask_use_mean,
    mask_100 = mask_use_mean * 100,
    mask_lag_100 = mask_lag * 100,
    date2 = as_date(date),
    date3 = floor_date(date2, "month")
  ) 

mod1 <- feols(stringency_index ~ new_deaths_covid_lag2 | location + date, data = dt21)
mod1
mod2 <- feols(mask_100 ~ stringency_index_lag2 | location + date, data = dt21)
mod2
mod3 <- feols(mask_100 ~ new_deaths_covid_lag2 | location + date, data = dt21)
mod3
mod4 <- feols(cum_excess_lead_2 ~ mask_100 | location + date, data = dt21)
mod4
mod5 <- feols(cum_excess_lead_2 ~ new_deaths_per_million | location + date, data = dt21)
mod5

mod6 <- feols(cum_excess_lead_2 ~ stringency_index | location + date, data = dt21)
mod6

etable(mod1, mod2, mod3, mod4,mod5,
       coefstat = "confint",
  vcov = "twoway"
)

adf_dt <- dt21 |>
  select(location, mask_use_mean, excess_mortality) |>
  drop_na(mask_use_mean) |>
  group_split(location)

adf_ps <- NULL
for (i in 1:24) {
  print(unique(adf_dt[[i]]$location))
  print(tseries::adf.test(adf_dt[[i]]$mask_use_mean))
  temp <- data.frame(
    pais = unique(adf_dt[[i]]$location),
    p_value = tseries::adf.test(adf_dt[[i]]$mask_use_mean)$p.value
  )
  adf_ps <- rbind(temp, adf_ps)
}
train_sec <- function(primary, secondary, na.rm = TRUE) {
  # Thanks Henry Holm for including the na.rm argument!
  from <- range(secondary, na.rm = na.rm)
  to   <- range(primary, na.rm = na.rm)
  # Forward transform for the data
  forward <- function(x) {
    rescale(x, from = from, to = to)
  }
  # Reverse transform for the secondary axis
  reverse <- function(x) {
    rescale(x, from = to, to = from)
  }
  list(fwd = forward, rev = reverse)
}
sec <- train_sec(
  dt21 |>
    filter(location == "Austria") |>
    pull(stringency_index),
  dt21 |>
    filter(location == "Austria") |>
    pull(new_deaths_per_million)
)
p1 <- dt21 |>
  filter(location == "Austria") |>
  ggplot(aes(date)) +
  geom_line(aes(y = mask_100, color = "Mask Usage")) +
  geom_line(aes(y = stringency_index, color = "Stringency Index")) +
  geom_vline(aes(xintercept = as.Date("2021-11-11")), linetype=2)+
  annotate(geom = "text", label = "Omicron",
           x = as.Date("2021-11-11"),
           y = 0.9*max(dt21 |>
                         filter(location == "Austria") |> 
                         pull(max(stringency_index))),
           hjust = 1.1, vjust = 1,size=8/.pt)+
  facet_wrap(~location, scales = "free") +
  scale_x_date( date_labels = "%b\n%Y" )+
  scale_y_continuous(limits = c(0,100))+
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(colour = "", x = "Date", y = "Mask Usage (%)\nStringency Index")

p2 <- dt21 |>
  filter(location == "Austria") |>
  ggplot(aes(date)) +
  
  geom_line(aes(y = new_deaths_per_million, color = "New Deaths per million")) +
  geom_line(aes(y = sec$fwd(stringency_index), color = "Stringency Index")) +
  geom_vline(aes(xintercept = as.Date("2021-11-11")), linetype=2)+
  annotate(geom = "text", label = "Omicron",
           x = as.Date("2021-11-11"),
           y = 0.95*max(dt21 |>
                         filter(location == "Austria") |> 
                         pull(max(new_deaths_per_million))),
           hjust = 1, vjust = 1,size=8/.pt)+
  scale_y_continuous(
    name = "New Deaths per million",
    labels = scales::comma,
    sec.axis = sec_axis(~ sec$rev(.),
      name = "Stringency Index"
    )
  ) +
  scale_x_date( date_labels = "%b\n%Y" )+
  facet_wrap(~location, scales = "free") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(colour = "", x = "Date", y = "")


sec <- train_sec(
  dt21 |>
    filter(location == "United Kingdom") |>
    pull(stringency_index),
  dt21 |>
    filter(location == "United Kingdom") |>
    pull(new_deaths_per_million)
)

p3 <- dt21 |>
  filter(location == "United Kingdom") |>
  ggplot(aes(date)) +
  geom_line(aes(y = mask_100, color = "Mask Usage")) +
  geom_line(aes(y = stringency_index, color = "Stringency Index")) +
  geom_vline(aes(xintercept = as.Date("2021-11-11")), linetype=2)+
  annotate(geom = "text", label = "Omicron",
           x = as.Date("2021-11-11"),
           y = 0.9*max(dt21 |>
                         filter(location == "Austria") |> 
                         pull(max(stringency_index))),
           hjust = 1.1, vjust = 1,size=8/.pt)+
  facet_wrap(~location, scales = "free") +
  scale_x_date( date_labels = "%b\n%Y" )+
  scale_y_continuous(limits = c(0,100))+
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(colour = "", x = "Date", y = "Mask Usage (%)\nStringency Index")


p4 <- dt21 |>
  filter(location == "United Kingdom") |>
  ggplot(aes(date)) +
  geom_line(aes(y = new_deaths_per_million, color = "New Deaths per million")) +
  geom_line(aes(y = sec$fwd(stringency_index), color = "Stringency Index")) +
  geom_vline(aes(xintercept = as.Date("2021-11-11")), linetype=2)+
  annotate(geom = "text", label = "Omicron",
           x = as.Date("2021-11-11"),
           y = 0.9*max(dt21 |>
                         filter(location == "United Kingdom") |> 
                         pull(max(new_deaths_per_million))),
           hjust = 1, vjust = 1,
           size=8/.pt)+
  scale_y_continuous(
    name = "New Deaths per million",
    labels = scales::comma,
    sec.axis = sec_axis(~ sec$rev(.),
      name = "Stringency Index"
    )
  ) +
  scale_x_date( date_labels = "%b\n%Y" )+
  facet_wrap(~location, scales = "free") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(colour = "", x = "Date", y = "")
library(patchwork)
((p1/p2) |(p3/p4)) + plot_annotation(tag_levels = "A")

