pacman::p_load(tidytable, lubridate, ggplot2, lme4, fixest, scales, slider)
masks20 <- arrow::read_parquet("masks20.parquet")
masks21 <- arrow::read_parquet("masks21.parquet")
# masks22 <- fread("data_download_file_best_masks_2022.csv")
owd <- arrow::read_parquet("owd.parquet") # all data owd

brazil <- bind_rows(masks20, masks21) |> filter(location_name == "Brazil")
brazil_owd <- owd |> filter(location=="Brazil")
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
  ) |> 
  group_by(location_id) |> 
  mutate(
    mask_7day = slide_index_dbl(
      mask_use_mean,                       # calculate on new_cases
      .i = date,       # indexed with date_onset 
      .f = ~mean(.x, na.rm = T),     # function is sum() with missing values removed
      .before = days(7),
      .complete = T)               # window is the DAY and 6 prior DAYS
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
    new_cases_lead_1 = lead(new_cases_per_million, n = 1),
    stringency_index_lead_2 = lead(stringency_index, n = 2),
    stringency_index_lag_2 = lag(stringency_index, n = 2),
    mask_7day_percent = mask_7day*100,
    excess_lead_1 = lead(excess_mortality),
    excess_lead_2 = lead(excess_mortality, n = 2),
    excess_lead_4 = lead(excess_mortality, n = 4),
    people_fully_vaccinated_per_hundred = replace_na(people_fully_vaccinated_per_hundred, 0) # 2020 didnt had any fully vaccinated
  )


dt21 <- dt |>
  filter(date < "2022-01-06") |>
  filter(!is.na(excess_mortality_cumulative_per_million)) |>
  mutate(
    date2 = as_date(date),
    date3 = floor_date(date2, "month")
  ) 

mod1 <- feols(excess_lead_2 ~ mask_7day_percent | location + date, data = dt21)
mod1
mod2 <- feols(cum_excess_lead_2 ~ mask_7day_percent | location + date, data = dt21)
mod2

mod3 <- feols(mask_7day_percent ~ stringency_index_lag_2 | location + date, data = dt21)


etable(mod1, mod2,mod3,
       coefstat = "confint",
  vcov = "twoway"
)


dt_plot <- dt21 |> 
  filter(date>="2020-03-01") |> 
  group_by(location) |> 
  mutate(mean_all = mean(mask_7day,na.rm=T),
         mean_50 = case_when(
           mean_all>0.5 ~ "high",
           mean_all>0.35 ~ "moderate",
           mean_all<=0.35 ~ "low")) |>
  ungroup()
  
p1 <- dt_plot |> 
  filter(mean_50 == "high") |> 
  ggplot(aes(date, y = mask_7day, group = location))+
  geom_line(alpha=0.2, linewidth=1)+
  geom_hline(aes(yintercept = 0.5), linetype="dashed")+
  scale_y_continuous(labels=scales::percent,
                     limits = c(0,1))+
  scale_x_date(breaks = c(as.Date("2020-06-01"),
                          as.Date("2021-01-01"),
                          as.Date("2021-06-01"),
                          as.Date("2022-01-01")),
               date_labels = "%b %d\n%Y")+
  labs(y= "Mask Usage - 7 days moving average",
       subtitle="Period average > 50%",
       x="")+
  theme_minimal()

p2 <- dt_plot |> 
  filter(mean_50 == "moderate") |> 
  ggplot(aes(date, y = mask_7day, group = location))+
  geom_line(alpha=0.2, linewidth=1)+
  geom_hline(aes(yintercept = 0.5), linetype="dashed")+
  scale_y_continuous(labels=scales::percent,
                     limits = c(0,1))+
  scale_x_date(breaks = c(as.Date("2020-06-01"),
                          as.Date("2021-01-01"),
                          as.Date("2021-06-01"),
                          as.Date("2022-01-01")),
               date_labels = "%b %d\n%Y")+
  labs(y= "",
       subtitle="Period average 36-50%",
       x="")+
  theme_minimal()

p3 <- dt_plot |> 
  filter(mean_50 == "low") |> 
  ggplot(aes(date, y = mask_7day, group = location))+
  geom_line(alpha=0.2, linewidth=1)+
  geom_hline(aes(yintercept = 0.5), linetype="dashed")+
  scale_y_continuous(labels=scales::percent,
                     limits = c(0,1))+
  scale_x_date(breaks = c(as.Date("2020-06-01"),
                          as.Date("2021-01-01"),
                          as.Date("2021-06-01"),
                          as.Date("2022-01-01")),
               date_labels = "%b %d\n%Y")+
  labs(y= "",
       subtitle="Period average ≤ 35%",
       x="")+
  theme_minimal()



geom_dotplot(aes(group=date,
                   y=mask_7day), 
               alpha=0.5,
               median.linewidth = 0.5,
               box.colour = "gray60",
               whisker.colour = "gray60",
               outliers = F)


p2 <- dt21 |> 
  filter(date>="2020-03-01") |> 
  ggplot(aes(date))+
  geom_boxplot(aes(group=date, y=excess_mortality), alpha=0.5, outliers=T)

p1/p2


# Brazil ####

# select only necessary variables
masks <- brazil |> 
  select(
    location_id,
    date,
    location_name,
    mandates_mean,
    mask_use_obs,
    mask_use_mean
  ) |> 
  group_by(location_id) |> 
  mutate(
    mask_7day = slide_index_dbl(
      mask_use_mean,                       # calculate on new_cases
      .i = date,       # indexed with date_onset 
      .f = ~mean(.x, na.rm = T),     # function is sum() with missing values removed
      .before = days(30),
      .complete = T)               # window is the DAY and 6 prior DAYS
  )
# select variables from owd data
owd1 <- brazil_owd |> 
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
    new_cases_lead_1 = lead(new_cases_per_million, n = 1),
    stringency_index_lead_2 = lead(stringency_index, n = 2),
    stringency_index_lag_2 = lag(stringency_index, n = 2),
    mask_7day_percent = mask_7day*100,
    excess_lead_1 = lead(excess_mortality),
    excess_lead_2 = lead(excess_mortality, n = 2),
    excess_lead_4 = lead(excess_mortality, n = 4),
    people_fully_vaccinated_per_hundred = replace_na(people_fully_vaccinated_per_hundred, 0) # 2020 didnt had any fully vaccinated
  )


dt21 <- dt |>
  filter(date < "2022-01-06") |>
  mutate(new_cases_per_million = na_if(new_cases_per_million,0)) |> 
  filter(!is.na(excess_mortality_cumulative_per_million)) |>
  mutate(
    date2 = as_date(date),
    date3 = floor_date(date2, "month")
  ) 


p4 <- dt21 |>
  filter(date >= "2020-03-01") |>
  ggplot(aes(x = date)) +
  
  # First geom_line for excess_mortality. Color is mapped inside aes().
  geom_line(aes(y = excess_mortality, color = "Excess Mortality"), linewidth = 1) +
  
  # Second geom_line for mask usage. Note we DIVIDE by the coefficient.
  geom_line(aes(y = mask_7day_percent , color = "Mask Usage"), linewidth = 1) +
  
  # Add the secondary axis. We MULTIPLY by the coefficient in the trans formula.
  scale_y_continuous(
    limits = c(0,100),
    # Name of the primary y-axis
    name = "Excess Mortality (%)",
    
    # Add the secondary y-axis
    sec.axis = sec_axis(~ . , name = "Mask Usage (%)")
  ) +
  
  # Define the colors for the lines and the legend title
  scale_color_manual(
    name = "Metric",
    values = c(
      "Excess Mortality" = "#E69F00", # Sky Blue
      "Mask Usage" = "#56B4E9"      # Orange
    )
  ) +
  
  # Formatting the x-axis (your original code was good)
  scale_x_date(
    breaks = as.Date(c("2020-06-01", "2021-01-01", "2021-06-01", "2022-01-01")),
    date_labels = "%b %d\n%Y"
  ) +
  
  # Optional: Match axis title colors to the line colors for clarity
  theme_minimal() +
  theme(
    axis.title.y.right = element_text(color = "#56B4E9"), # Sky Blue
    axis.text.y.right = element_text(color = "#56B4E9"),
    axis.title.y.left = element_text(color = "#E69F00"), # Orange
    axis.text.y.left = element_text(color = "#E69F00"),
    legend.position = "bottom"
  ) +
  labs(x = "Date", title = "Excess Mortality and Mask Usage Over Time-Brazil")
library(patchwork)
(p1 + p2 + p3) / p4+plot_annotation(tag_levels = "A")

dt21 |> 
  filter(date>="2020-03-01") |> 
  ggplot(aes(date))+
  geom_line(aes(y=mask_7day_percent),color="green")+
  geom_line(aes(y = excess_mortality))+
  scale_x_date(breaks = c(as.Date("2020-06-01"),
                          as.Date("2021-01-01"),
                          as.Date("2021-06-01"),
                          as.Date("2022-01-01")),
               date_labels = "%b %d\n%Y")+
  labs(y = "Mask usage - 30 days moving avarege", y ="")+
  theme_minimal()

  geom_dotplot(aes(group=date,
                   y=mask_7day), 
               alpha=0.5,
               median.linewidth = 0.5,
               box.colour = "gray60",
               whisker.colour = "gray60",
               outliers = F)+
  ggthemes::theme_clean()
