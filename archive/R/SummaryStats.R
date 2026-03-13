srdatwna <- read.csv(here::here("DataIn/SRinputfile.csv"))

# starwars %>% count(species)
sr <- srdatwna %>% 
  count(Name, Stream) %>% 
  rename("Num. Years" = n, "Stock" = Name, "LifeHist" = Stream) 
# Add a summary of start year and end year
# e.g. first row and last row of the count?

# kable(sr)

# write.csv(sr, here::here("DataOut/SR_timeseries_summary.csv"))