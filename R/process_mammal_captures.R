library(tidyverse)

## give more descriptive names to sites
sites <- tribble(~pop, ~SiteName,
                 "CS", "Carlo",
                 "FRN", "Field River North",
                 "FRS", "Field River South",
                 "KSE", "Kunnamuka Swamp East",
                 "MC", "Main Camp",
                 "SS", "South Site",
                 "WS", "Way Site"
)

abund_summ <- abund %>%
  select(-`Notomys alexis`, -`Sminthopsis youngsoni`) %>%
  left_join(sites) %>%
  drop_na(pop) %>%
  ## convert month-year character column to proper dates using lubridate
  separate(MonthYear, c("Month", "year"), "\\.") %>%
  mutate(Month = sapply(Month, function(x) which(month.abb == x))) %>%
  mutate(date = myd(paste(Month, year, sep = "-"), truncated = 1)) %>%
  ## summarise by site and date
  group_by(date, pop) %>%
  summarise(abund = mean(`Pseudomys hermannsburgensis`),
            .groups = "drop") 

write_csv(abund_summ, "data/hermannsburg_abund.csv")