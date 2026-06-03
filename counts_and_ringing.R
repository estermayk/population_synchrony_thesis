library(lme4)
library(insight)
library(dplyr)
library(ggplot2)
library(tidyverse)

adults <- read.csv("data/Adults.csv")
nestlings <- read.csv("data/Nestlings.csv")

#only looking at spring captures
adults <- adults[adults$season == "spring",]

length(adults$ring[duplicated(adults$ring)])

nestlingrings <- c(nestlings$ring)

?match

match(nestlingrings, adults$ring)

nestlingrings %in% adults$ring

recruits <- adults %>%
  filter(ring %in% nestlings$ring)

recruitsbyyear <- recruits %>%
  group_by(year) %>%
  summarise(countadults = n())

ggplot(recruitsbyyear, aes(x = factor(year), y = count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Recruits",
       x = "Year",
       y = "Recruits")

adultsbyyear <- adults %>% 
  group_by(year) %>%
  summarise(countadults = n())

adultsbyyear <- adultsbyyear %>%
  mutate(year = as.character(year))

recruitsbyyear <- recruitsbyyear %>%
  mutate(year = as.character(year))

recruitsandadults <- left_join(adultsbyyear, recruitsbyyear, by = "year")

recruitsplot <- recruitsandadults %>%
  pivot_longer(cols = c(count.x, count.y), 
               names_to = "Category", values_to = "Count") 

ggplot(recruitsplot, aes(x = factor(year), y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("skyblue", "lightcoral"), labels = c("Adults", "Recaptured Nestlings")) +
  theme_minimal() +
  labs(title = "Annual Counts of Adults and Recaptured Nestlings",
       x = "Year",
       y = "Count",
       fill = "Category")
