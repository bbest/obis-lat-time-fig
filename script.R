require(dplyr)
require(ggplot2)

load("obis.dat")

d <- data %>%
  filter(
    latitude <= 90 | latitude >= -90,
    !is.na(depth)) %>%
  mutate(
    zone = cut(depth, c(0, 20, 200, Inf), labels=c('< 20', '20 - 200', '> 200')), 
    band = round(latitude / 5) * 5) %>%
  group_by(band, zone) %>%
  summarize(
    records = n())

ggplot(d, aes(x = band, y = records, fill = zone)) +
  geom_bar(stat = 'identity') +
  scale_fill_brewer(palette='YlGnBu', name='zone (m)') +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x='latitude')
