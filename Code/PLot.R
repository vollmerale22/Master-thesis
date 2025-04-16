filtered_data <- data
filtered_data$date <- as.Date(filtered_data$date)

# Filter data to start from 2000-01-01
filtered_data <- filtered_data[filtered_data$date >= as.Date("2000-01-01"), ]
#rename the columns without the CPB at the end
colnames(filtered_data) <- gsub("CPB ", "", colnames(filtered_data))


plot_data <- filtered_data %>%
  select(date, 2:10) %>%  # Select date and columns 2-10
  pivot_longer(cols = -date, names_to = "Series", values_to = "Value") %>%
  # Clean the Series names by removing "CPB"
  mutate(Series = gsub("CPB", "", Series))

# Create the plot
ggplot(plot_data, aes(x = date, y = Value, color = Series)) +
  geom_line() +
  #make the lines a bit fatter
  geom_line(size = 0.5) +
  theme_minimal() +
  labs(title = "CPB World Trade Monitor ",
       x = "", 
       y = "Merchandise Trade Volume Index") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y",
               limits = c(as.Date("2000-01-01"), as.Date("2024-12-31")) # Set explicit limits
  ) +
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text = element_text(color = "black"),
        legend.background = element_rect(fill = "white", color = "grey90"))+
  coord_cartesian(expand = FALSE) # Remove padding around the plot
