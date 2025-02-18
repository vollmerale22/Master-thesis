library(readxl)
library(tidyr)
library(dplyr)
#library for wrap_plots
library(patchwork)
library(ggplot2)

# Read in the data (assumes your first row is column names):
neg <- read_csv("C:/Users/vollm/Desktop/Thesis - Imports EU/rmse_sel_1_n_60_reg_1_2_3_4_5_6_7_h_-1_2017-01-15_2024-07-15.csv")
now <- read_csv("C:/Users/vollm/Desktop/Thesis - Imports EU/rmse_sel_1_n_60_reg_1_2_3_4_5_6_7_h_0_2017-01-15_2024-07-15.csv")
pos <- read_csv("C:/Users/vollm/Desktop/Thesis - Imports EU/rmse_sel_1_n_60_reg_1_2_3_4_5_6_7_h_1_2017-01-15_2024-07-15.csv")

# Combine the datasets and add a Horizon column
data_combined <- bind_rows(
  neg %>% mutate(Horizon = "t-1"),
  now %>% mutate(Horizon = "t"),
  pos %>% mutate(Horizon = "t+1")
)
# Assign a name to the unnamed column with method names
colnames(data_combined)[1] <- "Method"  # Rename the first column to "Method"

# Reshape the data to long format for ggplot
data_tidy <- data_combined %>%
  filter(Method != "ar") %>%  # Exclude the "ar" row
  pivot_longer(cols = -c(Method, Horizon), names_to = "Scenario", values_to = "RMSE") %>%
  group_by(Horizon) %>%
  mutate(Standardised_RMSE = RMSE / RMSE[Method == "pred_ols"]) %>%  # Standardise relative to OLS
  ungroup()

# Define method groups for grouping in the plot
data_tidy <- data_tidy %>%
  mutate(Group = case_when(
    Method %in% c("pred_ms", "pred_qr") ~ "Trad.",
    Method %in% c("pred_rf", "pred_xgbt") ~ "ML tree",
    Method %in% c("pred_mrf", "pred_xgbl") ~ "ML reg."
  ),
  Method = factor(Method, levels = c("pred_ms", "pred_qr", "pred_rf", "pred_xgbt", "pred_mrf", "pred_xgbl")),
  Group = factor(Group, levels = c("Trad.", "ML tree", "ML reg."))
  )

# Define custom colours for each method
custom_colours <- c(
  "pred_ms" = "darkgrey", "pred_qr" = "lightgrey",
  "pred_rf" = "steelblue", "pred_xgbt" = "lightblue",
  "pred_mrf" = "darkred", "pred_xgbl" = "salmon"
)

# Define labels for the legend
custom_labels <- c(
  "pred_ms" = "Markov-switching",
  "pred_qr" = "Quantile reg.",
  "pred_rf" = "Random forest",
  "pred_xgbt" = "Gradient boosting",
  "pred_mrf" = "Macroeconomic random forest",
  "pred_xgbl" = "Gradient linear boosting"
)
data_tidy_filtered <- data_tidy %>% filter(Method != "pred_ols")
data_tidy_filtered <- data_tidy_filtered %>%
  mutate(Horizon = factor(Horizon, levels = c("t-1","t",  "t+1"))) 


# Filter data for each panel
panel_a <- data_tidy_filtered %>% filter(Scenario == "total")
panel_b <- data_tidy_filtered %>% filter(Scenario == "crisis")
panel_c <- data_tidy_filtered %>% filter(Scenario == "normal")

# Function to create each panel
create_panel <- function(data, title) {
  ggplot(data, aes(x = Group, y = Standardised_RMSE, fill = Method)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.7) +
    facet_wrap(~ Horizon, ncol = 3) +
    scale_fill_manual(values = custom_colours, labels = custom_labels) +
    labs(title = title, x = NULL, y = NULL, fill = NULL) +
    geom_hline(yintercept = 1.0, colour = "black", size = 0.6) +  # Add line at 1.0
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_line(colour = "grey85", size = 0.4),  # Light major grid
      panel.grid.minor = element_line(colour = "grey90", size = 0.2),  # Light minor grid
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      legend.key.size = unit(0.5, "cm")
    )
}

# Generate each panel
panel_a_plot <- create_panel(panel_a, "Panel A. Full sample (Jan. 2017 - Jul. 2024)")
panel_b_plot <- create_panel(panel_b, "Panel B. Crisis sample (Jan. 2020 - Apr. 2022)")
panel_c_plot <- create_panel(panel_c, "Panel C. Normal sample (Jan. 2017 - Dec. 2019)")

# Combine all panels using wrap_plots
final_plot <- wrap_plots(
  panel_a_plot + theme(legend.position = "none"),  # Remove legend from Panel A
  panel_b_plot + theme(legend.position = "none"),  # Remove legend from Panel B
  panel_c_plot + theme(legend.position = "none"),  # Remove legend from Panel C
  ncol = 1
) + theme(legend.position = "bottom")  # Add legend at the bottom
final_plot

