# ==============================================================================
# Allele and Genotype Distribution Visualizations in R
# Description: Generates publication-ready jitter plots, donut charts, and heatmaps 
#               to explore allele/genotype frequencies by Geographic & Genetic Groups.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Packages
# ------------------------------------------------------------------------------
required_packages <- c("tidyverse", "scales", "forcats")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(tidyverse)
library(scales)
library(forcats)

# ------------------------------------------------------------------------------
# 2. Data Preparation & Constants
# ------------------------------------------------------------------------------
# Path to input dataset
file_path <- "data/allele_genotype_data.xlsx"

# Import dataset (assumes 'd' is loaded from readxl or CSV)
# d <- readxl::read_excel(file_path)

# Mappings & Custom Ordering
geo_translation_map <- c(
  "Nom de la région1"       = "Region Name1",
  "Nom de la région2"       = "Region Name2",
  "Nom de la région3"       = "Region Name3",
  "Nom de la région4"       = "Region Name4",
  "Nom de la région5"       = "Region Name5",
  "Nom de la région6"       = "Region Name6"
)

custom_order_geo <- c("Region Name1", "Region Name2", "Region Name3", 
                      "Region Name4", "Region Name5", "Region Name6")

genetic_name_map <- c(
  "1"     = "Region Name1",
  "2"     = "Region Name2",
  "4"     = "Region Name3",
  "3"     = "Region Name4",
  "5" = "Region Name5",
  "NA"    = "NA"
)

custom_order_genetic <- c("Region Name1", "Region Name2", "Region Name3", "Region Name5", "Region Name5", "NA")

# ------------------------------------------------------------------------------
# 3. Jitter Plot: Alleles by Geographic Group (with Sample Counts)
# ------------------------------------------------------------------------------
# Calculate sample counts per geographic group
sample_counts_geo <- d %>%
  filter(!is.na(`Groupe géographique`)) %>%
  group_by(`Groupe géographique`) %>%
  summarise(n_ind = n(), .groups = 'drop') %>%
  mutate(translated_label = paste0(geo_translation_map[`Groupe géographique`], " (n=", n_ind, ")")) %>%
  deframe()

# Reshape and prepare data
data_alleles_geo_long <- d %>%
  select(`Allèle 1`, `Allèle 2`, `Groupe géographique`) %>%
  pivot_longer(cols = c(`Allèle 1`, `Allèle 2`), names_to = "Allele_Type", values_to = "Allele") %>%
  filter(!is.na(Allele)) %>%
  mutate(
    Geo_Group = factor(sample_counts_geo[`Groupe géographique`], 
                       levels = unname(sample_counts_geo[names(geo_translation_map)])),
    Allele_char = as.character(Allele),
    is_numeric = grepl("\\d", Allele_char),
    num_part = as.numeric(str_extract(Allele_char, "\\d+"))
  ) %>%
  arrange(is_numeric, num_part, Allele_char) %>%
  mutate(Allele = factor(Allele_char, levels = unique(Allele_char)))

p_jitter_geo <- ggplot(data_alleles_geo_long, aes(x = Geo_Group, y = Allele, color = Allele)) +
  geom_point(position = position_jitter(width = 0.3, height = 0.2), alpha = 0.7, size = 2) +
  scale_x_discrete(name = "Geographic Group") +
  labs(title = NULL, y = "Allele") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 16),
    legend.position = "none"
  )

print(p_jitter_geo)

# ------------------------------------------------------------------------------
# 4. Donut Chart: Allele Frequencies by Geographic Group
# ------------------------------------------------------------------------------
allele_prop_geo <- data_alleles_geo_long %>%
  group_by(`Groupe géographique`, Allele) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(`Groupe géographique`) %>%
  mutate(
    Proportion = Count / sum(Count),
    fraction = Proportion,
    percentage = paste0(round(fraction * 100, 1), "%"),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, -1)),
    label_y = ymin + (ymax - ymin) / 2
  )

p_donut_geo <- ggplot(allele_prop_geo, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Allele)) +
  geom_rect(colour = "black") +
  geom_segment(
    data = filter(allele_prop_geo, Proportion >= 0.05),
    aes(x = 4, xend = 4.2, y = label_y, yend = label_y),
    linewidth = 0.6, colour = "black"
  ) +
  geom_text(
    data = filter(allele_prop_geo, Proportion >= 0.05),
    aes(x = 4.3, y = label_y, label = percentage),
    hjust = 0, size = 4.5, fontface = "bold"
  ) +
  coord_polar(theta = "y") +
  facet_wrap(~ `Groupe géographique`, labeller = labeller(`Groupe géographique` = geo_translation_map)) +
  theme_void(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.position = "right",
    strip.text = element_text(size = 18, face = "bold")
  )

print(p_donut_geo)

# ------------------------------------------------------------------------------
# 5. Heatmap: Allele Frequencies by Geographic Group (%)
# ------------------------------------------------------------------------------
heatmap_geo_data <- d %>%
  select(`Allèle 1`, `Allèle 2`, `Groupe géographique`) %>%
  pivot_longer(cols = starts_with("Allèle"), values_to = "Allele") %>%
  filter(!is.na(Allele), !is.na(`Groupe géographique`)) %>%
  group_by(`Groupe géographique`, Allele) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(`Groupe géographique`) %>%
  mutate(freq = round((count / sum(count)) * 100)) %>%
  ungroup() %>%
  mutate(
    Geo_Group = factor(recode(`Groupe géographique`, !!!geo_translation_map), levels = custom_order_geo),
    Allele_char = as.character(Allele),
    num_part = as.numeric(str_extract(Allele_char, "\\d+")),
    Allele = fct_reorder(Allele_char, num_part, .na_rm = FALSE)
  )

p_heatmap_geo <- ggplot(heatmap_geo_data, aes(x = Geo_Group, y = Allele, fill = freq)) +
  geom_tile(color = "gray90") +
  geom_text(aes(label = freq), color = "black", size = 8, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#FFFFCC", "#FFEDA0", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#B10026"),
    name = "Allele Frequency (%)"
  ) +
  labs(x = "Geographic Group", y = "Allele") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14, color = "black"),
    axis.text.y = element_text(face = "bold", size = 14, color = "black"),
    legend.position = "top",
    legend.justification = "right",
    legend.direction = "horizontal",
    panel.grid = element_blank()
  ) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = unit(0.4, "cm")))

print(p_heatmap_geo)

# ------------------------------------------------------------------------------
# 6. Heatmap: Allele Frequencies by Genetic Group (%)
# ------------------------------------------------------------------------------
heatmap_gen_data <- d %>%
  select(`Allèle 1`, `Allèle 2`, `Groupe génétique`) %>%
  pivot_longer(cols = starts_with("Allèle"), values_to = "Allele") %>%
  filter(!is.na(Allele)) %>%
  mutate(group_char = tidyr::replace_na(as.character(`Groupe génétique`), "NA")) %>%
  group_by(group_char, Allele) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(group_char) %>%
  mutate(freq = round((count / sum(count)) * 100)) %>%
  ungroup() %>%
  mutate(
    Genetic_Group = factor(recode(group_char, !!!genetic_name_map), levels = custom_order_genetic),
    Allele_char = as.character(Allele),
    num_part = as.numeric(str_extract(Allele_char, "\\d+")),
    Allele = fct_reorder(Allele_char, num_part, .na_rm = FALSE)
  )

p_heatmap_gen <- ggplot(heatmap_gen_data, aes(x = Genetic_Group, y = Allele, fill = freq)) +
  geom_tile(color = "gray90") +
  geom_text(aes(label = freq), color = "black", size = 8, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#FFFFCC", "#FFEDA0", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#B10026"),
    name = "Allele Frequency (%)"
  ) +
  labs(x = "Genetic Group", y = "Allele") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14, color = "black"),
    axis.text.y = element_text(face = "bold", size = 14, color = "black"),
    legend.position = "top",
    legend.justification = "right",
    legend.direction = "horizontal",
    panel.grid = element_blank()
  ) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barheight = unit(0.4, "cm")))

print(p_heatmap_gen)
