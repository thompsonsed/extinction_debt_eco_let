# Table generation for some conceptual numbers

library(tidyverse)
source("preston.R")

# Generate the data for the summary table in the main text.
# Load the range for dispersal values across real maps (calculated by the algorithm as described)
# then calculate the minimum and maximum proportion of species remaining at each spatial and
# temporal scale, with varying levels of habitat fragmentation and habitat loss.
table_df <- read.csv(file.path("results", "dispersal_ranges.csv")) %>% select(-X) %>%
  mutate(a_max = a_max ^2,
         area = a_max * proportion_cover) %>% 
  group_by(a_max, area, map_file, proportion_cover, number_steps, mean_distance) %>%
  summarise(sigma_e = estimate_sigma_rayleigh(mean(mean_distance), n_steps = mean(number_steps)),
            effective_connectivity_sq = (area/a_max)*sigma_e^2) %>% 
  group_by(proportion_cover, a_max, area, number_steps) %>%
  summarise(min_sigma_e = min(sigma_e),
            max_sigma_e = min(16, max(sigma_e)),
            min_effective_connectivity = min(sqrt(effective_connectivity_sq)),
            max_effective_connectivity = max(sqrt(effective_connectivity_sq))) %>%
  filter(number_steps == 10000) %>%
  mutate(speciation_rate = 1/number_steps,
         contig_richness = S_contig(a_max, 0.00001, 16^2),
         inst_richness_upper = S_random(a_max, area, 0.00001, 16^2),
         inst_richness_lower = S_contig(area, 0.00001, 16^2),
         eq_richness_upper = S_random_equilibrium(a_max, area, 0.00001, max_sigma_e^2),
         eq_richness_lower = S_random_equilibrium(a_max, area, 0.00001, min_sigma_e^2),
         prop_inst_upper = inst_richness_upper/contig_richness,
         prop_inst_lower = inst_richness_lower/contig_richness,
         prop_eq_upper = eq_richness_upper/contig_richness,
         prop_eq_lower = eq_richness_lower/contig_richness)
