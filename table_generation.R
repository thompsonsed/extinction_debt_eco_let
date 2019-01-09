# Table generation for some conceptual numbers
 # TODO remove this 
library(tidyverse)
library(viridis)
library(akima)
source("preston.R")
source("plot_colours.R")
#################
# Data-Figure 6 #
#################
# Generate the data for the summary table styled figures in the main text.
# Load the range for dispersal values across real maps (calculated by the algorithm as described)
# then calculate the minimum and maximum proportion of species remaining at each spatial and
# temporal scale, with varying levels of habitat fragmentation and habitat loss.

table_df <- read.csv(file.path("results", "dispersal_ranges.csv")) %>% select(-X) %>%
  mutate(a_max = a_max ^2,
         area = a_max * proportion_cover) %>% 
  group_by(a_max, area, map_file, proportion_cover, number_steps, mean_distance) %>%
  summarise(sigma_e = estimate_sigma_rayleigh(mean(mean_distance), n_steps = mean(number_steps)),
            effective_connectivity_sq = (area/a_max)*sigma_e^2) %>%
  group_by(a_max, area, proportion_cover, number_steps) %>%
  summarise(min_sigma_e = min(sigma_e),
            max_sigma_e = max(sigma_e)) %>% 
  group_by(a_max, area, proportion_cover, number_steps, min_sigma_e, max_sigma_e) %>%
  expand(a_max, area, proportion_cover, number_steps,
         sigma_e = seq(min_sigma_e, max_sigma_e, 0.1)) %>%
  filter(number_steps == 10000) %>%
  mutate(speciation_rate = 1/number_steps,
         contig_richness = S_contig(a_max, speciation_rate, 16^2),
         inst_richness_upper = S_random(a_max, area, speciation_rate, 16^2),
         inst_richness_lower = S_contig(area, speciation_rate, 16^2),
         inst_richness = max_sigma_e,
         eq_richness = S_random_equilibrium(a_max, area, speciation_rate, sigma_e^2),
         prop_inst_upper = inst_richness_upper/contig_richness,
         prop_inst_lower = inst_richness_lower/contig_richness,
         prop_eq = eq_richness/contig_richness)
  
table_df <- read.csv(file.path("results", "dispersal_ranges.csv")) %>% select(-X) %>%
  mutate(a_max = a_max ^2,
         area = a_max * proportion_cover) %>% 
  group_by(a_max, area, map_file, proportion_cover, number_steps, mean_distance) %>%
  summarise(sigma_e = estimate_sigma_rayleigh(mean(mean_distance), n_steps = mean(number_steps)),
            effective_connectivity_sq = (area/a_max)*sigma_e^2) %>% 
  filter(number_steps == 10000) %>%
  mutate(speciation_rate = 1/number_steps,
         contig_richness = S_contig(a_max, 0.00001, 16^2),
         inst_richness_upper = S_random(a_max, area, 0.00001, 16^2),
         inst_richness_lower = S_contig(area, 0.00001, 16^2),
         eq_richness = S_random_equilibrium(a_max, area, 0.00001, sigma_e^2),
         prop_inst_upper = inst_richness_upper/contig_richness,
         prop_inst_lower = inst_richness_lower/contig_richness,
         prop_eq = eq_richness/contig_richness) %>%
  gather(key="key", value="proportion_richness", prop_inst_upper, prop_inst_lower,
         prop_eq) %>%
  mutate(time = ifelse(key %in% c("prop_inst_upper", "prop_inst_lower"), "Short-term", "Long-term"),
         rel_proportion_richness = proportion_richness/proportion_cover) %>%
  mutate(time=factor(time, levels=c("Short-term", "Long-term")),
         extinction_debt=1.0-proportion_richness)

interp_mat <- interp(x=table_df,
                     y=100*select_df$mean_richness/select_df$contig_richness,
                     z=select_df$estimated_sigma, duplicate="strip", extrap=TRUE, linear=TRUE,
                     xo=seq(2, 7, 0.05), yo=seq(0, 21, 0.1))
p6 <- table_df %>% ggplot() + 
  geom_raster(aes(fill=100*proportion_richness, y=sigma_e,
                x=time), interpolate=TRUE) + 
  theme_classic() + 
  facet_grid(as.factor(proportion_cover) ~ a_max)+
             # labeller=labeller(a_max=a_max_names_large, 
                               # proportion_cover = percent_cover_names2))+
  scale_y_continuous(element_blank())+ 
  scale_x_discrete("Spatiotemporal scale", position = "top")+
  # scale_fill_viridis("Percentage of\nspecies richness", 
                     # option="plasma", direction=-1, limits=c(0, 100), breaks=c(0, 100)) +
  theme(aspect.ratio = 1.5, legend.position = "bottom", 
        strip.placement = "outside")
        # strip.text.y = element_text(angle = 180),
        # axis.ticks.x=element_blank(),
        # axis.line.x=element_blank(),
        # strip.background = element_blank())
p6
p6 <- ggplot(table_df) + geom_tile(aes(fill=100*proportion_richness, y=as.factor(proportion_cover),
                                       x=time)) + 
  facet_grid(fragmentation ~ a_max, labeller=labeller(a_max=a_max_names_large),
             switch="y") +
  theme_classic()+ 
  scale_y_discrete(element_blank(), 
                   labels=c("20% habitat\nremaining", "40% habitat\nremaining",
                            "80% habitat\nremaining"))+ 
  scale_x_discrete("Spatiotemporal scale", position = "top")+
  scale_fill_viridis("Percentage of\nspecies richness", 
                     option="plasma", direction=-1, limits=c(0, 100)) + 
  theme(aspect.ratio = 1.5, legend.position = "bottom", strip.placement = "outside",
        strip.text.y = element_text(angle = 180),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        strip.background = element_blank())
pdf(file.path(figure_dir, "figure6.pdf"), 6.81, 6)
print(p6)
dev.off()
