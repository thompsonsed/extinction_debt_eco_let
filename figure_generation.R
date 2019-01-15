## Figure reproduction based on simulation outputs

library(tidyverse)
library(viridis)
library(lemon)
library(cowplot)
library(ggpubr)
# Install ggthemr with devtools::install_github('cttobin/ggthemr') - requires devtools package
library(magick) # Required for cowplot to plot images
library(ggthemr)
source("preston.R")
source("plot_colours.R")
##################
## Data reading ##
##################

# Set the input and output directories
# Folder containing the simulation results
data_dir <- "results" 
# Folder to output all figures generated here to
figure_dir <- "figures"
# Folder which holds the example map files
map_dir <- "maps"

dir.create(figure_dir)
dir.create(file.path(figure_dir, "appendices"))

#################
# Data-Figure 1 #
#################

# Read the selected example simulation data and manipulate to the desired format
example_sim_data <- read.csv(file.path(data_dir, "single_sim_results.csv")) %>% select(-X) %>% 
  mutate(time_since_pristine = 2501 - time) %>% unique() %>%
  group_by(map_file, speciation, time_since_pristine) %>%
  summarise(mean_richness = mean(richness),
            min_richness = quantile(richness)[1],
            max_richness = quantile(richness)[5],
            count=n())
example_sim_data$map_file <- factor(example_sim_data$map_file,
                                    levels = c("random_500_0.2_48.tif",
                                               "contrived_500_0.2_1024.tif",
                                               "map_500_0.2_1.tif"),
                                    labels = c("Random",
                                               "Clustered",
                                               "Real"))
# Find the mean richness prior clearing (should be 1536)
mean_prior_clearing <- (example_sim_data %>% ungroup %>%
  filter(time_since_pristine < 0 & time_since_pristine > -200) %>%
  summarise(mean_richness = mean(mean_richness)))$mean_richness
# Calculate the average species loss for each scenario
example_sim_averages <- example_sim_data %>% ungroup() %>%
  filter(time_since_pristine == -997499, speciation == "False") %>% 
  select(mean_richness, map_file) %>%
  mutate(mean_percentage_remaining = 100*mean_richness/mean_prior_clearing,
         mean_richness_lost = mean_prior_clearing - mean_richness,
         mean_percentage_lost = 100-mean_percentage_remaining)
# Generate the summary statistics for plotting the equilibrium lines
example_sim_summary <- example_sim_data%>% group_by(map_file, speciation) %>% 
  filter(speciation == "False") %>%
  summarise(contig_richness = mean_richness[which(time_since_pristine == 0)],
            immediate_richness_max = max_richness[which(time_since_pristine == 1)],
            immediate_richness_min = min_richness[which(time_since_pristine == 1)],
            eq_richness_max = max_richness[which(time_since_pristine == -997499)],
            eq_richness_min = min_richness[which(time_since_pristine == -997499)])
example_sim_summary_2 <- example_sim_summary %>%
  mutate(immediate_loss = contig_richness - immediate_richness_max,
         imm_or_ex = immediate_richness_max - immediate_richness_min,
         extinction_debt = immediate_richness_min - eq_richness_max,
         ex_or_remaining = eq_richness_max - eq_richness_min,
         remaining = eq_richness_min) %>% 
  select(map_file, immediate_loss, imm_or_ex, extinction_debt, ex_or_remaining, remaining) %>%
  gather(key=time, value=richness, -map_file) %>%
  mutate(time = factor(time, levels=c("immediate_loss", "imm_or_ex", "extinction_debt",
                                      "ex_or_remaining", "remaining")))
example_sim_summary <- example_sim_summary %>%
  select(-speciation) %>% gather(key=time, value=richness, -map_file)

###################
# Data-Figure 2,3 #
###################

# Read simulated data from main set of simulations
coalescence_data <- read.csv(file.path(data_dir, "coalescence.csv"), header=TRUE, sep=" ") %>%
  separate(file, into=c("type", "size", "proportion_cover"), sep="_", extra="drop", remove=FALSE) %>%
  mutate(type=ifelse(type == "contrived", "Clustered", ifelse(type=="map", "Real", "Random"))) %>%
  select(-job)
dispersal_data <- read.csv(file.path(data_dir, "dispersal.csv"), header=TRUE, sep=" ") %>%
  separate(file, into=c("type", "size", "proportion_cover"), sep="_", extra="drop", remove=FALSE) %>%
  mutate(type=ifelse(type == "contrived", "Clustered", ifelse(type=="map", "Real", "Random"))) %>%
  select(-c(seed, mean_dispersal, mean_distance_travelled, 
            mean_distance_travelled_100, mean_distance_travelled_1000, stdev_dispersal, 
            stdev_distance_travelled, stdev_distance_travelled_100, stdev_distance_travelled_1000))

main_df <- merge(coalescence_data, dispersal_data, 
                 by = c("file", "sigma", "type", "size", "proportion_cover")) %>% 
  mutate(size = as.numeric(size),
         proportion_cover = as.numeric(proportion_cover),
         area=size*size*proportion_cover)
# Add in data from contiguous landscape simulations
# Import the contiguous data and bind it to the same data frame
contiguous_df <- read.csv(file.path(data_dir, "contiguous_richness.csv"))
contiguous_dispersal <- read.csv(file.path(data_dir, "dispersals_contiguous.csv")) %>%
  select(-c(seed, mean_distance_travelled_100000, mean_dispersal, mean_distance_travelled, 
            mean_distance_travelled_100, mean_distance_travelled_1000, stdev_dispersal, 
            stdev_distance_travelled, stdev_distance_travelled_100, stdev_distance_travelled_1000))
contiguous_df$type <- "Contiguous"
contiguous_df$proportion_cover <- 1
contiguous_df$size <- contiguous_df$area^0.5
contiguous_df_combined <- merge(contiguous_df, contiguous_dispersal)
main_df <- rbind(main_df, contiguous_df_combined)

# Calculate the means for each combination of sigma, size, proportion cover
main_df <- main_df %>% group_by(file, type, sigma, size, proportion_cover, 
                                mean_distance_travelled_10000) %>% 
  summarise(richness_mean= mean(richness), richness_stdev= sd(richness)) %>%
  rename(richness = richness_mean)

# Calculate our sigma_e and effective connectivity metric 
main_df <- main_df %>% filter(proportion_cover != 0.5, proportion_cover != 0.25) %>% ungroup()  %>%
  mutate(area=size*size*proportion_cover,
         a_max = size^2,
         sigma_e = estimate_sigma_rayleigh(mean_distance_travelled_10000, 10000),
         effective_connectivity = sqrt(area/a_max) * sigma_e,
         preston_richness = S_contig(area, 0.0001, effective_connectivity^2)) %>% 
  filter(size<6000) %>% 
  mutate(type = factor(type, labels=c("Random", "Clustered", "Real", "Contiguous"),
                       levels=c("Random", "Clustered", "Real", "Contiguous")))

#################
# Data-Figure 3 #
#################

# Calculate the analytical approximations of extinction debt
analytical_approx_ED <- expand.grid(proportion_cover=c(0.1, 0.2, 0.4), sigma=c(8, 16, 32), 
                          a_max=10^seq(1, 9, 0.1)) %>% 
  mutate(area = proportion_cover * a_max,
         s_inst_lower = S_contig(area, 0.0001, sigma_sq = sigma^2),
         s_inst_upper = S_random(a_max, area, 0.0001, sigma^2),
         contig_diversity = S_contig(a_max, 0.0001, sigma^2),
         min_pc_inst = s_inst_lower* 100/contig_diversity,
         max_pc_inst = s_inst_upper * 100/contig_diversity)
# Calculate the minimum effective connectivity for each size
min_effective_connectivity <- main_df %>% ungroup%>% group_by(sigma) %>% 
  summarise(min_sigma_e = max(0.45, min(na.omit(estimate_sigma_rayleigh(mean_distance_travelled_10000,
                                                                        10000)))),
            max_sigma_e = max(na.omit(estimate_sigma_rayleigh(mean_distance_travelled_10000,
                                                              10000))))
analytical_approx_ED <- merge(min_effective_connectivity, analytical_approx_ED) %>% 
  mutate(min_equilibrium = S_random_equilibrium(a_max, area, 0.0001, min_sigma_e^2),
         max_equilibrium = S_random_equilibrium(a_max, area, 0.0001, max_sigma_e^2),
         min_pc_equilibrium = min_equilibrium*100/contig_diversity,
         max_pc_equilibrium = max_equilibrium*100/contig_diversity)

#################
# Data-Figure 4 #
#################

variable_size_data <- read.csv(file.path(data_dir, "variables_sizes.csv")) %>% select(-X) %>%
  filter(sigma == 16) %>%
  mutate(file=str_replace(file, pattern=".*/", replacement=""))  %>%
  group_by(file, sigma, type, speciation_rate) %>% 
  summarise(mean_richness = mean(richness),
            count=n()) 
variable_size_dispersal <- read.csv(file.path(data_dir, "variable_sizes_dispersal.csv")) %>% 
  select(-X) %>% filter(sigma == 16)
variable_size_data <- variable_size_data %>% full_join(variable_size_dispersal) %>% 
  mutate(sigma_e = estimate_sigma_rayleigh(mean_distance_10000, 10000),
         area = a_max * 0.2,
         effective_connectivity = sqrt(area/a_max) * sigma_e) %>% 
  filter(a_max <= 70000000, sigma==16, speciation_rate==0.0001) %>% 
  ungroup() %>%
  mutate(min_sigma_e = min(sigma_e),
         max_sigma_e = 16,
         preston_richness = S_random_equilibrium(a_max, area, 0.0001, sigma_e^2),
         contig_richness = S_contig(a_max, 0.0001, sigma^2),
         max_richness_inst = S_random(a_max, area, 0.0001, sigma^2),
         min_richness_inst = S_contig(area, 0.0001, sigma^2),
         max_richness_eq = S_random_equilibrium(a_max, area, 0.0001, max_sigma_e^2),
         min_richness_eq = S_random_equilibrium(a_max, area, 0.0001, min_sigma_e^2),
         max_pc_inst = 100*max_richness_inst / contig_richness,
         min_pc_inst = 100*min_richness_inst / contig_richness,
         min_pc_eq = 100*min_richness_eq/contig_richness,
         max_pc_eq = 100*max_richness_eq/contig_richness,
         type = as.factor(ifelse(type == "Contrived", "Clustered", as.character(type))))

# Calculate the bounds on the real landscapes for each size
real_bounds <- variable_size_data %>% filter(type == "Real") %>% 
  ungroup() %>%
  group_by(area) %>% 
  mutate(real_pc = 100*mean_richness/contig_richness,
         min_real_pc = min(real_pc),
         max_real_pc = max(real_pc))


#################
# Data-Figure 5 #
#################
# Generate the data for the summary table-styled figure in the main text.
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
         prop_eq_lower = eq_richness_lower/contig_richness) %>%
  gather(key="key", value="proportion_richness", prop_inst_upper, prop_inst_lower,
         prop_eq_upper, prop_eq_lower) %>%
  mutate(fragmentation=ifelse(key %in% c("prop_inst_upper", "prop_eq_upper"), 
                              "High effective\nconnectivity", 
                              "Low effective\nconnectivity"),
         time = ifelse(key %in% c("prop_inst_upper", "prop_inst_lower"), "Short-term", "Long-term"),
         rel_proportion_richness = proportion_richness/proportion_cover) %>%
  ungroup %>%
  mutate(time=factor(time, levels=c("Short-term", "Long-term")),
         fragmentation=factor(fragmentation, levels=c("Low effective\nconnectivity",
                                                      "High effective\nconnectivity")),
         proportion_cover = factor(proportion_cover, levels=c(0.8, 0.4, 0.2)),
         extinction_debt=1.0-proportion_richness)

#################
# Data-Figure 6 #
#################

zeroes_df <- data.frame(expand.grid(area=0.0, speciation_rate=c(1e-8, 1e-6, 1e-4))) %>%
  mutate(prop_cover =0.0,
         sigma=16, a_max=1000^2,
         z = ifelse(speciation_rate == 1e-8, 0.1, ifelse(speciation_rate == 1e-6, 0.2, 0.3)),
         a_max_str=as.character(sqrt(a_max)),
         pc_cover_str = as.character(0.0),
         best_case_richness = 0.0,
         contig_richness=S_contig(a_max, speciation_rate, sigma^2),
         best_case_instant = 0.0,
         worst_case_instant=0.0,
         contiguous_area_est=0.0,
         power_law_est = 0.0,
         best_case_pc = 0.0)

analytical_approx_best_case <- data.frame(expand.grid(prop_cover=seq(0.01, 0.99, 0.001),
                                                      speciation_rate=c(1e-8, 1e-6, 1e-4),
                                                      sigma=16, a_max=1000^2)) %>%
  mutate(area=a_max * prop_cover,
         z = ifelse(speciation_rate == 1e-8, 0.1, ifelse(speciation_rate == 1e-6, 0.2, 0.3)),
         a_max_str = as.character(sqrt(a_max)),
         pc_cover_str = as.character(prop_cover),
         best_case_richness = S_random_equilibrium(a_max, area, speciation_rate, sigma^2),
         contig_richness = S_contig(a_max, speciation_rate, sigma^2),
         best_case_instant = 100*S_random(a_max, area, speciation_rate, sigma^2)/contig_richness,
         worst_case_instant = 100*S_contig(area, speciation_rate, sigma^2)/contig_richness,
         contiguous_area_est = 100*S_contig(area, speciation_rate, sigma^2)/contig_richness,
         power_law_est = 100*power_law_estimation(a_max, area, z),
         best_case_pc = 100*best_case_richness/contig_richness) %>%
  rbind(zeroes_df)
###############
## Main text ##
###############

############
# Figure 1 #
############
p2 <- ggplot(example_sim_data %>% filter(time_since_pristine>=0)) +
  geom_smooth(data=example_sim_data %>%
                filter(time_since_pristine >= 5, time_since_pristine <=1001),
              aes(x=time_since_pristine, y=mean_richness, linetype=speciation),
              fill=NA, colour="black", size=0.75,
              method="glm", formula = y ~ poly(x,9))+
  geom_line(data=example_sim_data %>% filter(time_since_pristine <= 6, 
                                             time_since_pristine >=0),
            aes(x=time_since_pristine, y=mean_richness, linetype=speciation, group=speciation),
            size=0.75, colour="black")+
  theme_classic() +
  theme(legend.key.height=unit(2,"line")) + 
  geom_bar(data=example_sim_summary_2, 
           aes(x=1200, y=richness, fill=time), width=200, stat = "identity")+#, position="dodge") +
  geom_segment(data=example_sim_summary,
               aes(x=0, xend=1100, y=richness, yend = richness, colour=time),
               linetype="dotted", size=0.25) +
  scale_x_continuous("Generations since habitat loss", breaks=c(0, 200, 400, 600, 800, 1000),
                     limits=c(0, 1300), expand = c(0, 20))+
  scale_y_continuous("Species richness", expand = c(0, 0), limits = c(0, 1600)) + 
  scale_fill_manual("", labels=c("Immediate\nloss",
                                 "Immediate loss\nor extinction debt",
                                 "Extinction\ndebt",
                                 "Remaining or\nextinction debt",
                                 "Remaining"),
                    breaks=c("immediate_loss", "imm_or_ex", "extinction_debt",
                             "ex_or_remaining", "remaining"),
                    values=c("immediate_loss"=plot_colours[1],
                             "imm_or_ex"=plot_colours[3],
                             "extinction_debt"=plot_colours[7],
                             "ex_or_remaining"=plot_colours[9],
                             "remaining"=plot_colours[11])) +
  scale_colour_manual("", labels=c("Immediate\nloss",
                                   "Immediate loss\nor extinction debt",
                                   "Extinction\ndebt",
                                   "Remaining or\nextinction debt",
                                   "Remaining"),
                      breaks=c("contig_richness", "immediate_richness_max", "immediate_richness_min",
                               "eq_richness_max", "eq_richness_min"),
                      values=c("contig_richness"=plot_colours[1],
                               "immediate_richness_max"=plot_colours[3],
                               "immediate_richness_min"=plot_colours[3],
                               "eq_richness_max"=plot_colours[9],
                               "eq_richness_min"=plot_colours[11]), guide=FALSE)+
  guides(linetype=guide_legend(keywidth = 2, keyheight = 0.5))+
  theme(panel.spacing = unit(0.75, "lines"),
        legend.position = "bottom",
        legend.box="vertical",
        # legend.direction = "horizontal",
        legend.text.align = 0, axis.line=element_line(),
        legend.text=element_text(size=7),
        legend.key.width = unit(0.4,"cm"),
        legend.key.height = unit(0.4,"cm"),
        legend.spacing = unit(0.1, "cm"),
        legend.margin=margin())+
  facet_rep_grid(map_file~.) + 
  scale_linetype_manual("", breaks=c("True", "False"), labels=c("With\nspeciation", 
                                                                "Without\nspeciation"),
                        values=c("dashed", "solid")) 
# Arrange the maps on the plot as well
l <- get_legend(p2)
units <- c(-0.5,-1.25,-1.5,0.5)
im1 <- ggdraw() + 
  draw_image(file.path(map_dir, "random_500_0.2_48.png")) + 
  theme(plot.margin = unit(units, "lines"))
# im1
im2 <- ggdraw() + 
  draw_image(file.path(map_dir, "contrived_500_0.2_1024.png"))+ 
  theme(plot.margin = unit(units, "lines"))
im3 <- ggdraw() + 
  draw_image(file.path(map_dir, "map_500_0.2_1.png"))+ 
  theme(plot.margin = unit(units, "lines"))
gga1 <- ggarrange(im1, im2, im3, NULL, labels=c("a)", "b)", "c)", NA), nrow=4, ncol=1, 
                  heights=c(1, 1, 1, 0.25))
gga2 <- ggarrange(NULL, NULL, NULL, NULL, labels=c("d)", "e)", "f)", NA), nrow=4, ncol=1, 
                  heights=c(1, 1, 1, 0.25), label.x=0.3)
gga2b <- ggarrange(NULL, p2 + theme(legend.position = "none"), nrow=2, ncol=1, 
                   heights=c(0.01, 1))
gga3 <- ggarrange(gga2, gga2b, nrow=1, ncol=2,
                  widths=c(0.14, 1))
gga4 <- ggarrange(gga1, gga3, nrow=1, ncol=2,
                  widths=c(0.40, 1))
gga5 <- ggarrange(gga4, l, nrow=2, ncol=1, heights=c(1, 0.1))
pdf(file.path(figure_dir, "figure1.pdf"), 4.3, 6)
print(gga5)
dev.off()

############
# Figure 2 #
############
ggthemr("light")

p1 <- main_df %>% filter(sigma > 4) %>% 
  ggplot(aes(x=area,
             y=richness))+
  theme_classic() + scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  xlab(expression(paste("Area (", A[e], ")"))) + geom_point(aes(colour = type),
                                                                             alpha=0.4)+
  ylab(expression(paste("Species richness (", S, ")"))) + 
  scale_colour_discrete("Landscape type")+
  scale_shape_discrete(expression(sigma))+
  guides(colour = guide_legend(override.aes = list(alpha = 1),
                               nrow=2, byrow=TRUE,
                               title.position="top",
                               title.hjust = 0.5),
         shape = guide_legend(override.aes = list(alpha = 1),
                              nrow=2, byrow=TRUE, title.position="top",
                              title.hjust = 0.5))+
  theme(aspect.ratio = 1.0, plot.title = element_text(size = 10)) + 
  ggtitle("Unscaled species-area curve")
p2 <- main_df %>% filter(sigma > 4) %>% 
  ggplot(aes(x=area/effective_connectivity^2,
             y=richness/effective_connectivity^2))+
  theme_classic() + scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  xlab(expression(paste("Rescaled area (",A[e] / c[e]^2, ")"))) + geom_point(aes(colour = type),
                                                                             alpha=0.4)+
  ylab(expression(paste("Rescaled species richness (", S/c[e]^2, ")"))) + 
  scale_colour_discrete("Landscape type")+
  scale_shape_discrete(expression(sigma))+
  guides(colour = guide_legend(override.aes = list(alpha = 1),
                               nrow=2, byrow=TRUE,
                               title.position="top",
                               title.hjust = 0.5),
         shape = guide_legend(override.aes = list(alpha = 1),
                              nrow=2, byrow=TRUE, title.position="top",
                              title.hjust = 0.5))+

  theme(aspect.ratio = 1.0, plot.title = element_text(size = 10)) + 
  ggtitle("Rescaled species-area curve")
gga2 <- ggarrange(p1, p2, common.legend = TRUE, legend="bottom", nrow = 2, ncol=1,
                 labels = c("a)", "b)"), heights = c(0.95, 1.0))
pdf(file.path(figure_dir, "figure2.pdf"), 3.23, 7)
print(gga2)
dev.off()
ggthemr_reset()


############
# Figure 3 #
############
p <- ggplot(analytical_approx_ED, colour="black") + theme_classic() + 
  geom_ribbon(aes(x=a_max, ymin=max_pc_inst, ymax=100,
                  fill="Definite Immediate loss")) +
  geom_ribbon(aes(x=a_max, ymin=min_pc_inst, ymax=max_pc_inst,
                  fill="Possible Immediate loss")) +
  geom_ribbon(aes(x=a_max, ymin=max_pc_equilibrium, ymax=min_pc_inst,
                  fill="Definite long-term loss")) +
  geom_ribbon(aes(x=a_max, ymin=min_pc_equilibrium, ymax=max_pc_equilibrium,
                  fill="Possible remaining")) +
  geom_ribbon(aes(x=a_max, ymin=0, ymax=min_pc_equilibrium,
                  fill="Definite remaining")) +
  ylab("Percentage of species richness")+
  scale_x_log10(expression(paste("Total area (", A[max], ")")),limits=c(10^3, 10^9),
                breaks = scales::trans_breaks("log10", function(x) 10^x, n=4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_fill_manual("Scenario", labels=c("Immediate loss",
                                         "Immediate loss\nor extinction debt",
                                         "Extinction debt",
                                         "Remaining or\nextinction debt",
                                         "Remaining"),
                    breaks=c("Definite Immediate loss",
                             "Possible Immediate loss",
                             "Definite long-term loss",
                             "Possible remaining",
                             "Definite remaining"),
                    values=c("Definite Immediate loss"=plot_colours[1],
                             "Possible Immediate loss"=plot_colours[3],
                             "Definite long-term loss"=plot_colours[7],
                             "Possible remaining"=plot_colours[9],
                             "Definite remaining"=plot_colours[11]))+
  facet_grid(sigma~proportion_cover, labeller = labeller(proportion_cover=percent_cover_names,
                                                         sigma=sigma_names)) + 
  theme(legend.key.height=unit(2,"line")) + 
  theme(legend.key.width=unit(2,"line")) + 
  theme(aspect.ratio=1)
pdf(file.path(figure_dir, "figure3.pdf"), 6.7, 5)
print(p)
dev.off()

############
# Figure 4 #
############

# Highlight a few results on the plot
highlighted_sims <- variable_size_data %>% 
  filter(file %in% c("map_398_0.2_2.tif",
                     "random_398_0.2_2.tif",
                     "contrived_select_398_16.tif"))
px <- ggplot(variable_size_data, colour="black") +
  theme_classic() + 
  geom_ribbon(data=real_bounds,
              aes(x=a_max, ymin=min_real_pc, ymax=max_real_pc,
                  fill="Real landscapes"),  alpha=0.6)+
  geom_point(aes(x=a_max, y=100*mean_richness/contig_richness, colour=effective_connectivity))+
  
  geom_point(data=highlighted_sims,
             aes(x=a_max, y=100*mean_richness/contig_richness,
                 colour=effective_connectivity), shape=17, alpha=1, size=3)+
  
  geom_line(aes(x=a_max, y=max_pc_eq, linetype="Upper"),
            colour="grey20", size=1) +
  geom_line(aes(x=a_max, y=min_pc_eq, linetype="Lower"), colour="grey20", size=1) +
  ylab("Percentage of species\nrichness at equilibrium")+
  
  scale_x_log10(expression(paste("Total area (", A[max], ")")),limits=c(10^2.5, 10^8),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_fill_manual(NULL, label=c("Real\nlandscapes"), values=c("grey80"))+
  scale_linetype_manual("Theoretical\nlimits",breaks=c("Upper", "Lower"),
                        label=c("Upper", "Lower"), values=c("dotted", "dashed"))+
  guides(linetype=guide_legend(keywidth = 2.5))+
  scale_colour_viridis(expression(paste("Effective\nConnectivity (", c[e], ")")),
                       option="inferno", end=0.94)+
  theme( aspect.ratio = 1, legend.position = "left",
         legend.text=element_text(size=9),
         legend.title=element_text(size=10))
pdf(file.path(figure_dir, "figure4.pdf"), 5.57, 4)
print(px)
dev.off()
############
# Figure 5 #
############

p5 <- ggplot(table_df) + geom_tile(aes(fill=100*proportion_richness, y=fragmentation,
                                       x=time)) + 
  facet_grid(as.factor(proportion_cover) ~ a_max,
             labeller=labeller(a_max=a_max_names_large, 
                               "as.factor(proportion_cover)"=percent_cover_names2),
             switch="y") +
  theme_classic()+ 
  scale_y_discrete(element_blank())+
  # labels=c("20% habitat\nremaining", "40% habitat\nremaining",
  # "80% habitat\nremaining"))+ 
  scale_x_discrete("Spatiotemporal scale", position = "top")+
  scale_fill_viridis("Percentage of species\nrichness remaining", 
                     option="inferno", direction=-1, limits=c(0, 100), breaks=c(0, 100),
                     begin = 0.1, end=0.9) + 
  theme(aspect.ratio = 1, legend.position = "bottom", strip.placement = "outside",
        strip.text.y = element_text(angle = 180),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        strip.background = element_blank())
pdf(file.path(figure_dir, "figure5.pdf"), 6.81, 6)
print(p5)
dev.off()

############
# Figure 6 #
############
# Define the colours for the theme
ggthemr_light_colours = c("#785d37", "#62bba5", "#ffb84d", "#aaa488", "#b2432f", "#3a6589", 
                          "#9b5672", "#908150", "#373634")
p1 <- analytical_approx_best_case %>% 
  ggplot() + 
  geom_line(aes(x=prop_cover*100, y=power_law_est, colour=as.factor(z), linetype=as.factor(z))) + 
  geom_line(aes(x=prop_cover*100, y=best_case_pc, colour="Long-term\n(best case)",
                linetype="Long-term\n(best case)"))+
  theme_classic() + 
  scale_colour_manual(element_blank(), breaks=c(0.1, 0.2, 0.3, "Long-term\n(best case)"),
                      labels=c(expression(z == 0.1, z == 0.2, z == 0.3),
                               "Long-term\n(best case)"), 
                      values=c(ggthemr_light_colours[1], ggthemr_light_colours[3], 
                               ggthemr_light_colours[5], plot_colours[9]))+
  scale_linetype_manual(element_blank(), breaks=c(0.1, 0.2, 0.3, "Long-term\n(best case)"),
                        labels=c(expression(z == 0.1, z == 0.2, z == 0.3),
                                 "Long-term\n(best case)"), 
                        values=c("solid", "solid", "solid", "dotted"))+
  xlab("Percentage of habitat remaining")+
  ylab("Percentage of species\nrichness at equilibrium")+
  ylim(c(0, 100))+
  theme(aspect.ratio=1, legend.key.size = unit(1.5, 'lines'),
        legend.position = c(0.78, 0.28),
        legend.background = element_rect(size = 0.5, colour = 1),
        legend.title = element_blank(),
        legend.margin = margin(-0.25, 0.5, 0.25, 0.5, unit = "lines"),
        legend.text = element_text(size=8),
        plot.title = element_text(size=10)) + 
  ggtitle("Power law approach vs\nlong-term Preston function (best case)")

p2 <- analytical_approx_best_case %>% filter(speciation_rate == 10^-6) %>%
  ggplot() + 
  geom_line(aes(x=prop_cover*100, y=best_case_instant,
                  linetype="Short-term\n(best case)", colour="Short-term\n(best case)")) +
  geom_line(aes(x=prop_cover*100, y=worst_case_instant,
                linetype="Short-term\n(worst case)", colour="Short-term\n(worst case)")) +
  geom_line(aes(x=prop_cover*100, y=best_case_pc, linetype="Long-term\n(best case)",
                colour="Long-term\n(best case)"))+
  theme_classic() + 
  scale_linetype_manual(element_blank(),
                          breaks=c("Short-term\n(best case)", "Short-term\n(worst case)",
                                   "Long-term\n(best case)"),
                        values=c("Short-term\n(best case)"="solid",
                                 "Short-term\n(worst case)"="dashed", 
                                 "Long-term\n(best case)"="dotted"))+
  scale_colour_manual(element_blank(),
                      breaks=c("Short-term\n(best case)", "Short-term\n(worst case)",
                                   "Long-term\n(best case)"),
                      values=c("Short-term\n(best case)"=plot_colours[1],
                               "Short-term\n(worst case)"=plot_colours[3], 
                               "Long-term\n(best case)"=plot_colours[9]))+
  xlab("Percentage of habitat remaining")+
  ylab("Percentage of species\nrichness at equilibrium")+
  ylim(c(0, 100))+
  theme(aspect.ratio=1, legend.key.size = unit(1.5, 'lines'),
        plot.title = element_text(size=10),
        legend.margin = margin(-0.25, 0.5, 0.25, 0.5, unit = "lines"),
        legend.position = c(0.76, 0.22),
        legend.background = element_rect(size = 0.5, colour = 1),
        legend.title = element_text(size=10),
        legend.text = element_text(size=8)) + 
  ggtitle("Short-term Preston function vs\nlong-term Preston function (best case)")

gga1 <- ggarrange(p1, p2, labels = c("a)",
                                     "b)"), ncol=1, nrow=2)
pdf(file.path(figure_dir, "figure6.pdf"), 3.23, 7.0)
print(gga1)
dev.off()
ggthemr_reset()
################
## Appendices ##
################

##############
# Appendix 3 #
##############

# Figure 1
# The scaling collapse in the random scenario

p <- main_df %>% filter(type=="Random" | type == "Contiguous", sigma>2) %>% 
  mutate(scaled_sigma = sqrt(area/a_max)*sigma) %>%
  ggplot(aes(x=area/scaled_sigma^2, y=richness/scaled_sigma^2))+
  geom_point(aes(colour=sigma, shape=as.factor(100*proportion_cover))) +theme_classic() + 
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n=4),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n=4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  xlab(expression(A[max]/ sigma^2)) + 
  ylab(expression(S*A[max]/(A[e]*sigma^2))) + 
  scale_colour_viridis(expression(sigma), trans="log2")+
  geom_smooth(data= main_df %>% filter(type == "Contiguous", sigma>2) %>% 
                mutate(scaled_sigma = sqrt(area/a_max)*sigma),
              fill=NA,linetype="dotted", colour="grey")+
  scale_shape_discrete("Percentage\ncover")
pdf(file.path(figure_dir, "appendices", "appendix3_s1.pdf"), 6, 4)
print(p)
dev.off()

# Figure 2
# Lack of scaling collapse for random solution in real and contrived scenarios
ggthemr('light')
p1 <- main_df %>% filter(sigma > 4) %>% 
  ggplot(aes(x=area,
             y=richness))+
  theme_classic() + scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  xlab(expression(paste("Area (", A[e], ")"))) +
  geom_point(aes(colour = type, shape=as.factor(sigma)), alpha=0.4)+
  ylab(expression(paste("Species richness (", S, ")"))) + 
  scale_colour_discrete("Landscape type")+
  scale_shape_discrete(expression(sigma))+
  guides(colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = 1)),
         shape = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = 1)))+
  theme(aspect.ratio = 1.0, plot.title = element_text(size = 10)) + 
  ggtitle("Unscaled species-area curve")
p2 <-main_df %>% filter(sigma > 4) %>% ggplot()+
  theme_classic() + 
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n=4),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n=4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_point(aes(y = richness/(sigma^2),
                 x = area/sigma^2,
                 colour = type,
                 shape=as.factor(sigma)),alpha=0.4)+
  xlab(expression(paste("Scaled area (", A[e]/sigma^2, ")"))) + 
  ylab(expression(paste("Scaled species richness (", S/sigma^2, ")"))) +
  scale_colour_discrete("Landscape type", 
                        labels=c( "Random","Clustered", "Real", "Contiguous"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1)),
         shape = guide_legend(override.aes = list(alpha = 1)))+
  scale_shape_discrete(expression(sigma))+
  theme(aspect.ratio=1, plot.title = element_text(size = 10)) + 
  ggtitle(expression(paste("Scaling with ", sigma^2)))+
  guides(colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = 1)),
         shape = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = 1)))
p3 <- main_df %>% filter(sigma > 4) %>% ggplot()+
  theme_classic() + scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n=4),
                                  labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n=4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_point(aes(y = richness/(proportion_cover * sigma^2),
                 x = area/sigma^2,
                 colour = type,
                 shape=as.factor(sigma)),alpha=0.4)+
  xlab(expression(paste("Scaled area (", A[e]/omega^2, ")"))) + 
  ylab(expression(paste("Scaled species richness (", S/omega^2, ")"))) + 
  scale_colour_discrete("Landscape type", 
                        labels=c( "Random","Clustered", "Real", "Contiguous"))+
  scale_shape_discrete(expression(sigma)) + 
  theme(aspect.ratio=1, plot.title = element_text(size = 10)) + 
  ggtitle(expression(paste("Scaling with ", omega^2==(A[e]/A[max])*sigma^2)))+
  guides(colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = 1)),
         shape = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = 1)))
p4 <- main_df %>% filter(sigma > 4) %>% ggplot()+
  theme_classic() + scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n=4),
                                  labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n=4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_point(aes(y = richness/(effective_connectivity^2),
                 x = area/(effective_connectivity^2),
                 colour = type,
                 shape=as.factor(sigma)),alpha=0.4)+
  xlab(expression(paste("Scaled area (", A[e]/c[e]^2, ")"))) + 
  ylab(expression(paste("Scaled species richness (", S/c[e]^2, ")"))) + 
  scale_colour_discrete("Landscape type", 
                        labels=c( "Random","Clustered", "Real", "Contiguous"))+
  scale_shape_discrete(expression(sigma)) + 
  theme(aspect.ratio=1, plot.title = element_text(size = 10)) + 
  ggtitle(expression(paste("Scaling with ", c[e]^2==(A[e]/A[max])*sigma[e]^2)))+
  guides(colour = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = 1)),
         shape = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(alpha = 1)))
gga <- ggarrange(p1, p2, p3, p4, common.legend = TRUE, legend="bottom", nrow = 2, ncol=2,
                 labels = c("a)", "b)", "c)", "d)"), heights = c(0.95, 1.0), widths=c(0.92, 1.0))
pdf(file.path(figure_dir, "appendices", "appendix3_s2.pdf"), 6.7, 7.5)
print(gga)
dev.off()
ggthemr_reset()

# Figure 3
# The simulated vs analytical results to indicate accuracy of the method

# Determine the mean actual percentage error between two vectors
mape <- function (x, y){
  df <- data.frame(x, y)
  df <- df %>% rowwise() %>% mutate(prop_err = 100*((x-y)^2)^0.5/max(x, y))
  return(sum(df$prop_err)/length(x))
}
mapes_df <- main_df %>% ungroup() %>% filter(sigma > 4) %>% group_by(type) %>% 
  summarise(mape = mape(richness, preston_richness)) %>% rowwise() %>%
  mutate(
    str_mape =as.character(as.expression(substitute(paste(MPE==r, "%"),
                                                         list(r=round(mape, 2))))))
p <- ggplot(main_df %>% filter(type!="Contiguous", sigma > 4)) + 
  geom_point(aes(x=preston_richness, y=richness, colour=area), alpha=0.1) +
  geom_point(data=main_df %>% filter(type=="Contiguous", sigma > 4),
             aes(x=preston_richness, y=richness, colour=area), alpha=0.7)+
  geom_abline(intercept = 0, slope = 1, linetype="dotted", colour="grey")+
  theme_classic() + 
  scale_x_log10(expression("Theoretical richness (" ~ Psi ~ ")")) +
  scale_y_log10("Simulated richness") +
  facet_grid(.~type) + 
  scale_colour_viridis("Area", trans="log10", option="magma",
                       breaks = scales::trans_breaks("log10", function(x) 10^x, n=3),
                       labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_text(size=3,data=mapes_df, 
            aes(x = 50, y = 10000, label = str_mape, group=type),
            parse=TRUE, inherit.aes = FALSE, colour="black") + 
  theme(aspect.ratio = 1)
pdf(file.path(figure_dir, "appendices", "appendix3_s3.pdf"), 6.7, 3)    
print(p)
dev.off()

##############
# Appendix 4 #
##############

# Create dummy data for plotting speciation rate and sigma vs percent remaining
dummy_df <- data.frame(expand.grid(speciation_rate=10^seq(-10, -4, 0.1), proportion_cover=c(0.1, 0.2, 0.4), 
                                   a_max=c(100^2, 1000^2, 10000^2), sigma=2^seq(2, 5, 0.5))) %>%
  rowwise() %>%
  mutate(area=a_max * proportion_cover,
         a_max_str = as.character(sqrt(a_max)),
         pc_cover_str = as.character(proportion_cover),
         richness_inst_upper = S_random(a_max, area, speciation_rate, sigma^2),
         richness_inst_lower= S_contig(area, speciation_rate, sigma^2),
         richness_eq_upper = S_random_equilibrium(a_max, area, speciation_rate, (sigma*1.01)^2),
         richness_eq_lower = S_random_equilibrium(a_max, area, speciation_rate, 0.5^2),
         contig_richness = S_contig(a_max, speciation_rate, sigma^2),
         pc_remaining_eq_upper =100*richness_eq_upper/contig_richness,
         pc_remaining_eq_lower =100*richness_eq_lower/contig_richness,
         pc_remaining_inst_upper = 100*richness_inst_upper/contig_richness,
         pc_remaining_inst_lower = 100*richness_inst_lower/contig_richness)



p <- dummy_df %>% filter(speciation_rate == 1e-5) %>% ggplot(aes(x=sigma)) + 
  geom_ribbon(aes(ymin=0, ymax=pc_remaining_eq_lower,
                  fill="Remaining"))+
  geom_ribbon(aes(ymin=pc_remaining_eq_lower, ymax=pc_remaining_eq_upper,
                  fill="Remaining\nor extinction debt"))+
  geom_ribbon(aes(ymin=pc_remaining_eq_upper, ymax=pc_remaining_inst_lower,
                  fill="Extinction debt"))+
  geom_ribbon(aes(ymin=pc_remaining_inst_lower, ymax=pc_remaining_inst_upper,
                  fill="Immediate loss\nor extinction debt"))+
  geom_ribbon(aes(ymin=pc_remaining_inst_upper, ymax=100,
                  fill="Immediate loss"))+
  facet_grid(a_max_str ~ pc_cover_str, labeller = labeller(pc_cover_str=percent_cover_names,
                                                           a_max_str=a_max_names)) + 
  scale_x_continuous(expression(paste("Dispersal (", sigma, ")")), trans="log2",
                     breaks = scales::trans_breaks("log2", function(x) 2^x, n=4)
                     # labels = scales::trans_format("log2", scales::math_format(2^.x))
  ) + 
  scale_y_continuous("Percentage of species richness remaining")+
  scale_fill_manual("Scenario", labels=c("Immediate loss",
                                         "Immediate loss\nor extinction debt",
                                         "Extinction debt",
                                         "Remaining\nor extinction debt",
                                         "Remaining"),
                    breaks=c("Immediate loss",
                             "Immediate loss\nor extinction debt",
                             "Extinction debt",
                             "Remaining\nor extinction debt",
                             "Remaining"),
                    # values=c("#FF7F80", "#FFB0B0", "#FCFCFC", "#C8BFFF", "#A997FF"))+
                    values=c("Immediate loss"=plot_colours[1],
                             "Immediate loss\nor extinction debt"=plot_colours[3],
                             "Extinction debt"=plot_colours[7],
                             "Remaining\nor extinction debt"=plot_colours[9],
                             "Remaining"=plot_colours[11])) + 
  theme_classic() + theme(aspect.ratio = 1) + 
  theme(legend.key.height=unit(2,"line")) + 
  theme(legend.key.width=unit(2,"line"))

pdf(file.path(figure_dir, "appendices", "appendix4_s4.pdf"), 7, 5)
print(p)
dev.off()

p <- dummy_df %>% filter(sigma == 16, speciation_rate >= 10^-9) %>% ggplot() + 
  geom_ribbon(aes(x=speciation_rate, ymin=0, ymax=pc_remaining_eq_lower,
                  fill="Remaining"))+
  geom_ribbon(aes(x=speciation_rate, ymin=pc_remaining_eq_lower, ymax=pc_remaining_eq_upper,
                  fill="Remaining\nor extinction debt"))+
  geom_ribbon(aes(x=speciation_rate, ymin=pc_remaining_eq_upper, ymax=pc_remaining_inst_lower,
                  fill="Extinction debt"))+
  geom_ribbon(aes(x=speciation_rate, ymin=pc_remaining_inst_lower, ymax=pc_remaining_inst_upper,
                  fill="Immediate loss\nor extinction debt"))+
  geom_ribbon(aes(x=speciation_rate, ymin=pc_remaining_inst_upper, ymax=100,
                  fill="Immediate loss"))+
  facet_grid(a_max_str ~ pc_cover_str, labeller = labeller(pc_cover_str=percent_cover_names,
                                                           a_max_str=a_max_names)) + 
  scale_x_log10(expression(paste("Speciation rate (", nu, ")")),
                breaks = c(10^-8, 10^-6, 10^-4),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous("Percentage of species richness remaining")+
  scale_fill_manual("Scenario", labels=c("Immediate loss",
                                         "Immediate loss\nor extinction debt",
                                         "Extinction debt",
                                         "Remaining\nor extinction debt",
                                         "Remaining"),
                    breaks=c("Immediate loss",
                             "Immediate loss\nor extinction debt",
                             "Extinction debt",
                             "Remaining\nor extinction debt",
                             "Remaining"),
                    # values=c("#FF7F80", "#FFB0B0", "#FCFCFC", "#C8BFFF", "#A997FF"))+
                    values=c("Immediate loss"=plot_colours[1],
                             "Immediate loss\nor extinction debt"=plot_colours[3],
                             "Extinction debt"=plot_colours[7],
                             "Remaining\nor extinction debt"=plot_colours[9],
                             "Remaining"=plot_colours[11])) + 
  theme_classic() + theme(aspect.ratio = 1) + 
  theme(legend.key.height=unit(2,"line")) + 
  theme(legend.key.width=unit(2,"line"))

pdf(file.path(figure_dir, "appendices", "appendix4_s5.pdf"), 7, 5)
print(p)
dev.off()
