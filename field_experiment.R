#### run analysis models for field experiment ####
library(tidyverse)
library(lmerTest)
library(effects)
library(emmeans)
library(boot)
library(lubridate)

#### temp plot ####
temps <- read.csv("temps.csv") %>% 
  mutate(Time = ymd_hms(Time))

str(temps)

col_dates <- read.csv("col_dates.csv") %>% 
  mutate(col_date = paste(col_date, "00:00:00", sep = " ")) %>% 
  mutate(col_date = mdy_hms(col_date)) %>% 
  mutate(site = plyr::mapvalues(site, from = c("Green Bay", "Illinois", "Madison"), to = c("Green Bay, WI", "Casey, IL", "Madison, WI"))) %>% 
  mutate(site = factor(site, levels = c("Green Bay, WI", "Madison, WI", "Casey, IL")))

# This is the graph that Jon made in import_temp_data.R but I wanted it here so I could look at it when thinking more about the temperature data
temp_graph <- temps %>% filter(treatment != "air") %>%  separate(Time, into = c("date", "time"), sep = " ", remove = F) %>% mutate(date = ymd(date), Time = ymd_hms(Time)) %>% filter(date < "2020-05-11" & date > "2020-01-01") %>%   mutate(site = plyr::mapvalues(site, from = c("Green Bay", "Illinois", "Madison"), to = c("Green Bay, WI", "Casey, IL", "Madison, WI"))) %>% 
  mutate(site = factor(site, levels = c("Green Bay, WI", "Madison, WI", "Casey, IL"))) %>% 
  mutate(treatment = factor(treatment, levels = c("snow reduction", "control"))) %>% 
  mutate(treatment = plyr::mapvalues(treatment, from = c("snow reduction", "control"), to = c("Snow Reduction", "Control"))) %>% 
  ggplot(aes(x = Time, y = TempC, group = treatment, color = treatment)) +
  geom_hline(yintercept = 0) +
  geom_vline(aes(xintercept = col_date), data = col_dates, color = "lightgray", linetype = "dashed", size = 1) +
  stat_summary(geom = "line", size = 0.5) +
  facet_wrap(~site, ncol = 1) +
  theme(text = element_text(size = 14)) +
  theme_classic() +
  theme(legend.position = "top",
        text = element_text(size = 14),
        legend.title=element_blank()) +
  ylab("Temperature (°C)") +
  xlab("Date") +
  scale_color_manual(values = c("royalblue2", "gold3"))

jpeg(filename = "temp_graph.jpeg", width = 10, height = 15, res = 400, unit = "cm")
temp_graph
dev.off()




#### plots ####
germ <- read.csv("germ_resp_comb.csv")

temp <- read.csv("temp_summary.csv") %>% 
  mutate(site = factor(site, levels = c("Illinois", "Madison", "Green Bay"))) %>% 
  mutate(site = plyr::mapvalues(site, from = c("Illinois", "Madison", "Green Bay"), to = c("Casey, IL", "Madison, WI", "Green Bay, WI"))) %>% 
  mutate(type = plyr::mapvalues(type, from = c("cold_days", "cold_strat_days", "warm_days"), to = c("Temp <0C", "Temp 0-5C", "Temp >5C"))) %>% 
  mutate(type = factor(type, levels = c("Temp <0C", "Temp 0-5C", "Temp >5C"))) %>% 
  mutate(treatment = plyr::mapvalues(treatment, from = c("control", "snow reduction"), to = c("Control", "Snow Reduction"))) %>% 
  mutate(col_date = date(col_date))

col_dates <- read.csv("col_dates.csv") %>% 
  mutate(col_date = paste(col_date, "00:00:00", sep = " ")) %>% 
  mutate(col_date = mdy_hms(col_date)) %>% 
  mutate(site = plyr::mapvalues(site, from = c("Green Bay", "Illinois", "Madison"), to = c("Green Bay, WI", "Casey, IL", "Madison, WI"))) %>% 
  mutate(site = factor(site, levels = c("Green Bay, WI", "Madison, WI", "Casey, IL")))

temp_plot <- temp  %>% 
  ggplot(aes(x = as.Date(col_date), y = mean_val, color = type, group = paste(type, treatment), linetype = treatment)) +
  geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),position = position_dodge(5)) +
  geom_line(position = position_dodge(5), linewidth = 1) +
  geom_point(position = position_dodge(5), size = 2) +
  facet_wrap(~site) +
  theme_classic() +
  scale_color_manual(values = c("dodgerblue4", "cyan3", "darkorange")) +
  labs(y = "Accumulated\nDays", x = "Date") +
  theme(legend.title = element_blank(),
        text = element_text(size = 12),
        strip.text = element_blank(),
        legend.position = "bottom")

germ_resp_comb_summary <- germ %>% 
  mutate(site = plyr::mapvalues(site, from = c("GB", "MA", "IL"), to = c("Green Bay, WI", "Madison, WI", "Casey, IL"))) %>% 
  mutate(site = factor(site, levels = c("Casey, IL", "Madison, WI",  "Green Bay, WI"))) %>% 
  left_join(col_dates) %>% 
  group_by(site, round, treatment, col_date, species) %>% 
  summarize(mean_germ = mean(perc_germ, na.rm = T),
            se_germ = sqrt(var(perc_germ, na.rm = T))/sqrt(n()),
            mean_time = mean(days_to_germ, na.rm = T),
            se_time = sqrt(var(days_to_germ, na.rm = T))/sqrt(n())) %>% 
  mutate(treatment = plyr::mapvalues(treatment, from = c("control", "snow reduction"), to = c("Control", "Snow Reduction"))) 

germ_plot <- germ_resp_comb_summary %>% 
  mutate(species = factor(species, levels = c("TO", "DC"))) %>% 
  #filter(species == "DC") %>% 
  mutate(mean_germ = ifelse(round < 5 | species != "DC", mean_germ, NA)) %>% 
  ggplot(aes(x = as.Date(col_date), y = mean_germ, group = paste(treatment, species), linetype = treatment, color = species)) +
  geom_errorbar(aes(ymin = mean_germ - se_germ, ymax = mean_germ + se_germ),position = position_dodge(5)) +
  geom_line(position = position_dodge(5), linewidth = 1) +
  geom_point(position = position_dodge(5), size=1.5) +
  facet_wrap(~site) +
  theme_classic() +
  #scale_color_manual(values = c("black", "darkgray")) +
  labs(y = "Germination\nProportion", x = NULL) +
  theme(legend.title = element_blank(),
        text = element_text(size = 12),
        legend.position = "none",
        axis.text.x = element_blank()) +
  scale_color_manual(values = c("black", "gray"))


DC_germ_plot <- germ_resp_comb_summary %>% mutate(treatment = plyr::mapvalues(treatment, from = c("control", "snow reduction"), to = c("Control", "Snow Reduction"))) %>% 
  filter(species == "DC") %>% 
  mutate(mean_germ = ifelse(round < 5, mean_germ, NA)) %>% 
  ggplot(aes(x = as.Date(col_date), y = mean_germ, group = paste(treatment), shape = treatment, linetype = treatment)) +
  geom_errorbar(aes(ymin = mean_germ - se_germ, ymax = mean_germ + se_germ),position = position_dodge(5)) +
  geom_line(position = position_dodge(5)) +
  geom_point(position = position_dodge(5)) +
  facet_wrap(~site) +
  theme_classic() +
  #scale_color_manual(values = c("black", "darkgray")) +
  labs(y = "Desmodium canadense\nGermination Proportion", x = NULL) +
  theme(legend.title = element_blank(),
        text = element_text(size = 12),
        legend.position = "none",
        axis.text.x = element_blank())

TO_germ_plot <- germ_resp_comb_summary %>% mutate(treatment = plyr::mapvalues(treatment, from = c("control", "snow reduction"), to = c("Control", "Snow Reduction"))) %>% 
  filter(species == "TO") %>% 
  #mutate(mean_germ = ifelse(round < 5, mean_germ, NA)) %>% 
  ggplot(aes(x = as.Date(col_date), y = mean_germ, group = paste(treatment), shape = treatment, linetype = treatment)) +
  geom_errorbar(aes(ymin = mean_germ - se_germ, ymax = mean_germ + se_germ),position = position_dodge(5)) +
  geom_line(position = position_dodge(5)) +
  geom_point(position = position_dodge(5)) +
  facet_wrap(~site) +
  theme_classic() +
  #scale_color_manual(values = c("black", "darkgray")) +
  labs(y = "Tradescantia Ohioense\nGermination Proportion", x = NULL) +
  theme(legend.title = element_blank(),
        text = element_text(size = 12),
        legend.position = "none",
        axis.text.x = element_blank())

timing_plot <- germ_resp_comb_summary %>% 
  mutate(mean_time = ifelse(round < 5 | species != "DC", mean_time, NA)) %>% 
  mutate(species = factor(species, levels = c("TO", "DC"))) %>% 
  ggplot(aes(x = as.Date(col_date), y = mean_time, group = paste(treatment, species), linetype = treatment, color = species)) +
  geom_errorbar(aes(ymin = mean_time - se_time, ymax = mean_time + se_time),position = position_dodge(5)) +
  geom_line(position = position_dodge(5), linewidth = 1) +
  geom_point(position = position_dodge(5), size = 1.5) +
  facet_wrap(~site) +
  theme_classic() +
  #scale_color_manual(values = c("coral3", "deepskyblue3")) +
  labs(y = "Days to First\nGermination", x = NULL) +
  theme(legend.title = element_blank(),
        text = element_text(size = 12),
        legend.position = "none",
        axis.text.x = element_blank())+
  scale_color_manual(values = c("black", "gray"))

germ_plot <- ggpubr::ggarrange(germ_plot, temp_plot, nrow = 2, heights = c(1, 1.2))

jpeg(filename = "germ_graph.jpeg", width = 20, height = 12, res = 400, unit = "cm")
germ_plot
dev.off()


ggpubr::ggarrange(timing_plot, temp_plot, nrow = 2, heights = c(1, 1.2))

#model result plot
germ <- read.csv("germ_resp_comb.csv")

temp <- read.csv("temp_summary.csv")

temp <- temp %>% 
  select(-se_val) %>% 
  pivot_wider(names_from = type, values_from = mean_val) %>% 
  mutate(site = plyr::mapvalues(site, from = c("Green Bay", "Illinois", "Madison"), to = c("GB", "IL", "MA")))

mod_dat <- germ %>% 
  left_join(temp) %>% 
  mutate(round = as.factor(round)) %>% 
  mutate(logit_germ = log((perc_germ+.1)/(100-(perc_germ+.1))))

####model test ####

TO_mod <- lmer(logit_germ ~ site * treatment + cold_days + cold_strat_days + warm_days + (1|siterep) + (1|round:site), data = mod_dat %>% filter(species == "TO", perc_germ < 100))


summary(TO_mod)
anova(TO_mod)
plot(allEffects(TO_mod))
ranova(TO_mod)

emmeans(TO_mod, pairwise~site)



DC_mod <- lmer(logit_germ ~ site * treatment + cold_days + cold_strat_days + warm_days + (1|siterep) + (1|round:site), data = mod_dat %>% filter(species == "DC", perc_germ < 100, round %in% c("1", "2", "3", "4")))

summary(DC_mod)
anova(DC_mod)
plot(allEffects(DC_mod))

emmeans(DC_mod, pairwise~site)
emtrends(DC_mod, pairwise~site, var = "round")



#### soil temp trends models ####
#cold days
cold_mod <- lmer(cold_days ~ site * treatment + (1|rep) + (1|round:site), data = mod_dat)

summary(cold_mod)
anova(cold_mod)
plot(allEffects(cold_mod))

emmeans(cold_mod, pairwise~treatment|site)

#cold strat days
cold_strat_mod <- lmer(cold_strat_days ~ site * treatment + (1|round:site) + (1|rep), data = mod_dat)

summary(cold_strat_mod)
anova(cold_strat_mod)
plot(allEffects(cold_strat_mod))

emmeans(cold_strat_mod, pairwise~treatment|site)

#### TO bivariate plots ####

TO_trends_cold <- as.data.frame(allEffects(TO_mod, xlevels = list(cold_days = seq(0,75, 5)))[[1]]) %>% mutate(type = "cold_days") %>% rename("days" = cold_days)
TO_trends_cold_s <- as.data.frame(allEffects(TO_mod, xlevels = list(cold_strat_days = seq(5,65, 5)))[[2]]) %>% mutate(type = "cold_strat_days") %>% rename("days" = cold_strat_days)
TO_trends_warm <- as.data.frame(allEffects(TO_mod, xlevels = list(warm_days = seq(0,80, 5)))[[3]]) %>% mutate(type = "warm_days") %>% rename("days" = warm_days)

TO_trends <- bind_rows(TO_trends_cold, TO_trends_cold_s, TO_trends_warm) %>% mutate(type = plyr::mapvalues(type, from = c("cold_days", "cold_strat_days", "warm_days"), to = c("Temp <0C", "Temp 0-5C", "Temp >5C"))) %>%   mutate(type = factor(type, levels = c("Temp <0C", "Temp 0-5C", "Temp >5C")))

line_plot <- mod_dat %>% filter(species == "TO", perc_germ < 100) %>% pivot_longer(cols = c(cold_days, cold_strat_days, warm_days), names_to = "type", values_to = "days") %>% mutate(type = plyr::mapvalues(type, from = c("cold_days", "cold_strat_days", "warm_days"), to = c("Temp <0C", "Temp 0-5C", "Temp >5C"))) %>%   mutate(type = factor(type, levels = c("Temp <0C", "Temp 0-5C", "Temp >5C"))) %>% mutate(site = factor(site, levels = c("IL", "MA", "GB"))) %>% 
  mutate(site = plyr::mapvalues(site, from = c("IL", "MA", "GB"), to = c("Casey, IL", "Madison, WI", "Green Bay, WI"))) %>% rename("Site"=site) %>% 
  ggplot() +
  geom_point(aes(x = days, y = inv.logit(logit_germ), shape = Site), size = 1.5) +
  #stat_smooth(aes(x = days, y = inv.logit(logit_germ), group = Site), method = "glm", formula = "y ~ x") +
  geom_ribbon(aes(x = days, ymin = inv.logit(lower), ymax = inv.logit(upper)), data = TO_trends, alpha = 0.3) +
  geom_line(aes(x = days, y = inv.logit(fit), color = type), data = TO_trends, linewidth = 2) +
  facet_wrap(~type, scales = "free", strip.position = "bottom") +
  theme_classic() +
  scale_color_manual(values = c("dodgerblue4", "cyan3", "darkorange"), guide = "none") +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        text = element_text(size = 12),
        legend.position = "top") +
  labs(x = "Days", y = "Germination Proportion")

jpeg(filename = "line_graph.jpeg", width = 20, height = 8, res = 400, unit = "cm")
line_plot
dev.off()


