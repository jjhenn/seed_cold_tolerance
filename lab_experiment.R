# Seed Cold Tolerance
# LML
# Summer 2021

# Enter and organize data

library(tidyverse)
library(tidyr)
library(lmerTest)
library(broom)
library(gridExtra)
library(grid)

#Import data
mod_dat <- read.csv("seed_cold_tol_dat1.csv") %>% 
  mutate(trt_temp_planned = as.factor(trt_temp_planned))

trait_summary <- mod_dat %>% 
  select(spher2, sav, sc, perc, perc_cont, spp, trt_temp_planned) %>% 
  pivot_longer(cols = -c(spp, trt_temp_planned), names_to = "trait", values_to = "vals") %>% 
  group_by(spp, trt_temp_planned, trait) %>% 
  summarize(mean_val = mean(vals, na.rm = T), se_val = sd(vals, na.rm = T)/sqrt(n()), n = n())

#PCA of traits
pca_dat <- mod_dat %>% 
  select(mg_per_seed, sc, spher2, sav) %>% 
  distinct() %>% 
  na.omit()

pca <- prcomp(pca_dat, center = T, scale = T)
biplot(pca)

pc_points <- as.data.frame(pca$x)
pc_vectors <- as.data.frame(pca$rotation) 
row.names(pc_vectors) <- c("Seed Mass", "Seed Coat\nThickness", "Sphericity", "Surface Area:\nVolume")

trait_pca <- pc_points %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_segment(aes(x = 0, y= 0, xend = PC1 * 2, yend = PC2 * 2), data = pc_vectors, arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text(aes(x = PC1 * 3, y = PC2 * 2.5, label = row.names(pc_vectors)), data= pc_vectors) +
  theme_classic() +
  labs(x = "PC1 (55%)", y = "PC2 (28%)")


jpeg("trait_pca.jpeg", width = 5, height = 5, units = "in", res = 300)
trait_pca
dev.off()

#make models for each trait
mod_dat_long_cont <- mod_dat %>% filter(family != "Apiaceae") %>% 
  select(spp, rep, trt_temp_planned, rratio, "Seed Mass (mg)" = mg_per_seed , family, "Seed Coat Thickness (mm)" = sc , "Sphericity" = spher2, "Surface Area:Volume" = sav, coldStrat) %>% 
  pivot_longer(cols = -c(spp, rep, trt_temp_planned, rratio, family, coldStrat), names_to = "trait", values_to = "values") %>% 
  filter(!rep %in% c("1b", "2b", "3b")) %>% 
  mutate(type = paste(rep, trt_temp_planned))

mod_dat_long_cat <- mod_dat %>% filter(family != "Apiaceae") %>% 
  select(spp, rep, trt_temp_planned, rratio, "Family" = family,  "Cold Stratification" = coldStrat, "Maturation Timing (month)" = maturation) %>% mutate(`Cold Stratification` = as.factor(`Cold Stratification`)) %>% mutate(`Maturation Timing (month)` = factor(`Maturation Timing (month)`, levels = c("7", "8", "9", "10"))) %>% 
  pivot_longer(cols = -c(spp, rep, trt_temp_planned, rratio), names_to = "trait", values_to = "values") %>% 
  na.omit() 

#loop to fit models
anova_out_comb <- data.frame(NULL)
plot_out_comb <- data.frame(NULL)
for(i in unique(mod_dat_long_cont$trait)){
  fit <- lmer(rratio ~ trt_temp_planned*values + (1|spp:rep), data = mod_dat_long_cont %>% filter(trait == i))
  
  anova_out <- as.data.frame(anova(fit)) %>% mutate(trait = i)
  plot_out <- as.data.frame(allEffects(fit, xlevels = list(values = seq(min(mod_dat_long_cont$values[mod_dat_long_cont$trait == i], na.rm = T), max(mod_dat_long_cont$values[mod_dat_long_cont$trait == i], na.rm = T), 0.01)))[[1]]) %>% mutate(trait = i)
  
  anova_out_comb <- bind_rows(anova_out_comb, anova_out)
  plot_out_comb <- bind_rows(plot_out_comb, plot_out)
}


anova_text <- anova_out_comb %>% 
  mutate(symbol = ifelse(`Pr(>F)` < 0.01, "**", ifelse(`Pr(>F)` < 0.05, "*", "NS"))) %>% 
  mutate(type = rep(c("Treatment", "Trait", "Trait:Treatment"), 4)) %>% 
  mutate(comb = paste(type, symbol, sep = ":")) %>% 
  select(comb, trait, type) %>% 
  pivot_wider(names_from = type, values_from = comb) %>% 
  mutate(text = paste(Treatment, Trait, `Trait:Treatment`, sep = ", "))%>% 
  mutate(letter = plyr::mapvalues(trait, from = c("Seed Mass (mg)", "Seed Coat Thickness (mm)", "Sphericity", "Surface Area:Volume"), to = c("B.", "A.", "C.", "D.")))

plot_cont <- plot_out_comb %>% 
  ggplot(aes(x = values, y = fit)) +
  geom_ribbon(aes(x = values, ymin = lower, ymax = upper, group = trt_temp_planned), alpha = 0.1) +
  geom_point(aes(x = values, y = rratio, color = trt_temp_planned), data = mod_dat_long_cont, alpha = 0.2) +
  geom_line(aes(color = trt_temp_planned), linewidth = 1.5) +
  geom_text(aes(x = Inf, y = 2.5, label = text, hjust = 1), data = anova_text %>% filter(trait != "Family"), size = 2.5) +
  geom_text(aes(x = -Inf, y = 2.5, label = letter, hjust = -0.5), data = anova_text %>% filter(trait != "Family"), size = 4) +
  facet_wrap(~trait, scales = "free", strip.position = "bottom") +
  theme_classic() +
  theme(text = element_text(size = 12),
        strip.placement = "outside",
        strip.background = element_blank()) +
  scale_color_manual(values = c("royalblue3", "lightblue3"), name = "Treatment") +
  labs(y = " ", x = NULL)

#categorical traits
anova_out_comb <- data.frame(NULL)
plot_out_comb <- data.frame(NULL)
for(i in unique(mod_dat_long_cat$trait)){
  fit <- lmer(rratio ~ trt_temp_planned*values + (1|spp:rep), data = mod_dat_long_cat %>% filter(trait == i))
  
  anova_out <- as.data.frame(anova(fit)) %>% mutate(trait = i)
  plot_out <- as.data.frame(allEffects(fit)[[1]]) %>% mutate(trait = i)
  
  anova_out_comb <- bind_rows(anova_out_comb, anova_out)
  plot_out_comb <- bind_rows(plot_out_comb, plot_out)
}

anova_text <- anova_out_comb %>% 
  mutate(symbol = ifelse(`Pr(>F)` < 0.01, "**", ifelse(`Pr(>F)` < 0.05, "*", "NS"))) %>% 
  mutate(type = rep(c("Treatment", "Trait", "Trait:Treatment"), 3)) %>% 
  mutate(comb = paste(type, symbol, sep = ":")) %>% 
  select(comb, trait, type) %>% 
  pivot_wider(names_from = type, values_from = comb) %>% 
  mutate(text = paste(Treatment, Trait, `Trait:Treatment`, sep = ", ")) %>% 
  mutate(letter = ifelse(trait == "Cold Stratification", "E.", ifelse(trait == "Family", "G.", "F.")))

plot_cat <- plot_out_comb %>% filter(trait != "Family") %>% 
  ggplot(aes(x = values, y = fit)) +
  geom_jitter(aes(x = values, y = rratio, color = trt_temp_planned), data = mod_dat_long_cat  %>% filter(trait != "Family"), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), alpha = 0.1) +
  geom_errorbar(aes(x = values, ymin = lower, ymax = upper, group = trt_temp_planned), position = position_dodge(0.5), width = 0, linewidth = 1) +
  geom_point(aes(color = trt_temp_planned), size = 2, position = position_dodge(0.5)) +
  geom_text(aes(x = Inf, y = 2.5, label = text, hjust = 1), data = anova_text %>% filter(trait != "Family"), size = 2.5) +
  geom_text(aes(x = -Inf, y = 2.5, label = letter, hjust = -0.5), data = anova_text %>% filter(trait != "Family"), size = 4) +
  facet_wrap(~trait, scales = "free", strip.position = "bottom", ncol = 2) +
  theme_classic() +
  theme(text = element_text(size = 12),
        strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 0)) +
  scale_color_manual(values = c("royalblue3", "lightblue3"), name = "Treatment") +
  labs(y = "Log Response Ratio", x = NULL)

plot_family <- plot_out_comb %>% filter(trait == "Family") %>% 
  ggplot(aes(x = values, y = fit)) +
  geom_jitter(aes(x = values, y = rratio, color = trt_temp_planned), data = mod_dat_long_cat  %>% filter(trait == "Family"), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), alpha = 0.1) +
  geom_errorbar(aes(x = values, ymin = lower, ymax = upper, group = trt_temp_planned), position = position_dodge(0.5), width = 0, linewidth = 1) +
  geom_point(aes(color = trt_temp_planned), size = 2, position = position_dodge(0.5)) +
  geom_text(aes(x = Inf, y = 2.5, label = text, hjust = 1), data = anova_text %>% filter(trait == "Family"), size = 2.5) +
  geom_text(aes(x = -Inf, y = 2.5, label = letter, hjust = -0.5), data = anova_text %>% filter(trait == "Family"), size = 4) +
  facet_wrap(~trait, scales = "free", strip.position = "bottom", ncol = 2) +
  theme_classic() +
  theme(text = element_text(size = 12),
        strip.placement = "outside",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0)) +
  scale_color_manual(values = c("royalblue3", "lightblue3"), name = "Treatment") +
  labs(y = " ", x = NULL) +
  scale_y_continuous(breaks = c(-6, -4.5, -3, -1.5, 0, 1.5, 3))

trait_fig <- ggpubr::ggarrange(plot_cont, plot_cat, plot_family, ncol = 1, common.legend = T, heights = c(2.5,1.2,2))

jpeg("trait_fig.jpeg", width = 5.6, height = 8.5, units = "in", res = 300)
trait_fig
dev.off()

