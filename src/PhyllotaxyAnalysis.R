options(contrasts=c('contr.sum', 'contr.poly'))
library(lme4)
library(tidyverse)
library(readxl)
library(cowplot)
library(MoMAColors)
library(BiocIO) # read GFF file
library(GenomicRanges) 
library(Gviz)
library(ggrepel)
library(scales)
library(viridis)
library(stringi)
library(car)
# paths assume wd = src/
colors <- rep(moma.colors('VanGogh', 5), 2)
colors2 <- moma.colors('VanGogh', 2)

normalize360 <- function(x)
{
  return (x %% 360)
}

normalize180 <- function(x)
{
  x <- normalize360(x)
  return (abs(x - 180))
}

calcHeritability <- function(data, trait)
{
  # Filter out NAs for this trait and the reference genotype
  dfh <- filter(data, !is.na({{trait}}), PI_num!='PI_656058')
  dfh <- group_by(dfh, PI_num, JS_ID)
  dfh <- summarize(dfh, n = n(), trait = {{trait}})
  dfh <- ungroup(dfh)
  dfh <- filter(dfh, n>1)
  print(length(unique(dfh$PI_num)))
  fit <- lmer(trait ~ (1|PI_num), data = dfh)
  vc <- as.data.frame(formatVC(VarCorr(fit)))
  sigma <- as.numeric(vc$Std.Dev.)
  h <- sigma[1]^2/(sigma[1]^2 + sigma[2]^2)
  resids <- residuals(fit)
  print(plot(resids))
  qqPlot(resids, main = deparse(substitute(trait)))
  return (h)
}

# Function not requiring JS_ID column for dataframes already reduced to 
calcHeritability2 <- function(data, trait)
{
  # Filter out NAs for this trait and the reference genotype
  dfh <- filter(data, !is.na({{trait}}), PI_num!='PI_656058')
  # Randomly add back in 2 reps of ref genotype
  dfh_ref <- filter(data, !is.na({{trait}}), PI_num =='PI_656058')
  r <- sample(1:8, size = 2, replace = FALSE)
  dfh <- add_row(dfh, dfh_ref[(r[1]), ])
  dfh <- add_row(dfh, dfh_ref[(r[2]), ])
  dfh <- group_by(dfh, PI_num)
  dfh <- summarize(dfh, n = n(), trait = {{trait}})
  dfh <- ungroup(dfh)
  dfh <- filter(dfh, n>1)
  print(length(unique(dfh$PI_num)))
  fit <- lmer(trait ~ (1|PI_num), data = dfh)
  vc <- as.data.frame(formatVC(VarCorr(fit)))
  sigma <- as.numeric(vc$Std.Dev.)
  h <- sigma[1]^2/(sigma[1]^2 + (sigma[2]^2)/2)
  resids <- residuals(fit)
  print(plot(resids))
  qqPlot(resids, main = deparse(substitute(trait)))
  
  hist <- ggplot(data, aes({{trait}})) + 
    geom_histogram(binwidth = 2) + 
    theme_minimal()
  print(hist)
  
  return (h)
}

getSupport <- function(path, trait)
{
  files <- list.files(path = path, pattern = {trait}, full.names = TRUE)
  
  if (length(files) >= 1)
  {
    currFile <- files[1]
    signals <- read.csv(currFile, header = FALSE, skip = 1, col.names = c('SNP', 'CHROM', 'POS', 'REF', 'ALT', 'EFFECT', 'SE', 'P')) %>% 
      as_tibble() %>%
      rowwise() %>%
      mutate(CHROM = as.numeric(str_remove(CHROM, 'Chr'))) %>%
      select(c(SNP:POS, EFFECT:P))
    # print(tbl_sum(signals))
    
    for (i in 2:length(files))
    {
      currFile <- files[i]
      currFile <- read.csv(currFile, header = FALSE, skip = 1, col.names = c('SNP', 'CHROM', 'POS', 'REF', 'ALT', 'EFFECT', 'SE', 'P')) %>%
        as_tibble() %>%
        rowwise() %>%
        mutate(CHROM = as.numeric(str_remove(CHROM, 'Chr'))) %>%
        select(c(SNP:POS, EFFECT:P))
      # print(tbl_sum(currFile))
      signals <- bind_rows(signals, currFile)
      # print(tbl_sum(signals))
    }
    
    # print(length(unique(signals$SNP)))
    
    signals <- signals %>%
      group_by(SNP, CHROM, POS) %>%
      summarise(EFFECT = mean(EFFECT), SE = mean(SE), P = mean(P), {{trait}} := n()/100)
    # print(tbl_sum(signals))
    # write.csv(signals, file = paste0(path, 'Z', trait, 'signals.csv'), quote = FALSE, row.names = FALSE, col.names = TRUE)
    
    return(signals)
  }
  else 
  {
    print(paste0('File not found.', trait))
  }
}

plotRMIPManhattan <- function(data, trait, title=NULL, colorVec = colors)
{
  t <- 'a'
  if (is.null(title))
  {
    t <- deparse(substitute(trait))
  }
  else
  {
    t <- title
  }
  data_cum <- data %>%
    group_by(.data[['CHROM']]) %>%
    summarise(max_bp = max(POS)) %>%
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
    select(c(CHROM, bp_add))
  df <- data
  df <- inner_join(df, data_cum, by = 'CHROM') %>%
    mutate(loc = (POS + bp_add))
  axis_set <- df %>%
    ungroup()%>%
    summarise(center = mean(loc), .by = 'CHROM')
  
  manhattan <- ggplot(df, aes(loc, {{trait}}, color = factor(CHROM))) +
    geom_point(size = 3, alpha = 1) +
    geom_hline(yintercept = 0.1, linetype = 2, color = 'black') +
    scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
    scale_y_continuous(limits = c(0, 0.3)) +
    scale_size_continuous(range = c(1, 8)) +
    scale_color_manual(values = colorVec) +
    ylab('RMIP') +
    xlab(NULL + unique(axis_set$CHROM)) +
    labs(title = t) +
    labs(x = 'Chromosome') +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0), 
                                     vjust = 0.5, hjust = 0.5),
          axis.text.y = element_text(size = 11, color = 'black', vjust = 0, hjust = 0.5),
          legend.text = element_text(size = 11, color = 'black'),
          plot.title = element_text(size = 11, color = 'black', vjust = 0, hjust = 0.5),
          text = element_text(size = 11, color = 'black'),
          legend.position = 'none',
          line = element_line(color = 'black', linewidth = 1),
          panel.grid = element_blank())  
  print(manhattan)
  return(manhattan)
}

phiData1 <- read.csv('../Data/processedData.csv', header = TRUE)
outliers_remove <- read.csv('../Data/outliers_remove.csv', header = TRUE)
outliers_remove$plant.date <- stri_join(outliers_remove$plant_num, outliers_remove$img_date, sep = '.')
outliers_keep <- read.csv('../Data/outliers_keep.csv', header = TRUE)
outliers_keep <- unique(outliers_keep$plant_num)

# First row is empty so remove it
phiData1 <- phiData1[2:length(phiData1$plant_name),]
phiData1 <- filter(phiData1, voxel_len == '6cm')
phiData1$plant.date <- stri_join(phiData1$plant_num, phiData1$img_date, sep = '.')

reconstructionReliability <- phiData1 %>%
  select(c(plant_num, PI_num, img_date, phi_0, phi_1, phi_2, phi_3, phi_4, phi_5)) %>%
  rowwise() %>%
  mutate(varphi1 = normalize360(normalize360(phi_1) - normalize360(phi_0)),
         varphi2 = normalize360(normalize360(phi_2) - normalize360(phi_1)),
         varphi3 = normalize360(normalize360(phi_3) - normalize360(phi_2)),
         varphi4 = normalize360(normalize360(phi_4) - normalize360(phi_3)),
         varphi5 = normalize360(normalize360(phi_5) - normalize360(phi_4))) %>%
  pivot_longer(starts_with('varphi'), 
               names_to = 'varphiNum',
               values_to = 'val', 
               names_prefix = 'varphi') %>%
  pivot_wider(id_cols = c(plant_num, varphiNum), names_from = img_date, values_from = val)
colnames(reconstructionReliability) <- c("plant_num", "varphiNum", 'day1', 'day2', 'day3')

# Check reliability if we don't remove outliers
reconstructionReliabilityPlantLvl <- reconstructionReliability %>%
  group_by(plant_num) %>%
  summarise(cor1.2 = cor(day1, day2, use = 'na.or.complete'),
            cor1.3 = cor(day1, day3, use = 'na.or.complete'),
            cor2.3 = cor(day2, day3, use = 'na.or.complete'))
conjugate1.3 <- filter(reconstructionReliabilityPlantLvl, cor1.3 < 0)
conjugate1.3 <- conjugate1.3$plant_num

conjugate1.2 <- filter(reconstructionReliabilityPlantLvl, cor1.2 < 0)
conjugate1.2 <- conjugate1.2$plant_num

conjugate2.3 <- filter(reconstructionReliabilityPlantLvl, cor2.3 < 0)
conjugate2.3 <- conjugate2.3$plant_num

reconstructionReliability <- reconstructionReliability %>%
  rowwise() %>%
  mutate(day3.1 = case_when(plant_num %in% conjugate1.3 ~ 360 - day3, 
                            .default = day3), 
         day2.1 = case_when(plant_num %in% conjugate1.2 ~ 360 - day2,
                            .default = day2),
         day3.2 = case_when(plant_num %in% conjugate2.3 ~ 360 - day3, 
                            .default = day3))

cor(reconstructionReliability$day1, reconstructionReliability$day3.1, use = 'na.or.complete')^2
cor(reconstructionReliability$day1, reconstructionReliability$day2.1, use = 'na.or.complete')^2
cor(reconstructionReliability$day2, reconstructionReliability$day3.2, use = 'na.or.complete')^2

phe_orig <- phiData1 %>%
  rowwise() %>%
  select(c(plant_num, PI_num, img_date, phi_0, phi_1, phi_2, phi_3, phi_4)) %>%
  rowwise() %>%
  mutate(varphi1 = normalize360(normalize360(phi_1) - normalize360(phi_0)),
         varphi2 = normalize360(normalize360(phi_2) - normalize360(phi_1)),
         varphi3 = normalize360(normalize360(phi_3) - normalize360(phi_2)),
         varphi4 = normalize360(normalize360(phi_4) - normalize360(phi_3))) %>%
  pivot_longer(starts_with('varphi'), 
               names_to = 'varphiNum',
               values_to = 'val', 
               names_prefix = 'varphi') %>% 
  select(plant_num, PI_num, varphiNum, val) %>%
  mutate(val_normalized = normalize180(val)) %>%
  group_by(plant_num) %>%
  summarise(med = median(val_normalized, na.rm = TRUE))

# Remove known outliers due to tillers, missing leaves, lodging
phiData1 <- filter(phiData1, !(plant.date %in% outliers_remove$plant.date))

# Calculate normalized equivalent of each angle (in range of 0 - 359)
# R modulo operation always gives positive result (not true in every language)
phiData1$phi_0Norm <- normalize360(phiData1$phi_0)
phiData1$phi_1Norm <- normalize360(phiData1$phi_1)
phiData1$phi_2Norm <- normalize360(phiData1$phi_2)
phiData1$phi_3Norm <- normalize360(phiData1$phi_3)
phiData1$phi_4Norm <- normalize360(phiData1$phi_4)
phiData1$phi_5Norm <- normalize360(phiData1$phi_5)
phiData1$phi_6Norm <- normalize360(phiData1$phi_6)
phiData1$phi_7Norm <- normalize360(phiData1$phi_7)
phiData1$phi_8Norm <- normalize360(phiData1$phi_8)
phiData1$phi_9Norm <- normalize360(phiData1$phi_9)
phiData1$phi_10Norm <- normalize360(phiData1$phi_10)
phiData1$phi_11Norm <- normalize360(phiData1$phi_11)
phiData1$phi_12Norm <- normalize360(phiData1$phi_12)
phiData1$phi_13Norm <- normalize360(phiData1$phi_13)
phiData1$phi_14Norm <- normalize360(phiData1$phi_14)

phi6.wide <- phiData1 %>%
  rowwise() %>%
  mutate(varphi1 = normalize360(phi_1Norm - phi_0Norm),
         varphi2 = normalize360(phi_2Norm - phi_1Norm),
         varphi3 = normalize360(phi_3Norm - phi_2Norm),
         varphi4 = normalize360(phi_4Norm - phi_3Norm),
         varphi5 = normalize360(phi_5Norm - phi_4Norm)) %>%
  select(c(plant_num, img_date, contains('varphi'))) %>%
  pivot_longer(contains('varphi'), names_to = 'varphiNum', values_to = 'val') %>%
  pivot_wider(id_cols = c(plant_num, varphiNum), names_from = img_date, values_from = val)
colnames(phi6.wide) <- c('plant_num', 'varphiNum', 'day1', 'day2', 'day3')

# Check reconstruction reliability
phi6.pl <- phi6.wide %>%
  group_by(plant_num) %>%
  summarise(cor1.3 = cor(day1, day3, use = 'na.or.complete'), 
            cor1.2 = cor(day1, day2, use = 'na.or.complete'),
            cor2.3 = cor(day2, day3, use = 'na.or.complete'))

conjugate1.3 <- filter(phi6.pl, cor1.3 < 0)
conjugate1.3 <- conjugate1.3$plant_num

conjugate1.2 <- filter(phi6.pl, cor1.2 < 0)
conjugate1.2 <- conjugate1.2$plant_num

conjugate2.3 <- filter(phi6.pl, cor2.3 < 0)
conjugate2.3 <- conjugate2.3$plant_num

phi6.wide <- phi6.wide %>%
  rowwise() %>%
  mutate(day3.1 = case_when(plant_num %in% conjugate1.3 ~ 360 - day3, .default = day3),
         day2.1 = case_when(plant_num %in% conjugate1.2 ~ 360 - day2, .default = day2),
         day3.2 = case_when(plant_num %in% conjugate2.3 ~ 360 - day3, .default = day3))

cor(phi6.wide$day1, phi6.wide$day3.1, use = 'na.or.complete')^2
cor(phi6.wide$day1, phi6.wide$day2.1, use = 'na.or.complete')^2
cor(phi6.wide$day2, phi6.wide$day3.2, use = 'na.or.complete')^2

phi6.wide <- phiData1 %>%
  rowwise() %>%
  mutate(varphi1 = normalize360(phi_1Norm - phi_0Norm),
         varphi2 = normalize360(phi_2Norm - phi_1Norm),
         varphi3 = normalize360(phi_3Norm - phi_2Norm),
         varphi4 = normalize360(phi_4Norm - phi_3Norm),
         varphi5 = normalize360(phi_5Norm - phi_4Norm)) %>%
  select(c(plant_num, img_date, contains('varphi'))) %>%
  pivot_longer(contains('varphi'), names_to = 'varphiNum', values_to = 'val') %>%
  filter(val >= 90 & val <= 270) %>%
  pivot_wider(id_cols = c(plant_num, varphiNum), names_from = img_date, values_from = val)
colnames(phi6.wide) <- c("plant_num", "varphiNum", 'day1', 'day2', 'day3')

phi6.pl <- phi6.wide %>%
  group_by(plant_num) %>%
  summarise(cor1.3 = cor(day1, day3, use = 'na.or.complete'), 
            cor1.2 = cor(day1, day2, use = 'na.or.complete'),
            cor2.3 = cor(day2, day3, use = 'na.or.complete'))

conjugate1.3 <- filter(phi6.pl, cor1.3 < 0)
conjugate1.3 <- conjugate1.3$plant_num

conjugate1.2 <- filter(phi6.pl, cor1.2 < 0)
conjugate1.2 <- conjugate1.2$plant_num

conjugate2.3 <- filter(phi6.pl, cor2.3 < 0)
conjugate2.3 <- conjugate2.3$plant_num

phi6.wide <- phi6.wide %>%
  rowwise() %>%
  mutate(day3.1 = case_when(plant_num %in% conjugate1.3 ~ 360 - day3, .default = day3),
         day2.1 = case_when(plant_num %in% conjugate1.2 ~ 360 - day2, .default = day2),
         day3.2 = case_when(plant_num %in% conjugate2.3 ~ 360 - day3, .default = day3))

cor(phi6.wide$day1, phi6.wide$day3.1, use = 'na.or.complete')^2
r2Day1.2 <- round(cor(phi6.wide$day1, phi6.wide$day2.1, use = 'na.or.complete')^2, digits = 2)
cor(phi6.wide$day2, phi6.wide$day3.2, use = 'na.or.complete')^2

reconstructionReliabilityPlot <- ggplot(phi6.wide, aes(day1, day2.1)) +
  geom_point(color = colors2[2]) +
  scale_x_continuous(breaks = c(90, 135, 180, 225, 270),
                     labels = c(expression('90'*degree), expression('135'*degree),
                                expression('180'*degree), expression('225'*degree),
                                expression('270'*degree)),
                     limits = c(90, 270)) +
  scale_y_continuous(breaks = c(90, 135, 180, 225, 270),
                     labels = c(expression('90'*degree), expression('135'*degree),
                                expression('180'*degree), expression('225'*degree),
                                expression('270'*degree)),
                     limits = c(90, 270)) +
  labs(x = 'Timepoint 1 Reconstruction',
       y = 'Timepoint 2 Reconstruction',
       title = bquote('R'^2 == .(r2Day1.2))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0),
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 11, color = 'black', vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11, color = 'black'),
        plot.title = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0),
                                  vjust = 0.5, hjust = 0.5),
        text = element_text(size = 11, color = 'black'),
        legend.position = 'bottom',
        line = element_line(color = 'black', linewidth = 1),
        panel.grid = element_blank())
reconstructionReliabilityPlot

r2Day1.3 <- round(cor(phi6.wide$day1, phi6.wide$day3.1, use = 'complete.obs')^2, digits = 2)
p1.3 <- ggplot(phi6.wide, aes(day1, day3.1)) + 
  geom_point(color = colors2[2]) + 
  scale_x_continuous(breaks = c(90, 135, 180, 225, 270), 
                     labels = c(expression('90'*degree), expression('135'*degree),
                                expression('180'*degree), expression('225'*degree), 
                                expression('270'*degree)), 
                     limits = c(90, 270)) + 
  scale_y_continuous(breaks = c(90, 135, 180, 225, 270), 
                     labels = c(expression('90'*degree), expression('135'*degree),
                                expression('180'*degree), expression('225'*degree), 
                                expression('270'*degree)),
                     limits = c(90, 270)) + 
  labs(x = 'Timepoint 1 Reconstruction', 
       y = 'Timepoint 3 Reconstruction', 
       title = bquote('R'^2 == .(r2Day1.3))) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0), 
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 11, color = 'black', vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11, color = 'black'),
        plot.title = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0), 
                                  vjust = 0.5, hjust = 0.5),
        text = element_text(size = 11, color = 'black'),
        legend.position = 'bottom',
        line = element_line(color = 'black', linewidth = 1),
        panel.grid = element_blank())
p1.3

r2Day2.3 <- round(cor(phi6.wide$day2, phi6.wide$day3.2, use = 'complete.obs')^2, digits = 2)
p2.3 <- ggplot(phi6.wide, aes(day2, day3.2)) + 
  geom_point(color = colors2[2]) + 
  scale_x_continuous(breaks = c(90, 135, 180, 225, 270), 
                     labels = c(expression('90'*degree), expression('135'*degree),
                                expression('180'*degree), expression('225'*degree), 
                                expression('270'*degree)), 
                     limits = c(90, 270)) + 
  scale_y_continuous(breaks = c(90, 135, 180, 225, 270), 
                     labels = c(expression('90'*degree), expression('135'*degree),
                                expression('180'*degree), expression('225'*degree), 
                                expression('270'*degree)),
                     limits = c(90, 270)) + 
  labs(x = 'Timepoint 2 Reconstruction', 
       y = 'Timepoint 3 Reconstruction', 
       title = bquote('R'^2 == .(r2Day2.3))) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0), 
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 11, color = 'black', vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11, color = 'black'),
        plot.title = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0), 
                                  vjust = 0.5, hjust = 0.5),
        text = element_text(size = 11, color = 'black'),
        legend.position = 'bottom',
        line = element_line(color = 'black', linewidth = 1),
        panel.grid = element_blank())
p2.3

figs1 <- plot_grid(p2.3, p1.3, nrow = 1, labels = 'AUTO')
figs1
# ggsave('../Images/figS1.png', plot = figs1, width = 6.5, height = 3.5, units = 'in', dpi = 1000, bg = 'white')

# Normalize to the absolute value of the difference from 180 degrees: the choice between two conjugate angles depends only on the direction of the measurement, which likely has little biological relevance
phiData1$dif_phi1_phi0 <- normalize180(phiData1$phi_1Norm - phiData1$phi_0Norm)  
phiData1$dif_phi2_phi1 <- normalize180(phiData1$phi_2Norm - phiData1$phi_1Norm)  
phiData1$dif_phi3_phi2 <- normalize180(phiData1$phi_3Norm - phiData1$phi_2Norm)  
phiData1$dif_phi4_phi3 <- normalize180(phiData1$phi_4Norm - phiData1$phi_3Norm)  
phiData1$dif_phi5_phi4 <- normalize180(phiData1$phi_5Norm - phiData1$phi_4Norm)  
phiData1$dif_phi6_phi5 <- normalize180(phiData1$phi_6Norm - phiData1$phi_5Norm)  
phiData1$dif_phi7_phi6 <- normalize180(phiData1$phi_7Norm - phiData1$phi_6Norm)  
phiData1$dif_phi8_phi7 <- normalize180(phiData1$phi_8Norm - phiData1$phi_7Norm)  
phiData1$dif_phi9_phi8 <- normalize180(phiData1$phi_9Norm - phiData1$phi_8Norm) 
phiData1$dif_phi10_phi9 <- normalize180(phiData1$phi_10Norm - phiData1$phi_9Norm)
phiData1$dif_phi11_phi10 <- normalize180(phiData1$phi_11Norm - phiData1$phi_10Norm)
phiData1$dif_phi12_phi11 <- normalize180(phiData1$phi_12Norm - phiData1$phi_11Norm)
phiData1$dif_phi13_phi12 <- normalize180(phiData1$phi_13Norm - phiData1$phi_12Norm)
phiData1$dif_phi14_phi13 <- normalize180(phiData1$phi_14Norm - phiData1$phi_13Norm)

# Pivot 
phiData1_long <- pivot_longer(phiData1, cols = starts_with('dif'), names_to = 'angle', names_prefix = 'dif_', 
                              values_to = 'deg', values_drop_na = TRUE)

# Filter to first 4 leaves
phi4 <- filter(phiData1_long, angle %in% c('phi1_phi0', 'phi2_phi1', 'phi3_phi2', 'phi4_phi3'))
# # Plot distribution of angles in first 4 leaves
# hist4 <- ggplot(data = phi4, aes(deg, fill = 'red')) +
#   geom_histogram(binwidth = 2) +
#   theme_minimal() +
#   theme(panel.background = element_blank(), legend.position = 'none', 
#         text = element_text(family  = 'Times New Roman', size = 16, colour = 'black', ), axis.line =
#           element_line(color = 'black')) + 
#   labs(x = 'Difference of Angles Between Sequential Leaves in Lower Canopy from 180 Degrees', y = 'Frequency')
# print(hist4)
# # Calculate quantiles
# quant <- quantile(phi4$deg, probs = c(0.25,0.75))
# iqr <- (quant[2] - quant[1])
# t1.5 <- 1.5*iqr
# t1.75 <- 1.75*iqr
# t2 <- 2*iqr
# 
# hist4t <- ggplot(data = phi4, aes(deg)) +
#   geom_histogram(binwidth = 5) +
#   geom_vline(aes(xintercept = t1.5)) + 
#   geom_vline(aes(xintercept = t1.75)) +
#   geom_vline(aes(xintercept = t2)) +
#   theme_minimal()
# print(hist4t)
# 
# #Remove obs w/ deg > 100 unless the plant_num is in outliers_keep
# phi4_100 <- filter(phi4, (deg <= 100)|(plant_num %in% outliers_keep))
# 
# hist4_100 <- ggplot(data = phi4_100, aes(deg)) +
#   geom_histogram(binwidth = 5) +
#   theme_minimal()
# print(hist4_100)

# Remove obs w/ deg > 90. Strict cutoff. 
phi4_90t <- filter(phi4, deg <= 90)
# phi4_90t_2 <- filter(phi4_90t, angle == 'phi2_phi1')
# hist4_90t <- ggplot(data = phi4_90t, aes(deg)) + 
#   geom_histogram(binwidth = 2) +
#   theme_minimal()
# print(hist4_90t)

# Simple mean, median for each plant - take median across all obs remaining for first 4 leaves
phi4_90t_s <- summarise(phi4_90t, 
                        med = median(deg), 
                        am = mean(deg), n_obs = n(),  
                        cmed3 = median(deg[angle %in% c('phi1_phi0', 'phi2_phi1', 'phi3_phi2')], na.rm = TRUE), 
                        cam2 = mean(deg[angle %in% c('phi1_phi0', 'phi2_phi1')], na.rm = TRUE), 
                        cam3 = mean(deg[angle %in% c('phi1_phi0', 'phi2_phi1', 'phi3_phi2')], na.rm = TRUE), 
                        med1 = median(deg[angle == 'phi1_phi0'], na.rm = TRUE),
                        med2 = median(deg[angle == 'phi2_phi1'], na.rm = TRUE),
                        med3 = median(deg[angle == 'phi3_phi2'], na.rm = TRUE), 
                        med4 = median(deg[angle == 'phi4_phi3'], na.rm = TRUE),
                        am1 = mean(deg[angle == 'phi1_phi0'], na.rm = TRUE), 
                        am2 = mean(deg[angle == 'phi2_phi1'], na.rm = TRUE),
                        am3 = mean(deg[angle == 'phi3_phi2'], na.rm = TRUE), 
                        am4 = mean(deg[angle == 'phi4_phi3'], na.rm = TRUE),
                        #cmed3_pm = median(med1, med2, med3, na.rm = TRUE), 
                        #cmed4_pm = median(med1, med2, med3, med4, na.rm = TRUE),
                        #cam2_pm = mean(med1, med2, na.rm = TRUE), 
                        #cam3_pm = mean(med1, med2, med3, na.rm = TRUE), 
                        #cam4_pm = mean(med1, med2, med3, med4, na.rm = TRUE),
                        .by = plant_num:PI_num)
phi4_90t_s$plant_num[314] <- 49

phi4_90t_s <- mutate(rowwise(phi4_90t_s), 
                     cmed3_pmed = median(med1, med2, med3, na.rm = TRUE),
                     cmed4_pmed = median(med1, med2, med3, med4, na.rm = TRUE),
                     cmed2_pam = median(am1, am2, na.rm = TRUE),
                     cmed3_pam = median(am1, am2, am3, na.rm = TRUE),
                     cmed4_pam = median(am1, am2, am3, am4, na.rm = TRUE),
                     cam2_pmed = mean(c(med1, med2), na.rm =TRUE),
                     cam3_pmed = mean(c(med1, med2, med3), na.rm = TRUE), 
                     cam4_pmed = mean(c(med1, med2, med3, med4), na.rm = TRUE),
                     cam2_pam = mean(c(am1, am2), na.rm = TRUE),
                     cam3_pam = mean(c(am1, am2, am3), na.rm = TRUE),
                     cam4_pam = mean(c(am1, am2, am3, am4), na.rm = TRUE))


# hist4_90t_s_med <- ggplot(phi4_90t_s, aes(med)) + 
#   geom_histogram(binwidth = 2) +
#   theme_minimal() 
# print(hist4_90t_s_med)
# 
# hist4_90t_s_am <- ggplot(phi4_90t_s, aes(am)) +
#   geom_histogram(binwidth = 2) +
#   theme_minimal()
# print(hist4_90t_s_am)
# 
# hist4_90t_s_cmed3 <- ggplot(phi4_90t_s, aes(cmed3)) +
#   geom_histogram(binwidth = 2) +
#   theme_minimal()
# print(hist4_90t_s_cmed3)
# 
# hist4_90t_s_cam2 <- ggplot(phi4_90t_s, aes(cam2)) +
#   geom_histogram(binwidth = 2) +
#   theme_minimal()
# print(hist4_90t_s_cam2)
# 
# hist4_90t_s_cam3 <- ggplot(phi4_90t_s, aes(cam3)) + 
#   geom_histogram(binwidth = 2) +
#   theme_minimal()
# print(hist4_90t_s_cam3)
# 
# hist4_90t_s_med1 <- ggplot(phi4_90t_s, aes(med1)) + 
#   geom_histogram(binwidth = 2) +
#   theme_minimal()
# print(hist4_90t_s_med1)
# 
# hist4_90t_s_med2 <- ggplot(phi4_90t_s, aes(med2)) + 
#   geom_histogram(binwidth = 2) +
#   theme_minimal()
# print(hist4_90t_s_med2)
# 
# hist4_90t_s_med3 <- ggplot(phi4_90t_s, aes(med3)) + 
#   geom_histogram(binwidth = 2) +
#   theme_minimal()
# print(hist4_90t_s_med3)
# 
# hist4_90t_s_med4 <- ggplot(phi4_90t_s, aes(med4)) + 
#   geom_histogram(binwidth = 2) +
#   theme_minimal()
# print(hist4_90t_s_med4)
# 
# 
# hist4_90t_s_am1 <- ggplot(phi4_90t_s, aes(am1)) + 
#   geom_histogram(binwidth = 2) +
#   theme_minimal()
# print(hist4_90t_s_am1)
# 
# hist4_90t_s_am2 <- ggplot(phi4_90t_s, aes(am2)) + 
#   geom_histogram(binwidth = 2) +
#   theme_minimal()
# print(hist4_90t_s_am2)
# 
# hist4_90t_s_am3 <- ggplot(phi4_90t_s, aes(am3)) + 
#   geom_histogram(binwidth = 2) +
#   theme_minimal()
# print(hist4_90t_s_am3)
# 
# hist4_90t_s_am4 <- ggplot(phi4_90t_s, aes(am4)) + 
#   geom_histogram(binwidth = 2) +
#   theme_minimal()
# print(hist4_90t_s_am4)

h_med <- calcHeritability2(phi4_90t_s, med)
h_am <- calcHeritability2(phi4_90t_s, am)
h_cmed3 <- calcHeritability2(phi4_90t_s, cmed3)
h_cam2 <- calcHeritability2(phi4_90t_s, cam2)
h_cam3 <- calcHeritability2(phi4_90t_s, cam3)
h_med1 <- calcHeritability2(phi4_90t_s, med1)
h_med2 <- calcHeritability2(phi4_90t_s, med2)
h_med3 <- calcHeritability2(phi4_90t_s, med3)
h_med4 <- calcHeritability2(phi4_90t_s, med4)
h_am1 <- calcHeritability2(phi4_90t_s, am1)
h_am2 <- calcHeritability2(phi4_90t_s, am2)
h_am3 <- calcHeritability2(phi4_90t_s, am3)
h_am4 <- calcHeritability2(phi4_90t_s, am4)
h_cmed3_pmed <- calcHeritability2(phi4_90t_s, cmed3_pmed)
h_cmed4_pmed <- calcHeritability2(phi4_90t_s, cmed4_pmed)
h_cmed2_pam <- calcHeritability2(phi4_90t_s, cmed2_pam)
h_cmed3_pam <- calcHeritability2(phi4_90t_s, cmed3_pam)
h_cmed4_pam <- calcHeritability2(phi4_90t_s, cmed4_pam)
h_cam2_pmed <- calcHeritability2(phi4_90t_s, cam2_pmed)
h_cam3_pmed <- calcHeritability2(phi4_90t_s, cam3_pmed)
h_cam4_pmed <- calcHeritability2(phi4_90t_s, cam4_pmed)
h_cam2_pam <- calcHeritability2(phi4_90t_s, cam2_pam)
h_cam3_pam <- calcHeritability2(phi4_90t_s, cam3_pam)
h_cam4_pam <- calcHeritability2(phi4_90t_s, cam4_pam)

geno_reps <- phi4_90t_s %>%
  ungroup() %>% 
  group_by(PI_num) %>%
  summarise(n = n())

## BLUE GLM
glm_fit <- glm(deg ~ factor(plant_num), family = Gamma(link = 'identity'), data = phi4_90t)
phi4_90t$blue_glm <- glm_fit$fitted.values
phi4_blue <- summarise(phi4_90t, blue_glm = mean(blue_glm), .by = c(plant_num, PI_num))

phi4_90t_s <- full_join(phi4_90t_s, phi4_blue, join_by(plant_num, PI_num), keep = FALSE)
h_blueGLM <- calcHeritability2(phi4_blue, blue_glm)

hist_a <- ggplot() + 
  geom_histogram(aes(med), data = phe_orig, binwidth = 2, fill = colors2[1]) +
  geom_histogram(aes(med), data = phi4_90t_s, binwidth = 2, fill = colors2[2]) + 
  scale_x_continuous(name = expression(Median~of~Phi[1-4]),
                     breaks = c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180), 
                     labels = c(expression('0'*degree), expression('15'*degree), expression('30'*degree),
                                expression('45'*degree), expression('60'*degree), expression('75'*degree),
                                expression('90'*degree), expression('105'*degree), expression('120'*degree),
                                expression('135'*degree), expression('150'*degree), expression('165'*degree),
                                expression('180'*degree))) +
  labs(x = expression(Median~of~Phi[1-4]), y = 'Frequency') +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0), 
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 11, color = 'black', vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11, color = 'black'),
        text = element_text(size = 11, color = 'black'),
        legend.position = 'bottom',
        line = element_line(color = 'black', linewidth = 1),
        panel.grid = element_blank()) 
hist_a
ggsave('../Images/hist_a.svg', plot = hist_a, width = 6.5, height = 2.72, units = 'in', bg = 'white', dpi = 1000)

# Export phenotype data
keep_cols <- c('PI_num', 'plant_num', 'am', 'am2', 'am3', 'blue_glm', 'cam3', 'cam3_pmed', 'cam4_pam', 'cam4_pmed', 'cmed3', 'med', 'med2', 'med3')
pheno_gwas2 <- phi4_90t_s[, which((names(phi4_90t_s) %in% keep_cols)==TRUE)]
pheno_gwas2$PI_num <- sub('PI_', 'PI', pheno_gwas2$PI_num)
write.csv(pheno_gwas2, file = '../Data/phenotypeData_490t.csv', row.names = FALSE) #column names manually changed in exported phenotype data to end with underscores

phi4_90t_s %>%
  select(PI_num, cam3_pam) %>%
  dplyr::rename(cam3pam_ = cam3_pam) %>%
  rowwise() %>%
  mutate(PI_num = str_remove(PI_num, '_')) %>%
  write.csv('../Data/phenotypeData_cam3pam.csv', row.names = FALSE, quote = FALSE)
# Get list of unique PI_nums appearing in phenotype data
genotypes <- pheno_gwas2$PI_num
genotypes <- unique(genotypes)
write.csv(genotypes, file = '../Data/genotypes.csv', quote = FALSE, row.names = FALSE)

# Check reliability vs. manual method
manualDataNew <- read.csv('../Data/leafLevelTraits.csv') %>%
  dplyr::rename(plantNum = plantID)

reconstructionDataNew <- read_excel('../Data/angles.xlsx', skip = 1)[, c(1, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35)]
colnames(reconstructionDataNew) <- c('plantID', 'phi1', 'phi2', 'phi3', 'phi4', 'phi5', 'phi6', 'phi7', 'phi8', 'phi9', 'phi10', 'phi11')

reconstructionDataNew <- reconstructionDataNew %>%
  rowwise() %>%
  mutate(imageRun = str_split_i(plantID, '-', 9) %>%
           str_split_i('_', 1),
         plantNum = str_split_i(plantID, '-', 8))

secondRunReconstructions <- filter(reconstructionDataNew, imageRun=='2')
firstRunReconstructions <- filter(reconstructionDataNew, imageRun=='03')

reconstructionDataNew <- firstRunReconstructions %>%
  filter(plantNum=='189_2024') %>%
  bind_rows(secondRunReconstructions) %>%
  rowwise() %>%
  mutate(plantNum = str_remove(plantNum, '_2024') %>%
           as.numeric()) %>%
  pivot_longer(starts_with('phi'), names_to = 'leafNumber', values_to = 'phi', names_prefix = 'phi') %>%
  mutate(leafNumber = as.numeric(leafNumber))
# filter(imageRun=='03') %>%
# select(c(plantID, starts_with('phi'))) %>%
# mutate(plantNum = str_split_i(plantID, '_', 2) %>%
#          str_split_i('-', 3) %>%
#          as.numeric()) %>%
# ungroup() %>%
# pivot_longer(starts_with('phi'), names_to = 'leafNumber', values_to = 'phi', names_prefix = 'phi') %>%
# mutate(leafNumber = as.numeric(leafNumber))

val_raw <- full_join(manualDataNew, reconstructionDataNew, join_by(plantNum==plantNum, leafNumber==leafNumber)) %>%
  filter(!is.na(compassDirectionNikee)|!is.na(compassDirectionJensina)) %>%
  rowwise() %>%
  mutate(phiNorm = normalize360(phi)) %>%
  select(plantNum, leafNumber, phiNorm, compassDirectionNikee, compassDirectionJensina) %>%
  filter(leafNumber <= 6) %>%
  pivot_wider(id_cols = plantNum, names_from = leafNumber, values_from = c(contains('compassDirection'), phiNorm)) %>%
  rowwise() %>%
  mutate(varphi_1_n = normalize360(compassDirectionNikee_2 - compassDirectionNikee_1),
         varphi_2_n = normalize360(compassDirectionNikee_3 - compassDirectionNikee_2),
         varphi_3_n = normalize360(compassDirectionNikee_4 - compassDirectionNikee_3),
         varphi_4_n = normalize360(compassDirectionNikee_5 - compassDirectionNikee_4),
         varphi_5_n = normalize360(compassDirectionNikee_6 - compassDirectionNikee_5),
         varphi_1_j = normalize360(compassDirectionJensina_2 - compassDirectionJensina_1),
         varphi_2_j = normalize360(compassDirectionJensina_3 - compassDirectionJensina_2),
         varphi_3_j = normalize360(compassDirectionJensina_4 - compassDirectionJensina_3),
         varphi_4_j = normalize360(compassDirectionJensina_5 - compassDirectionJensina_4),
         varphi_5_j = normalize360(compassDirectionJensina_6 - compassDirectionJensina_5),
         varphi_1_r = normalize360(phiNorm_2 - phiNorm_1),
         varphi_2_r = normalize360(phiNorm_3 - phiNorm_2),
         varphi_3_r = normalize360(phiNorm_4 - phiNorm_3),
         varphi_4_r = normalize360(phiNorm_5 - phiNorm_4),
         varphi_5_r = normalize360(phiNorm_6 - phiNorm_5)) %>%
  select(c(plantNum, contains('varphi'))) %>%
  pivot_longer(contains('varphi'), names_to = 'col', values_to = 'varphi', names_prefix = 'varphi_') %>%
  rowwise() %>%
  mutate(method = str_split_i(col, '_', 2),
         leafNumber = str_split_i(col, '_', 1)) %>%
  select(plantNum, method, leafNumber, varphi) %>%
  pivot_wider(id_cols = c(plantNum, leafNumber),
              values_from = varphi,
              names_from = method) %>%
  group_by(plantNum) %>%
  mutate(plantCorr.RJ = cor(j, r, use = 'complete.obs'),
         plantCorr.RN = cor(n, r, use = 'complete.obs'),
         plantCorr.JN = cor(j, n, use = 'complete.obs')) %>%
  rowwise() %>%
  mutate(rConjugate.J = case_when(plantCorr.RJ < 0 ~ 360 - r, .default = r),
         rConjugate.N = case_when(plantCorr.RN < 0 ~ 360 - r, .default = r),
         nConjugate.J = case_when(plantCorr.JN < 0 ~ 360 - n, .default = n),
         PHI_J = normalize180(j),
         PHI_N = normalize180(n),
         PHI_R = normalize180(r))

r2J <- round(cor(val_raw$j, val_raw$r, use = 'complete.obs')^2, digits = 4)
r2N <- round(cor(val_raw$n, val_raw$r, use = 'complete.obs')^2, digits = 4)

r2ConjugateJ <- round(cor(val_raw$j, val_raw$rConjugate.J, use = 'complete.obs')^2, digits = 4)
r2ConjugateN <- round(cor(val_raw$j, val_raw$rConjugate.N, use = 'complete.obs')^2, digits = 4)

PHI_r2J <- round(cor(val_raw$PHI_J, val_raw$PHI_R, use = 'complete.obs')^2, digits = 4)
PHI_r2N <- round(cor(val_raw$PHI_N, val_raw$PHI_R, use = 'complete.obs')^2, digits = 4)

val_filt <- filter(val_raw, PHI_R <= 90)
val_filtJ <- filter(val_filt, PHI_J <= 90)
val_filtN <- filter(val_filt, PHI_N <= 90)

r2Jfilt <- round(cor(val_filtJ$j, val_filtJ$r, use = 'complete.obs')^2, digits = 4)
r2Nfilt <- round(cor(val_filtN$n, val_filtN$r, use = 'complete.obs')^2, digits = 4)

r2ConjugateJfilt <- round(cor(val_filtJ$j, val_filtJ$rConjugate.J, use = 'complete.obs')^2, digits = 4)
r2ConjugateNfilt <- round(cor(val_filtN$n, val_filtN$rConjugate.N, use = 'complete.obs')^2, digits = 4)

conjugateNfiltPlot <- ggplot(val_filtN, aes(n, rConjugate.N, use = 'complete.obs')) +
  geom_point() +
  scale_x_continuous(limits = c(90, 270)) +
  scale_y_continuous(limits = c(90, 270)) +
  labs(x = 'Nikee',
       y = 'Reconstruction',
       subtitle = paste0('R-squared = ', r2ConjugateNfilt))

combinedValidation <- val_filtJ %>%
  dplyr::rename(manual = j,
         reconstruction = rConjugate.J) %>%
  select(plantNum, leafNumber, manual, reconstruction) %>%
  mutate(person = 'Individual 1')
combinedValidation2 <- val_filtN %>%
  dplyr::rename(manual = n,
         reconstruction = rConjugate.N) %>%
  select(plantNum, leafNumber, manual, reconstruction) %>%
  mutate(person = 'Individual 2')
combinedValidation <- bind_rows(combinedValidation2, combinedValidation)
combinedR2 <- round(cor(combinedValidation$manual, combinedValidation$reconstruction, use = 'complete.obs')^2, digits = 2)

combinedPlot <- ggplot(combinedValidation, aes(manual, reconstruction, color = person)) +
  geom_point() +
  scale_x_continuous(breaks = c(90, 135, 180, 225, 270),
                     labels = c(expression('90'*degree), expression('135'*degree),
                                expression('180'*degree), expression('225'*degree),
                                expression('270'*degree)),
                     limits = c(90, 270)) +
  scale_y_continuous(breaks = c(90, 135, 180, 225, 270),
                     labels = c(expression('90'*degree), expression('135'*degree),
                                expression('180'*degree), expression('225'*degree),
                                expression('270'*degree)),
                     limits = c(90, 270)) +
  scale_color_manual(values = colors2) +
  labs(x = 'Manual', y = 'Reconstruction', color = NULL,
       title = bquote('R'^2 == .(combinedR2))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0),
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 11, color = 'black', vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11, color = 'black'),
        plot.title = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0),
                                  vjust = 0.5, hjust = 0.5),
        text = element_text(size = 11, color = 'black'),
        legend.position = 'none',
        line = element_line(color = 'black', linewidth = 1),
        panel.grid = element_blank())
combinedPlot

manualReliability <- filter(val_raw, PHI_J <=90 & PHI_N <= 90)
r2jnfilt <- round(cor(manualReliability$j, manualReliability$nConjugate.J, use = 'complete.obs')^2, digits = 2)
manualReliabilityPlot <- ggplot(manualReliability, aes(j, nConjugate.J)) +
  geom_point(color = colors2[2]) +
  scale_x_continuous(breaks = c(90, 135, 180, 225, 270),
                     labels = c(expression('90'*degree), expression('135'*degree),
                                expression('180'*degree), expression('225'*degree),
                                expression('270'*degree)),
                     limits = c(90, 270)) +
  scale_y_continuous(breaks = c(90, 135, 180, 225, 270),
                     labels = c(expression('90'*degree), expression('135'*degree),
                                expression('180'*degree), expression('225'*degree),
                                expression('270'*degree)),
                     limits = c(90, 270)) +
  labs(x = 'Individual 1 (Manual)', y = 'Individual 2 (Manual)',
       title = bquote('R'^2 == .(r2jnfilt))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0),
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 11, color = 'black', vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11, color = 'black'),
        plot.title = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0),
                                  vjust = 0.5, hjust = 0.5),
        text = element_text(size = 11, color = 'black'),
        legend.position = 'none',
        line = element_line(color = 'black', linewidth = 1),
        panel.grid = element_blank())
reliabilityFig <- plot_grid(manualReliabilityPlot, reconstructionReliabilityPlot, combinedPlot,
                            nrow = 1, labels = 'AUTO')
reliabilityFig
# ggsave('../Images/reliabilityPlots.svg', width = 6.5, height = 2.5, units = 'in', dpi = 1000, bg = 'white')
 
#Summarize GWAS results
path <- '../Data/signals/'
phenotypeData <- read.csv(paste0(path, 'phenotypes.csv'))
traits <- colnames(phenotypeData)[-1]
df_rmip <- getSupport(path, traits[1])
for (j in 2:length(traits))
{
  df_rmip <- bind_rows(df_rmip, getSupport(path, traits[j]))
}

cam3pamRMIP <- getSupport(paste0(path, '/cam3pam/'), 'cam3pam_')
df_rmip <- bind_rows(df_rmip, cam3pamRMIP)
traits <- c(traits, 'cam3pam_')

rmip_summary <- df_rmip %>%
  pivot_longer(any_of(traits), names_to = 'trait', 
               values_to = 'RMIP') %>%
  filter(!is.na(RMIP)) %>%
  group_by(trait) %>%
  summarise(maxRMIP = max(RMIP, na.rm = TRUE), 
            numHits = n())
rmip0.02 <- df_rmip %>%
  pivot_longer(any_of(traits), names_to = 'trait', 
               values_to = 'RMIP') %>%
  filter(!is.na(RMIP)) %>%
  filter(RMIP > 0.02) %>% 
  arrange(CHROM, POS)
rmip0.02SNPsummary <- rmip0.02 %>%
  group_by(SNP) %>%
  summarise(nTraits = n(), 
            maxRMIP = max(RMIP, na.rm = TRUE))

manhattan_med <- plotRMIPManhattan(df_rmip, med_, colorVec = colors, title = '')
# ggsave('../Images/manhattan_workflow.png', plot = manhattan_med, width = 2.3, height = 1.82, units = 'in', dpi=1000)

manhattan_am <- plotRMIPManhattan(df_rmip, am_, colorVec = colors, title = expression(Mean~of~Phi[1-4]))

manhattan_cmed3 <- plotRMIPManhattan(df_rmip, cmed3_, colorVec = colors, title = expression(Median~of~Phi[1-3]))

manhattan_cam3 <- plotRMIPManhattan(df_rmip, cam3_, colorVec = colors, title = expression(Mean~of~Phi[1-3]))

manhattan_med2 <- plotRMIPManhattan(df_rmip, med2_, colorVec = colors, title = expression(Median~of~Phi[2]))

manhattan_med3 <- plotRMIPManhattan(df_rmip, med3_, colorVec = colors, title = expression(Median~of~Phi[3]))

manhattan_am2 <- plotRMIPManhattan(df_rmip, am2_, colorVec = colors, title = expression(Mean~of~Phi[2]))

manhattan_am3 <- plotRMIPManhattan(df_rmip, am3_, colorVec = colors, title = expression(Mean~of~Phi[3]))

manhattan_cam3pmed <- plotRMIPManhattan(df_rmip, cam3pmed_, colorVec = colors, title = expression(Mean~of~Median~Phi[1-3]))

manhattan_cam4pmed <- plotRMIPManhattan(df_rmip, cam4pmed_, colorVec = colors, title = expression(Mean~of~Median~Phi[1-4]))

manhattan_cam4pam <- plotRMIPManhattan(df_rmip, cam4pam_, colorVec = colors, title = expression(Mean~of~Mean~Phi[1-4]))

manhattan_blueglm <- plotRMIPManhattan(df_rmip, blueglm, colorVec = colors, title = expression(atop(Mean~GLM~BLUE, Fitted~Phi[1-4])))

manhattan_cam3pam <- plotRMIPManhattan(df_rmip, cam3pam_, colorVec = colors, title = expression(Mean~of~Mean~Phi[1-3]))

sManhattan <- plot_grid(manhattan_am, manhattan_cam4pmed, manhattan_cmed3, 
                        manhattan_med2, manhattan_am2, manhattan_med3, 
                        manhattan_cam3, manhattan_am3, manhattan_cam3pmed,
                        manhattan_cam3pam, manhattan_cam4pam, manhattan_blueglm, 
                        nrow = 4, ncol = 3, labels = 'AUTO')
sManhattan
# ggsave('../Images/supplementalManhattanPlots.svg', plot = sManhattan, width = 6.5, height = 9, units = 'in', dpi = 1000, 
       # bg = 'white')

hist_med <- ggplot(phenotypeData, aes(med_)) +
  geom_histogram(binwidth = 2, fill = colors2[2]) + 
  scale_x_continuous(name = expression(Median~of~Phi[1-4]),
                     breaks = c(0, 15, 30, 45, 60, 75, 90), 
                     labels = c(expression('0'*degree), expression('15'*degree), expression('30'*degree),
                                expression('45'*degree), expression('60'*degree), 
                                expression('75'*degree), expression('90'*degree))) +
  labs(x = expression(Median~of~Phi[1-4]), y = 'Frequency') +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0), 
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 11, color = 'black', vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11, color = 'black'),
        text = element_text(size = 11, color = 'black'),
        legend.position = 'bottom',
        line = element_line(color = 'black', linewidth = 1),
        panel.grid = element_blank()) 
hist_med
ggsave('../Images/histogram_workflow.png', plot = hist_med, width = 2.3, height = 1.82, units = 'in', dpi = 1000)

fig6top <- plot_grid(hist_med, manhattan_med, nrow = 1, labels = 'AUTO')
fig6top

PATH_SIGSNPSVCF <- '../Data/sigSNPSMed.tsv'
sigSNPsVCF <- read.table(PATH_SIGSNPSVCF, sep = '\t', comment.char = ' ')
colnames(sigSNPsVCF) <- sigSNPsVCF[1, ]

allele_state <- sigSNPsVCF[2:4,] %>%
  as_tibble() %>%
  pivot_longer(starts_with('PI'), values_to = 'state', names_to = 'genotype')

med <- phenotypeData %>%
  select(PI_num, med_)

allele_state_med <- full_join(allele_state, med, join_by(genotype == PI_num), keep = FALSE, 
                              relationship = 'many-to-many') %>%
  rowwise() %>%
  mutate(alleleState1 = str_sub(state, 1, 1) %>%
           as.numeric(),
         alleleState2 = str_sub(state, 3, 3) %>%
           as.numeric()) %>%
  filter(alleleState1==alleleState2) %>%
  mutate(alleleStateFactor = case_when(alleleState1 == 1 ~ 'Minor',
                                       .default = 'Major'))

allele_state_summary <- allele_state_med %>%
  group_by(ID, genotype, alleleStateFactor) %>%
  summarise(med = mean(med_, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(ID, alleleStateFactor) %>%
  summarise(genotypes = n())

allele_state5.1 <- allele_state_med %>%
  filter(str_detect(ID, 'Chr05_12109370'))
fit5.1 <- lm(med_ ~ alleleStateFactor, data = allele_state5.1)
anova(fit5.1)

box5.1 <- allele_state5.1 %>%
  ggplot(aes(alleleStateFactor, med_, fill = alleleStateFactor)) +
  geom_boxplot(color = 'black') +
  scale_fill_manual(values = colors2) +
  scale_y_continuous(name = expression(Median~of~Phi[1-4]),
                     breaks = c(0, 15, 30, 45, 60, 75, 90), 
                     labels = c(expression('0'*degree), expression('15'*degree), expression('30'*degree),
                                expression('45'*degree), expression('60'*degree), 
                                expression('75'*degree), expression('90'*degree))) +
  scale_x_discrete(breaks = waiver(),
                   labels = c(paste0('Major', '\n', 'n=201'), 
                              paste0('Minor', '\n', 'n=8'))) +
  labs(x = NULL, y = expression(Median~of~Phi[1-4]), 
       fill = '') +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0), 
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 11, color = 'black', vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11, color = 'black'),
        text = element_text(size = 11, color = 'black'),
        legend.position = 'none',
        line = element_line(color = 'black', linewidth = 1),
        panel.grid = element_blank())
box5.1

allele_state5.2 <- allele_state_med %>%
  filter(str_detect(ID, 'Chr05_65733791'))
fit5.2 <- lm(med_ ~ alleleStateFactor, data = allele_state5.2)
anova(fit5.2)

box5.2 <- allele_state5.2 %>%
  ggplot(aes(alleleStateFactor, med_, fill = alleleStateFactor)) +
  geom_boxplot(color = 'black') +
  scale_fill_manual(values = colors2) +
  scale_y_continuous(name = expression(Median~of~Phi[1-4]),
                     breaks = c(0, 15, 30, 45, 60, 75, 90), 
                     labels = c(expression('0'*degree), expression('15'*degree), expression('30'*degree),
                                expression('45'*degree), expression('60'*degree), 
                                expression('75'*degree), expression('90'*degree))) +
  scale_x_discrete(breaks = waiver(),
                   labels = c(paste0('Major', '\n', 'n=193'), 
                              paste0('Minor', '\n', 'n=13'))) +
  labs(x = NULL, y = expression(Median~of~Phi[1-4]), 
       fill = '') +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0), 
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 11, color = 'black', vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11, color = 'black'),
        text = element_text(size = 11, color = 'black'),
        legend.position = 'none',
        line = element_line(color = 'black', linewidth = 1),
        panel.grid = element_blank())
box5.2

allele_state6 <- allele_state_med %>%
  filter(str_detect(ID, 'Chr06_41390777'))
fit6 <- lm(med_ ~ alleleStateFactor, data = allele_state6)
anova(fit6)

box6 <- allele_state6 %>%
  ggplot(aes(alleleStateFactor, med_, fill = alleleStateFactor)) +
  geom_boxplot(color = 'black') +
  scale_fill_manual(values = colors2) +
  scale_y_continuous(name = expression(Median~of~Phi[1-4]),
                     breaks = c(0, 15, 30, 45, 60, 75, 90), 
                     labels = c(expression('0'*degree), expression('15'*degree), expression('30'*degree),
                                expression('45'*degree), expression('60'*degree), 
                                expression('75'*degree), expression('90'*degree))) +
  scale_x_discrete(breaks = waiver(),
                   labels = c(paste0('Major', '\n', 'n=188'), 
                              paste0('Minor', '\n', 'n=14'))) +
  labs(x = NULL, y = expression(Median~of~Phi[1-4]), 
       fill = '') +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0), 
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 11, color = 'black', vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11, color = 'black'),
        text = element_text(size = 11, color = 'black'),
        legend.position = 'none',
        line = element_line(color = 'black', linewidth = 1),
        panel.grid = element_blank())
box6

fig6boxes <- plot_grid(box5.1, box5.2, box6, nrow = 1, labels = c('C', 'D', 'E'))
fig6boxes

gff <- '/Users/jensinadavis/davisjensina@gmail.com - Google Drive/My Drive/Schnable-Lab/Sbicolor_454_v3.1.1.gene_exons.gff3'
gene_models <- import(gff, format = 'gff3')
seq_names <- unique(gene_models@seqinfo@seqnames)

WINDOW_SIZE <- 100000

chrom5.1site <- 12109370
chrom5.1start <- chrom5.1site - WINDOW_SIZE
chrom5.1end <- chrom5.1site + WINDOW_SIZE
chrom5.1_reg <- GRanges(seqnames = 'Chr05', 
                        ranges = IRanges(start = (chrom5.1start), end = chrom5.1end))
chrom5.1_reg_genes <- subsetByOverlaps(gene_models, chrom5.1_reg)
genes5.1_df <- as.data.frame(chrom5.1_reg_genes)
genes5.1_df <- genes5.1_df %>%
  as_tibble() %>% 
  filter(str_detect(type, 'gene')) %>%
  rowwise() %>%
  mutate(midpt = (start + end)/2)

chrom5.2site <- 65733791
chrom5.2start <- chrom5.2site - WINDOW_SIZE
chrom5.2end <- chrom5.2site + WINDOW_SIZE
chrom5.2_reg <- GRanges(seqnames = 'Chr05', 
                        ranges = IRanges(start = (chrom5.2start), end = chrom5.2end))
chrom5.2_reg_genes <- subsetByOverlaps(gene_models, chrom5.2_reg)
genes5.2_df <- as.data.frame(chrom5.2_reg_genes)
genes5.2_df <- genes5.2_df %>%
  as_tibble() %>% 
  filter(str_detect(type, 'gene')) %>%
  rowwise() %>%
  mutate(midpt = (start + end)/2)


chrom6site <- 41390777
chrom6start <- chrom6site - WINDOW_SIZE
chrom6end <- chrom6site + WINDOW_SIZE
chrom6_reg <- GRanges(seqnames = 'Chr06', 
                      ranges = IRanges(start = (chrom6start), end = chrom6end))
chrom6_reg_genes <- subsetByOverlaps(gene_models, chrom6_reg)
genes6_df <- as.data.frame(chrom6_reg_genes)
genes6_df <- genes6_df %>%
  as_tibble() %>% 
  filter(str_detect(type, 'gene')) %>%
  rowwise() %>%
  mutate(midpt = (start + end)/2)

ld_colnames <- c('CHROM', 'POS', 'SNP', 'R2')
ld5.1 <- read.table('../Data/LD_Chr05_12109370.ld', header = TRUE)[, 4:7]
colnames(ld5.1) <- ld_colnames
ld5.1 <- ld5.1 %>%
  arrange(POS) %>%
  filter(POS %in% chrom5.1start:chrom5.1end)
for(i in 1:length(ld5.1$SNP))
{
  if(i!=1)
  {
    ld5.1$xmin[i] <- (ld5.1$POS[i - 1] + ld5.1$POS[i])/2
  }
  
  if(i!=length(ld5.1$SNP))
  {
    ld5.1$xmax[i] <- (ld5.1$POS[i + 1] + ld5.1$POS[i])/2
  }
  
  if(i==1)
  {
    ld5.1$xmin[i] <- ld5.1$POS[i]
  }
  if(i==length(ld5.1$SNP))
  {
    ld5.1$xmax[i] <- ld5.1$POS[i]
  }
}

ld5.2 <- read.table('../Data/LD_Chr05_65733791.ld', header = TRUE)[, 4:7]
colnames(ld5.2) <- ld_colnames
ld5.2 <- ld5.2 %>%
  arrange(POS) %>%
  filter(POS %in% chrom5.2start:chrom5.2end)
for(i in 1:length(ld5.2$SNP))
{
  if(i!=1)
  {
    ld5.2$xmin[i] <- (ld5.2$POS[i - 1] + ld5.2$POS[i])/2
  }
  
  if(i!=length(ld5.2$SNP))
  {
    ld5.2$xmax[i] <- (ld5.2$POS[i + 1] + ld5.2$POS[i])/2
  }
  
  if(i==1)
  {
    ld5.2$xmin[i] <- ld5.2$POS[i]
  }
  if(i==length(ld5.2$SNP))
  {
    ld5.2$xmax[i] <- ld5.2$POS[i]
  }
}

ld6 <- read.table('../Data/LD_Chr06_41390777.ld', header = TRUE)[, 4:7]
colnames(ld6) <- ld_colnames
ld6 <- ld6 %>%
  arrange(POS) %>%
  filter(POS %in% chrom6start:chrom6end)
for(i in 1:length(ld6$SNP))
{
  if(i!=1)
  {
    ld6$xmin[i] <- (ld6$POS[i - 1] + ld6$POS[i])/2
  }
  
  if(i!=length(ld6$SNP))
  {
    ld6$xmax[i] <- (ld6$POS[i + 1] + ld6$POS[i])/2
  }
  
  if(i==1)
  {
    ld6$xmin[i] <- ld6$POS[i]
  }
  if(i==length(ld6$SNP))
  {
    ld6$xmax[i] <- ld6$POS[i]
  }
}

p5.1reg <- ggplot() + 
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 0.1), data = genes5.1_df, fill = colors2[1]) +
  geom_point(aes(chrom5.1site, 0)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -0.1, ymax = 0, fill = R2), data = ld5.1) +
  geom_hline(yintercept = 0) +
  geom_text_repel(aes(x = midpt, y = 0, label = Name), data = genes5.1_df, direction = 'both', 
                  box.padding = 0.6, min.segment.length = 0, size = 3.881, max.overlaps = 100) +
  scale_x_continuous(limits = c(chrom5.1start, chrom5.1end)) +
  scale_y_continuous(limits = c(-0.35, 0.35)) +
  scale_fill_viridis(direction = -1, option = 'mako') +
  guides(fill = guide_colourbar(barwidth = 12,
                                barheight = 1)) +
  labs(x = NULL, y = NULL, fill = expression(LD~R^2)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        text = element_text(size = 11, color = 'black'),
        legend.text = element_text(size = 11, color = 'black', hjust = 0.5),
        line = element_line(color = 'black', linewidth = 1),
        panel.grid = element_blank(), 
        legend.position = 'bottom')
p5.1reg

ld_legend <- get_plot_component(p5.1reg, 'guide-box-bottom', return_all = TRUE)

p5.1reg <- ggplot() + 
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 0.1), data = genes5.1_df, fill = colors2[1]) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -0.1, ymax = 0, fill = R2), data = ld5.1) +
  geom_hline(yintercept = 0) +
  geom_point(aes(chrom5.1site, 0)) +
  geom_text_repel(aes(x = midpt, y = 0, label = Name), data = genes5.1_df, direction = 'both', 
                  box.padding = 0.6, min.segment.length = 0, size = 3.881, max.overlaps = 100) +
  scale_x_continuous(limits = c(chrom5.1start, chrom5.1end)) +
  scale_y_continuous(limits = c(-0.35, 0.35)) +
  scale_fill_viridis(direction = -1, option = 'mako') +
  guides(fill = guide_colourbar(barwidth = 12,
                                barheight = 1)) +
  labs(x = NULL, y = NULL, fill = expression(LD~R^2)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        text = element_text(size = 11, color = 'black'),
        legend.text = element_text(size = 11, color = 'black', hjust = 0.5),
        line = element_line(color = 'black', linewidth = 1),
        panel.grid = element_blank(), 
        legend.position = 'none')
p5.1reg

p5.1reg <- ggplot() + 
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 0.1), data = genes5.1_df, fill = colors2[1]) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -0.1, ymax = 0, fill = R2), data = ld5.1) +
  geom_hline(yintercept = 0) +
  geom_point(aes(chrom5.1site, 0), size = 3) +
  geom_text_repel(aes(x = midpt, y = 0, label = Name), data = genes5.1_df, direction = 'both', 
                  box.padding = 0.6, min.segment.length = 0, size = 3.881, max.overlaps = 100) +
  scale_x_continuous(limits = c(chrom5.1start, chrom5.1end)) +
  scale_y_continuous(limits = c(-0.35, 0.35)) +
  scale_fill_viridis(direction = -1, option = 'mako') +
  guides(fill = guide_colourbar(barwidth = 12,
                                barheight = 1)) +
  labs(x = NULL, y = NULL, fill = expression(LD~R^2)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        text = element_text(size = 11, color = 'black'),
        legend.text = element_text(size = 11, color = 'black', hjust = 0.5),
        line = element_line(color = 'black', linewidth = 1),
        panel.grid = element_blank(), 
        legend.position = 'none')
p5.1reg

p5.2reg <- ggplot() + 
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 0.1), data = genes5.2_df, fill = colors2[2]) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -0.1, ymax = 0, fill = R2), data = ld5.2) +
  geom_hline(yintercept = 0) +
  geom_point(aes(chrom5.2site, 0), size = 3) +
  geom_text_repel(aes(x = midpt, y = 0, label = Name), data = genes5.2_df, direction = 'both', 
                  box.padding = 0.6, min.segment.length = 0, size = 3.881, max.overlaps = 100) +
  scale_x_continuous(limits = c(chrom5.2start, chrom5.2end)) +
  scale_y_continuous(limits = c(-0.35, 0.35)) +
  scale_fill_viridis(direction = -1, option = 'mako') +
  guides(fill = guide_colourbar(barwidth = 12,
                                barheight = 1)) +
  labs(x = NULL, y = NULL, fill = expression(LD~R^2)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        text = element_text(size = 11, color = 'black'),
        legend.text = element_text(size = 11, color = 'black', hjust = 0.5),
        line = element_line(color = 'black', linewidth = 1),
        panel.grid = element_blank(), 
        legend.position = 'none')
p5.2reg

p6reg <- ggplot() + 
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 0.1), data = genes6_df, fill = colors[3]) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -0.1, ymax = 0, fill = R2), data = ld6) +
  geom_hline(yintercept = 0) +
  geom_point(aes(chrom6site, 0), size = 3) +
  geom_text_repel(aes(x = midpt, y = 0, label = Name), data = genes6_df, direction = 'both', 
                  box.padding = 0.6, min.segment.length = 0, size = 3.881, max.overlaps = 100) +
  scale_x_continuous(limits = c(chrom6start, chrom6end)) +
  scale_y_continuous(limits = c(-0.35, 0.35)) +
  scale_fill_viridis(direction = -1, option = 'mako') +
  guides(fill = guide_colourbar(barwidth = 12,
                                barheight = 1)) +
  labs(x = NULL, y = NULL, fill = expression(LD~R^2)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        text = element_text(size = 11, color = 'black'),
        legend.text = element_text(size = 11, color = 'black', hjust = 0.5),
        line = element_line(color = 'black', linewidth = 1),
        panel.grid = element_blank(), 
        legend.position = 'none')
p6reg

fig6regions <- plot_grid(p5.1reg, p5.2reg, p6reg, nrow = 1, labels = c('F', 'G', 'H'))
fig6regions <- plot_grid(fig6regions, ld_legend, ncol = 1, labels = NULL, rel_heights = c(0.8, 0.2))
fig6regions

fig6 <- plot_grid(fig6top, fig6boxes, fig6regions, ncol = 1, labels = NULL, rel_heights = c(0.4, 0.3, 0.35))
fig6
ggsave('../Images/figure6.svg', plot = fig6, width = 6.5, height = 9.45, units = 'in', dpi = 300, bg = 'white')
