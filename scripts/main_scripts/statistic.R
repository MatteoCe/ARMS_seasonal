#!/usr/bin/Rscript

# load list of R packages used in the analyses

library(phyloseq)
library(tidyverse)
library(lattice)
library(data.table)
library(mgcv)
library(fitdistrplus)
library(corrplot)
library(olsrr)
library(ncf)
library(Metrics)
library(gratia)
library(RColorBrewer)
library(patchwork)
library(cowplot)
library(plotrix)

### DATASET PREPARATION

# load the otu table and sample data

count_tab <- read.table("clust_swarm_mod.counts.tsv",
                        header = TRUE,
                        row.names = 1,
                        check.names = FALSE)

sample_info <- read.table("sample_info_env.tsv",
                          header = TRUE,
                          row.names = 1,
                          check.names = FALSE,
                          sep = "\t",
                          stringsAsFactors = TRUE)

# Some samples may have been removed during the bioinformatic analyses. Thus,
# the sample info tab will be filetered to remove them and match those remaining
# in the count tab

sample_info_tab <- sample_info[rownames(sample_info) %in% colnames(count_tab), ]

# format the sample data as a "sampledata" phyloseq object

sampledata <- sample_data(sample_info_tab)

# format the otu table as a "OTU_table" phyloseq object

OTU <- otu_table(count_tab, taxa_are_rows = TRUE)

# merge the two objects into a single "phyloseq" object

physeq <- phyloseq(OTU, sampledata)

# calculate the number of mOTUs from phyloseq object as vectors and Shannon

rich <- estimate_richness(physeq, split = TRUE, measures = c("Observed", "Shannon"))

# store number of reads in a variable

total_reads <- sample_sums(physeq)

# create a dataframe with mOTU richness, all the environmental variables and 
# other useful factors

# set vector of most important field's names

fields <- c("arms_code",
            "id_site",
            "protection",
            "protection_level",
            "months",
            "depth",
            "lat",
            "lon",
            "chl_mean",
            "k490_mean",
            "par_mean",
            "npp_nsi",
            "npp_mean",
            "npp_total",
            "sss_mean",
            "sss_range",
            "sst_mean",
            "sst_range",
            "ECOREGION",
            "PROVINCE",
            "Lat_Zone")

fields_labels <- c("Sample Code",
                   "Random Effect group",
                   "Protection type",
                   "Protection level",
                   "# Months of deployment",
                   "Seafloor depth",
                   "Latitude",
                   "Longitude",
                   "Chlorophyll mean",
                   "Diffuse Attenuation Coef.",
                   "Photosynthetically Active Radiation",
                   "NSI",
                   "Net Primary Productivity mean",
                   "Total Net Primary Productivity",
                   "Sea Surface Salinity mean",
                   "Sea Surface Salinity range",
                   "Sea Surface Temperature mean",
                   "Sea Surface Temperature range",
                   "Ecoregion title",
                   "Province title",
                   "Latitudinal zone")

fields_df <- data.frame(field = fields, label = fields_labels)

# create a data.frame with the sample information uploaded earlier and formatted
# as phyloseq object, the total number of reads and the number of mOTUs. 

tab <- cbind(as.data.frame(rownames(rich)),
             sample_data(physeq)[, fields],
             total_reads,
             rich$Observed,
             rich$Shannon)

# Reset column names

colnames(tab) <- c("samples",
                   fields,
                   "total_reads",
                   "Observed",
                   "Shannon")

# change the measure unit of NPP variables by reporting the grams insteas of
# milligrams

tab$npp_mean <- tab$npp_mean / 1000
tab$npp_total <- tab$npp_total / 1000

### DATASET EXPLORATION

# Dotchart:

# check the presence of outliers for the response and explanatory variables
# print Supplementary Material S3.1
png("plots/S3.1_dotchart.png", units = "cm", width = 30, height = 30, res = 300)

op <- par(mfrow = c(5, 3), mar = c(3, 3, 3, 1))
dotchart(tab$Observed, main = "#OTUs")
dotchart(tab$chl_mean, main = "CHL - Chlorophyll mean")
dotchart(tab$k490_mean, main = "DA - Diffuse Attenuation coefficient mean")
dotchart(tab$par_mean, main = "PAR - Photosynthetically Active Radiation mean")
dotchart(tab$npp_nsi, main = "NSI - Normalized Seasonality Index")
dotchart(tab$npp_mean, main = "NPP_m - Net Primary Production mean")
dotchart(tab$npp_total, main = "NPP_t - Total Net Primary Production")
dotchart(tab$sss_mean, main = "SSS - Sea Surface Salinity mean")
dotchart(tab$sss_range, main = "SSS - Sea Surface Salinity range")
dotchart(tab$sst_mean, main = "SST - Sea Surface Temperature mean")
dotchart(tab$sst_range, main = "SST - Sea Surface Temperature range")
dotchart(tab$months, main = "Months of deployment")
dotchart(tab$depth, main = "Depth")
dotchart(abs(tab$lat), main = "Latitude")
dotchart(tab$total_reads, main = "Number of reads")
par(op)

dev.off()

# Define the variables relevant for model testing:

model_vars <- c("chl_mean",
                "k490_mean",
                "par_mean",
                "npp_nsi",
                "npp_mean",
                "npp_total",
                "sss_mean",
                "sss_range",
                "sst_mean",
                "sst_range",
                "months",
                "depth")

# Pearson Correlation #OTUs vs environmental variables:

# plot each environmental variable against the #OTUs and calculate Pearson's rho
# print Supplementary Material S3.2
png("plots/S3.2_Pearson_correlation.png",
                             units = "cm",
                             width = 25,
                             height = 30,
                             res = 300)

rows <- ceiling(length(model_vars) / 3)
cols <- length(model_vars) / rows

op <- par(mfrow = c(rows, cols))

for (var in model_vars) {

  pear_corr <- round(cor(tab[, var], tab$Observed, method = "pearson"), 2)
  plot(tab[, var], tab$Observed, 
       main = pear_corr,
       xlab = fields_df[fields_df$field == var, 2],
       ylab = "#OTUs")
  lines(lowess(tab[, var], tab$Observed), col = "red")

}

par(op)

dev.off()

# Distribution family of the response variable (Not included in the paper):

# the AIC is used to check which distribution has the lower akaike value if fit 
# to the real data, using the fitdistrplus package

distributions <- c("nbinom", "norm", "pois", "gamma", "geom")

# a vector with all the possible values of the response variables (from minimum
# values to maximum) is created and the density of different common
# distributions are computed using the main characteristics of those
# distributions, but calculated on the real observed dataset

# A distribution, with the characteristics of the real dataset, but with all 
# possible observations is created and can be compared to the histogram of the
# values in the real dataset
resp_var <- tab$Observed
name <- "# OTUs"

fits <- lapply(distributions, function(dist) {
  fitdist(resp_var, dist)
})

# compute the AIC for each distribution
AIC_values <- sapply(fits, function(fit) {
  k <- length(fit$estimate)
  logLik <- logLik(fit)
  2*k - 2*logLik
})

# combine the distribution names and AIC values into a data frame
results <- data.frame(distribution = distributions, AIC = AIC_values)

# before continuing the break parameter of the hist command can be checked and
# chosen on behalf of better representativity of the real distribution of the 
# data. Here 30 is set.
br <- 30

# set minimum and maximum values of the response variable chosen
min_set <- min(resp_var)
max_set <- max(resp_var)

# define a theoric distribution with all values from minimum and maximum of 
# those observed 
x <- min_set:max_set

# mean and variance
mean_tab <- mean(resp_var)
var_tab <- var(resp_var)

# k value for negative binomial
size_tab <- (mean_tab^2) / (var_tab - mean_tab)

# standard deviation for normal distribution
sdv_tab <- sqrt(var_tab)

# for gamma distribution
shape_tab <- (mean_tab^2) / var_tab

# calculate density for a negative binomial distribution
y_nbinom <- dnbinom(x, size = size_tab, mu = mean_tab)

# calculate density for a normal distribution
y_norm <- dnorm(x, mean = mean_tab, sd = sdv_tab)

# calculate density for a Poisson distribution
y_pois <- dpois(x, lambda = mean_tab)

# calculate density for a Gamma distribution
y_gamma <- dgamma(x, shape_tab)

# calculate density for a Geometric distribution (k value of Geometric 
# distribution = 1 of Negative binomial and variance as quadratic function of 
# the mean)
y_geom <- dnbinom(x, size = 1, mu = mean_tab)

# plot the response variable with all possible theoretic distributions
png("plots/Distribution Families.png",
                  units = "cm",
                  width = 30,
                  height = 30,
                  res = 300)

op <- par(mfrow = c(3, 2))

# plot the histogram of the real response variable dataset
hist(resp_var, breaks = br, main = name)

plot(x, y_nbinom, main = "Density Plot of Negative Binomial Distribution",
     type = "h",
     lwd = 2,
     xlab = "x",
     ylab = "Density")
text(400, .004, paste0("AIC=", round(results[1, 2], 3)), cex = 1.5)

plot(x, y_norm, main = "Density Plot of Normal Distribution",
	 type = "h",
	 lwd = 2,
	 xlab = "x",
	 ylab = "Density")
text(400, .0025, paste0("AIC=", round(results[2, 2], 3)), cex = 1.5)

plot(x, y_pois, main = "Density Plot of Poisson Distribution",
	 type = "h",
	 lwd = 2,
	 xlab = "x",
	 ylab = "Density")
text(400, .025, paste0("AIC=", round(results[3, 2], 3)), cex = 1.5)

plot(x, y_gamma, main = "Density Plot of Gamma Distribution",
	 type = "h",
	 lwd = 2,
	 xlab = "x",
	 ylab = "Density")
text(400, .000000125, paste0("AIC=", round(results[4, 2], 3)), cex = 1.5)

plot(x, y_geom, main = "Density Plot of Geometric Distribution",
	 type = "h",
	 lwd = 2,
	 xlab = "x",
	 ylab = "Density")
text(400, .0045, paste0("AIC=", round(results[5, 2], 3)), cex = 1.5)

par(op)

dev.off()

# Boxplots of number of OTUs index grouped by Ecoregion Spalding et al. 2007:
# Part of Figure 1

tab2 <- tab %>% group_by(ECOREGION) %>% mutate(Obs_median = median(Observed))
tab2_ord <- tab2[order(tab2$Obs_median), ]
tab2_ord$ECOREGION <- factor(tab2_ord$ECOREGION, 
                             levels = unique(tab2_ord$ECOREGION))

# Simple plot with only OTU number
p1_simp <- ggplot(tab2_ord) +
                 geom_boxplot(aes(ECOREGION, Observed, fill = Lat_Zone)) +
                 theme_half_open() +
                 ggtitle("# OTUs") +
                 theme(axis.text.x = element_text(angle = 45,
                                                  hjust = 1,
                                                  vjust = 1,
                                                  size = 11,
                                                  face = "bold"),
                       plot.title = element_text(hjust = 0.5,
                                                 face = "bold"),
                       panel.grid.major.x = element_blank()) +
                 theme(legend.position = "none") +
                 xlab(NULL) +
                 ylab(NULL) +
                 scale_y_continuous(position = "right") +
                 scale_fill_manual(values = c("#98D3FF",
                                              "#83EC87",
                                            "#FF9898"),
                                   name = ("Latitudinal \n    zone"))

p1_simp_fin <- cowplot::plot_grid(plot.new(), p1_simp, plot.new(),
                                  align = "vh",
                                  axis = "t",
                                  ncol = 3,
                                  nrow = 1,
                                  rel_heights = c(1, 1, 1),
                                  rel_widths = c(0.16, 1, 0.01))

# print Part of Figure 1
ggsave("plots/Figure1_OTU_Ecoregion.png",
                           plot = p1_simp_fin,
                           device = png(),
                           path = NULL,
                           scale = 1,
                           width = 7,
                           height = 10,
                           dpi = 300,
                           limitsize = TRUE)

dev.off()

# Schematic plot of seasonality dynamics for map:
# Part of Figure 1

t <- seq(0, 6, .1)
y <- sin(t)
y1 <- rep(.5, 61)
y2 <- y + (abs(min(y)) + .2)
y3 <- dpois(seq(0:60), lambda = 7)
seas <- data.frame(t, y1, y2, y3)

p1 <- ggplot(seas,aes(t, y1)) +
             theme_classic() +
             theme(axis.title = element_text(size = 13,
                                             face = "bold"),
                                             axis.text.y = element_blank(),
                                             axis.ticks.y = element_blank()) +
             ylab("Net Primary Production") +
             xlab("") +
             geom_area(fill = "#e8e1ef",
                       color = "black") +
             scale_y_continuous(limits = c(0, 2.5)) +
             scale_x_continuous(breaks = seq(0, 6, .5),
                                labels = NULL)
p2 <- ggplot(seas, aes(t, y2)) +
             theme_classic() +
             theme(axis.title = element_text(size = 13, face = "bold"),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank()) +
             ylab("") +
             xlab("Time of year") +
             geom_area(fill = "#a486c1", color = "black") +
             scale_y_continuous(limits = c(0, 3.5)) +
             scale_x_continuous(breaks = seq(0, 6, .5),
                                labels = NULL)

p3 <- ggplot(seas, aes(t, y3)) +
             theme_classic() +
             theme(axis.title = element_text(size = 13,
                                             face="bold"),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank()) +
             ylab("") +
             xlab("") +
             geom_area(fill = "#3f007d",
                     color = "black") +
             scale_y_continuous(limits = c(0, .15)) +
             scale_x_continuous(breaks = seq(0, 6, .5),
                                labels = NULL)

seasonality <- cowplot::plot_grid(p1, p2, p3,
                                  align = "vh",
                                  axis = "t",
                                  ncol = 3,
                                  nrow = 1,
                                  rel_heights = c(1, 1, 1),
                                  rel_widths = c(1, 1, 1))

# print Part of Figure 1
ggsave("plots/Figure1_Schematic_Seasonality.png", plot = seasonality,
                                    device = png(),
                                    path = NULL,
                                    scale = 3,
                                    width = 3,
                                    height = 1,
                                    dpi = 200,
                                    limitsize = TRUE)

dev.off()

# Collinearity testing:

# Due to the problematic of including highly collinear variables in GAM models,
# the uncollinear variables are selected by calculating the Variance Inflating 
# Factor in a linear model including all of them, and removing the variable with
# higher VIF, until all variables in the model have a VIF < 3. 

model_vars <- c("chl_mean",
                "k490_mean",
                "par_mean",
                "npp_nsi",
                "npp_mean",
                "npp_total",
                "sss_mean",
                "sss_range",
                "sst_mean",
                "sst_range",
                "months",
                "depth")

tab_corr_mat <- cor(tab[, model_vars])

png("plots/S3.3_Correlation_plot.png", units = "cm",
                                             width = 30,
                                             height = 30,
                                             res = 300)
corrplot(tab_corr_mat, method = "number",
                       type = "lower",
                       tl.col = "black")

dev.off()

# testing variance inflating factors, each variable showing the highest score of
# VIFs is removed and the model refitted without it

# The process is performed untill all variables have a VIF < 3

model_vars_mod <- model_vars

for (n in seq(1:length(model_vars_mod))) {

  model_vars_fun <- paste("Observed ~", paste(model_vars_mod, collapse = "+"))

  lm <- lm(model_vars_fun, data = tab)

  lm_data <- ols_coll_diag(lm)
  lm_data_df <- lm_data$vif_t %>% dplyr::filter(!(Variables %in% c("npp_nsi", 
                                                                   "npp_total")))

  lm_data_df_sel <- lm_data_df[lm_data_df$VIF == max(lm_data_df$VIF),]

  if (lm_data_df_sel$VIF <= 3) {

    # print part of values from Table 2
    file <- file("prints/VIFs_final.txt", "w")
    sink(file)
    print(lm_data)
    sink()
    close(file)
    sink()

    vif_tolerance_df <- lm_data$vif_t

    Sys.sleep(1)

    break

  } else {

    var_rem <- lm_data_df_sel$Variables

    model_vars_mod <- model_vars_mod[!(model_vars_mod %in% var_rem)]

  }

}

# Plot each environmental variable chosen by the VIF selection process with the
# #OTUs and calculate the Pearson's rho

model_vars_labels <- c("Chlorophyll mean",
                       "Diffuse Attenuation Coef.",
                       "Photosynthetically Active Radiation",
                       "NSI",
                       "Net Primary Productivity mean",
                       "Total Net Primary Productivity",
                       "Sea Surface Salinity mean",
                       "Sea Surface Salinity range",
                       "Sea Surface Temperature mean",
                       "Sea Surface Temperature range",
                       "# Months of deployment",
                       "Seafloor depth")

model_vars_df <- data.frame(var = model_vars, label = model_vars_labels)

vars_df <- model_vars_df[model_vars_df$var %in% model_vars_mod, ]

png("plots/Pearson_correlation_uncollinear.png", 
                          units = "cm",
                          width = 30,
                          height = 30,
                          res = 300)

rows <- ceiling(length(model_vars_mod)/2)
cols <- length(model_vars_mod) - rows

op <- par(mfrow = c(rows, cols))

for (var in model_vars_mod) {

  pear_corr <- round(cor(tab[, var],
                         tab$Observed,
                         method = "pearson"),
                         2)
  plot(tab[, var], tab$Observed,
       main = pear_corr,
       xlab = vars_df[vars_df$var == var, 2],
       ylab = "#OTUs")
  lines(lowess(tab[, var], tab$Observed),
        col = "red")

}

par(op)

dev.off()

# Spatial autocorrelation:

# The presence of spatial auto-correlation will be inspected with a plot, using
# the package "ncf", that will show what type of correlation (or not) is present
# at differing distances between the samples

corr_spaz <- spline.correlog(x = tab[, "lon"],
                             y = tab[, "lat"],
                             z = tab[, "Observed"],
                             latlon = TRUE)

png("plots/Spatial autocorrelation Dataset.png", units = "cm",
                                           width = 30,
                                           height = 20,
                                           res = 300)
plot(corr_spaz)

dev.off()

# The wiggly form of the plot at higher distances reflects the sparsely
# distributed dataset we have, with most of the samples from the Mediterranean,
# Red and Northern Seas, and sparse samples in the pacific and antarctica. This
# shows more negative and positive correlation without a clear sense, but the
# most important part is the high correlation at reduced distances. This is
# mostly due to the fact that for each X and Y (latitude and longitude
# combinations) we have more than one sample, often in three replicates. This
# can be inspected by looking at the factor level "site", based on the lat lon
# and the number of ARMS deployed in that location in a particular year

# In order to account for this auto-correlation, a random effect using that
# factor can be used. The effectiveness of the random effect can be assessed by
# checking the AIC of the model without it with the one that has it, and then
# comparing the autocorrelation left in the pearson residuals of the two models.
# The model, including the variables our main hypothesis, will be used as a 
# fixed effect.

GAMM_norand <- gam(Observed ~ s(npp_nsi) + s(npp_total),
                   method = "REML",
                   family = nb(theta = NULL,
                               link = "log"),
                               data = tab)

GAMM_rand <- gam(Observed ~ s(npp_nsi) + s(npp_total) +
                 s(id_site, bs = "re"),
                 method = "REML",
                 family = nb(theta = NULL,
                             link = "log"),
                             data = tab)

corr_spaz_norand <- spline.correlog(x = tab[, "lon"],
                                    y = tab[, "lat"],
                                    z = resid(GAMM_norand,
                                              type = "pearson"),
                                              latlon = TRUE)

corr_spaz_rand <- spline.correlog(x = tab[, "lon"],
                                  y = tab[, "lat"],
                                  z = resid(GAMM_rand,
                                            type = "pearson"),
                                            latlon = TRUE)

sink("prints/AIC_randomEffect.txt")
AIC(GAMM_norand, GAMM_rand)
sink()

# print Supplementary Material S3.4
png("plots/S3.4_SpatialAutocorrelation_RandomEffect.png", units = "cm",
                                                            width = 30,
                                                            height = 30,
                                                            res = 300)

op <- par(mfrow = c(2, 2))

plot(corr_spaz_norand, main = "Spatial autocorrelation without Random Effect")
text(x = 16500, y = .85, cex = 1.5, labels = paste("AIC:", round(AIC(GAMM_norand), 3)))

lin.pred_norand <- napredict(GAMM_norand$na.action,
                             GAMM_norand$linear.predictors)
plot(lin.pred_norand, resid(GAMM_norand), 
                      ylab = ("residuals"),
                      xlab = ("Linear predictor"),
                      main = "Residuals vs. Linear predictor \n without Random Effect")
lines(lowess(lin.pred_norand, resid(GAMM_norand)), col = "red")

plot(corr_spaz_rand, main = "Spatial autocorrelation with Random Effect")
text(x = 16500, y = .85, cex = 1.5, labels = paste("AIC:", round(AIC(GAMM_rand), 3)))

lin.pred_rand <- napredict(GAMM_rand$na.action, 
                           GAMM_rand$linear.predictors)
plot(lin.pred_rand, resid(GAMM_rand),
                    ylab = ("residuals"),
                    xlab = ("Linear predictor"),
                    main = "Residuals vs. Linear predictor \n with Random Effect")
lines(lowess(lin.pred_rand, resid(GAMM_rand)), col = "red")

par(op)

dev.off()

# Thus, the base model will include a smoother for both total primary 
# productivity and the normalized seasonality index, plus a random effect 
# smoother for the site factor

### Model selection

# A backword selection is performed testing which explanatory variables, of
# those with VIFs lower than 3, can be included in the model. However, all of 
# the additional uncollinear variables resulted statistically significant at the
# first round of the selection, although with p-values relatively high.  

# Another model will be computed including exclusively the smoothers for the NSI
# and total NPP.

# The fitness of the model will be evaluated by comparic the AIC with the base
# model, together with the info from the gam.check command. This will be
# computed using the gam command of the mgcv package following Pedersen et al
# 2019.

GAMM0 <- gam(Observed ~ s(npp_nsi) +
                        s(npp_total) +
                        s(months) +
                        s(depth) +
                        s(sss_range) +
                        s(id_site, bs = "re"),
                                   method = "REML",
                                   family = nb(theta = NULL,
                                               link = "log"),
                                   data = tab)

# print part of values from Table 2
sink("prints/Table2_Summary statistics HGAM model_uncollinear.txt")
summary(GAMM0)
sink()

# fit the model with only the variables of the main hypothesis: NSI and total 
# NPP

GAMM1 <- gam(Observed ~ s(npp_nsi) +
                        s(npp_total) +
                        s(id_site, bs = "re"),
                                   method = "REML",
                                   family = nb(theta = NULL,
                                               link = "log"),
                                   data = tab)

# print part of values from Table 2
sink("prints/Table2_Summary statistics HGAM model_hypothesis.txt")
summary(GAMM1)
sink()

# print Supplementary Material S3.5
png("plots/S3.5_NSI-TotalNPP uncollinear smoothers plot.png", units = "cm",
                                    width = 20,
                                    height = 15,
                                    res = 300)
draw(GAMM0)
dev.off()

# print Supplementary Material S3.4
png("plots/S3.6_NSI-TotalNPP mgcv diagnostic.png", units = "cm",
                                                          width = 20,
                                                          height = 20,
                                                          res = 300)
gam.check(GAMM1)
dev.off()

# Calculate the AIC values of each model

models_list <- list()

models_list[[1]] <- GAMM0
models_list[[2]] <- GAMM1

coefs <- list()

for (i in seq(1, 2)) {

	model <- models_list[[i]]

	sum <- summary(model)

	edf <- sum(round(sum$edf,2))
	dev.expl <- round(sum$dev.expl*100,2)
	aic <- AIC(model)
	par_names <- c("edfs","Dev. expl. %","AIC")
	pars_gam <- cbind(edf,dev.expl,aic)
	colnames(pars_gam) <- par_names

	coefs[[i]] <- pars_gam


}

param_full_models <- as.data.frame(do.call(rbind,coefs))
prev_cols <- colnames(param_full_models)

rownames(param_full_models) <- c("GAMM0", "GAMM1")

param_full_models$FixEff <- c("NSI + NPP total + All non-collinear","NSI + NPP total")

# print part of values from Table 3
sink("prints/Table3_Final_models_AIC.txt")
param_full_models
sink()

# plot smoothers

# create dummy data.frame for the legend, and use the same limit for the colors
# colors of the linear predictor in the next plot
df <- data.frame(x = seq(1:100), 
                 value = c(50, sample(x=seq(51, 299), size = 98), 300))

# create fake plot for the legend
pleg2 <- ggplot(df, aes(x = x, y)) +
         geom_bar(aes(fill = value), stat = "identity") + 
         scale_fill_gradient(low = c(sort(heat.colors(2), 
                                           decreasing = TRUE)[2]),
                                high = c(sort(heat.colors(2), 
                                           decreasing = TRUE)[1]),
                           name = "OTU\nnumber\n",
                           labels = seq(from = 50, to = 300, by = 50),
                           n.breaks = length(seq(from = 50, to = 300, by = 50)),
                           limits = c(50, 300))


# extract the legend
leg2 <- ggpubr::get_legend(pleg2)

GAMM0_models_effects <- wrap_elements(panel = 

~ {

  op <- par(mfrow = c(1, 1), mar = c(5, 5, 3, 1))

  vis.gam(GAMM0, type = "response",
                 plot.type = "contour",
                 main = "NSI + NPP total + All non-collinear",
                 xlab = "NSI",
                 ylab = "Total NPP",
                 zlim = c(0, 350))
  par(op)

})

GAMM1_models_effects <- wrap_elements(panel = 

~ {

  op <- par(mfrow = c(1, 1), mar = c(5, 5, 3, 1))

  vis.gam(GAMM1, type = "response",
                 plot.type = "contour",
                 main = "NSI + NPP total",
                 xlab = "NSI",
                 ylab = "",
                 zlim = c(0, 350),
                 yaxt = "n")

  par(op)

})

dev.off()

layout <- c(area(1, 1, 3, 8),
            area(1, 9, 3, 16),
            area(2, 17, 2, 18))

models_effects <- GAMM0_models_effects + 
                  GAMM1_models_effects + 
                  leg2 + 
                  plot_layout(design = layout)

# print Figure 2

ggsave("plots/Figure2_NSI_NPP models effects.jpg", plot = models_effects, 
                                                   path = NULL, 
                                                   width = 12, 
                                                   height = 6,
                                                   limitsize = TRUE)

# Different GAMs will be fitted, including a single variable between those in 
# the dataset that are usually included as explanatory variables of richness
# such as Sea Surface Temperature, Chlorophyll and Net primary productivity alone

coefs <- list()

mod_res <- data.frame(vars = c("npp_nsi",
                               "npp_total",
                               "sst_mean",
                               "chl_mean",
                               "npp_mean"),
                      names = c("NSI",
                                "Total NPP",
                                "SST mean",
                                "CHL mean",
                                "NPP mean"))

mod_expr <- list("1" = "degree of seasonality",
                 "2" = expression("g" %*% "m" ^-2 %*% "year" ^-1),
                 "3" = expression("C" ^o),
                 "4" = expression("mg" %*% "m" ^3),
                 "5" = expression("g" %*% "m" ^-2 %*% "day" ^-1))

# print Figure 3
png("plots/Figure3_Alternative GAMs smoothers.png", units = "cm",
                                      width = 25,
                                      height = 15,
                                      res = 300)
op <- par(mfrow = c(2, 3))

for (i in seq(1, nrow(mod_res))) {

    tab_name <- mod_res[i, ]$vars

    formula <- substitute(Observed ~ s(x) + 
                                     s(id_site, bs = "re"),
                          list(x = as.name(tab_name)))

    model <- gam(formula, family = nb(theta = NULL, link = "log"),
                          method = "REML",
                          data = tab)

    sum <- summary(model)

    # extract estimated degrees of freedom and p value for the fixed effect
    chi.sq <- round(sum$chi.sq, 2)[1]
    pv <- round(sum$s.table[1, 4], 3)

    if (pv == 0) {pv <- "<2e-16 ***"
        } else if (pv > 0.01 && pv < 0.5) {
          paste(pv, "*")
        } else {}

    par_names <- c("chi-square", "p-value")
    pars_gam <- cbind(chi.sq, pv)
    colnames(pars_gam) <- par_names

    coefs[[i]] <- pars_gam

    if (i == 1 || i == 4) {

      plot(model, 
         select = 1,
         main = mod_res[i, ]$names,
         xlab = mod_expr[[i]],
         ylab = "Partial effect on linear predictor")

    } else {

      plot(model, 
           select = 1,
           main = mod_res[i, ]$names,
           xlab = mod_expr[[i]],
           ylab = "")

    }

}

grobTable <- do.call(rbind, coefs)
rownames(grobTable) <- mod_res$names

plot.new()
addtable2plot(0, 0, grobTable,
              xjust = .1,
              yjust = 1,
              bty = "n",
              display.rownames = TRUE,
              hlines = TRUE,
              vlines = FALSE,
              cex = 1.2,
              ypad = 2,
              xpad = .1)

par(op)
dev.off()

# Cross validation

# The prediction performance of the main model is compared with HGAM models 
# including environmental variables commonly used in literature to describe 
# benthic metazoan richness following a cross-validation approach. For each
# ecoregion in Spalding et al., two sites are randomly selected and all the
# different models are fitted on this "training" dataset. The obtained model is 
# used to calculate the prediction of the response variable on the validation 
# dataset, which is obtained by selecting the samples that are not included in 
# the training dataset. The proportion of training vs validation datasets sample
# is approximately 50% vs 50%, with a certain variability given by the uneven
# number of replicates in each site. However, this allows to keep approximately
# the same number of sites in different ecoregions of the world, and limit the
# weight of highly sampled regions in the dataset (e.g. Red Sea)

# The prediction performance is evaluated by calculating the RMSE (Root Mean Squared Error)
# between the actual and the predicted value of number of OTUs, using the package
# Metrics. The linearity of NPP total, included in the model with NSI, if assessed
# graphycally by temporarily modifying the value of NSI in eahc validation dataset to
# different levels of seasonality (from 0.3 to 0.8) and plotting each prediction 
# as a line.

# The random selection of vsites is performed 100 times, and the id_sites 
# selected for each trial is stored in a tsv table (Supplementary Material 2.4)

# Cross-validation parameters

# The number of trials
ntrial <- 100
# The region onto which perform the selection of sites
region <- "ECOREGION"
# The number of sites to sample in each region previously defined
sample_size <- 2

# This provides the number of sites that will be retained each trial, allowing
# to construct the data.frame that will store the sites names
training_set_size <- tab %>% 
                     group_by(.[[region]], id_site) %>% 
                     summarise() %>% 
                     slice_sample(n = sample_size) %>% 
                     nrow()

# Create empty data.frame to store the sites names
training_sets <- data.frame(matrix(ncol = ntrial, 
                                   nrow = training_set_size), 
                            row.names = seq(1:training_set_size))

colnames(training_sets) <- seq(1, 100)

# Create a data.frame including model name, the ame of the variables included in 
# each model and the formula used for the HGAM model
models_info <- data.frame(name = c("NSI",
                                   "CHL mean",
                                   "NPP mean",
                                   "Total NPP",
                                   "SST mean",
                                   "NSI + NPP Total"),
                          variables = c("npp_nsi", 
                                        "chl_mean",
                                        "npp_mean",
                                        "npp_total",
                                        "sst_mean",
                                        "npp_nsi|npp_total"),
                          formula = c("Observed ~ s(npp_nsi) + s(id_site, bs = 're')",
                                      "Observed ~ s(chl_mean) + s(id_site, bs = 're')",
                                      "Observed ~ s(npp_mean) + s(id_site, bs = 're')",
                                      "Observed ~ s(npp_total) + s(id_site, bs = 're')",
                                      "Observed ~ s(sst_mean) + s(id_site, bs = 're')",
                                      "Observed ~ s(npp_nsi) + s(npp_total) + s(id_site, bs = 're')"))

# To produce the RMSE plots, store the error values in a data.frame
rmse_df <- data.frame(matrix(ncol = ntrial, 
                             nrow = nrow(models_info)),
                      row.names = models_info$name)

# Now run the cross-validation.

# two empty lists are created to store the predictions of the model including 
# only NPP Total, and another to store the predictions of the model "NSI + NPP total"
# at fixed levels of seasonality
preds_NPP_total <- list()
preds_fixed_nsi_GAMM1 <- list()

for (i in 1:ntrial) {

    # create training and validation datasets (for the training dataset, 
    # 2 sites for each region (approximately 1/2 of original number of samples):
    suppressMessages(

    # Randomply select 2 sites for each region
    sites_sel <- tab %>% group_by(.[[region]], id_site) %>% 
                         summarise() %>% 
                         slice_sample(n = sample_size) %>% 
                         ungroup() %>% 
                         dplyr::select(id_site) %>% 
                         c(.) %>% unlist() %>% unname() %>% as.character()

    )

    # Assign the data to the correct sets
    training <- tab %>% dplyr::filter(id_site %in% sites_sel)
    validation <- tab %>% dplyr::filter(!(samples %in% training[, "samples", drop = TRUE]))

    # register sites names of training set in a dataframe
    training_sets[, i] <- as.character(unique(training$id_site))

    # Use the training dataset for the RMSE calculation, after fitting the HGAM
    # model. This is done in a loop for each model in models_info
    for (y in seq(1, nrow(models_info))) {

      # get informations for the current model in distinct objects
      name <- models_info[y, "name"]
      variables <- stringr::str_split(models_info[y, "variables"], "\\|") %>% unlist()
      formula <- models_info[y, "formula"]

      # fit the model
      GAMM <- gam(as.formula(formula),
                  method = "REML",
                  family = nb(theta = NULL,
                              link = "log"),
                  data = training)
      
      # use the model to predict the OTU number from the validation dataset
      prediction <- c(predict.gam(GAMM, 
                      newdata = validation[, c(variables, "id_site"), 
                                           drop = FALSE],
                      newdata.guaranteed = TRUE,
                      se.fit = TRUE,
                      type = c("response")))

      # calculate the RMSE
      rmse_df[name, i] <- rmse(validation$Observed, prediction$fit)

      # If the model is the one of the main hypothesis, run the same prediction
      # on a validation dataset with fixed values of NSI 
      if (name == "NSI + NPP Total") {

        # fixed NSI smoother for NPP mean
        for (j in seq(.3, .8, by = .1)) {

          # create "copy" of the validation dataset
          validation_nsi_mod <- validation

          # fix the values for the seasonality index of the validation dataset 
          # to the values chosen by the loop arguments (from 0.3 to 0.8). These
          # "moderate" values were chosen to restrain the prediction of the
          # model to less extreme values, not necessarily more represented
          validation_nsi_mod$npp_nsi <- j

          prediction <- c(predict.gam(GAMM, 
                                      newdata = validation_nsi_mod[, c(variables, "id_site"), 
                                                                   drop = FALSE],
                                      newdata.guaranteed = TRUE,
                                      se.fit = TRUE,
                                      type = c("response")))

          res_fix <- data.frame(prediction$fit, prediction$se.fit,
                                validation_nsi_mod[, variables, drop = FALSE])

          preds_fixed_nsi_GAMM1[[paste0(i, "_", j)]] <- res_fix

          }

        } else if (name == "Total NPP") {

          preds_NPP_total[[i]] <- data.frame(prediction$fit, 
                                             prediction$se.fit,
                                             validation[, variables, drop = FALSE])

        }

    }

}

# store the id_site names for each training set as a tsv file for reproducibility 
# of figures. Validation datasets can be recreated by selecting the id_site in 
# the original tab that were not included in each training set
write.table(training_sets, file = "prints/Cross-validation_training-sets.tsv",
                         col.names = NA,
                         row.names = TRUE,
                         quote = FALSE,
                         sep = "\t",
                         fileEncoding = "UTF-8")

# Store the RMSE values in a two column data.frame with the RMSE for each trial 
# and the corresponding HGAM model name as factor variable
rmse_tidy <- stack(as.data.frame(t(rmse_df)))
colnames(rmse_tidy) <- c("RMSE", "MODEL")

# For each model, calculate the median of the RMSE value and sort the model by
# decreasing RMSE median
rmse_tidy <- rmse_tidy %>% 
             group_by(MODEL) %>% 
             mutate(med = median(RMSE))
rmse_tidy_ord <- rmse_tidy[order(-rmse_tidy$med), ]
rmse_tidy_ord$MODEL <- factor(rmse_tidy_ord$MODEL, 
                              levels = unique(rmse_tidy_ord$MODEL))

# print Violin plots for the RMSE of each trial in each MODEL
violins <- ggplot(rmse_tidy_ord) +
                  geom_violin(aes(MODEL, RMSE), 
                              fill = '#A4A4A4', 
                              alpha = .5) +
                  scale_y_continuous(limits = quantile(rmse_tidy_ord$RMSE, c(.01, .99))) +
                  theme_minimal() +
                  xlab("Fixed effects") +
                  ylab("RMSE") +
                  geom_boxplot(aes(MODEL, RMSE), width = .05) +
                  theme(panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        plot.background = element_rect(fill = "white",
                                                       color = "white"))

# Save Figure 4
ggsave("plots/Figure4_RMSE models.png", 
       plot = violins, 
       device = png(), 
       path = NULL, 
       scale = 1, 
       width = 9, 
       height = 6, 
       dpi = 300, 
       limitsize = TRUE)

# Now run the smooth plots showing consistency of linearity of NPP total in the
# cross validation process. First, rbind all dataframes stored in the lists with
# predictions for NPP total and NPP total + NSI at fixed seasonality levels

names(preds_NPP_total) <- c(1:length(preds_NPP_total))
NPP_total_smooth <- data.frame(rbindlist(preds_NPP_total, 
                                         idcol = "name"))

names(preds_fixed_nsi_GAMM1) <- c(1:length(preds_fixed_nsi_GAMM1))
GAMM1_fix_nsi_smooth <- data.frame(rbindlist(preds_fixed_nsi_GAMM1, 
                                             idcol = "name"))

# Foth the data conserning the model with total NPP only, create two new 
# variables that store the approximate 95% confidence interval for each prediction
# This will allow geom_smooth to plot a gam smooth of all confidence intervals
# for all trials
NPP_total_smooth$uprP <- NPP_total_smooth$prediction.fit + (2 * NPP_total_smooth$prediction.se.fit)
NPP_total_smooth$lwrP <- NPP_total_smooth$prediction.fit - (2 * NPP_total_smooth$prediction.se.fit)

# Plot the predictions for each trial of the NPP total model and the gam smooth
# of all confidence intervals
NPP_total_predict <- ggplot(data = NPP_total_smooth) +
                     theme_classic() +
                     geom_smooth(aes(npp_total, prediction.fit,
                                     group = name),
                                method = gam,
                                formula = y ~ s(x),
                                se = FALSE,
                                color = "grey",
                                linewidth = .5) +
                     geom_smooth(aes(npp_total, uprP),
                                method = gam,
                                formula = y ~ s(x),
                                se = FALSE,
                                color = "black",
                                linewidth = .5,
                                linetype = "dashed") +
                     geom_smooth(aes(npp_total, lwrP),
                                method = gam,
                                formula = y ~ s(x),
                                se = FALSE,
                                color = "black",
                                linewidth = .5,
                                linetype = "dashed") +
                     ylab("Predicted #OTUs") +
                     xlab("Total NPP") +
                     labs(title = models_info[4, "name"]) +
                     theme(plot.title = element_text(size = 11)) +
                     scale_y_continuous(position = "left",
                                        limits = c(0, 500))

# Similarly, for the hypothesis HGAM, plot all predictions at fixed levels of
# seasonality
GAMM1_fix_nsi_predict <- ggplot(data = GAMM1_fix_nsi_smooth) +
                         theme_classic() +
                         geom_smooth(aes(npp_total, prediction.fit,
                                         colour = factor(npp_nsi),
                                         group = name),
                                    method = gam, 
                                    formula = y ~ s(x),
                                    se = FALSE,
                                    linewidth = .05,
                                    alpha = .95) +
                         geom_line(aes(npp_total, prediction.fit,
                                        group = factor(npp_nsi)),
                                    colour = "black",
                                    alpha = 1,
                                    stat = "smooth",
                                    method = gam,
                                    formula = y ~ s(x),
                                    linewidth = 1.5,
                                    linetype = "solid") +
                         geom_smooth(aes(npp_total, prediction.fit,
                                        color = factor(npp_nsi),
                                        fill = factor(npp_nsi)),
                                    method = gam,
                                    formula = y ~ s(x),
                                    se = FALSE,
                                    linewidth = 1,
                                    linetype = "solid") +
                         ylab("") +
                         xlab("Total NPP") +
                         labs(title = models_info[6, "name"]) +
                         theme(legend.position = "none",
                               axis.text.y = element_blank(),
                               plot.title = element_text(size = 11)) +
                         scale_y_continuous(position = "left",
                                            limits = c(0, 500))

# The legend will be created from a dummy dataset that allows to change the
# legend shapes to square instead of lines, for enhanced visibility
legend_df <- GAMM1_fix_nsi_smooth %>% 
             dplyr::group_by(npp_nsi) %>% 
             slice_sample(n = 1) %>% 
             dplyr::ungroup()

pleg <- ggplot(legend_df) + 
        geom_bar(aes(x = npp_total, 
                     fill = factor(npp_nsi)),
                     colour = "black") +
        scale_colour_manual(values = c("#F8766D", 
                                       "#B79F00", 
                                       "#00BA38", 
                                       "#00BFC4", 
                                       "#619CFF", 
                                       "#F564E3")) +
        theme_minimal() +
        guides(fill = guide_legend(title = "Fixed NSI \n   value"))

# extract the legend and combine all models in a single one with cowplot
leg <- ggpubr::get_legend(pleg)

# Included plot.new on first column to add more margin for the tilted x labels
smooth_plots <- cowplot::plot_grid(NPP_total_predict, 
                                   GAMM1_fix_nsi_predict, 
                                   leg,
                                   align = "hv",
                                   axis = "b",
                                   ncol = 3,
                                   nrow = 1,
                                   rel_heights = c(1, 1, 1),
                                   rel_widths = c(1, 1, .3))

# print Figure 5
ggsave("plots/Figure5_NPP predictions.png", 
       plot = smooth_plots, 
       device = png(), 
       path = NULL, 
       scale = 1, 
       width = 11, 
       height = 6, 
       dpi = 300, 
       limitsize = TRUE)

# Missing rich regions

# The completeness of the dataset is given by the diverse environments 
# investigated by the different projects that conducted the experiments in very
# different areas of the world coastal regions. This is confirmed by the high 
# dataset variability in both seasonality index (spanning from 0.3 to 0.9) and 
# total NPP, given by the presence of samples obtained from regions with eutrophic
# dynamism (Baltic Sea), tropical environments (Red sea, Hawaii), polar and 
# temperate areas. Nonetheless, other, peculiar environments go inevitably
# missed in the context of global ecological studies, especially in the case of 
# new emerging techniques and sampling methods. This issue is here represented 
# by the reduced representativity of "low seasonality - high productivity" 
# environments, which, are notably less distributed in the world compared to 
# other combinations of NSI and magnitude. This is, of course, to consider 
# keeping in mind the range of NPP total and seasonality values here observed, 
# i.e, the minimum and maximum values of both variables, which should not be considered 
# as distinct variables, but intrinsecally connected, to the point where, for
# example, specific extreme values of NPP magnitude, could not, or very rarely,
# be observed concurrently with specific degress of seasonality. For example,
# NPP magnitude tends to reach very high values in highly seasonal environments,
# but the same extreme values are rarely, or never, observed in low seasonal 
# environments. The relationship between the two is thus relatively strict, 
# due to the fenology of planktonic communities. Nonetheless, a certain 
# variability of NPP magnitude is still observed, at varying degrees of 
# seasonality, with a reduced range of NPP magnitude values at lower 
# seasonality, in respect to areas with higher seasonality.

# The model, as shown in Figure 2, shows a negative relationship between
# NPP magnitude or NSI and #OTUs. The concurrent effect of both variables defines
# the combination of environmental variables producing the highest number of OTUs
# as the "low seasonality - low productivity", with decreasing OTUs at increasing
# seasonality and NPP magnitude. However, as mentioned before, the absence of 
# high confidence in the "low seasonality - high productivity" plot area, given
# by the absence of samples from those regions, impede us from confirming at the
# statistical level the consistent negative relationship between magnitude and
# OTU number, also recurring at Figure 5 in almost all cross validation trials.

model_95confidence <- wrap_elements(panel = 

~ {

  op <- par(mfrow = c(1, 1), 
            mar = c(1.4, 2.5, 1.8, .46), 
            mgp = c(2, .6, .1), 
            bty = "n")

  vis.gam(GAMM1, 
          type = "response", 
          too.far = 0.05, 
          plot.type = "contour", 
          n.grid = 100, 
          contour.col = "transparent", 
          xlab = "", 
          ylab = "", 
          xaxt = "n",
          yaxt = "n",
          zlim = c(0, 350),
          main = "")
  axis(side = 2, las = 2, cex.axis = .8, tck = -.01)
  axis(side = 1, las = 1, cex.axis = .8, tck = -.01)

  par(op)

  })

# Here, the density of points falling on a 100 bins plot of total NPP and NSI 
# plot is shown, from randomly generated points (10 thousenad for ecoregion at 
# 0.04 degrees of minimum distance between them) in a 0.04 degrees buffer of the
# world coastlines is shown. This is plotted alongside the linear predictor, in
# the form of a contour plot, of the HGAM model GAMM1, with 95% confidence 
# interval values of the predictor shown, indicating the lack of 95% confidence 
# of the model for areas including "low seasonality - high productivity" values.

path_coastal_points <- "~/lavori_PhD/ARMS_antartide/manuscript_june2021/new_analyses_revision_Feb2024/total/environmental_variables/scripts/qgis/additional_layers/"

points_coastal <- read.csv(paste0(path_coastal_points, "10t_random_points_sampled_nsi_npp_total.csv"), 
                           header = TRUE)

coastal_model <- ggplot(points_coastal, aes(x = npp_nsi1, y = npp_total1/1000)) + 
                 stat_bin2d(bins = 100) + 
                 scale_y_continuous(limits = c(min(tab$npp_total), 270), 
                                    n.breaks = 6) +
                 scale_x_continuous(limits = c(min(tab$npp_nsi), max(tab$npp_nsi)), 
                                    n.breaks = 8) +
                 scale_fill_gradientn(limits = c(0, 500), 
                                      #breaks = seq(10, 500, by=50), 
                                      colours = brewer.pal(n = 6, 
                                                           name = "Blues")[2:6], 
                                      na.value = "transparent") +
                 theme_classic() + 
                 theme(legend.position = "NULL", 
                       axis.text.y = element_blank()) +
                 xlab("") + 
                 ylab("") +
                 annotate("rect", 
                          xmin = min(tab$npp_nsi), 
                          xmax = .4, 
                          ymin = 75, 
                          ymax = 270,
                          alpha = .1, 
                          fill = "red")

dev.off()

# extract legend in a separate plot
coastal_model_leg <- ggplot(points_coastal, aes(x = npp_nsi1, y = npp_total1/1000)) + 
                     stat_bin2d(bins = 100) + 
                     scale_y_continuous(limits = c(min(tab$npp_total), 270), 
                                        n.breaks = 6) +
                     scale_x_continuous(limits = c(min(tab$npp_nsi), max(tab$npp_nsi)), 
                                        n.breaks = 8) +
                     scale_fill_gradientn(limits = c(0, 500), 
                                          n.breaks = 3, 
                                          colours = brewer.pal(n = 6, 
                                                               name = "Blues")[2:6], 
                                          na.value = "transparent",
                                          name = "Point\ncount",
                                          labels=c(1, 250, 500)) +
                     theme(legend.key.height = unit(.3, "cm"))

leg3 <- ggpubr::get_legend(coastal_model_leg)

layout <- c(area(1, 1, 6, 8),
            area(1, 7, 1, 8),
            area(1, 9, 6, 16),
            area(1, 15, 1, 16))

plot <- model_95confidence + leg2 + coastal_model + leg3 + plot_layout(design = layout)

# An alternative interpretation of the model's effects can be drawn from examining
# which ecoregions in the world host these kind of environments, thus providing
# an insight on what the model would look like if those areas were included. For
# this reason,  

ggsave("plots/Figure6_Seasonality and Magnitude missing rich areas.jpg", plot = plot, path = NULL, width = 10, height = 6, dpi = 600, limitsize = TRUE)

# The percentage of high energy availability points depicted in Figure 6 (c), 
# is here also reported in a table, grouped by REALM instead of ecoregion. 
stats_energy_availability <- readr::read_delim(file = paste0(path_coastal_points, "energy_availability_ecoregions.tsv"), 
                                               delim = "\t")

stats <- stats_energy_availability %>% dplyr::group_by(REALM) %>% 
                                       dplyr::summarise(n_filt = sum(NUMPOINTS, 
                                                                     na.rm = TRUE), 
                                                        n_total = sum(Total_NUMP, 
                                                                      na.rm = TRUE)) %>% 
                                       dplyr::mutate(perc_points_realm = (n_filt / n_total) * 100) %>% 
                                       dplyr::arrange(., -perc_points_realm, 
                                                      by_group = FALSE)

readr::write_delim(stats, file = "prints/Table4_energy_availability_realms_stats", 
                   delim = "\t")

# The distribution of total NPP and NSI values is here depicted from the whole
# dataset. Dashed rectangles show the extreme of NPP and NSI combination values:

# set the size for contour line around points

NSI_NPP_extremes <- ggplot(tab, aes(x = npp_nsi, 
                                    y = npp_total)) + 
                    geom_point(aes(fill = "black", 
                                   size = Observed,
                                   stroke = 1.5)) + 
                    geom_point(aes(color = ECOREGION, 
                                   size = Observed)) +
                    theme_minimal() + 
                    xlab("NSI") + 
                    ylab("Total NPP") + 
                    annotate("rect", 
                             xmin = min(tab$npp_nsi), 
                             xmax = .4, 
                             ymin = 75, 
                             ymax = 270,
                             color = "black",
                             lwd = .5,
                             alpha = 0,
                             linetype = 2) +
                    annotate(geom = "text", 
                             x = .45, 
                             y = 270, 
                             label = "S3.8",
                             color = "black") + 
                    annotate("rect", 
                             xmin = .7, 
                             xmax = max(tab$npp_nsi), 
                             ymin = 75, 
                             ymax = 270,
                             color = "black",
                             lwd = .5,
                             alpha = 0,
                             linetype = 2) +
                    annotate(geom = "text", 
                             x = .65, 
                             y = 250, 
                             label = "S3.10",
                             color = "black") + 
                    annotate("rect", 
                             xmin = min(tab$npp_nsi), 
                             xmax = .4, 
                             ymin = 0, 
                             ymax = 30,
                             color = "black",
                             lwd = .5,
                             alpha = 0,
                             linetype = 2) +
                    annotate(geom = "text", 
                             x = .45, 
                             y = 5, 
                             label = "S3.9",
                             color = "black") + 
                    annotate("rect", 
                             xmin = .7, 
                             xmax = max(tab$npp_nsi), 
                             ymin = 0, 
                             ymax = 30,
                             color = "black",
                             lwd = .5,
                             alpha = 0,
                             linetype = 2) +
                    annotate(geom = "text", 
                             x = .65, 
                             y = 5, 
                             label = "S3.11",
                             color = "black") 

ggsave("plots/S3.7_Extreme values of total NPP and NSI.jpg", 
       plot = NSI_NPP_extremes, 
       path = NULL, 
       width = 9, 
       height = 9, 
       dpi = 600, 
       limitsize = TRUE)

