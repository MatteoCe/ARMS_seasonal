#!/usr/bin/Rscript

library(phyloseq)
library(tidyverse)
library(lattice)
library(cowplot)
library(data.table)
library(mgcv)
library(fitdistrplus)
library(corrplot)
library(olsrr)
library(ncf)
library(Metrics)
library(gratia)

# DATASET PREPARATION

# load the otu table and sample data

count_tab <- read.table("clust_swarm_mod.counts.tsv", 
						header=T, 
						row.names=1, 
						check.names=F)

sample_info <- read.table("sample_info_env.tsv",
						  header=T,
						  row.names=1,
						  check.names=F,
						  sep ='\t',
						  stringsAsFactors=T)

# as some samples may have been removed during the bioinformatic analyses, the sample info tab will be filetered to
# remove them and match the count tab

sample_info_tab <- sample_info[rownames(sample_info) %in% colnames(count_tab),]

# format the sample data as a phyloseq object

sampledata <- sample_data(sample_info_tab)

# format the otu table as a phyloseq object

OTU <- otu_table(count_tab, taxa_are_rows = TRUE)

# merge the two objects into a single phyloseq object

physeq <- phyloseq(OTU, sampledata)

# calculate the number of OTUs from phyloseq object as vectors

rich <- estimate_richness(physeq, split = TRUE, measures = "Observed")

# store number of reads

total_reads <- sample_sums(physeq)

# create a dataframe with OTU richness, all the environmental variables and other useful factors

# set vector of most important field's names

fields <- c("id_site",
			"protection",
			"protection_level",
			"months",
			"depth",
			"lat",
			"lon",
			"chl_mean",
			"da_mean",
			"nsi",
			"par_mean",
			"prim_mean",
			"sal_mean",
			"sal_range",
			"temp_mean",
			"temp_range",
			"ECOREGION",
			"Lat_Zone")

tab <- cbind(as.data.frame(rownames(rich)),
			 sample_data(physeq)[,fields],
			 total_reads,
			 rich$Observed)

colnames(tab) <- c("samples",
						fields,
						"total_reads",
						"Observed")

### DATASET EXPLORATION

# check the presence of outliers for the response and explanatory variables
png("dotchart.png", units="cm", width=30, height=30, res=300)

op<- par(mfrow=c(5,3),mar=c(3,3,3,1))
dotchart(tab$Observed,main="#OTUs")
dotchart(tab$nsi,main="NSI - Normalized Seasonality Index")
dotchart(tab$chl_mean,main="CHL - Chlorophyll mean")
dotchart(tab$da_mean,main="DA - Diffuse Attenuation coefficient mean")
dotchart(tab$par_mean,main="PAR - Photosynthetically Active Radiation mean")
dotchart(tab$prim_mean,main="NPP - Net Primary Production mean")
dotchart(tab$sal_mean,main="SSS - Sea Surface Salinity mean")
dotchart(tab$sal_range,main="SSS - Sea Surface Salinity range")
dotchart(tab$temp_mean,main="SST - Sea Surface Temperature mean")
dotchart(tab$temp_range,main="SST - Sea Surface Temperature range")
dotchart(tab$months,main="Months of deployment")
dotchart(tab$depth,main="Depth")
dotchart(abs(tab$lat),main="Latitude")
dotchart(tab$total_reads,main="Number of reads")
par(op)

dev.off()

# the dotcharts show the presence of observations with extreme values of both SSS mean and range, which
# derive from samples that were located in areas (the Black and Baltic seas) naturally showing these values.
# Something similar can be said for both chlorophyll mean and primary productivity mean, but the observed
# values are not "extremes", rather higher then the mean and less "sampled" conditions.

# Salinity transformation:

# in order to test effectively the influence of "extreme" salinity values for those two variables, including
# conditions of lower salinity mean and higher salinity range, these two were log-transformed. However, as
# salinity mean shows extremes with values lower than the mean, simply applying the logarithm
# is not a valuable option. The distribution was first reflected, the log-transformation applied and the reflected again,
# due to left-skewed nature of the data.

# reflecting the distribution (to all values, the maximum was subtracted and then turned to absolute positive
# values). A constant of 1 is applied.
sal_ref <- abs(tab$sal_mean - (max(tab$sal_mean) + 1))

# apply the logarithm
sal_ref_log <- log(sal_ref)

# reflect the values again and include it in the correlation table
sal_ref_log_ref <- abs(sal_ref_log -(max(sal_ref_log) +1))
tab$l.ref.sal_mean <- sal_ref_log_ref

# apply the logarithm to the salinity range
tab$l.sal_range <- log(tab$sal_range)

# The effect of the transformation for salinity mean is depicted in this dotchart plot
png("SSS mean transformation process.png", units="cm", width=20, height=30, res=300)

op<- par(mfrow=c(5,1),mar=c(3,1,2,1))
layout(matrix(c(1,2,3,4,5),ncol=1),heights=c(1.5,3,3,3,3))
plot.new()
text(0.5,0.5,"SSS mean transformation process",cex=1.3,font=2)
dotchart(tab$sal_mean,main="Original distribution")
dotchart(sal_ref,main="Reflected and turned to absolute")
dotchart(sal_ref_log,main="Log-transformation")
dotchart(sal_ref_log_ref,main="Reflected and turned to absolute")
par(op)

dev.off()

# the next plot will compare the relationship between the response variable (OTU richness) and the SSS - related variables
# before and after the transformation processes.
png("SSS transformation vs number of OTUs.png", units="cm", width=30, height=30, res=300)

op<- par(mfrow=c(3,2),mar=c(3,3,1,3))
layout(matrix(c(1,2,3,1,4,5),ncol=2),heights=c(1,3,3))
plot.new()
text(0.5,0.5,"Effect of transformation on the relationship SSS vs #OTUs",cex=1.5,font=2)
plot(tab$sal_mean,tab$Observed,main="SSS mean")
plot(tab$l.ref.sal_mean,tab$Observed,main="SSS mean transformed")
plot(tab$sal_range,tab$Observed,main="SSS range")
plot(tab$l.sal_range,tab$Observed,main="SSS range transformed")
par(op)

dev.off()

# check the distribution of the response variable

# here two different methods are implemented: in the first case, the AIC is used to check which distribution has the
# lower akaike value if fit to the real data, using the fitdistrplus package

distributions <- c("nbinom", "norm", "pois", "gamma", "geom")

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



# a vector with all the possible values
# of the response variables (from minimum values to maximum) is created and the density of different
# common distributions are computed using the main characteristics of those distributions, but calculated
# on the real observed dataset.

# A distribution, with the characteristics of the real dataset, but with all possible
# observations is created and can be compared to the histogram of the values in the real dataset
resp_var <- tab$Observed
name <- "# OTUs"

# before continuing the break parameter of the hist command can be checked and chosen on behalf of
# better representativity of the real distribution of the data. Here 20 is set, but 5 can work as well
br <- 10

# set minimum and maximum values of the response variable chosen
min_set <- min(resp_var)
max_set <- max(resp_var)

# define a theoric distribution with all values from minimum and maximum of those observed 
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

# calculate density for a Geometric distribution (k value of Geometric distribution =1 of Negative binomial and variance as quadratic function of the mean)
y_geom <- dnbinom(x, size = 1, mu = mean_tab)

# plot the response variable with all possible theoretic distributions
png("Distribution Families.png", units="cm", width=30, height=30, res=300)

op <- par(mfrow=c(3,2))

# plot the histogram of the real response variable dataset
hist(resp_var, breaks=br, main=name)

plot(x, y_nbinom, main="Density Plot of Negative Binomial Distribution",
	 type = "h",
	 lwd = 2,
	 xlab = "x",
	 ylab = "Density")
text(400,.004,paste0("AIC=",round(results[1,2],3)),cex=1.5)

plot(x, y_norm, main="Density Plot of Normal Distribution",
	 type = "h",
	 lwd = 2,
	 xlab = "x",
	 ylab = "Density")
text(400,.0025,paste0("AIC=",round(results[2,2],3)),cex=1.5)

plot(x, y_pois, main="Density Plot of Poisson Distribution",
	 type = "h",
	 lwd = 2,
	 xlab = "x",
	 ylab = "Density")
text(400,.025,paste0("AIC=",round(results[3,2],3)),cex=1.5)

plot(x, y_gamma, main="Density Plot of Gamma Distribution",
	 type = "h",
	 lwd = 2,
	 xlab = "x",
	 ylab = "Density")
text(400,.000000125,paste0("AIC=",round(results[4,2],3)),cex=1.5)

plot(x, y_geom, main="Density Plot of Geometric Distribution",
	 type = "h",
	 lwd = 2,
	 xlab = "x",
	 ylab = "Density")
text(400,.0045,paste0("AIC=",round(results[5,2],3)),cex=1.5)

par(op)

dev.off()

# plot each environmental variable against the #OTUs and calculate the Pearson's rho
png("Response vs explanatory Pearson's rho.png", units="cm", width=30, height=30, res=300)

op <- par(mfrow=c(3,3))

pear_corr <- round(cor(tab$nsi,tab$Observed,method="pearson"),2)
plot(tab$nsi,tab$Observed, main=pear_corr,xlab="NSI",ylab="#OTUs")
lines(lowess(tab$nsi, tab$Observed), col = "red")

pear_corr <- round(cor(tab$da_mean,tab$Observed,method="pearson"),2)
plot(tab$da_mean,tab$Observed, main=pear_corr,xlab="Diffuse Attenuation Coef.",ylab="#OTUs")
lines(lowess(tab$da_mean, tab$Observed), col = "red")

pear_corr <- round(cor(tab$par_mean,tab$Observed,method="pearson"),2)
plot(tab$par_mean,tab$Observed, main=pear_corr,xlab="Photosynthetically Active Radiation",ylab="#OTUs")
lines(lowess(tab$par_mean, tab$Observed), col = "red")

pear_corr <- round(cor(tab$chl_mean,tab$Observed,method="pearson"),2)
plot(tab$chl_mean,tab$Observed, main=pear_corr,xlab="Chlorophyll mean",ylab="#OTUs")
lines(lowess(tab$chl_mean, tab$Observed), col = "red")

pear_corr <- round(cor(tab$prim_mean,tab$Observed,method="pearson"),2)
plot(tab$prim_mean,tab$Observed, main=pear_corr,xlab="Net Primary Productivity mean",ylab="#OTUs")
lines(lowess(tab$prim_mean, tab$Observed), col = "red")

pear_corr <- round(cor(tab$l.ref.sal_mean,tab$Observed,method="pearson"),2)
plot(tab$l.ref.sal_mean,tab$Observed, main=pear_corr,xlab="Sea Surface Salinity mean",ylab="#OTUs")
lines(lowess(tab$l.ref.sal_mean, tab$Observed), col = "red")

pear_corr <- round(cor(tab$temp_mean,tab$Observed,method="pearson"),2)
plot(tab$temp_mean,tab$Observed, main=pear_corr,xlab="Sea Surface Temperature mean",ylab="#OTUs")
lines(lowess(tab$temp_mean, tab$Observed), col = "red")

pear_corr <- round(cor(tab$l.sal_range,tab$Observed,method="pearson"),2)
plot(tab$l.sal_range,tab$Observed, main=pear_corr,xlab="Sea Surface Salinity range",ylab="#OTUs")
lines(lowess(tab$l.sal_range, tab$Observed), col = "red")

pear_corr <- round(cor(tab$temp_range,tab$Observed,method="pearson"),2)
plot(tab$temp_range,tab$Observed, main=pear_corr,xlab="Sea Surface Temperature range",ylab="#OTUs")
lines(lowess(tab$temp_range, tab$Observed), col = "red")

par(op)

dev.off()

# the effect of management protection on species richness and diversity is assessed using data from www.protectedplanet.net
# this is also assessed by looking at the variability 
png("Protection effect on diversity global.png", units="cm", width=25, height=30, res=300)

op <- par(mfrow=c(2,1),mar=c(5,5,2,1))

boxplot(Observed ~ protection, data=tab,
		xlab="",
		ylab="Number of OTUs",
		main="Global effect of protection status on diversity metrics",
		xaxt="n")
axis(1, at=c(1,2),
	 cex.axis=1.3,
	 labels=c("Protection","No protection"))

boxplot(Observed ~ protection_level, data=tab,
		xlab="Protection status - IUCN Category",
		ylab="Number of OTUs",
		main="",
		cex.axis=1.3,
		cex.lab=1.5)

par(op)

dev.off()

# the same is observed at lower spatial scale, based on the REALM regionalization of Spalding et al. 2007
png("Protection effect on richness at the Ecoregion scale.png", units="cm", width=30, height=30, res=300)

bwplot(Observed ~ protection | ECOREGION, data=tab)

dev.off()

# no variability issues are detected.

# on the other hand, a high variability can be detected between samples collected at the same "site", factoring group
# created using a buffer resembling the spatial resolution of the environmental layers used in the analyses.
# From the following boxplots, it appears evident that the variance in diversity or richness is greater for some sites
# and lower in others. This variability does not appear to reflect the latitude per se, rather it appears greater
# with increasing values of OTU richness, possibly due to the type of distribution of this data (negative binomial
# with variance increasing with mean), or simply to the methodology, such as DNA metabarcoding, with increasing
# "uncertainty" analyzing samples with a higher number of species.

# here, boxplots are produced on a version of the dataset that has been ordered by Latitudinal Zone, showing a trend of
# increasing variability per site with an higher mean of the number of OTUs, rather than the latitude per se.

tab_ord <- tab[order(tab$Lat_Zone),]
tab_ord$id_site <- factor(tab_ord$id_site, levels=unique(tab_ord$id_site))

ggplot(tab_ord) +
	geom_boxplot(aes(id_site,Observed)) +
	theme_minimal() +
    theme(axis.text.x = element_text(angle = 45,
									 hjust = 1,
									 vjust = 1,
									 size= 11,
									 face="bold")) +
	xlab("Site/year") + ylab("# OTUs")

### Richness

# Boxplots of number of OTUs index grouped by Ecoregion Spalding et al. 2007

tab2 <- tab %>% group_by(ECOREGION) %>% mutate(Obs_median = median(Observed))
tab2_ord <- tab2[order(tab2$Obs_median),]
tab2_ord$ECOREGION <- factor(tab2_ord$ECOREGION,levels=unique(tab2_ord$ECOREGION))

pleg <- ggplot(tab2_ord) +
		geom_boxplot(aes(ECOREGION,Observed,fill=Lat_Zone)) +
		theme(legend.text=element_text(size=11,
									   face="bold"),
			  legend.title=element_text(size=12,
										face="plain")) +
		scale_fill_manual(values=c("#98D3FF","#83EC87","#FF9898"),
						  name=("Latitudinal \n    zone"))

leg <- ggpubr::get_legend(pleg)

# Simple plot with only OTU number

p1_simp <- ggplot(tab2_ord) +
		geom_boxplot(aes(ECOREGION,Observed,fill=Lat_Zone)) +
		theme_half_open() +
		ggtitle("# OTUs") +
		theme(axis.text.x = element_text(angle = 45,
										 hjust = 1,
										 vjust = 1,
										 size= 11,
										 face="bold"),
			  plot.title = element_text(hjust=0.5,
										face="bold"),
			  panel.grid.major.x = element_blank(),
			  legend.text=element_text(size=11,
									   face="bold"),
			  legend.title=element_text(size=12,
										face="plain")) +
		xlab(NULL) +
		ylab(NULL) +
		scale_y_continuous(position = "right") +
		scale_fill_manual(values=c("#98D3FF","#83EC87","#FF9898"),
						  name=("Latitudinal \n    zone"))

p1_simp_fin <- cowplot::plot_grid(plot.new(), p1_simp, plot.new(),
								align = "vh",
								axis="t",
								ncol = 3,
								nrow = 1,
								rel_heights = c(1, 1, 1),
								rel_widths = c(0.16, 1, 0.01))

ggsave("Number of OTUs per Ecoregion.png", plot = p1_simp_fin, device = png(), path = NULL,scale = 1, width = 7.5, height = 9, dpi = 300, limitsize = TRUE)
dev.off()

# Schematic plot of seasonality dynamics for map

# draw dataframe
t=seq(0,6,0.1)
y=sin(t)
y1=rep(0.5,61)
y2=y+(abs(min(y))+0.2)
y3 <- dpois(seq(0:60), lambda = 7)
seas <- data.frame(t,y1,y2,y3)

p1 <- ggplot(seas,aes(t,y1)) + theme_classic() + theme(axis.title=element_text(size=13,face="bold"), axis.text.y=element_blank(), axis.ticks.y = element_blank()) + ylab("Net Primary Production") + xlab("") + geom_area(fill="#00802b",color="black") + scale_y_continuous(limits = c(0,2.5)) + scale_x_continuous(breaks=seq(0,6,.5),labels=c("J","F","M","A","M","J","J","A","S","O","N","D","J"))
p2 <- ggplot(seas,aes(t,y2)) + theme_classic() + theme(axis.title=element_text(size=13,face="bold"), axis.text.y=element_blank(), axis.ticks.y = element_blank()) + ylab("") + xlab("Month of year") + geom_area(fill="#00802b",color="black") + scale_y_continuous(limits = c(0,3.5)) + scale_x_continuous(breaks=seq(0,6,.5),labels=c("J","F","M","A","M","J","J","A","S","O","N","D","J"))
p3 <- ggplot(seas,aes(t,y3)) + theme_classic() + theme(axis.title=element_text(size=13,face="bold"), axis.text.y=element_blank(), axis.ticks.y = element_blank()) + ylab("") + xlab("") + geom_area(fill="#00802b",color="black") + scale_y_continuous(limits = c(0,0.15)) + scale_x_continuous(breaks=seq(0,6,.5),labels=c("J","F","M","A","M","J","J","A","S","O","N","D","J"))

seasonality <- cowplot::plot_grid(p1,p2,p3, align = "vh", axis="t", ncol = 3, nrow = 1, rel_heights = c(1, 1, 1), rel_widths = c(1, 1, 1))

ggsave("Schematic Seasonality.png", plot = seasonality, device = png(), path = NULL,scale = 3, width = 3, height = 1, dpi = 200, limitsize = TRUE)
dev.off()

# Tests


# Prepare a subset of the dataset with only response and explanatory variables.
# Two subset tables are created, one with the original names and another with labels for plots and such

tab_corr <- tab[,c("Observed",
				   "nsi",
				   "da_mean",
				   "par_mean",
				   "chl_mean",
				   "prim_mean",
				   "l.ref.sal_mean",
				   "temp_mean",
				   "l.sal_range",
				   "temp_range",
				   "depth",
				   "months",
				   "lat")]

tab_corr_names <- tab_corr

colnames(tab_corr_names) <- c("#OTUs",
							"NSI",
							"DA mean",
							"PAR mean",
							"CHL mean",
							"NPP mean",
							"SSS mean",
							"SST mean",
							"SSS range",
							"SST range",
							"Depth",
							"Months",
							"Latitude")

# To prepare the base model, the presence of spatial auto-correlation can be evaluated using the ncf package.
# The presence of spatial auto-correlation will be inspected with a plot, that will show what type of correlation
# (or not) is present at differing distances between the samples.

corr_spaz <- spline.correlog(x=tab[,"lon"],y=tab[,"lat"],z=tab[,"Observed"],latlon=TRUE)

png("Spatial autocorrelation Dataset.png", units="cm", width=30, height=20, res=300)
plot(corr_spaz)

dev.off()

# The wiggly form of the plot at higher distances reflects the sparsely distributed dataset we have, with most
# of the samples from the mediterranean, red and northern seas, and sparse samples in the pacific and antarctica.
# This shows more negative and positive correlation without a clear sense, but the most important part is the
# high correlation at reduced distances. This is mostly due to the fact that for each X and Y (latitude and
# longitude combinations) we have more than one sample, often in three replicates. This can be inspected by
# looking at the factor level "site", based on the lat lon and the number of ARMS deployed in that location in
# a particular year.

# In order to account for this auto-correlation, a random effect using that factor can be used.
# The effectiveness of the random effect can be assessed by checking the AIC of the model without it with the one that
# has it, and then comparing the autocorrelation left in the pearson residuals of the two models.
# A tensor product, representing our main hypothesis, will be used as a fixed effect.

GAMM_norand <- gam(Observed ~ te(nsi, prim_mean),method="REML",family=nb(theta=NULL,link="log"),data=tab)
GAMM_rand <- gam(Observed ~ te(nsi, prim_mean) + s(id_site, bs="re"),method="REML",family=nb(theta=NULL,link="log"),data=tab)

corr_spaz_norand <- spline.correlog(x=tab[,"lon"],y=tab[,"lat"],z=resid(GAMM_norand,type="pearson"),latlon=TRUE)
corr_spaz_rand <- spline.correlog(x=tab[,"lon"],y=tab[,"lat"],z=resid(GAMM_rand,type="pearson"),latlon=TRUE)

sink("AIC_randomEffect.txt")
AIC(GAMM_norand,GAMM_rand)
sink()

png("Spatial Autocorrelation Effect with Random model.png", units="cm", width=30, height=30, res=300)

op <- par(mfrow=c(2,2))

plot(corr_spaz_norand,main="Spatial autocorrelation without Random Effect")

lin.pred_norand <- napredict(GAMM_norand$na.action, GAMM_norand$linear.predictors)
plot(lin.pred_norand, resid(GAMM_norand), ylab=("residuals"), xlab=("Linear predictor"), main="Residuals vs. Linear predictor \n without Random Effect")
lines(lowess(lin.pred_norand, resid(GAMM_norand)), col = "red")

plot(corr_spaz_rand,main="Spatial autocorrelation with Random Effect")

lin.pred_rand <- napredict(GAMM_rand$na.action, GAMM_rand$linear.predictors)
plot(lin.pred_rand, resid(GAMM_rand), ylab=("residuals"), xlab=("Linear predictor"), main="Residuals vs. Linear predictor \n with Random Effect")
lines(lowess(lin.pred_rand, resid(GAMM_rand)), col = "red")

par(op)

dev.off()



# Thus, the base model will include a tensor product "an interaction of each pair of basis functions for each
# marginal term, and a penalty for each marginal term that penalizes the average wiggliness in that term"
# between primary productivity mean and the normalized seasonality index, and a random effect smoother for
# the site factor.


# Hypothesis testing part 1:
# FULL MODEL

# A Generalized Linear Model assuming a negative binomial distribution of the response variable is fitted using the
# environmental variables that contribute the least to the collinearity. The VIFs were thus less then 5, a bit
# higher than the VIF values reported as acceptable by many papers, but lower then ohers
# (see Zuur https://besjournals.onlinelibrary.wiley.com/doi/pdf/10.1111/j.2041-210X.2009.00001.x).

# Considering the established non-normality of our response variable, the first approach is to apply
# a generalized linear model that assume a negative-binomial distribution of the response variable, typical
# of count data. Using the package MASS it is possible to fit it, and in this case the link between the response
# variable and the systematic part, a function of the explanatory variables, is the logarithm.

# In order to test whether the seasonality is a good environmental descriptor of the observed richness (#OTUs)
# different GLM models will be fitted, one for each variable, and the estimates and p values will be compared
# between them. This will be done using both the summary and drop1 output, the last one for a Chisquare test
# of log-likelihood.

# use the subset tab as made before the glms and gams

tab_corr_mat <- cor(tab_corr_names)

png("Correlation plot Explanatory Vars.png", units="cm", width=30, height=30, res=300)
corrplot(tab_corr_mat, method="number", type="lower", tl.col="black")

dev.off()

# testing variance inflating factors, each variable showing the highest score of VIFs is removed and the model
# refitted without it

# The process is performed untill all variables have a VIF < 3

M1 <- lm(Observed~nsi+prim_mean+da_mean+chl_mean+par_mean+l.ref.sal_mean+temp_mean+l.sal_range+temp_range+months+depth,data=tab)

ols_coll_diag(M1)

M2 <- lm(Observed~nsi+prim_mean+da_mean+chl_mean+par_mean+l.ref.sal_mean+l.sal_range+temp_range+months+depth,data=tab)

ols_coll_diag(M2)

M3 <- lm(Observed~nsi+prim_mean+da_mean+par_mean+l.ref.sal_mean+l.sal_range+temp_range+months+depth,data=tab)

ols_coll_diag(M3)

M4 <- lm(Observed~nsi+prim_mean+par_mean+l.ref.sal_mean+l.sal_range+temp_range+months+depth,data=tab)

ols_coll_diag(M4)

M5 <- lm(Observed~nsi+prim_mean+par_mean+l.sal_range+temp_range+months+depth,data=tab)

ols_coll_diag(M5)

M6 <- lm(Observed~nsi+prim_mean+l.sal_range+temp_range+months+depth,data=tab)

sink("VIFs_final.txt")
ols_coll_diag(M6)
sink()

# Plot each environmental variable of the one chosen by the VIF selection process
# against the #OTUs and calculate the Pearson's rho

png("Response vs uncollinear explanatory Pearson's rho.png", units="cm", width=30, height=30, res=300)
op <- par(mfrow=c(3,2))

pear_corr <- round(cor(tab$nsi,tab$Observed,method="pearson"),2)
plot(tab$nsi,tab$Observed, main=pear_corr,xlab="NSI",ylab="#OTUs")
lines(lowess(tab$nsi, tab$Observed), col = "red")

pear_corr <- round(cor(tab$prim_mean,tab$Observed,method="pearson"),2)
plot(tab$prim_mean,tab$Observed, main=pear_corr,xlab="Net Primary Productivity mean",ylab="#OTUs")
lines(lowess(tab$prim_mean, tab$Observed), col = "red")

pear_corr <- round(cor(tab$l.sal_range,tab$Observed,method="pearson"),2)
plot(tab$l.sal_range,tab$Observed, main=pear_corr,xlab="Sea Surface Salinity range",ylab="#OTUs")
lines(lowess(tab$l.sal_range, tab$Observed), col = "red")

pear_corr <- round(cor(tab$temp_range,tab$Observed,method="pearson"),2)
plot(tab$temp_range,tab$Observed, main=pear_corr,xlab="Sea Surface Temperature range",ylab="#OTUs")
lines(lowess(tab$temp_range, tab$Observed), col = "red")

pear_corr <- round(cor(tab$months,tab$Observed,method="pearson"),2)
plot(tab$months,tab$Observed, main=pear_corr,xlab="Months of deployment",ylab="#OTUs")
lines(lowess(tab$months, tab$Observed), col = "red")

pear_corr <- round(cor(tab$depth,tab$Observed,method="pearson"),2)
plot(tab$depth,tab$Observed, main=pear_corr,xlab="Depth",ylab="#OTUs")
lines(lowess(tab$depth, tab$Observed), col = "red")

par(op)

dev.off()

# MODEL selection

# A forward selection is performed testing wether other explanatory variables can be included apart from the
# main hypothesis which refers to the main effect being the availability of NPP throughout the year.
# Considering that, as the primary productivity can become more avaiable with less seasonality, the NSI index
# will be included in an interaction term with the mean Primary Productivity, testing whether with more
# availability through time, the amount of NPP might have an effect.

# Based on the information gathered from the above GAM models with single explanatory variables, from the least
# collinear variables the one explaining the more variance is chosen to add the the base model, and the
# fitness of this model will be evaluated by comparic the AIC with the base model, together with the info
# from the gam.check

# This will be computed using the gam command of the mgcv package following Pedersen et al 2019.

GAMM0 <- gam(Observed ~ te(nsi, prim_mean) + s(months) + s(depth) + s(temp_range) + s(l.sal_range) + s(id_site, bs="re"),method="REML",family=nb(theta=NULL,link="log"),data=tab)

sink("Summary statistics First HGAM model.txt")
summary(GAMM0)
sink()

GAMM1 <- gam(Observed ~ te(nsi, prim_mean) + s(temp_range) + s(l.sal_range) + s(id_site, bs="re"),method="REML",family=nb(theta=NULL,link="log"),data=tab)

sink("Summary statistics Second HGAM model.txt")
summary(GAMM1)
sink()

png("NSI-NPP mean Physical extremes.png", units="cm", width=20, height=20, res=300)
gam.check(GAMM1)
dev.off()

png("NSI-NPP mean Physical extremes_smoothers.png", units="cm", width=20, height=20, res=300)
draw(GAMM1)
dev.off()

GAMM2 <- gam(Observed ~ te(nsi, prim_mean) + s(id_site, bs="re"),method="REML",family=nb(theta=NULL,link="log"),data=tab)

models_list <- list()

models_list[[1]] <- GAMM0
models_list[[2]] <- GAMM1
models_list[[3]] <- GAMM2

coefs <- list()

for (i in seq(1,3)) {
	
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

rownames(param_full_models) <- c("GAMM0","GAMM1","GAMM2")

param_full_models$FixEff <- c("NSI*NPP mean + All non-collinear","NSI*NPP mean + Physical extremes","NSI*NPP mean")

sink("Final_models_AIC.txt")
param_full_models
sink()

# plot smoothers

png("NSI_NPP models effects.png", units="cm", width=30, height=10, res=300)

op <- par(mfrow=c(1,3),mar=c(5,5,3,1))
vis.gam(GAMM0,type="response",plot.type="contour",main="NSI*NPP mean + All non-collinear",xlab="",ylab="NPP mean")
vis.gam(GAMM1,type="response",plot.type="contour",main="NSI*NPP mean + Physical extremes",xlab="NSI",ylab="",yaxt="n")
vis.gam(GAMM2,type="response",plot.type="contour",main="NSI*NPP mean",xlab="",ylab="",yaxt="n")
par(op)

dev.off()


# Testing part 2:

# Similarly, GAMs are fitted for each var. The results are again a table and a plot. The first includes
# the estimated degrees of freedom of the smoother for each var, the AIC, the % variance explained and the
# estimated theta parameter. The table is order by ascending values of estimated degrees of freedom. The plot
# includes the smoothers of each var with estimated degrees of freedom and % variance explained.

coefs <- list()
models_gam_list <- list()

mod_res <- as.data.frame(cbind(c("nsi","chl_mean","prim_mean","temp_mean"),c("NSI","CHL mean","NPP mean","SST mean")))
colnames(mod_res) <- c("vars","names")

png("Alternative GAMs smoothers.png", units="cm", width=20, height=20, res=300)
op <- par(mfrow=c(2,2))

for (i in seq(1,nrow(mod_res))) {
	
	tab_name <- mod_res[i,]$vars
	off_name <- mod_res[i,]$names
	
	formula <- substitute(Observed ~ s(x) + s(id_site, bs="re"), list(x=as.name(tab_name)))
	
	## This is just a temporary model to estimate the theta value using the function nb (for negative binomial distributions)
	## The theta will be extracted and the model refitted using the negbin function, similarly to the
	## function glm.nb used earlier for the GLMs
	#temp_model <- gam(formula,family=nb(theta=NULL,link="log"),data=tab)
	## extract the family component of the model
	#f <- family(temp_model)
	## get the theta value
	#k <- f$getTheta(trans=T)
	## refit the model using the negbin function and specifying the theta
	model <- gam(formula,family=nb(theta=NULL,link="log"),method="REML",data=tab)
	
	models_gam_list[[paste0("GAM",i)]] <- model
	
	sum <- summary(model)
	
	# extract estimated degrees of freedom and p value for the fixed effect
	edf <- round(sum$edf,2)[1]
	pv <- round(sum$s.table[1,4],3)
	
	if (pv == 0) {pv <- "<2e-16"
		} else {
			
		}
	
	par_names <- c("edfs","p-value")
	pars_gam <- cbind(edf,pv)
	colnames(pars_gam) <- par_names
	
	coefs[[i]] <- pars_gam
	
	plot(model, select=1, main=paste("df:",edf,"\n","p-value",pv),xlab=off_name,ylab="#OTUs")
	
	}

par(op)
dev.off()


# to compare accuracy of models following Simon et al 2022 on Random Forest interpretation
# REMEMBER WHY I SET K=8 IN TEMP RANGE e SSS RANGE FOR THE FULL MODEL BECAUSE THE TRAINING DATASET DOESN'T HAVE ENOUGH OBSERVATIONS 

# Cross-Validation on different models with SST mean, NPP mean, CHL mean, NSI and the GAMM3 and GAMM1 from the previous
# model selection steps. In the GAMM1 (here GAMMd), the smoothers were 

ntrial <- 100

# create the dataframe storing the RMSE values

res <- data.frame(matrix(0,ncol=ntrial,nrow=6),row.names=c("SST mean","NSI","Tensor product \n NSI and NPP mean","Full Model","CHL mean","NPP mean"))

colnames(res) <- c(1:ntrial)

moda <- list()
modb <- list()
modc <- list()
modd <- list()
mode <- list()
modf <- list()

modc_fix <- list()
modd_fix <- list()

modela <- list()
modelb <- list()
modelc <- list()
modeld <- list()
modele <- list()
modelf <- list()

for (i in 1:ntrial) {

# create training and validation datasets:

data_set_size <- floor(nrow(tab)/1.5)
# Generate a random sample of "data_set_size" indexes
indexes <- sample(1:nrow(tab), size = data_set_size)
# Assign the data to the correct sets
training <- tab[indexes,]
validation1 <- tab[-indexes,]

# run the different models
GAMMa <- gam(Observed ~ s(temp_mean) + s(id_site, bs="re"),method="REML",family=nb(theta=NULL,link="log"),data=training)
GAMMb <- gam(Observed ~ s(nsi) + s(id_site, bs="re"),method="REML",family=nb(theta=NULL,link="log"),data=training)
GAMMc <- gam(Observed ~ te(nsi,prim_mean) + s(id_site, bs="re"),method="REML",family=nb(theta=NULL,link="log"),data=training)
GAMMd <- gam(Observed ~ te(nsi,prim_mean) + s(temp_range, k=8) + s(l.sal_range, k=8) + s(id_site, bs="re"),method="REML",family=nb(theta=NULL,link="log"),data=training)
GAMMe <- gam(Observed ~ s(chl_mean) + s(id_site, bs="re"),method="REML",family=nb(theta=NULL,link="log"),data=training)
GAMMf <- gam(Observed ~ s(prim_mean) + s(id_site, bs="re"),method="REML",family=nb(theta=NULL,link="log"),data=training)

# store the models in different lists
modela[[i]] <- GAMMa
modelb[[i]] <- GAMMb
modelc[[i]] <- GAMMc
modeld[[i]] <- GAMMd
modele[[i]] <- GAMMe
modelf[[i]] <- GAMMf

# the prediction is performed excluding the random effect on the prediction
catPreda=c(predict.gam(GAMMa,newdata=validation1[,"temp_mean",drop=FALSE],exclude='s(id_site)',newdata.guaranteed = TRUE,se.fit = TRUE,type=c("response")))
catPredb=c(predict.gam(GAMMb,newdata=validation1[,"nsi",drop=FALSE],exclude='s(id_site)',newdata.guaranteed = TRUE,se.fit = TRUE,type=c("response")))
catPredc=c(predict.gam(GAMMc,newdata=validation1[,c("nsi","prim_mean")],exclude='s(id_site)',newdata.guaranteed = TRUE,se.fit = TRUE,type=c("response")))
catPredd=c(predict.gam(GAMMd,newdata=validation1[,c("nsi","prim_mean","temp_range","l.sal_range")],exclude='s(id_site)',newdata.guaranteed = TRUE,se.fit = TRUE,type=c("response")))
catPrede=c(predict.gam(GAMMe,newdata=validation1[,"chl_mean",drop=FALSE],exclude='s(id_site)',newdata.guaranteed = TRUE,se.fit = TRUE,type=c("response")))
catPredf=c(predict.gam(GAMMf,newdata=validation1[,"prim_mean",drop=FALSE],exclude='s(id_site)',newdata.guaranteed = TRUE,se.fit = TRUE,type=c("response")))

# calculate the RMSE "Root-mean-square deviation"
rmsea <- rmse(validation1$Observed,catPreda$fit)
rmseb <- rmse(validation1$Observed,catPredb$fit)
rmsec <- rmse(validation1$Observed,catPredc$fit)
rmsed <- rmse(validation1$Observed,catPredd$fit)
rmsee <- rmse(validation1$Observed,catPrede$fit)
rmsef <- rmse(validation1$Observed,catPredf$fit)

# store them as a dataframe for comparison
res[,i] <- data.frame(c(rmsea,rmseb,rmsec,rmsed,rmsee,rmsef))

resPata <- data.frame(catPreda$fit,catPreda$se.fit,validation1[,"temp_mean",drop=FALSE])
resPatb <- data.frame(catPredb$fit,catPredb$se.fit,validation1[,"nsi",drop=FALSE])
resPatc <- data.frame(catPredc$fit,catPredc$se.fit,validation1[,c("nsi","prim_mean")])
resPatd <- data.frame(catPredd$fit,catPredd$se.fit,validation1[,c("nsi","prim_mean","temp_range","l.sal_range")])
resPate <- data.frame(catPrede$fit,catPrede$se.fit,validation1[,"chl_mean",drop=FALSE])
resPatf <- data.frame(catPredf$fit,catPredf$se.fit,validation1[,"prim_mean",drop=FALSE])

# aggregate all predictions in different lists
moda[[i]] <- resPata
modb[[i]] <- resPatb
modc[[i]] <- resPatc
modd[[i]] <- resPatd
mode[[i]] <- resPate
modf[[i]] <- resPatf

# fixed NSI smoother for NPP mean

for (j in seq(0.4,0.8,by=0.1)) {

# create "copy" of the validation dataset
validation2 <- validation1

# fix the values for the seasonality index of the validation dataset to the values chosen by the loop arguments
# (from 0.4 to 0.8). These "moderate" values were chosen to restrain the prediction of the model to less extreme
# values, not necessarily more represented.
validation2$nsi <- j

catPredc_nsifix=c(predict.gam(GAMMc,newdata=validation2[,c("nsi","prim_mean")],exclude='s(id_site)',newdata.guaranteed = TRUE,se.fit = TRUE,type=c("response")))
catPredd_nsifix=c(predict.gam(GAMMd,newdata=validation2[,c("nsi","prim_mean","temp_range","l.sal_range")],exclude=c('s(id_site)','s(l.sal_range)','s(temp_range)'),newdata.guaranteed = TRUE,se.fit = TRUE,type=c("response")))

resPatc_fix <- data.frame(catPredc_nsifix$fit,catPredc_nsifix$se.fit,validation2[,c("nsi","prim_mean")])
resPatd_fix <- data.frame(catPredd_nsifix$fit,catPredd_nsifix$se.fit,validation2[,c("nsi","prim_mean","temp_range","l.sal_range")])
modc_fix[[paste0(i,"_",j)]] <- resPatc_fix
modd_fix[[paste0(i,"_",j)]] <- resPatd_fix

}

}


names(moda) <- c(1:length(moda))
moda.df <- data.frame(rbindlist(moda,idcol="name"))

names(modb) <- c(1:length(modb))
modb.df <- data.frame(rbindlist(modb,idcol="name"))

names(modc) <- c(1:length(modc))
modc.df <- data.frame(rbindlist(modc,idcol="name"))

names(modd) <- c(1:length(modd))
modd.df <- data.frame(rbindlist(modd,idcol="name"))

names(mode) <- c(1:length(mode))
mode.df <- data.frame(rbindlist(mode,idcol="name"))

names(modf) <- c(1:length(modf))
modf.df <- data.frame(rbindlist(modf,idcol="name"))


res2 <- stack(as.data.frame(t(res)))
colnames(res2) <- c("value","col")

# Add factor for number of variables in the models
res2$nvar <- NA
res2[res2$col == "Tensor product \n NSI and NPP mean" | res2$col == "Full Model",]$nvar <- "Multiple"
res2[is.na(res2$nvar),]$nvar <- "Single"

# Group and calculate median for ordering purposes
res2 <- res2 %>% group_by(col) %>% mutate(med = median(value))
res2_ord <- res2[order(-res2$med),]
res2_ord$col <- factor(res2_ord$col,levels=unique(res2_ord$col))
res2_ord$nvar <- factor(res2_ord$nvar,levels=c("Single","Multiple"))

# other version

violins <- ggplot(res2_ord) +
				facet_grid(~ nvar, scales = "free", space='free') +
				geom_violin(aes(col,value), fill='#A4A4A4', alpha=0.5) +
				scale_y_continuous(limits = quantile(res2_ord$value, c(.01, .99))) +
				theme_minimal() +
				xlab("Fixed effects") +
				ylab("RMSE") +
				geom_boxplot(aes(col,value), width=0.05) +
				theme(panel.grid.major.x = element_blank(),
					  strip.text.x = element_text(size = 13),
					  axis.title.x = element_text(size = 13),
					  axis.title.y = element_text(size = 13),
					  axis.text.x=element_text(size=11),
					  plot.background = element_rect(fill = "white"),
					  panel.background = element_rect(fill = NA,
													  color = "black"))

ggsave("RMSE models.png", plot = violins, device = png(), path = NULL,scale = 1, width = 9, height = 8, dpi = 300, limitsize = TRUE)
dev.off()

# Fixed NSI effect of NPP mean

# the smooth of the model with only NPP mean (to show non linearity and overall decrease of rich with high NPP mean)

modf.df$uprP <- modf.df$catPredf.fit + (2 * modf.df$catPredf.se.fit)
modf.df$lwrP <- modf.df$catPredf.fit - (2 * modf.df$catPredf.se.fit)

p_nppModel <- ggplot(data=modf.df) +
					theme_classic() +
					geom_smooth(aes(prim_mean,catPredf.fit,
									group=name),
								method=gam, formula=y~s(x),
								se=FALSE,
								color="grey",
								size=.5) +
					geom_smooth(aes(prim_mean,uprP),
								se=FALSE,
								color="black",
								size=.5,
								linetype= "dashed") +
					geom_smooth(aes(prim_mean,lwrP),
								se=FALSE,
								color="black",
								size=.5,
								linetype= "dashed") +
					ylab("Predicted #OTUs") +
					xlab("") +
					labs(title="\nNPP mean model") +
					theme(plot.title = element_text(size=11)) +
					scale_y_continuous(position = "left",
									   limits = c(0,800))


names(modc_fix) <- c(1:length(modc_fix))
modc_fix.df <- data.frame(rbindlist(modc_fix,idcol="name"))

names(modd_fix) <- c(1:length(modd_fix))
modd_fix.df <- data.frame(rbindlist(modd_fix,idcol="name"))

modc_fix.df$uprP <- modc_fix.df$catPredc_nsifix.fit + (2 * modc_fix.df$catPredc_nsifix.se.fit)
modc_fix.df$lwrP <- modc_fix.df$catPredc_nsifix.fit - (2 * modc_fix.df$catPredc_nsifix.se.fit)



p_c_fix <- ggplot(data=modc_fix.df) +
			theme_classic() +
			geom_smooth(aes(prim_mean,catPredc_nsifix.fit,
							colour=factor(nsi),
							group=name),
						method=gam, formula=y~te(x),
						se=FALSE,size=.05,alpha=0.95) +
			geom_smooth(aes(prim_mean,catPredc_nsifix.fit,
							color=factor(nsi),
							fill=factor(nsi)),
						method=gam, formula=y~te(x),
						se=FALSE,size=1.7,
						linetype="dashed") +
			ylab("") +
			xlab("NPP mean") +
			labs(title="Fixed NSI\nTensor product NSI*NPP mean") +
			theme(legend.position="none",
				  axis.text.y=element_blank(),
				  plot.title = element_text(size=11)) +
			scale_y_continuous(position = "left",
							   limits = c(0,800))

p_d_fix <- ggplot(data=modd_fix.df) +
			theme_classic() +
			geom_smooth(aes(prim_mean,catPredd_nsifix.fit,
							colour=factor(nsi),
							group=name),
						method=gam, formula=y~te(x),
						se=FALSE,size=.05,alpha=0.95) +
			geom_smooth(aes(prim_mean,catPredd_nsifix.fit,
							color=factor(nsi),
							fill=factor(nsi)),
						method=gam, formula=y~te(x),
						se=FALSE,size=1.7,
						linetype="dashed") +
			ylab("") +
			xlab("") +
			labs(title="Fixed NSI\nTensor product NSI*NPP mean + Physical extremes") +
			theme(legend.position="none",
				  axis.text.y=element_blank(),
				  plot.title = element_text(size=11)) +
			scale_y_continuous(position = "left",
							   limits = c(0,800))


pleg <- ggplot(data=modc_fix.df) +
			theme_minimal() +
			geom_smooth(aes(prim_mean,catPredc_nsifix.fit,
							color=factor(nsi)),
						se=FALSE,
						size=.5) +
			scale_colour_manual(values=c("#F8766D","#B79F00","#00BA38","#00BFC4","#F564E3"),
							  name="Fixed NSI \n   value")

leg <- ggpubr::get_legend(pleg)

# Included plot.new on first column to add more margin for the tilted x labels
smooth_plots <- cowplot::plot_grid(p_nppModel, p_c_fix, p_d_fix, leg,
								align = "hv",
								axis="b",
								ncol = 4,
								nrow = 1,
								rel_heights = c(1, 1, 1, 1),
								rel_widths = c(1, 1, 1, 0.3))

ggsave("NPP predictions.png", plot = smooth_plots, device = png(), path = NULL,scale = 1, width = 12, height = 4, dpi = 300, limitsize = TRUE)
dev.off()

