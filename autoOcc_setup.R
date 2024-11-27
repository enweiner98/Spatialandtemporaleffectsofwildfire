#### AUTO_OCC SET UP ####
# redoing statistical analysis beginning on 07/27/2023

#### LIBRARIES ####
library(tidyverse)
library(caret)
library(stats)
library(autoOcc)
library(beepr)
library(gt)

#### COVARIATE DATA ####
var_cov <- read.csv("var_cov_091423.csv")
# 110623 version with recoded tsf/burn status
con_cov <- read.csv("con_cov_083123.csv")
post_varCov <- read.csv("post_varCov_100223.csv")
# new covariate list:
    # burn_status: burned_post, unburned_post, prefire (ref level = burned_post)
    # time since fire: -22:22
    # ndvi: mean vegetation greenness at each site for a 500m (radius) buffer
    # ndvi_het: standard deviation of ndvi at each site for a 500m buffer
    # rdnbr: mean change in aboveground biomass at each site for a 500m buffer
    # rdnbrQuad: quadratic term for rdnbr
    # rdnbr_het: standard deviation of rdnbr at each site for a 500m buffer
    # dist_perim: distance to fire perimeter in kilometers
    # dist_perimQuad: quadratic term for distance to perimeter 
    # imperv_500: mean impervious surface at each site for a 500m buffer
    # microsite attractant (onroad): camera on natural features (0) or trail/trough (1)
    # county: Orange or Santa Cruz

#### TRANSFORMING COVARIATES ####
# centering and scaling covariates

# constant covariates
sub_conCov <- subset(con_cov,
                     select = c("dist_perim", "dist_perimQuad", "rdnbr_abs", "rdnbrQuad", "rdnbr_het", "imperv_500"))  # creating a subset of values

sub_cc_cs <- preProcess(sub_conCov, 
                            method = c("center", "scale"))  # estimating a transformation from training data

conCov_scaled <- predict(sub_cc_cs, sub_conCov)  # creating df with transformed version of covariates 
conCov_scaled <- round(conCov_scaled, 3)  # round to 2 decimal places

conCov_scaled <- conCov_scaled %>% 
  add_column(site = NA,
             county = NA,
             burn_status = NA,
             onroad = NA)
conCov_scaled$site <- con_cov$site
conCov_scaled$county <- con_cov$county
conCov_scaled$burn_status <- con_cov$burn_status
conCov_scaled$onroad <- con_cov$onroad

conCov_scaled <- conCov_scaled %>% 
  select(site, county, burn_status, onroad, dist_perim, dist_perimQuad, rdnbr_abs, rdnbrQuad, rdnbr_het, imperv_500)

conCov_scaled$county_bin <- ifelse(conCov_scaled$county == "Orange", 0, 1)
view(conCov_scaled)

# variable covariates 
sub_varCov <- subset(var_cov, 
                   select = c("ndvi_abs", "ndvi_het"))  # create a training dataset 

sub_vc_cs <- preProcess(sub_varCov, 
                          method = c("center", "scale"))  # estimate transformation from training data

varCov_scaled <- predict(sub_vc_cs, sub_varCov) # apply transformation 
varCov_scaled <- round(varCov_scaled, 3)

varCov_scaled <- varCov_scaled %>% 
  add_column(site = NA,
             season = NA,
             burn_status = NA,
             tsf = NA)
varCov_scaled$site <- var_cov$site
varCov_scaled$season <- var_cov$season
varCov_scaled$burn_status <- var_cov$burned
varCov_scaled$tsf <- var_cov$tsf

varCov_scaled <- varCov_scaled %>% 
  select(site, season, burn_status, tsf, ndvi_abs, ndvi_het)

view(varCov_scaled)

# post-fire variable covariates 
sub_pVC <- subset(post_varCov, 
                     select = c("rdnbr_abs", "rdnbrQuad", "rdnbr_het"))  # create a training dataset 

sub_pVC_cs <- preProcess(sub_pVC, 
                        method = c("center", "scale"))  # estimate transformation from training data

pVC_scaled <- predict(sub_pVC_cs, sub_pVC) # apply transformation 
pVC_scaled  <- round(pVC_scaled, 3)

pVC_scaled  <- pVC_scaled  %>% 
  add_column(site = NA,
             season = NA,
             burn_status = NA,
             dist_perim = NA,
             dist_perimQuad = NA)
pVC_scaled$site <- post_varCov$site
pVC_scaled$season <- post_varCov$season
pVC_scaled$burn_status <- post_varCov$burned
pVC_scaled$dist_perim <- post_varCov$dist_perim
pVC_scaled$dist_perimQuad <- post_varCov$dist_perimQuad

pVC_scaled <- pVC_scaled %>% 
  select(site, season, burn_status, dist_perim, dist_perimQuad, rdnbr_abs, rdnbrQuad, rdnbr_het)

view(pVC_scaled)

#### COVARIATE FRAME SET UP ####
# spatial analysis covariate frame: no change over time
covFrame_S <- list(
  burn_status = factor(conCov_scaled$burn_status,
                       levels = c("unburned", "burned")),
  dist_perim = conCov_scaled$dist_perim,
  dist_perimQuad = conCov_scaled$dist_perimQuad,
  rdnbr = conCov_scaled$rdnbr_abs,
  rdnbrQuad = conCov_scaled$rdnbrQuad,
  rdnbr_het = conCov_scaled$rdnbr_het,
  imperv = conCov_scaled$imperv_500,
  onroad = conCov_scaled$onroad,
  county = factor(conCov_scaled$county)
)
view(covFrame_S)

# spatial analysis covariate frame 2/3
  # covFrame_2: change over time for 2017_1-2022_4
  # covFrame_3: change over time for 2020_3-2022_4
covFrame_S3 <- list(
  burn_status = matrix(
    pVC_scaled$burn_status,
    nrow = length(unique(pVC_scaled$site)),
    ncol = length(unique(pVC_scaled$season))
  ),
  dist_perim = matrix(
    pVC_scaled$dist_perim,
    nrow = length(unique(pVC_scaled$site)),
    ncol = length(unique(pVC_scaled$season))
  ),
  rdnbr = matrix(
    pVC_scaled$rdnbr_abs,
    nrow = length(unique(pVC_scaled$site)),
    ncol = length(unique(pVC_scaled$season))
  ),
  rdnbrQuad = matrix(
    pVC_scaled$rdnbrQuad,
    nrow = length(unique(pVC_scaled$site)),
    ncol = length(unique(pVC_scaled$season))
  ),
  rdnbr_het = matrix(
    pVC_scaled$rdnbr_het,
    nrow = length(unique(pVC_scaled$site)),
    ncol = length(unique(pVC_scaled$season))
  ),
  imperv = conCov_scaled$imperv_500,
  onroad = conCov_scaled$onroad,
  county = conCov_scaled$county_bin
  )
)

# creating factors and assigning levels
covFrame_S3$burn_status <- as.data.frame(covFrame_S3$burn_status)

covFrame_S3$burn_status <- as.data.frame(
  lapply(
    covFrame_S3$burn_status,
    function(x) factor(x, levels = c("unburned", "burned")
    )
  )
)

# temporal analysis covariate frame
covFrame_T <- list(
  burn_status = matrix(
    varCov_scaled$burn_status,
    nrow = length(unique(varCov_scaled$site)),
    ncol = length(unique(varCov_scaled$season))
  ),
  tsf = matrix(
    varCov_scaled$tsf,
    nrow = length(unique(varCov_scaled$site)),
    ncol = length(unique(varCov_scaled$season))
  ),
  ndvi = matrix(
    varCov_scaled$ndvi_abs,
    nrow = length(unique(varCov_scaled$site)),
    ncol = length(unique(varCov_scaled$season))
  ),
#  ndvi_het = matrix(
#    varCov_scaled$ndvi_het,
#    nrow = length(unique(varCov_scaled$site)),
#    ncol = length(unique(varCov_scaled$season))
#  ),
  imperv = conCov_scaled$imperv_500,
  onroad = conCov_scaled$onroad,
  county = conCov_scaled$county_bin
)
view(varCov_scaled)

# creating factors and assigning levels
covFrame_T$burn_status <- as.data.frame(covFrame_T$burn_status)

covFrame_T$burn_status <- as.data.frame(
  lapply(
    covFrame_T$burn_status,
    function(x) factor(x, levels = c("burned_post", "unburned_post", "prefire")
    )
  )
)

#### COVARIATES ####
#### correlation matrices ####
# continuous temporal covariates: time since fire, ndvi, impervious surface
compare_T <- data.frame(
  tsf = var_cov$tsf,
  ndvi = var_cov$ndvi_abs,
  imperv = con_cov$imperv_500
)
cor(compare_T) 
# tsf:ndvi = -0.42
# ndvi:imperv = -0.29
# imperv:tsf = 0.16

compare_S <- data.frame(
  dist_perim = post_varCov$dist_perim,
  rdnbr = post_varCov$rdnbr_abs,
  rdnbrQ = post_varCov$rdnbrQuad,
  rdnbr_het = post_varCov$rdnbr_het,
  imperv = con_cov$imperv_500,
  onroad = con_cov$onroad
)
cor(compare_S)
# dp:rdnbr = -0.55
# dp:rdnbrHet = -0.37
# dp:imperv = -0.04
# dp:onroad = -0.06
# rdnbr:rdnbrHet = 0.89
# rdnbr:imperv = -0.057
# rdnbr:onroad = 0.14
# rdnbrHet:imperv = 0.008
# rdnbrHet:onroad = 0.20
# imperv:onroad = 0.008

compare_scaled <- data.frame(
  dist_perim = pVC_scaled$dist_perim,
  rdnbr = pVC_scaled$rdnbr_abs,
  rdnbrQuad = pVC_scaled$rdnbrQuad,
  rdnbr_het = pVC_scaled$rdnbr_het,
  imperv = conCov_scaled$imperv_500,
  onroad = con_cov$onroad
)
cor(compare_scaled) # didn't make a difference 

# rdnbr v. rdnbr heterogeneity
ggplot(compare_S, aes(rdnbr, rdnbr_het)) +
  geom_point() + 
  labs(x = "Fire Severity (Mean RdNBR)", y = "Fire Heterogeneity (SD RdNBR)") +
  theme_classic()

# rdnbr v. distance to perimeter
ggplot(compare_S, aes(rdnbr, dist_perim)) +
  geom_point() +
  labs(x = "Fire Severity (Mean RdNBR)", y = "Distance to Fire Perimeter (km)") +
  theme_classic() 

# tsf v. ndvi
ggplot(compare_T, aes(tsf, ndvi)) + 
  geom_point() +
  labs(x = "Time Since Fire", y = "NDVI") +
  theme_classic()

# comparing fire heterogeneity and ndvi heterogeneity
comp <- read.csv("fh_nh.csv")

ggplot(comp, aes(rdnbr_het, ndvi_het_2)) +
  geom_point() +
  theme_classic() +
  labs(x = "Fire Severity", y = "Vegetation Heterogeneity")


#### OCCUPANCY DATA ####
all_occu <- read.csv("allOccu_weekly_062323.csv")

all_occu <- dplyr::distinct(all_occu)  # remove duplicate sites/seasons

all_occu <- split(all_occu,
                  all_occu$Species)  # create a list split by species 

#### SUMMARY STATS ####
sum_stat <- data.frame(
  Species = c("Coyote", "Gray fox", "Bobcat", "Mountain lion", "Mule deer", "Striped skunk", "Opossum", "Raccoon", "Rabbit spp.", "CA ground squirrel"),
  t_det = c(1234, 1078, 1259, 740, 2367, 545, 190, 64, 361, 161),
  s_det = c(391, 370, 369, 195, 714, 205, 28, 4, 154, 61),
  t_occu = c(0.362, 0.391, 0.535, 0.510, 0.701, 0.301, 0.289, 0.236, 0.121, 0.070),
  t_det = c(0.488, 0.526, 0.397, 0.236, 0.532, 0.349, 0.178, 0.067, 0.507, 0.343),
  s_occu = c(0.468, 0.413, 0.533, 0.492, 0.693, 0.322, 0.128, NA, 0.162, 0.091),
  s_det = c(0.465, 0.510, 0.409, 0.234, 0.546, 0.359, 0.219, NA, 0.516, 0.345)
) 
sum_stat <- sum_stat %>% 
  arrange(Species)

sum_stat %>% 
  gt() %>% 
  cols_align(align = "left",
             columns = "Species") %>% 
  cols_align(align = "center",
             columns = c("s_det", "s_occu", "s_det.1", "t_det", "t_occu", "t_det.1")) %>% 
  tab_spanner(label = "Spatial Model", 
              columns = c("s_det", "s_occu", "s_det.1")) %>% 
  tab_spanner(label = "Temporal Model", 
              columns = c("t_det", "t_occu", "t_det.1")) %>% 
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_spanners(spanners = everything())) %>% 
  cols_label(t_det = "No. Detections",
             t_occu = "Overall Occupancy",
             t_det.1 = "Overall Detection",
             s_det = "No. Detections",
             s_occu = "Overall Occupancy",
             s_det.1 = "Overall Detection") %>% 
  cols_width(everything() ~ px(120))
  


