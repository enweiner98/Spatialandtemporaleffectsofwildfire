#### COYOTE FINAL ANALYSIS #####
# single-species, multi-factor models
# approach: temporal & spatial global models

#### LIBRARIES ####
library(tidyverse)
library(caret)
library(stats)
library(autoOcc)
library(beepr)
library(gt)
library(greekLetters)
library(cowplot)

#### OCCUPANCY DATA ####
# temporal
all_occu <- read.csv("allOccu_weekly_092723.csv") # complete pre/post dataset

all_occu <- dplyr::distinct(all_occu)  # remove duplicate sites/seasons

all_occu <- split(all_occu,
                  all_occu$Species)  # create a list split by species 

# spatial 
post_occu <- read.csv("postOccu_weekly_092723.csv") # complete post only dataset

post_occu <- dplyr::distinct(post_occu)

post_occu <- split(post_occu,
                   post_occu$Species)

# spatial v2
postOccu2 <- read.csv("postOccu_weekly_092923.csv")

postOccu2 <- dplyr::distinct(postOccu2)

postOccu2 <- split(postOccu2,
                   postOccu2$Species)

# spatial v3
postOccu3 <- read.csv("postOccu_weekly_100223.csv")

postOccu3 <- dplyr::distinct(postOccu3)

postOccu3 <- split(postOccu3,
                   postOccu3$Species)

#### COYOTE DATAFRAMES ####
coy_all <- all_occu$coyote  # all coyote detections
coy_post <- post_occu$coyote  # post-fire coyote detections
coy_post2 <- postOccu2$coyote
coy_post3 <- postOccu3$coyote

#### summary statistics ####
# temporal: total detections
coy_all %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 1234 total detections

# temporal: detections at burned v. unburned sites 
coy_bs$burn_status <- var_cov$burned

bs_count <- coy_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 370 detections 
# unburned: 452 detections 
# prefire: 422 detections

# temporal: total sites 
tot_sites <- coy_all %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
view(tot_sites)

# spatial: total detections
coy_post3 %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 391 total detections

# spatial: detections by burn status
coy_bs <- coy_post3
coy_bs$burn_status <- post_varCov$burned

bs_count <- coy_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 300 detections 
# unburned: 91 detections 

# spatial: total sites 
tot_sites <- coy_post3 %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

## creating dataframes for autoOcc
# post-fire only
coyote_S <- format_y(  # NAs for Santa Cruz sites from 2017-2020_3
  x = coy_post,
  site_column = "Site",
  time_column = "Season",
  history_column = "Week"
)

coyote_S2 <- format_y(  # NAs for all sites except Canyon 2 sites from 2017-2020_3
  x = coy_post2,
  site_column = "Site",
  time_column = "Season",
  history_column = "Week"
)

coyote_S3 <- format_y(  # detection histories for all sites begin at 2020_3
  x = coy_post3,
  site_column = "Site",
  time_column = "Season",
  history_column = "Week"
)

# includes pre-fire and post-fire data
coyote_T <- format_y(        
  x = coy_all,
  site_column = "Site",
  time_column = "Season",
  history_columns = "Week" 
)

#### SPATIAL ANALYSIS ####
#### null model ####
coy_null_S <- auto_occ(
  formula = ~1~1,
  y = coyote_S3
)

# overall expected occupancy: 0.468
(
  intercept_preds_psi <- predict(
    coy_null_S,
    type = "psi"
  )
)

# overall expected detection: 0.465
(
  intercept_preds_rho <- predict(
    coy_null_S, 
    type = "rho"
  )
)

#### prediction 1: fire heterogeneity ####
coyS_het <- auto_occ(
  formula = ~onroad + burn_status + rdnbr_het + burn_status:rdnbr_het
  ~burn_status + rdnbr_het + burn_status:rdnbr_het,
  y = coyote_S3,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(coyS_het)

#### prediction 2: fire severity ####
coyS_rdnbr <- auto_occ(
  formula = ~onroad + burn_status + rdnbr + rdnbrQuad
  ~burn_status + rdnbr + rdnbrQuad,
  y = coyote_S3,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(coyS_rdnbr)

#### prediction 3: distance to perimeter ####
coyS_dp <- auto_occ(
  formula = ~onroad + burn_status + dist_perim + burn_status:dist_perim
  ~burn_status + dist_perim + burn_status:dist_perim,
  y = coyote_S3,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(coyS_dp)

#### global model ####
coyS_global <- auto_occ(
  formula = ~onroad + burn_status + rdnbr_het + rdnbr + rdnbrQuad + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim
  ~burn_status + rdnbr_het + rdnbr + rdnbrQuad + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim,
  y = coyote_S3,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(coyS_global)

#### 2 factor: severity ####
coyS_sev_dp <- auto_occ(
  formula = ~onroad + burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim
  ~burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim,
  y = coyote_S3,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(coyS_sev_dp)

#### 2 factor: heterogeneity ####
coyS_het_dp <- auto_occ(
  formula = ~onroad + burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim
  ~burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim,
  y = coyote_S3,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(coyS_het_dp)

#### SPATIAL MODEL SELECTION ####
# competing predictions framework
compare_models(list(coy_null_S, coyS_het, coyS_rdnbr, coyS_dp, coyS_sev_dp, coyS_het_dp),
               digits = 2)  # best fit model: m4 - distance to fire perimeter

# global model framework
compare_models(list(coy_null_S, coyS_global),
               digits = 2) # best fit model: m2 - global

# 2-factor model framework
compare_models(list(coy_null_S, coyS_sev_dp, coyS_het_dp),
               digits = 2) # best fit model: m2 - severity + dp

beep()

#### SPATIAL MODEL PREDICTIONS ####
# burn_status, rdnbr_het, rdnbr, rdnbrQuad, dist_perim, dist_perimQuad 

#### fire heterogeneity ####
coy_rdnbrHet_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
                       ),
   
  rdnbr_het = rdnbrHet,
  rdnbr = 0, 
  rdnbrQuad = 0,
  dist_perim = 0,
  onroad = 0
)

# scaling
coy_rdnbrHet_scaled <- coy_rdnbrHet_real
coy_rdnbrHet_scaled$rdnbr_het <- (
  coy_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)
view(coy_rdnbrHet_scaled)

# the model prediction
rdnbrHet_coy <- predict(
  object = coyS_global,
  type = "psi",
  newdata = coy_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_coy <- data.frame(
  rdnbrHet_coy,
  coy_rdnbrHet_real
)

# plot 
ggplot(rdnbrHet_coy, aes(rdnbr_het, estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = "white", color = "gray") +
  geom_line(color = "blue") +
  annotate("rect", xmin=-90, xmax=1900, ymin=0.37, ymax=0.57, alpha=0.25, fill="gray") +
  labs(x = "Fire Heterogeneity", y = "Occupancy") + 
  ylim(0,1) +
  theme_classic(18)

#### burn_status*fire heterogeneity ####
rdnbrHet <- seq(20, 1750, 10)

# burned: lowest = 43.8, highest = 1749.5, median = 360.3, mean = 652.9
# unburned: lowest = 22.9, highest = 1365.1, median = 79.8, mean = 187.4

post_varCov %>% 
  filter(burned == "burned") %>% 
  summarize(quantile = quantile(rdnbr_het))  

ggplot(post_varCov, aes(burned, rdnbr_het)) + 
  geom_boxplot() + 
  theme_classic() +
  labs(x = "Burn Status", y = "Fire Heterogeneity")

coy_rdnbrHet_real <- data.frame(
  burn_status = c(
    rep("burned", 174),
    rep("unburned", 174)
  ),
  rdnbr_het = c(
    rdnbrHet,
    rdnbrHet
  ),
  rdnbr = 0,
  rdnbrQuad = 0,
  dist_perim = 0,
  onroad = 0
)

# scaling
coy_rdnbrHet_scaled <- coy_rdnbrHet_real
coy_rdnbrHet_scaled$rdnbr_het <- (
  coy_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_coy <- predict(
  object = coyS_het_dp,
  type = "psi",
  newdata = coy_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_coy <- data.frame(
  rdnbrHet_coy,
  coy_rdnbrHet_real
)

rdnbrHet_coy <- split(rdnbrHet_coy, rdnbrHet_coy$burn_status)

#plot
plot(1~1, type = "n", xlab = "Fire Heterogeneity", ylab = "Detection",
     xlim = c(20, 1750), ylim = c(0, 1), bty = "l", las = 1, cex.lab=1.5, cex.axis=1.25)

#burned
polygon(
  x = c(rdnbrHet_coy$burned$rdnbr_het, rev(rdnbrHet_coy$burned$rdnbr_het)),
  y = c(rdnbrHet_coy$burned$lower, rev(rdnbrHet_coy$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = rdnbrHet_coy$burned$rdnbr_het, y = rdnbrHet_coy$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(rdnbrHet_coy$unburned$rdnbr_het, rev(rdnbrHet_coy$unburned$rdnbr_het)),
  y = c(rdnbrHet_coy$unburned$lower, rev(rdnbrHet_coy$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = rdnbrHet_coy$unburned$rdnbr_het, y = rdnbrHet_coy$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("bottomleft", c("Unburned", "Burned"), lwd = 2, lty = c(2, 1), col = c("darkgreen", "darkorange"), bty = "n", cex=1.25)

c1 <- recordPlot()
c2 <- recordPlot()

c1 # occupancy
c2 # detection

#### fire severity ####
rdnbr_seq <- seq(-90, 1900, 10)
rdnbrQ <- seq(-0.60, 4, 0.1)

post_varCov %>% 
  filter(burned == "unburned") %>% 
  summarize(quantile = quantile(rdnbr_abs))  
# burned: -10 to 1911
# unburned: -94 to 630

ggplot(post_varCov, aes(burned, rdnbr_abs)) + 
  geom_boxplot() + 
  theme_classic() +
  labs(x = "Burn Status", y = "Fire Severity")

coy_rdnbr_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("unburned", "burned")
                       ),
  rdnbr_het = 0,
  rdnbr = rdnbr_seq,
  rdnbrQuad = 0,
  dist_perim = 0,
  onroad = 0
)

# scaling
coy_rdnbr_scaled <- coy_rdnbr_real
coy_rdnbr_scaled$rdnbr <- (
  coy_rdnbr_scaled$rdnbr - mean(post_varCov$rdnbr_abs)
) / sd(post_varCov$rdnbr_abs)

# the model prediction
rdnbr_coy <- predict(
  object = coyS_sev_dp,
  type = "rho",
  newdata = coy_rdnbr_scaled
)

# add on covariate data
rdnbr_coy <- data.frame(
  rdnbr_coy,
  coy_rdnbr_real
)

c4 <- ggplot(rdnbr_coy, aes(x = rdnbr, y = estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  geom_line(color = "darkorange") +
  labs(x = "Fire Severity (RdNBR)", y = "Detection", color = "County") +
  ylim(0,1) +
  scale_color_manual(values = c("darkorange", "darkgreen"),
                     labels = c("Orange", "Santa Cruz")) +
  theme_classic(18)

c3 # occupancy
c4 # detection

#### distance to fire perimeter ####
dp_seq <- seq(-6.2, 9.1, 0.1)
dpQuad <- seq(-0.50, 3.4, 0.1)

coy_dp_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  rdnbr_het = 0,
  rdnbr = 0,
  rdnbrQuad = 0,
  dist_perim = dp_seq,
  onroad = 0
)

# scaling
coy_dp_scaled <- coy_dp_real
coy_dp_scaled$dist_perim <- (
  coy_dp_scaled$dist_perim - mean(post_varCov$dist_perim)
) / sd(post_varCov$dist_perim)

# the model prediction
dp_coy <- predict(
  object = coyS_global,
  type = "psi",
  newdata = coy_dp_real
)

# add on covariate data
dp_coy <- data.frame(
  dp_coy,
  coy_dp_real
)

ggplot(dp_coy, aes(x = dist_perim, y = estimate)) +
  geom_line(color="darkorange") + 
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  scale_color_manual(values = c("darkorange", "darkgreen"),
                     labels = c("Orange", "Santa Cruz")) +
  labs(x = "Distance to Fire Perimeter (km)", y = "Occupancy", color = "County") + 
  ylim(0,1) +
  theme_classic(18)

#### burn_status*dist_perim ####
dp_seq <- seq(-6.3, 9.1, 0.1)

coy_dp_real <- data.frame(
  burn_status = c(
    rep("burned", 7),
    rep("unburned", 10)
  ),
  rdnbr_het = 0,
  rdnbr = 0,
  rdnbrQuad = 0,
  dist_perim = c(
    -6:0,
    0:9
  ),
  onroad = 0
)

# scaling
coy_dp_scaled <- coy_dp_real
coy_dp_scaled$dist_perim <- (
  coy_dp_scaled$dist_perim - mean(post_varCov$dist_perim)
) / sd(post_varCov$dist_perim)

# the model prediction
dp_coy <- predict(
  object = coyS_dp,
  type = "rho",
  newdata = coy_dp_real
)

dp_coy$dist_perim <- c(-6:0, 0:9)
dp_coy$burn_status <- coy_dp_real$burn_status

dp_coy <- split(dp_coy, coy_dp_real$burn_status)

#plot
plot(1~1, type = "n", xlab = "Distance to Fire Perimeter (km)", ylab = "Detection",
     xlim = c(-5, 5), ylim = c(0, 1), bty = "l", las = 1, cex.axis = 1.25, cex.lab=1.5)

#burned
polygon(
  x = c(dp_coy$burned$dist_perim, rev(dp_coy$burned$dist_perim)),
  y = c(dp_coy$burned$lower, rev(dp_coy$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = dp_coy$burned$dist_perim, y = dp_coy$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(dp_coy$unburned$dist_perim, rev(dp_coy$unburned$dist_perim)),
  y = c(dp_coy$unburned$lower, rev(dp_coy$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = dp_coy$unburned$dist_perim, y = dp_coy$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topright", c("Burned", "Unburned", "Fire Perimeter"), lwd = 2, lty = c(1, 2, 1), col = c("darkorange", "darkgreen", "black"), bty = "n", cex=1.25)

c5 <- recordPlot()
c6 <- recordPlot()

c5 # occupancy
c6 # detection

#### TEMPORAL ANALYSIS ####

#### null model ####
coy_null_T <- auto_occ(
  formula = ~1~1,
  y = coyote_T
)

# overall expected occupancy: 0.362
(
  intercept_preds_psi <- predict(
    coy_null_T,
    type = "psi"
  )
)

# overall expected detection: 0.488
(
  intercept_preds_rho <- predict(
    coy_null_T, 
    type = "rho"
  )
)

#### global model ####
coyT_global <- auto_occ(
  formula = ~onroad + burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi
      ~burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi,
  y = coyote_T,
  det_covs = covFrame_T,
  occ_covs = covFrame_T
)
summary(coyT_global)

# without burn status*heterogeneity AIC: 4556.477
# with burn status*heterogeneity AIC: 4560.477 

#### TEMPORAL MODEL SELECTION ####
compare_models(list(coy_null_T, coyT_global), digits = 2)

#### TEMPORAL MODEL PREDICTIONS ####
# formula = ~burn_status * tsf + ndvi + ndvi_het

#### burn status ####
burned_seq <- c("prefire", "burned_post", "unburned_post")
season_levels <- levels(cov_frame$season$V1)

coyT_burn_real <- data.frame(
  burn_status = factor(burned_seq, 
                       levels = burned_seq),
  tsf = -3,
  ndvi = 0, 
  ndvi_het = 0,
  onroad = 0
)

# the model prediction
burn_coyT <- predict(
  object = coyT_global,
  type = "psi",
  newdata = coyT_burn_real
)

# add on the covariate data
burn_coyT <- data.frame(
  burn_coyT,
  coyT_burn_real
)

# plot
ggplot(burn_coyT, aes(x = factor(burn_status,
                                levels = c("prefire", "burned_post", "unburned_post")), 
                     y = estimate,
                     color = burn_status)) +
  #geom_rect(aes(xmin=0, xmax=4, ymin=0.47, ymax=0.51), fill = "gray", color = NA, alpha = 0.25) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, show.legend = FALSE) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_vline(aes(xintercept = 1.5), color = "gray40", linetype = "longdash", size = 1) +
  labs(x = "Burn Status", y = "Detection") +
  scale_x_discrete(labels = c("Pre-Fire", "Burned", "Unburned")) +
  scale_color_manual(values = c("darkgreen", "darkorange", "darkgreen")) +
  ylim(0, 1) +
  theme_classic(18)

c7 # occupancy
c8 # detection

#### burn_status*tsf ####
coyT_int <- data.frame(
  burn_status = c(
    rep("burned_post", 23),
    rep("unburned_post", 23),
    rep("prefire", 23)
  ),
  tsf = c(
    0:22, 
    0:22,
    -22:0
  ),
  ndvi = 0,
  ndvi_het = 0,
  onroad = 0
)
coyT_int$burn_status <- factor(coyT_int$burn_status,
                              levels = c("burned_post", "unburned_post", "prefire"))

coyT_int_pred <- predict(coyT_global,
                        type = 'psi',
                        newdata = coyT_int,
                        interval = "confidence")

coyT_int_pred$tsf <- c(0:22, 0:22, -22:0)
coyT_int_pred$burn_status <- coyT_int$burn_status

coyT_int_pred <- split(coyT_int_pred, coyT_int_pred$burn_status)

#plot
plot(1~1, type = "n", xlab = "Seasons since fire", ylab = "Occupancy",
     xlim = c(-22, 23), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(coyT_int_pred$prefire$tsf, rev(coyT_int_pred$prefire$tsf)),
  y = c(coyT_int_pred$prefire$lower, rev(coyT_int_pred$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = coyT_int_pred$prefire$tsf, y = coyT_int_pred$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(coyT_int_pred$burned_post$tsf, rev(coyT_int_pred$burned_post$tsf)),
  y = c(coyT_int_pred$burned_post$lower, rev(coyT_int_pred$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = coyT_int_pred$burned_post$tsf, y = coyT_int_pred$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(coyT_int_pred$unburned_post$tsf, rev(coyT_int_pred$unburned_post$tsf)),
  y = c(coyT_int_pred$unburned_post$lower, rev(coyT_int_pred$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = coyT_int_pred$unburned_post$tsf, y = coyT_int_pred$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topright", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

c9 <- recordPlot()
c10 <- recordPlot()

c9 # occupancy
c10 # detection

#### ndvi ####
ndvi_seq <- seq(0, 0.30, 0.01)

coyT_ndvi_real <- data.frame(
  burn_status = factor("burned_post",
                       levels = c("prefire", "burned_post", "unburned_post")
  ),
  tsf = -3,
  ndvi = ndvi_seq,
  ndvi_het = 0,
  onroad = 0
)

# scaling
coyT_ndvi_scaled <- coyT_ndvi_real
coyT_ndvi_scaled$ndvi <- (
  coyT_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_coyT <- predict(
  object = coyT_global,
  type = "rho",
  newdata = coyT_ndvi_scaled
)

# add on covariate data
ndvi_coyT <- data.frame(
  ndvi_coyT,
  coyT_ndvi_real
)

# plot
ggplot(ndvi_coyT, aes(x = ndvi, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  scale_x_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)) +
  labs(x = "NDVI", y = "Occupancy") + 
  ylim(0,1) +
  theme_classic(18)

#### burn_status*ndvi ####
ndvi_seq <- seq(0, 0.35, 0.01)

coyT_ndvi_real <- data.frame(
  burn_status = c(
    rep("prefire", 36),
    rep("burned_post", 36),
    rep("unburned_post", 36)
  ),
  tsf = -3,
  ndvi = c(
    ndvi_seq,
    ndvi_seq,
    ndvi_seq
  ),
  onroad = 0
)

# scaling
coyT_ndvi_scaled <- coyT_ndvi_real
coyT_ndvi_scaled$ndvi <- (
  coyT_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)
view(coyT_ndvi_scaled)

# the model prediction
ndvi_coyT <- predict(
  object = coyT_global,
  type = "psi",
  newdata = coyT_ndvi_scaled
)

# add on covariate data
ndvi_coyT <- data.frame(
  ndvi_coyT,
  coyT_ndvi_real
)

ndvi_coyT <- split(ndvi_coyT, ndvi_coyT$burn_status)

#plot
plot(1~1, type = "n", xlab = "Vegetation Biomass", ylab = "Occupancy", xlim = c(0, 0.3), ylim = c(0,1), bty = "l", las = 1, cex.lab=1.5, cex.axis=1.25)

#prefire
polygon(
  x = c(ndvi_coyT$prefire$ndvi, rev(ndvi_coyT$prefire$ndvi)),
  y = c(ndvi_coyT$prefire$lower, rev(ndvi_coyT$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = ndvi_coyT$prefire$ndvi, y = ndvi_coyT$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(ndvi_coyT$burned_post$ndvi, rev(ndvi_coyT$burned_post$ndvi)),
  y = c(ndvi_coyT$burned_post$lower, rev(ndvi_coyT$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = ndvi_coyT$burned_post$ndvi, y = ndvi_coyT$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(ndvi_coyT$unburned_post$ndvi, rev(ndvi_coyT$unburned_post$ndvi)),
  y = c(ndvi_coyT$unburned_post$lower, rev(ndvi_coyT$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = ndvi_coyT$unburned_post$ndvi, y = ndvi_coyT$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topright", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.15)

c11 <- recordPlot()
c12 <- recordPlot()

c11 # occupancy
c12 # detection

#### not using - heterogeneity ####
ndviHet_seq <- seq(0, 0.15, 0.01)

coyT_ndviHet_real <- data.frame(
  burn_status = factor("burned_post",
                       levels = c("prefire", "burned_post", "unburned_post")
  ),
  tsf = -3,
  ndvi = 0,
  ndvi_het = ndviHet_seq,
  onroad = 0
)

# scaling
coyT_ndviHet_scaled <- coyT_ndviHet_real
coyT_ndviHet_scaled$ndvi_het <- (
  coyT_ndviHet_scaled$ndvi_het - mean(var_cov$ndvi_het)
) / sd(var_cov$ndvi_het)
view(coyT_ndviHet_scaled)

# the model prediction
ndviHet_coyT <- predict(
  object = coyT_global,
  type = "psi",
  newdata = coyT_ndviHet_scaled
)

# add on covariate data
ndviHet_coyT <- data.frame(
  ndviHet_coyT,
  coyT_ndviHet_real
)

# plot
ggplot(ndviHet_coyT, aes(x = ndvi_het, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Vegetation Heterogeneity", y = "Occupancy Probability", title = "Coyote") + 
  ylim(0,1) +
  theme_classic(18)

#### not using - burn_status*heterogeneity ####
ndviHet_seq <- seq(0, 0.15, 0.01)

ggplot(var_cov, aes(burned, ndvi_het)) + 
  geom_boxplot() +
  theme_classic() +
  labs(x = "Burn Status", y = "Vegetation Heterogeneity")

coyT_ndviHet_real <- data.frame(
  burn_status = c(
    rep("prefire", 16),
    rep("burned_post", 16),
    rep("unburned_post", 16)
  ),
  tsf = -3,
  ndvi = 0,
  ndvi_het = c(
    ndviHet_seq,
    ndviHet_seq,
    ndviHet_seq
  ),
  onroad = 0
)

# scaling
coyT_ndviHet_scaled <- coyT_ndviHet_real
coyT_ndviHet_scaled$ndvi_het <- (
  coyT_ndviHet_scaled$ndvi_het - mean(var_cov$ndvi_het)
) / sd(var_cov$ndvi_het)
view(coyT_ndviHet_scaled)

# the model prediction
ndviHet_coyT <- predict(
  object = coyT_global,
  type = "psi",
  newdata = coyT_ndviHet_scaled
)

# add on covariate data
ndviHet_coyT <- data.frame(
  ndviHet_coyT,
  coyT_ndviHet_real
)

ndviHet_coyT <- split(ndviHet_coyT, ndviHet_coyT$burn_status)

plot(1~1, type = "n", xlab = "Vegetation Heterogeneity", ylab = "Occupancy",
     xlim = c(0, 0.15), ylim = c(0, 1), bty = "l", las = 1)

#prefire
polygon(
  x = c(ndviHet_coyT$prefire$ndvi_het, rev(ndviHet_coyT$prefire$ndvi_het)),
  y = c(ndviHet_coyT$prefire$lower, rev(ndviHet_coyT$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = ndviHet_coyT$prefire$ndvi_het, y = ndviHet_coyT$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)


#burned
polygon(
  x = c(ndviHet_coyT$burned_post$ndvi_het, rev(ndviHet_coyT$burned_post$ndvi_het)),
  y = c(ndviHet_coyT$burned_post$lower, rev(ndviHet_coyT$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = ndviHet_coyT$burned_post$ndvi_het, y = ndviHet_coyT$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(ndviHet_coyT$unburned_post$ndvi_het, rev(ndviHet_coyT$unburned_post$ndvi_het)),
  y = c(ndviHet_coyT$unburned_post$lower, rev(ndviHet_coyT$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = ndviHet_coyT$unburned_post$ndvi_het, y = ndviHet_coyT$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("topleft", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n")

#### grid arrange ####
# work in progress - need to figure out how to scale the plots to fit properly
coy_grid <- plot_grid(c1, c3, c5, c2, c4, c6, 
                      ncol = 3, nrow = 2)

coy_grid
# mess with ggsave, customizable height and width
ggsave("coy_grid.jpg", coy_grid, width = 20, height = 14)
?ggsave

#### SUMMARY TABLES ####

#### spatial ####
# setup
param_h <- c(
  paste(greeks("psi"), "- Intercept"), paste(greeks("psi"), "- Burn Status (Burned)"), paste(greeks("psi"), "- Fire Heterogeneity"), 
  paste(greeks("psi"), "- Burned x Heterogeneity"), paste(greeks("psi"), "-", greeks("theta")), paste(greeks("rho"), "- Intercept"), paste(greeks("rho"), "- Microsite Attractant"),
  paste(greeks("rho"), "- Burn Status (Burned)"),  paste(greeks("rho"), "- Fire Heterogeneity"),  paste(greeks("rho"), "- Burned x Heterogeneity")
  )
param_r <- c(
  paste(greeks("psi"), "- Intercept"), paste(greeks("psi"), "- Burn Status (Burned)"), paste(greeks("psi"), "- Fire Severity (RdNBR)"), paste(greeks("psi"), "- RdNBR^2"),
  paste(greeks("psi"), "-", greeks("theta")), paste(greeks("rho"), "- Intercept"), paste(greeks("rho"), "- Microsite Attractant"), paste(greeks("rho"), "- Burn Status (Burned)"),
  paste(greeks("rho"), "- Fire Severity (RdNBR)"), paste(greeks("rho"), "- RdNBR^2")
  )
param_d <- c(
  paste(greeks("psi"), "- Intercept"), paste(greeks("psi"), "- Burn Status (Burned)"), paste(greeks("psi"), "- Distance to Fire Perimeter"), paste(greeks("psi"), "- Burned x Distance to Perimeter"),
  paste(greeks("psi"), "-", greeks("theta")), paste(greeks("rho"), "- Intercept"), paste(greeks("rho"), "- Microsite Attractant"), 
  paste(greeks("rho"), "- Burn Status (Burned)"), paste(greeks("rho"), "- Distance to Fire Perimeter"), paste(greeks("rho"), "- Burned x Distance to Perimeter")
  )

estimate_h <- round(coyS_het@estimates$Est, 3)
estimate_r <- round(coyS_rdnbr@estimates$Est, 3)
estimate_d <- round(coyS_dp@estimates$Est, 3)

SE_h <- round(coyS_het@estimates$SE, 3)
SE_r <- round(coyS_rdnbr@estimates$SE, 3)
SE_d <- round(coyS_dp@estimates$SE, 3)

p_h <- round(coyS_het@estimates$p, 3)
p_r <- round(coyS_rdnbr@estimates$p, 3)
p_d <- round(coyS_dp@estimates$p, 3)

# fire heterogeneity 
coyHet_tab <- bind_cols(param_h, estimate_h, SE_h, p_h)
coyHet_tab <- coyHet_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
coyHet_tab <- coyHet_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

ch <- coyHet_tab %>% 
  gt() %>% 
  cols_align(align = "left", 
             columns = "Parameter") %>%
  tab_style(style = cell_text(weight = "normal"),
            locations = cells_column_labels(c("Parameter",
                                              "Estimate",
                                              "SE",
                                              "p"))) %>% 
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(
              columns = p,
              rows = p < 0.05
            )) %>% 
  tab_spanner(label = md("Fire Heterogeneity Model (AIC = 1498.5)"),
              columns = everything())

ch %>% 
  gtsave(filename = "ch.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# fire severity
coyR_tab <- bind_cols(param_r, estimate_r, SE_r, p_r)
coyR_tab <- coyR_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
coyR_tab <- coyR_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

cr <- coyR_tab %>% 
  gt() %>% 
  cols_align(align = "left", 
             columns = "Parameter") %>%
  tab_style(style = cell_text(weight = "normal"),
            locations = cells_column_labels(c("Parameter",
                                              "Estimate",
                                              "SE",
                                              "p"))) %>% 
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(
              columns = p,
              rows = p < 0.05
            )) %>% 
  tab_spanner(label = md("Fire Severity Model (AIC = 1504.6)"),
              columns = everything())

cr %>%  
  gtsave(filename = "cr.html",
             path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# distance to fire perimeter
coyD_tab <- bind_cols(param_d, estimate_d, SE_d, p_d)
coyD_tab <- coyD_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
coyD_tab <- coyD_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

cd <- coyD_tab %>% 
  gt() %>% 
  cols_align(align = "left", 
             columns = "Parameter") %>%
  tab_style(style = cell_text(weight = "normal"),
            locations = cells_column_labels(c("Parameter",
                                              "Estimate",
                                              "SE",
                                              "p"))) %>% 
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(
              columns = p,
              rows = p < 0.05
            )) %>% 
  tab_spanner(label = md("Distance to Fire Perimeter Model (AIC = 1473.2)"),
              columns = everything())

cd %>% 
  gtsave(filename = "cd.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

#### temporal ####
parameter <- coyT_global@estimates$parameter
parameter <- c(
  paste(greeks("psi"), "- Intercept"), paste(greeks("psi"), "- Burn Status (Unburned)"), 
  paste(greeks("psi"), "- Burn Status (Pre-Fire)"), paste(greeks("psi"), "- Time Since Fire"),
  paste(greeks("psi"), "- NDVI"), paste(greeks("psi"), "- Unburned x Time Since Fire"), 
  paste(greeks("psi"), "- Pre-Fire x Time Since Fire"), paste(greeks("psi"), "- Unburned x NDVI"),
  paste(greeks("psi"), "- Pre-Fire x NDVI"), paste(greeks("psi"), "-", greeks("theta")), 
  paste(greeks("rho"), "- Intercept"), paste(greeks("rho"), "- Microsite Attractant"), paste(greeks("rho"), "- Burn Status (Unburned)"), 
  paste(greeks("rho"), "- Burn Status (Pre-Fire)"), paste(greeks("rho"), "- Time Since Fire"),
  paste(greeks("rho"), "- NDVI"), paste(greeks("rho"), "- Unburned x Time Since Fire"), 
  paste(greeks("rho"), "- Pre-Fire x Time Since Fire"), paste(greeks("rho"), "- Unburned x NDVI"),
  paste(greeks("rho"), "- Pre-Fire x NDVI")
  )

estimate <- coyT_global@estimates$Est
estimate <- round(estimate, 3)

SE <- coyT_global@estimates$SE
SE <- round(SE, 3)

p <- coyT_global@estimates$p
p <- round(p, 3)

coyGlobal_tab <- bind_cols(parameter, estimate, SE, p)
coyGlobal_tab <- coyGlobal_tab %>% 
  rename("Parameter" = ...1, 
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
coyGlobal_tab <- coyGlobal_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

cg <- coyGlobal_tab %>% 
  gt() %>% 
  cols_align(align = "left", 
             columns = "Parameter") %>%
  tab_style(style = cell_text(weight = "normal"),
            locations = cells_column_labels(c("Parameter",
                                              "Estimate",
                                              "SE",
                                              "p"))) %>% 
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(
              columns = p,
              rows = p < 0.05
            )) %>% 
  tab_spanner(label = md("Global Model"),
              columns = everything())
 
cg %>% 
  gtsave(filename = "cg.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")
