#### MULE DEER FINAL ANALYSIS #####
# single-species, multi-factor models
# approach: temporal & spatial global models

#### LIBRARIES ####
library(tidyverse)
library(caret)
library(stats)
library(autoOcc)
library(beepr)

#### OCCUPANCY DATA ####
# temporal
all_occu <- read.csv("allOccu_weekly_092723.csv") # complete pre/post dataset

all_occu <- dplyr::distinct(all_occu)  # remove duplicate sites/seasons

all_occu <- split(all_occu,
                  all_occu$Species)  # create a list split by species 

# spatial 
post_occu <- read.csv("postOccu_weekly_100223.csv")  # post only dataset v3

post_occu <- dplyr::distinct(post_occu)

post_occu <- split(post_occu,
                   post_occu$Species)

#### MULE DEER DATAFRAMES ####
md_all <- all_occu$mule_deer  # all mule deer detections
md_post <- post_occu$mule_deer  # post-fire mule deer detections

## summary statistics
# temporal: total detections
md_all %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 2367 total detections

# detections at burned v. unburned sites 
md_bs <- md_all
md_bs$burn_status <- var_cov$burned

bs_count <- md_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 576 detections 
# unburned: 1046 detections 
# prefire: 745 detections

# temporal: total sites 
tot_sites <- md_all %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 73 sites 

# spatial: total detections
md_post %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE))) # 714 total detections

# detections at burned v. unburned sites 
md_bs <- md_post
md_bs$burn_status <- post_varCov$burned

bs_count <- md_bs %>% 
  group_by(burn_status) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 

rowSums(bs_count[,2:6]) 
# burned: 438 detections 
# unburned: 276 detections 

# spatial: total sites 
tot_sites <- md_post %>% 
  group_by(Site) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

tot_sites$sum <- rowSums(tot_sites[,2:6])

tot_sites$status <- ifelse(tot_sites$sum > 0, 1, 0)
sum(tot_sites$status) # detected at 62 sites 

## creating dataframes for autoOcc
# post-fire only
muleDeer_S <- format_y(
  x = md_post,
  site_column = "Site",
  time_column = "Season",
  history_column = "Week"
)

# includes pre-fire and post-fire data
muleDeer_T <- format_y(        
  x = md_all,
  site_column = "Site",
  time_column = "Season",
  history_columns = "Week" 
)

#### SPATIAL ANALYSIS ####
#### null model ####
mdS_null <- auto_occ(
  formula = ~1~1,
  y = muleDeer_S
)

# overall expected occupancy: 0.693
(
  intercept_preds_psi <- predict(
    mdS_null,
    type = "psi"
  )
)

# overall expected detection: 0.546
(
  intercept_preds_rho <- predict(
    mdS_null, 
    type = "rho"
  )
)

#### p1: fire heterogeneity ####
mdS_het <- auto_occ(
  formula = ~onroad + burn_status * rdnbr_het
  ~burn_status * rdnbr_het,
  y = muleDeer_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(mdS_het)

#### p2: fire severity ####
mdS_rdnbr <- auto_occ(
  formula = ~onroad + burn_status + rdnbr + rdnbrQuad
  ~burn_status + rdnbr + rdnbrQuad,
  y = muleDeer_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(mdS_rdnbr)

#### p3: distance to fire perimeter ####
mdS_dp <- auto_occ(
  formula = ~onroad + burn_status * dist_perim
  ~burn_status * dist_perim,
  y = muleDeer_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(mdS_dp)

#### 2 factor: severity ####
mdS_sev_dp <- auto_occ(
  formula = ~onroad + burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim
  ~burn_status + rdnbr + rdnbrQuad + dist_perim + burn_status:dist_perim,
  y = muleDeer_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(mdS_sev_dp)

#### 2 factor: heterogeneity ####
mdS_het_dp <- auto_occ(
  formula = ~onroad + burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim
  ~burn_status + rdnbr_het + dist_perim + burn_status:rdnbr_het + burn_status:dist_perim,
  y = muleDeer_S,
  det_covs = covFrame_S3,
  occ_covs = covFrame_S3
)
summary(mdS_het_dp)

#### SPATIAL MODEL SELECTION ####
compare_models(list(mdS_null, mdS_het, mdS_rdnbr, mdS_dp, mdS_sev_dp, mdS_het_dp),
               digits = 2) # best fit model: distance to perimeter

#### SPATIAL MODEL PREDICTIONS ####

#### burn status ####
md_bs <- data.frame(
  county = 1,
  burn_status = factor(c("burned", "unburned"),
                       levels = c("unburned", "burned")),
  dist_perim = 0,
  onroad = 0
)

bs_md <- predict(
  object = mdS_dp,
  type = "psi",
  newdata = md_bs
)

# add on covariate data
bs_md <- data.frame(
  bs_md,
  md_bs
)

ggplot(bs_md, aes(x = factor(burn_status,
                               levels = c("unburned", "burned")), 
                    y = estimate,
                    color = burn_status)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, show.legend = FALSE) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_vline(aes(xintercept = 1.5), color = "gray40", linetype = "longdash", size = 1) +
  labs(x = "Burn Status", y = "Occupancy") +
  scale_x_discrete(labels = c("Burned", "Unburned")) +
  scale_color_manual(values = c("darkorange", "darkgreen")) +
  ylim(0, 1) +
  theme_classic(18)

#### fire heterogeneity ####
md_rdnbrHet_real <- data.frame(
  county = 1,
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  
  rdnbr_het = rdnbrHet,
  onroad = 0
)

# scaling
md_rdnbrHet_scaled <- md_rdnbrHet_real
md_rdnbrHet_scaled$rdnbr_het <- (
  md_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_md <- predict(
  object = mdS_het,
  type = "psi",
  newdata = md_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_md <- data.frame(
  rdnbrHet_md,
  md_rdnbrHet_real
)

# plot 
ggplot(rdnbrHet_md, aes(rdnbr_het, estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Fire Heterogeneity", y = "Occupancy") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*fire heterogeneity ####
md_rdnbrHet_real <- data.frame(
  county = 1,
  burn_status = c(
    rep("burned", 174),
    rep("unburned", 174)
  ),
  rdnbr_het = c(
    rdnbrHet,
    rdnbrHet
  ),
  onroad = 0
)

# scaling
md_rdnbrHet_scaled <- md_rdnbrHet_real
md_rdnbrHet_scaled$rdnbr_het <- (
  md_rdnbrHet_scaled$rdnbr_het - mean(post_varCov$rdnbr_het)
) / sd(post_varCov$rdnbr_het)

# the model prediction
rdnbrHet_md <- predict(
  object = mdS_het,
  type = "psi",
  newdata = md_rdnbrHet_scaled
)

# add on covariate data
rdnbrHet_md <- data.frame(
  rdnbrHet_md,
  md_rdnbrHet_real
)

rdnbrHet_md <- split(rdnbrHet_md, rdnbrHet_md$burn_status)

#plot
plot(1~1, type = "n", xlab = "Fire Heterogeneity", ylab = "Occupancy",
     xlim = c(20, 1750), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(rdnbrHet_md$burned$rdnbr_het, rev(rdnbrHet_md$burned$rdnbr_het)),
  y = c(rdnbrHet_md$burned$lower, rev(rdnbrHet_md$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = rdnbrHet_md$burned$rdnbr_het, y = rdnbrHet_md$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(rdnbrHet_md$unburned$rdnbr_het, rev(rdnbrHet_md$unburned$rdnbr_het)),
  y = c(rdnbrHet_md$unburned$lower, rev(rdnbrHet_md$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = rdnbrHet_md$unburned$rdnbr_het, y = rdnbrHet_md$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("bottomleft", c("Unburned", "Burned"), lwd = 2, lty = c(2, 1), col = c("darkgreen", "darkorange"), bty = "n", cex=1.25)

md1 <- recordPlot()
md2 <- recordPlot()

md1 # occupancy
md2 # detection

#### fire severity ####
county_seq <- c("Orange", "Santa_Cruz")

md_rdnbr_real <- data.frame(
  county = 1,
  burn_status = factor("burned",
                       levels = c("unburned", "burned")
  ),
  rdnbr = rdnbr_seq,
  rdnbrQuad = 0,
  dist_perim= 0,
  onroad = 0
)

# scaling
md_rdnbr_scaled <- md_rdnbr_real
md_rdnbr_scaled$rdnbr <- (
  md_rdnbr_scaled$rdnbr - mean(post_varCov$rdnbr_abs)
) / sd(post_varCov$rdnbr_abs)

# the model prediction
rdnbr_md <- predict(
  object = mdS_sev_dp,
  type = "rho",
  newdata = md_rdnbr_scaled
)

# add on covariate data
rdnbr_md <- data.frame(
  rdnbr_md,
  md_rdnbr_real
)

md4 <- ggplot(rdnbr_md, aes(x = rdnbr, y = estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  geom_line(color = "darkorange") +
  labs(x = "Fire Severity (RdNBR)", y = "Detection") +
  ylim(0,1) +
  theme_classic(18)

md3 # occupancy
md4 # detection

#### distance to fire perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)

md_dp_real <- data.frame(
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  dist_perim = dp_seq,
  rdnbr = 0,
  rdnbrQuad = 0,
  onroad = 0
)

# scaling
md_dp_scaled <- md_dp_real
md_dp_scaled$dist_perim <- (
  md_dp_scaled$dist_perim - mean(post_varCov$dist_perim)
) / sd(post_varCov$dist_perim)

# the model prediction
dp_md <- predict(
  object = mdS_sev_dp,
  type = "rho",
  newdata = md_dp_real
)

# add on covariate data
dp_md <- data.frame(
  dp_md,
  md_dp_real
)

ggplot(dp_md, aes(x = dist_perim, y = estimate)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = NA, color = "gray40", linetype = 2) +
  geom_line(color = "darkorange") +
  labs(x = "Distance to Fire Perimeter (km)", y = "Detection") + 
  ylim(0,1) +
  theme_classic(18)

#### distance to fire perimeter ####
dp_seq <- seq(-6.2, 9.1, 0.1)
dpQuad <- seq(-0.50, 3.4, 0.1)

md_dp_real <- data.frame(
  county = county_seq,
  burn_status = factor("burned",
                       levels = c("burned", "unburned")
  ),
  dist_perim = dp_seq,
  rdnbr = 0, 
  rdnbrQuad = 0,
  onroad = 0
)

# scaling
md_dp_scaled <- md_dp_real
md_dp_scaled$dist_perim <- (
  md_dp_scaled$dist_perim - mean(post_varCov$dist_perim)
) / sd(post_varCov$dist_perim)

# the model prediction
dp_md <- predict(
  object = mdS_sev_dp,
  type = "psi",
  newdata = md_dp_real
)

# add on covariate data
dp_md <- data.frame(
  dp_md,
  md_dp_real
)

ggplot(dp_md, aes(x = dist_perim, y = estimate, color=county)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper, color=county), fill = NA, linetype = 2) +
  scale_color_manual(values = c("darkorange", "darkgreen"),
                     labels = c("Orange", "Santa Cruz")) +
  labs(x = "Distance to Fire Perimeter (km)", y = "Occupancy", color = "County") + 
  ylim(0,1) +
  theme_classic(18)

#### burn status*distance to perimeter ####
dp_seq <- seq(-6.3, 9.1, 0.1)

md_dp_real <- data.frame(
  burn_status = c(
    rep("burned", 7),
    rep("unburned", 10)
  ),
  dist_perim = c(
    -6:0,
    0:9
  ),
  rdnbr = 0,
  rdnbrQuad = 0,
  onroad = 0
)

# the model prediction
dp_md <- predict(
  object = mdS_sev_dp,
  type = "rho",
  newdata = md_dp_real
)

dp_md$dist_perim <- c(-6:0, 0:9)
dp_md$burn_status <- md_dp_real$burn_status

dp_md <- split(dp_md, md_dp_real$burn_status)

#plot
plot(1~1, type = "n", xlab = "Distance to Fire Perimeter (km)", ylab = "Detection",
     xlim = c(-5, 5), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#burned
polygon(
  x = c(dp_md$burned$dist_perim, rev(dp_md$burned$dist_perim)),
  y = c(dp_md$burned$lower, rev(dp_md$burned$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = dp_md$burned$dist_perim, y = dp_md$burned$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned
polygon(
  x = c(dp_md$unburned$dist_perim, rev(dp_md$unburned$dist_perim)),
  y = c(dp_md$unburned$lower, rev(dp_md$unburned$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = dp_md$unburned$dist_perim, y = dp_md$unburned$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topright", c("Burned", "Unburned", "Fire Perimeter"), lwd = 2, lty = c(1, 2, 1), col = c("darkorange", "darkgreen", "black"), bty = "n", cex=1.25)

md5 <- recordPlot()
md6 <- recordPlot()

md5 # occupancy
md6 # detection

#### TEMPORAL ANALYSIS ####

#### null model ####
mdT_null <- auto_occ(
  formula = ~1~1,
  y = muleDeer_T
)

# overall expected occupancy: 0.701
(
  intercept_preds_psi <- predict(
    mdT_null,
    type = "psi"
  )
)

# overall expected detection: 0.532
(
  intercept_preds_rho <- predict(
    mdT_null, 
    type = "rho"
  )
)

#### global model ####
mdT_global <- auto_occ(
  formula = ~onroad + county + burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi
  ~county + burn_status + tsf + ndvi + burn_status:tsf + burn_status:ndvi,
  y = muleDeer_T,
  det_covs = covFrame_T,
  occ_covs = covFrame_T
)
summary(mdT_global)

#### TEMPORAL MODEL SELECTION ####
compare_models(list(mdT_null, mdT_global),
               digits = 2)

#### TEMPORAL MODEL PREDICTIONS ####

#### burn status ####
burned_seq <- c("prefire", "burned_post", "unburned_post")
season_levels <- levels(cov_frame$season$V1)

md_burn_real <- data.frame(
  county = 1,
  burn_status = factor(burned_seq, 
                       levels = burned_seq),
  tsf = -3,
  ndvi = 0, 
  ndvi_het = 0,
  onroad = 0
)

# the model prediction
burn_md <- predict(
  object = mdT_global,
  type = "psi",
  newdata = md_burn_real
)

# add on the covariate data
burn_md <- data.frame(
  burn_md,
  md_burn_real
)

# plot
md7 <- ggplot(burn_md, aes(x = factor(burn_status,
                               levels = c("prefire", "burned_post", "unburned_post")), 
                    y = estimate,
                    color = burn_status)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, show.legend = FALSE) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_vline(aes(xintercept = 1.5), color = "gray40", linetype = "longdash", size = 1) +
  labs(x = "Burn Status", y = "Occupancy") +
  scale_x_discrete(labels = c("Pre-Fire", "Burned", "Unburned")) +
  scale_color_manual(values = c("darkgreen", "darkorange", "darkgreen")) +
  ylim(0, 1) +
  theme_classic(18)

md7 # occupancy
md8 # detection

#### burn_status*tsf ####
md_int <- data.frame(
  county = 1,
  burn_status = c(
    rep("prefire", 23),
    rep("burned_post", 23),
    rep("unburned_post", 23)
  ),
  tsf = c(
    -22:0,
    0:22, 
    0:22
  ),
  ndvi = 0,
  onroad = 0
)
md_int$burn_status <- factor(md_int$burn_status,
                             levels = c("prefire", "burned_post", "unburned_post"))

md_int_pred <- predict(mdT_global,
                       type = 'rho',
                       newdata = md_int,
                       interval = "confidence")

md_int_pred$tsf <- c(-22:0, 0:22, 0:22)
md_int_pred$burn_status <- md_int$burn_status

md_int_pred <- split(md_int_pred, md_int_pred$burn_status)

#plot
plot(1~1, type = "n", xlab = "Seasons Since Fire", ylab = "Detection",
     xlim = c(-22, 23), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(md_int_pred$prefire$tsf, rev(md_int_pred$prefire$tsf)),
  y = c(md_int_pred$prefire$lower, rev(md_int_pred$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = md_int_pred$prefire$tsf, y = md_int_pred$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(md_int_pred$burned_post$tsf, rev(md_int_pred$burned_post$tsf)),
  y = c(md_int_pred$burned_post$lower, rev(md_int_pred$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = md_int_pred$burned_post$tsf, y = md_int_pred$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(md_int_pred$unburned_post$tsf, rev(md_int_pred$unburned_post$tsf)),
  y = c(md_int_pred$unburned_post$lower, rev(md_int_pred$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = md_int_pred$unburned_post$tsf, y = md_int_pred$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

abline(v = 0, lwd = 4)
legend("topright", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.25)

md9 <- recordPlot() # occupancy
md10 <- recordPlot() # detection

#### ndvi ####
ndvi_seq <- seq(0, 0.30, 0.01)

md_ndvi_real <- data.frame(
  county = 1,
  burn_status = factor("burned_post",
                       levels = c("prefire", "burned_post", "unburned_post")
  ),
  tsf = -3,
  ndvi = ndvi_seq,
  ndvi_het = 0,
  onroad = 0
)

# scaling
md_ndvi_scaled <- md_ndvi_real
md_ndvi_scaled$ndvi <- (
  md_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_md <- predict(
  object = mdT_global,
  type = "rho",
  newdata = md_ndvi_scaled
)

# add on covariate data
ndvi_md <- data.frame(
  ndvi_md,
  md_ndvi_real
)

# plot
ggplot(ndvi_md, aes(x = ndvi, y = estimate)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.30) +
  labs(x = "Vegetation Greenness", y = "Detection", title = "Mule Deer") + 
  ylim(0,1) +
  theme_classic(18)

#### burn_status*ndvi ####
ndvi_seq <- seq(0, 0.35, 0.01)

md_ndvi_real <- data.frame(
  county = 1,
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
md_ndvi_scaled <- md_ndvi_real
md_ndvi_scaled$ndvi <- (
  md_ndvi_scaled$ndvi - mean(var_cov$ndvi_abs)
) / sd(var_cov$ndvi_abs)

# the model prediction
ndvi_md <- predict(
  object = mdT_global,
  type = "rho",
  newdata = md_ndvi_scaled
)

# add on covariate data
ndvi_md <- data.frame(
  ndvi_md,
  md_ndvi_real
)

ndvi_md <- split(ndvi_md, ndvi_md$burn_status)

#plot
plot(1~1, type = "n", xlab = "Vegetation Biomass", ylab = "Detection",
     xlim = c(0, 0.3), ylim = c(0, 1), bty = "l", las = 1, cex.axis=1.25, cex.lab=1.5)

#prefire
polygon(
  x = c(ndvi_md$prefire$ndvi, rev(ndvi_md$prefire$ndvi)),
  y = c(ndvi_md$prefire$lower, rev(ndvi_md$prefire$upper)),
  border = NA,
  col = scales::alpha("gray40", 0.3)
)
lines(x = ndvi_md$prefire$ndvi, y = ndvi_md$prefire$estimate, lwd = 3, col = "gray40")
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$lower, lwd = 2, lty = 2)
#lines(x = my_pred2$prefire$x, y = my_pred2$prefire$upper, lwd = 2, lty = 2)

#burned_post
polygon(
  x = c(ndvi_md$burned_post$ndvi, rev(ndvi_md$burned_post$ndvi)),
  y = c(ndvi_md$burned_post$lower, rev(ndvi_md$burned_post$upper)),
  border = NA,
  col = scales::alpha("darkorange", 0.3)
)
lines(x = ndvi_md$burned_post$ndvi, y = ndvi_md$burned_post$estimate, lwd = 3, col = "darkorange")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$lower, lwd = 2, lty = 2, col = "red")
#lines(x = my_pred2$burned_post$x, y = my_pred2$burned_post$upper, lwd = 2, lty = 2, col = "red")

#unburned_post
polygon(
  x = c(ndvi_md$unburned_post$ndvi, rev(ndvi_md$unburned_post$ndvi)),
  y = c(ndvi_md$unburned_post$lower, rev(ndvi_md$unburned_post$upper)),
  border = NA,
  col = scales::alpha("darkgreen", 0.3)
)
lines(x = ndvi_md$unburned_post$ndvi, y = ndvi_md$unburned_post$estimate, lwd = 3, lty = 2, col = "darkgreen")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$lower, lwd = 2, lty = 2, col = "gray40")
#lines(x = my_pred2$unburned_post$x, y = my_pred2$unburned_post$upper, lwd = 2, lty = 2, col = "gray40")

legend("bottomright", c("Pre-Fire", "Unburned", "Burned"), lwd = 2, lty = c(1, 2, 1), col = c("gray40", "darkgreen", "darkorange"), bty = "n", cex=1.5)

md11 <- recordPlot() # occupancy
md12 <- recordPlot() # detection

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

estimate_h <- round(mdS_het@estimates$Est, 3)
estimate_r <- round(mdS_rdnbr@estimates$Est, 3)
estimate_d <- round(mdS_dp@estimates$Est, 3)

SE_h <- round(mdS_het@estimates$SE, 3)
SE_r <- round(mdS_rdnbr@estimates$SE, 3)
SE_d <- round(mdS_dp@estimates$SE, 3)

p_h <- round(mdS_het@estimates$p, 3)
p_r <- round(mdS_rdnbr@estimates$p, 3)
p_d <- round(mdS_dp@estimates$p, 3)

# fire heterogeneity 
mdHet_tab <- bind_cols(param_h, estimate_h, SE_h, p_h)
mdHet_tab <- mdHet_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
mdHet_tab <- mdHet_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

md_h <- mdHet_tab %>% 
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
  tab_spanner(label = md("Fire Heterogeneity Model (AIC = 2220.2)"),
              columns = everything())

md_h %>% 
  gtsave(filename = "md_h.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# fire severity
mdR_tab <- bind_cols(param_r, estimate_r, SE_r, p_r)
mdR_tab <- mdR_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
mdR_tab <- mdR_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

md_r <- mdR_tab %>% 
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
  tab_spanner(label = md("Fire Severity Model (AIC = 2206.7)"),
              columns = everything())

md_r %>%  
  gtsave(filename = "md_r.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

# distance to fire perimeter
mdD_tab <- bind_cols(param_d, estimate_d, SE_d, p_d)
mdD_tab <- mdD_tab %>% 
  rename("Parameter" = ...1,
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
mdD_tab <- mdD_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

md_d <- mdD_tab %>% 
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
  tab_spanner(label = md("Distance to Fire Perimeter Model (AIC = 2206.3)"),
              columns = everything())

md_d %>% 
  gtsave(filename = "md_d.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")

#### temporal ####
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

estimate <- mdT_global@estimates$Est
estimate <- round(estimate, 3)

SE <- mdT_global@estimates$SE
SE <- round(SE, 3)

p <- mdT_global@estimates$p
p <- round(p, 3)

mdGlobal_tab <- bind_cols(parameter, estimate, SE, p)
mdGlobal_tab <- mdGlobal_tab %>% 
  rename("Parameter" = ...1, 
         "Estimate" = ...2,
         "SE" = ...3,
         "p" = ...4)
mdGlobal_tab <- mdGlobal_tab %>% 
  arrange(desc(abs(Estimate))) %>% 
  arrange(p)

md_g <- mdGlobal_tab %>% 
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

md_g %>% 
  gtsave(filename = "md_g.html",
         path = "/Users/Erin/Documents/School/CSU Long Beach/Thesis/Thesis")


