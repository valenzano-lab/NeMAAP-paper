Ne_surv <- read.csv("rstan/MAAP/multipleNe.csv")
Ne_surv <- Ne_surv[-c(1)]

# subset the df for the numeric values only
surv2 <- Ne_surv[,4:53]
age <- c(1:50)

# produce the cumulative product for survival across all rows
cumsurv2 <- t(apply(surv2, 1, 
                    function(x) cumprod(as.numeric(x)))
)

# I want to generate the median lifespan for
median.lifespan <- apply(cumsurv2, 1,
                         function(x) which.min(abs(x-0.5))
)

# Generate new column names for new data frame
repr <- ifelse(Ne_surv[,3] == "sexual", "sex", "asex")
Ne <- Ne_surv[,2]
mode <- Ne_surv[,1]
melif <- median.lifespan

# New data frame
df <- data.frame(Ne, repr, mode, melif)

# Subset df by mode
df_ap <- subset(df, mode == 'AP')
plot(log(df_ap$Ne), df_ap$melif)

# subset df_ap by repr
df_ap_sex <- df[(df$mode == 'AP') & (df$repr == 'sex'),]
plot(log(df_ap_sex$Ne), df_ap_sex$melif)

# a simple plot of median lifespan as a function of Ne AND repr
plot(log(df_ap$Ne), df_ap$melif, col= ifelse(df_ap$repr == 'sex', "red", "blue"), pch=20, bty="n")
# Note that I am using the log of population size - it's a magnitude measure, so it makes sense to transform this way

write.csv(df,file="~/rstan/MAAP/Nesurv.csv")

x <- df_ap_sex$Ne
y <- df_ap_sex$melif
N <- nrow(df_ap_sex)

data_surv <- list(x = x,
                  y = y,
                  N = N)

model_A <- stan(file = "rstan/MAAP/maap_Nesurv01.stan",
                data = data_surv, 
                iter = 1000,
                warmup = 500)

traceplot(model_A, 
          inc_warmup=TRUE)

#####RESOURCES#####
# https://widdowquinn.github.io/Teaching-Stan-Hierarchical-Modelling/10-partial_pooling_varying_slope_and_intercept.html
# https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html
# https://vasishth.github.io/bayescogsci/book/ch-complexstan.html
##### ##### ##### #

## Let's improve the model - First, let's standardize the inputs
df_Ne <- read.csv("./rstan/MAAP/Nesurv.csv", header=TRUE)
# I will separate the dataset based on AP and MA
df_Ne_AP <- df_Ne[which(df_Ne$mode == "AP"),]
df_Ne_MA <- df_Ne[which(df_Ne$mode == "MA"),]

df_Ne_AP_S <- df_Ne[(df_Ne$mode == 'AP') & (df_Ne$repr == 'sex'),]
df_Ne_MA_S <- df_Ne[(df_Ne$mode == 'MA') & (df_Ne$repr == 'sex'),]

# Then, I generate two new data frames for the AP and MA models
melif_sap <- (df_Ne_AP_S$melif - mean(df_Ne_AP_S$melif))/sd(df_Ne_AP_S$melif)
df_Ne_AP_S$melif_st <- melif_sap

melif_sma <- (df_Ne_MA_S$melif - mean(df_Ne_MA_S$melif))/sd(df_Ne_MA_S$melif)
df_Ne_MA_S$melif_st <- melif_sma

### Below I estimate intercept and slope for the AP, sexual model
x <- df_Ne_AP_S$Ne
y <- df_Ne_AP_S$melif_st
N <- nrow(df_Ne_AP_S)

data_surv <- list(x = x,
                  y = y,
                  N = N)

model_AP <- stan(file = "rstan/MAAP/maap_Nesurv01.stan",
                data = data_surv, 
                iter = 1000,
                warmup = 500)


traceplot(model_AP, 
          inc_warmup=TRUE)

post_AP <- extract(model_AP)

#################################################################


# plot(log(df_Ne_MA_S$Ne), df_Ne_MA_S$melif_st, col="brown", frame.plot = FALSE)
# abline(mean(post$beta), mean(post$alpha), #abline requires first intercept, then slope
#        col= "navy", lwd=3)

plot(log(df_Ne_AP_S$Ne), df_Ne_AP_S$melif_st, col="brown", frame.plot = FALSE, ylim=c(-2.5,2.5), ylab="standardized median lifespan", xlab="log(Ne)")

# plot(df_Ne_MA_S$Ne, df_Ne_MA_S$melif_st, col="brown", frame.plot = FALSE, ylim=c(-2.5,2.5), ylab="standardized median lifespan", xlab="log(Ne)", log="x")

for (i in 1:50)
  curve(post_AP$beta[i] + post_AP$alpha[i]*(x-mean(df_Ne_AP_S$melif_st)), 
         col= alpha("black", 0.3), add=TRUE)


### Below I estimate intercept and slope for the MA, sexual model
x <- df_Ne_MA_S$Ne
y <- df_Ne_MA_S$melif_st
N <- nrow(df_Ne_MA_S)

data_surv <- list(x = x,
                  y = y,
                  N = N)

model_MA <- stan(file = "rstan/MAAP/maap_Nesurv01.stan",
                 data = data_surv, 
                 iter = 1000,
                 warmup = 500)


traceplot(model_MA, 
          inc_warmup=TRUE)

post_MA <- extract(model_MA)

#################################################################

# plot(log(df_Ne_MA_S$Ne), df_Ne_MA_S$melif_st, col="brown", frame.plot = FALSE, ylim=c(-2.5,2.5), ylab="standardized median lifespan", xlab="log(Ne)")
# 
# for (i in 1:50)
#   curve(post_MA$beta[i] + post_MA$alpha[i]*(x-mean(df_Ne_MA_S$melif_st)), 
#         col= alpha("black", 0.3), add=TRUE)

#################################################################

# par(mfrow=c(1,2))
# 
# plot(log(df_Ne_AP_S$Ne), df_Ne_AP_S$melif_st, col="brown", frame.plot = FALSE, ylim=c(-2.5,2.5), ylab="standardized median lifespan", xlab="log(Ne)", main="Antagonistic Pleiotropy")
# abline(mean(post_AP$beta), mean(post_AP$alpha), #abline requires first intercept, then slope
#        col= "navy", lwd=3)
# 
# plot(log(df_Ne_MA_S$Ne), df_Ne_MA_S$melif_st, col="brown", frame.plot = FALSE, ylim=c(-2.5,2.5), ylab="standardized median lifespan", xlab="log(Ne)", main="Mutation Accumulation")
# abline(mean(post_MA$beta), mean(post_MA$alpha), #abline requires first intercept, then slope
#        col= "navy", lwd=3)

#################################################################

# par(mfrow=c(1,2))
# 
# plot(log(df_Ne_AP_S$Ne), df_Ne_AP_S$melif_st, col="brown", frame.plot = FALSE, ylim=c(-2.5,2.5), ylab="standardized median lifespan", xlab="log(Ne)", main="Antagonistic Pleiotropy")
# for (i in 1:50)
#   curve(post_AP$beta[i] + post_AP$alpha[i]*(x-mean(df_Ne_AP_S$melif_st)), 
#         col= alpha("black", 0.3), add=TRUE)
# 
# plot(log(df_Ne_MA_S$Ne), df_Ne_MA_S$melif_st, col="brown", frame.plot = FALSE, ylim=c(-2.5,2.5), ylab="standardized median lifespan", xlab="log(Ne)", main="Mutation Accumulation")
# for (i in 1:50)
#   curve(post_MA$beta[i] + post_MA$alpha[i]*(x-mean(df_Ne_MA_S$melif_st)), 
#         col= alpha("black", 0.3), add=TRUE)
# dev.off()

################################################################################
png(file="./rstan/MAAP/figures/Figure4.png",
    width=900, height=600)

tiff(file="./rstan/MAAP/figures/Figure4.tiff",
    width=900, height=600)

jpeg(file="./rstan/MAAP/figures/Figure4.jpeg",
     width=900, height=600)

png(file="./rstan/MAAP/figures/Figure4.png",
    width=1800, height=1200)

pdf(file="./rstan/MAAP/figures/Figure4.pdf",
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")

par(mfrow=c(1,2))

plot(log(df_Ne_MA_S$Ne), df_Ne_MA_S$melif_st, 
     col=alpha("brown", 0.5), 
     pch=19,
     cex=1.5,
     cex.lab=1.5,
     cex.axis=1.5,
     cex.main=1.5,
     frame.plot = FALSE, 
     ylim=c(-2.5,2.5),
     ylab="standardized median lifespan", 
     xlab="log(Ne)", 
     main="Mutation Accumulation")

for (i in 1:50)
  curve(post_MA$beta[i] + post_MA$alpha[i]*(x-mean(df_Ne_MA_S$melif_st)), 
        col= alpha("black", 0.3), add=TRUE)

plot(log(df_Ne_AP_S$Ne), df_Ne_AP_S$melif_st, 
     col=alpha("brown", 0.5), 
     pch=19, 
     cex=1.5,
     cex.lab=1.5,
     cex.axis=1.5,
     cex.main=1.5,
     frame.plot = FALSE, 
     ylim=c(-2.5,2.5), 
     ylab="standardized median lifespan", xlab="log(Ne)", 
     main="Antagonistic Pleiotropy")

for (i in 1:50)
  curve(post_AP$beta[i] + post_AP$alpha[i]*(x-mean(df_Ne_AP_S$melif_st)), 
        col= alpha("black", 0.3), add=TRUE)

#par(fig = c(0.02,0.3, 0.45, 0.95), new = T)
par(fig = c(0.22,0.48, 0.08, 0.55), new = T)
plot(density(post_MA$alpha), 
     frame.plot=FALSE, 
     col="navy", lwd=3,
     cex=1.5,
     cex.lab=1.5,
     cex.axis=1,
     main = "", 
     ylab="", 
     xlab="", 
     yaxt="n", 
     xlim=c(0.6,1.4))

abline(v=1.0, col="grey", lwd=2)
abline(v=median(post_MA$alpha), col="dark red", lwd=2, lty=2)

par(fig = c(0.72,0.98, 0.08, 0.55), new = T)
plot(density(post_AP$alpha), 
     frame.plot=FALSE, 
     col="navy", 
     cex=1.5,
     cex.lab=1.5,
     cex.axis=1,
     lwd=3,  
     main = "", 
     ylab="", 
     xlab="", 
     yaxt="n", 
     xlim=c(0.6,1.4))

abline(v=1.0, col="grey", lwd=2)
abline(v=median(post_AP$alpha), col="dark red", lwd=2, lty=2)

dev.off()

#######################################################################
####BELOW, THE ANALYSIS FOR THE ASEXUALLY REPRODUCING POPULATIONS######
#######################################################################

df_Ne_AP_A <- df_Ne[(df_Ne$mode == 'AP') & (df_Ne$repr == 'asex'),]
df_Ne_MA_A <- df_Ne[(df_Ne$mode == 'MA') & (df_Ne$repr == 'asex'),]

# Then, I generate two new data frames for the AP and MA models
melif_sap_ <- (df_Ne_AP_A$melif - mean(df_Ne_AP_A$melif))/sd(df_Ne_AP_A$melif)
df_Ne_AP_A$melif_st <- melif_sap_

melif_sma_ <- (df_Ne_MA_A$melif - mean(df_Ne_MA_A$melif))/sd(df_Ne_MA_A$melif)
df_Ne_MA_A$melif_st <- melif_sma_

### Below I estimate intercept and slope for the AP, asexual model
x <- df_Ne_AP_A$Ne
y <- df_Ne_AP_A$melif_st
N <- nrow(df_Ne_AP_A)

data_surv <- list(x = x,
                  y = y,
                  N = N)

model_AP_ <- stan(file = "rstan/MAAP/maap_Nesurv01.stan",
                 data = data_surv, 
                 iter = 1000,
                 warmup = 500)


traceplot(model_AP_, 
          inc_warmup=TRUE)

post_AP_ <- extract(model_AP_)
#################################################################

### Below I estimate intercept and slope for the MA, asexual model
x <- df_Ne_MA_A$Ne
y <- df_Ne_MA_A$melif_st
N <- nrow(df_Ne_MA_A)

data_surv <- list(x = x,
                  y = y,
                  N = N)

model_MA_ <- stan(file = "rstan/MAAP/maap_Nesurv01.stan",
                 data = data_surv, 
                 iter = 1000,
                 warmup = 500)


traceplot(model_MA_, 
          inc_warmup=TRUE)

post_MA_ <- extract(model_MA_)
#################################################################

png(file="./rstan/MAAP/figures/SupplFigure4.png",
    width=900, height=600)

tiff(file="./rstan/MAAP/figures/SupplFigure4.tiff",
     width=900, height=600)

jpeg(file="./rstan/MAAP/figures/SupplFigure4.jpeg",
     width=900, height=600)

pdf(file="./rstan/MAAP/figures/SupplFigure4.pdf",
    width = 8, height = 5, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")

par(mfrow=c(1,2))

plot(log(df_Ne_MA_A$Ne), df_Ne_MA_A$melif_st, 
     col=alpha("brown", 0.5), 
     pch=19, 
     cex.lab=1.5,
     cex.axis=1.5,
     cex.main=1.5,
     frame.plot = FALSE, 
     ylim=c(-2.5,2.5), 
     ylab="standardized median lifespan", 
     xlab="log(Ne)", 
     main="Mutation Accumulation")

for (i in 1:50)
  curve(post_MA_$beta[i] + post_MA_$alpha[i]*(x-mean(df_Ne_MA_A$melif_st)), 
        col= alpha("black", 0.3), add=TRUE)

plot(log(df_Ne_AP_A$Ne), df_Ne_AP_A$melif_st, 
     col=alpha("brown", 0.5), 
     pch=19, 
     cex.lab=1.5,
     cex.axis=1.5,
     cex.main=1.5,
     frame.plot = FALSE, 
     ylim=c(-2.5,2.5), 
     ylab="standardized median lifespan", 
     xlab="log(Ne)", 
     main="Antagonistic Pleiotropy")

for (i in 1:50)
  curve(post_AP_$beta[i] + post_AP_$alpha[i]*(x-mean(df_Ne_AP_A$melif_st)), 
        col= alpha("black", 0.3), add=TRUE)

#par(fig = c(0.02,0.3, 0.45, 0.95), new = T)
par(fig = c(0.22,0.48, 0.08, 0.55), new = T)

plot(density(post_MA_$alpha), 
     frame.plot=FALSE, 
     col="navy", 
     cex=1.5,
     cex.lab=1.5,
     cex.axis=1,
     lwd=3,  
     main = "", 
     ylab="", 
     xlab="", 
     yaxt="n", 
     xlim=c(0.6,1.4))

abline(v=1.0, col="grey", lwd=2)
abline(v=median(post_MA_$alpha), col="dark red", lwd=2, lty=2)

par(fig = c(0.72,0.98, 0.08, 0.55), new = T)

plot(density(post_AP_$alpha), 
     frame.plot=FALSE, 
     col="navy",
     cex=1.5,
     cex.lab=1.5,
     cex.axis=1,
     lwd=3,  
     main = "", 
     ylab="", 
     xlab="", 
     yaxt="n", 
     xlim=c(0.6,1.4))

abline(v=1.0, col="grey", lwd=2)
abline(v=median(post_AP_$alpha), col="dark red", lwd=2, lty=2)

dev.off()


