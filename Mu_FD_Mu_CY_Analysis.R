################
## 03/2020 ##
# G. Montazeaud #
##################################################
# Code for the study "Multifaceted functional diversity for multifaceted crop yield: towards ecological assembly rules for varietal mixtures"
###################################################


#########################
### Packages loading  ###
#########################
library(RCurl)
library(stringr)
library(Hmisc)
library(RColorBrewer)
library(corrplot)
library(MASS)
library(MuMIn)
library(parallel)
library(glmulti)
library(leaps)
library(agricolae)
library(cluster)
library(gplots)


#########################
### Data loading  ###
#########################

## Trait indices (CWM and FD). Mixtures' CWMs and FDs are computed based on single-variety trait values (see Material & Methods in the main text)
trait_data <- read.csv(text = getURL("https://raw.githubusercontent.com/germmtz/multifaceted_FD_multifaceted_YLD/master/CWM_FD.csv"), sep=";",header=T)


## Performance variables (biomass yield, grain yield, etc). Mixtures' absolute and relative performance are reported with prefix "RAW_" and "RYT_", respectively. Only absolute performances are reported for single-variety plots.
perf_data <- read.csv(text = getURL("https://raw.githubusercontent.com/germmtz/multifaceted_FD_multifaceted_YLD/master/RAW_RYT.csv"), sep=";",header=T)


################
################
# TABLE summarizing trait variation in single-variety plots (Table 1) #
################
################

# Selecting  single-variety data
traitM <- trait_data[which(trait_data$assoc=="M"),]
row.names(traitM) <- traitM$focal

# Retaining  traits of interest. We just retain "CWM" variables which, in single-variety plots, correspond to the trait mean over replicated  measurements. All traits are retained except "filling time" which is not used in the manuscript.
traitM <- traitM[,c(4:22)]

# Computing summary statistics
min <- apply(traitM,2, min)
max <- apply(traitM,2, max)
mean <-apply(traitM,2, mean)
sd <-apply(traitM,2, sd)

# Pooling summary statistics within a single data-frame
table <- as.data.frame(do.call(cbind, list(min,max,mean,sd)))
colnames(table) <- c("min", "max","mean","sd")
table[,"CV"] <- table[,"sd"]/table[,"mean"]
table_char <- apply(table,2,function(x) {sprintf("%.2f",round(x,2))})
row.names(table_char) <-  str_sub(row.names(table),5)

## Outputting the summary table
write.csv(table_char, file="functionnal_diversity.csv",row.name=T)


######################
#### PAIWISE CORRELATIONS BETWEEN MONOCULTURE TRAITS (Supplementary Figure 5) ######
#####################

## Renaming variables
colnames(traitM) <- c(":Angle[aer]", ":Angle[rac]",":Diam[sem]",":SRL[sem]",":RTD[sem]",":RBI[sem]",":RLD[sem]",":Diam[adv]",":SRL[adv]",":RTD[adv]",":RBI[adv]",":RLD[adv]","Till. nb.","Ear. bio.","SLA","Leaf N","Height","Heading","Maturity")

## computing correlation matrix (pearson coefficients)
cormat <- rcorr(as.matrix(traitM), type="pearson")

## Outputting correlation heatmap
png("Correlations_monoc_traits.png", height=1500, width=1500, res=190)
cols <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(cormat$r, p.mat=cormat$P, sig.level=0.05, insig="blank", method="color", type="upper", order="hclust", col=cols(200), tl.col="black", tl.srt=45, diag=FALSE, tl.cex=1.2, cl.cex = 1.2,  number.cex=1.2,  addCoefasPercent = T, addCoef.col="black", mar=c(0,0,1,0))
mtext("Trait correlations (Pearson's r)",side = 3,line=1, font=2, cex=1.5)
dev.off()


################
################
#  RYT  DISTRIBUTION (Figure 2) #
################
################

## Retaining only data from mixtures
perfP <- droplevels(perf_data[which(perf_data$assoc=="P"),])

## Plotting RYT distributions for Grain Yield (GY), Biomass Yield (BY), and Grain Protein Content (GPC)

png("RYT_distribution.png", width = 1100, height = 400, res=150)

par(mfrow=c(1,3), mar=c(4,5,6,2))


##
hist(perfP[,"RYT_BY"], main="Biomass Yield (BY)", xlab=expression("RYT"["BY"]), ylab="Nb of observations", col="grey", border="white", axes=F, xlim=c(0,2), ylim=c(0,40), breaks=seq(0,2,by=0.05), cex.lab=1.3, cex.main=1.5)
axis(2, las=1, at=seq(0,40, by=10), cex.axis=1.3)
axis(1, at=seq(0,2,by=0.5), tick=T, xlab="", cex.axis=1.3, cex.lab=1.5)
abline(v=1, lty=2, lwd=2)
#mtext("Mixture (relative)", side=3, line=0.5, cex=1.2)
text(0.40,38, labels = bquote(atop(NA,atop(textstyle(mu==.(1.04)), textstyle(sigma==.(0.14))))), cex=1.3)
mtext(expression(bold(a)), side=3, line=3, at=-0.5, cex=1.1)
text(1.75,35,expression(bold("*")), cex=3)
##

##
hist(perfP[,"RYT_GY"], main="Grain Yield (GY)", xlab=expression("RYT"["GY"]), ylab="Nb of observations", col="grey", border="white", axes=F, xlim=c(0,2), ylim=c(0,40),breaks=seq(0,2,by=0.05), cex.lab=1.3, cex.main=1.5)
axis(2, las=1, seq(0,40, by=10), cex.axis=1.3)
axis(1, at=seq(0,2,by=0.5), tick=T, xlab="RYT", cex.axis=1.3)
abline(v=1, lty=2, lwd=2)
#mtext("Mixture (relative)", side=3, line=0.5, cex=1.2)
text(0.40,38, labels = bquote(atop(NA,atop(textstyle(mu==.(1.04)), textstyle(sigma==.(0.16))))), cex=1.3)
mtext(expression(bold(b)), side=3, line=3, at=-0.5, cex=1.1)
text(1.75,35,expression(bold("*")), cex=3)
##


##
hist(perfP[,"RYT_GPC"], main="Grain Protein\nConcentration (GPC)", xlab=expression("RYT"["GPC"]), ylab="Nb of observations", col="grey", border="white", axes=F, xlim=c(0,2), ylim=c(0,120),breaks=seq(0,2,by=0.05), cex.lab=1.3, cex.main=1.5)
axis(2, las=1, seq(0,120, by=20), cex.axis=1.3)
axis(1, at=seq(0,2,by=0.5), tick=T, xlab="RYT", cex.axis=1.3)
abline(v=1, lty=2, lwd=2)
#mtext("Mixture (relative)", side=3, line=0.5, cex=1.2)
text(0.40,112, labels = bquote(atop(NA,atop(textstyle(mu==.("1.00")), textstyle(sigma==.(0.04))))), cex=1.3)
mtext(expression(bold(c)), side=3, line=3, at=-0.5, cex=1.1)
##

dev.off()



#################
### RYT Multivariate analysis (Main text - Figure 3) ##
#################

######  GENERAL PARAMETERS ####
N <- 10L # number of models to be retained among the best models
w <- 0.2 # barplot width space
conversion <- 0.005 # conversion factor to define colors on a rgb scale between [0,1]
cols <- c(rgb(0,150*conversion,50*conversion),rgb(150*conversion,75*conversion,0),rgb(0,50*conversion,150*conversion),rgb(1,1,1))  # vector of colors to differenciate aboveground traits (green), belowground traits (brown), and phenological traits (blue)
hgt <- 600 # image height
wdt <- 850 # image width
rsl <- 120 # image resolution

# retaining only mixture trait data
traitP <- droplevels(trait_data[which(trait_data$assoc=="P"),])

## Merging perf and trait data
dataP <- merge(perfP[,c(1:3,17,21,25)], traitP, by=c("genotype_1","genotype_2","assoc"))


##
#### A. BIOMASS YIELD RYT ###
##

## Subsetting and scaling the data set 
data_reg <- as.data.frame(scale(dataP[,c(5,7:44)], center=T, scale = T))

## Running model selection from the full_model with all traits
full_model <- lm(RYT_BY ~ ., data=data_reg)
gl_multi <- glmulti(full_model, level = 1, crit = aicc, method="l", report = T)

## Retrieving model-averaged estimates based on the top N models and their unconditional sampling variance
effect_sizes <- coef(gl_multi, select=N, alphaIC=0.05, icmethod = "Standard")[seq(nrow(coef(gl_multi, select=N, alphaIC=0.05, icmethod = "Standard")),1,-1),]
effect_sizes <- as.data.frame(effect_sizes[-grep("(Intercept)",row.names(effect_sizes)),])
effect_sizes <- effect_sizes[order(gsub("CWM|FD","",row.names(effect_sizes)), decreasing = F),]
effect_sizes <- effect_sizes[order(effect_sizes$Importance, decreasing = T),]


## Retrieving AICcs and Akaike weights for the top N models. Note that the weights are not computed based on the top N models here, but based on the whole range of models evaluated by the glmulti() function. These weights will be computed again in the next step so that their sum is equal to 1.
weights <- weightable(gl_multi)
weights <- weights[c(1:N),]


## Loop to build a table with several information on the top N models (AICc, delta AICc, adjusted R², conditional estimates of predictor's effect, and re-computed weights)
var_names <- row.names(effect_sizes)

info_model <- data.frame(AICc=NA, delta=NA, R2_adj=NA)
for (i in 1:nrow(weights)) {
  info_model[i,"AICc"] <- weights[i,"aicc"] # AICc
  info_model[i,"delta"] <- weights[i,"aicc"]-min(weights[,"aicc"]) # delta AICc
  info_model[i,"R2_adj"] <- summary(gl_multi@objects[[i]])$adj.r.squared # Adjusted r²
  for (j in 1:length(var_names)){
    info_model[i,var_names[j]] <- coef(gl_multi@objects[[i]])[var_names[j]] # estimates
  }
}

Lik <- exp((-1/2)*info_model$delta)
info_model$weight <- Lik/sum(Lik) # weights

# Re-ordering columns
info_model <- info_model[,c(2,ncol(info_model),3:(ncol(info_model)-1))]

# Outputting model informations
write.csv(info_model, file="model_selection_mixture_RYT_BY.csv", row.names=F)

######
### Plotting results
########


## Creating an artifical "y" variables to plot results
effect_sizes$y <- rev(seq(0.5+w,(0.5+w)+(nrow(effect_sizes)-1)*(1+w),(1+w)))
#effect_sizes$y <- rev(seq(0.5+w,(0.5+w)+(nrow(effect_sizes)-1)*(1+w),(1+w)))

## Computing the upper and lower bounds of confidence intervals for model-averaged estimates
effect_sizes$ub <- effect_sizes[,"Estimate"]+effect_sizes[,"+/- (alpha=0.05)"]
effect_sizes$lb <- effect_sizes[,"Estimate"]-effect_sizes[,"+/- (alpha=0.05)"]

## Defining specific vector colors to match the traits retained in the model selection & model averaging procedure
colvec <- cols[c(3,1,2,2,1,2,2,3,1,1,1,2)]
bgvec <- cols[c(3,1,4,2,4,4,4,4,1,1,1,2)]
densvec <- c(-9,-9,20,-9,20,20,20,20,-9,-9,-9,-9)
anglevec <- rep(45,12)

png("Model_selection_mixture_RYT_BY.png", height=hgt, width=wdt, res=rsl)

par(mfrow=c(1,2))
par(mar=c(4,4,5,1), oma=c(0,3,1,1), lwd=3)

# 1st plot: estimates and confidence intervales
plot(y~Estimate, data=effect_sizes, xlim=c(-1,1), ylim=c(0,nrow(effect_sizes)+2), xlab="Standardized estimates", ylab="", axes=F, cex=1.3, type="n", cex.lab=1.3, cex.axis=1.3)

abline(v=0, lty=2)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$ub,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$lb,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
points(effect_sizes$Estimate, effect_sizes$y, pch=21, bg=bgvec, col=colvec, cex=1.6, lwd=1.5)

axis(2, las=1, labels=c("Heading", "Leaf N", expression("RBI"["sem"]), expression("RTD"["adv"]), "Till. nb.", expression("Diam"["adv"]), expression("Diam"["sem"]),"Maturity", expression("Angle"["aer"]), "Ear. bio.", "Till. nb.", expression("Diam"["sem"])), at=effect_sizes$y, tick=F, cex.lab=1.3, cex.axis=1.3) 

axis(1, at=seq(-1,1,by=0.5), labels=c("-1","-0.5","0","0.5","1"), cex.lab=1.3, cex.axis=1.3)

mtext(expression(bold("Biomass Yield RYT (RYT"["BY"]*")")), side=3, line=3, at=1.3, font=2, cex=1.7)
R2 <- round(mean(info_model$R2_adj),2)
mtext(bquote(bar(R²[adj])*" = "*.(R2)), side=3, line=0.8, at=1.3, font=2, cex=1.5)
mtext("a", side=3, line=3.2, at=-2, font=2, cex=1.5)


## 2nb plot: variable importance
par(mar=c(4,1,5,1))
barplot(rev(effect_sizes[,"Importance"]),horiz=T, col=rev(colvec), bg=rev(bgvec), xlab="Relative variable importance", space=w, width=1,density=rev(densvec), angle=rev(anglevec), ylim=c(0,nrow(effect_sizes)+2), cex.lab=1.3, cex.axis=1.3)

dev.off()


##
#### B. GRAIN YIELD ###
##


## Subsetting and scaling the data set 
data_reg <- as.data.frame(scale(dataP[,c(4,7:44)], center=T, scale = T))

## Running model selection from the full_model with all traits
full_model <- lm(RYT_GY ~ ., data=data_reg)
gl_multi <- glmulti(full_model, level = 1, crit = aicc, method="l", report = T)

## Retrieving model-averaged estimates based on the top N models and their unconditional sampling variance
effect_sizes <- coef(gl_multi, select=N, alphaIC=0.05, icmethod = "Standard")[seq(nrow(coef(gl_multi, select=N, alphaIC=0.05, icmethod = "Standard")),1,-1),]
effect_sizes <- as.data.frame(effect_sizes[-grep("(Intercept)",row.names(effect_sizes)),])
effect_sizes <- effect_sizes[order(gsub("CWM|FD","",row.names(effect_sizes)), decreasing = F),]
effect_sizes <- effect_sizes[order(effect_sizes$Importance, decreasing = T),]

## Retrieving the AICcs and Akaike weights of the top N models. Note that the weights are not computed based on the top N models here, but based on the whole range of models evaluated by the glmulti() function. These weights will be computed again in the next step so that their sum is equal to 1.
weights <- weightable(gl_multi)
weights <- weights[c(1:N),]


## Loop to build a table with several information on the top N models (AICc, delta AICc, adjusted R², conditional estimates of predictor's effect, and re-computed weights)
var_names <- row.names(effect_sizes)

info_model <- data.frame(AICc=NA, delta=NA, R2_adj=NA)
for (i in 1:nrow(weights)) {
  info_model[i,"AICc"] <- weights[i,"aicc"] # AICc
  info_model[i,"delta"] <- weights[i,"aicc"]-min(weights[,"aicc"]) # delta AICc
  info_model[i,"R2_adj"] <- summary(gl_multi@objects[[i]])$adj.r.squared # Adjusted r²
  for (j in 1:length(var_names)){
    info_model[i,var_names[j]] <- coef(gl_multi@objects[[i]])[var_names[j]] # estimates
  }
}

Lik <- exp((-1/2)*info_model$delta)
info_model$weight <- Lik/sum(Lik) # weights

# Re-ordering columns
info_model <- info_model[,c(2,ncol(info_model),3:(ncol(info_model)-1))]

# Outputting model informations
write.csv(info_model, file="model_selection_mixture_RYT_GY.csv", row.names=F)

######
### Plotting results
########


## Creating an artifical "y" variables to plot results
effect_sizes$y <- rev(seq(0.5+w,(0.5+w)+(nrow(effect_sizes)-1)*(1+w),(1+w)))

## Computing the upper and lower bounds of confidence intervals for model-averaged estimates
effect_sizes$ub <- effect_sizes[,"Estimate"]+effect_sizes[,"+/- (alpha=0.05)"]
effect_sizes$lb <- effect_sizes[,"Estimate"]-effect_sizes[,"+/- (alpha=0.05)"]

## Defining specific vector colors to match the traits retained in the model selection & model averaging procedure
colvec <- cols[c(1,2,2,1,1,2,3,2,2,2,1,3,1)]
bgvec <- cols[c(1,4,2,4,1,4,3,4,4,4,1,4,4)]
densvec <- c(-9,20,-9,20,-9,20,-9,20,20,20,-9,20,20)
anglevec <- rep(45,13)

png("Model_selection_mixture_RYT_GY.png", height=hgt, width=wdt, res=rsl)

par(mfrow=c(1,2))
par(mar=c(4,4,5,1), oma=c(0,3,1,1), lwd=3)

# 1st plot: estimates and confidence intervales
plot(y~Estimate, data=effect_sizes, xlim=c(-1,1), ylim=c(0,nrow(effect_sizes)+2), xlab="Standardized estimates", ylab="", axes=F, cex=1.3, type="n", cex.lab=1.3, cex.axis=1.3)
abline(v=0, lty=2)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$ub,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$lb,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
points(effect_sizes$Estimate, effect_sizes$y, pch=21, bg=bgvec, col=colvec, cex=1.6, lwd=1.5)
axis(2, las=1, labels=c("Ear. bio.", expression("RBI"["sem"]), expression("RTD"["adv"]),"Till. nb.", "Till. nb.", expression("Diam"["sem"]), "Heading", expression("Diam"["adv"]), expression("RTD"["adv"]), expression("SRL"["adv"]), "Leaf N", "Maturity", "SLA"), at=effect_sizes$y, tick=F, cex.lab=1.3, cex.axis=1.3) 

axis(1, at=seq(-1,1,by=0.5), labels=c("-1","-0.5","0","0.5","1"), cex.lab=1.3, cex.axis=1.3)


mtext(expression(bold("Grain Yield RYT (RYT"["GY"]*")")), side=3, line=3, at=1.3, font=2, cex=1.7)
R2 <- round(mean(info_model$R2_adj),2)
mtext(bquote(bar(R²[adj])*" = "*.(R2)), side=3, line=0.8, at=1.3, font=2, cex=1.5)
mtext("b", side=3, line=3.2, at=-2, font=2, cex=1.5)

## 2nb plot: variable importance
par(mar=c(4,1,5,1))
barplot(rev(effect_sizes[,"Importance"]),horiz=T, col=rev(colvec), bg=rev(bgvec), xlab="Relative variable importance", space=w, width=1, density=rev(densvec), angle=rev(anglevec),ylim=c(0,nrow(effect_sizes)+2), cex.lab=1.3, cex.axis=1.3)


dev.off()



##
#### C. GRAIN PROTEIN CONTENT ###
##


## Subsetting and scaling the data set 
data_reg <- as.data.frame(scale(dataP[,c(6:44)], center=T, scale = T))

## Running model selection from the full_model with all traits
full_model <- lm(RYT_GPC ~ ., data=data_reg)
gl_multi <- glmulti(full_model, level = 1, crit = aicc, method="l", report = T)

## Retrieving model-averaged estimates based on the top N models and their unconditional sampling variance
effect_sizes <- coef(gl_multi, select=N, alphaIC=0.05, icmethod = "Standard")[seq(nrow(coef(gl_multi, select=N, alphaIC=0.05, icmethod = "Standard")),1,-1),]
effect_sizes <- as.data.frame(effect_sizes[-grep("(Intercept)",row.names(effect_sizes)),])
effect_sizes <- effect_sizes[order(gsub("CWM|FD","",row.names(effect_sizes)), decreasing = F),]
effect_sizes <- effect_sizes[order(effect_sizes$Importance, decreasing = T),]

## Retrieving the AICcs and Akaike weights of the top N models. Note that the weights are not computed based on the top N models here, but based on the whole range of models evaluated by the glmulti() function. These weights will be computed again in the next step so that their sum is equal to 1.
weights <- weightable(gl_multi)
weights <- weights[c(1:N),]


## Loop to build a table with several information on the top N models (AICc, delta AICc, adjusted R², conditional estimates of predictor's effect, and re-computed weights)
var_names <- row.names(effect_sizes)

info_model <- data.frame(AICc=NA, delta=NA, R2_adj=NA)
for (i in 1:nrow(weights)) {
  info_model[i,"AICc"] <- weights[i,"aicc"] # AICc
  info_model[i,"delta"] <- weights[i,"aicc"]-min(weights[,"aicc"]) # delta AICc
  info_model[i,"R2_adj"] <- summary(gl_multi@objects[[i]])$adj.r.squared # Adjusted r²
  for (j in 1:length(var_names)){
    info_model[i,var_names[j]] <- coef(gl_multi@objects[[i]])[var_names[j]] # estimates
  }
}

Lik <- exp((-1/2)*info_model$delta)
info_model$weight <- Lik/sum(Lik) # weights

# Re-ordering columns
info_model <- info_model[,c(2,ncol(info_model),3:(ncol(info_model)-1))]

# Outputting model informations
write.csv(info_model, file="model_selection_mixture_RYT_GPC.csv", row.names=F)

######
### Plotting results
########


## Creating an artifical "y" variables to plot results
effect_sizes$y <- rev(seq(0.5+w,(0.5+w)+(nrow(effect_sizes)-1)*(1+w),(1+w)))

## Computing the upper and lower bounds of confidence intervals for model-averaged estimates
effect_sizes$ub <- effect_sizes[,"Estimate"]+effect_sizes[,"+/- (alpha=0.05)"]
effect_sizes$lb <- effect_sizes[,"Estimate"]-effect_sizes[,"+/- (alpha=0.05)"]

## Defining specific vector colors to match the traits retained in the model selection & model averaging procedure
colvec <- cols[c(1,2,3,1,2,3,2,2,2,2,1,3,1,2)]
bgvec <- cols[c(1,2,3,1,4,3,2,4,4,2,4,4,4,4)]
densvec <- c(-9,-9,-9,-9,20,-9,-9,20,20,-9,20,20,20,20)
anglevec <- rep(45,14)


png("Model_selection_mixture_RYT_GPC.png", height=hgt, width=wdt, res=rsl)

par(mfrow=c(1,2))
par(mar=c(4,4,5,1), oma=c(0,3,1,1), lwd=3)

# 1st plot: estimates and confidence intervales
plot(y~Estimate, data=effect_sizes, xlim=c(-1,1), ylim=c(0,nrow(effect_sizes)+2), xlab="Standardized estimates", ylab="", axes=F, cex=1.3, type="n", cex.lab=1.3, cex.axis=1.3)
abline(v=0, lty=2)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$ub,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$lb,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
points(effect_sizes$Estimate, effect_sizes$y, pch=21, bg=bgvec, col=colvec, cex=1.6, lwd=1.5)
axis(2, las=1, labels=c(expression("Angle"["aer"]),expression("Angle"["root"]),"Maturity", "Till. nb.", expression("Diam"["sem"]),"Heading",expression("RLD"["sem"]),expression("RTD"["adv"]), expression("SRL"["adv"]),expression("SRL"["sem"]),expression("Angle"["aer"]), "Heading","Height", expression("SRL"["sem"])), at=effect_sizes$y, tick=F, cex.lab=1.3, cex.axis=1.3) 
axis(1, at=seq(-1,1,by=0.5), labels=c("-1","-0.5","0","0.5","1"), cex.lab=1.3, cex.axis=1.3)


mtext(expression(bold("Grain Protein Concentration RYT (RYT"["GPC"]*")")), side=3, line=3, at=1.3, font=2, cex=1.7)
R2 <- round(mean(info_model$R2_adj),2)
mtext(bquote(bar(R²[adj])*" = "*.(R2)), side=3, line=0.8, at=1.3, font=2, cex=1.5)
mtext("c", side=3, line=3.2, at=-2, font=2, cex=1.5)


## 2nb plot: variable importance
par(mar=c(4,1,5,1))
barplot(rev(effect_sizes[,"Importance"]),horiz=T, col=rev(colvec), bg=rev(bgvec), xlab="Relative variable importance", space=w, width=1, density=rev(densvec), angle=rev(anglevec),ylim=c(0,nrow(effect_sizes)+2), cex.lab=1.3, cex.axis=1.3)

dev.off()



#################
### Analysis of Spike Number per m² (SNb) and Thousand Kernel Weight (TKW) (Supplementary Figure 7) ##
#################


## Merging yield and trait data
dataP <- merge(perfP[,c(1:3,20,23)], traitP, by=c("genotype_1","genotype_2","assoc"))


##
#### A. Spike Number per m² (SNb) ###
##

## Subsetting and scaling the data set 
data_reg <- as.data.frame(scale(dataP[,c(4,6:43)], center=T, scale = T))

## Running model selection from the full_model with all traits
full_model <- lm(RYT_SNb ~ ., data=data_reg)
gl_multi <- glmulti(full_model, level = 1, crit = aicc, method="l", report = T)

## Retrieving model-averaged estimates based on the top N models and their unconditional sampling variance
effect_sizes <- coef(gl_multi, select=N, alphaIC=0.05, icmethod = "Standard")[seq(nrow(coef(gl_multi, select=N, alphaIC=0.05, icmethod = "Standard")),1,-1),]
effect_sizes <- as.data.frame(effect_sizes[-grep("(Intercept)",row.names(effect_sizes)),])
effect_sizes <- effect_sizes[order(gsub("CWM|FD","",row.names(effect_sizes)), decreasing = F),]
effect_sizes <- effect_sizes[order(effect_sizes$Importance, decreasing = T),]


## Retrieving the AICcs and Akaike weights of the top N models. Note that the weights are not computed based on the top N models here, but based on the whole range of models evaluated by the glmulti() function. These weights will be computed again in the next step so that their sum is equal to 1.
weights <- weightable(gl_multi)
weights <- weights[c(1:N),]


## Loop to build a table with several information on the top N models (AICc, delta AICc, adjusted R², conditional estimates of predictor's effect, and re-computed weights)
var_names <- row.names(effect_sizes)

info_model <- data.frame(AICc=NA, delta=NA, R2_adj=NA)
for (i in 1:nrow(weights)) {
  info_model[i,"AICc"] <- weights[i,"aicc"] # AICc
  info_model[i,"delta"] <- weights[i,"aicc"]-min(weights[,"aicc"]) # delta AICc
  info_model[i,"R2_adj"] <- summary(gl_multi@objects[[i]])$adj.r.squared # Adjusted r²
  for (j in 1:length(var_names)){
    info_model[i,var_names[j]] <- coef(gl_multi@objects[[i]])[var_names[j]] # estimates
  }
}

Lik <- exp((-1/2)*info_model$delta)
info_model$weight <- Lik/sum(Lik) # weights

# Re-ordering columns
info_model <- info_model[,c(2,ncol(info_model),3:(ncol(info_model)-1))]

# Outputting model informations
write.csv(info_model, file="model_selection_mixture_RYT_SNb.csv", row.names=F)

######
### Plotting results
########


## Creating an artifical "y" variables to plot results
effect_sizes$y <- rev(seq(0.5+w,(0.5+w)+(nrow(effect_sizes)-1)*(1+w),(1+w)))
#effect_sizes$y <- rev(seq(0.5+w,(0.5+w)+(nrow(effect_sizes)-1)*(1+w),(1+w)))

## Computing the upper and lower bounds of confidence intervals for model-averaged estimates
effect_sizes$ub <- effect_sizes[,"Estimate"]+effect_sizes[,"+/- (alpha=0.05)"]
effect_sizes$lb <- effect_sizes[,"Estimate"]-effect_sizes[,"+/- (alpha=0.05)"]

## Defining specific vector colors to match the traits retained in the model selection & model averaging procedure
colvec <- cols[c(3,3,2,2,1,2,1,1,1,2,1)]
bgvec <- cols[c(3,4,4,2,1,2,1,4,1,4,1)]
densvec <- c(-9,20,20,-9,-9,-9,-9,20,-9,20,-9)
anglevec <- rep(45,11)

png("Model_selection_mixture_RYT_SNb.png", height=hgt, width=wdt, res=rsl)

par(mfrow=c(1,2))
par(mar=c(4,4,5,1), oma=c(0,3,1,1), lwd=3)

# 1st plot: estimates and confidence intervales
plot(y~Estimate, data=effect_sizes, xlim=c(-1,1), ylim=c(0,nrow(effect_sizes)+2), xlab="Standardized estimates", ylab="", axes=F, cex=1.3, type="n", cex.lab=1.3, cex.axis=1.3)

abline(v=0, lty=2)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$ub,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$lb,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
points(effect_sizes$Estimate, effect_sizes$y, pch=21, bg=bgvec, col=colvec, cex=1.6, lwd=1.5)

axis(2, las=1, labels=c("Heading", "Maturity", expression("RBI"["sem"]), expression("RTD"["adv"]),  expression("Angle"["aer"]), expression("Diam"["sem"]), "Height", "Till. nb.", "Leaf N", expression("RLD"["adv"]), "SLA"), at=effect_sizes$y, tick=F, cex.lab=1.3, cex.axis=1.3) 

axis(1, at=seq(-1,1,by=0.5), labels=c("-1","-0.5","0","0.5","1"), cex.lab=1.3, cex.axis=1.3)

mtext(expression(bold("Spike nb per m² RYT (RYT"["SNb"]*")")), side=3, line=3, at=1.3, font=2, cex=1.7)
R2 <- as.character(format(round(mean(info_model$R2_adj),2), nsmall=2))
mtext(bquote(bar(R²[adj])*" = "*.(R2)), side=3, line=0.8, at=1.3, font=2, cex=1.5)
mtext("a", side=3, line=3.2, at=-2, font=2, cex=1.5)


## 2nb plot: variable importance
par(mar=c(4,1,5,1))
barplot(rev(effect_sizes[,"Importance"]),horiz=T, col=rev(colvec), bg=rev(bgvec), xlab="Relative variable importance", space=w, width=1,density=rev(densvec), angle=rev(anglevec), ylim=c(0,nrow(effect_sizes)+2), cex.lab=1.3, cex.axis=1.3)

dev.off()





##
#### B. Thousand Kernel Weight (TKW)  ###
##


## Subsetting and scaling the data set 
data_reg <- as.data.frame(scale(dataP[,c(5:43)], center=T, scale = T))

## Running model selection from the full_model with all traits
full_model <- lm(RYT_TKW ~ ., data=data_reg)
gl_multi <- glmulti(full_model, level = 1, crit = aicc, method="l", report = T)

## Retrieving model-averaged estimates based on the top N models and their unconditional sampling variance
effect_sizes <- coef(gl_multi, select=N, alphaIC=0.05, icmethod = "Standard")[seq(nrow(coef(gl_multi, select=N, alphaIC=0.05, icmethod = "Standard")),1,-1),]
effect_sizes <- as.data.frame(effect_sizes[-grep("(Intercept)",row.names(effect_sizes)),])
effect_sizes <- effect_sizes[order(gsub("CWM|FD","",row.names(effect_sizes)), decreasing = F),]
effect_sizes <- effect_sizes[order(effect_sizes$Importance, decreasing = T),]

## Retrieving the AICcs and Akaike weights of the top N models. Note that the weights are not computed based on the top N models here, but based on the whole range of models evaluated by the glmulti() function. These weights will be computed again in the next step so that their sum is equal to 1.
weights <- weightable(gl_multi)
weights <- weights[c(1:N),]


## Loop to build a table with several information on the top N models (AICc, delta AICc, adjusted R², conditional estimates of predictor's effect, and re-computed weights)
var_names <- row.names(effect_sizes)

info_model <- data.frame(AICc=NA, delta=NA, R2_adj=NA)
for (i in 1:nrow(weights)) {
  info_model[i,"AICc"] <- weights[i,"aicc"] # AICc
  info_model[i,"delta"] <- weights[i,"aicc"]-min(weights[,"aicc"]) # delta AICc
  info_model[i,"R2_adj"] <- summary(gl_multi@objects[[i]])$adj.r.squared # Adjusted r²
  for (j in 1:length(var_names)){
    info_model[i,var_names[j]] <- coef(gl_multi@objects[[i]])[var_names[j]] # estimates
  }
}

Lik <- exp((-1/2)*info_model$delta)
info_model$weight <- Lik/sum(Lik) # weights

# Re-ordering columns
info_model <- info_model[,c(2,ncol(info_model),3:(ncol(info_model)-1))]

# Outputting model informations
write.csv(info_model, file="model_selection_mixture_RYT_TKW.csv", row.names=F)

######
### Plotting results
########


## Creating an artifical "y" variables to plot results
effect_sizes$y <- rev(seq(0.5+w,(0.5+w)+(nrow(effect_sizes)-1)*(1+w),(1+w)))

## Computing the upper and lower bounds of confidence intervals for model-averaged estimates
effect_sizes$ub <- effect_sizes[,"Estimate"]+effect_sizes[,"+/- (alpha=0.05)"]
effect_sizes$lb <- effect_sizes[,"Estimate"]-effect_sizes[,"+/- (alpha=0.05)"]

## Defining specific vector colors to match the traits retained in the model selection & model averaging procedure
colvec <- cols[c(2,2,3,2,2,1,1,3,1,2,1,1,2,2,2)]
bgvec <- cols[c(2,4,3,2,4,1,1,3,4,4,4,1,4,2,2)]
densvec <- c(-9,20,-9,-9,20,-9,-9,-9,20,20,20,-9,20,-9,-9)
anglevec <- rep(45,15)

png("Model_selection_mixture_RYT_TKW.png", height=hgt, width=wdt, res=rsl)

par(mfrow=c(1,2))
par(mar=c(4,4,5,1), oma=c(0,3,1,1), lwd=3)

# 1st plot: estimates and confidence intervales
plot(y~Estimate, data=effect_sizes, xlim=c(-1,1), ylim=c(0,nrow(effect_sizes)+2.5), xlab="Standardized estimates", ylab="", axes=F, cex=1.3, type="n", cex.lab=1.3, cex.axis=1.3)
abline(v=0, lty=2)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$ub,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$lb,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
points(effect_sizes$Estimate, effect_sizes$y, pch=21, bg=bgvec, col=colvec, cex=1.6, lwd=1.5)
axis(2, las=1, labels=c(expression("Angle"["root"]), expression("Diam"["sem"]), "Maturity", expression("RLD"["sem"]), expression("RTD"["adv"]),expression("Angle"["aer"]), "Ear. bio.", "Heading", "Ear. bio.",  expression("RTD"["sem"]), "Height","SLA", expression("RLD"["sem"]), expression("RLD"["adv"]), expression("RBI"["sem"])), at=effect_sizes$y, tick=F, cex.lab=1.3, cex.axis=1.3) 

axis(1, at=seq(-1,1,by=0.5), labels=c("-1","-0.5","0","0.5","1"), cex.lab=1.3, cex.axis=1.3)


mtext(expression(bold("Thousand kernel weight RYT (RYT"["TKW"]*")")), side=3, line=3, at=1.3, font=2, cex=1.7)
R2 <- as.character(format(round(mean(info_model$R2_adj),2), nsmall=2))
mtext(bquote(bar(R²[adj])*" = "*.(R2)), side=3, line=0.8, at=1.3, font=2, cex=1.5)
mtext("b", side=3, line=3.2, at=-2, font=2, cex=1.5)

## 2nb plot: variable importance
par(mar=c(4,1,5,1))
barplot(rev(effect_sizes[,"Importance"]),horiz=T, col=rev(colvec), bg=rev(bgvec), xlab="Relative variable importance", space=w, width=1, density=rev(densvec), angle=rev(anglevec),ylim=c(0,nrow(effect_sizes)+2.5), cex.lab=1.3, cex.axis=1.3)


dev.off()


#################
### Analysis of Grain Protein Deviation (GPD) (Supplementary Figure 8) ##
#################


## Merging yield and trait data
dataP <- merge(perfP[,c(1:3,29)], traitP, by=c("genotype_1","genotype_2","assoc"))

## Subsetting and scaling the data set 
data_reg <- as.data.frame(scale(dataP[,c(4:42)], center=T, scale = T))

## Running model selection from the full_model with all traits
full_model <- lm(RYT_GPD ~ ., data=data_reg)
gl_multi <- glmulti(full_model, level = 1, crit = aicc, method="l", report = T)

## Retrieving model-averaged estimates based on the top N models and their unconditional sampling variance
effect_sizes <- coef(gl_multi, select=N, alphaIC=0.05, icmethod = "Standard")[seq(nrow(coef(gl_multi, select=N, alphaIC=0.05, icmethod = "Standard")),1,-1),]
effect_sizes <- as.data.frame(effect_sizes[-grep("(Intercept)",row.names(effect_sizes)),])
effect_sizes <- effect_sizes[order(gsub("CWM|FD","",row.names(effect_sizes)), decreasing = F),]
effect_sizes <- effect_sizes[order(effect_sizes$Importance, decreasing = T),]


## Retrieving the AICcs and Akaike weights of the top N models. Note that the weights are not computed based on the top N models here, but based on the whole range of models evaluated by the glmulti() function. These weights will be computed again in the next step so that their sum is equal to 1.
weights <- weightable(gl_multi)
weights <- weights[c(1:N),]


## Loop to build a table with several information on the top N models (AICc, delta AICc, adjusted R², conditional estimates of predictor's effect, and re-computed weights)
var_names <- row.names(effect_sizes)

info_model <- data.frame(AICc=NA, delta=NA, R2_adj=NA)
for (i in 1:nrow(weights)) {
  info_model[i,"AICc"] <- weights[i,"aicc"] # AICc
  info_model[i,"delta"] <- weights[i,"aicc"]-min(weights[,"aicc"]) # delta AICc
  info_model[i,"R2_adj"] <- summary(gl_multi@objects[[i]])$adj.r.squared # Adjusted r²
  for (j in 1:length(var_names)){
    info_model[i,var_names[j]] <- coef(gl_multi@objects[[i]])[var_names[j]] # estimates
  }
}

Lik <- exp((-1/2)*info_model$delta)
info_model$weight <- Lik/sum(Lik) # weights

# Re-ordering columns
info_model <- info_model[,c(2,ncol(info_model),3:(ncol(info_model)-1))]

# Outputting model informations
write.csv(info_model, file="model_selection_mixture_RYT_GPD.csv", row.names=F)

######
### Plotting results
########


## Creating an artifical "y" variables to plot results
effect_sizes$y <- rev(seq(0.5+w,(0.5+w)+(nrow(effect_sizes)-1)*(1+w),(1+w)))
#effect_sizes$y <- rev(seq(0.5+w,(0.5+w)+(nrow(effect_sizes)-1)*(1+w),(1+w)))

## Computing the upper and lower bounds of confidence intervals for model-averaged estimates
effect_sizes$ub <- effect_sizes[,"Estimate"]+effect_sizes[,"+/- (alpha=0.05)"]
effect_sizes$lb <- effect_sizes[,"Estimate"]-effect_sizes[,"+/- (alpha=0.05)"]

## Defining specific vector colors to match the traits retained in the model selection & model averaging procedure
colvec <- cols[c(1,1,3,1,1,2,2,2,1,2,2,1)]
bgvec <- cols[c(1,4,3,4,1,4,2,2,1,2,2,4)]
densvec <- c(-9,20,-9,20,-9,20,-9,-9,-9,-9,-9,20)
anglevec <- rep(45,12)

png("Model_selection_mixture_RYT_GPD.png", height=hgt, width=wdt, res=rsl)

par(mfrow=c(1,2))
par(mar=c(4,4,5,1), oma=c(0,3,1,1), lwd=3)

# 1st plot: estimates and confidence intervales
plot(y~Estimate, data=effect_sizes, xlim=c(-1,1), ylim=c(0,nrow(effect_sizes)+2), xlab="Standardized estimates", ylab="", axes=F, cex=1.3, type="n", cex.lab=1.3, cex.axis=1.3)

abline(v=0, lty=2)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$ub,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
arrows(effect_sizes$Estimate,effect_sizes$y,effect_sizes$lb,effect_sizes$y, code=2, angle=90, length=0.03, col=colvec, lwd=2.5)
points(effect_sizes$Estimate, effect_sizes$y, pch=21, bg=bgvec, col=colvec, cex=1.6, lwd=1.5)

axis(2, las=1, labels=c(expression("Angle"["aer"]), "Height", "Maturity", "Till. nb.", "Till. nb.",expression("Diam"["adv"]), expression("SRL"["sem"]),  expression("Angle"["root"]), "SLA", expression("RBI"["sem"]), expression("Diam"["sem"]), "Ear. bio."), at=effect_sizes$y, tick=F, cex.lab=1.3, cex.axis=1.3) 

axis(1, at=seq(-1,1,by=0.5), labels=c("-1","-0.5","0","0.5","1"), cex.lab=1.3, cex.axis=1.3)

mtext(expression(bold("Grain protein deviation RYT (RYT"["GPD"]*")")), side=3, line=3, at=1.3, font=2, cex=1.7)
R2 <- as.character(format(round(mean(info_model$R2_adj),2), nsmall=2))
mtext(bquote(bar(R²[adj])*" = "*.(R2)), side=3, line=0.8, at=1.3, font=2, cex=1.5)
mtext("b", side=3, line=3.2, at=-2, font=2, cex=1.5)


## 2nb plot: variable importance
par(mar=c(4,1,5,1))
barplot(rev(effect_sizes[,"Importance"]),horiz=T, col=rev(colvec), bg=rev(bgvec), xlab="Relative variable importance", space=w, width=1,density=rev(densvec), angle=rev(anglevec), ylim=c(0,nrow(effect_sizes)+2), cex.lab=1.3, cex.axis=1.3)

dev.off()


#########################################
############ Plot legend (main text - Figure 3, Supplementary Figures 7 and 8) ################
#########################################

png("Model_selection_legend.png", height=200, width=1000, res=150)
par(mar=c(2,2,2,2), oma=c(1,3,1,1))
plot.new()
xlegend = -0.2
ylegend =3


legend(x=xlegend,y=ylegend,legend=c(" "," "), pch=c(19,NA),fill=c(NA,"black"), bty="n", col=c("black","black"), border=c(NA,"black"), xjust=0.2, x.intersp = c(2,0.2), pt.cex = 1.5, cex=1.5, xpd=NA)
legend(xlegend+0.02,ylegend-1.2, legend="CWM", bty="n", xpd=NA, cex=1.4)

legend(x=xlegend+0.3,y=ylegend,legend=c(" "," "), pch=c(21,NA),fill=c(NA,"black"), bty="n", col=c("black","black"), border=c(NA,NA), xjust=0.2, x.intersp = c(2,0.2), pt.cex = 1.5, cex=1.5, xpd=NA, density=c(0,20))
legend(xlegend+0.02+0.3,ylegend-1.2, legend=expression("F"["D"]), bty="n", xpd=NA, cex=1.4)

legend(x=xlegend+0.5,y=ylegend,legend=c(" "," "), fill=c(cols[1],cols[2]), bty="n", col=c("black","black"), border=c(NA,NA), xjust=0, cex=1.5, xpd=NA)
legend(x=xlegend+0.5,y=ylegend,legend=c("Aboveground","Belowground"), fill=c(NA,NA), bty="n", col=c(NA,NA), border=c(NA,NA), xjust=0,cex=1.4, xpd=NA)
# 
legend(x=xlegend+0.95,y=ylegend,legend=c(" "), fill=c(cols[3]), bty="n", col=c("black"), border=c(NA), xjust=0, cex=1.5, xpd=NA)
legend(x=xlegend+0.95,y=ylegend,legend=c("Phenology"), fill=c(NA), bty="n", col=c(NA), border=c(NA), xjust=0,cex=1.4, xpd=NA)


dev.off()




#########################################
##  Hierarchical classification of genotypes (Figure 4 and 5 in the manuscript) 
#########################################

## retaining only single-variety trait data and standardizing the genotype x trait matrix
trait <- droplevels(trait_data[which(trait_data$assoc=="M"), c(1,grep("CWM",colnames(trait_data)))])
rownames(trait) <- trait$genotype_1
trait_classif <- scale(trait[,-which(colnames(trait)=="genotype_1")])


## Clustering genotypes with the Ward criterion
classif <- agnes(trait_classif, method="ward")    

## plotting the height of tree branches
png("Functional_groups_branches_height.png", width = 2000, height = 2000, res=300)
classif2 <- as.hclust(classif)
plot(rev(classif2$height), type="h", ylab="branch height")
dev.off()


## plotting the clustering tree with a delimitation of the 3 functional groups (Supplementary Figure 4)
png("Functional_groups_tree.png", width = 5000, height = 2500, res=300)
par(mar=c(2,5,3,2))
plot(classif, xlab="", ask=F, which.plots=2, cex=0.5, main="Functional Groups", las=1, ylab="Height", cex.lab=2, cex.axis=2, cex.main=4)
classes <- cutree(classif, k=3)
classes <- as.factor(classes)

## inverting groups "2" and "3" to be in line with the name of the groups in the previous version of the manuscript (the input file as been reformated between the two versions)
levels(classes)[levels(classes)=="2"] <- "4"
levels(classes)[levels(classes)=="3"] <- "2"
levels(classes)[levels(classes)=="4"] <- "3"
classes = factor(classes,levels(classes)[c(1,3,2)])

rect.hclust(classif, k=3, border="red")
text(35, 1.9,"1", cex=6, col="red")
text(96, 1.9,"3", cex=6, col="red")
text(154, 1.9,"2", cex=6, col="red")
dev.off()


## Adding the grouping factor to the original genotype x trait matrix, and computing trait means and variances per functional groups
trait_classif <- cbind.data.frame(trait_classif, classes)
colnames(trait_classif)[ncol(trait_classif)] <- "grp"
trait_means <- aggregate(.~grp, data=trait_classif, FUN= function(x) {round(mean(x),3)})
trait_variance <- aggregate(.~grp, data=trait_classif, FUN= function(x) {round(var(x),3)} )


## Ploting heatmap representing the mean trait values per functional group (Figure 4)

# re-framing the table of trait means
t_trait_means <- t(trait_means[,-1])
colnames(t_trait_means) <- c("1","2","3")

# Defining a colour gradient
colfunc <- colorRampPalette(c("#0f9b0f","#000000"))

# plotting heatmap
png("Functional_group_mean_heatmap.png",width = 2400, height = 2000, res=400)
hmmeans <- heatmap.2(as.matrix(t_trait_means),Rowv = T, Colv = F, dendrogram = "none", trace="none", col = colfunc(4), scale = "row", distfun = dist, hclustfun = hclust,xlab = "", ylab = "", cexRow = 1.3, cexCol = 2, srtRow = 0, srtCol = 0, density.info = "none",margins=c(5,15), key.ylab=NA, key.xlab=NA, key.title="", keysize=1, labRow =  c(expression(Angle["aer"]), expression(Angle["root"]), expression("Diam"["sem"]), expression("SRL"["sem"]),expression("RTD"["sem"]), expression("RBI"["sem"]),  expression("RLD"["sem"]), expression("Diam"["adv"]),expression("SRL"["adv"]),expression("RTD"["adv"]), expression("RBI"["adv"]),expression("RLD"["adv"]), "Till. nb.","Ear. bio.", "SLA", "Leaf N","Height","Heading", "Maturity"), labCol="",main = "Standardized\ntrait means", key=T)
dev.off()


##### Merging clustering information & RYT data
for (i in c(1:nrow(perfP))) {
  geno1 <- as.character(perfP[i,"genotype_1"])
  geno2 <- as.character(perfP[i,"genotype_2"])
  perfP[i,"fun_grp"] <- paste(min(as.numeric(trait_classif[geno1,"grp"]), as.numeric(trait_classif[geno2, "grp"])),max(as.numeric(trait_classif[geno1,"grp"]), as.numeric(trait_classif[geno2, "grp"])), sep="-")
}
perfP$fun_grp <- as.factor(perfP$fun_grp)


## Plotting RYT per functional associations (main text - Figure 5)
png("Functional_groups.png", height=400, width=1300, res=110)

# Defining some graphical parameters
par(mar=c(4,6,6,2), oma=c(1,1,1,1), mfrow=c(1,3))
size_axis<- 1.7
size_lab <- 1.8
size_main <-1.6
size_letter <- 1.5

## A. Biomass yield RYT
data_reg <- perfP[,which(colnames(perfP)%in%c("RYT_BY", "fun_grp"))]
by(data_reg, data_reg$fun_grp, FUN = function(x) {t.test(x[,1], mu=1)$p.value})
mod <- lm(RYT_BY ~ fun_grp, data=data_reg)
multi.test <- LSD.test(mod,"fun_grp")
multi.test <- multi.test$groups[match(levels(data_reg$fun_grp),row.names(multi.test$groups)),]

boxplot(RYT_BY~fun_grp, data=data_reg, ylab="", xlab="Functional association",cex.axis=size_axis, cex=2, cex.lab=size_lab, las=1, col="grey75", ylim=c(0.5,1.7))
abline(h=1, lty=2, lwd=2)
mtext(side=3, line=2, "Biomass yield RYT", font=2, cex=size_main, at=3.5)
mtext(side=2, line=4, expression("RYT"["BY"]),cex=1.5)
mtext(side=3, line=4, "a", font=2, cex=size_letter, at=-1.5)
axis(3, as.character(multi.test$groups), at=c(1:6), tick=F, line=-3, font=2, cex.axis=size_axis)
text(c(4,5),c(1.105,1.13),"*",cex=3, font=2)


## B. Grain yield RYT
data_reg <- perfP[,which(colnames(perfP)%in%c("RYT_GY", "fun_grp"))]
by(data_reg, data_reg$fun_grp, FUN = function(x) {t.test(x[,1], mu=1)$p.value})
mod <- lm(RYT_GY ~ fun_grp, data=data_reg)
multi.test <- LSD.test(mod,"fun_grp")
multi.test <- multi.test$groups[match(levels(data_reg$fun_grp),row.names(multi.test$groups)),]

boxplot(RYT_GY~fun_grp, data=data_reg, ylab="", xlab="Functional association", cex.axis=size_axis, cex=2, cex.lab=size_lab, las=1, col="grey75", ylim=c(0.5,1.7))
abline(h=1, lty=2, lwd=2)
mtext(side=3, line=2, "Grain yield RYT", font=2, cex=size_main, at=3.5)
mtext(side=2, line=4, expression("RYT"["GY"]),cex=1.5)
mtext(side=3, line=4, "b", font=2, cex=size_letter, at=-1.5)
axis(3, as.character(multi.test$groups), at=c(1:6), tick=F, line=-3, font=2, cex.axis=size_axis)
text(c(1,5),c(1.085,1.13),"*",cex=3, font=2)


## C. Grain protein concentration RYT
data_reg <- perfP[,which(colnames(perfP)%in%c("RYT_GPC", "fun_grp"))]
by(data_reg, data_reg$fun_grp, FUN = function(x) {t.test(x[,1], mu=1)$p.value})
mod <- lm(RYT_GPC ~ fun_grp, data=data_reg)
multi.test <- LSD.test(mod,"fun_grp")
multi.test <- multi.test$groups[match(levels(data_reg$fun_grp),row.names(multi.test$groups)),]

boxplot(RYT_GPC~fun_grp, data=data_reg, ylab="", xlab="Functional association", cex.axis=size_axis, cex=2, cex.lab=size_lab, las=1, col="grey75", ylim=c(0.85,1.15))
abline(h=1, lty=2, lwd=2)
mtext(side=3, line=2, "Grain protein\nconcentration RYT", font=2, cex=size_main, at=3.5)
mtext(side=2, line=4, expression("RYT"["GPC"]),cex=1.5)
mtext(side=3, line=4, "c", font=2, cex=size_letter, at=-1.5)
axis(3, as.character(multi.test$groups), at=c(1:6), tick=F, line=-3, font=2, cex.axis=size_axis)
text(c(1),c(0.978),"*",cex=3, font=2)

dev.off()



#########################################
##  Table with means and variances of trait values per functional group (Supplementary table 5) ##
#########################################


## Merging trait data and clustering information
trait_clust <- cbind(traitM,classes)
colnames(trait_clust)[ncol(trait_clust)] <- "grp"
Mmeans <- aggregate(.~grp, data=trait_clust, FUN= function(x) {round(mean(x),2)})
Mvar <- aggregate(.~grp, data=trait_clust, FUN= function(x) {round(sd(x),2)} )


tot <- Mmeans
for (t in colnames(Mmeans)[2:ncol(Mmeans)]) {
  sub_trait_clust <- data.frame(trait=trait_clust [,t], grp=trait_clust $grp)
  mod <- lm(trait~grp, data=sub_trait_clust)
  multi.test <- LSD.test(mod,"grp")
  multi.test <- multi.test$groups[match(levels(trait_clust $grp),row.names(multi.test$groups)),]
  tot[,which(colnames(tot)==t)] <- multi.test[,2]
}

Mmeans <- t(Mmeans)
Mvar <- t(Mvar)
tot <- t(tot)

for (i in 2:nrow(Mmeans)) {
  for (j in 1:ncol(Mmeans)) {
    Mmeans[i,j] <- paste(Mmeans[i,j], " (", Mvar[i,j],")",tot[i,j], sep="")
  }
}
write.csv(Mmeans, file="functional_group_means_var.csv", row.names=T)
