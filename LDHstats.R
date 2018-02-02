

# Load packages

library(readr)
library(car)
library(ggplot2)
library(dplyr)
library(lsmeans)
library(multcompView)

# ggplot theme
theme_bw()

# Load CSV file

library(readr)
LDHdata <- read_csv("~/Documents/Davis/Burbot Experiment/BURBOT-PROTEIN-ASSAYS/2018-jan-LDH-summary.csv", 
                    col_types = cols(Family = col_factor(levels = 
                                                           c("Ma", "Mb", "Mc", "Md", "RV")), 
                                     `Feeding Strategy` = col_factor(levels = c("Cannibal", "Planktivore"))))
View(LDHdata)

# Briefly check data for weirdness:
str(LDHdata) #check each column's data type
summary(LDHdata) #easy way to check that data is entered correctly and that distributions make sense
head(LDHdata)


# Check spread of LDH values without looking at treatment. Everything looks reasonable. :)
p=ggplot(data=LDHdata, aes(x=Sample, y=`Cold LDH activity`))
p+geom_point()

p=ggplot(data=LDHdata, aes(x=Sample, y=`Warm LDH activity`))
p+geom_point()


# Calculate Q10 while data is in this format

## T1=14°C T2=24°C 
LDHdataQ10 <- LDHdata %>%
  mutate(Q10=((`Warm LDH activity`/`Cold LDH activity`)^(10/(24-14)))) %>% 
  rename(LDH.cold=`Cold LDH activity`, LDH.warm=`Warm LDH activity`, 
                       feedstrat=`Feeding Strategy`)

##Look at Q10 between cannibals and planktivores to see if they have different temperature sensititivty. 

pairs(LDHdataQ10)

pd <-position_dodge(0.8)
p=ggplot(data=LDHdataQ10, aes(x=Family, y=Q10))
p+geom_boxplot(aes(fill=feedstrat), position=pd) + 
  labs(x = "Family", y="Q10") + geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5), linetype="solid", color="darkgrey") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.title = element_blank(), 
    axis.text = element_text(size=12, color="black"), 
    title = element_text(size=14),
    legend.text = element_text(size=12), 
    legend.box.background = element_rect(color="black"), 
    legend.position = c(0.9, 0.91))

#Stats on LDH Activity at enviornmental temp (14°C)

##Boxplot of LDH activity at 14°C
pd <-position_dodge(0.8)
p=ggplot(data=LDHdataQ10, aes(x=Family, y=LDH.cold))
p+geom_boxplot(aes(fill=feedstrat), position=pd) + 
  labs(x = "Family", y="LDH Activity (U/mg protein)") + geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5), linetype="solid", color="darkgrey") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.title = element_blank(), 
    axis.text = element_text(size=12, color="black"), 
    title = element_text(size=14),
    legend.text = element_text(size=12), 
    legend.box.background = element_rect(color="black"), 
    legend.position = c(0.9, 0.91))
  
##Linear model
LDHmod = lm(LDH.cold ~ feedstrat*Family, data=LDHdataQ10)
plot(LDHmod) #Residuals look fine. Data is normally distributed. 
summary(LDHmod)
Anova(LDHmod) #Both feeding strategy and family explain variation in the data! And their interaction borderline explains variation in the data.
anova(LDHmod)

##Test variance within feeding strategies and families 
op <- par(mfrow = c(2,2))
plot(LDHmod, which=1)
residuals <- resid(LDHmod)
hist(residuals, xlab="Residuals", main="")
plot(LDHdataQ10$feedstrat, residuals, xlab = "Feeding Strategy", ylab = "Residuals")
plot(LDHdataQ10$Family, residuals, xlab = "Family", ylab="Residuals")
par(op) 
###Heterogeneity by family?

## 3 Models to include heterogeneity 
library(nlme)
f1 <- formula(LDH.cold ~ feedstrat*Family)
LDHmod.gls0 <- gls(f1, data=LDHdataQ10)
LDHmod.gls1 <- gls(LDH.cold ~ feedstrat*Family, data=LDHdataQ10, 
                   weights = varIdent(form=~1 | feedstrat))
LDHmod.gls2 <- gls(LDH.cold ~ feedstrat*Family, data=LDHdataQ10, 
                   weights = varIdent(form=~1 | Family))
LDHmod.gls3 <- gls(LDH.cold ~ feedstrat*Family, data=LDHdataQ10, 
                   weights = varIdent(form=~1 | feedstrat*Family))
anova(LDHmod.gls0, LDHmod.gls1, LDHmod.gls2, LDHmod.gls3)
#AIC values for gls models are larger --> No evidence that using a gls model is necessary. 

##Linear model with family nested in feeding strategy
LDHmod.nested = lm(LDH.cold ~ feedstrat/Family, data=LDHdataQ10)
plot(LDHmod.nested)
Anova(LDHmod.nested)
summary(LDHmod.nested)
#Actually, probs don't want a nested model because we're also interested in seeing if each family follows the same pattern


#_ _ _ _ _ _ _ _
# 1/26/18
# Mandy Frazier

## Tukey Post Hoc Test to see if the means are significantly different from each other. 

TukeyTest = lsmeans(LDHmod, tukey ~ feedstrat:Family)
cld(TukeyTest)

TukeyTest = lsmeans(LDHmod, tukey ~ feedstrat)
cld(TukeyTest)

TukeyTest = lsmeans(LDHmod, tukey ~ Family)
cld(TukeyTest)


## In general, cannibals tend to have higher LDH activity, but there is a strong familial effect. 

### If you look at the mass boxplot and the LDH boxplot, they look VERY similar. Perhaps LDH is simply associated with mass. 

LDHdataQ10mass <- LDHdataQ10 %>%
  mutate(mass=burbotmasterdata$`Whole Fish Mass (g)`)

ggplot(LDHdataQ10mass, aes(x=mass, y=LDH.cold, group= feedstrat, color=feedstrat)) +
  geom_smooth(method='lm', se=F) + geom_point() +
  labs(x = "Mass (g)", y = "LDH Activity (U / mg protein)", color= "Feeding Strategy") +
  theme_bw() +
  theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  axis.text = element_text(size=14), 
  title = element_text(size=14),
  strip.text.x = element_text(size=14), 
  legend.text = element_text(size=12), 
  legend.position = c(0.8, 0.15), 
  legend.box.background = element_rect(color="black"))
 

ggplot(LDHdataQ10mass, aes(x=mass, y=LDH.cold, group= Family, color=Family)) +
  geom_smooth(method='lm', se=F) + geom_point() +
  labs(x = "Mass (g)", y = "LDH Activity (U / mg protein)", color= "Feeding Strategy") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size=14), 
    title = element_text(size=14),
    strip.text.x = element_text(size=14), 
    legend.text = element_text(size=12), 
    legend.position = c(0.8, 0.15), 
    legend.box.background = element_rect(color="black"))





massmodel = lm(LDH.cold~mass*feedstrat, LDHdataQ10mass)
summary(massmodel)

##Strong evidence that mass influences LDH 

### Higher LDH levels indicate that that individual has a higher capacity for anaerobic capacity. 

# - - - - - 
# 1/31/18
# Mandy Frazier

# Mass has an effect on LDH activity, so it should probably be incorporated as a part of the linear model for LDH. Let's try: 

LDHmod.mass = lm(LDH.cold ~ feedstrat*Family*mass, data=LDHdataQ10mass)
plot(LDHmod.mass) #What's going on with the "Residuals vs Leverage" plot?? Something is wacky. 
summary(LDHmod.mass)
Anova(LDHmod.mass) #If this is correct, it appears that the interaction between feeding strategy, family, and mass are the only thing that significantly explains LDH activity in this model. How do I know which model is better?
#We don't really want to do a 3-way ANOVA (which this is doing). We want to do an ANCOVA because mass is a co-variate. 

#ANCOVA

LDHmod = lm(LDH.cold ~ feedstrat*Family, data=LDHdataQ10)
LDHmodmass = lm(LDH.cold ~ mass*Family, data=LDHdataQ10mass)



LDHmod.mass.ancova = lm(LDH.cold ~ mass*feedstrat+feedstrat*Family, data=LDHdataQ10mass)
plot(LDHmod.mass.ancova)

LDHmod.mass.ancova2 = lm(LDH.cold ~ mass+feedstrat*Family, data=LDHdataQ10mass)
plot(LDHmod.mass.ancova2)

LDHmod.mass.ancova3 = lm(LDH.cold ~ mass + mass:feedstrat + feedstrat*Family, data=LDHdataQ10mass)
plot(LDHmod.mass.ancova3)

LDHmod.mass = lm(LDH.cold ~ mass, data=LDHdataQ10mass)


anova(LDHmod.mass.ancova2, LDHmod.mass.ancova3)


#Test variance within feeding strategies and families:
op <- par(mfrow = c(2,2))
plot(LDHmod.mass.ancova, which=1)
residuals <- resid(LDHmod.mass.ancova)
hist(residuals, xlab="Residuals", main="")
plot(LDHdataQ10mass$feedstrat, residuals, xlab = "Feeding Strategy", ylab = "Residuals")
plot(LDHdataQ10mass$Family, residuals, xlab = "Family", ylab="Residuals")
par(op) 

library(nlme)
f1.mass <- formula(LDH.cold ~ feedstrat*Family)
LDHmod.gls0.mass <- gls(f1, data=LDHdataQ10mass)
LDHmod.gls1.mass <- gls(LDH.cold ~ feedstrat*Family*mass, data=LDHdataQ10mass, 
                   weights = varIdent(form=~1 | feedstrat))
LDHmod.gls2.mass <- gls(LDH.cold ~ feedstrat*Family*mass, data=LDHdataQ10mass, 
                   weights = varIdent(form=~1 | Family))
LDHmod.gls3.mass <- gls(LDH.cold ~ feedstrat*Family*mass, data=LDHdataQ10mass, 
                   weights = varIdent(form=~1 | feedstrat*Family))
library(nlme)
f1.mass <- formula(LDH.cold ~ feedstrat*Family*mass)
LDHmod.gls0.mass <- gls(f1, data=LDHdataQ10)
LDHmod.gls1.mass <- gls(LDH.cold ~ feedstrat*Family*mass, data=LDHdataQ10mass, 
                   weights = varIdent(form=~1 | feedstrat))
LDHmod.gls2.mass <- gls(LDH.cold ~ feedstrat*Family*mass, data=LDHdataQ10mass, 
                   weights = varIdent(form=~1 | Family))
LDHmod.gls3.mass <- gls(LDH.cold ~ feedstrat*Family*mass, data=LDHdataQ10mass, 
                   weights = varIdent(form=~1 | feedstrat*Family))
LDHmod.gls4.mass <- gls(LDH.cold ~ feedstrat*Family*mass, data=LDHdataQ10mass, 
                   weights = varIdent(form=~1 | feedstrat*mass))
LDHmod.gls5.mass <- gls(LDH.cold ~ feedstrat*Family*mass, data=LDHdataQ10mass, 
                   weights = varIdent(form=~1 | Family*mass))
LDHmod.gls6.mass <- gls(LDH.cold ~ feedstrat*Family*mass, data=LDHdataQ10mass, 
                   weights = varIdent(form=~1 | feedstrat*Family*mass))
anova(LDHmod.gls0.mass, LDHmod.gls1.mass, LDHmod.gls2.mass, LDHmod.gls3.mass, LDHmod.gls4.mass, LDHmod.gls5.mass, LDHmod.gls6.mass)

#The main issue is that we have two variables (feeding strategy and mass) that are colinear and we don't know which variable is driving the other. We can't include both mass and feeding strategy in the model because of their colinearity. The model that includes feeding strategy and family suits the data, but so does the model that includes mass and family. So we're not sure if model or feeding strategy is driving the LDH activity. 
 

