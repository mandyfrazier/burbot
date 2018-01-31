#### Plotting LDH activity ####

## Load packages

library(readr)
library(car)
library(ggplot2)
library(dplyr)
library(lsmeans)
library(multcompView)

## ggplot theme
theme_bw()

## Load CSV file

library(readr)
LDHdata <- read_csv("~/Documents/Davis/Burbot Experiment/BURBOT-PROTEIN-ASSAYS/2018-jan-LDH-summary.csv", 
                    col_types = cols(Family = col_factor(levels = 
                                                           c("Ma", "Mb", "Mc", "Md", "RV")), 
                                     `Feeding Strategy` = col_factor(levels = c("Cannibal", "Planktivore"))))
View(LDHdata)

## Briefly check data for weirdness:
str(LDHdata) ##check each column's data type
summary(LDHdata) ##easy way to check that data is entered correctly and that distributions make sense
head(LDHdata)


## Check spread of LDH values without looking at treatment. Everything looks reasonable. :)
p=ggplot(data=LDHdata, aes(x=Sample, y=`Cold LDH activity`))
p+geom_point()

p=ggplot(data=LDHdata, aes(x=Sample, y=`Warm LDH activity`))
p+geom_point()


## Calculate Q10 while data is in this format

# T1=14°C T2=24°C 
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

##Test variance within feeding strategies and families 
op <- par(mfrow = c(2,2))
plot(LDHmod, which=1)
residuals <- resid(LDHmod)
hist(residuals, xlab="Residuals", main="")
plot(LDHdataQ10$feedstrat, residuals, xlab = "Feeding Strategy", ylab = "Residuals")
plot(LDHdataQ10$Family, residuals, xlab = "Family", ylab="Residuals")
par(op) 
###Heterogeneity by family?

#### 3 Models to include heterogeneity 
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


_ _ _ _ _ _ _ _
## 1/26/18
## Mandy Frazier

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
 
massmodel = lm(LDH.cold~mass*feedstrat, LDHdataQ10mass)
summary(massmodel)

##Strong evidence that mass influences LDH 

### Higher LDH levels indicate that that individual has a higher capacity for anaerobic capacity. 








