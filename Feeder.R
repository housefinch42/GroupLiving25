#Copyright (c) 2024 [Marissa Mae Langager]

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

library(dplyr)
library(ordinal)
library(car)
library(emmeans)
library(lattice)
library(RVAideMemoire)
library(ggplot2)
library(survival)
library(survminer)
library(Rcpp)
library(coxme)

feed <- read.csv("Feeder_Behavior_22.csv")
View(feed)

str(feed$high.lvl.nov) #integer
feed[,'high.lvl.norm'] <- factor(feed[,'high.lvl.norm']) 

feed[,'high.lvl.nov'] <- factor(feed[,'high.lvl.nov'])
class(feed$high.lvl.nov)

feed[,'soc.trt'] <- factor(feed[,'soc.trt'])
class(feed$soc.trt)

feed[,'inf.feed'] <- factor(feed[,'inf.feed'])

#releveling predictor variables
feed <- mutate(feed,
               soc.trt=factor(soc.trt, levels = c("s", "fl")), 
               inf.feed=factor(inf.feed, levels = c("c", "i")))


#Fisher's Exact Test - normal feeder
fisher.test(feed$high.lvl.norm, feed$trt1)

#removing additional cagemates
#checking for possible effects of pseudoreplication
psr <-feed[-c(18,2,3,20,21,22,5,6,24,8,26,9,28,12,13,30,31,32,33,35,36,37,39,16,40),]
View(psr)

fisher.test(psr$high.lvl.norm, psr$trt1) #p=0.0039, still see sig diff without extra cagemates

#normal feeder CLMMs
#feeder1 <- clmm(high.lvl.norm ~ soc.trt+inf.feed + (1|cage.id), data = feed3)
##Anova.clmm(feeder1, type = 2) #jacked up, Hessian singularity
#don't do this. use chi-squared test on the differences in proportions of "successes"
#does random effect matter since there were no differences among any group other than the single-housed diseased birds?

#CLMM - novel feeder
nov1 <- clmm(high.lvl.nov ~ soc.trt*inf.feed + (1|cage.id), data = feed)
summary(nov1)
Anova.clmm(nov1, type = 3) #interaction not significant, remove from model

#subsetted data
novpsr <- clmm(high.lvl.nov ~ soc.trt*inf.feed + (1|cage.id), data = psr)
summary(novpsr)
Anova.clmm(novpsr, type = 3)

#taking out interaction, which isn't significant, effect social treatment p=0.0017
nov2 <- clmm(high.lvl.nov ~ soc.trt + inf.feed + (1|cage.id), data = feed)
summary(nov2)
Anova.clmm(nov2, type = 2)

#subsetted data
novpsr2 <- clmm(high.lvl.nov ~ soc.trt + inf.feed + (1|cage.id), data = psr)
summary(novpsr2)
Anova.clmm(novpsr2, type = 2)

###COX PROPORTIONAL HAZARDS###
survfeed <- read.csv("survfeed.csv")
View(survfeed)

#checking for possible effects of pseudoreplication
#subsetted data, removing additional cagemates
survpsr <- survfeed[-c(3,4,5,8,9,10,14,15,16,24,25,26,32,33,34,36,37,38,39,47,48,49,52,54,55),]
View(survpsr)

#normal feeder
normfeed1 <- coxme(Surv(time.norm, status.norm) ~ soc.trt * inf.feed + (1|cage.id), data = survfeed)
summary(normfeed1)

#subsetted data
normpsr1 <- coxme(Surv(time.norm, status.norm) ~ soc.trt*inf.feed + (1|cage.id), data = survpsr)
summary(normpsr1)

#interaction not significant, additive model
normfeed2 <- coxme(Surv(time.norm, status.norm) ~ soc.trt + inf.feed + (1|cage.id), data = survfeed)
summary(normfeed2)
Anova(normfeed2, type = 2)

#subsetted data
normpsr2 <- coxme(Surv(time.norm, status.norm) ~ soc.trt + inf.feed + (1|cage.id), data = survpsr)
summary(normpsr2)
Anova(normpsr2, type = 2)

#novel feeder
nov.cox1 <- coxme(Surv(time.nov, status.nov) ~ soc.trt * inf.feed + (1|cage.id), data = survfeed)
#error message, have tried to decode, but not sure what it means

#subsetted data
novpsr1 <- coxme(Surv(time.nov, status.nov) ~ soc.trt * inf.feed  + (1|cage.id), data = survpsr)
#well, same error message

#additive model works
nov.cox2 <- coxme(Surv(time.nov, status.nov) ~ soc.trt + inf.feed + (1|cage.id), data = survfeed)
summary(nov.cox2)
Anova(nov.cox2, type = 2)

#subsetted data
novpsr2 <- coxme(Surv(time.nov, status.nov) ~ soc.trt + inf.feed + (1|cage.id), data = survpsr)
summary(novpsr2)
Anova(novpsr2, type = 2)

####cox prop haz for prevalence###
prevfeed <- read.csv("survfeed_prevhighlow.csv")
View(prevfeed)

#subsetted data
prevgroup <- prevfeed[-c(41:56),]
View(prevgroup)

normprev1 <- coxme(Surv(time.norm, status.norm) ~ prev + (1|cage.id), data = prevgroup)
summary(normprev1)
Anova(normprev1, type = 3)

novprev1 <- coxme(Surv(time.nov, status.nov) ~ prev + (1|cage.id), data = prevgroup)
summary(novprev1)
Anova(novprev1, type = 3)

###FIGS###
#stacked bar plots#

feedprop <- read.csv("propfeed.csv")

str(feedprop$resp)
feedprop[,'resp'] <- factor(feedprop[,'resp']) #change resp to factor
class(feedprop$resp)

feedprop$trt.2[feedprop$trt=="fl_d"]="group-housed diseased"
feedprop$trt.2[feedprop$trt=="fl_u"]="group-housed healthy"
feedprop$trt.2[feedprop$trt=="s_d"]="single-housed diseased"
feedprop$trt.2[feedprop$trt=="s_u"]="single-housed healthy"

feedprop<- mutate(feedprop, 
                  trt.2=factor(trt.2,
                               levels = c("single-housed healthy", "single-housed diseased", "group-housed healthy", "group-housed diseased")))
#normal feeder#
normfeed <- ggplot(data = feedprop, aes(x=trt.2, y=prop.norm, fill=resp))+
  geom_bar(stat = "identity", color="black")+
  scale_fill_brewer(palette = "PuBu")+
  ylab("Proportion Birds")+
  xlab("Treatment")+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12.5),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(t=7, r=0, b=0, l=0)),
        axis.title.y = element_text(margin = margin(t=0, r=5, b=0, l=0)),
        axis.line=element_line(), 
        legend.position = "none")
normfeed

#novel feeder#
novfeed <- ggplot(data = feedprop, aes(x=trt.2, y=prop.nov, fill=resp))+
  geom_bar(stat = "identity", color="black")+
  scale_fill_brewer(palette = "PuBu")+
  guides(fill = guide_legend(title = "Novel Response"))+
  ylab("Proportion Birds")+
  xlab("Treatment")+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12.5),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(t=7, r=0, b=0, l=0)),
        axis.title.y = element_text(margin = margin(t=0, r=5, b=0, l=0)),
        axis.line=element_line(),
        legend.position = "none")
novfeed

###Decay Plots###
#Normal Feeder#
normfeed <- survfit(Surv(time.norm, status.norm) ~ soc.trt+inf.feed, data = survfeed)

feedplot1 <-ggsurvplot(normfeed, 
                       pval = FALSE, conf.int = TRUE,
                       linetype = "strata",
                       ggtheme = theme_bw(), 
                       palette = c("#63C5DA", "#0a1172", "#BC544B", "#420C09"))



feedplot1

#Novel Feeder#
nov.surv1 <- survfit(Surv(time.nov, status.nov) ~ soc.trt+inf.feed, data = survfeed)

feedplot2 <- ggsurvplot(nov.surv1, 
                        pval = FALSE, conf.int = TRUE,
                        linetype = "strata", 
                        ggtheme = theme_bw(),
                        palette = c("#63c5da", "#0a1172", "#bc544b", "#420c09"))
feedplot2

#prevalence figs#
prevfeed <- read.csv("survfeed_prevhighlow.csv")

str(prevfeed$prev)
prevfeed[,'prev'] <- factor(prevfeed[,'prev']) #change prev to factor
class(prevfeed$prev)

prevfeed$prev2[prevfeed$prev=="high"]="high disease prevalence group"
prevfeed$prev2[prevfeed$prev=="low"]="low disease prevalence group"
prevfeed$prev2[prevfeed$prev=="diseased"]="diseased single-housed"
prevfeed$prev2[prevfeed$prev=="healthy"]="healthy single-housed"

prevfeed<- mutate(prevfeed, 
                  prev2=factor(prev2,
                               levels = c("healthy single-housed", "diseased single-housed", "low disease prevalence group", "high disease prevalence group")))

normprev <- survfit(Surv(time.norm, status.norm) ~ prev2, data = prevfeed)

normprev2 <- survfit(Surv(time.norm, status.norm) ~ cage.id, data = prevfeed)

prevplot1 <-ggsurvplot(normprev, 
                       pval = FALSE, conf.int = TRUE,
                       linetype = "strata",
                       ggtheme = theme_bw(), 
                       palette = c("#63c5da", "#0a1172", "#bc544b", "#420c09"), 
                       xlim = c(0,2500))

prevplot1


prevplot2 <- ggsurvplot(normprev2, 
                        pval = FALSE, conf.int = FALSE, 
                        linetype = "strata", 
                        ggtheme = theme_bw(),
                        palette = c("#ff0000","#ffa500", "#ffff00", "#008000", "#0000ff", "#4b0082", "#ee82ee", "#000000ff" ), 
                        xlim = c (0, 2500))
prevplot2

novprev <- survfit(Surv(time.nov, status.nov) ~ prev2, data = prevfeed)

novprev2 <- survfit(Surv(time.nov, status.nov) ~ cage.id, data = prevfeed)

novplot1 <- ggsurvplot(novprev, 
                       pval = FALSE, conf.int = TRUE, 
                       linetype = "strata", 
                       ggtheme = theme_bw(), 
                       palette = c("#63c5da", "#0a1172", "#bc544b", "#420c09"))
novplot1

novplot2 <- ggsurvplot(novprev2, 
                        pval = FALSE, conf.int = FALSE, 
                        linetype = "strata", 
                        ggtheme = theme_bw(),
                        palette = c("#ff0000","#ffa500", "#ffff00", "#008000", "#0000ff", "#4b0082", "#ee82ee", "#000000ff" ), 
                        xlim = c (0, 3500))
novplot2

############# Who eats first? ###############

flockfeed <- survfeed %>% filter(soc.trt == "fl")
View(flockfeed)

flockfeed <- flockfeed %>%
  group_by(cage.id) %>%
  mutate(norm.rank = rank(time.norm),
         nov.rank = rank(time.nov),
         norm.eatFirst = ifelse(norm.rank == 1, "y", "n"),
         nov.eatFirst = ifelse(nov.rank == 1, "y", "n"))  # use -score for descending order

head(flockfeed$norm.rank)

ggplot(data = flockfeed, aes(x = inf.feed, y = norm.rank))+
  geom_jitter(width = 0.1)

feedinf <- flockfeed %>% filter(inf.feed == "i")
feedc <- flockfeed %>% filter(inf.feed == "c")

inf.norm <- feedinf$norm.rank

c.norm <- feedc$norm.rank

#Mann Whitney U test of ranks
wilcox.test(inf.norm, c.norm, paired = FALSE) #No difference in ranks

#### Novel feeder
ggplot(data = flockfeed, aes(x = nov.rank))+
  geom_bar(aes(fill = inf.feed), position = "dodge")

inf.nov <- feedinf$nov.rank

c.nov <- feedc$nov.rank

#Mann Whitney U test of ranks
wilcox.test(inf.nov, c.nov, paired = FALSE) #Significant difference in ranks

