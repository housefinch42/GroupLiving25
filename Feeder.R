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
levels(behav2$soc.trt)
levels(behav2$inf.pred)

#Fisher's Exact Test - normal feeder
fisher.test(feed$high.lvl.norm, feed$trt1)

#normal feeder CLMMs
#feeder1 <- clmm(high.lvl.norm ~ soc.trt+inf.feed + (1|cage.id), data = feed3)
##Anova.clmm(feeder1, type = 2) #jacked up, Hessian singularity
#don't do this. use chi-squared test on the differences in proportions of "successes"
#does random effect matter since there were no differences among any group other than the single-housed diseased birds?

#CLMM - novel feeder
nov1 <- clmm(high.lvl.nov ~ soc.trt*inf.feed + (1|cage.id), data = feed)
summary(nov1)
Anova.clmm(nov1, type = 3) #interaction not significant, remove from model

#taking out interaction, which isn't significant, effect social treatment p=0.0017
nov2 <- clmm(high.lvl.nov ~ soc.trt + inf.feed + (1|cage.id), data = feed)
summary(nov2)
Anova.clmm(nov2, type = 2)



###COX PROPORTIONAL HAZARDS###
survfeed <- read.csv("survfeed.csv")

#normal feeder
normfeed1 <- coxme(Surv(time.norm, status.norm) ~ soc.trt * inf.feed + (1|cage.id), data = survfeed)
summary(normfeed1)

#interaction not significant, additive model
normfeed2 <- coxme(Surv(time.norm, status.norm) ~ soc.trt + inf.feed + (1|cage.id), data = survfeed)
summary(normfeed2)
Anova(normfeed2, type = 2)

#novel feeder
nov.cox1 <- coxme(Surv(time.nov, status.nov) ~ soc.trt * inf.feed + (1|cage.id), data = survfeed)
#error message, have tried to decode, but not sure what it means

#additive model works
nov.cox2 <- coxme(Surv(time.nov, status.nov) ~ soc.trt + inf.feed + (1|cage.id), data = survfeed)
summary(nov.cox2)
Anova(nov.cox2, type = 2)

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
