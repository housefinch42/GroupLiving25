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


#CLMMS-quartiles
behav2 <- read.csv("Pred_Behav.csv") #load predator dataset 
behav2$tot.hawk <- behav2$sec.hawk.end-behav2$sec.hawk.start #adding column to denote total time hawk was "flying"

behav2$tot.assay.time <- (behav2$pred.end - behav2$pred.start) #total time of assay (start of predator playback to end of recording)

#dividing response to predator flyover into reaction time quartiles
behav2$quartile <- ntile(behav2$sec.to.resp, 4)
View(behav2)

behav2[,'quartile'] <- factor(behav2[,'quartile']) #making quartile a factor
class(behav2$quartile)

#original models using response times split into 6 separate bins 
#(too many 0s, causing error in clmm)
#pred1 <- clmm (resp.6bin ~ soc.trt + inf.pred + soc.trt*inf.pred + (1|cage.id), data = behav2) 
#getting error message, too many zeroes among response factors
#model using quartiles instead of response - much less likelihood of excess zeros in each bin
#interactive model - quartile ~ social group treatment*disease status during predator assay

behav2[,'inf.pred'] <- factor(behav2[,'inf.pred'])
behav2[,'soc.trt'] <- factor(behav2[,'soc.trt'])

str(behav2$inf.pred)
str(behav2$soc.trt)

#reordering levels from alphabetical
behav2 <- mutate(behav2,
                 soc.trt=factor(soc.trt, levels = c("s", "fl")), 
                 inf.pred=factor(inf.pred, levels = c("c", "i")))
levels(behav2$soc.trt)
levels(behav2$inf.pred)

#model with quartile as response variable
pred1 <- clmm(quartile ~ soc.trt*inf.pred + (1|cage.id), data = behav2) 
summary(pred1)
Anova.clmm(pred1, type = 3) #type 3 where interaction present

#####IMMOBILITY##### 
behav2$prop.imm <- (behav2$sec.imm.pred/behav2$tot.assay.time)

imod1 <- lmer(prop.imm ~ soc.trt*inf.pred + (1|cage.id), 
              data = behav2, weights = tot.assay.time)
summary(imod1)
Anova(imod1, type = 3) 

emmeans(imod1, pairwise ~ soc.trt*inf.pred)

######COX PROPORTIONAL HAZARDS######
pred <- read.csv("survpred.csv")
View(pred)

predcox1 <- coxme(Surv(time, status) ~ soc.trt * inf.pred + (1|cage.id), data = pred)
summary(predcox1)
Anova(predcox1, type = 3) 
#interaction not significant, remove from model

predcox2 <- coxme(Surv(time, status) ~ soc.trt + inf.pred + (1|cage.id), data = pred)
summary(predcox2)
Anova(predcox2, type = 2)#USE

#FIGS
###Stacked bar plot###
propbehav <- read.csv("proppred.csv")

propbehav$prop.quart <- propbehav$num.pred.quart/propbehav$total.bird

str(propbehav$quart)
propbehav[,'quart'] <- factor(propbehav[,'quart']) #change resp to factor
class(propbehav$quart)

#renaming treatments
propbehav$trt.2[propbehav$trt=="fl_d"]="group-housed diseased"
propbehav$trt.2[propbehav$trt=="fl_u"]="group-housed healthy"
propbehav$trt.2[propbehav$trt=="s_d"]="single-housed diseased"
propbehav$trt.2[propbehav$trt=="s_u"]="single-housed healthy"

#reordering levels from alphabetical
propbehav <- mutate(propbehav, 
                    trt.2=factor(trt.2,
                                 levels = c("single-housed healthy", "single-housed diseased", "group-housed healthy", "group-housed diseased")))

predplot <- ggplot(data = propbehav, aes(x=trt.2, y=prop.quart, fill=quart))+
  geom_bar(stat = "identity", color="black")+
  scale_fill_brewer(palette = "BuPu")+
  guides(fill=guide_legend(title = "Pred. Resp."))+
  ylab("Proportion Birds")+
  xlab("Treatment")+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 12.5),
        panel.grid = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t=7,r=0, b=0, l=0)),
        axis.title.y = element_text(margin = margin(t=0,r=5,b=0,l=0)),
        axis.line = element_line(),
        legend.box.background = element_rect(colour = "black"))


predplot

###Immobility Fig###
behav2$comb.trt <- paste0(behav2$soc.trt, ".", behav2$inf.pred) #new column with combined treatments
View(behav2)

behav2$trt.2[behav2$comb.trt=="fl.i"]="group-housed diseased"
behav2$trt.2[behav2$comb.trt=="fl.c"]="group-housed healthy"
behav2$trt.2[behav2$comb.trt=="s.i"]="single-housed diseased"
behav2$trt.2[behav2$comb.trt=="s.c"]="single-housed healthy"

behav2[,'trt.2'] <- factor(behav2[,'trt.2'])

behav2 <- mutate(behav2, 
                 trt.2=factor(trt.2,
                              levels = c("single-housed healthy", "single-housed diseased", "group-housed healthy", "group-housed diseased")))

immfigall1 <- ggplot(data = behav2, aes(x=trt.2, y=prop.imm, color=comb.trt))+
  geom_boxplot(outlier.shape = NA)+
  theme_light()+
  scale_color_manual(values = c("s.c"="#BC544B","s.i"="#420C09","fl.c"="#63C5DA","fl.i"="#0a1172"))+
  scale_y_continuous(limits = c(0,1))+
  geom_jitter(position = position_jitterdodge(jitter.width=0.6, jitter.height = 0, dodge.width=0.75))+
  guides(fill = guide_legend(title = "Social Treatment"))+
  ylab("Proportion Time Spent Immobile")+
  xlab("Treatment")

immfigall1

###survival plot###
fit2 <- survfit(Surv(time,status) ~ soc.trt + inf.pred, data = pred)

predsurv <- ggsurvplot(fit2,
                       pval = FALSE, conf.int = TRUE,
                       linetype = "strata",
                       ggtheme = theme_bw(), 
                       palette = c("#63c5da", "#0a1172", "#bc544b", "#420c09"))

predsurv
