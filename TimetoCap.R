library(ggplot2)
library(dplyr)
library(car)
library(effects)
library(lme4)
library(emmeans)


time7 <- read.csv("d7ttc.csv")

View(time7)

time7 %>% count(inf.cap, sort = TRUE)

time7[,'inf.cap'] <- factor(time7[,'inf.cap'])
time7[,'soc.trt'] <- factor(time7[,'soc.trt'])

time7 <- mutate(time7,
                 soc.trt=factor(soc.trt, levels = c("s", "fl")), 
                 inf.cap=factor(inf.cap, levels = c("c", "i")))
levels(time7$soc.trt)
levels(time7$inf.cap)

#d7/8 models - gamma distribution - capture time ~ social treatment*disease status 
#random effect of cage id
d7glmm1 <- glmer(time ~ soc.trt*inf.cap+(1|cage.id), data = time7,
                 family = Gamma(link = "log"))

summary(d7glmm1)
Anova(d7glmm1, type = 3) #type 3 anova, significant interaction

plot(allEffects(d7glmm1))


#d15/16
time15 <- read.csv("d15ttc.csv")
View(time15)
time15 %>% count(comb.trt, sort = TRUE)

time15[,'inf.cap'] <- factor(time15[,'inf.cap'])
time15[,'soc.trt'] <- factor(time15[,'soc.trt'])

levels(time15$soc.trt)
levels(time15$inf.cap)

time15 <- mutate(time15,
                 soc.trt=factor(soc.trt, levels = c("s", "fl")), 
                 inf.cap=factor(inf.cap, levels = c("c", "i")))

levels(time15$soc.trt)
levels(time15$inf.cap)

#d15/16 interactive model - gamma distribution - capture time ~ social behavior*disease status
#random effect cage id
d15glmm1 <- glmer(time ~ soc.trt*inf.cap+(1|cage.id), data = time15,
                 family = Gamma(link = "log"))
summary(d15glmm1)
Anova(d15glmm1, type = 3) #interaction not significant, take out of model

#d15/16 additive model - gamma distribution - capture time ~ social behavior+disease status
#random effect cage id
d15glmm2 <- glmer(time ~ soc.trt + inf.cap + (1|cage.id), data = time15, 
                  family = Gamma(link = "log"))
summary(d15glmm2)
Anova(d15glmm2, type = 2)
plot(allEffects(d15glmm2))

#plotting
time7$trt.2[time7$comb.trt=="fl.i"]="group-housed diseased"
time7$trt.2[time7$comb.trt=="fl.c"]="group-housed healthy"
time7$trt.2[time7$comb.trt=="s.i"]="single-housed diseased"
time7$trt.2[time7$comb.trt=="s.c"]="single-housed healthy"

time7 <- mutate(time7, 
                    trt.2=factor(trt.2,
                                 levels = c("single-housed healthy", "single-housed diseased", "group-housed healthy", "group-housed diseased")))

ttc7 <- ggplot(data = time7, aes(x=trt.2, y=time, color = comb.trt))+
  geom_boxplot()+
  scale_colour_manual(values = c("fl.i" = "#0A1172","fl.c" = "#63C5DA", "s.i" = "#420C09", "s.c" = "#BC544B"))+
  ylab("Time to Capture (s)")+
  xlab("Treatment")+
  geom_jitter()+
  theme_bw()+
  theme(legend.position = "none", 
        text = element_text(colour = "black"), 
        axis.text.y = element_text(colour = "black"), 
        axis.text.x = element_text(colour = "black"))+
  ylim(0,40)
  

ttc7 #32 points, wee

#supplement w/all infected separate
time15$trt.2[time15$comb.trt=="fl.i"]="group-housed diseased"
time15$trt.2[time15$comb.trt=="fl.c"]="group-housed healthy"
time15$trt.2[time15$comb.trt=="s.i"]="single-housed diseased"
time15$trt.2[time15$comb.trt=="s.c"]="single-housed healthy"

time15 <- mutate(time15, 
                trt.2=factor(trt.2,
                             levels = c("single-housed healthy", "single-housed diseased", "group-housed healthy", "group-housed diseased")))

ttc2 <- ggplot(data = time15, aes(x=trt.2, y=time, color=comb.trt))+
  geom_boxplot()+
  scale_colour_manual(values = c("fl.i" = "#0A1172","fl.c" = "#63C5DA", "s.i" = "#420C09", "s.c" = "#BC544B"))+
  ylab("Time to Capture (s)")+
  xlab("Treatment")+
  geom_jitter()+
  theme_bw()+
  theme(legend.position = "none", 
        text = element_text(colour = "black"), 
        axis.text.y = element_text(colour = "black"), 
        axis.text.x = element_text(colour = "black"))+
  ylim(0,40)

ttc2 #it looks like the birds without any non-diseased birds in their flock were caught faster - very interesting! 

install.packages("ggpubr")
library(ggpubr)
 
ggarrange(ttc7, ttc2)


#dataframe with only flock infected birds #summary stats
flinf7 <-time7[which(time7$comb.trt=="fl.i"),]
View(flinf7)
range(flinf7$time)

#dataframe with only flock healthy birds #summary stats
flnot7 <-time7[which(time7$comb.trt=="fl.c"),]
range(flnot7$time)

#dataframe with only single healthy birds #summary stats
snot7 <-time7[which(time7$comb.trt=="s.c"),]
range(snot7$time)

#dataframe with only single infected birds #summary stats
sinf7 <-time7[which(time7$comb.trt=="s.i"),]
range(sinf7$time)

#dataframe with only flock infected #summary stats
flinf15 <-time15[which(time15$comb.trt=="fl.i"),]
range(flinf15$time)

#dataframe with only flock healthy #summary stats
flnot15 <-time15[which(time15$comb.trt=="fl.c"),]
range(flnot15$time)

#dataframe with only single healthy #summary stats
snot15 <-time15[which(time15$comb.trt=="s.c"),]
range(snot15$time)

#dataframe with only single infected #summary stats
sinf15 <-time15[which(time15$comb.trt=="s.i"),]
range(sinf15$time)