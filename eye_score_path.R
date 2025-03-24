library(ggplot2)
library(car)
library(dplyr)
library(effects)
library(scales)
library(lme4)


#full dataset
flock_eye <- read.csv("eye_score_all.csv")
View(flock_eye)


#converting full dataset into proper classes
str(flock_eye)

flock_eye$bird.id <- as.factor(flock_eye$bird.id)
flock_eye$flock.trt <- as.factor(flock_eye$flock.trt)
flock_eye$assign.trt <- as.factor(flock_eye$assign.trt)
flock_eye$mg.status <- as.factor(flock_eye$mg.status)
flock_eye$PID.code <- as.factor(flock_eye$PID.code)

flock_eye <- mutate(flock_eye, 
                 flock.trt=factor(flock.trt,
                              levels = c("s", "fl")))
levels(flock_eye$flock.trt)

#Nice Labels
flock_eye$PID2[flock_eye$PID.code=="a"]="Day -4"
flock_eye$PID2[flock_eye$PID.code=="b"]="Day 7/8"
flock_eye$PID2[flock_eye$PID.code=="c"]="Day 15/16"

View(flock_eye)

str(flock_eye$PID2)
flock_eye$PID2 <- as.factor(flock_eye$PID2)


#models
#dataset with only birds inoculated with MG
infonly <-flock_eye[which(flock_eye$assign.trt=="i"),]
View(infonly)

infonly2 <- infonly[-which(infonly$PID=="-4"),] #dataset removing pre-inoculation eye score (reduce 0s)


infonly4 <- infonly[which(infonly$PID2=="Day 7/8"),]
mean(infonly4$comb.score)
sd(infonly4$comb.score)

infonly3 <-infonly[which(infonly$PID2=="Day 15/16"),]
View(infonly3)
mean(infonly3$comb.score)
sd(infonly3$comb.score)

mod2 <- lmer(comb.score ~ flock.trt + (1|bird.id), data = infonly2)
summary(mod2)
Anova(mod2, type = 2)


plot(mod2)
plot(allEffects(mod2))
resmod2 <- resid(mod2)
hist(resmod2)


eyeplot <- ggplot(data = infonly2, aes(x=PID2, y=comb.score, color=flock.trt))+
  scale_color_manual(values = c("#BC544B", "#63C5DA"))+
  geom_boxplot(outlier.shape = NA)+
  ylab("Eye Score")+
  xlab("Post-inoculation Day")+
  scale_x_discrete(limits=c("Day 7/8","Day 15/16"))+
  geom_jitter(position = position_jitterdodge(jitter.width=0.6, jitter.height = 0, dodge.width=0.75))

eyeplot

#path
flock_path <- read.csv("path_load_all.csv")
View(flock_path)

str(flock_path)
flock_path$bird.id <- as.factor(flock_path$bird.id)
flock_path$flock.trt <- as.factor(flock_path$flock.trt)
flock_path$assign.trt <- as.factor(flock_path$assign.trt)
flock_path$mg.status.path <- as.factor(flock_path$mg.status.path)
flock_path$PID.code <- as.factor(flock_path$PID.code)

flock_path <- mutate(flock_path, 
                    flock.trt=factor(flock.trt,
                                     levels = c("s", "fl")))

infonlypath <-flock_path[which(flock_path$assign.trt=="i"),]
View(infonlypath)

pathmod1 <- lmer(path.load ~ flock.trt + (1|bird.id), data = infonlypath)
summary(pathmod1)
Anova(pathmod1, type = 2)

plot(allEffects(pathmod1))
plot(pathmod1)

infonlypath$PID2[infonlypath$PID.code=="a"]="Day 7/8"

infonlypath$PID2[infonlypath$PID.code=="b"]="Day 15/16"

pathplot1 <- ggplot(data = infonlypath, aes(x=PID2, y=path.load, color=flock.trt))+
  scale_colour_manual(values = c("fl" = "#63C5DA", "s" = "#BC544B"), labels = c("fl" = "flock-housed", "s" = "single-housed"))+
  geom_boxplot(outlier.shape = NA)+
  ylab("Pathogen Load")+
  xlab("Post-inoculation Day")+
  scale_x_discrete(limits=c("Day 7/8", "Day 15/16"))+
  geom_jitter(position = position_jitterdodge(jitter.width=0.6, jitter.height = 0, dodge.width=0.75))

pathplot1

ggarrange(eyeplot, pathplot1)



