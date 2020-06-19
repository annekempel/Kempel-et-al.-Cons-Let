##17/16/2020
##Anne Kempel
##Code for Kempel et al., Nationwide re-visitation reveals thousands of local extinctions across the ranges of
## 713 threatened and rare plant species

#Dataframe needed: Data_Kempel et al.txt
#(Listened to Pink Floyd - The Wall)


library(multcomp)
library(data.table)
library(base)
library(plyr)  
library(ggplot2)

#load dataset
redlist<-read.table("Data_Kempel et al.txt", check.names = FALSE)
str(redlist)

#Prepair dataframe
#Log-transform Max_height
redlist$Maxheight_log<-log(redlist$Maxheight_cm+1)

#Scale variables and convert to factors
#scale(x, center = TRUE, scale = TRUE)
redlist$logyeardif_s<-scale(log(redlist$yeardif), center=TRUE, scale =TRUE)
redlist$C_s<-scale(as.numeric(redlist$C), center=TRUE, scale =TRUE) 
redlist$R_s<-scale(as.numeric(redlist$R), center=TRUE, scale =TRUE)
redlist$Kn_s<-scale(redlist$Kn, center=TRUE, scale =TRUE)
redlist$Nn_s<-scale(redlist$Nn, center=TRUE, scale =TRUE)
redlist$Fn_s<-scale(redlist$Fn, center=TRUE, scale =TRUE)
redlist$Tn_s<-scale(redlist$Tn, center=TRUE, scale =TRUE)
redlist$Ln_s<-scale(redlist$Ln, center=TRUE, scale =TRUE)
redlist$Rn_s<-scale(redlist$Rn, center=TRUE, scale =TRUE)
redlist$Maxheight_s<-scale(redlist$Maxheight_log, center=TRUE, scale =TRUE)
redlist$C_s<-as.factor(redlist$C_s)
redlist$R_s<-as.factor(redlist$R_s)
redlist$yeardif_s<-scale(redlist$yeardif, center=TRUE, scale =TRUE)
redlist$xy_coded<-as.factor(redlist$xy_coded)
redlist$Species_ID<-as.factor(redlist$Species_ID)
redlist$N_value <- as.factor(redlist$N_value)
redlist$T_value <- as.factor(redlist$T_value)
redlist$K_value <- as.factor(redlist$K_value)
redlist$L_value <- as.factor(redlist$L_value)
redlist$F_value <- as.factor(redlist$F_value)
redlist$R_value <- as.factor(redlist$R_value)
redlist$H_Rivers_and_Lakes <-as.factor(redlist$H_Rivers_and_Lakes )
redlist$H_Shorelines       <-as.factor(redlist$H_Shorelines       )
redlist$H_Bogs_and_mires   <-as.factor(redlist$H_Bogs_and_mires   )
redlist$H_Rocks_and_debris <-as.factor(redlist$H_Rocks_and_debris )  
redlist$H_Dry_meadows      <-as.factor(redlist$H_Dry_meadows      )
redlist$H_Rich_meadows     <-as.factor(redlist$H_Rich_meadows     )
redlist$H_Alpine_pastures  <-as.factor(redlist$H_Alpine_pastures  )
redlist$H_Herbaceous_fringe<-as.factor(redlist$H_Herbaceous_fringe) 
redlist$H_Shrubs_hedges    <-as.factor(redlist$H_Shrubs_hedges    )
redlist$H_Forests          <-as.factor(redlist$H_Forests          )
redlist$H_Ruderal_areas    <-as.factor(redlist$H_Ruderal_areas    )
redlist$H_Cropfields       <-as.factor(redlist$H_Cropfields       ) 

str(redlist)


####Model 1: IUCN####
#The bobyqa addition helps to get rid of the convergence warning
#yeardif scaled is added as a covariate. 
#Grid cell ID in random term is coded, data protection because of rare species

library(lme4)
m<-glmer(presence_01~
            + yeardif_s
            +IUCNorder
            +(1|xy_coded)+(1|Species_ID)
            ,data=redlist, family="binomial",
            glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))   

m1<-glmer(presence_01~
          + yeardif_s
        # +IUCNorder
         +(1|xy_coded)+(1|Species_ID)
         ,data=redlist, family="binomial",
         glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 

summary(m)
anova(m,m1)

stepba(m)  
str(redlist)


#####Models 2: habitat types####
#We run one model for each habitat type 
m1<-glmer(presence_01~yeardif_s+H_Shorelines +(1|xy_coded)+(1|Species_ID)
          ,data=redlist[!is.na(redlist$H_Shorelines),], family="binomial")
drop1(m1, test="Chisq") #gives the significancies (P, Chi2, Df) of the factors
summary(m1)  #estimates

m2<-glmer(presence_01~yeardif_s+H_Bogs_and_mires +(1|xy_coded)+(1|Species_ID)
          ,data=redlist[!is.na(redlist$H_Bogs_and_mires),], family="binomial")
drop1(m2, test="Chisq")
summary(m2)
m3<-glmer(presence_01~yeardif_s+H_Rocks_and_debris +(1|xy_coded)+(1|Species_ID)
          ,data=redlist[!is.na(redlist$H_Rocks_and_debris),], family="binomial")
drop1(m3, test="Chisq")
summary(m3)
m4<-glmer(presence_01~yeardif_s+H_Dry_meadows +(1|xy_coded)+(1|Species_ID)
          ,data=redlist[!is.na(redlist$H_Dry_meadows),], family="binomial")
drop1(m4, test="Chisq")
summary(m4)
m5<-glmer(presence_01~yeardif_s+H_Rich_meadows +(1|xy_coded)+(1|Species_ID)
          ,data=redlist[!is.na(redlist$H_Rich_meadows),], family="binomial")
drop1(m5, test="Chisq")
summary(m5)
m6<-glmer(presence_01~yeardif_s+H_Alpine_pastures +(1|xy_coded)+(1|Species_ID)
          ,data=redlist[!is.na(redlist$H_Alpine_pastures),], family="binomial")
drop1(m6, test="Chisq")
summary(m6)

m7<-glmer(presence_01~yeardif_s+H_Herbaceous_fringe +(1|xy_coded)+(1|Species_ID)
          ,data=redlist[!is.na(redlist$H_Herbaceous_fringe),], family="binomial")
drop1(m7, test="Chisq")
summary(m7)
m8<-glmer(presence_01~yeardif_s+H_Shrubs_hedges +(1|xy_coded)+(1|Species_ID)
          ,data=redlist[!is.na(redlist$H_Shrubs_hedges),], family="binomial")
drop1(m8, test="Chisq")
summary(m8)

m9<-glmer(presence_01~yeardif_s+H_Forests +(1|xy_coded)+(1|Species_ID)
          ,data=redlist[!is.na(redlist$H_Forests),], family="binomial")
drop1(m9, test="Chisq")
summary(m9)

m10<-glmer(presence_01~yeardif_s+H_Ruderal_areas +(1|xy_coded)+(1|Species_ID)
           ,data=redlist[!is.na(redlist$H_Ruderal_areas),], family="binomial")
drop1(m10, test="Chisq")
summary(m10)

m11<-glmer(presence_01~yeardif_s+H_Cropfields +(1|xy_coded)+(1|Species_ID)
           ,data=redlist[!is.na(redlist$H_Cropfields),], family="binomial")
drop1(m11, test="Chisq")
summary(m11)
m12<-glmer(presence_01~yeardif_s+H_Rivers_and_Lakes +(1|xy_coded)+(1|Species_ID)
           ,data=redlist[!is.na(redlist$H_Rivers_and_Lakes),], family="binomial")
drop1(m12, test="Chisq")
summary(m12)

#Get estimates
eff1<-data.frame(effect("H_Shorelines", m1))
eff2<-effect("H_Bogs_and_mires", m2)
eff2<-data.frame(eff2)
eff3<-effect("H_Rocks_and_debris", m3)
eff3<-data.frame(eff3)
eff4<-effect("H_Dry_meadows", m4)
eff4<-data.frame(eff4)
eff5<-effect("H_Rich_meadows", m5)
eff5<-data.frame(eff5)
eff6<-effect("H_Alpine_pastures", m6)
eff6<-data.frame(eff6)
eff7<-effect("H_Herbaceous_fringe", m7)
eff7<-data.frame(eff7)
eff8<-effect("H_Shrubs_hedges", m8)
eff8<-data.frame(eff8)
eff9<-effect("H_Forests", m9)
eff9<-data.frame(eff9)
eff10<-effect("H_Ruderal_areas", m10)
eff10<-data.frame(eff10)
eff11<-effect("H_Cropfields", m11)
eff11<-data.frame(eff11)
eff12<-effect("H_Rivers_and_Lakes", m12)
eff12<-data.frame(eff12)

eff1<-eff1[c(0,2),-1]
eff2<-eff2[c(0,2),-1]
eff3<-eff3[c(0,2),-1]
eff4<-eff4[c(0,2),-1]
eff5<-eff5[c(0,2),-1]
eff6<-eff6[c(0,2),-1]
eff7<-eff7[c(0,2),-1]
eff8<-eff8[c(0,2),-1]
eff9<-eff9[c(0,2),-1]
eff10<-eff10[c(0,2),-1]
eff11<-eff11[c(0,2),-1]
eff12<-eff12[c(0,2),-1]

#Habibitat types
habitat<-c("Crop fields and vineyards","Shorelines", "Bogs and mires", "Rocks and debris", 
           "Dry meadows and pastures", "Rich meadows and pastures", "Alpine pastures",
           "Herbaceous fringe", "Shrubs and hedges", "Forests", "Ruderal areas", "Rivers and lakes")

#Combine
dat<-rbind(eff11, eff1, eff2, eff3, eff4, eff5, eff6, eff7, eff8, eff9, eff10, eff12)
dat2<-cbind(dat, habitat)
dat2$fit<-dat2$fit*100
dat2$upper<-dat2$upper*100
dat2$lower<-dat2$lower*100
dat2$se<-dat2$se*100
dat2$ex_fit<-100-dat2$fit
dat2$ex_lower<-100-dat2$lower
dat2$ex_upper<-100-dat2$upper

#Code for Figure 4. 
dat2$habitat<-factor(dat2$habitat, levels =dat2$habitat[order(dat2$fit)])
ggplot(dat2, aes(x=habitat, y=ex_fit))+
  ylim(0,80)+
  geom_errorbar(aes(ymin=ex_lower, ymax=ex_upper), width=.0) +
  geom_point(shape=22, size=6, fill=c("#FFC107", "#0288d1", "#0288d1", "#FFC107", "#b8e186", "#b8e186", "#b8e186", "#b8e186", "#008837",
                                      "#008837", "#FFC107", "#0288d1")) +
  coord_flip()+
  theme_bw()+
  theme(axis.line = element_line(size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),  
        text = element_text(size=18),
        axis.ticks.y = element_blank(),
        axis.line.y= element_blank(),
        axis.text = element_text(size=16))+
  ylab("Proportion of extinct populations")+
  xlab("")+
  geom_text(x = 1, y = 78, label="66", size =4, colour="darkgrey")+
  geom_text(x = 2, y = 78, label="73", size =4,colour="darkgrey")+
  geom_text(x = 3, y = 78, label="95", size =4,colour="darkgrey")+
  geom_text(x = 4, y = 78, label="63", size =4,colour="darkgrey")+
  geom_text(x = 5, y = 78, label="107", size =4,colour="darkgrey")+
  geom_text(x = 6, y = 78, label="133", size =4,colour="darkgrey")+
  geom_text(x = 7, y = 78, label="100", size =4,colour="darkgrey")+
  geom_text(x = 8, y = 78, label="17", size =4,colour="darkgrey")+
  geom_text(x = 9, y = 78, label="113", size =4,colour="darkgrey")+
  geom_text(x = 10, y = 78, label="38", size =4,colour="darkgrey")+
  geom_text(x = 11, y = 78, label="74", size =4,colour="darkgrey")+
  geom_text(x = 12, y = 78, label="84", size =4,colour="darkgrey")



#______________________

####Code for Figure 2: IUCN Order####
#use yeardiff_s as covariate
redlist_n<-redlist[!redlist$IUCNorder=="1_DD",]  #exclude species with DD (data defficient)

m1<-glmer(presence_01~
            +yeardif_s
          +IUCNorder
          +(1|xy_coded)+(1|Species_ID)
          ,data=redlist_n, family="binomial")
summary(m1)
allEffects(m1)

plot(allEffects(m1))
anova(m2,m1)

eff<-effect("IUCNorder", m1)
eff<-data.frame(eff)

eff$extinction<-1-eff$fit
eff$ex_lower<-1-eff$lower
eff$ex_upper<-1-eff$upper
eff$N_cat<-c(8,94,214,268,110,17)

#Plot
par(mfrow=c(1,1))
par(mar=c(5,6,4,1)+.1)
x1<-c("RE", "CR", "EN", "VU", "NT", "LC")
x<-c(1,2,3,4,5,6)
y1<-c("0", "25", "50", "75", "100")

# http://colorbrewer2.org/#type=sequential&scheme=Reds&n=5     # colors from this website
plot(x,eff$extinction, pch= 22,  bg= c("#a50f15", "#de2d26", "#fcae91", "#fee5d9", "grey81", "grey93") ,  cex=2.3,  ylim=c( 0,1.0), xlim=c(1, 7), ylab = "",  xaxt="n", 
     bty="n", xlab="IUCN threat category", yaxt="n", cex.lab=1.6)
axis(side=1,at=c(1,2,3,4,5,6), labels =x1, cex.axis=1.3 ) #specify tick marks to be drawn at 1,2,3,4
axis(side=2, at=c(0, 0.25, 0.5, 0.75, 1), labels = y1, cex.axis=1.3)
title(ylab="Proportion of extinct populations", line = 4, cex.lab = 1.6)
segments(x, eff$ex_lower,x,  eff$ex_upper)
text(x, 0, eff$N_cat, cex=0.8)


####Code for Figure 2: yeardif (now yeardif is not scaled) ####
head(redlist)
m1<-glmer(presence_01~
          +yeardifnew
          +IUCNorder
          +(1|xy_coded)+(1|Species_ID)
          ,data=redlist, family="binomial")
summary(m1)
plot(allEffects(m1))
allEffects(m1)
eff<-effect("yeardifnew", m1, xlevels=(list(yeardifnew=c(9,15,20,30,40,50))))
eff<-data.frame(eff)
eff$fit<-eff$fit*100
eff$lower<-eff$lower*100
eff$upper<-eff$upper*100
eff$ex_fit<-100-eff$fit
eff$ex_lower<-100-eff$lower
eff$ex_upper<-100-eff$upper
eff
y1<-c("0%", "25%", "50%", "75%", "100%")

str(redlist)
levels(redlist$IUCNorder)
#Make a variable for the jitter
redlist$Graphpresence<-as.factor(redlist$presence_01)
levels(redlist$Graphpresence)
levels(redlist$Graphpresence)[levels(redlist$Graphpresence)=="0"] <- "100"
levels(redlist$Graphpresence)[levels(redlist$Graphpresence)=="1"] <- "0"
redlist$Graphpresence<-as.numeric(as.character(redlist$Graphpresence))
hist(redlist$Graphpresence) 

y<-ggplot(data=eff, aes(x=yeardifnew, y=ex_fit))+
  xlim(9,50) +  ylim(0, 100) +
  geom_line(size=1.2)+
  theme_bw() +
  theme(axis.line = element_line(colour = "black",size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(text = element_text(size=18))+
  xlab("Years since last record")+
  ylab("Proportion of extinct populations")+
  theme(legend.position="bottom")+
  geom_ribbon(aes(ymin=eff$ex_lower, ymax=eff$ex_upper), linetype=2, alpha=0.1)+
  geom_jitter(data=redlist, aes(x=yeardifnew, y=Graphpresence) ,alpha=0.2, height=5, width =0.5, colour="grey")
y


#__________________________________________________

#### Model 3: Traits and indicator values ####
#Grimes C and R, height, all indicator values and yeardif_s. 
str(redlist)
m<-glmer(presence_01~
           +yeardif_s
         +C_s
         +R_s
         + Maxheight_s
         +Kn_s
         + Nn_s          
         + Fn_s          
         + Tn_s           
         + Ln_s            
         + Rn_s 
         +(1|xy_coded)+(1|Species_ID)
         ,data=redlist[!is.na(redlist$Ln_s),], family="binomial",
         glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 
summary(m)

m2<-glmer(presence_01~yeardif_s+C_s +R_s+ Maxheight_s+Kn_s+ Nn_s + Fn_s + Tn_s + Rn_s 
          +(1|xy_coded)+(1|Species_ID),data=redlist[!is.na(redlist$Maxheight_s),], family="binomial"
          ,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 
summary(m2)

#remove height
m3<-glmer(presence_01~yeardif_s+C_s +R_s+ Kn_s+ Nn_s + Fn_s + Tn_s + Rn_s 
          +(1|xy_coded)+(1|Species_ID),data=redlist[!is.na(redlist$Maxheight_s),], family="binomial"
          ,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 
anova(m2,m3)

m3<-glmer(presence_01~yeardif_s+C_s +R_s+ Kn_s+ Nn_s + Fn_s + Tn_s + Rn_s 
          +(1|xy_coded)+(1|Species_ID),data=redlist, family="binomial"
          ,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 
summary(m3)

drop1(m3)

#C_s
m4<-glmer(presence_01~yeardif_s+R_s+ Kn_s+ Nn_s + Fn_s + Tn_s + Rn_s 
          +(1|xy_coded)+(1|Species_ID),data=redlist, family="binomial"
          ,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 
anova(m3,m4)

#Final model m4
summary(m4)

#Kn
m5<-glmer(presence_01~yeardif_s+R_s+  Nn_s + Fn_s + Tn_s + Rn_s 
          +(1|xy_coded)+(1|Species_ID),data=redlist, family="binomial"
          ,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 
anova(m4,m5)
#Fn
m6<-glmer(presence_01~yeardif_s+R_s+ Kn_s+ Nn_s +  Tn_s + Rn_s 
          +(1|xy_coded)+(1|Species_ID),data=redlist, family="binomial"
          ,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 
anova(m4,m6)

#Tn
m7<-glmer(presence_01~yeardif_s+R_s+ Kn_s+ Nn_s + Fn_s +  Rn_s 
          +(1|xy_coded)+(1|Species_ID),data=redlist, family="binomial"
          ,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 
anova(m4,m7)

#Rs
m8<-glmer(presence_01~yeardif_s+ Kn_s+ Nn_s + Fn_s + Tn_s + Rn_s 
          +(1|xy_coded)+(1|Species_ID),data=redlist, family="binomial"
          ,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 
anova(m4,m8)

#Nn
m9<-glmer(presence_01~yeardif_s+R_s+ Kn_s+  Fn_s + Tn_s + Rn_s 
          +(1|xy_coded)+(1|Species_ID),data=redlist, family="binomial"
          ,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))  
anova(m4,m9)

#Rn
m4_<-glmer(presence_01~yeardif_s+R_s+ Kn_s+ Nn_s + Fn_s + Tn_s + Rn_s 
           +(1|xy_coded)+(1|Species_ID),data=redlist[!is.na(redlist$Rn_s),], family="binomial"
           ,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 

m10<-glmer(presence_01~yeardif_s+R_s+ Kn_s+ Nn_s + Fn_s + Tn_s 
           +(1|xy_coded)+(1|Species_ID),data=redlist[!is.na(redlist$Rn_s),], family="binomial"
           ,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 
anova(m4_,m10)


#### Code for the Figure S2 (species characteristics and indicator values) #### 
#Fn - scale everything except Fn
mFn<-glmer(presence_01~
                  +yeardif_s
                  +R_s
                  +Kn_s
                  +Nn_s         
                  +Fn          
                  +Tn_s           
                  +Rn_s 
                  +(1|xy_coded)+(1|Species_ID)
                  ,data=redlist, family="binomial",
                  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))   
       
eff_F<-data.frame(effect("Fn", mFn))
eff_F$fit<-eff_F$fit*100
eff_F$lower<-eff_F$lower*100
eff_F$upper<-eff_F$upper*100
eff_F$ex_fit<-100-eff_F$fit
eff_F$ex_lower<-100-eff_F$lower
eff_F$eff_upper<-100-eff_F$upper
       
#Plot F, only continous  
Fn<-ggplot(data=eff_F, aes(x=Fn, y=ex_fit))+
    xlim(0.9,5.1) +
    geom_line(size=1.2)+
    theme_bw() +
    theme(axis.line = element_line(colour = "black",size = 1),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),   
               axis.title.y=element_blank())+
         scale_y_continuous(limits=c(0,69), breaks = c(0,10,20,30,40,50)) +  
         theme(text = element_text(size=16))+
         xlab("Moisture level")+
         ylab("Proportio of local extinctions")+
         theme(legend.position="bottom")+
         geom_ribbon(aes(ymin=eff_F$ex_lower, ymax=eff_F$eff_upper), linetype=2, alpha=0.1)
Fn
 
#add the categorical values as well:
mFn_k<-glmer(presence_01~
                      +yeardif_s
                    +R_s
                    +Kn_s
                    + Nn_s          
                    + F_value          
                    + Tn_s           
                    + Rn_s 
                    +(1|xy_coded)+(1|Species_ID)
                    ,data=redlist, family="binomial",
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))  
eff_Fk<-data.frame(effect("F_value", mFn_k))
eff_Fk$fit<-eff_Fk$fit*100
eff_Fk$lower<-eff_Fk$lower*100
eff_Fk$upper<-eff_Fk$upper*100
eff_Fk$ex_fit<-100-eff_Fk$fit
eff_Fk$ex_lower<-100-eff_Fk$lower
eff_Fk$ex_upper<-100-eff_Fk$upper
str(eff_Fk)      
eff_Fk$F_value<-as.numeric(as.character(eff_Fk$F_value))
       
#plot F, for continous variables and categorical variables
Fnew<-ggplot (data=eff_F, aes(x=Fn, y=ex_fit))+
         geom_line(size=1.2) +
         #xlim(0.9,5.1) +
         theme_bw() +
         theme(axis.line = element_line(colour = "black",size = 1),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),   
               axis.title.y=element_blank())+
         scale_y_continuous(limits=c(0,76), breaks = c(0,10,20,30,40,50, 60)) +  
         theme(text = element_text(size=16))+
         xlab("Moisture level")+
         #ylab("Proportio of extinct populations")+
         theme(legend.position="bottom")+
         geom_ribbon(data=eff_F, aes(x=Fn, ymin=ex_lower, ymax=eff_upper), linetype=2, alpha=0.1)+
         geom_errorbar(data=eff_Fk, aes(x=F_value, ymin=ex_lower, ymax=ex_upper), width=.0, colour="darkgrey") +
         geom_point(data=eff_Fk, aes(x=F_value, y=ex_fit), shape=22, size=4, colour="darkgrey")
Fnew
       
#********************************      
#Plot Nn
mNn<-glmer(presence_01~
                    +yeardif_s
                  +R_s
                  +Kn_s
                  + Nn          
                  + Fn_s          
                  + Tn_s           
                  + Rn_s 
                  +(1|xy_coded)+(1|Species_ID)
                  ,data=redlist, family="binomial",
                  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))    
eff_N<-effect("Nn", mNn)
eff_N<-data.frame(eff_N)
eff_N$fit<-eff_N$fit*100
eff_N$lower<-eff_N$lower*100 
eff_N$upper<-eff_N$upper*100
eff_N$ex_fit<-100-eff_N$fit
eff_N$ex_lower<-100-eff_N$lower
eff_N$ex_upper<-100-eff_N$upper
       
#get categorical estimates
m19_k<-glmer(presence_01~
                      +yeardif_s
                    +R_s
                    +Kn_s
                    + N_value          
                    + Fn_s          
                    + Tn_s           
                    + Rn_s 
                    +(1|xy_coded)+(1|Species_ID)
                    ,data=redlist, family="binomial",
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))    
eff_Nk<-effect("N_value", m19_k)
eff_Nk<-data.frame(eff_Nk)
eff_Nk$fit<-eff_Nk$fit*100
eff_Nk$lower<-eff_Nk$lower*100
eff_Nk$upper<-eff_Nk$upper*100
eff_Nk$ex_fit<-100-eff_Nk$fit
eff_Nk$ex_lower<-100-eff_Nk$lower
eff_Nk$ex_upper<-100-eff_Nk$upper
eff_Nk$N_value<-as.numeric(as.character(eff_Nk$N_value))
eff_Nk
       
#combined plot
Nn<-ggplot(data=eff_N, aes(x=Nn, y=ex_fit))+
         xlim(0.9,5.1) +
         geom_line(size=1.2)+
         theme_bw() +
         theme(axis.line = element_line(colour = "black",size = 1),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               axis.title.y=element_blank(),    
               text = element_text(size=16))+
         scale_y_continuous(limits=c(0,76), breaks = c(0,10,20,30,40,50,60)) +   
         xlab("Nutrient level")+
         # ylab("Proportion of extinct populations")+
         theme(legend.position="bottom")+
         geom_ribbon(aes(ymin=eff_N$ex_lower, ymax=eff_N$ex_upper), linetype=2, alpha=0.1)  +
         geom_errorbar(data=eff_Nk, aes(x=N_value, ymin=ex_lower, ymax=ex_upper), width=.0, colour="darkgrey") +
         geom_point(data=eff_Nk, aes(x=N_value, y=ex_fit), shape=22, size=4, colour="darkgrey")
Nn
       
#Add Line for model predictions when highest N=5 are excluded
#Nn without N=5
mN5<-glmer(presence_01~
                    +yeardif_s
                  +R_s
                  +Kn_s
                  + Nn          
                  + Fn_s          
                  + Tn_s           
                  + Rn_s 
                  +(1|xy_coded)+(1|Species_ID)
                  ,data=redlist[!redlist$Nn=="5",], family="binomial",
                  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 
eff_Nex<-effect("Nn", mN5)
eff_Nex<-data.frame(eff_Nex)
eff_Nex$fit<-eff_Nex$fit*100
eff_Nex$lower<-eff_Nex$lower*100
eff_Nex$upper<-eff_Nex$upper*100
eff_Nex$ex_fit<-100-eff_Nex$fit
eff_Nex$ex_lower<-100-eff_Nex$lower
eff_Nex$ex_upper<-100-eff_Nex$upper
eff_Nex
       
#Make another graph excluding N=5
Nnex<-ggplot(data=eff_Nex, aes(x=Nn, y=ex_fit))+
         xlim(0.9,5.1) +
         geom_line(size=1.2)+
         theme_bw() +
         theme(axis.line = element_line(colour = "black",size = 1),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               axis.title.y=element_blank(),    
               text = element_text(size=16))+
         scale_y_continuous(limits=c(0,73), breaks = c(0,10,20,30,40,50,60)) +   
         xlab("Nutrient level")+
         # ylab("Proportion of local extinctions")+
         theme(legend.position="bottom")+
         geom_ribbon(data = eff_Nex[eff_Nex$Nn < 5, ],aes(ymin=eff_Nex$ex_lower, ymax=eff_Nex$ex_upper), linetype=2, alpha=0.1)  +
         geom_errorbar(data=eff_Nk, aes(x=N_value, ymin=ex_lower, ymax=ex_upper), width=.0, colour="darkgrey") +
         geom_point(data=eff_Nk, aes(x=N_value, y=ex_fit), shape=22, size=4, colour="darkgrey")
Nnex
       
       
#*****************
#Plot Rn
mRn<-glmer(presence_01~
                    +yeardif_s
                  +R_s
                  +Kn_s
                  + Nn_s          
                  + Fn_s          
                  + Tn_s           
                  + Rn 
                  +(1|xy_coded)+(1|Species_ID)
                  ,data=redlist, family="binomial",
                  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))    
eff_R<-effect("Rn", mRn)
eff_R<-data.frame(eff_R)
eff_R$fit<-eff_R$fit*100
eff_R$lower<-eff_R$lower*100
eff_R$upper<-eff_R$upper*100
eff_R$ex_fit<-100-eff_R$fit
eff_R$ex_lower<-100-eff_R$lower
eff_R$ex_upper<-100-eff_R$upper
       
#categorical
mRn_k<-glmer(presence_01~
                      +yeardif_s
                    +R_s
                    +Kn_s
                    + Nn_s          
                    + Fn_s          
                    + Tn_s           
                    + R_value
                    +(1|xy_coded)+(1|Species_ID)
                    ,data=redlist, family="binomial",
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))    
eff_Rk<-effect("R_value", mRn_k)
eff_Rk<-data.frame(eff_Rk)
eff_Rk$fit<-eff_Rk$fit*100
eff_Rk$lower<-eff_Rk$lower*100
eff_Rk$upper<-eff_Rk$upper*100
eff_Rk$ex_fit<-100-eff_Rk$fit
eff_Rk$ex_lower<-100-eff_Rk$lower
eff_Rk$ex_upper<-100-eff_Rk$upper
eff_Rk$R_value<-as.numeric(as.character(eff_Rk$R_value))
       
#Plot
Rn<-ggplot(data=eff_R, aes(x=Rn, y=ex_fit))+
         xlim(0.9,5.1) + 
         geom_line(size=1.2)+
         theme_bw() +
         theme(axis.line = element_line(colour = "black",size = 1),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank()) +
         scale_y_continuous(limits=c(0,76), breaks = c(0,10,20,30,40,50,60)) +  
         theme(text = element_text(size=16),
               axis.title.y=element_blank()
               #  axis.text.y=element_blank()
         )+
         xlab("Reaction level")+
         ylab("Proportion of extinct populations")+
         theme(legend.position="bottom")+
         geom_ribbon(aes(ymin=eff_R$ex_lower, ymax=eff_R$ex_upper), linetype=2, alpha=0.1)+
         geom_errorbar(data=eff_Rk, aes(x=R_value, ymin=ex_lower, ymax=ex_upper), width=.0, colour="darkgrey") +
         geom_point(data=eff_Rk, aes(x=R_value, y=ex_fit), shape=22, size=4, colour="darkgrey")
Rn
       
#**********************
mTn<-glmer(presence_01~
                    +yeardif_s
                  +R_s
                  +Kn_s
                  + Nn_s          
                  + Fn_s          
                  + Tn           
                  + Rn_s 
                  +(1|xy_coded)+(1|Species_ID)
                  ,data=redlist, family="binomial",
                  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
eff_T<-effect("Tn", mTn)
eff_T<-data.frame(eff_T)
eff_T$fit<-eff_T$fit*100
eff_T$lower<-eff_T$lower*100
eff_T$upper<-eff_T$upper*100
eff_T$ex_fit<-100-eff_T$fit
eff_T$ex_lower<-100-eff_T$lower
eff_T$ex_upper<-100-eff_T$upper
       
#categorical
mTn_k<-glmer(presence_01~
                      +yeardif_s
                    +R_s
                    +Kn_s
                    + Nn_s          
                    + Fn_s          
                    + T_value          
                    + Rn_s 
                    +(1|xy_coded)+(1|Species_ID)
                    ,data=redlist, family="binomial",
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
eff_Tk<-effect("T_value", mTn_k)
eff_Tk<-data.frame(eff_Tk)
eff_Tk$fit<-eff_Tk$fit*100
eff_Tk$lower<-eff_Tk$lower*100
eff_Tk$upper<-eff_Tk$upper*100
eff_Tk$ex_fit<-100-eff_Tk$fit
eff_Tk$ex_lower<-100-eff_Tk$lower
eff_Tk$ex_upper<-100-eff_Tk$upper
eff_Tk$T_value<-as.numeric(as.character(eff_Tk$T_value))
       
#combined graph
Tn<-ggplot(data=eff_T, aes(x=Tn, y=ex_fit))+
         xlim(0.9,5.1) + 
         geom_line(size=1.2)+
         theme_bw() +
         theme(axis.line = element_line(colour = "black",size = 1),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               axis.title.y=element_blank()) +
         theme(text = element_text(size=16))+
         scale_y_continuous(limits=c(0,73), breaks = c(0,10,20,30,40,50,60)) +  
         xlab("Temperature level")+
         ylab("Proportion of extinct populations")+
         theme(legend.position="bottom")+
         geom_ribbon(aes(ymin=eff_T$ex_lower, ymax=eff_T$ex_upper), linetype=2, alpha=0.1)+
         geom_errorbar(data=eff_Tk, aes(x=T_value, ymin=ex_lower, ymax=ex_upper), width=.0, colour="darkgrey") +
         geom_point(data=eff_Tk, aes(x=T_value, y=ex_fit), shape=22, size=4, colour="darkgrey")
Tn
       

#**********************
#Kn
mKn<-glmer(presence_01~
                  +yeardif_s
                  +R_s
                  +Kn
                  + Nn_s          
                  + Fn_s          
                  + Tn_s          
                  + Rn_s 
                  +(1|xy_coded)+(1|Species_ID)
                  ,data=redlist, family="binomial",
                  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
eff_K<-effect("Kn", mKn)
eff_K<-data.frame(eff_K)
eff_K$fit<-eff_K$fit*100
eff_K$lower<-eff_K$lower*100
eff_K$upper<-eff_K$upper*100
eff_K$ex_fit<-100-eff_K$fit
eff_K$ex_lower<-100-eff_K$lower
eff_K$ex_upper<-100-eff_K$upper
       
#categorical
mKn_k<-glmer(presence_01~
                      +yeardif_s
                    +R_s
                    +K_value
                    + Nn_s          
                    + Fn_s          
                    + Tn_s          
                    + Rn_s 
                    +(1|xy_coded)+(1|Species_ID)
                    ,data=redlist, family="binomial",
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
eff_Kk<-effect("K_value", mKn_k)
eff_Kk<-data.frame(eff_Kk)
eff_Kk$fit<-eff_Kk$fit*100
eff_Kk$lower<-eff_Kk$lower*100
eff_Kk$upper<-eff_Kk$upper*100
eff_Kk$ex_fit<-100-eff_Kk$fit
eff_Kk$ex_lower<-100-eff_Kk$lower
eff_Kk$ex_upper<-100-eff_Kk$upper
eff_Kk$K_value<-as.numeric(as.character(eff_Kk$K_value))
eff_Kk       
#combined plot
Kn<-ggplot(data=eff_K, aes(x=Kn, y=ex_fit))+
         xlim(0.9,5.1) +
         geom_line(size=1.2)+
         theme_bw() +
         theme(axis.line = element_line(colour = "black",size = 1),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               axis.title.y=element_blank()) +
         theme(text = element_text(size=16))+
         scale_y_continuous(limits=c(0,73), breaks = c(0,10,20,30,40,50,60)) +  
         xlab("Continentality level")+
         ylab("Percentage of confirmed populations")+
         theme(legend.position="bottom")+
         geom_ribbon(aes(ymin=eff_K$ex_lower, ymax=eff_K$ex_upper), linetype=2, alpha=0.1)+
         geom_errorbar(data=eff_Kk, aes(x=K_value, ymin=ex_lower, ymax=ex_upper), width=.0, colour="darkgrey") +
         geom_point(data=eff_Kk, aes(x=K_value, y=ex_fit), shape=22, size=4, colour="darkgrey")
Kn
       
#Grime R Graph as ggplot Graph - use R non-scaled, as a factor.
m23<-glmer(presence_01~
                    +yeardif_s
                  +R
                  +Kn_s
                  + Nn_s          
                  + Fn_s          
                  + Tn_s           
                  + Rn_s 
                  +(1|xy_coded)+(1|Species_ID)
                  ,data=redlist, family="binomial",
                  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))  
       
eff_Ru<-effect("R", m23)
eff_Ru<-data.frame(eff_Ru)
eff_Ru$fit<-eff_Ru$fit*100
eff_Ru$lower<-eff_Ru$lower*100
eff_Ru$upper<-eff_Ru$upper*100
eff_Ru$ex_fit<-100-eff_Ru$fit
eff_Ru$ex_lower<-100-eff_Ru$lower
eff_Ru$ex_upper<-100-eff_Ru$upper
       
R<-ggplot(eff_Ru, aes(x=R, y=ex_fit))+
         geom_errorbar(aes(ymin=ex_lower, ymax=ex_upper), width=.0) +
         geom_point(shape=22, size=4) +
         theme_bw()+
         theme(axis.line = element_line(colour = "black",size = 1),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               text = element_text(size=16),
               axis.text.x = element_text(size=14), 
               axis.title.y=element_blank())+
         scale_y_continuous(limits=c(0,73), breaks = c(0,10,20,30,40,50, 60)) +  
         xlab("Degree of ruderality")
       #+   ylab("Percentage of confirmed populations")
R
       
       
#**********
#Combine in one graph: load packages
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggpubr)
       
library(grid)
library(ggplot2)
library(gridBase)
library(gridExtra)
library(ggpubr)
       
       
figure<-ggarrange(Fnew, Nn ,  Rn , Tn, Kn, R,
                         labels = c("", "" , ""),
                         ncol = 3, nrow = 2)
figure
annotate_figure(figure,
                       left =text_grob("Proportion of extinct populations", rot = 90, 
                                       face="bold", size=18),
                       fig.lab="", fig.lab.face="bold")
       
       
#### ####
       
      
  
       
  


