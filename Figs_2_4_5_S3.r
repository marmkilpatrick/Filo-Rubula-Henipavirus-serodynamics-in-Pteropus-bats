#Figures 1 (near bottom), 2, 4, 5, S1(near bottom), S2(near bottom), S3, S4, S5
lapply(c("tidyverse","ggthemes","viridis","mixtools","mclust","mgcv","gratia"),require,character.only=T) #load packages

setwd("c:/marm/research/other/nipah/manyviruses/")
q1=read.csv("bat_serology_processed.csv")
q1$date=as.Date(q1$date)
q1=q1[order(q1$date),]
q1$antibody[q1$antibody=="Nipah"]="Henipavirus"
q1$antibody=factor(q1$antibody,levels=c("Henipavirus","Filovirus","Rubulavirus"))
q1$Infect=0;q1$Infect[q1$seropos]=1 #make numerical seropos column
q1$Rlocation=ifelse(q1$location=="chakhoria","chakhoria",
                    ifelse(q1$location=="ramnagar","ramnagar","faridpur"))

#Fraction of bats recaptured during each sampling period
q1FnpwNiv=q1[q1$Rlocation=="faridpur"& q1$antibody=="Henipavirus" & 
               !q1$age_class=="preweaned"&!q1$microchip_id=="NA_NA_NA",]
q1FnpwNiv$microchip_id[q1FnpwNiv$microchip_id=="NA_NA_NA"]=NA
q1FnpwNiv$recap=as.integer(duplicated(q1FnpwNiv$microchip_id)) #recap column
q1FnpwNivG = q1FnpwNiv %>% group_by(sheet_name) %>%  
  summarize(recapfrac= mean(recap),N=n(),Nrecap=sum(recap),date=mean(date))
# ggplot(q1FnpwNivG,aes(x=date,y=recapfrac,weight=N))+theme_few()+geom_point()+
#   stat_smooth(method="glm",method.args=list(family="binomial"))
mean(q1FnpwNivG$recapfrac)
summary(glm(Infect~recap+age_class,data=q1FnpwNiv,family="binomial")) #no diff in NiV serop for recaptures

#Mixture models to categorize seropositivity
NiVd=normalmixEM(log10(q1$luminex[q1$antibody=="Henipavirus"]));summary(NiVd)
source("cutoffcalc.r") #Function for finding intersection of positive and negative distributions to calculate cutoff
NiVcut=max(find_normal_intersections(m1=NiVd$mu[1],s1=NiVd$sigma[1],m2=NiVd$mu[2],s2=NiVd$sigma[2]))
NiVdmcL=Mclust(log10(q1$luminex[q1$antibody=="Henipavirus"]),G=1:2);NiVdmcL$BIC #NiVdmcL$parameters;NiVdmcL$classification
NiVcut2=max(find_normal_intersections(m1=NiVdmcL$parameters$mean[1],s1=sqrt(NiVdmcL$parameters$variance$sigmasq[1]),
                                      m2=NiVdmcL$parameters$mean[2],s2=sqrt(NiVdmcL$parameters$variance$sigmasq[2]) ))

FiVd=normalmixEM(log10(q1$luminex[q1$antibody=="Filovirus"]),maxit=5000)
FiVcut=max(find_normal_intersections(m1=FiVd$mu[1],s1=FiVd$sigma[1],m2=FiVd$mu[2],s2=FiVd$sigma[2]));10^(FiVcut)
FiVdmcL=Mclust(log10(q1$luminex[q1$antibody=="Filovirus"]),G=1:2);FiVdmcL$BIC[2,2]-min(FiVdmcL$BIC[-4])
FiVcut2=max(find_normal_intersections(m1=FiVdmcL$parameters$mean[1],s1=sqrt(FiVdmcL$parameters$variance$sigmasq[1]),
                                      m2=FiVdmcL$parameters$mean[2],s2=sqrt(FiVdmcL$parameters$variance$sigmasq[2]) ))

RuVd=normalmixEM(log10(q1$luminex[q1$antibody=="Rubulavirus"]))
RuVcut=min(find_normal_intersections(m1=RuVd$mu[1],s1=RuVd$sigma[1],m2=RuVd$mu[2],s2=RuVd$sigma[2]));10^RuVcut
RuVdmcL=Mclust(log10(q1$luminex[q1$antibody=="Rubulavirus"]),G=1:2);RuVdmcL$BIC[2,2]-min(RuVdmcL$BIC[-4])
RuVcut2=max(find_normal_intersections(m1=RuVdmcL$parameters$mean[1],s1=sqrt(RuVdmcL$parameters$variance$sigmasq[1]),
                                      m2=RuVdmcL$parameters$mean[2],s2=sqrt(RuVdmcL$parameters$variance$sigmasq[1]) ))

cutoffs = data.frame(antibody = c("Henipavirus","Filovirus","Rubulavirus"),
                     cutoff = c(3,FiVcut,RuVcut))
cutoffs$antibody=factor(cutoffs$antibody,levels=c("Henipavirus","Filovirus","Rubulavirus"))

#Figure S3
LMIF=seq(1,4.5,by=0.1)
df1=data.frame(antibody=rep(c("Henipavirus","Filovirus","Rubulavirus"),each=2*length(LMIF)),
               LMIF=rep(LMIF,2),
               density=c(dnorm(LMIF,NiVd$mu[1],NiVd$sigma[1])*NiVd$lambda[1],
                         dnorm(LMIF,NiVd$mu[2],NiVd$sigma[2])*NiVd$lambda[2],
                         dnorm(LMIF,FiVd$mu[1],FiVd$sigma[1])*FiVd$lambda[1],
                         dnorm(LMIF,FiVd$mu[2],FiVd$sigma[2])*FiVd$lambda[2],
                         dnorm(LMIF,RuVd$mu[1],RuVd$sigma[1])*RuVd$lambda[1],
                         dnorm(LMIF,RuVd$mu[2],RuVd$sigma[2])*RuVd$lambda[2]),
               group=rep(LETTERS[1:6],each=length(LMIF)))
df1$antibody=factor(df1$antibody,levels=c("Henipavirus","Filovirus","Rubulavirus"))

#Fig S3
FigS3=ggplot(q1,aes(x=log10(luminex),color=antibody)) + theme_few()+
  geom_density()+
  geom_vline(data=cutoffs,aes(xintercept=cutoff,color=antibody),linetype="dashed")+
  geom_line(data=df1,aes(x=LMIF,y=density,color=antibody,group=group),linetype="dotted")+
  scale_color_viridis_d(option="C", begin=0.2, end=0.8) +
  facet_wrap(~antibody,nrow=3)+
  theme(text=element_text(size=20),legend.position="None")+
  labs(x=expression(log[10]~Luminex~MFI));FigS3
ggsave("Figure S3.pdf",plot=FigS3,width=8,height=6)

diagCutOffs=q1 %>% group_by(antibody) %>% 
  summarize(minP=min(luminex[Infect==1]),maxN=max(luminex[Infect==0]))
diagCutOffs$LogCutOff=with(diagCutOffs,log10((minP+maxN)/2));diagCutOffs

q1$jday=yday(q1$date)
q1$wk=as.numeric(strftime(q1$date,format="%V")) #make week column
q1$yr=as.numeric(format(q1$date,"%Y")) #make year column
q1$week=q1$wk+(q1$yr-2006)*52 #make running week, starting w/ 2006
q2=q1[!q1$age_class=="preweaned"&!q1$location%in%c("chakhoria", "ramnagar"),] #subset to Faridpur
q2=q2 %>% arrange(date) %>% mutate(sess = cumsum(c(TRUE, diff(date) > 15)))
q2$pdate=as.Date(ymd("2007-01-01") + q2$jday) #yday("2007-06-01") #make wrapping date for plot
q2$pdate[q2$jday<152]=as.Date(ymd("2008-01-01")+q2$jday[q2$jday<152])

#Data for Figure 2 seasonal prev over time for juveniles
qA = data.frame(q2 %>% group_by(sess,yr,age_class,antibody) %>% 
                  summarize(N=n(),prev=mean(Infect),pdate=mean(pdate),jday=mean(jday)))

for (i in 1:nrow(qA)) {
  qA$L95[i]=binom.test(qA$prev[i]*qA$N[i],qA$N[i])$conf.int[1]
  qA$U95[i]=binom.test(qA$prev[i]*qA$N[i],qA$N[i])$conf.int[2] }

qAj=qA[qA$age_class=="juvenile",] #Juvenile data
qAj$jday[qAj$jday<152]=qAj$jday[qAj$jday<152]+365 #Make julian day starting in June

g0=gam(prev~antibody+s(jday),weights=N,gamma=3,data=qAj,family="binomial");summary(g0)
g1=gam(prev~antibody+s(jday,by=antibody),weights=N,gamma=3,data=qAj,family="binomial");summary(g1)
AIC(g0,g1); AIC(g0)-AIC(g1)
pdays=seq(range(qAj$jday)[1],range(qAj$jday)[2],by=1)
gamder=derivatives(g1,data=data.frame(jday=rep(pdays,3),
                                      antibody=rep(c("Henipavirus","Filovirus","Rubulavirus"),each=length(pdays))))
draw(gamder)
gamder$pdate=as.Date(ymd("2007-01-01") + gamder$jday) #yday("2007-06-01") #make wrapping date for plot
#gamder$pdate[gamder$jday<152]=as.Date(ymd("2008-01-01")+gamder$jday[gamder$jday<152])
# qAj = merge(qAj,slopes[,c("jday",".derivative",".lower_ci",".upper_ci","antibody")],by=c("jday","antibody"))
gamder$sig=ifelse(sign(gamder$.lower_ci)==sign(gamder$.upper_ci),"red","black")
gamder$prev=predict(g1,newdata=gamder,type="response")
gamder$prevlink=predict(g1,newdata=gamder,type="link")
gamder$SE=predict(g1,newdata=gamder,type="link",se.fit=T)$se
gamder$L95=plogis(gamder$prevlink-1.96*gamder$SE);gamder$U95=plogis(gamder$prevlink+1.96*gamder$SE)
gamder$antibody=factor(gamder$antibody,levels=c("Henipavirus","Filovirus","Rubulavirus"))
  
#Figure 2 Longitudinal juvenile serodynamics
Fig2=ggplot(qAj,aes(x=pdate,y=prev))+theme_few()+
  facet_wrap(~antibody,nrow=3)+
  scale_y_continuous(labels=scales::percent)+
  geom_ribbon(data=gamder,aes(x=pdate,y=prev,ymin=L95,ymax=U95,fill=antibody,color=antibody))+
  geom_path(data=gamder,aes(x=pdate,y=prev),color=gamder$sig)+
  geom_pointrange(aes(ymin=L95,ymax=U95),size=.9)+
  labs(x="Date",y="Prevalence")+
  scale_x_date(date_breaks="1 month",date_labels="%b",expand=c(0,0))+
  scale_color_viridis_d(option="C", begin=0.2, end=0.8) +
  scale_fill_viridis_d(option="C", begin=0.2, end=0.8)+
  theme(text=element_text(size=20),legend.position="None")+
  coord_cartesian(x=as.Date(c("2007-06-01","2008-05-30")),ylim=c(0,1),expand = c(0.0))+
  labs(x="",y="Seroprevalence");Fig2
ggsave("Figure 2.pdf",plot=Fig2,width=16,height=12)
ggsave("Figure 2.png",plot=Fig2,width=16,height=12)

#Figure 3 is separate R file

#Figure 4: mom-pup MFIs
p1=q1[q1$age_class=="preweaned",c("unique_id","luminex","mother","antibody")] #subset to Faridpur pups
p1=rename(p1,pup_luminex=luminex,pup_id=unique_id)
m1=q1[!q1$age_class=="preweaned",] #subset to Faridpur moms
m1=merge(m1,p1,by.x=c("antibody","unique_id"),by.y=c("antibody","mother"))
m1$momlogMFI=log10(m1$luminex);m1$puplogMFI=log10(m1$pup_luminex)

#Figure 4 Mom - pup MFI
Fig4=ggplot(m1,aes(x=momlogMFI,y=puplogMFI,fill=antibody))+theme_few()+
  geom_point(size=4,shape=21,alpha=.8)+
  facet_wrap(~antibody,ncol=3)+
  scale_fill_viridis_d(option="C", begin=0.2, end=0.8)+
  theme(text=element_text(size=25),legend.position="None")+
  geom_abline(slope=1,intercept=0,linetype="dashed",linewidth=2)+
  labs(x=expression("Mother log"["10"]~"MFI"),
       y=expression("Pup log"["10"]~"MFI"));Fig4
ggsave("Figure 4.pdf",plot=Fig4,width=10,height=8)
ggsave("Figure 4.png",plot=Fig4,width=10,height=8)

#Figure 5 One year studies in Chakhoria and Ramnagar
q3=q1[!q1$age_class=="preweaned"&q1$location%in%c("chakhoria", "ramnagar"),] #subset
q3=q3 %>% arrange(date) %>% mutate(sess = cumsum(c(TRUE, diff(date) > 15)))
qB = data.frame(q3 %>% group_by(sess,age_class,antibody,location) %>% 
                  summarize(N=n(),prev=mean(Infect),Npos=sum(Infect),mdate=mean(date),jday=mean(jday)))
qB$ageloc=with(qB,paste(age_class,location,sep=", ") )
qB=qB[!(qB$sess==1&qB$age_class=="juvenile"),] #remove juvs from first session b/c from previous year

for (i in 1:nrow(qB)) {
  qB$L95[i]=with(qB,binom.test(Npos[i],N[i])$conf.int[1])
  qB$U95[i]=with(qB,binom.test(Npos[i],N[i])$conf.int[2]) }


qB$jday[qB$jday<152]=qB$jday[qB$jday<152]+365 #Make julian day starting in June

# g2=gam(prev~antibody+ageloc+s(jday),weights=N,gamma=3,data=qB,family="binomial");summary(g2)
# g3=gam(prev~antibody+ageloc+s(jday,by=interaction(antibody,ageloc)),weights=N,gamma=3,data=qB,family="binomial");summary(g3)
# AIC(g2,g3); AIC(g2)-AIC(g3)
# pdays=seq(range(qB$jday)[1],range(qB$jday)[2],by=1)
# ndChRa=expand.grid(jday=pdays,antibody=unique(qB$antibody),ageloc=unique(qB$ageloc))
# ndChRa$`interaction(antibody, ageloc)` <- interaction(ndChRa$antibody, ndChRa$ageloc)
# gamder2=cbind(ndChRa[,!names(ndChRa) %in% c("jday", "interaction(antibody, ageloc)")],
#               derivatives(g3,data=ndChRa))
# gamder2$pdate=as.Date(ymd("2010-01-01") + gamder2$jday) #yday("2007-06-01") #make wrapping date for plot
# gamder2$sig=ifelse(sign(gamder2$.lower_ci)==sign(gamder2$.upper_ci),"red","black")
# 
# gamder2$prev=predict(g3,newdata=gamder2,type="response")
# gamder2$prevlink=predict(g3,newdata=gamder2,type="link")
# gamder2$SE=predict(g3,newdata=gamder2,type="link",se.fit=T)$se
# gamder2$L95=plogis(gamder2$prevlink-1.96*gamder2$SE);gamder2$U95=plogis(gamder2$prevlink+1.96*gamder2$SE)
# gamder2$antibody=factor(gamder2$antibody,levels=c("Henipavirus","Filovirus","Rubulavirus"))

Fig5=ggplot(qB,aes(x=mdate,y=prev))+theme_few()+
  facet_grid(ageloc~antibody)+
  scale_y_continuous(labels=scales::percent)+
  # stat_smooth(aes(color=antibody,fill=antibody))+
  stat_smooth(aes(color=antibody,fill=antibody,weight=N),method="gam",formula = y ~ s(x, bs = "tp"),method.args = list(family = binomial, gamma = 1.25))+
  # geom_ribbon(data=gamder2,aes(x=pdate,y=prev,ymin=L95,ymax=U95,fill=antibody,color=antibody,alpha=.1))+
  # geom_path(data=gamder2,aes(x=pdate,y=prev),color=gamder2$sig)+
  geom_pointrange(aes(ymin=L95,ymax=U95),size=.9)+
  labs(x="Date",y="Prevalence")+
  scale_x_date(date_breaks="1 month",date_labels="%b",expand=c(0,0))+
  scale_color_viridis_d(option="C", begin=0.2, end=0.8) +
  scale_fill_viridis_d(option="C", begin=0.2, end=0.8)+
  theme(text=element_text(size=20),legend.position="None",
        axis.text.x=element_text(angle=45,vjust=0.25,hjust=0))+
  coord_cartesian(x=as.Date(c("2010-04-15","2011-05-14")),ylim=c(0,1),expand = T)+#c(0.001))+
  labs(x="",y="Seroprevalence");Fig5

ggsave("Figure 5.pdf",plot=Fig5,width=12,height=9)
ggsave("Figure 5.png",plot=Fig5,width=12,height=9)

#Figure S4
q1=read.csv("bat_serology_processed.csv")
q1$collection_date=as.Date(q1$collection_date)
q1$wk=as.numeric(strftime(q1$collection_date,format="%V")) #make week column
q1$yr=as.numeric(format(q1$collection_date,"%Y")) #make year column
q1$week=q1$wk+(q1$yr-2006)*52 #make running week, starting w/ 2006 same as old Rweek
q1=q1[!q1$age_class=="preweaned"&!q1$location%in%c("chakhoria", "ramnagar"),] #subset to Ruv,Ad+Juv,Faridpur
q1$Sero=0;q1$Sero[q1$seropos]=1 #make numerical seropos column
q1$antibody[q1$antibody=="Nipah"]="Henipavirus"
q1$antibody=factor(q1$antibody,levels=c("Henipavirus","Filovirus","Rubulavirus"))

#Make aggregated data frame to plot prev over time for adults & juveniles
qA = q1 %>% group_by(antibody,age_class,week,collection_date) %>% 
  summarize(count=n(),npos=sum(Sero),prev=mean(Sero),SE=(prev*(1-prev)/count)^.5 )
for (i in 1:nrow(qA)) { qA$L95[i] = binom.test(qA$count[i]*qA$prev[i],qA$count[i])$conf.int[1]
qA$U95[i]= binom.test(qA$count[i]*qA$prev[i],qA$count[i])$conf.int[2] }

ud=unique(qA$week)#unique time points
UA=unique(qA$antibody)
qA$Change=0;#qA$ChangeJ=0
for (i in 2:length(ud)) {
  for (j in 1:length(UA)) {
    datapair=qA[qA$antibody==UA[j]&qA$week%in%ud[i:(i-1)]&qA$age_class=="adult",]
    f1=glm(prev~week,family=binomial,data=datapair,weights=count)
    pval=summary(f1)$coefficients["week","Pr(>|z|)"]
    qA$Change[qA$antibody==UA[j]&qA$week==ud[i]&qA$age_class=="adult"]=
      sign(coef(f1)["week"])*ifelse(pval<0.05,1,0)  } }

qA$Change[qA$Change==-1]="Sig. decrease"
qA$Change[qA$Change==0]="No detec. change"
qA$Change[qA$Change==1]="Sig. increase"
qA$Change=relevel(factor(qA$Change),ref="Sig. decrease")
qA$collection_date_start=as.Date(c(NA,qA$collection_date[1:(nrow(qA)-1)]))
qA$prev_start=c(NA,qA$prev[1:(nrow(qA)-1)])
qA$prev_start[qA$week==82]=NA #remove segment for first week for all viruses

FigS4=ggplot(qA[qA$age_class=="adult",],aes(x=collection_date,y=prev))+theme_few()+
  facet_wrap(~antibody,nrow=3)+
  geom_segment(aes(xend = collection_date_start,yend = prev_start,color=antibody,
                   linetype = Change),linewidth=.65)+
  geom_pointrange(aes(ymin=L95,ymax=U95,shape=Change,fill=antibody),
                  position=position_dodge(width=10),size=1.5)+
  scale_fill_viridis_d(option="C", begin=0.2, end=0.8,guide = "none")+
  scale_color_viridis_d(option="C", begin=0.2, end=0.8,guide = "none")+
  scale_shape_manual(values=c(22,21,24))+
  scale_linetype_manual(values=c(2,3,1))+ #3,1,2
  scale_x_date(date_labels="%b-%y",date_breaks  ="6 month",
               limits=c(as.Date("2007-07-11"), as.Date("2012-12-25")),expand=c(0.02,0.01))+
  labs(x="Month-Year",y="Seroprevalence",color=NULL)+lims(y=c(0,1))+
  theme(text=element_text(size=18),legend.position="bottom",
        panel.grid.major.x =  element_line(colour = "gray"));FigS4

ggsave("Figure S4.pdf",plot=FigS4,width=12,height=8)
ggsave("Figure S4.png",plot=FigS4,width=12,height=8)

#Figure S5 Count data plot--------------------------
#Decline in counts over time
m=read.csv("niVFarLong3.csv");m$Date=as.Date(m$Date,"%m/%d/%Y") #serop data
m$WeekCC=m$Rweek-132 #weeks since counts were started
f1=lm(N~WeekCC,data=m);summary(f1)
#f1=lm(N~Date,data=m);summary(f1)
FigS5=ggplot(m,aes(x=Date,y=N))+theme_few()+
  geom_point(size=2)+
  stat_smooth(method="lm",color="black",linewidth=.5,alpha=.2)+
  scale_x_date(date_labels="%b-%y",date_breaks  ="6 month",
               limits=c(as.Date("2008-07-19"), as.Date("2012-12-25")),expand=c(0.02,0.01))+
  theme(text=element_text(size=18))+
  labs(x="",y="Colony size");FigS5
ggsave("Figure S5.pdf",plot=FigS5,width=12,height=8)

#--Fig S1--------------------------- From create_demo_plot
processed_data=read.csv("bat_serology_processed.csv")
processed_data$collection_date=as.Date(processed_data$collection_date)
bts <- processed_data |>
  filter(antibody == "Nipah", !location %in% c("chakhoria", "ramnagar")) |>
  mutate(repro = case_when(pregnant ~ "pregnant", has_pup | lactating ~ "lactating", TRUE ~ NA_character_)) |>
  mutate(
    status = if_else(!is.na(repro) & sex == "female" & age_class == "adult", paste(age_class, sex, repro, sep = ", "), paste(age_class, sex, sep = ", "), ),
    year = year(collection_date),
    age_class = fct_recode(age_class, Juveniles = "juvenile", Juveniles = "preweaned", Adults = "adult")
  ) |>
  mutate(mo = strftime(collection_date, format = "%b '%y")) |>
  arrange(mo, age_class, status)

bts_count <- bts |>
  mutate(status = fct_recode(status, 
                             "adult, female, pregnant/lactating" = "adult, female, lactating",
                             "adult, female, pregnant/lactating" = "adult, female, pregnant",
                             "juvenile" = "juvenile, male", "juvenile" = "juvenile, female",
                             "preweaned" = "preweaned, male", "preweaned" = "preweaned, female")) |> 
  count(mo, collection_date, age_class, status) |>
  arrange(collection_date) |> 
  mutate(n = if_else(age_class == "Juveniles", -n, n)) |> 
  mutate(status = forcats::fct_relevel(status, c(
    "adult, male", "adult, female", "adult, female, pregnant/lactating",
    "juvenile", "preweaned"
  ))) |> 
  mutate(mo = fct_reorder(mo, collection_date))

v <- viridis::viridis(5)

FigS1 <- ggplot(bts_count, aes(x = collection_date, y = n, fill = status)) + #formerly fig1_demo_plot
  geom_col(col = "black", lwd = 0.25) +
  scale_fill_manual(values = c(v[1], v[3], v[5], v[2], v[4])) +
  facet_grid(age_class~., scales = "free_y") +
  ylab("Count") +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%b '%y", name = "", limits = as.POSIXct(c("2007-01-01", "2013-03-01")), expand = c(0,0)) +
  # scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(labels = abs, breaks =c(-80, -60, -40, -20, 0, 20, 40, 60, 80), expand = c(0,0)) +
  geom_hline(yintercept = 0, size = 0.5, color = "black") +
  # noamtools::theme_nr +
  theme(
    text = element_text(color = "black"), 
    legend.position = "bottom", legend.background = element_rect(fill="white"),
    legend.title = element_blank(),
    strip.background = element_blank(), 
    strip.text = element_text(size = 15, margin = margin(b = 5), face = "bold"),
    axis.text = element_text(size = 8, color = "black"), 
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    panel.background = element_blank(), 
    panel.grid.major.x =  element_line(color="grey", size  = 0.25),
    panel.grid.minor.x =  element_line(color="grey", size  = 0.25),
    panel.grid.major.y = element_line(color="grey", size  = 0.25),
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_line(size = 0.25, color = "grey"),
    legend.text = element_text(size = 8, hjust = 0, margin = margin(l = -3, r = 10)), legend.key.size = unit(0.27, "cm"),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-15,-10,-10,-10),
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(15, 15, 15, 15, "pt"),
    panel.ontop = FALSE);FigS1

#---------------Figure S2
bts <- processed_data |>
  filter(antibody == "Nipah", location %in% c("chakhoria", "ramnagar")) |>
  mutate(repro = case_when(pregnant ~ "pregnant", has_pup | lactating ~ "lactating", TRUE ~ NA_character_)) |>
  mutate(
    status = if_else(!is.na(repro) & sex == "female" & age_class == "adult", paste(age_class, sex, repro, sep = ", "), paste(age_class, sex, sep = ", "), ),
    year = year(collection_date),
    age_class = fct_recode(age_class, Juveniles = "juvenile", Juveniles = "preweaned", Adults = "adult")
  ) |>
  mutate(mo = as.character(collection_date, format = "%Y-%b"),
         location = stringi::stri_trans_totitle(location)) |>
  arrange(mo, age_class, status)

bts_count <- bts |>
  count(mo, collection_date, age_class, status, location) |>
  arrange(collection_date) |>
  mutate(status = forcats::fct_relevel(status, c(
    "adult, female", "adult, female, pregnant", "adult, female, lactating", "adult, male",
    "preweaned, male", "preweaned, female", "juvenile, male", "juvenile, female"
  )))

v <- viridis::viridis(8)

FigS2 <- ggplot(bts_count, aes(x = forcats::fct_rev(as.factor(collection_date)), y = n, fill = status)) + #si_fig1_monthly_demo_plot
  geom_col(col = "black", lwd = 0.25) +
  coord_flip() +
  scale_fill_manual(values = c(v[1], v[3], v[5], v[7], v[2], v[4], v[6], v[8])) +
  facet_grid(location~age_class, scales = "free_y") +
  ylab("Count") +
  scale_x_discrete(labels = rev(unique(bts_count$mo))) +
  theme(
    text = element_text(color = "black"), 
    legend.position = "right", legend.background = element_rect(fill="white"),
    legend.title = element_blank(),
    strip.background = element_blank(), 
    strip.text = element_text(size = 15, margin = margin(b = 5), face = "bold"),
    strip.text.y = element_text(margin = margin(l = 8)),
    axis.text = element_text(size = 8, color = "black"), 
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    panel.background = element_blank(), panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color="grey", size  = 0.25),
    panel.grid.minor.x = element_blank(),
    axis.ticks = element_line(size = 0.25),
    legend.text = element_text(size = 10, hjust = 0, margin = margin(l = 5)), legend.key.size = unit(0.5, "cm"),
    panel.ontop = TRUE);FigS2


#--Figure 1 Venn diagram------------ from create_venn_plots()
sero_summary <- processed_data |>
  select(age_class, unique_id, antibody, seropos) |>
  spread(antibody, seropos) |>
  filter(!(is.na(Nipah) | is.na(Filovirus) | is.na(Rubulavirus))) |> 
  mutate(
    None = !Nipah & !Filovirus & !Rubulavirus,
    Nipah1 = Nipah & !Filovirus & !Rubulavirus,
    Filovirus1 = !Nipah & Filovirus & !Rubulavirus,
    Nipah_Filovirus = Nipah & Filovirus & !Rubulavirus,
    Rubulavirus1 = !Nipah & !Filovirus & Rubulavirus,
    Nipah_Rubulavirus = Nipah & !Filovirus & Rubulavirus,
    Filovirus_Rubulavirus = !Nipah & Filovirus & Rubulavirus,
    All = Nipah & Filovirus & Rubulavirus,
    N = 1,
  ) |>
  select(-unique_id, Nipah, Filovirus, Rubulavirus) |>
  group_by(age_class) |>
  summarize_all(sum)


sero_pct = sero_summary |>
  mutate(across(Nipah:All, ~ .x / N)) |>
  select(age_class, Nipah:All)

labels = structure(as.list(sero_summary$N), .Names = sero_summary$age_class)


x <- structure(unlist(sero_summary[1, 5:11]),
               .Names = c("A", "B", "A&B", "C", "A&C", "B&C", "A&B&C")
)
fit <- eulerr::euler(x)
a1 <- plot(fit,
           labels = c("Henipavirus", "Filovirus", "Rubulavirus"),
           fill = c("#9569BF", "#DA7FA0", "#FBC07B"),
           edges = "grey", quantities = list(type = "percent"), lty = 1, main = paste0("Adults (n = ", labels$adult, ")")
)

x <- structure(unlist(sero_summary[2, 5:11]),
               .Names = c("A", "B", "A&B", "C", "A&C", "B&C", "A&B&C")
)
fit <- eulerr::euler(x)
a2 <- plot(fit,
           labels = c("Henipavirus", "Filovirus", "Rubulavirus"),
           fill = c("#9569BF", "#DA7FA0", "#FBC07B"),
           edges = "grey", quantities = list(type = "percent"), lty = 1, main = paste0("Juveniles (n = ", labels$juvenile, ")")
)

x <- structure(unlist(sero_summary[3, 5:11]),
               .Names = c("A", "B", "A&B", "C", "A&C", "B&C", "A&B&C")
)
fit <- eulerr::euler(x)
a3 <- plot(fit,
           labels = c("Henipavirus", "Filovirus", "Rubulavirus"),
           fill = c("#9569BF", "#DA7FA0", "#FBC07B"),
           edges = "grey", quantities = list(type = "percent"), lty = 1, main = paste0("Pre-weaned (n = ", labels$preweaned, ")")
)

fig2_venn_plots <- gridExtra::arrangeGrob(a1, a2, a3, nrow = 1)
ggsave("Figure 1.png", fig2_venn_plots, width = 11.5, height = 3, dpi = 300, bg = "white")

#----Figure 1 analysis probit model
probit_data <- processed_data |>  #From prep_data_for_probit
  mutate(seropos = as.integer(seropos), 
         male = as.integer(sex == "male"),
         juvenile = as.integer(age_class == "juvenile"),
         preweaned = as.integer(age_class == "preweaned")
         #location = as.integer(as.factor(if_else(location %in% c("chakhoria", "ramnagar"), location, "faridpur"))),
  ) |> 
  select(unique_id, antibody, male, juvenile, preweaned, seropos, location) |> 
  spread(antibody, seropos) |> 
  filter(!(is.na(Nipah) | is.na(Filovirus) | is.na(Rubulavirus)))

probit_prepped_data <- lst(x = cbind(probit_data$juvenile, probit_data$preweaned, probit_data$male),
                           y = cbind(probit_data$Nipah, probit_data$Filovirus, probit_data$Rubulavirus),
                           #Z = as.matrix(probit_data$location),
                           #N_1 = ncol(r),
                           #M_1 = length(unique(r))
                           N = nrow(x),
                           K = ncol(x),
                           D = ncol(y))

probit_model_compiled = cmdstanr::cmdstan_model("multiprobit.stan")

probit_model_fit <- probit_model_compiled$sample(
  data = probit_prepped_data,
  seed = 1, iter_warmup = 1000, iter_sampling = 1000,
  chains = 4, parallel_chains = 4)

draws <- probit_model_fit$draws(variables = "Omega") |> 
  tidybayes::gather_draws(Omega[i, j])

probit_model_summary <- draws |> 
  group_by(i, j) |> 
  ggdist::mean_hdi(.value, .width = 0.95) |> 
  mutate_at(c("i","j"), \(x) c("Nipah", "Filovirus", "Rubulavirus")[x]) |> 
  rename(correlation = .value)

probit_model_summary

