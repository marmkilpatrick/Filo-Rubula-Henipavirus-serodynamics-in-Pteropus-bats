#Figures 2, 4, 5, S3
#This code creates the juvenile seroprevalence plot, Fig 2
lapply(c("tidyverse","ggthemes","viridis","brglm2"),require,character.only=T) #load packages
setwd("c:/marm/research/other/nipah/manyviruses/")
q1=read.csv("bat_serology_processed.csv")
q1$date=as.Date(q1$date)
q1$antibody[q1$antibody=="Nipah"]="Henipavirus"
q1$antibody=factor(q1$antibody,levels=c("Henipavirus","Filovirus","Rubulavirus"))
q1$Infect=0;q1$Infect[q1$seropos]=1 #make numerical seropos column
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

#Plot seasonal prev over time for juveniles
qA = data.frame(q2 %>% group_by(sess,yr,age_class,antibody) %>% 
                  summarize(N=n(),prev=mean(Infect),pdate=mean(pdate)))

for (i in 1:nrow(qA)) {
  qA$L95[i]=binom.test(qA$prev[i]*qA$N[i],qA$N[i])$conf.int[1]
  qA$U95[i]=binom.test(qA$prev[i]*qA$N[i],qA$N[i])$conf.int[2] }

Fig2=ggplot(qA[qA$age_class=="juvenile",],aes(x=pdate,y=prev))+theme_few()+
  facet_wrap(~antibody,nrow=3)+
  scale_y_continuous(labels=scales::percent)+
  stat_smooth(aes(color=antibody,fill=antibody))+
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

#Figure 4: mom-pup MFIs
p1=q1[q1$age_class=="preweaned",c("unique_id","luminex","mother","antibody")] #subset to Faridpur pups
p1=rename(p1,pup_luminex=luminex,pup_id=unique_id)
m1=q1[!q1$age_class=="preweaned",] #subset to Faridpur moms
m1=merge(m1,p1,by.x=c("antibody","unique_id"),by.y=c("antibody","mother"))
m1$momlogMFI=log10(m1$luminex);m1$puplogMFI=log10(m1$pup_luminex)

Fig4=ggplot(m1,aes(x=momlogMFI,y=puplogMFI,fill=antibody))+theme_few()+
  # geom_hline(data=diagCutOffs,aes(yintercept=LogCutOff),color="gray")+
  # geom_vline(data=diagCutOffs,aes(xintercept=LogCutOff),color="gray")+
  geom_point(size=4,shape=21,alpha=.8)+
  facet_wrap(~antibody,ncol=3)+
  scale_fill_viridis_d(option="C", begin=0.2, end=0.8)+
  theme(text=element_text(size=25),legend.position="None")+
  geom_abline(slope=1,intercept=0,linetype="dashed",linewidth=2)+
  labs(x=expression("Mother log"["10"]~"MFI"),
       y=expression("Pup log"["10"]~"MFI"));Fig4
ggsave("Figure 4.pdf",plot=Fig4,width=8,height=6)
ggsave("Figure 4.png",plot=Fig4,width=8,height=6)

#Figure 5 One year studies in Chakhoria and Ramnagar
q3=q1[!q1$age_class=="preweaned"&q1$location%in%c("chakhoria", "ramnagar"),] #subset
q3=q3 %>% arrange(date) %>% mutate(sess = cumsum(c(TRUE, diff(date) > 15)))

qB = data.frame(q3 %>% group_by(sess,age_class,antibody,location) %>% 
                  summarize(N=n(),prev=mean(Infect),Npos=sum(Infect),mdate=mean(date)))
qB$ageloc=with(qB,paste(age_class,location,sep=", ") )
for (i in 1:nrow(qB)) {
  qB$L95[i]=with(qB,binom.test(Npos[i],N[i])$conf.int[1])
  qB$U95[i]=with(qB,binom.test(Npos[i],N[i])$conf.int[2]) }

qB=qB[!(qB$sess==1&qB$age_class=="juvenile"),] #remove juvs from first session b/c from previous year
Fig5=ggplot(qB,aes(x=mdate,y=prev))+theme_few()+
  facet_grid(ageloc~antibody)+
  scale_y_continuous(labels=scales::percent)+
  stat_smooth(aes(color=antibody,fill=antibody))+
  geom_pointrange(aes(ymin=L95,ymax=U95),size=.9)+
  labs(x="Date",y="Prevalence")+
  scale_x_date(date_breaks="1 month",date_labels="%b",expand=c(0,0))+
  scale_color_viridis_d(option="C", begin=0.2, end=0.8) +
  scale_fill_viridis_d(option="C", begin=0.2, end=0.8)+
  theme(text=element_text(size=15),legend.position="None",
        axis.text.x=element_text(angle=45,vjust=0.25,hjust=0))+
  coord_cartesian(x=as.Date(c("2010-04-15","2011-05-14")),ylim=c(0,1),expand = c(0.001))+
  labs(x="",y="Seroprevalence");Fig5

ggsave("Figure 5.pdf",plot=Fig5,width=16,height=12)
ggsave("Figure 5.png",plot=Fig5,width=16,height=12)

#Figure S3
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

FigS3=ggplot(qA[qA$age_class=="adult",],aes(x=collection_date,y=prev))+theme_few()+
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
        panel.grid.major.x =  element_line(colour = "gray"));FigS3

ggsave("Figure S3.pdf",plot=FigS3,width=12,height=8)
ggsave("Figure S3.png",plot=FigS3,width=12,height=8)
