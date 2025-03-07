#Rubulavirus model with DDT transmission, recrudescence
#This is final model for paper
lapply(c("tidyverse","ggthemes","bbmle","beepr","viridis"),require,character.only=T) #load packages
setwd("c:/marm/research/other/nipah/manyviruses/")
q1=read.csv("bat_serology_processed.csv")
q1$collection_date=as.Date(q1$collection_date)
q1$wk=as.numeric(strftime(q1$collection_date,format="%V")) #make week column
q1$yr=as.numeric(format(q1$collection_date,"%Y")) #make year column
q1$week=q1$wk+(q1$yr-2006)*52 #make running week, starting w/ 2006 same as old Rweek
q1=q1[q1$antibody=="Rubulavirus"&!q1$age_class=="preweaned"&
        !q1$location%in%c("chakhoria", "ramnagar"),] #subset to RuV,Ad+Juv,Faridpur
q1$Age="Adult";q1$Age[q1$age_class=="juvenile"]="Juv"
q1$NiV=0;q1$NiV[q1$seropos]=1 #make numerical seropos column; used NiV to use code from NiV model

#Grouped data for prev over time for adults & juveniles
qA = q1 %>% group_by(week,age_class,collection_date) %>% 
  summarize(count=n(),npos=sum(NiV),prev=mean(NiV),SE=(prev*(1-prev)/count)^.5 )
for (i in 1:nrow(qA)) { qA$L95[i] = binom.test(qA$count[i]*qA$prev[i],qA$count[i])$conf.int[1]
  qA$U95[i]= binom.test(qA$count[i]*qA$prev[i],qA$count[i])$conf.int[2] }
qA$Age=qA$age_class #Age column used below

#model function
NiV= function(Bjj,Bja,Baj,Baa,Ias,r,SA_A,MAlr,LAA,LAJ,gt){ #Model which simulated dynamics given parameter estimates
  SA=plogis(SA_A)^(1/52);SJ=((1.02-SA^52)/.5)^(1/52) #SA is weekly survival; SJ is set to give correct decline in counts 1.02
  st=52;tend=st+max(qA$week);#starting and ending weeks - start on week 52 to create juvs from previous year
  A=rep(NA,tend);Rw=A;Rw[]=0;Nj=A;wc=A;wc[1:(st)]=1:52;nJ=Nj;#initializes variables
  #A - #adults, Nj-# juvs; nJ-# new juveniles born that week; wc-week
  Ast=265;A[1:(st+1)]=c(Ast*SA^(1:17),Ast*SA^17*SA^(1:11)*1.012^(1:11),Ast*1.012^11*SA^(29:(st+1))) #adult initial conditions for year before data
  IasL=plogis(Ias)#Initial antibody seroprevalence; using logit transform to constrain IasL to be 0 to 1
  Sa=(1-IasL-.005)*A;Ia=(.005)*A;Ra=(IasL)*A; # # of adult (a) susceptibles, infecteds, recovereds
#  Sa=max((1-IasL-.005),0)*A;Ia=(.005)*A;Ra=max((IasL-.005),0)*A; # constraining all to be +
  nJ[]=0;Nj[]=0;nJ[c(19:27)]=Ast*.5/9;
  Nji=c(90*SJ^(1:17),90*SJ^17+1:11*Ast*SA^17*.5/9*SJ^51);Nji=c(Nji,Nji[28]*SJ^(1:24));Nj[1:(st)]=Nji  #Juv initial conditions for year before data
  Mj=IasL*Nj
  Sj=Nj-Mj;Ij=0.005*Nj;Rj=0*Nj ## # of juveniles (j) susceptibles, infecteds, recovereds
  g=plogis(gt)#g=0.25;
  w=1#gamma = recovery rate = 1/ infectious period; w = week counter
  for (t in st:tend) {#loop over weeks of data
    if(w>18&w<28) {nJ[t]= A[t-9]*0.5*(1/9)} #births joining juvenile pop
    #    w/in age trans                     b/w age trans    recovery mortality    births            Loss of mat antib  Mature into adults              Recrudescence         Loss/gain immunity
    dSj=max(-exp(Bjj)*Sj[t]*Ij[t]-exp(Baj)*Sj[t]*Ia[t]        -(1-SJ)*Sj[t]+nJ[t]*Sa[t-5]/A[t-5]+plogis(MAlr)*Mj[t]     -nJ[t-51]*Sj[t]/Nj[t]*(SJ)^51                          +exp(LAJ)*Rj[t],
            -Sj[t])
    dIj=min(+exp(Bjj)*Sj[t]*Ij[t]+exp(Baj)*Sj[t]*Ia[t]        -(g-(1-SJ))*Ij[t]-(1-SJ)*Ij[t]                            -nJ[t-51]*Ij[t]/Nj[t]*(SJ)^51,
            Sj[t])
    dRj=+                                                       g*SJ*Ij[t]-(1-SJ)*Rj[t]                                 -nJ[t-51]*Rj[t]/Nj[t]*(SJ)^51                          -exp(LAJ)*Rj[t]
    dMj=                                                              -(1-SJ)*Mj[t]+nJ[t]*Ra[t-5]/A[t-5]-plogis(MAlr)*Mj[t]
    dSa=max(-exp(Baa)*Sa[t]*Ia[t]-exp(Bja)*Sa[t]*Ij[t]          -(1-SA)*Sa[t]+                                     nJ[t-51]*Sj[t]/Nj[t]*(SJ)^51                        +exp(LAA)*Ra[t],
            -Sa[t])
    dIa=min(+exp(Baa)*Sa[t]*Ia[t]+exp(Bja)*Sa[t]*Ij[t]  -(g-(1-SA))*Ia[t]-(1-SA)*Ia[t]+                                                                 exp(r)*Ra[t],
            Sa[t])
    dRa=+                                                      g*SA*Ia[t]-(1-SA)*Ra[t]+                                  nJ[t-51]*(Rj[t]+Ij[t])/Nj[t]*(SJ)^51-exp(r)*Ra[t]    -exp(LAA)*Ra[t]
    Sj[t+1]=Sj[t]+dSj;Ij[t+1]=Ij[t]+dIj;Rj[t+1]=Rj[t]+dRj;Mj[t+1]=Mj[t]+dMj
    Sa[t+1]=Sa[t]+dSa;Ia[t+1]=Ia[t]+dIa;Ra[t+1]=Ra[t]+dRa
    Nj[t+1]=Sj[t+1]+Ij[t+1]+Rj[t+1]+Mj[t+1];A[t+1]=Sa[t+1]+Ia[t+1]+Ra[t+1]
    wc[t]=w;w=w+1;Rw[t]=t-st;if(w==53) {w=1}} #Rw is week of actual data
    prev=data.frame(wc=c(NA,wc),w=c(Rw,tend-st+1),t=1:(tend+1),Aprev=Ra/A,Jprev=(Rj+Mj)/Nj,Sa,Ia,Ra,Sj,Ij,Rj)
  return(prev)
}

lkf=function(prev) {   #This function estimates the neg log likelihood using data in qA data.frame 
  ud=unique(qA$week)#unique time points
  lkA=matrix(0,(length(ud)-1),1);#likelihoods for Adult data
  lkJ=matrix(0,(length(ud)-1),1);#likelihoods for juvenile data
  for (i in 1:length(ud)) {
    rowA=qA$count[qA$week==ud[i]&qA$Age=="adult"]#total adults on week
    ppA=qA$npos[qA$week==ud[i]&qA$Age=="adult"]#observed # pos adult
    npA=rowA-ppA#observed # neg adult
    rowJ=qA$count[qA$week==ud[i]&qA$Age=="juvenile"]#total juvs on week
    ppJ=qA$npos[qA$week==ud[i]&qA$Age=="juvenile"]#observed # pos Juv
    npJ=rowJ-ppJ#observed # neg Juv
    lkA[i]=(prev$Aprev[prev$w==ud[i]]^ppA)*(1-prev$Aprev[prev$w==ud[i]])^npA #likelihood of adult data
    lkJ[i]=(prev$Jprev[prev$w==ud[i]]^ppJ)*(1-prev$Jprev[prev$w==ud[i]])^npJ #likelihood of juvenile data
  }
  tkl=-sum(log(lkA),log(lkJ)) #sum neg likelihood from both
  return(tkl)			}

##unpack parameters, call model function and call likelihood function 
ff2=function(Bjj,Bja,Baj,Baa,Ias,r,SA_A,MAlr,LAA,LAJ,gt) {#Function for fitting with mle2
  prev=  NiV(Bjj,Bja,Baj,Baa,Ias,r,SA_A,MAlr,LAA,LAJ,gt)
  B=       c(Bjj,Bja,Baj,Baa,Ias,r,SA_A,MAlr,LAA,LAJ,gt)
  lke=lkf(prev)
  print(c(B,lke))
  return(lke)  }

#-----------Fit model to data
#B=trans coeffjuv-juv juv->ad adu->juv  adu-adu   initial adult seroprevalence recrudescence ad surv   Mat. Antibod loss rate     Adult antibod loss rate
#                Bjj  Bja          Baj         Baa    Ias                                  r     SA_A      MAlr                        LAA

#This produces best fitting model and starting a bit away from MLEs provides SEs for all but 4 params
fp2 <- mle2(ff2,start=list(        #Fit model with starting values
  Bjj=-3.55,Bja=-4.28,Baj=-12.,Baa=-12.33,Ias=6.,
  r=-13,SA_A=.96,MAlr=-0.35,LAA=-11.1,LAJ=-4.,gt=8.),
  control=list(maxit=1e4)) #864.99
summary(fp2) #params + SE, stats, if hessian invertable + deviance (-2Log Likelihood)
873.5-864.99

#untransformed params
utp=data.frame(utv=c(exp(coef(fp2)[1:4]),plogis(coef(fp2)["Ias"]),exp(coef(fp2)["r"]),plogis(coef(fp2)["SA_A"]),
  plogis(coef(fp2)["MAlr"]),exp(coef(fp2)["LAA"]),exp(coef(fp2)["LAJ"]),plogis(coef(fp2)["gt"])));utp #param values untransformed

# #profiling - takes quite a while; individual profiles were assembled and saved: see below------------
# SEfp=summary(fp2)@coef[,"Std. Error"] #Define std.err vector for profiling
# SEfp[is.nan(SEfp)]=.5 #rough estimate to use for starting profiling
# #Profile parameters individually
# prof.Baj = profile(fitted=fp2,which="Baj",tol.newmin=0.1,try_harder=T,maxsteps=200,zmax=3,
#                     std.err=SEfp);beep(3) #tol.newmin - dev difference = 2*LL diff; Z = dev^.5
# write.csv(prof.Baj,file="prof.Baj.csv") 
# confint(prof.Baj, method="quad") #confidence intervals
# prof.r = profile(fitted=fp2,which="r",tol.newmin=0.1,try_harder=T,maxsteps=200,zmax=3,
#                    std.err=SEfp);beep(3) #tol.newmin - dev difference = 2*LL diff; Z = dev^.5
# write.csv(prof.r,file="prof.r.csv") #confidence intervals
# confint(prof.r, method="quad") #confidence intervals
# prof.SA_A = profile(fitted=fp2,which="SA_A",tol.newmin=0.1,try_harder=T,maxsteps=200,zmax=3,
#                  std.err=SEfp);beep(3) #tol.newmin - dev difference = 2*LL diff; Z = dev^.5
# write.csv(prof.SA_A,file="prof.SA_A.csv") #confidence intervals
# confint(prof.SA_A, method="quad") #confidence intervals
# prof.LAA = profile(fitted=fp2,which="LAA",tol.newmin=0.1,try_harder=T,maxsteps=200,zmax=3,
#                    std.err=SEfp);beep(3) #tol.newmin - dev difference = 2*LL diff; Z = dev^.5
# write.csv(prof.LAA,file="prof.LAA.csv") 
# confint(prof.LAA, method="quad") #confidence intervals
# plot(prof.LAA) #plot parameter profiles
# prof.LAJ = profile(fitted=fp2,which="LAJ",tol.newmin=0.1,try_harder=T,maxsteps=200,zmax=3,
#                   std.err=SEfp);beep(3) #tol.newmin - dev difference = 2*LL diff; Z = dev^.5
# write.csv(prof.LAJ,file="prof.LAJ.csv") 
# confint(prof.LAJ, method="quad") #confidence intervals
# plot(prof.LAJ) #plot parameter profiles
# prof.gt = profile(fitted=fp2,which="gt",tol.newmin=0.1,try_harder=T,maxsteps=200,zmax=3,
#                    std.err=SEfp);beep(3) #tol.newmin - dev difference = 2*LL diff; Z = dev^.5
# write.csv(prof.gt,file="prof.gt.csv") 
# confint(prof.gt, method="quad") #confidence intervals

#Simulate fitted model---------------------------------
prevM=do.call(NiV,as.list(coef(fp2))) #call model function w/ fitted params
prevM$iPrev_A=prevM$Ia/(prevM$Sa+prevM$Ia+prevM$Ra)
prevM$iPrev_J=prevM$Ij/(prevM$Sj+prevM$Ij+prevM$Rj)
prevL=pivot_longer(prevM,cols=c(Aprev,Jprev,iPrev_A,iPrev_J),names_to="age_class",values_to="prev") #pivot for plotting
prevL$type="Seroprevalence";prevL$type[grep("i",prevL$age_class)]="Infection Prevalence"
prevL$age_class=fct_recode(prevL$age_class,adult="Aprev",juvenile="Jprev",
                           adult="iPrev_A",juvenile="iPrev_J") #rename ages for plotting to match data
prevL$Rep=0;prevL$Rep[prevL$wc%in%c(19:27)]=1 #makes green reproduction rectangles when new juvs added to pop
prevL$collection_date=as.Date("2007-07-25")+(prevL$w-82)*7

ggplot(qA,aes(x=collection_date,y=prev))+geom_line(alpha=.1)+theme_few()+
  facet_wrap(~age_class,nrow=2)+
  geom_line(data=prevL,aes(x=collection_date,y=prev,color=type))+ #add fitted model
  geom_line(data=prevL,aes(x=collection_date,y=Rep),col="gray",alpha=.3)+ #add time of new juvs
  geom_pointrange(aes(ymin=L95,ymax=U95),position=position_dodge(width=10))+
  scale_color_manual(values=c(2,1))+
  scale_x_date(date_labels="%b-%y",date_breaks  ="6 month",
               limits = c(as.Date("2007-07-11"), as.Date("2012-11-25")),expand=c(0.01,0.01))+
  labs(x="Time",y="Seroprevalence or Infection prevalence",color="")+
  annotate(geom = "text", x=as.Date("2010-02-10"),y=.08,label="Infec prev",color="red")+
  theme(text=element_text(size=20))

#---------------- #create many simulations using stored parameter profiles
vcov_adj=vcov(fp2)
RUV_SE=read.csv("RUV_param_SE.csv") #load stored parameter profiles
for (i in 1:nrow(vcov_adj)) {vcov_adj[i,i]=RUV_SE$SE[i]^2} #Use stored SEs

prevLs=NULL;nsim=200;prevLm=prevL;s=1 #run 1k simulations to get ~50-100 sims w/in 30 logLik of MLE
for (i in 1:nsim) { #takes 2-3 min
  cfd=rnorm(length(diag(vcov_adj)),mean=coef(fp2),sd=diag(vcov_adj)^.5 )
  prevM=do.call(NiV,as.list(cfd));lkf(prevM)-fp2@min #call model function w/ fitted params
  if(!is.nan(lkf(prevM)) ) { if(lkf(prevM)-fp2@min<30) { #if likelihood is within X of MLE, include simulation
    prevM$iPrev_A=prevM$Ia/(prevM$Sa+prevM$Ia+prevM$Ra)
    prevM$iPrev_J=prevM$Ij/(prevM$Sj+prevM$Ij+prevM$Rj)
    #prevL=pivot_longer(prevM,cols=c(Aprev,Jprev),names_to="age_class",values_to="prev") #pivot for plotting
    prevLd=pivot_longer(prevM,cols=c(Aprev,Jprev,iPrev_A,iPrev_J),names_to="age_class",values_to="prev") #pivot for plotting
    prevLd$type="Seroprevalence";prevLd$type[grep("i",prevLd$age_class)]="Infection Prevalence"
    prevLd$age_class=fct_recode(prevLd$age_class,adult="Aprev",juvenile="Jprev",
                                adult="iPrev_A",juvenile="iPrev_J") #rename ages for plotting to match data
    prevLd$collection_date=as.Date("2007-07-25")+(prevLd$w-82)*7
    prevLd$sim=s;s=s+1
    prevLs=bind_rows(prevLs,prevLd)}}
};beep(3);s

#Store data frames since other code uses similar names
qA_RUV=qA
prevLm_RUV=prevLm
prevLs_RUV=prevLs

RUVp=ggplot(qA_RUV,aes(x=collection_date,y=prev))+theme_few()+
  facet_wrap(~age_class,nrow=2,strip.position="right")+
  geom_line(alpha=1,color=viridis(1,option="C", begin=0.8, end=0.8))+
  geom_line(data=prevL,aes(x=collection_date,y=Rep),col="gray",alpha=.5,linewidth=.5)+ #add time of new juvs
  geom_line(data=prevLm_RUV,aes(x=collection_date,y=prev,color=type),linewidth=1.25)+ #add fitted model  
  geom_pointrange(aes(ymin=L95,ymax=U95),position=position_dodge(width=10),
                  color=viridis(1,option="C", begin=0.8, end=0.8))+
  scale_color_manual(values=c(1,viridis(1,option="C", begin=0.8, end=0.8)))+
  scale_x_date(date_labels="%m/%y",date_breaks  ="6 month",
               limits = c(as.Date("2007-07-11"), as.Date("2012-11-25")),expand=c(0.01,0.01))+
  labs(x="",y="",color="")+
  ggtitle("B. Rubulavirus",)+
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size=20));RUVp

#Max infection prev in period w/ data
max(prevLm_RUV$prev[prevLm_RUV$type=="Infection Prevalence"&
                      prevLm_RUV$collection_date>as.Date("2007-07-11")])

for (i in 1:min(s,50)) { #add 50 simulations
  RUVp=RUVp+geom_line(data=prevLs_RUV[prevLs_RUV$sim==i,],aes(x=collection_date,y=prev,color=type),alpha=.2) };RUVp

