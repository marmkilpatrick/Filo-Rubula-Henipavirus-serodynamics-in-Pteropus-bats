#This is final model for paper
#Filovirus model with DDT transmission, loss of immunity
#This model doesn't include parameters for recrudescence because more complex model estimated it to be zero
lapply(c("tidyverse","ggthemes","bbmle","beepr","viridis"),require,character.only=T) #load packages
setwd("c:/marm/research/other/nipah/manyviruses/")
q1=read.csv("bat_serology_processed.csv")
q1$collection_date=as.Date(q1$collection_date)
q1$wk=as.numeric(strftime(q1$collection_date,format="%V")) #make week column
q1$yr=as.numeric(format(q1$collection_date,"%Y")) #make year column
q1$week=q1$wk+(q1$yr-2006)*52 #make running week, starting w/ 2006
q1=q1[q1$antibody=="Filovirus"&!q1$age_class=="preweaned"&
        !q1$location%in%c("chakhoria", "ramnagar"),] #subset to FiV,Ad+Juv,Faridpur
q1$Age="Adult";q1$Age[q1$age_class=="juvenile"]="Juv"
q1$NiV=0;q1$NiV[q1$seropos]=1 #make numerical seropos column
#Grouped data for prev over time for adults & juveniles
qA = q1 %>% group_by(antibody,week,age_class,collection_date) %>% 
  summarize(count=n(),prev=mean(NiV),SE=(prev*(1-prev)/count)^.5,
            pp=sum(NiV),np=count-pp)
#95% CIs for prev
for (i in 1:nrow(qA)) { qA$L95[i] = binom.test(qA$count[i]*qA$prev[i],qA$count[i])$conf.int[1]
  qA$U95[i]= binom.test(qA$count[i]*qA$prev[i],qA$count[i])$conf.int[2] }

#model function
NiV= function(Bjj,Bja,Baj,Baa,SA_A,LAA,LAJ){ #Model which simulated dynamics given parameter estimates
  SA=plogis(SA_A)^(1/52);SJ=((1.02-SA^52)/.5)^(1/52) #SA is weekly survival; SJ is set to give correct decline in counts 1.02
  st=52;tend=st+max(qA$week);#starting and ending weeks - start on week 52 to create juvs from previous year
  A=rep(NA,tend);Rw=A;Rw[]=0;Nj=A;wc=A;wc[1:(st)]=1:52;nJ=Nj;#initializes variables
  #A - #adults, Nj-# juvs; nJ-# new juveniles born that week; wc-week
  Ast=265;A[1:(st+1)]=c(Ast*SA^(1:17),Ast*SA^17*SA^(1:11)*1.012^(1:11),Ast*1.012^11*SA^(29:(st+1))) #adult initial conditions for year before data
  IasL=plogis(-11.3997644) #Using fit value from more complex model; essentially 0
  Sa=(1-IasL-.005)*A;Ia=(.005)*A;Ra=(IasL-.005)*A; # # of adult (a) susceptibles, infecteds, recovereds
#  Sa=max((1-IasL-.005),0)*A;Ia=(.005)*A;Ra=max((IasL-.005),0)*A; # constraining all to be +
  nJ[]=0;Nj[]=0;nJ[c(19:27)]=Ast*.5/9;
  Nji=c(90*SJ^(1:17),90*SJ^17+1:11*Ast*SA^17*.5/9*SJ^51);Nji=c(Nji,Nji[28]*SJ^(1:24));Nj[1:(st)]=Nji  #Juv initial conditions for year before data
  Mj=IasL*Nj
  Sj=Nj-Mj;Ij=0.005*Nj;Rj=0*Nj ## # of juveniles (j) susceptibles, infecteds, recovereds
  g=10.1032418;r=-28.5;MAlr=11.94
  w=1#gamma = recovery rate = 1/ infectious period; w = week counter
  for (t in st:tend) {#loop over weeks of data
    if(w>18&w<28) {nJ[t]= A[t-9]*0.5*(1/9)} #births joining juvenile pop
    #    w/in age trans                     b/w age trans    recovery mortality    births            Loss of mat antib  Mature into adults              Recrudescence         Loss/gain immunity
    dSj=max(-exp(Bjj)*Sj[t]*Ij[t]-exp(Baj)*Sj[t]*Ia[t]        -(1-SJ)*Sj[t]+nJ[t]*Sa[t-5]/A[t-5]+plogis(MAlr)*Mj[t]     -nJ[t-51]*Sj[t]/Nj[t]*(SJ)^51                                      +exp(LAJ)*Rj[t],
            -Sj[t])
    dIj=min(+exp(Bjj)*Sj[t]*Ij[t]+exp(Baj)*Sj[t]*Ia[t]-plogis(g)*Ij[t]-(1-SJ)*Ij[t]                            -nJ[t-51]*Ij[t]/Nj[t]*(SJ)^51,
            Sj[t])
    dRj=+                                                       plogis(g)*SJ*Ij[t]-(1-SJ)*Rj[t]                                 -nJ[t-51]*Rj[t]/Nj[t]*(SJ)^51                          -exp(LAJ)*Rj[t]
    dMj=                                                              -(1-SJ)*Mj[t]+nJ[t]*Ra[t-5]/A[t-5]-plogis(MAlr)*Mj[t]
    dSa=max(-exp(Baa)*Sa[t]*Ia[t]-exp(Bja)*Sa[t]*Ij[t]          -(1-SA)*Sa[t]+                                     nJ[t-51]*Sj[t]/Nj[t]*(SJ)^51                        +exp(LAA)*Ra[t],
            -Sa[t])
    dIa=min(+exp(Baa)*Sa[t]*Ia[t]+exp(Bja)*Sa[t]*Ij[t]  -plogis(g)*Ia[t]-(1-SA)*Ia[t]+                                                                 exp(r)*Ra[t],
            Sa[t])
    dRa=+                                                      plogis(g)*SA*Ia[t]-(1-SA)*Ra[t]+                                  nJ[t-51]*(Rj[t]+Ij[t])/Nj[t]*(SJ)^51-exp(r)*Ra[t]    -exp(LAA)*Ra[t]
    Sj[t+1]=Sj[t]+dSj;Ij[t+1]=Ij[t]+dIj;Rj[t+1]=Rj[t]+dRj;Mj[t+1]=Mj[t]+dMj
    Sa[t+1]=Sa[t]+dSa;Ia[t+1]=Ia[t]+dIa;Ra[t+1]=Ra[t]+dRa
    Nj[t+1]=Sj[t+1]+Ij[t+1]+Rj[t+1]+Mj[t+1];A[t+1]=Sa[t+1]+Ia[t+1]+Ra[t+1]
    wc[t]=w;w=w+1;Rw[t]=t-st;if(w==53) {w=1}} #Rw is week of actual data
    prev=data.frame(wc=c(NA,wc),w=c(Rw,tend-st+1),t=1:(tend+1),Aprev=Ra/A,Jprev=(Rj+Mj)/Nj,Sa,Ia,Ra,Sj,Ij,Rj,N=A+Nj)
  return(prev)
}

lkf=function(prev) {   #This function estimates the neg log likelihood using data in data.frame qA
  ud=unique(qA$week)#unique time points
  lkA=matrix(0,(length(ud)-1),1);#likelihoods for Adult data
  lkJ=matrix(0,(length(ud)-1),1);#likelihoods for juvenile data
  for (i in 1:length(ud)) {
    ppA=qA$pp[qA$week==ud[i]&qA$age_class=="adult"]
    npA=qA$np[qA$week==ud[i]&qA$age_class=="adult"]
    ppJ=qA$pp[qA$week==ud[i]&qA$age_class=="juvenile"]
    npJ=qA$np[qA$week==ud[i]&qA$age_class=="juvenile"]
    lkA[i]=(prev$Aprev[prev$w==ud[i]]^ppA)*(1-prev$Aprev[prev$w==ud[i]])^npA #likelihood of adult data
    lkJ[i]=(prev$Jprev[prev$w==ud[i]]^ppJ)*(1-prev$Jprev[prev$w==ud[i]])^npJ #likelihood of juvenile data
  }
  tkl=-sum(log(lkA),log(lkJ)) #sum likelihood from both
  return(tkl)			}

##unpack parameters, call model function and call likelihood function 
ff2=function(Bjj,Bja,Baj,Baa,SA_A,LAA,LAJ) {#Function for fitting with mle2
  prev=  NiV(Bjj,Bja,Baj,Baa,SA_A,LAA,LAJ)
  B=       c(Bjj,Bja,Baj,Baa,SA_A,LAA,LAJ)
  lke=lkf(prev)
  print(c(B,lke))
  return(lke)  }

#-----------Fit model to data
#B=trans coeffjuv-juv juv->ad adu->juv  adu-adu   initial adult seroprevalence recrudescence ad surv   Mat. Antibod loss rate     Adult antibod loss rate
#                Bjj  Bja          Baj         Baa    Ias                                  r     SA_A      MAlr                        LAA
fp2 <- mle2(ff2,start=list(#FiV starting/ending values
  Bjj=-3.81,Bja=-6.35,Baj=-7.98,Baa=-4.84,SA_A=1.26,
  LAA=-4.49,LAJ=-4.69)) #fit model 1134.86

1156.7 -1134.86 #Delta AIC FDT - DDT;

show(fp2)
coef(summary(fp2))
utp=data.frame(utv=c(exp(coef(fp2)[1:4]),plogis(coef(fp2)["SA_A"]),exp(coef(fp2)["LAA"]),
  exp(coef(fp2)["LAJ"])));utp #param values untransformed

# #profiling - takes 20-30 min - for saved values see below
# SEfp=summary(fp2)@coef[,"Std. Error"] #Define std.err vector for profiling
# SEfp[]=0.05 #starting scale for profiling
# pfit2 <- profile(fp2,tol.newmin=0.1,try_harder=T,maxsteps=200,zmax=3,
#                  std.err=SEfp);beep(3) #tol.newmin - dev difference = 2*LL diff; Z = dev^.5
# plot(pfit2,show.points=T) #plot parameter profiles
# FIV_CI=confint(pfit2, method="quad") #confidence intervals
# FIV_MLE=coef(fp2)
# #Store profiling results
# write.csv(FIV_MLE,file="FIV_MLE.csv") #MLE estimates
# write.csv(FIV_CI,file="FIV_CI.csv") #confidence intervals
# write.csv(pfit2,file="profile_FIV.csv") 
# # ------------end of profiling

#Simulate fitted model---------------------------------
prevM=do.call(NiV,as.list(coef(fp2))) #call model function w/ fitted params
prevM$iPrev_A=prevM$Ia/(prevM$Sa+prevM$Ia+prevM$Ra)
prevM$iPrev_J=prevM$Ij/(prevM$Sj+prevM$Ij+prevM$Rj)
prevL=pivot_longer(prevM,cols=c(Aprev,Jprev,iPrev_A,iPrev_J),names_to="age_class",values_to="prev") #pivot for plotting
prevL$type="Seroprevalence";prevL$type[grep("i",prevL$age_class)]="Infection Prevalence"
prevL$age_class=fct_recode(prevL$age_class,adult="Aprev",juvenile="Jprev",
                           adult="iPrev_A",juvenile="iPrev_J") #rename ages for plotting to match data
prevL$Rep=0;prevL$Rep[prevL$wc%in%c(19:27)]=1 #makes green reproduction rectangles when new juvs added to pop
prevL$collection_date=as.Date("2007-07-25")+(prevL$w-82)*7 #make date column for plotting

ggplot(qA,aes(x=collection_date,y=prev))+geom_line(alpha=.1)+theme_few()+
  facet_wrap(~age_class,nrow=2)+
  geom_line(data=prevL,aes(x=collection_date,y=prev,color=type))+ #add fitted model
  geom_line(data=prevL,aes(x=collection_date,y=Rep),col="gray",alpha=.6,linewidth=.5)+ #add time of new juvs
  geom_pointrange(aes(ymin=prev-SE,ymax=prev+SE))+scale_color_manual(values=c(2,1))+
  scale_x_date(date_labels="%b-%y",date_breaks  ="6 month",
               limits = c(as.Date("2007-07-11"), as.Date("2012-11-25")),expand=c(0.01,0.01))+
  labs(x="",y="Seroprevalence or Infection prevalence",color="")+
  annotate(geom = "text", x=as.Date("2009-12-20"),y=.08,label="Infec prev",color="red")+
  theme(text=element_text(size=20))

#------------------------------ #create many simulations
vcov_adj=vcov(fp2)
FIV_CI = read.csv("FIV_CI.csv") #read in stored profiling 95% CIs
for (i in 1:nrow(vcov_adj)) {vcov_adj[i,i]=((FIV_CI[FIV_CI$X==rownames(vcov_adj)[i],"X97.5.."]-
                                            FIV_CI[FIV_CI$X==rownames(vcov_adj)[i],"X2.5.." ])/3.92)^2}

prevLs=NULL;nsim=300;prevLm=prevL;s=1 #run 10k simulations to get ~50-100 sims w/in 30 logLik of MLE
for (i in 1:nsim) { #takes 1-2 minutes
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
qA_FIV=qA
prevLm_FIV=prevLm
prevLs_FIV=prevLs

FIVp=ggplot(qA_FIV,aes(x=collection_date,y=prev))+theme_few()+
  facet_wrap(~age_class,nrow=2,strip.position="right")+
  geom_line(alpha=1,color=viridis(1,option="C", begin=0.5, end=0.8))+
  geom_line(data=prevL,aes(x=collection_date,y=Rep),col="gray",alpha=.5,linewidth=.5)+ #add time of new juvs
  geom_line(data=prevLm_FIV,aes(x=collection_date,y=prev,color=type),linewidth=1.25)+ #add fitted model  
  geom_pointrange(aes(ymin=L95,ymax=U95),position=position_dodge(width=10),
                  color=viridis(1,option="C", begin=0.5, end=0.8))+
  scale_color_manual(values=c(1,viridis(1,option="C", begin=0.5, end=0.8)))+
  scale_x_date(date_labels="%m/%y",date_breaks  ="6 month",
               limits = c(as.Date("2007-07-11"), as.Date("2012-11-25")),expand=c(0.01,0.01))+
  labs(x="",y="Sero- or Infection prevalence",color="")+
  ggtitle("A. Filovirus",)+
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size=20));FIVp

max(prevLm_FIV$prev[prevLm_FIV$type=="Infection Prevalence"&
                      prevLm_FIV$collection_date>as.Date("2007-07-11")])

for (i in 1:min(s,50)) { #add simulations
  FIVp=FIVp+geom_line(data=prevLs_FIV[prevLs_FIV$sim==i,],aes(x=collection_date,y=prev,color=type),alpha=.1) };FIVp

