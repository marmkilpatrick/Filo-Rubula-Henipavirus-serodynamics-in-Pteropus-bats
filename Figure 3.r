#This code makes Figure 3 and requires code for RUV and FIV (see below)
#The code for NiV replots Fig 4 from Epstein et al 2020 PNAS; it uses fitted param estimates and 95CI
# to show model with params accounting for uncertainty; params.csv are from from model fitting code NIV_SIRR_DDT.R
lapply(c("tidyverse","ggthemes","beepr","ggpubr","viridis","patchwork"),require,character.only=T) #load packages
setwd("c:/marm/research/other/nipah/manyviruses/")
s=read.csv("params.csv") #parameter values from fitting
m=read.csv("niVFarLong3.csv");m$Date=as.Date(m$Date,"%m/%d/%Y") #serop data
m=rename(m,Npop=N)
m=dplyr::select(m,-c("Adult.neg","Juv.neg","Njuv","Nadult"))
m$t=m$Rweek+52 #make week that matches model column t
mL=pivot_longer(m, cols  = c(Adult.tot,Juv.tot,Adult.pos,Juv.pos), 
                names_to = c("age_class",".value"), names_sep = "\\.")
mL$prev=mL$pos/mL$tot
for (i in 1:nrow(mL)) { mL$L95[i] = binom.test(mL$pos[i],mL$tot[i])$conf.int[1]
  mL$U95[i]= binom.test(mL$pos[i],mL$tot[i])$conf.int[2] }
mL$age_class=case_match(mL$age_class,"Adult"~"adult","Juv"~"juvenile")

pv=data.frame(matrix(NA,nrow=1,ncol=9));names(pv)=s$Param;pv[]=s$MLE;pvMLE=pv
st=52;tend=st+max(m$Rweek)+10
source("NIVmodelfunction.r") #load NiV function 
dm=with(pv,NiV(Bjj,Bja,Baj,Baa,Ias,r,SA_A,MAlr,LAA))
dm$Date=as.Date("2006-02-03") +7*(dm$t-57) #add date column for plotting
dm$A.Iprev=with(dm,Ia/A);dm$J.Iprev=with(dm,Ij/Nj)
dm=rename(dm,A.Sprev=Aprev,J.Sprev=Jprev)

dmL=data.frame(pivot_longer(dm, cols  = c(A.Sprev, J.Sprev, A.Iprev, J.Iprev), 
                names_to = c("age_class",".value"), names_sep = "\\.") )

dmL2=data.frame(pivot_longer(dmL,col=c(Sprev,Iprev),names_to="type",values_to="prev"))
dmL2$age_class=case_match(dmL2$age_class,"A"~"adult","J"~"juvenile")
dmL2$type=case_match(dmL2$type,"Sprev"~"Seroprevalence","Iprev"~"Infection Prevalence")
dmL2$Rep=0;dmL2$Rep[dmL2$wc%in%c(19:27)]=1 #when new juvs added to pop

dms=NULL;set.seed(3)
for (i in 1:50) {
  pv[]=rnorm(length(pv),s$MLE,(s$upper95ci-s$MLE)/2)
  dmT=with(pv,NiV(Bjj,Bja,Baj,Baa,Ias,r,SA_A,MAlr,LAA))
  dmT$nsim=i
  dms=bind_rows(dms,dmT)
}

dms$Date=as.Date("2006-02-03") +7*(dms$t-57) #add date column for plotting
dms$A.Iprev=with(dms,Ia/A);dms$J.Iprev=with(dms,Ij/Nj)
dms=rename(dms,A.Sprev=Aprev,J.Sprev=Jprev)

dmLs=data.frame(pivot_longer(dms, cols  = c(A.Sprev, J.Sprev, A.Iprev, J.Iprev), 
                             names_to = c("age_class",".value"), names_sep = "\\.") )

dmL2s=data.frame(pivot_longer(dmLs,col=c(Sprev,Iprev),names_to="type",values_to="prev"))
dmL2s$age_class=case_match(dmL2s$age_class,"A"~"adult","J"~"juvenile")
dmL2s$type=case_match(dmL2s$type,"Sprev"~"Seroprevalence","Iprev"~"Infection Prevalence")

NIVpleg=ggplot(mL,aes(x=Date,y=prev))+theme_few()+
  facet_wrap(~age_class,nrow=2,strip.position="right")+
  geom_pointrange(aes(ymin=L95,ymax=U95),position=position_dodge(width=10),
                  color=viridis(1,option="C", begin=0.2, end=0.8))+geom_line()+
  geom_line(data=dmL2,aes(x=Date,y=prev,color=type),linewidth=1.25)+ #fitted model
  geom_line(data=dmL2,aes(x=Date,y=Rep),col="gray",alpha=.5,linewidth=.5)+ #add time of new juvs
  scale_color_manual(values=c(1,viridis(1,option="C", begin=0.2, end=0.8)))+
  scale_x_date(date_labels="%m/%y",date_breaks  ="6 month",
               limits = c(as.Date("2007-07-11"), as.Date("2012-11-25")),expand=c(0.01,0.01))+
  ggtitle("C. Henipavirus")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="",y="Sero- or Infection prevalence",color="")+
  theme(text=element_text(size=20));NIVpleg

max(dmL2$prev[dmL2$type=="Infection Prevalence"&dmL2$Date>as.Date("2007-07-11")])

for (i in 1:50) { #add simulations
  NIVpleg=NIVpleg+geom_line(data=dmL2s[dmL2s$nsim==i,],aes(x=Date,y=prev,color=type),alpha=0.1) };NIVpleg

#Use Patchwork to make multipanel fig w legend
design = "111222
          3334##"
#Need to run these scripts to create FIV and RUV panels before running next line: 
#FIV_Final.r
#RUV_Final.r

#Then:
fig3=FIVp+RUVp+NIVpleg+guide_area()+plot_layout(design=design, guides = "collect");fig3
ggsave("Figure 3.pdf",plot=fig3,width=16,height=12) 
ggsave("Figure 3.png",plot=fig3,width=16,height=12) 

