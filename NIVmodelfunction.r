NiV= function(Bjj,Bja,Baj,Baa,Ias,r,SA_A,MAlr,LAA){ #Model which simulated dynamics given parameter estimates
  SA=SA_A^(1/52);SJ=((1.02-SA^52)/.5)^(1/52) #SA is weekly survival; SJ is set to give correct decline in counts 1.02
  st=52;tend=st+max(m$Rweek)+10;#starting and ending weeks - start on week 52 to create juvs from previous year
  A=rep(NA,tend);Rw=A;Rw[]=0;Nj=A;wc=A;wc[1:(st)]=1:52;nJ=Nj;#initializes variables
  #A - #adults, Nj-# juvs; nJ-# new juveniles born that week; wc-week
  Ast=265;A[1:(st+1)]=c(Ast*SA^(1:17),Ast*SA^17*SA^(1:11)*1.012^(1:11),Ast*1.012^11*SA^(29:(st+1))) #adult initial conditions for year before data
  IasL=1/(1+exp(-Ias))#Initial antibody seroprevalence; using logit transform to constrain IasL to be 0 to 1
  Sa=(1-IasL-.005)*A;Ia=(.005)*A;Ra=(IasL-.005)*A; # # of adult (a) susceptibles, infecteds, recovereds
  nJ[]=0;Nj[]=0;nJ[c(18:28)]=Ast*.5/9;
  Nji=c(90*SJ^(1:17),90*SJ^17+1:11*Ast*SA^17*.5/9*SJ^51);Nji=c(Nji,Nji[28]*SJ^(1:24));Nj[1:(st)]=Nji  #Juv initial conditions for year before data
  Mj=IasL*Nj
  Sj=Nj-Mj;Ij=0.005*Nj;Rj=0*Nj ## # of juveniles (j) susceptibles, infecteds, recovereds
  g=1;w=1#gamma = recovery rate = 1/ infectious period; w = week counter
  for (t in st:tend) {#loop over weeks of data
    if(w>18&w<28) {nJ[t]= A[t-9]*0.5*(1/9)} #births joining juvenile pop
    #    w/in age trans                     b/w age trans    recovery mortality    births            Loss of mat antib  Mature into adults              Recrudescence
    dSj=-exp(Bjj)*Sj[t]*Ij[t]-exp(Baj)*Sj[t]*Ia[t]        -(1-SJ)*Sj[t]+nJ[t]*Sa[t-5]/A[t-5]+MAlr*Mj[t]     -nJ[t-51]*Sj[t]/Nj[t]*(SJ)^51
    dIj=+exp(Bjj)*Sj[t]*Ij[t]+exp(Baj)*Sj[t]*Ia[t]-(g-(1-SA))*Ij[t]-(1-SJ)*Ij[t]                            -nJ[t-51]*Ij[t]/Nj[t]*(SJ)^51
    dRj=+                                                       g*SJ*Ij[t]-(1-SJ)*Rj[t]                                 -nJ[t-51]*Rj[t]/Nj[t]*(SJ)^51
    dMj=                                                              -(1-SJ)*Mj[t]+nJ[t]*Ra[t-5]/A[t-5]-MAlr*Mj[t]
    dSa=max(-exp(Baa)*Sa[t]*Ia[t]-exp(Bja)*Sa[t]*Ij[t],-Sa[t])          -(1-SA)*Sa[t]+                                     nJ[t-51]*Sj[t]/Nj[t]*(SJ)^51+exp(LAA)*Ra[t]
    dIa=min(+exp(Baa)*Sa[t]*Ia[t]+exp(Bja)*Sa[t]*Ij[t],Sa[t])  -(g-(1-SA))*Ia[t]-(1-SA)*Ia[t]+                                                                 exp(r)*Ra[t]
    dRa=+                                                      g*SA*Ia[t]-(1-SA)*Ra[t]+                                  nJ[t-51]*(Rj[t]+Ij[t])/Nj[t]*(SJ)^51-exp(r)*Ra[t]-exp(LAA)*Ra[t]
    Sj[t+1]=Sj[t]+dSj;Ij[t+1]=Ij[t]+dIj;Rj[t+1]=Rj[t]+dRj;Mj[t+1]=Mj[t]+dMj
    Sa[t+1]=Sa[t]+dSa;Ia[t+1]=Ia[t]+dIa;Ra[t+1]=Ra[t]+dRa
    Nj[t+1]=Sj[t+1]+Ij[t+1]+Rj[t+1]+Mj[t+1];A[t+1]=Sa[t+1]+Ia[t+1]+Ra[t+1]
    wc[t]=w;w=w+1;Rw[t]=t-st;if(w==53) {w=1}} #Rw is week of actual data
  #   prev=data.frame(w=c(Rw,tend-st+1),Aprev=Ra/A,Jprev=(Rj+Mj)/Nj)
  prev=data.frame(wc=c(wc,NA),w=c(Rw,tend-st+1),t=1:(tend+1),Aprev=Ra/A,Jprev=(Rj+Mj)/Nj,Sa,Ia,Ra,Sj,Ij,Rj,Mj,Nj,A)
  return(prev)
}
