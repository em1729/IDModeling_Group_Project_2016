library (deSolve)

#SEIR shared host model for Fair A attendees including both age groups

SEIR2host=function(t,state,parameters){
  with(as.list(c(state,parameters)),{

    dSc = -cSH * betaSH * (Sc)* (Is/Ns)
    dEc = cSH * betaSH * (Sc) * (Is/Ns)  - alpha * Ec 
    dIc = alpha * Ec - gamma * Ic
    
    
    dSa = -cSH * (Sa) * betaSH * (Is/Ns)
    dEa = cSH * (Sa) * betaSH * (Is/Ns) - alpha * Ea
    dIa = alpha * Ea - gamma * Ia
    
    dSs = -betaSS * Ss *  (Is/Ns)
    dIs = betaSS * Ss  * (Is/Ns) - gamma * Is
   
        list(c(dSc, dEc, dIc, dSa, dEa, dIa, dSs, dIs))#put in order by equation
  })
}


# Inital states and parameters provided by the paper

times=seq(0,9) 

#children parameters

Nc=6468; Ec=0; Ic=0; Sc=Nc*.9-Ec-Ic;  # IMc = 0.1 preexisting immunity to H3N2v among children


#adult parameters

Na=8442; Ea=0; Ia=0; Sa=Na*.5-Ea-Ia; # IMa = 0.5 preeisitng immunity to H3N2v among adults 


#swine parameters

Ns= 208; 
Is=1; 
Ss=Ns-Is;  # number of suseptible swine   
RoSS=2;  # RoSS = swine to swine H3N2pM reproductive number
gamma=1/5;   # gamma = swine and human recovery rate 
betaSS= RoSS * gamma  # beta!! transmission! from swine to swine the paper says its force of infection but its really the transmission rate



state=c(Sc=Sc, Ec=Ec, Ic=Ic, Sa=Sa, Ea=Ea, Ia=Ia, Ss=Ss, Is=Is); #same order as above in the list

parameters=c(betaSH=.024,alpha=1/2,gamma,  
             cSH=5, betaSS);

# betaSH = transmission probabilty of H3N2v, swine to human
# cSH = mean duration of contact between humans and swine for all general attendees at Fair A,  (5min/d)
# alpha = rate of exposed humans to infected human = 1/incubation period 



#simulation
simSEIR2host=ode(y=state,times=times,func=SEIR2host,parms=parameters) 

View(simSEIR2host)

SC=simSEIR2host[,'Sc']/Nc
EC=simSEIR2host[,'Ec']/Nc
IC=simSEIR2host[,'Ic']/Nc
SA=simSEIR2host[,'Sa']/Na
EA=simSEIR2host[,'Ea']/Na
IA=simSEIR2host[,'Ia']/Na
SS=simSEIR2host[,'Ss']/Ns
IS=simSEIR2host[,'Is']/Ns

SC
EC
IC
SA
EA
IA
SS
IS

#calculating incidence cases among younger age group

# Element at first row (record # at time 0 days)and second column (Sc)
simSEIR2host[1,2] # 5821.2
# Element at 10th row (record # at time 9 days)and second column (Sc)
simSEIR2host[10,2] # 5740.012
## Total number of infections among person aged <20 years at the fair
(simSEIR2host[1,2]) - (simSEIR2host[10,2])  #81.18825  # the paper got 80 

(81.18825/6468)* 100 # = 1.25 % Paper got 1.25%

#calculating incidence cases among older age group

# Element at first row (record # at time 0 days)and second column (Sa)
simSEIR2host[1,5] # 4221
# Element at 10th row (record # at time 9 days)and second column (Sa)
simSEIR2host[10,5] # 4174.445
## Total number of infections among person aged >=20 years at the fair
(simSEIR2host[1,5]) - (simSEIR2host[10,5]) #58.87027  # the paper got 58

(58.87027/8442) * 100  # = 0.69% Paper got 0.55% 

#calculating incidence cases among swine

# Element at first row (record # at time 0 days)and second column (Sa)
simSEIR2host[1,8] # 207
# Element at 10th row (record # at time 9 days)and second column (Sa)
simSEIR2host[10,8] # 197.5322
## Total number of infections among person aged >=20 years at the fair
(simSEIR2host[1,8]) - (simSEIR2host[10,8]) # 9.467



#Graphs

#To plot simulated proportion incident influenza A (H3N2) virus infections among humans and swine 

plot(simSEIR2host[,'time'],IC, ylab='Incidence',xlab='Duration of Fair (days)', ylim = c(0,0.035), type='b',lwd=2,  col="blue", main = "Proportion of Infectious Incidence of Swine Flu at Fair A", cex.main=0.85, font.main=14)
lines(simSEIR2host[,'time'],IA,type='b',lwd=2,col='olivedrab4')
lines(simSEIR2host[,'time'],IS,type='b',lwd=2,col='deeppink')
legend('topleft',legend = (c('Humans<20yrs','Humans>=20yrs', 'Swine')),lty=1, lwd=2, col=c('blue', 'olivedrab4', 'deeppink'), bty='n', cex=.95)



#To plot simulated # incident cases of influenza A (H3N2) virus infections among humans and swine 

plot(simSEIR2host[,'time'],simSEIR2host[,'Ic'], ylab='No. of Incidence Cases',xlab='Duration of Fair (days)', ylim = c(0,42), type='b',lwd=2,  col="blue", main = "# of Infectious cases of Swine Flu at Fair A",cex.main=0.85, font.main=14)
lines(simSEIR2host[,'time'],simSEIR2host[,'Ia'],type='b',lwd=2,col='olivedrab4')
lines(simSEIR2host[,'time'],simSEIR2host[,'Is'],type='b',lwd=2,col='deeppink')
legend('topleft',legend = (c('Humans<20yrs','Humans>=20yrs', 'Swine')),lty=1, lwd=2, col=c('blue', 'olivedrab4', 'deeppink'), bty='n', cex=.95)

###########################################################################################
##########################################################################################
library (deSolve)
#Sensitivity Analysis

SEIR2host0.75=function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    dSc = -cSH * betaSH * (Sc)* (Is/Ns)
    dEc = cSH * betaSH * (Sc) * (Is/Ns)  - alpha * Ec 
    dIc = alpha * Ec - gamma * Ic
    
    dSa = -cSH * (Sa) * betaSH * (Is/Ns)
    dEa = cSH * (Sa) * betaSH * (Is/Ns) - alpha * Ea
    dIa = alpha * Ea - gamma * Ia
    
    dSs = -betaSS * Ss *  (Is/Ns)
    dIs = betaSS * Ss  * (Is/Ns) - gamma * Is
    
    list(c(dSc, dEc, dIc, dSa, dEa, dIa, dSs, dIs))#put in order by equation
  })
}


# Inital states and parameters provided by the paper

times=seq(0,9) 

#children parameters

Nc=6468; Ec=0; Ic=0; Sc=Nc*.9-Ec-Ic;  # IMc = 0.1 preexisting immunity to H3N2v among children


#adult parameters

Na=8442; Ea=0; Ia=0; Sa=Na*.5-Ea-Ia; # IMa = 0.5 preeisitng immunity to H3N2v among adults 


#swine parameters

Ns= 208; 
Is=1; 
Ss=Ns-Is;  # number of suseptible swine   
RoSS=2;  # RoSS = swine to swine H3N2pM reproductive number
gamma=1/5;   # gamma = swine and human recovery rate 
betaSS= RoSS * gamma  # beta!! transmission! from swine to swine the paper says its force of infection but its really the transmission rate



state=c(Sc=Sc, Ec=Ec, Ic=Ic, Sa=Sa, Ea=Ea, Ia=Ia, Ss=Ss, Is=Is); #same order as above in the list

parameters=c(betaSH=.017,alpha=1/2,gamma,  
             cSH=5, betaSS);

# betaSH = transmission probabilty of H3N2v, swine to human
# cSH = mean duration of contact between humans and swine for all general attendees at Fair A,  (5min/d)
# alpha = rate of exposed humans to infected human = 1/incubation period 



#simulation
simSEIR2host0.75=ode(y=state,times=times,func=SEIR2host0.75,parms=parameters) 

View(simSEIR2host0.75)

SC=simSEIR2host0.75[,'Sc']/Nc
EC=simSEIR2host0.75[,'Ec']/Nc
IC=simSEIR2host0.75[,'Ic']/Nc
SA=simSEIR2host0.75[,'Sa']/Na
EA=ssimSEIR2host0.75[,'Ea']/Na
IA=simSEIR2host0.75[,'Ia']/Na
SS=simSEIR2host0.75[,'Ss']/Ns
IS=simSEIR2host0.75[,'Is']/Ns

SC
EC
IC
SA
EA
IA
SS
IS


#calculating incidence cases among younger age group

# Element at first row (record # at time 0 days)and second column (Sc)
simSEIR2host0.75[1,2] 
# Element at 10th row (record # at time 9 days)and second column (Sc)
simSEIR2host0.75[10,2] 
## Total number of infections among person aged <20 years at the fair
 (simSEIR2host0.75[1,2]) - (simSEIR2host0.75[10,2]) # =57.62
57.62602/6468 *100  # = 0.89%  paper got 0.9%


#calculating incidence cases among older age group

# Element at first row (record # at time 0 days)and second column (Sa)
simSEIR2host0.75[1,5] 
# Element at 10th row (record # at time 9 days)and second column (Sa)
simSEIR2host0.75[10,5] 
## Total number of infections among person aged >=20 years at the fair
(simSEIR2host0.75[1,5]) - (simSEIR2host0.75[10,5]) # = 41.78
41.7851/8442 * 100 # = 0.49% same as paper(0.5%)

#Graphs

#To plot simulated proportion incident influenza A (H3N2) virus infections among humans and swine 

plot(simSEIR2host0.75[,'time'],IC, ylab='Incidence',xlab='Duration of Fair (days)', ylim = c(0,0.035), type='b',lwd=2,  col="blue", main = "Proportion of Infectious Incidence of Swine Flu at Fair A", sub = "Sensitivity Analysis (Beta=0.017)", cex.main=0.85, font.main=14, font.sub=3)
lines(simSEIR2host0.75[,'time'],IA,type='b',lwd=2,col='olivedrab4')
lines(simSEIR2host0.75[,'time'],IS,type='b',lwd=2,col='deeppink')
legend('topleft',legend = (c('Humans<20yrs','Humans>=20yrs', 'Swine')),lty=1, lwd=2, col=c('blue', 'olivedrab4', 'deeppink'), bty='n', cex=.95)



#To plot simulated # incident cases of influenza A (H3N2) virus infections among humans and swine 

plot(simSEIR2host0.75[,'time'],simSEIR2host0.75[,'Ic'], ylab='No. of Incidence Cases',xlab='Duration of Fair (days)', ylim = c(0,42), type='b',lwd=2,  col="blue", main = "# of Infectious cases of Swine Flu at Fair A",cex.main=0.85, font.main=14, sub = "Sensitivity Analysis (Beta=0.017)", font.sub=3)
lines(simSEIR2host0.75[,'time'],simSEIR2host0.75[,'Ia'],type='b',lwd=2,col='olivedrab4')
lines(simSEIR2host0.75[,'time'],simSEIR2host0.75[,'Is'],type='b',lwd=2,col='deeppink')
legend('topleft',legend = (c('Humans<20yrs','Humans>=20yrs', 'Swine')),lty=1, lwd=2, col=c('blue', 'olivedrab4', 'deeppink'), bty='n', cex=.95)



#######################################################################################
#######################################################################################

#Adding higher swine contact for <20 age group

SEIR2hostCSHc=function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    dSc = -cSHc * betaSH * (Sc)* (Is/Ns)
    dEc = cSHc * betaSH * (Sc) * (Is/Ns)  - alpha * Ec 
    dIc = alpha * Ec - gamma * Ic
    
    dSa = -cSHa * (Sa) * betaSH * (Is/Ns)
    dEa = cSHa * (Sa) * betaSH * (Is/Ns) - alpha * Ea
    dIa = alpha * Ea - gamma * Ia
    
    dSs = -betaSS * Ss *  (Is/Ns)
    dIs = betaSS * Ss  * (Is/Ns) - gamma * Is
    
    list(c(dSc, dEc, dIc, dSa, dEa, dIa, dSs, dIs))#put in order by equation
  })
}


# Inital states and parameters provided by the paper

times=seq(0,9) 

#children parameters

Nc=6468; Ec=0; Ic=0; Sc=Nc*.9-Ec-Ic;  # IMc = 0.1 preexisting immunity to H3N2v among children


#adult parameters

Na=8442; Ea=0; Ia=0; Sa=Na*.5-Ea-Ia; # IMa = 0.5 preeisitng immunity to H3N2v among adults 


#swine parameters

Ns= 208; 
Is=1; 
Ss=Ns-Is;  # number of suseptible swine   
RoSS=2;  # RoSS = swine to swine H3N2pM reproductive number
gamma=1/5;   # gamma = swine and human recovery rate 
betaSS= RoSS * gamma  # beta!! transmission! from swine to swine the paper says its force of infection but its really the transmission rate



state=c(Sc=Sc, Ec=Ec, Ic=Ic, Sa=Sa, Ea=Ea, Ia=Ia, Ss=Ss, Is=Is); #same order as above in the list

parameters=c(betaSH=.024,alpha=1/2,gamma,  
             cSHa=5, cSHc=8, betaSS);

# betaSH = transmission probabilty of H3N2v, swine to human
# cSHa = mean duration of contact between humans and swine for adult attendees at Fair A,  (5min/d)
# cSHc = mean duration of contact between humans and swine for children attendees at Fair A,  
          #(8min/d) # 8 minutes seems like a good number. Tried 10 but # of cases seemed to high
# alpha = rate of exposed humans to infected human = 1/incubation period 



#simulation
simSEIR2hostcSHc=ode(y=state,times=times,func=SEIR2hostCSHc,parms=parameters) 

View(simSEIR2hostcSHc)

SC=simSEIR2hostcSHc[,'Sc']/Nc
EC=simSEIR2hostcSHc[,'Ec']/Nc
IC=simSEIR2hostcSHc[,'Ic']/Nc
SA=simSEIR2hostcSHc[,'Sa']/Na
EA=simSEIR2hostcSHc[,'Ea']/Na
IA=simSEIR2hostcSHc[,'Ia']/Na
SS=simSEIR2hostcSHc[,'Ss']/Ns
IS=simSEIR2hostcSHc[,'Is']/Ns

SC
EC
IC
SA
EA
IA
SS
IS


#calculating incidence cases among younger age group

# Element at first row (record # at time 0 days)and second column (Sc)
simSEIR2hostcSHc[1,2] # 5821.2
# Element at twenty-first row (record # at time 9 days)and second column (Sc)
simSEIR2hostcSHc[10,2] # 5691.843
## Total number of infections among person aged <20 years at the fair

(simSEIR2hostcSHc[1,2]) - (simSEIR2hostcSHc[10,2]) # = 129.35 cases

(129.35/6468)* 100  # = 1.9%

#calculating incidence cases among older age group

# Element at first row (record # at time 0 days)and second column (Sa)
simSEIR2hostcSHc[1,5] # 4221
# Element at twenty-first row (record # at time 9 days)and second column (Sa)
simSEIR2hostcSHc[10,5] # 4162.13
## Total number of infections among person aged >=20 years at the fair
(simSEIR2hostcSHc[1,5]) - (simSEIR2hostcSHc[10,5]) #58.87027 

(58.87027/8442) * 100 # = 0.69% 

#To plot simulated proportion incident influenza A (H3N2) virus infections among humans and swine 

plot(simSEIR2hostcSHc[,'time'],IC, ylab='Incidence',xlab='Duration of Fair (days)', ylim = c(0,0.035), type='b',lwd=2,  col="blue", main = "Proportion of Infectious Incidence of Swine Flu at Fair A", cex.main=0.85, font.main=14, sub = "With greater contact duration for <20 yr olds", font.sub=3)
lines(simSEIR2hostcSHc[,'time'],IA,type='b',lwd=2,col='olivedrab4')
lines(simSEIR2hostcSHc[,'time'],IS,type='b',lwd=2,col='deeppink')
legend('topleft',legend = (c('Humans<20yrs','Humans>=20yrs', 'Swine')),lty=1, lwd=2, col=c('blue', 'olivedrab4', 'deeppink'), bty='n', cex=.95)



#To plot simulated # incident cases of influenza A (H3N2) virus infections among humans and swine 

plot(simSEIR2hostcSHc[,'time'],simSEIR2hostcSHc[,'Ic'], ylab='No. of Incidence Cases',xlab='Duration of Fair (days)', ylim = c(0,42), type='b',lwd=2,  col="blue", main = "# of Infectious cases of Swine Flu at Fair A",cex.main=0.85, font.main=14, sub = "With greater contact duration for <20 yr olds", font.sub=3)
lines(simSEIR2hostcSHc[,'time'],simSEIR2hostcSHc[,'Ia'],type='b',lwd=2,col='olivedrab4')
lines(simSEIR2hostcSHc[,'time'],simSEIR2hostcSHc[,'Is'],type='b',lwd=2,col='deeppink')
legend('topleft',legend = (c('Humans<20yrs','Humans>=20yrs', 'Swine')),lty=1, lwd=2, col=c('blue', 'olivedrab4', 'deeppink'), bty='n', cex=.95)

