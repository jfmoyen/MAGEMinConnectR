#######################################################
#### Demonstration of the use of MAGEMin connector
#### Fractionation
#######################################################

library(MAGEMinConnectR)

MAGEMin_setup(JULIA_HOME="C:\\Program Files\\Julia-1.11.5\\bin")

# MAGEMin oxides - here for the ig dataset
# We will use MAGEMin output which returns all oxides, so it is easier to
# start with a longer list...
Xoxides2 <- c("SiO2", "Al2O3", "CaO", "MgO", "FeO",
              "K2O", "Na2O", "TiO2", "O", "Cr2O3","H2O")

rock = c(SiO2 = 57.28,
         Al2O3 = 16.5,
         CaO = 7.12,
         MgO = 5.42,
         FeO = 8.48/1.111,
         K2O = 0.76,
         Na2O = 3.52,
         TiO2 = 0.67,
         O = EvaluateO(8.48/1.111,0.3), # Recommended for basaltic andesite, Middlemost 1989
         Cr203 = 0,
         H2O = 5.0
)  #ATAC-4

T_range <- seq(1050,600,by=-25)
c0 <- rock
liqs <- as.data.frame(t(c(c0,Temp=NA,FF=1,remLiq=1)))
pb <- progress::progress_bar$new(total=length(T_range))

for(Temp in T_range){
  # cat(Temp,"\n")
  out <- MAGEMin(Xoxides2,X = c0,
                 Pkbar = 4.5,TC=Temp,
                 showMinimizationResults=F)
  Cl <- getSSComp(out,"liq","wt") * 100
  FF <- getPhProp(out)["liq"]

  if(is.na(FF)){
    liqs <- rbind(liqs,
                  c(rep(NA,length(c0)),Temp,0,0
                  )
    )
    break()
  }else{
    remLiq <- FF * tail(liqs,1)[,"remLiq"]

    liqs <- rbind(liqs,c(Cl,Temp=Temp,FF=FF,remLiq=remLiq))
    c0 <- Cl
  }

  pb$tick()
}

print(liqs)
