#### Demonstration of the use of MAGEMin connector ####
#######################################################

library(MAGEMinConnectR)

#### Simple case - one compo, manually defined ####

MAGEMin_setup(JULIA_HOME="C:\\Program Files\\Julia-1.11.5\\bin")

# MAGEMin oxides - here for the ig dataset
Xoxides <- c("SiO2", "Al2O3", "CaO", "MgO", "FeO",
             "K2O", "Na2O", "TiO2", "O", "H2O")

rock = c(SiO2 = 57.28,
         Al2O3 = 16.5,
         CaO = 7.12,
         MgO = 5.42,
         FeO = 8.48/1.111,
         K2O = 0.76,
         Na2O = 3.52,
         TiO2 = 0.67,
         O = EvaluateO(8.48/1.111,0.3), # Recommended for basaltic andesite, Middlemost 1989
         H2O = 5.0 # julia needs float to be explicitely given as decimals (5.0 not 5)
         )  #ATAC-4

# The first minimization always takes a bit longer
# because julia compiles the code the first time...
out <- MAGEMin(Xoxides,X = rock,
               Pkbar = 4.5,TC=900,
               showResults=T)

# Phase Proportions
getPhProp(out)

# Phase composition
getSSComp(out,"liq")
getSSComp(out,"pl","emmol")

#### With a buffer ####

MAGEMin_setup(JULIA_HOME="C:\\Program Files\\Julia-1.11.5\\bin",extraMAGEMinParams = 'buffer=\"qfm\"')

# MAGEMin oxides - here for the ig dataset
Xoxides <- c("SiO2", "Al2O3", "CaO", "MgO", "FeO",
             "K2O", "Na2O", "TiO2", "O", "H2O")

rock = c(SiO2 = 57.28,
         Al2O3 = 16.5,
         CaO = 7.12,
         MgO = 5.42,
         FeO = 8.48/1.111,
         K2O = 0.76,
         Na2O = 3.52,
         TiO2 = 0.67,
         O = 10, # Excess
         H2O = 5.0
)  #ATAC-4

out <- MAGEMin(Xoxides,X = rock,
               Pkbar = 4.5,TC=1000,
               extraMAGEMinParams = 'B=2.0',
               showResults=T)

# Phase Proportions
getPhProp(out)

# Phase composition
getSSComp(out,"liq")
getSSComp(out,"pl","emmol")


#### Fractionation ####

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
                 showResults=F)
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

liqs

#### With tidyverse functions ####
library(tidyverse)
library(readxl)

MAGEMin_setup(JULIA_HOME="C:\\Program Files\\Julia-1.11.5\\bin")

# Here I use one of GCDkit's datasets, but you may want to use your own
# Pay attention in this case to FeO/Fe2O3/FeOt etc., also check H2O in your file
atac <- read_xlsx("D:\\R\\win-library\\4.3\\GCDkit\\Test_data\\atacazo.xlsx")

pb <- progress::progress_bar$new(total=nrow(atac))
atac %>%  mutate(FeO = Fe2O3/1.111,
                O =  EvaluateO(FeO,0.3),
                H2O = 5) %>% # No H2O given, add manually (excess)
        rowwise() %>%
        mutate(
          MAGEMin =list(
            {
              pb$tick()
              MAGEMin(X = across({{Xoxides}}), Pkbar=5,TC=850,Xoxides = Xoxides)
              }
          )
        ) %>%
        {.} -> results

# Extract useable results
results %>%
  rowwise() %>% #Rowwise is important because getPhProp() is not vectorized
  mutate(phprop = list(getPhProp(MAGEMin)) ) %>%
  unnest_wider(phprop,names_sep="_",names_repair = "unique") %>%
  rowwise() %>%
  mutate(plComp_em = list(getSSComp(MAGEMin,"pl","emmol")) ) %>%
  unnest_wider(plComp_em,names_sep="_",names_repair = "unique") %>%
  {.} -> rr

rr %>% ggplot()+
  geom_point(aes(x=SiO2,y=100*plComp_em_an/(plComp_em_an+plComp_em_ab)))+
  ylab("plag An%")


