#######################################################
#### Demonstration of the use of MAGEMin connector
#### Tidyverse to process solidus/liquidus temp
#######################################################

library(MAGEMinConnectR)

library(tidyverse)
library(readxl)

MAGEMin_setup(JULIA_HOME="C:\\Program Files\\Julia-1.11.5\\bin")

# Here I use one of GCDkit's datasets, but you may want to use your own
# Pay attention in this case to FeO/Fe2O3/FeOt etc., also check H2O in your file
atac <- read_xlsx("D:\\R\\win-library\\4.3\\GCDkit\\Test_data\\atacazo.xlsx")

## The following will take a while !
# You are doing on average 5-10 minimizations (twice) for each of 109 samples
# (we will find out how many exactly latter on)
# So something like 1-2000, at 0.2-0.5 s each on a modern machine this means 5-10 mn
# For a shorter calculation,adjust N_sample

N_sample <- nrow(atac)
# N_sample <- 5

pb <- progress::progress_bar$new(total=N_sample*2)
atac %>%
  slice_head(n=N_sample)  %>%
  mutate(FeO = Fe2O3/1.111,
                 O =  EvaluateO(FeO,0.3),
                 H2O = 5) %>% # No H2O given, add manually (excess)
  rowwise() %>%
  mutate(
    solidus =list(
      {
        pb$tick()
        findSolidus(X = across({{Xoxides}}),
                    Pkbar=5,Xoxides = Xoxides,
                    printPhaseBoundaryResults =F,Tlow=650,Thigh=1000)
      }
    ),
    liquidus =list(
      {
        pb$tick()
        findLiquidus(X = across({{Xoxides}}),
                    Pkbar=5,Xoxides = Xoxides,
                    printPhaseBoundaryResults =F,Tlow=900,Thigh=1400)
      }
    )
  ) %>%
  {.} -> results

# Extract useable results
results %>%
  mutate(T_solidus = solidus$T_C,
         solidus_found = solidus$solidus_found,
         solidus_minimizations = solidus$iterations,
         T_liquidus = liquidus$T_C,
         liquidus_found = liquidus$liquidus_found,
         liquidus_minimizations = liquidus$iterations) %>%
  # Liquidus phase
         rowwise() %>% #Rowwise is important because getPhProp() is not vectorized
         mutate(phprop = list(getPhProp(liquidus)) ) %>%
  {.} -> rr

# Total minimizations done:
sum(rr$liquidus_minimizations)+sum(rr$solidus_minimizations)

rr %>% ggplot()+
  geom_point(aes(x=SiO2,y=T_solidus,shape=Volcano,colour="Solidus"))+
  geom_point(aes(x=SiO2,y=T_liquidus,shape=Volcano,colour="Liquidus"))+
  scale_color_manual(values=c(Solidus="blue",Liquidus="red"))+
  guides(shape = guide_legend(override.aes = list( colour="black") ),
         color = guide_legend(title=""))+
  labs(y="Temperature")
