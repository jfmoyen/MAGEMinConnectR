#######################################################
#### Demonstration of the use of MAGEMin connector
#### Tidyverse
#######################################################

library(MAGEMinConnectR)

library(tidyverse)
library(readxl)

MAGEMin_setup(JULIA_HOME="C:\\Program Files\\Julia-1.11.5\\bin")

# Here I use one of GCDkit's datasets, but you may want to use your own
# Pay attention in this case to FeO/Fe2O3/FeOt etc., also check H2O in your file
atac <- read_xlsx("D:\\R\\win-library\\4.3\\GCDkit\\Test_data\\atacazo.xlsx")

## The following may take a while !
# You are doing on average 1 minimization for each of 109 samples
#  at 0.2-0.5 s each on a modern machine this means about a minute
# For a shorter calculation,adjust N_sample

N_sample <- nrow(atac)
# N_sample <- 10
pb <- progress::progress_bar$new(total=N_sample)
atac %>%  slice_head(n=N_sample)  %>%
  mutate(FeO = Fe2O3/1.111,
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
