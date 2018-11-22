library(data.table)
library(tidyverse)
library(ggplot2)

#### Definitions of equations to calculate pH from measured
#### Salinity (Sal), Temperature in K (Tem) and mCP absorption ratio (Rspec = A434/A578)

#### Mosley et al. 2004

pHT.Mosley <- function(Sal, Tem, Rspec,
                      e1 = 0.00691, e2 = 2.222, e3 = 0.1331) {
  (1245.69/Tem) + 4.4572353 - (0.3238*(Sal^0.5)) + (0.0807*Sal) - (0.01157*(Sal^1.5)) + (0.000694*(Sal^2)) +
    log10( (Rspec-e1) / (e2-Rspec*e3) ) 
  }


#### Mueller and Rehder (2018)

pHT.Mueller <- function(Sal, Tem, Rspec){
  
  #first set of coefficients defines pK2e2 = f(Sal, Tem)
  
    1.08071477e+03                      -
    1.35394946e-01  *Sal^0.5            -   
    1.98063716e+02  *Sal^1.5            +
    6.31924397e+01  *Sal^2              -
    5.18141866e+00  *Sal^2.5            -
    2.66457425e+04  *Tem^-1             +
    5.08796578e+03  *Sal^1.5 * Tem^-1   -
    1.62454827e+03  *Sal^2 * Tem^-1     +
    1.33276788e+02  *Sal^2.5 * Tem^-1   -
    1.89671212e+02  *log(Tem)           +
    3.49038762e+01  *Sal^1.5 * log(Tem) -
    1.11336508e+01  *Sal^2 * log(Tem)   +
    9.12761930e-01  *Sal^2.5 * log(Tem) +
    3.27430677e-01  *Tem              -
    7.51448528e-04  *Sal^0.5 * Tem      +
    3.94838229e-04  *Sal * Tem          -
    6.00237876e-02  *Sal^1.5 * Tem      +
    1.90997693e-02  *Sal^2 * Tem        -
    1.56396488e-03  *Sal^2.5 * Tem      +
    
  #second set of coefficients includes the definition of mCP absorptivity ratios e1 and e3/e3
  #as determined by Liu et al. (2011) and defines the log-term calculation 
    
    
    log10(
      (Rspec -
         (-0.007762 + 4.5174e-5*Tem)) /
        (1 - (Rspec *  (- 0.020813 + 2.60262e-4*Tem + 1.0436e-4*(Sal-35))))
    )
}


#### create data frame with a range of combined theoretical Sal and Rspec values and
#### calculate pHT values incl. deviations according to Mosley and Mueller 

df <- data.table(expand.grid(
  Sal = seq(0,40,0.1),
  Tem = 298.15,
  RSpec = seq(0.2,1.5,0.1)
))



df <- df %>% 
  mutate(pHT.Mueller = pHT.Mueller(Sal, Tem, RSpec)) %>% 
  mutate(pHT.Mosley = pHT.Mosley(Sal, Tem, RSpec)) %>% 
  mutate(dpHT = pHT.Mueller - pHT.Mosley) 


#### plot results and write data file

df %>% 
ggplot(aes(Sal, RSpec, fill=dpHT, col="grey50"))+
  geom_raster(interpolate = TRUE)+
  scale_fill_gradient2(low = "blue", high = "red", mid="white", name=expression(Delta~pHT))+
  labs(x="Salinity", y="R ratio")+
  scale_x_continuous(expand = c(0,0), breaks = seq(0,40,5))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,5,0.1))+
  theme_test()
ggsave("dpH_Mueller_Mosley_Hovmoeller.jpg")

df %>% 
  ggplot(aes(Sal, dpHT, col=as.factor(RSpec))) +
  geom_line()+
  scale_color_viridis_d(name="R ratio")+
  labs(x="Salinity", name=expression(Delta~pHT))+
  theme_bw()
ggsave("dpH_Mueller_Mosley_linegraph.jpg")


write.csv(df, "dpH_Mueller_Mosley.csv", row.names = FALSE)





