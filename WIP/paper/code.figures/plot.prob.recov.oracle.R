# libraries 
####### libraries #######

library(ggplot2)
library(scales)
library(latex2exp)
library(gridExtra)
library(tidyverse)


####### Asymptotic regimes #######


log.nb.inactive1 <-function(n){
  return(log(n/2))
}

log.nb.inactive2 <-function(n){
  return(log(n))
}

log.nb.inactive3 <-function(n){
  return(2*log(n))
}

log.nb.inactive4 <-function(n){
  return(n/20)
}

log.nb.inactive5<- function(n){
  return(log(log(n)))
}

log.nb.inactive6 <-function(n){
  return(0.5*log(n))
}

log.nb.inactive7 <-function(n){
  return(log(1.5*n))
}

log.nb.active <- function(n){
  log.s <- log(3*log(n))
  return(log.s)
}

log.nb.active.j <- function(n){
  log.sj <- log(0.5)+log.nb.active(n)
  return(log.sj)
}

bmin <- function(n,fun.log.nb.inactive,fun.log.nb.active,constant=FALSE){
bmin.val <- sqrt(2*(fun.log.nb.inactive(n)+fun.log.nb.active(n))/n)
  return(1/2)
}


ratio.scen1<-function(n,
                      fun.log.nb.inactive=log.nb.inactive7,
                      fun.log.nb.active=log.nb.active,
                      fun.log.nb.inactive.j2=log.nb.inactive6,
                      fun.log.nb.active.j=log.nb.active.j,
                      fun.bmin=bmin){
  fun.log.nb.inactive.j1 <- function(n){
    return(fun.log.nb.inactive(n)+log(1-exp(fun.log.nb.inactive.j2(n)-fun.log.nb.inactive(n))))
  }
  return(exp(fun.log.nb.inactive.j1(n)-fun.log.nb.inactive(n))*exp(-(1/8)*n*(1.3*bmin(n,fun.log.nb.inactive,fun.log.nb.active)-bmin(n,fun.log.nb.inactive,fun.log.nb.active)))+exp(fun.log.nb.inactive.j2(n)-fun.log.nb.inactive(n)))
}

ratio.scen2<-function(n,
                      fun.log.nb.inactive=log.nb.inactive4,
                      fun.log.nb.active=log.nb.active,
                      fun.log.nb.inactive.j2=log.nb.inactive3,
                      fun.log.nb.active.j=log.nb.active.j,
                      fun.bmin=bmin){
  fun.log.nb.inactive.j1 <- function(n){
    return(fun.log.nb.inactive(n)+log(1-exp(fun.log.nb.inactive.j2(n)-fun.log.nb.inactive(n))))
  }
  return(exp(fun.log.nb.inactive.j1(n)-fun.log.nb.inactive(n))*exp(-(1/8)*n*(1.3*bmin(n,fun.log.nb.inactive,fun.log.nb.active)-bmin(n,fun.log.nb.inactive,fun.log.nb.active)))+exp(fun.log.nb.inactive.j2(n)-fun.log.nb.inactive(n)))
}


ratio.scen3<-function(n,
    fun.log.nb.inactive=log.nb.inactive2,
    fun.log.nb.active=log.nb.active,
    fun.log.nb.inactive.j2=log.nb.inactive5,
    fun.log.nb.active.j=log.nb.active.j,
    fun.bmin=bmin){
  fun.log.nb.inactive.j1 <- function(n){
    return(fun.log.nb.inactive(n)+log(1-exp(fun.log.nb.inactive.j2(n)-fun.log.nb.inactive(n))))
  }
  return(exp(fun.log.nb.inactive.j1(n)-fun.log.nb.inactive(n))*exp(-(1/8)*n*(1.3*bmin(n,fun.log.nb.inactive,fun.log.nb.active)-bmin(n,fun.log.nb.inactive,fun.log.nb.active)))+exp(fun.log.nb.active.j(n)-fun.log.nb.inactive(n)))
}


ratio.scen4<-function(n,
                      fun.log.nb.inactive=log.nb.inactive2,
                      fun.log.nb.active=log.nb.active,
                      fun.log.nb.inactive.j2=log.nb.inactive1,
                      fun.log.nb.active.j=log.nb.active.j,
                      fun.bmin=bmin){
  fun.log.nb.inactive.j1 <- function(n){
    return(fun.log.nb.inactive(n)+log(1-exp(fun.log.nb.inactive.j2(n)-fun.log.nb.inactive(n))))
  }
  return(exp(fun.log.nb.inactive.j1(n)-fun.log.nb.inactive(n))*exp(-(1/8)*n*(1*bmin(n,fun.log.nb.inactive,fun.log.nb.active)-bmin(n,fun.log.nb.inactive,fun.log.nb.active)))+exp(fun.log.nb.inactive.j2(n)-fun.log.nb.inactive(n)))
}

#### Conv rate values, same betamin across blocks

df.samebmin <- data.frame(n = seq(1000,10^6,200))

df.samebmin$ratio.scen1 <- sapply(df.samebmin$n,ratio.scen1)
df.samebmin$ratio.scen2 <- sapply(df.samebmin$n,ratio.scen2)
df.samebmin$ratio.scen3 <- sapply(df.samebmin$n,ratio.scen3)
df.samebmin$ratio.scen4 <- sapply(df.samebmin$n,ratio.scen4)

ggplot(df.samebmin,aes(x=n)) +
  geom_line(aes(y=ratio.scen1, linetype = '1_')) +
  geom_line(aes(y=ratio.scen2, linetype = '2_')) +
  geom_line(aes(y=ratio.scen3 , linetype = '3_')) +
  #geom_line(aes(y=ratio.scen4 , linetype = '4')) +
  scale_x_log10(labels = scales::comma) +
  scale_linetype_manual('Example',values=c("1_"=1, "2_"=2,"3_"=3,"4_"=4),labels=c(1,2,3,4))+
  ylab('Ratio of oracle convergence rates') +
  coord_cartesian(xlim = c(10^3,10^5))+
  theme_light(base_size = 8)+
  #ggtitle('Ratio of smallest signal recoverable')+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.margin=margin(0,0,0,0), 
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.position="bottom")

ggsave(paste0('p.ratio.orconvrate.pdf'), width = 60, height = 72, units='mm')


ggplot(df.samebmin,aes(x=n)) +
  geom_line(aes(y=ratio.scen1, linetype = '1_')) +
  geom_line(aes(y=ratio.scen2, linetype = '2_')) +
  #geom_line(aes(y=ratio.scen3 , linetype = '3_')) +
  #geom_line(aes(y=ratio.scen4 , linetype = '4')) +
  scale_x_log10(labels = scales::comma) +
  scale_linetype_manual('Example',values=c("1_"=1, "2_"=2),labels=c(1,2))+
  ylab('Ratio of oracle convergence rates') +
  coord_cartesian(xlim = c(10^3,10^5))+
  theme_light(base_size = 12)+
  #ggtitle('Ratio of smallest signal recoverable')+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.margin=margin(0,0,0,0), 
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.position="bottom")

ggsave(paste0('p.ratio.orconvrate.posterA1.pdf'), width = 110, height = 100, units='mm')


ggplot(df.samebmin,aes(x=n)) +
  geom_line(aes(y=ratio.scen1, linetype = '1_')) +
  geom_line(aes(y=ratio.scen2, linetype = '2_')) +
  #geom_line(aes(y=ratio.scen3 , linetype = '3_')) +
  #geom_line(aes(y=ratio.scen4 , linetype = '4')) +
  scale_x_log10(labels = scales::comma) +
  scale_linetype_manual('Example',values=c("1_"=1, "2_"=2),labels=c(1,2))+
  ylab('Ratio of oracle convergence rates') +
  coord_cartesian(xlim = c(10^3,10^5))+
  theme_light(base_size = 12)+
  #ggtitle('Ratio of smallest signal recoverable')+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.margin=margin(0,0,0,0), 
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.position="bottom")

ggsave(paste0('p.ratio.orconvrate.posterA0.pdf'), width = 150, height = 120, units='mm')
############# BACK 

####### Probability of model recovery with and without external information ######


# prob.recov.bound <- function(n,betamin,fun.log.nb.inactive,fun.log.nb.active){
#   return(2*exp(-(n*betamin^2/8-max(fun.log.nb.inactive(n),fun.log.nb.active(n)))))
# }
# 
# #### Standard convergence rates  
# conv.rate.std.scen1 <- function(n,
#                                 fun.log.nb.inactive=log.nb.inactive7,
#                                 fun.log.nb.active=log.nb.active,
#                                 fun.bmin = bmin){
#   betamin <- fun.bmin(n,fun.log.nb.inactive,fun.log.nb.active,constant=FALSE)
#   return(prob.recov.bound(n,betamin,fun.log.nb.inactive,fun.log.nb.active))
# }
# 
# conv.rate.std.scen2 <- function(n,
#                                 fun.log.nb.inactive=log.nb.inactive4,
#                                 fun.log.nb.active=log.nb.active,
#                                 fun.bmin = bmin){
#   betamin <- fun.bmin(n,fun.log.nb.inactive,fun.log.nb.active,constant=FALSE)
#   return(prob.recov.bound(n,betamin,fun.log.nb.inactive,fun.log.nb.active))
# }
# 
# conv.rate.std.scen3 <- function(n,
#                                 fun.log.nb.inactive=log.nb.inactive2,
#                                 fun.log.nb.active=log.nb.active,
#                                 fun.bmin = bmin){
#   betamin <- fun.bmin(n,fun.log.nb.inactive,fun.log.nb.active,constant=FALSE)
#   return(prob.recov.bound(n,betamin,fun.log.nb.inactive,fun.log.nb.active))
# }
# 
# conv.rate.std.scen4 <- function(n,
#                                 fun.log.nb.inactive=log.nb.inactive2,
#                                 fun.log.nb.active=log.nb.active,
#                                 fun.bmin = bmin){
#   betamin <- fun.bmin(n,fun.log.nb.inactive,fun.log.nb.active,constant=FALSE)
#   return(prob.recov.bound(n,betamin,fun.log.nb.inactive,fun.log.nb.active))
# }

# #### Block convergence rates, same betamin across blocks
# 
# conv.rate.block.scen1 <- function(n,
#                                   fun.log.nb.inactive=log.nb.inactive7,
#                                   fun.log.nb.active=log.nb.active,
#                                   fun.log.nb.inactive.j2=log.nb.inactive6,
#                                   fun.log.nb.active.j=log.nb.active.j,
#                                   fun.bmin=bmin){
#   betamin <- fun.bmin(n,fun.log.nb.inactive,fun.log.nb.active,constant=FALSE)
#   
#   prob.block2 <- prob.recov.bound(n,betamin,fun.log.nb.inactive.j2,fun.log.nb.active.j)
#   
#   fun.log.nb.inactive.j1 <- function(n){
#     return(fun.log.nb.inactive(n)+log(1-exp(fun.log.nb.inactive.j2(n)-fun.log.nb.inactive(n))))
#   }
#   
#   
#   prob.block1 <- prob.recov.bound(n,betamin,fun.log.nb.inactive.j2,fun.log.nb.active.j)
#   return( prob.block1+prob.block2)
# }
# 
# 
# conv.rate.block.scen2 <- function(n,
#                                   fun.log.nb.inactive=log.nb.inactive4,
#                                   fun.log.nb.active=log.nb.active,
#                                   fun.log.nb.inactive.j2=log.nb.inactive3,
#                                   fun.log.nb.active.j=log.nb.active.j,
#                                   fun.bmin=bmin){
#   betamin <- fun.bmin(n,fun.log.nb.inactive,fun.log.nb.active,constant=FALSE)
#   
#   prob.block2 <- prob.recov.bound(n,betamin,fun.log.nb.inactive.j2,fun.log.nb.active.j)
#   
#   fun.log.nb.inactive.j1 <- function(n){
#     return(fun.log.nb.inactive(n)+log(1-exp(fun.log.nb.inactive.j2(n)-fun.log.nb.inactive(n))))
#   }
#   
#   
#   prob.block1 <- prob.recov.bound(n,betamin,fun.log.nb.inactive.j2,fun.log.nb.active.j)
#   return( prob.block1+prob.block2)
# }
# 
# conv.rate.block.scen3 <- function(n,
#                                   fun.log.nb.inactive=log.nb.inactive2,
#                                   fun.log.nb.active=log.nb.active,
#                                   fun.log.nb.inactive.j2=log.nb.inactive5,
#                                   fun.log.nb.active.j=log.nb.active.j,
#                                   fun.bmin=bmin){
#   betamin <- fun.bmin(n,fun.log.nb.inactive,fun.log.nb.active,constant=FALSE)
#   
#   prob.block2 <- prob.recov.bound(n,betamin,fun.log.nb.inactive.j2,fun.log.nb.active.j)
#   
#   fun.log.nb.inactive.j1 <- function(n){
#     return(fun.log.nb.inactive(n)+log(1-exp(fun.log.nb.inactive.j2(n)-fun.log.nb.inactive(n))))
#   }
#   
#   
#   prob.block1 <- prob.recov.bound(n,betamin,fun.log.nb.inactive.j2,fun.log.nb.active.j)
#   return( prob.block1+prob.block2)
# }
# 
# 
# conv.rate.block.scen4 <- function(n,
#                                   fun.log.nb.inactive=log.nb.inactive2,
#                                   fun.log.nb.active=log.nb.active,
#                                   fun.log.nb.inactive.j2=log.nb.inactive1,
#                                   fun.log.nb.active.j=log.nb.active.j,
#                                   fun.bmin=bmin){
#   betamin <- fun.bmin(n,fun.log.nb.inactive,fun.log.nb.active,constant=FALSE)
#   
#   prob.block2 <- prob.recov.bound(n,betamin,fun.log.nb.inactive.j2,fun.log.nb.active.j)
#   
#   fun.log.nb.inactive.j1 <- function(n){
#     return(fun.log.nb.inactive(n)+log(1-exp(fun.log.nb.inactive.j2(n)-fun.log.nb.inactive(n))))
#   }
#   
#   
#   prob.block1 <- prob.recov.bound(n,betamin,fun.log.nb.inactive.j2,fun.log.nb.active.j)
#   return( prob.block1+prob.block2)
# }
# 
# 
# df.samebmin$std.scen1 <- sapply(df.samebmin$n,conv.rate.std.scen1)
# df.samebmin$std.scen2 <- sapply(df.samebmin$n,conv.rate.std.scen2)
# df.samebmin$std.scen3 <- sapply(df.samebmin$n,conv.rate.std.scen3)
# df.samebmin$std.scen4 <- sapply(df.samebmin$n,conv.rate.std.scen4)
# 
# 
# df.samebmin$block.scen1 <- sapply(df.samebmin$n,conv.rate.block.scen1)
# df.samebmin$block.scen2 <- sapply(df.samebmin$n,conv.rate.block.scen2)
# df.samebmin$block.scen3 <- sapply(df.samebmin$n,conv.rate.block.scen3)
# df.samebmin$block.scen4 <- sapply(df.samebmin$n,conv.rate.block.scen4)
# 
# df.samebmin$bmin.j2 <- 'same'
# 
# 
# ggplot(df.samebmin,aes(x=n)) +
#   geom_line(aes(y=block.scen1/std.scen1, linetype = '1')) +
#   geom_line(aes(y=block.scen2/std.scen2, linetype = '2')) +
#   geom_line(aes(y=block.scen3/std.scen3 , linetype = '3')) +
#   geom_line(aes(y=block.scen4/std.scen4 , linetype = '4')) +
#   scale_x_log10() +
#   scale_linetype_manual(values=c(1,2,3,4))+
#   ylab(TeX(r"(Probability)")) +
#   coord_cartesian(xlim = c(10^3,10^4)
#                   #ylim = c(0,1)
#                   )+
#   theme_light(base_size = 13)+
#   theme(panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),legend.text.align = 0)


