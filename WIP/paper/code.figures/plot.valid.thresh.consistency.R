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

bmin <- function(n){
  #bmin <- sqrt(8*fun.log.nb.inactive(n)/n*(2*log(log(n))/3))
  return(1/10)
}

bmin.block1 <- function(n){
  #bmin <- sqrt(8*fun.log.nb.inactive(n)/n*(2*log(log(n))/3))
  return(2/3)
}

####### lower bounds of valid ranges for thresholds #######

lower.block1.scen1 <- function(n){
  return(sqrt(2*(log.nb.inactive7(n)+log(1-exp(log.nb.inactive6(n)-log.nb.inactive7(n))))/n))
}

lower.block2.scen1<- function(n){
  return(sqrt(2*log.nb.inactive6(n)/n))
}

lower.std.scen1<- function(n){
  return(sqrt(2*log.nb.inactive7(n)/n))
}

lower.block1.scen2 <- function(n){
  return(sqrt(2*(log.nb.inactive4(n)+log(1-exp(log.nb.inactive3(n)-log.nb.inactive4(n))))/n))
}

lower.block2.scen2 <- function(n){
  return(sqrt(2*log.nb.inactive3(n)/n))
}

lower.std.scen2 <- function(n){
  return(sqrt(2*log.nb.inactive4(n)/n))
}

lower.block1.scen3 <- function(n){
  return(sqrt(2*(log.nb.inactive2(n)+log(1-exp(log.nb.inactive5(n)-log.nb.inactive2(n))))/n))
}

lower.block2.scen3<- function(n){
  return(sqrt(2*log.nb.inactive5(n)/n))
}

lower.std.scen3<- function(n){
  return(sqrt(2*log.nb.inactive2(n)/n))
}


lower.block12.scen4<- function(n){
  return(sqrt(2*log.nb.inactive1(n)/n))
}


lower.std.scen4 <- function(n){
  return(sqrt(2*log.nb.inactive2(n)/n))
}


####### upper bounds of valid ranges for thresholds #######

upper.std<- function(n){
  return(bmin(n)-sqrt(2*log.nb.active(n)/n))
}

upper.block1<- function(n){
  return(bmin.block1(n)-sqrt(2*(log.nb.active(n)+log(0.5))/n))
}

upper.block2<- function(n){
  return(bmin(n)-sqrt(2*(log.nb.active(n)+log(0.5))/n))
}

####### Bounds values #######

n <- seq(800,10^5,100)

df.thresh.scen1.std <- data.frame(n = n, scenario='Example 1', est='')
df.thresh.scen1.block1 <- data.frame(n = n, scenario='Example 1', est='1')
df.thresh.scen1.block2 <- data.frame(n = n, scenario='Example 1', est='2')

df.thresh.scen2.std <- data.frame(n = n, scenario='Example 2', est='')
df.thresh.scen2.block1 <- data.frame(n = n, scenario='Example 2', est='1')
df.thresh.scen2.block2 <- data.frame(n = n, scenario='Example 2', est='2')

df.thresh.scen3.std <- data.frame(n = n, scenario='Example 3', est='')
df.thresh.scen3.block1 <- data.frame(n = n, scenario='Example 3', est='1')
df.thresh.scen3.block2 <- data.frame(n = n, scenario='Example 3', est='2')

df.thresh.scen4.std <- data.frame(n = n, scenario='Example 4', est='')
df.thresh.scen4.block1 <- data.frame(n = n, scenario='Example 4', est='1')
df.thresh.scen4.block2 <- data.frame(n = n, scenario='Example 4', est='2')

df.thresh.scen1.std$lower <- sapply(n,lower.std.scen1)
df.thresh.scen2.std$lower <- sapply(n,lower.std.scen2)
df.thresh.scen3.std$lower <- sapply(n,lower.std.scen3)
df.thresh.scen4.std$lower <- sapply(n,lower.std.scen4)

df.thresh.scen1.std$upper <- sapply(n,upper.std)
df.thresh.scen2.std$upper <- sapply(n,upper.std)
df.thresh.scen3.std$upper <- sapply(n,upper.std)
df.thresh.scen4.std$upper <- sapply(n,upper.std)

df.thresh.scen1.block1$lower <- sapply(n,lower.block1.scen1)
df.thresh.scen2.block1$lower <- sapply(n,lower.block1.scen2)
df.thresh.scen3.block1$lower <- sapply(n,lower.block1.scen3)
df.thresh.scen4.block1$lower <- sapply(n,lower.block12.scen4)

df.thresh.scen1.block1$upper <- sapply(n,upper.block1)
df.thresh.scen2.block1$upper <- sapply(n,upper.block1)
df.thresh.scen3.block1$upper <- sapply(n,upper.block1)
df.thresh.scen4.block1$upper <- sapply(n,upper.block1)

df.thresh.scen1.block2$lower <- sapply(n,lower.block2.scen1)
df.thresh.scen2.block2$lower <- sapply(n,lower.block2.scen2)
df.thresh.scen3.block2$lower <- sapply(n,lower.block2.scen3)
df.thresh.scen4.block2$lower <- sapply(n,lower.block12.scen4)

df.thresh.scen1.block2$upper <- sapply(n,upper.block2)
df.thresh.scen2.block2$upper <- sapply(n,upper.block2)
df.thresh.scen3.block2$upper <- sapply(n,upper.block2)
df.thresh.scen4.block2$upper <- sapply(n,upper.block2)

df.thresh.allscen <- rbind(
  df.thresh.scen1.std,
  df.thresh.scen1.block1,
  df.thresh.scen1.block2,
  df.thresh.scen2.std,
  df.thresh.scen2.block1,
  df.thresh.scen2.block2,
  df.thresh.scen3.std,
  df.thresh.scen3.block1,
  df.thresh.scen3.block2,
  df.thresh.scen4.std,
  df.thresh.scen4.block1,
  df.thresh.scen4.block2
)

df.thresh.allscen$valid <- as.character(df.thresh.allscen$upper > df.thresh.allscen$lower)

####### plot of valid ranges for thresholds #######

 
p1<- df.thresh.allscen %>% filter(scenario %in% c('Example 1','Example 2')) %>%
  ggplot(aes(x=n))+
  facet_grid(est~scenario,labeller = label_bquote(tau [.(est)] ), scales='free') +
  geom_line(aes(y=lower,linetype='1_'),size=0.5) +
  geom_line(aes(y=upper,linetype='2_'),size=0.5) +
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=valid),alpha=0.3) +
  scale_x_log10(labels =scales::comma) +
  scale_fill_manual(
    values=c('FALSE'='red','TRUE'="lightgreen"),
    guide="none")+
  scale_linetype_manual('', values=c("1_"=2,"2_"=1),labels=c('Lower bound','Upper bound'))+
  scale_y_continuous(position = "left",n.breaks = 4) +
  coord_cartesian(xlim = c(1*10^3,1*10^4))+
  theme_light()+
  xlab('n') +
  theme(axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),legend.text.align = 0, legend.position = 'bottom',
        strip.background = element_rect(colour="transparent", fill="transparent"),strip.text = element_text(colour = 'black'),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),panel.spacing = unit(1.2, "lines"))

p1

p2<- df.thresh.allscen %>% filter(scenario %in% c('Example 3','Example 4')) %>%
  ggplot(aes(x=n))+
  facet_grid(est~scenario,labeller = label_bquote(tau [.(est)] ), switch='y', scales='free') +
  geom_line(aes(y=lower,linetype='1_'),size=0.5) +
  geom_line(aes(y=upper,linetype='2_'),size=0.5) +
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=valid),alpha=0.3) +
  scale_x_log10(labels =scales::comma) +
  scale_fill_manual(
    values=c('FALSE'='red','TRUE'="lightgreen"),
    guide="none")+
  scale_linetype_manual('', values=c("1_"=2,"2_"=1),labels=c('Lower bound','Upper bound'))+
  scale_y_continuous(position = "right",n.breaks = 4) +
  coord_cartesian(xlim = c(1*10^3,1*10^4),clip="off")+
  theme_light()+
  xlab('n') +
  theme(axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),legend.text.align = 0, legend.position = 'bottom',
        strip.background = element_rect(colour="transparent", fill="transparent"),strip.text = element_text(colour = 'black'),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),panel.spacing = unit(1.2, "lines"))

p2

#### one by one
i=0
for(scen in c('Example 1','Example 2','Example 3','Example 4')){
  i=i+1
  for(est2 in c('','1','2')){
    p<- df.thresh.allscen %>% filter(scenario ==scen, est==est2) %>% 
      ggplot(aes(x=n))+
      geom_line(aes(y=lower,linetype='1_'),size=0.5) +
      geom_line(aes(y=upper,linetype='2_'),size=0.5) +
      geom_ribbon(aes(ymin=lower,ymax=upper,fill=valid),alpha=0.3) +
      scale_x_log10(labels =scales::comma) +
      scale_fill_manual(
        values=c('FALSE'='red','TRUE'="blue"),
        guide="none")+
      scale_linetype_manual('', values=c("1_"=2,"2_"=1),labels=c('Lower bound','Upper bound'),guide='none')+
      scale_y_continuous(position = "left",n.breaks = 4) +
      coord_cartesian(xlim = c(1*10^3,1*10^4))+
      theme_light(base_size=5)+
      xlab('n') +
      #ylab(bquote(tau[.(est2)])) +
      #labs(title=bquote(.(scen)~'-'~tau[.(est2)])) +
      theme(axis.title.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank())
    ggsave(paste0('p.scen',i,'tau',est2,'.pdf'), width = 33, height = 30, units='mm')
    assign(paste0('p.scen',i,'tau',est2),p)
  }
}



#### one by one - poster A1
i=0
for(scen in c('Example 1','Example 2','Example 3','Example 4')){
  i=i+1
  for(est2 in c('','1','2')){
    p<- df.thresh.allscen %>% filter(scenario ==scen, est==est2) %>% 
      ggplot(aes(x=n))+
      geom_line(aes(y=lower,linetype='1_'),size=0.5) +
      geom_line(aes(y=upper,linetype='2_'),size=0.5) +
      geom_ribbon(aes(ymin=lower,ymax=upper,fill=valid),alpha=0.3) +
      scale_x_log10(labels =scales::comma) +
      scale_fill_manual(
        values=c('FALSE'='red','TRUE'="blue"),
        guide="none")+
      scale_linetype_manual('', values=c("1_"=2,"2_"=1),labels=c('Lower bound','Upper bound'),guide='none')+
      scale_y_continuous(position = "left",n.breaks = 4) +
      coord_cartesian(xlim = c(1*10^3,1*10^4))+
      theme_light(base_size=10)+
      xlab('n') +
      #ylab(bquote(tau[.(est2)])) +
      #labs(title=bquote(.(scen)~'-'~tau[.(est2)])) +
      theme(axis.title.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank())
    ggsave(paste0('p.scen',i,'tau',est2,'.pdf'), width = 41, height = 41, units='mm')
    assign(paste0('p.scen',i,'tau',est2),p)
  }
}

#### one by one - poster A0
i=0
for(scen in c('Example 1','Example 2','Example 3','Example 4')){
  i=i+1
  for(est2 in c('','1','2')){
    p<- df.thresh.allscen %>% filter(scenario ==scen, est==est2) %>% 
      ggplot(aes(x=n))+
      geom_line(aes(y=lower,linetype='1_'),size=0.5) +
      geom_line(aes(y=upper,linetype='2_'),size=0.5) +
      geom_ribbon(aes(ymin=lower,ymax=upper,fill=valid),alpha=0.3) +
      scale_x_log10(labels =scales::comma) +
      scale_fill_manual(
        values=c('FALSE'='red','TRUE'="blue"),
        guide="none")+
      scale_linetype_manual('', values=c("1_"=2,"2_"=1),labels=c('Lower bound','Upper bound'),guide='none')+
      scale_y_continuous(position = "left",n.breaks = 4) +
      coord_cartesian(xlim = c(1*10^3,1*10^4))+
      theme_light(base_size=12)+
      xlab('n') +
      #ylab(bquote(tau[.(est2)])) +
      #labs(title=bquote(.(scen)~'-'~tau[.(est2)])) +
      theme(axis.title.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank())
    ggsave(paste0('p.scen',i,'tau',est2,'.pdf'), width = 55, height = 55, units='mm')
    assign(paste0('p.scen',i,'tau',est2),p)
  }
}



#### one by one
i=0
for(scen in c('Example 1','Example 2','Example 3','Example 4')){
  i=i+1
  for(est2 in c('','1','2')){
    p<- df.thresh.allscen %>% filter(scenario ==scen, est==est2) %>% 
      ggplot(aes(x=n))+
      geom_line(aes(y=lower,linetype='1_'),size=0.5) +
      geom_line(aes(y=upper,linetype='2_'),size=0.5) +
      geom_ribbon(aes(ymin=lower,ymax=upper,fill=valid),alpha=0.3) +
      scale_x_log10(labels =scales::comma) +
      scale_fill_manual(
        values=c('FALSE'='red','TRUE'="lightgreen"),
        guide="none")+
      scale_linetype_manual('', values=c("1_"=2,"2_"=1),labels=c('Lower bound','Upper bound'),guide='none')+
      scale_y_continuous(position = "left",n.breaks = 4) +
      coord_cartesian(xlim = c(1*10^3,1*10^4))+
      theme_light(base_size=8)+
      xlab('n') +
      #ylab(bquote(tau[.(est2)])) +
      #labs(title=bquote(.(scen)~'-'~tau[.(est2)])) +
      theme(axis.title.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank())
    ggsave(paste0('p60mm.scen',i,'tau',est2,'.pdf'), width = 60, height = 18, units='mm')
    assign(paste0('p.scen',i,'tau',est2),p)
  }
}













