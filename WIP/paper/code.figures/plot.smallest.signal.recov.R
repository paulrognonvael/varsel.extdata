# libraries 
####### libraries #######

library(ggplot2)
library(scales)
library(latex2exp)
library(gridExtra)
library(tidyverse)


path = '~/GitHub/varsel.extdata/WIP/paper/code.figures/'

####### Plot ratio of betamin lower bounds #######


log.nb.active <- function(n){
  log.s <- log(3*log(n))
  return(log.s)
}


smallest.std.scen1 <- function(n){
  return(sqrt(log(1.5*n)/n) + sqrt(log.nb.active(n)/n))
}


smallest.std.scen2 <- function(n){
  return(1/4 + sqrt(log.nb.active(n)/n))
}

smallest.std.scen3 <- function(n){
  return(sqrt(log(n)/n) + sqrt(log.nb.active(n)/n))
}

smallest.std.scen4 <- function(n){
  return(sqrt(log(n)/n) + sqrt(log.nb.active(n)/n))
}

smallest.block.scen1 <- function(n){
  return(sqrt(0.5*log(n)/n) + sqrt((log(0.5)+log.nb.active(n))/n))
}


smallest.block.scen2 <- function(n){
  return(sqrt(2*log(n)/n) + sqrt((log(0.5)+log.nb.active(n))/n))
}

smallest.block.scen3 <- function(n){
  return(sqrt(log(log(n))/n) + sqrt((log(0.5)+log.nb.active(n))/n))
}

smallest.block.scen4 <- function(n){
  return(sqrt((0.5+log(n))/n) + sqrt((log(0.5)+log.nb.active(n))/n))
}

ratio.scen1 <- function(n){
  return(smallest.block.scen1(n)/smallest.std.scen1(n))
}

ratio.scen2 <- function(n){
  return(smallest.block.scen2(n)/smallest.std.scen2(n))
}


ratio.scen3 <- function(n){
  return(smallest.block.scen3(n)/smallest.std.scen3(n))
} 


ratio.scen4 <- function(n){
  return(smallest.block.scen4(n)/smallest.std.scen4(n))
}

p.sc1<-ggplot(data.frame(n = seq(300,1.4*10^5,1000)),aes(x=n)) +
  geom_function(fun = smallest.std.scen1, aes(col='black'), linewidth = 0.5) +
  geom_function(fun = smallest.block.scen1, aes(col='darkgrey'), linewidth = 0.5) +
  scale_x_log10(breaks=c(300,10^3,10^4,10^5),labels = scales::comma) +
  scale_color_manual(values=c('black', 'darkgrey'),
                     labels=c(bquote({tilde(beta)^'*'}['min']),bquote({tilde(beta)^'*'}['min,2'])))+
  coord_cartesian(expand = FALSE, ylim=c(0,0.25))+
  theme_light(base_size = 9)+
  ylab(bquote({tilde(beta)^'*'}['min']~','~{tilde(beta)^'*'}['min,2']))+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,0,-5),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.key.spacing.x = unit(0, 'pt'),
        legend.key.spacing.y = unit(0, 'pt'),
        legend.key.width = unit(14, 'pt'),
        legend.key.height = unit(9, 'pt'),
        legend.text=element_text(size=9),
  )
p.sc1
ggsave(paste0(path,'p.betamin.sc1.ex.pdf'),width = 70, height = 60, units='mm')

p.sc2<-ggplot(data.frame(n = seq(300,1.4*10^5,1000)),aes(x=n)) +
  geom_function(fun = smallest.std.scen2, aes(col='black'), linewidth = 0.5) +
  geom_function(fun = smallest.block.scen2, aes(col='darkgrey'), linewidth = 0.5) +
  scale_x_log10(breaks=c(300,10^3,10^4,10^5),labels = scales::comma) +
  scale_color_manual(values=c('black', 'darkgrey'),
                     labels=c(bquote({tilde(beta)^'*'}['min']),bquote({tilde(beta)^'*'}['min,2'])))+
  coord_cartesian(expand = FALSE, ylim=c(0,0.4))+
  theme_light(base_size = 9)+
  ylab(bquote({tilde(beta)^'*'}['min']~','~{tilde(beta)^'*'}['min,2']))+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,0,-5),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.key.spacing.x = unit(0, 'pt'),
        legend.key.spacing.y = unit(0, 'pt'),
        legend.key.width = unit(14, 'pt'),
        legend.key.height = unit(9, 'pt'),
        legend.text=element_text(size=9),
  )
p.sc2
ggsave(paste0(path,'p.betamin.sc2.ex.pdf'),width = 70, height = 60, units='mm')

p.sc3<-ggplot(data.frame(n = seq(300,1.4*10^5,1000)),aes(x=n)) +
  geom_function(fun = smallest.std.scen3, aes(col='black'), linewidth = 0.5) +
  geom_function(fun = smallest.block.scen3, aes(col='darkgrey'), linewidth = 0.5) +
  scale_x_log10(breaks=c(300,10^3,10^4,10^5),labels = scales::comma) +
  scale_color_manual(values=c('black', 'darkgrey'),
                     labels=c(bquote({tilde(beta)^'*'}['min']),bquote({tilde(beta)^'*'}['min,2'])))+
  coord_cartesian(expand = FALSE, ylim=c(0,0.25))+
  theme_light(base_size = 9)+
  ylab(bquote({tilde(beta)^'*'}['min']~','~{tilde(beta)^'*'}['min,2']))+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,0,-5),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.key.spacing.x = unit(0, 'pt'),
        legend.key.spacing.y = unit(0, 'pt'),
        legend.key.width = unit(14, 'pt'),
        legend.key.height = unit(9, 'pt'),
        legend.text=element_text(size=9),
  )
p.sc3
ggsave(paste0(path,'p.betamin.sc3.ex.pdf'),width = 70, height = 60, units='mm')

p.sc4<-ggplot(data.frame(n = seq(300,1.4*10^5,1000)),aes(x=n)) +
  geom_function(fun = smallest.std.scen4, aes(col='black'), linewidth = 0.5) +
  geom_function(fun = smallest.block.scen4, aes(col='darkgrey'), linewidth = 0.5) +
  scale_x_log10(breaks=c(300,10^3,10^4,10^5),labels = scales::comma) +
  scale_color_manual(values=c('black', 'darkgrey'),
                     labels=c(bquote({tilde(beta)^'*'}['min']),bquote({tilde(beta)^'*'}['min,2'])))+
  coord_cartesian(expand = FALSE, ylim=c(0,0.25))+
  theme_light(base_size = 9)+
  ylab(bquote({tilde(beta)^'*'}['min']~','~{tilde(beta)^'*'}['min,2']))+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,0,-5),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.key.spacing.x = unit(0, 'pt'),
        legend.key.spacing.y = unit(0, 'pt'),
        legend.key.width = unit(14, 'pt'),
        legend.key.height = unit(9, 'pt'),
        legend.text=element_text(size=9),
  )
p.sc4
ggsave(paste0(path,'p.betamin.sc4.ex.pdf'),width = 70, height = 60, units='mm')


p1bis<-ggplot(data.frame(n = seq(100,10^6,1000)),aes(x=n)) +
  geom_function(fun = ratio.scen1, aes(linetype = "1_")) +
  geom_function(fun = ratio.scen2, aes(linetype = "2_")) +
  geom_function(fun = ratio.scen3, aes(linetype = "3_")) +
  geom_function(fun = ratio.scen4, aes(linetype = "4_")) +
  scale_x_log10(labels = scales::comma) +
  scale_linetype_manual('Example',values=c("1_"=1, "2_"=2,"3_"=3,"4_"=4),labels=c(1,2,3,4))+
  ylab('Ratio of smallest signals recoverable') +
  coord_cartesian(xlim = c(10^3,10^6),ylim=c(0,1))+
  theme_light(base_size = 8)+
  #ggtitle('Ratio of smallest signal recoverable')+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.margin=margin(0,0,0,0), 
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.position="bottom")
p1bis
ggsave(paste0('p.ratio.smallestbeta.pdf'), width = 60, height = 72, units='mm')



p1bis<-ggplot(data.frame(n = seq(100,10^6,1000)),aes(x=n)) +
  geom_function(fun = ratio.scen1, aes(linetype = "1_")) +
  geom_function(fun = ratio.scen2, aes(linetype = "2_")) +
  #geom_function(fun = ratio.scen3, aes(linetype = "3_")) +
  geom_function(fun = ratio.scen4, aes(linetype = "3_")) +
  scale_x_log10(labels = scales::comma) +
  scale_linetype_manual('Example',values=c("1_"=1, "2_"=2,"3_"=3),labels=c(1,2,3))+
  ylab('Ratio of smallest signals recoverable') +
  coord_cartesian(xlim = c(10^3,10^5),ylim=c(0,1))+
  theme_light(base_size = 12)+
  #ggtitle('Ratio of smallest signal recoverable')+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.margin=margin(0,0,0,0), 
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.position="bottom")
p1bis
ggsave(paste0('p.ratio.smallestbeta.posterA1.pdf'), width = 110, height = 100, units='mm')


p1bis<-ggplot(data.frame(n = seq(100,10^6,1000)),aes(x=n)) +
  geom_function(fun = ratio.scen1, aes(linetype = "1_")) +
  geom_function(fun = ratio.scen2, aes(linetype = "2_")) +
  #geom_function(fun = ratio.scen3, aes(linetype = "3_")) +
  geom_function(fun = ratio.scen4, aes(linetype = "3_")) +
  scale_x_log10(labels = scales::comma) +
  scale_linetype_manual('Example',values=c("1_"=1, "2_"=2,"3_"=3),labels=c(1,2,3))+
  ylab('Ratio of smallest signals recoverable') +
  coord_cartesian(xlim = c(10^3,10^5),ylim=c(0,1))+
  theme_light(base_size = 12)+
  #ggtitle('Ratio of smallest signal recoverable')+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.margin=margin(0,0,0,0), 
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.position="bottom")
p1bis
ggsave(paste0('p.ratio.smallestbeta.posterA0.pdf'), width = 150, height = 120, units='mm')



