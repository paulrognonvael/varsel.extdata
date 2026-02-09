# libraries 
####### libraries #######

library(ggplot2)
library(scales)
library(latex2exp)
library(gridExtra)
library(tidyverse)


path = '~/GitHub/varsel.extdata/WIP/paper/code.figures/'

####### Plot smallest betamin allowed #######


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