summarize.sim <- function(scenario,path){
  require(tidyverse)
  
  setwd(paste0(path,'/scenario',scenario))
  
  res.varn <- data.frame()
  time.varn <- data.frame()
  
  nb.active<- function(n){
    return(3*log(n))
  }
  
  nb.inactive<- function(n,scenario){
    if(scenario%in%c(1,2,3,4,5))  return(1.5*n)
  }
  
  
  for (n  in c(#20,40,
    60,80,seq(100,700,100))){
    temp<-read.csv(paste0('sim.result.scenario',scenario,'.n',n,".csv"))
    summary.sim.result <- temp %>% group_by(method) %>% 
      summarise(prob.recovery = mean(recovery),
                mean.size.sel =mean(est.size), 
                mean.FDR = mean(FDR),
                mean.power = mean(power),
                mean.mse = mean(est.mse))
    summary.sim.result$n <- rep(n, nrow(summary.sim.result))
    res.varn<- rbind(res.varn,summary.sim.result)
    
    summary.sim.time <- temp %>% summarise(median.time = mean(time.l0))
    summary.sim.time$n <- rep(n, nrow(summary.sim.time))
    summary.sim.time$p <- rep(nb.inactive(n,scenario)+nb.active(n), nrow(summary.sim.time))
    
    time.varn <- rbind(time.varn,summary.sim.time)
  }
  
  method.sel <- c(#"S.EB.b", 
    "S.A.b", 
    #"S.EB", 
    "S.A",
    'EBIC',"lasso.cv", "scad.cv")
  
  res.varn <- res.varn %>% 
    filter(method%in%method.sel)
  
  res.varn$method <- factor(res.varn$method, levels = method.sel)
  
  
  #### Prob recovery ####
  res.varn %>%
    ggplot() + 
    geom_line(aes(x=n, y=prob.recovery, col=method, linetype=method), size=0.5) + 
    # coord_cartesian(ylim=c(0,1)) +
    geom_point(aes(x=n,y=prob.recovery,shape=method),size=1.5) +
    theme_light(base_size = 8)+
    scale_color_manual(values=c(#'black',
      'black',
      #'darkgrey',
      'darkgrey',
      'lightgrey','lightgrey', 'lightgrey'),
      labels=c(#bquote(hat(S)^'EB,b'),
        'Tran-s-ell0',
        #bquote(hat(S)^'EB'),
        bquote(hat(S)^'A'),
        'EBIC','LASSO', 'SCAD'))+
    scale_linetype_manual(values=c(#1,
      1,
      #1,
      1,
      2,1,1), 
      labels=c(#bquote(hat(S)^'EB,b'),
        'Tran-s-ell0',
        #bquote(hat(S)^'EB'), 
        bquote(hat(S)^'A'),
        'EBIC','LASSO', 'SCAD'))+
    scale_shape_manual(values=c(#NA,
      NA,
      #NA,
      NA,
      NA,3,4), 
      labels=c(#bquote(hat(S)^'EB,b'),
        'Tran-s-ell0',
        #bquote(hat(S)^'EB'),
        bquote(hat(S)^'A'),
        'EBIC','LASSO', 'SCAD'))+
    ylab('Prob. recovery') +
    scale_x_continuous(labels = ~paste0('(',.,',',round(nb.inactive(.,scenario)+nb.active(.)), ")"), name = "(n,p)") +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.margin=margin(0,0,0,0), 
          legend.box.margin=margin(-5,-10,-5,-30),
          legend.position="bottom",
          legend.title = element_blank(),
          legend.key.spacing.x = unit(0, 'pt'),
          legend.key.spacing.y = unit(0, 'pt'),
          legend.key.width = unit(14, 'pt'),
          legend.key.height = unit(9, 'pt'),
          legend.text=element_text(size=8),
    ) +
    guides(color=guide_legend(nrow=2,byrow=TRUE))
  
  ggsave(paste0('p.TL.probrecov.ex',scenario,'.pdf'),width = 52, height = 60, units='mm')
  
  
  
  # res.varn %>% 
  #   filter(method%in%c('baseline1','method1.1', 'EBIC')) %>%
  #   ggplot() + 
  #   geom_line(aes(x=n, y=prob.recovery, col=method, linetype=method), size=0.7) + 
  #   coord_cartesian(ylim=c(0,1)) +
  #   theme_light(base_size = 12)+
  #   scale_color_manual(values=c('black','grey','grey'),
  #                      labels=c(bquote(hat(S)^'EB,b'),bquote(hat(S)^'EB'),'EBIC'))+
  #   scale_linetype_manual(values=c(1,4,5), 
  #                         labels=c(bquote(hat(S)^'EB,b'),bquote(hat(S)^'EB'),'EBIC'))+
  #   ylab('Prob. recovery') +
  #   theme(panel.grid.minor.x = element_blank(),
  #         panel.grid.minor.y = element_blank(),
  #         legend.margin=margin(0,0,0,0), 
  #         legend.box.margin=margin(-5,-10,-5,-30),
  #         legend.position="bottom",
  #         legend.title = element_blank(),
  #         legend.key.spacing.x = unit(0, 'pt'),
  #         legend.key.spacing.y = unit(0, 'pt'),
  #         legend.key.width = unit(10, 'pt'),
  #         legend.text=element_text(size=12)
  #   ) +
  #   guides(color=guide_legend(nrow=1,byrow=TRUE))
  # 
  # ggsave(paste0('p.TL.probrecov.ex1.slogn.posterA1.pdf'),width = 110, height = 100, units='mm')
  # 
  # 
  # res.varn %>% 
  #   filter(method%in%c('baseline1','method1.1', 'EBIC')) %>%
  #   ggplot() + 
  #   geom_line(aes(x=n, y=prob.recovery, col=method, linetype=method), size=0.7) + 
  #   coord_cartesian(ylim=c(0,1)) +
  #   theme_light(base_size = 12)+
  #   scale_color_manual(values=c('black','grey','grey'),
  #                      labels=c(bquote(hat(S)^'EB,b'),bquote(hat(S)^'EB'),'EBIC'))+
  #   scale_linetype_manual(values=c(1,4,5), 
  #                         labels=c(bquote(hat(S)^'EB,b'),bquote(hat(S)^'EB'),'EBIC'))+
  #   ylab('Prob. recovery') +
  #   theme(panel.grid.minor.x = element_blank(),
  #         panel.grid.minor.y = element_blank(),
  #         legend.margin=margin(0,0,0,0), 
  #         legend.box.margin=margin(-5,-10,-5,-30),
  #         legend.position="bottom",
  #         legend.title = element_blank(),
  #         legend.key.spacing.x = unit(0, 'pt'),
  #         legend.key.spacing.y = unit(0, 'pt'),
  #         legend.key.width = unit(10, 'pt'),
  #         legend.text=element_text(size=12)
  #   ) +
  #   guides(color=guide_legend(nrow=1,byrow=TRUE))
  # 
  # ggsave(paste0('p.TL.probrecov.ex1.slogn.posterA0.pdf'),width = 150, height = 120, units='mm')
  
  
  
  #### FDR ####
  
  res.varn %>%
    ggplot() + 
    geom_line(aes(x=n, y=mean.FDR, col=method, linetype=method), size=0.5) + 
    geom_point(aes(x=n,y=mean.FDR,shape=method),size=1.5) +
    # coord_cartesian(ylim=c(0,1)) +
    theme_light(base_size = 8)+
    scale_color_manual(values=c(#'black',
      'black',
      #'darkgrey',
      'darkgrey',
      'lightgrey','lightgrey', 'lightgrey'),
      labels=c(#bquote(hat(S)^'EB,b'),
        'Tran-s-ell0',
        #bquote(hat(S)^'EB'),
        bquote(hat(S)^'A'),
        'EBIC','LASSO', 'SCAD'))+
    scale_linetype_manual(values=c(#1,
      1,
      #1,
      1,
      2,1,1), 
      labels=c(#bquote(hat(S)^'EB,b'),
        'Tran-s-ell0',
        #bquote(hat(S)^'EB'), 
        bquote(hat(S)^'A'),
        'EBIC','LASSO', 'SCAD'))+
    scale_shape_manual(values=c(#NA,
      NA,
      #NA,
      NA,
      NA,3,4), 
      labels=c(#bquote(hat(S)^'EB,b'),
        'Tran-s-ell0',
        #bquote(hat(S)^'EB'),
        bquote(hat(S)^'A'),
        'EBIC','LASSO', 'SCAD'))+
    ylab('FDR') +
    scale_x_continuous(labels = ~paste0('(',.,',',round(nb.inactive(.,scenario)+nb.active(.)), ")"), name = "(n,p)") +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.margin=margin(0,0,0,0), 
          legend.box.margin=margin(-5,-10,-5,-30),
          legend.position="bottom",
          legend.title = element_blank(),
          legend.key.spacing.x = unit(0, 'pt'),
          legend.key.spacing.y = unit(0, 'pt'),
          legend.key.width = unit(14, 'pt'),          
          legend.key.height = unit(9, 'pt'),           
          legend.text=element_text(size=8),
    ) +
    guides(color=guide_legend(nrow=2,byrow=TRUE))
  
  ggsave(paste0('p.TL.FDR.ex',scenario,'.pdf'),width = 52, height = 60, units='mm')
  
  
  #### Power ####
  
  res.varn %>%
    ggplot() + 
    geom_line(aes(x=n, y=mean.power, col=method, linetype=method), size=0.5) + 
    geom_point(aes(x=n,y=mean.power,shape=method),size=1.5) +
    # coord_cartesian(ylim=c(0,1)) +
    theme_light(base_size = 8)+
    scale_color_manual(values=c(#'black',
      'black',
      #'darkgrey',
      'darkgrey',
      'lightgrey','lightgrey', 'lightgrey'),
      labels=c(#bquote(hat(S)^'EB,b'),
        'Tran-s-ell0',
        #bquote(hat(S)^'EB'),
        bquote(hat(S)^'A'),
        'EBIC','LASSO', 'SCAD'))+
    scale_linetype_manual(values=c(#1,
      1,
      #1,
      1,
      2,1,1), 
      labels=c(#bquote(hat(S)^'EB,b'),
        'Tran-s-ell0',
        #bquote(hat(S)^'EB'), 
        bquote(hat(S)^'A'),
        'EBIC','LASSO', 'SCAD'))+
    scale_shape_manual(values=c(#NA,
      NA,
      #NA,
      NA,
      NA,3,4), 
      labels=c(#bquote(hat(S)^'EB,b'),
        'Tran-s-ell0',
        #bquote(hat(S)^'EB'),
        bquote(hat(S)^'A'),
        'EBIC','LASSO', 'SCAD'))+
    ylab('Power') +
    scale_x_continuous(labels = ~paste0('(',.,',',round(nb.inactive(.,scenario)+nb.active(.)), ")"), name = "(n,p)") +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.margin=margin(0,0,0,0), 
          legend.box.margin=margin(-5,-10,-5,-30),
          legend.position="bottom",
          legend.title = element_blank(),
          legend.key.spacing.x = unit(0, 'pt'),
          legend.key.spacing.y = unit(0, 'pt'),
          legend.key.width = unit(14, 'pt'),
          legend.key.height = unit(9, 'pt'),          
          legend.text=element_text(size=8),
    ) +
    guides(color=guide_legend(nrow=2,byrow=TRUE))
  
  ggsave(paste0('p.TL.power.ex',scenario,'.pdf'),width = 52, height = 60, units='mm')
  
  #### Estimation mse ####
  
  res.varn %>%
    ggplot() + 
    geom_line(aes(x=n, y=mean.mse, col=method, linetype=method), size=0.5) + 
    geom_point(aes(x=n,y=mean.mse,shape=method),size=1.5) +
    # coord_cartesian(xlim=c(60,700)) +
    theme_light(base_size = 8)+
    scale_y_log10() +
    scale_color_manual(values=c(#'black',
      'black',
      #'darkgrey',
      'darkgrey',
      'lightgrey','lightgrey', 'lightgrey'),
      labels=c(#bquote(hat(S)^'EB,b'),                
        'Tran-s-ell0',                
        #bquote(hat(S)^'EB'),                
        bquote(hat(S)^'A'),                
        'EBIC','LASSO', 'SCAD'))+
    scale_linetype_manual(values=c(#1,
      1,
      #1,
      1,
      2,1,1), 
      labels=c(#bquote(hat(S)^'EB,b'),
        'Tran-s-ell0',
        #bquote(hat(S)^'EB'), 
        bquote(hat(S)^'A'),
        'EBIC','LASSO', 'SCAD'))+
    scale_shape_manual(values=c(#NA,
      NA,
      #NA,
      NA,
      NA,3,4), 
      labels=c(#bquote(hat(S)^'EB,b'),
        'Tran-s-ell0',
        #bquote(hat(S)^'EB'),
        bquote(hat(S)^'A'),
        'EBIC','LASSO', 'SCAD'))+
    ylab('MSE') +
    scale_x_continuous(labels = ~paste0('(',.,',',round(nb.inactive(.,scenario)+nb.active(.)), ")"), name = "(n,p)") +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.margin=margin(0,0,0,0), 
          legend.box.margin=margin(-5,-10,-5,-30),
          legend.position="bottom",
          legend.title = element_blank(),
          legend.key.spacing.x = unit(0, 'pt'),
          legend.key.spacing.y = unit(0, 'pt'),
          legend.key.width = unit(14, 'pt'),
          legend.key.height = unit(9, 'pt'),           
          legend.text=element_text(size=8),
    ) +
    guides(color=guide_legend(nrow=2,byrow=TRUE))
  
  ggsave(paste0('p.TL.est.mse.ex',scenario,'.pdf'),width = 52, height = 60, units='mm')
  
  
  #### L0 sol computation time ####
  
  #if(sceario==1){y.lab <- paste0('p  (=',bquote(frac( .(3), .(2)) ),'(n+ln(n)))')  }
  
  time.varn %>% 
    ggplot() + 
    geom_line(aes(x=n, y=median.time), size=0.5) + 
    # coord_cartesian(ylim=c(0,1)) +
    theme_light(base_size = 8)+
    scale_x_continuous(labels = ~paste0('(',.,',',round(nb.inactive(.,scenario)+nb.active(.)), ")"), name = "(n,p)") +
    ylab(expression(paste('Computation time for informed \u2113'[0], ' in secs') )) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())
  
  ggsave(paste0('p.TL.time.ex',scenario,'.pdf'),width = 52, height = 60, units='mm')
}

library(tidyverse)

path = "C:/Users/Rognon/Documents/GitHub/varsel.extdata/simulations/transfer_learning/"
summarize.sim(1,path)
summarize.sim(2,path)
summarize.sim(3,path)
summarize.sim(4,path)
summarize.sim(5,path)



