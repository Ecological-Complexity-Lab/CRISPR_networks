library(pomp)
library(tidyverse)

on_Midway <- function(){ifelse(Sys.getenv('USER')=='pilosofs',T,F)}

if (on_Midway()){
  run <- Sys.getenv('SLURM_ARRAY_TASK_ID')
}

if (length(commandArgs(trailingOnly=TRUE))==0) {
  args <- c('12499', 'imm_density_124_231')
} else {
  args <- commandArgs(trailingOnly=TRUE)
}
seed <- args[1]
network_density_file <- args[2]
base_name <- paste('mu1e-7_initialDiffDp1_S10P15_R-',seed,sep='')
print(base_name)
print(network_density_file)

if(on_Midway()){
  M_seq <- read_csv(paste('/project2/pascualmm/CRISPR/mu1e-7_initialDiffDp1_S10P15_R-/',base_name,'/',base_name,'_',network_density_file,'.csv',sep=''))
  M_seq <- M_seq$d
  stop_time <- length(M_seq)
  print(stop_time)
}



# library(ggplot2)


CRISPR_changing_M <- function(t, M_t, N_t, V_t){
  #### this code follows Aaron King's tutorial on https://kingaa.github.io/clim-dis/stochsim/gillespie.html ####
  #### section: The Gillespie algorithm in pomp  ####
  
  ## 1. the model -------
  ## define a pomp object with statenames, parameters, rate matrix, 
  ## how statenames change, interdependence of equations
  ## time steps, and measurement models
  CRISPR <- pomp(data=data.frame(
    time=seq(0,1,by=0.2),
    reportN=NA, reportV=NA),
    times="time",t0=0,
    rprocess=gillespie.sim(
      v=cbind(
        bacteria_growth=c(1,0), # Events that refer to the events in rate.fun, in the same order. length and order of each vector corresponds to length and order of statenames
        bacteria_competition=c(-1,0), # death due to competition
        not_immune_interaction = c(-1,50),
        immune_interaction = c(-1,50),
        viral_infection = c(0,-1),
        viral_deactivation = c(0,-1)
      ),
      rate.fun=function(j,x,t,params,...) {
        # if(x['N']>params["K"]){print(x['N'])};
        switch(j, #events
               params["r"]*x["N"],
               (params["r"]*x["N"]^2)/params["K"],
               (1-params['q'])*(1-params['M'])*params['phi']*x['N']*x['V'], # Interaction of non-immune host
               params['p']*params['M']*params['phi']*x['N']*x['V'], # Immune interaction
               params['phi']*x['N']*x['V'], # Viral infection
               params['m']*x['V']) # Viral deatctivation
      }
    ),
    rmeasure=Csnippet("reportN=N;reportV=V;"), # Measurement model. If there is no measurement model, just put the state name (e.g., N).
    statenames=c("N","V"),
    paramnames=c("K","r","beta","phi","p","q","mu","m",'M')
  )
  # 2. parameters and starting values
  paramSets<-c(K=10^5.5,r=1,beta=50,phi=10^-7,p=10^-5,q=10^-5,mu=10^-7,m=0.1,M=M_t)
  # 3. running the actual simulations
  simState<-simulate(CRISPR,
                     seed = 123,
                     nsims=1, 
                     params=c(N.0 = N_t, V.0 = V_t, paramSets),as.data.frame=TRUE) #.0 MEANS STARTING VALUE
  simState$time <- simState$time+t
  return(simState)
}

G <- NULL
N_t = 2*10^5
V_t = 2*10^5
for (t in 1:stop_time){
  # M_seq[t] <- 0.1 # Fix the density to test
  M_t <- ifelse(M_seq[t]>runif(1),1,0)
  print(paste('Time: ',t,' | Immm. density: ',round(M_seq[t],5), ' | M_t: ',M_t,' | N_t:', N_t,' | V_t:', V_t,sep=''))
  x <- CRISPR_changing_M(t, M_t, N_t, V_t)
  N_t <- x$N[nrow(x)]
  V_t <- x$V[nrow(x)]
  print(paste('N_t+1:',N_t,'V_t+1:',V_t))
  x$M_t <- M_t
  x$density <- M_seq[t]
  x$run <- run
  x$seed <- seed
  x$discrete_time <- t
  G <- rbind(G, x[-c(2,3)])
  write_csv(x, paste('/project2/pascualmm/CRISPR/Gillespie/','G_',seed,'_',run,'_',t,'.csv',sep=''))
}
write_csv(G, paste('/project2/pascualmm/CRISPR/Gillespie/','G_all_time_series_',seed,'_',run,'.csv',sep=''))

# Skip plotting on Midway
if (!on_Midway()){
  library(plotly)
  
  G <- read_csv('/home/shai/Documents/CRISPR_data/G_all_time_series_12499_5.csv')  
  hr_seq <- 1:max(G$time)
  label_seq <- pretty(hr_seq, n=10)
  label_seq <- subset(label_seq, label_seq<max(hr_seq))
  plt_M <- ggplot(G, aes(time, density))+geom_line()+
      theme(axis.title.x = element_blank())+
      scale_x_continuous(breaks=label_seq)
  
  # Divide by 10^5 for the scales of the labels
  p <- ggplot(G, aes(x=time))
  p <- p+geom_line(aes(y=N/10^5, color='N'))+labs(y='N (*10^5)')
  p <- p+geom_line(aes(y=V/10^7, color='V'))
  p <- p + scale_y_continuous(sec.axis = sec_axis(~.*10^2, name = "V (*10^5)"))
  p <- p + scale_colour_manual(values = c("blue", "red"))+
    theme(legend.position = c(0.8,0.8), legend.title = element_blank())+
    labs(x="Time")
  p2 <- ggplot(G, aes(N/10^5,V/10^5))+geom_path()+
    # geom_point(color='purple')+
    geom_vline(xintercept = 10^5.5/10^5)+
    labs(x="N (*10^5)", y="V (*10^5)")
  
  png('/home/shai/Documents/CRISPR/Presentations/Gillespie_reduced_model_example.png', 1920, 1080, res=200)
  plot_grid(p,p2,align='vh',nrow = 1,ncol=2, scale = 0.9)
  dev.off()
  # toplot <- G %>% select(time, N, V) %>% gather(key='variable', value='value', -time)
  # plt_abundances <- toplot %>% ggplot(aes(time, y=value, color=variable))+
  #   geom_line()+
  #   facet_wrap(~variable, scales='free')+
  #   theme(legend.position = c(0.7,0.8), legend.title = element_blank(), legend.direction = 'horizontal')+
  #   scale_x_continuous(breaks=label_seq)

Fig <- plot_grid(plt_abundances, plt_M, nrow=2, align='vh')



# plotly
G <- G[-nrow(G),]
G$time_discrete <- floor(G$time)
# plot_colors <- tibble(time_discrete=unique(G$time_discrete), col=rep(gg_color_hue(300), each=6))
plot_colors <- tibble(time_discrete=unique(G$time_discrete), col=gg_color_hue(max(G$time_discrete)))
G %<>% left_join(plot_colors)
p <- plot_ly(G, x = ~time, y = ~V, z = ~N, type = 'scatter3d', mode = 'lines',
             opacity = 1, line = list(width = 6, color = ~col, reverscale = FALSE))
Sys.setenv("plotly_username"="shainova")
Sys.setenv("plotly_api_key"="DzHyp7FUFxL6vpd4e4fn")
chart_link = api_create(p, filename="Temporal_Gillespie_1")
chart_link

}