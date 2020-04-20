# !diagnostics off

# parameters of the model; used by the functions, so need to initialize here
beta=50;phi=10^-7;p=10^-5;q=10^-5;m=0.1;spacer_len=10

# A function to identify if analysis is run on the UChicago's Midway HPC. Change it to fit your HPC.
on_Midway <- function(){ifelse(Sys.getenv('USER')=='pilosofs',T,F)}

# If run with an sbatch pipeline take arguments from there. Otherwise those specified here.
if (length(commandArgs(trailingOnly=TRUE))==0) {
  # args <- c('mu5e-7_initialDiffDp1_S50P15_R-13997','5*10^-7','15',F, F)
  args <- c('mu1e-7_initialDiffDp1_S10P15_R-12499','1*10^-7','15',F, T)
} else {
  args <- commandArgs(trailingOnly=TRUE)
}
base_name <- args[1]
mu <- eval(parse(text = args[2]))
protospacer_len <- eval(parse(text = args[3]))
make_plots <- as.logical(args[4])

if(on_Midway()){system('module load gcc/6.1')}
if(!on_Midway()){setwd(paste('data/',base_name,sep=''))}

dir.create('figures')

library(ggtree)
library(ape)
library(treeio)
library(igraph)
library(tidyverse)
library(magrittr)
library(bipartite)
library(cowplot)
library(grid)
library(infomapecology)

if(check_infomap()==F){install_infomap()} # Install infomap

# Functions ---------------------------------------------------------------
prep.packages <- function(package.list, verbose=T) {
  loaded = package.list %in% .packages()
  if ( all(loaded) ) return(invisible())
  
  package.list = package.list[!loaded]
  installed = package.list %in% .packages(TRUE)
  if ( !all(installed) ) install.packages(package.list[!installed], repos="http://cran.rstudio.com/")
  for ( p in package.list ){
    print(paste("Loading package:",p))
    if (verbose){
      library(p,character.only=TRUE)
    } else {
      suppressMessages(library(p,character.only=TRUE))  
    }
  }
}

returnnull <- function(x) if (is.null(x)){'none'} else {x}

notify <- function(x){
  print(paste('[',Sys.time(),'] ',x,sep=''))
}

record_data <- function(x){
  notify(paste('Recording ',deparse(substitute(x)),sep=''))
  write_csv(x, paste(base_name,'_',deparse(substitute(x)),'.csv',sep=''), col_names = T)
}

make_name <- function(x,hr){
  hr_str <- str_pad(hr, width = 4, side = 'left', pad = '0')
  paste('data/',x,'_',hr_str,'.csv',sep='')
}

make_png <- function(p, method='ggsave'){
  if (method=='ggsave'){
    ggplot2::ggsave(paste('figures/',base_name,'_',deparse(substitute(p)),'.png',sep=''), p, device = 'png', width = 32, height = 18, units = 'cm', dpi = 200)
  } else {
    png(paste('figures/',base_name,'_',deparse(substitute(p)),'.png',sep=''),1920,1080,res=150)
    print(p)
    dev.off()
  }
}

make_svg <- function(p){
  svg(paste('figures/',base_name,'_',deparse(substitute(p)),'.svg',sep=''),12.8,8)
  print(p)
  dev.off()
}

make_png_svg <- function(p){
  make_png(p)
  make_svg(p)
}

standard_plot <- function(p){
  p+scale_x_continuous(breaks=label_seq)+
    geom_vline(xintercept=BDRs$start, col='#27AE60',size=1.2)+
    geom_vline(xintercept=VDRs$start, col='purple',size=1.2)+
    labs(x='Time')
}

gg_color_hue <- function(n, hue_min = 10, hue_max = 280, tune1 = 62, tune2 = 100) {
  hues = seq(hue_min, hue_max, length=n+1)
  hcl(h=hues, l=tune1, c=tune2)[1:n]
}

list_to_matrix <- function(l, directed=F, bipartite=T){
  # Only deal with first 3 columns.
  g <- graph.data.frame(l,directed = directed)
  if(bipartite){
    V(g)$type <- V(g)$name %in% as.data.frame(l)[,1] # As.data.frame is necessary because if l is a tibble then l[,1] does not work
    output_mat <- as_incidence_matrix(g, names = T, attr = 'w', sparse = F)
  } else {
    output_mat <- as_adjacency_matrix(g, names = T, sparse = F, attr = 'w')
  }
  # print(dim(output_mat))
  # if(any(rowSums(output_mat)==0)){stop('One or more rows sum to 0')}
  # if(any(colSums(output_mat)==0)){stop('One or more columns sum to 0')}
  return(output_mat)
}


create_networks_hr <- function(virus_data,bacteria_data,hr){
  virus_data_hr <- virus_data %>% filter(timesOfRecord==hr)
  bacteria_data_hr <- bacteria_data %>% filter(timesOfRecord==hr)
  virus_abund_hr <- virus_data_hr %>% select(label, density) %>% rename(V_ID=label) %>% mutate(V_ID=paste('V_',str_pad(V_ID, 4, 'left', '0'),sep=''))
  bacteria_abund_hr <- bacteria_data_hr %>% select(label, density) %>% rename(B_ID=label) %>% mutate(B_ID=paste('B_',str_pad(B_ID, 4, 'left', '0'),sep=''))
  
  # Create networks in a list form
  
  ## Virus-protospacer network
  virus_ps_hr_list <- virus_data_hr %>% select(-time, -timesOfRecord, -density) %>% 
    gather(key='PSidx', value='value', -label) %>% 
    rename(V_ID=label, PS=value) %>% 
    arrange(V_ID,PS) %>% 
    mutate(w=1) %>% 
    mutate(V_ID=paste('V_',str_pad(V_ID, 4, 'left', '0'),sep='')) %>% 
    mutate(PS=paste('P_',str_pad(PS, 4, 'left', '0'),sep='')) %>% 
    select(V_ID,PS,w,PSidx)
  
  ## Bacteria-spacer network
  bacteria_sp_hr_list <- bacteria_data_hr %>% 
    select(-time, -timesOfRecord, -density) %>% 
    gather(key='SPidx', value='value', -label) %>% 
    filter(value!=-1) %>%
    rename(B_ID=label, SP=value) %>% 
    arrange(B_ID,SP) %>% 
    mutate(w=1) %>%
    # mutate(w=ifelse(SP==-1,0,1)) %>% 
    mutate(B_ID=paste('B_',str_pad(B_ID, 4, 'left', '0'),sep='')) %>% 
    mutate(SP=paste('P_',str_pad(SP, 4, 'left', '0'),sep='')) %>% 
    select(B_ID,SP,w,SPidx)
 
  # Bacteria without spacers: This is needed for NEUTRAL SCENARIOS in which
  # interactions are not related to spacer-protospacer acquisition. So bacteria
  # can actually survive without acquiting protospacers.
  bacteria_no_sp <- bacteria_data_hr %>% 
    select(-time, -timesOfRecord, -density) %>% 
    gather(key='SPidx', value='value', -label) %>% 
    filter(value==-1) %>% 
    group_by(label) %>% 
    summarise(no_spacers=n()) %>% filter(no_spacers==spacer_len) %>%
    rename(B_ID=label) %>% 
    mutate(B_ID=paste('B_',str_pad(B_ID, 4, 'left', '0'),sep=''))
  
  # Add to the bacteria-spacer network those bacteria without spacers.
  suppressWarnings(
    bacteria_sp_hr_list <- bind_rows(bacteria_sp_hr_list, 
                                     as.tibble(expand.grid(B_ID=bacteria_no_sp$B_ID, SP=unique(bacteria_sp_hr_list$SP), w=0)))
  )
  ## Immunity network as a list
  x <- virus_ps_hr_list %>% select(V_ID, PS)
  y <- bacteria_sp_hr_list %>% 
    filter(w!=0) %>% # MUST include only the interactions usign filter(w!=0), otherwise the immunity network will have false edges that should not exist.
    select(B_ID, SP) 
  immunity_list <- x %>% inner_join(y, by = c('PS'='SP')) %>% 
    arrange(B_ID, V_ID) %>% 
    group_by(V_ID, B_ID) %>% count() %>%
    rename(w=n)
  
  # Define nodes in the system
  viruses_hr <- sort(unique(virus_abund_hr$V_ID))
  bacteria_hr <- sort(unique(bacteria_abund_hr$B_ID))
  spacers_hr <- sort(union(bacteria_sp_hr_list$SP, virus_ps_hr_list$PS))
  nodes <- tibble(nodeID=1:(length(viruses_hr)+length(bacteria_hr)+length(spacers_hr)),
                  nodeName=c(viruses_hr, bacteria_hr, spacers_hr),
                  type=c(rep(1,length(viruses_hr)), rep(2,length(bacteria_hr)), rep(3,length(spacers_hr)))
  )
  nodes$type <- as.integer(nodes$type)
  
  # Check networks
  if (nrow(bacteria_sp_hr_list)==0) {bacteria_sp_hr_list <- NULL}
  if (nrow(immunity_list)==0) {immunity_list <- NULL}
  
  ## Infection network
  # To get this we first need to create the immunity nertwork as a matrix
  if (!is.null(immunity_list)){
    immunity_matrix <- list_to_matrix(immunity_list) # transform to a matrix format
  } else {
    immunity_matrix <- matrix(0, nrow=length(bacteria_hr), ncol=length(viruses_hr), dimnames = list(bacteria_hr, viruses_hr))
  }

  # Add all the nodes which were not included in the list because they have degree of zero
  # missing_bacteria <- setdiff(unique(bacteria_sp_hr_list$B_ID),rownames(immunity_matrix)) 
  # immunity_matrix <- rbind(immunity_matrix, matrix(0, 
  #                                                    ncol = ncol(immunity_matrix),
  #                                                    nrow = length(missing_bacteria),
  #                                                    dimnames = list(missing_bacteria, colnames(immunity_matrix))))
  # 
  
  infection_matrix <- 1*(immunity_matrix==0) # Binary infection network is the "negative (photography-wise)" of the immunity network
  if (any(infection_matrix==1)){ # If there is at least one non-zero interaction in the infection network
    N_T <- sum(bacteria_abund_hr$density)
    # This will produce a matrix with values: (N_i*V_j)/N_T
    A <- crossprod(matrix(bacteria_abund_hr$density, nrow=1,
                          ncol=length(bacteria_abund_hr$density),
                          dimnames=list(NULL, bacteria_abund_hr$B_ID)),
                   matrix(virus_abund_hr$density, 
                          nrow=1,ncol=length(virus_abund_hr$density)                        ,
                          dimnames=list(NULL, virus_abund_hr$V_ID)))
    A <- A/N_T
    # This will remove all cells that are 0.
    A <- A[rownames(infection_matrix),colnames(infection_matrix)]
    infection_matrix <- infection_matrix*A 
    # Transform the infection network to a list
    g <- graph.incidence(t(infection_matrix), directed = F, weighted = T) # Need to transpose so the "from" would be viruses and the "to" would be bacteria
    infection_list <- as.tibble(igraph::as_data_frame(g, what = 'edges'))
    names(infection_list) <- c('V_ID', 'B_ID', 'w')
  } else { # All interactions are 0
    infection_list <- NULL
  }
  
  # Calculate the effective mutation matrix. These lines here save a double for loop
  # and selecting bacteria and virus abundances at each step.
  # First produce a matrix with N_i*V_j values:
  A <- crossprod(matrix(bacteria_abund_hr$density, nrow=1,  
                        ncol=length(bacteria_abund_hr$density),
                        dimnames=list(NULL, bacteria_abund_hr$B_ID)),
                 matrix(virus_abund_hr$density, 
                        nrow=1,ncol=length(virus_abund_hr$density)                        ,
                        dimnames=list(NULL, virus_abund_hr$V_ID)))
  A <- A[rownames(immunity_matrix),colnames(immunity_matrix)]
  M0 <- 1*(immunity_matrix==0)
  M0 <- M0*A*beta*mu*phi*(1-q) # This is from Childs et al 2012 Suppl Info, page 4. The multiplication by M0 thakes only the no matches
  M1 <- 1*(immunity_matrix!=0)
  M1 <- M1*A*beta*mu*phi*p # This is from Childs et al 2012 Suppl Info, page 4. The multiplication by M1 thakes only the  matches
  mutation_matrix <- M0+M1

  return(list(hr=hr,
              virus_ps_list=virus_ps_hr_list,
              bacteria_sp_list=bacteria_sp_hr_list,
              immunity_list=immunity_list,
              immunity_matrix=immunity_matrix,
              infection_list=infection_list,
              infection_matrix=infection_matrix,
              mutation_matrix=mutation_matrix,
              bacteria_no_sp=bacteria_no_sp,
              virus_abund_hr=virus_abund_hr,
              bacteria_abund_hr=bacteria_abund_hr,
              nodes=nodes
  ))
}

get_regimes <- function (phage_time_series, d2_threshold=0.001, do_smoothing = T, make_plots=F) {
  if (do_smoothing){print('Finding regimes with smoothing...')} else {print('Finding regimes without smoothing...')}
  # Need to use the relative abundance of phages, because the total abundaces
  # create huge numbers in the derivatives.
  x <- phage_time_series %>% 
    group_by(timesOfRecord) %>% 
    summarise(a=sum(Pdensity)) %>% 
    rename(time=timesOfRecord) %>% 
    mutate(rel_abund=a/max(a, na.rm = T))
  # Find 1st and 2nd derivatives of the time series
  d1 <- diff(x$rel_abund)/diff(x$time) 
  d1_df <- tibble(time=1:length(d1), d1=d1)
  d2 <- diff(d1_df$d1)/diff(d1_df$time) 
  d2_df <- tibble(time=1:length(d2), d2=d2)
  # Define time series points as regimes. A point is defined in a regime if its 2nd derivative is lower than a threshold, d2_threshold
  suppressMessages(
  regimes <- x %>% left_join(d1_df) %>% left_join(d2_df) %>% 
    mutate(regime=ifelse(abs(d2)<=d2_threshold,'BDR','VDR'))
  )
  # Find sequences of regimes/no regimes.
  regime_seq <- rle(regimes$regime)
  regime_seq <- tibble(len=regime_seq$lengths, regime=regime_seq$values) # len is the length of the sequence
  # This allows to see points where dynamics switch from regimes to no regimes
  regime_seq$switch=cumsum(regime_seq$len) # cumsum gives the point in the time series
  
  # Smoothing eliminates small blips where there is a few points of 'VDR'
  # surrounded by many 'BDR'. This threshold is decided by the distribution of
  # the lengths of no regime sequences. Sequences with length smaller than that
  # threshold are defined as yes.
  if (do_smoothing){
    smoothing_threshold <- quantile(subset(regime_seq, regime=='VDR')$len, 0.75)
    smoothing <- subset(regime_seq, len<=smoothing_threshold & regime=='VDR') %>% mutate(start=switch-len)
    # smoothing$start[1] <- 1
    for (i in 1:nrow(smoothing)){
      # print(smoothing$start[i]:smoothing$switch[i])
      regimes$regime[smoothing$start[i]:smoothing$switch[i]] <- 'BDR'
    }
    # After smoothing need to define regime sequences again
    regime_seq <- rle(regimes$regime)
    regime_seq <- tibble(len=regime_seq$lengths, regime=regime_seq$values)
    regime_seq$switch=cumsum(regime_seq$len)
  }
  # To be a regime, a sequences has to be larger than the largest no-regime section,
  # which is the virus outbreak with the largest duration.
  max_virus_cycle <- max(subset(regime_seq, regime=='VDR')$len)
  
  regime_end <- which(regime_seq$len>max_virus_cycle)
  regime_start <- which(regime_seq$len>max_virus_cycle)-1
  regime_end <- regime_seq[regime_end,]$switch
  regime_start <- regime_seq[regime_start,]$switch
  if (length(regime_start)<length(regime_end)){
    regime_end <- regime_end[-1]
  }
  
  regimes_df <- rbind(
    tibble(start=1, end=regime_start[1]-1, regime_type='VDR'), # Add the first no regime period
    tibble(start=regime_start, end=regime_end, regime_type='BDR'),
    tibble(start=regime_end+1, end=c(regime_start-1, stop_time)[-1], regime_type='VDR')
  )
  regimes_df$duration=regimes_df$end-regimes_df$start
  regimes_df %<>% arrange(start)
  
  # Plot some diagnostics
  if (make_plots){
    plt_diagnostics <- 
      regimes %>% 
      gather(key='key', value='value', -time, -regime) %>% 
      mutate(key=factor(key,levels = c('a','rel_abund','d1','d2'))) %>% 
      ggplot()+
      geom_line(aes(time, value))+
      geom_point(aes(time, value, color=regime), size=0.5)+
      scale_color_manual(values = c('#27AE60','purple'))+
      scale_x_continuous(breaks = seq(0,max(x$time),500))+
      # scale_x_continuous(limits=c(800,1000))+
      facet_wrap(~key, scales='free', labeller = as_labeller(c(`a` = "Virus abundance",
                                                               `d1` = "First derivative",
                                                               `d2` = "Second derivative",
                                                               `rel_abund` = "Virus relative abundance")))+
      geom_vline(xintercept=regime_start, col='#27AE60')+
      geom_vline(xintercept=regime_end, col='purple')+
      theme(legend.position = 'none')+labs(x='Time')
  } else {
    plt_diagnostics <- NULL
  }
  return(list(regimes_df=regimes_df,
              plt_diagnostics=plt_diagnostics))
}



print_pattern <- function(parent, child, parent_death, child_birth, child_death, node, i){
  sprintf("(%s:%f, %s:%f)%s",
          parent, parent_death-child_birth,
          child, child_death-child_birth,
          paste(node, i , sep = "_"))
}

nodes_dataframe_to_one_root <- function(nodes, parent, children, parent_death) {
  for(i in 1:nrow(children)) {
    newChildren <- nodes %>% filter(parent_id == children$id[i]) %>% arrange(desc(creation_time))
    if (i == 1) {
      tempParent<-parent
      parent_death <- parent_death
    }else{
      tempParent<-out
      parent_death<-children$creation_time[i-1]
    }
    if (nrow(newChildren) > 0) {
      child <- nodes_dataframe_to_one_root(nodes, children$id[i], newChildren, children$death[i])
      child_birth <- children$creation_time[i]
      child_death <- newChildren$creation_time[nrow(newChildren)]
    }else{
      child <- children$id[i]
      child_birth <- children$creation_time[i]
      child_death <- children$death[i]
    }
    #print(list(tempParent, child, parent_death, child_birth, child_death))
    out<-print_pattern(tempParent, child, parent_death, child_birth, child_death, parent, nrow(children)-i)
  } 
  return(out)
}

nodes_dataframe_to_newick <- function(nodes) {
  root <- nodes %>% filter(is.na(parent_id))
  
  stopifnot(nrow(root) == 1)
  
  children <- nodes %>% filter(parent_id == root$id) %>% arrange(desc(creation_time))
  head(children)
  out<-nodes_dataframe_to_one_root(nodes, root$id[1], children, root$death[1])
  return(paste(out,  ":", children$creation_time[nrow(children)] - root$creation_time[1], ";",sep = ""))
}


# Function to test for phylogenetic signal in modules
test_PD_modules<- function(tree, module_object, node_start_letter){
  # Phylogenetic signal analysis
  D <- ape::cophenetic.phylo(tree) # Phyloegentic distance
  D <- matrix_to_list_unipartite(D, directed = T) # Use directed to make sure that the from column has all the nodes (need it for joining later)
  D <- D$edge_list
  
  # Difference between tree and matrix
  nodes_in_modules <- module_object$modules %>%
    filter(str_starts(node_name, node_start_letter)) %>%
    distinct(node_name) %>%
    mutate(node_name=str_replace_all(node_name, pattern = '\\.', ''))
  nodes_in_modules <- nodes_in_modules$node_name
  nodes_in_tree <- tree$tip.label
  # print(setdiff(nodes_in_modules, nodes_in_tree)) # In modules but not in tree
  # print(setdiff(nodes_in_tree, nodes_in_modules)) # In tree but not in modules
  # Overlapping nodes:
  overlapping <- intersect(nodes_in_tree, nodes_in_modules)
  
  # Observed modules
  M_obs <- module_object$modules %>%
    filter(str_starts(node_name, node_start_letter)) %>%
    mutate(node_name=str_replace_all(node_name, pattern = '\\.', '')) %>%
    filter(node_name %in% overlapping) %>%
    rename(m=module_level1) %>%
    select(node_name, m)
  
  #Mean PDistance between hosts within modules
  D_obs <- M_obs %>%
    inner_join(D, by=c('node_name'='from')) %>% # join PD distances
    rename(d=weight) %>%
    arrange(m, node_name) %>%
    group_by(m) %>% # Per module
    filter(to %in% node_name) %>% #Host pairs within a module
    summarise(d_mean=mean(d), mod_size=n())
  D_obs_mean <- mean(D_obs$d_mean)
  
  # print('Observed network:')
  # print(D_obs)
  
  #Shuffle to create permuted modules of the same size,
  #and recalculate the meand PD within modules. The shuffling permutes the ID of the strains.
  D_perm <- NULL
  nperm <- 500
  for (i in 1:nperm){
    # print(i)
    D_perm %<>% bind_rows(
      M_obs %>%
        mutate(node_name=sample(node_name, replace = F)) %>%
        inner_join(D, by=c('node_name'='from')) %>% # join PD distances
        rename(d=weight) %>%
        arrange(m, node_name) %>%
        group_by(m) %>% # Per module
        filter(to %in% node_name) %>% #Host pairs within a module
        summarise(d_mean=mean(d)) %>% # Calculate mean PD within modules
        mutate(run=i)
    )
  }
  
  # Null hypothesis is that the permuted distance is smaller than the observed for
  # each module (i.e., no signal). If we reject this hypothesis then there is
  # phylogenetic signal because the observed PD beteween hosts within each module
  # would be smaller than expected by chance (closely related hosts share a module).
  
  # Plot the means
  plt_across_modules <- 
    D_perm %>% group_by(run) %>%
    summarise(D_perm_mean=mean(d_mean)) %>%
    ggplot(aes(x=D_perm_mean))+geom_histogram()+geom_vline(xintercept = D_obs_mean)
  
  result_across_moduels <-
    D_perm %>% group_by(run) %>%
    summarise(D_perm=mean(d_mean)) %>%
    mutate(test=D_perm<D_obs_mean) %>%
    summarise(pvalue=sum(test)/nperm) %>%
    mutate(res=ifelse(pvalue<0.05,'Signal','No signal'))
  
  # This can also be tested per module
  plt_within_modules <-
    D_perm %>%
    full_join(D_obs, by='m') %>%
    rename(d_perm=d_mean.x, d_obs=d_mean.y) %>%
    ggplot(aes(x=d_perm))+
    geom_histogram()+
    facet_wrap(~m)+
    geom_vline(data = D_obs, aes(xintercept = d_mean))
  
  result_within_moduels <-
    D_perm %>%
    full_join(D_obs, by='m') %>%
    rename(d_perm=d_mean.x, d_obs=d_mean.y) %>%
    mutate(test=d_perm<d_obs) %>%
    group_by(m) %>%
    summarise(pvalue=sum(test)/nperm) %>%
    mutate(Signif=ifelse(pvalue<0.05,'Signal','No signal'),
           Signif_Bonferroni=ifelse(pvalue<0.05/nrow(D_obs),'Signal','No signal')) # Need to divide by number of modules for Bonferroni correction
  
  out <- list(D_obs=D_obs,
              D_obs_mean=D_obs_mean,
              plt_across_modules=plt_across_modules,
              plt_within_modules=plt_within_modules,
              result_across_moduels=result_across_moduels,
              result_within_moduels=result_within_moduels,
              nodes_in_modules=nodes_in_modules,
              nodes_in_tree=nodes_in_tree,
              overlapping=overlapping)
  return(out)
}

# Inititalize -------------------------------------------------------------

virus_data <- read_delim(paste(base_name,'_data-phage.txt',sep=''), delim=' ', col_names = T)
bacteria_data <- read_delim(paste(base_name,'_data-bact.txt',sep=''), delim=' ', col_names = T)
bacteria_abundance <- read_delim(paste(base_name,'_Bacteria-abundance.txt',sep=''), delim = ' ')
phage_abundance <- read_delim(paste(base_name,'_Phage-abundance.txt',sep=''), delim = ' ')

stop_time <- min(max(virus_data$timesOfRecord), max(phage_abundance$timesOfRecord))
hr_seq <- seq(1, stop_time, 1)

print(paste('-------- Working simulation:',base_name,' | stop time: ',stop_time,' | mu: ',mu,' | protospacers: ', protospacer_len,' | make plots: ',make_plots,'---------'))

virus_data %<>% filter(timesOfRecord<=stop_time)
bacteria_data %<>% filter(timesOfRecord<=stop_time)
bacteria_abundance %<>% filter(timesOfRecord<=stop_time)
phage_abundance %<>% filter(timesOfRecord<=stop_time)
 
regimes_df <- get_regimes(phage_time_series = phage_abundance, do_smoothing = T)$regimes_df
if (regimes_df[nrow(regimes_df),]$duration==1){regimes_df <- regimes_df[-nrow(regimes_df),]} # remove the end of simulation spurious effect
record_data(regimes_df)
BDRs <- subset(regimes_df, regime_type=='BDR')
VDRs <- subset(regimes_df, regime_type=='VDR')

# Vectors with the time points of VDRs and BDRs
VDR_hrs <- unlist(apply(VDRs, MARGIN = 1, FUN = function(x) seq(x[1],x[2])))
BDR_hrs <- unlist(apply(BDRs, MARGIN = 1, FUN = function(x) seq(x[1],x[2])))

# This is for the x axis labels when plotting
label_seq <- pretty(hr_seq, n=10)
label_seq <- subset(label_seq, label_seq<stop_time)
regimes_seq <- tibble(hr=c(VDR_hrs,BDR_hrs), regime_type=c(rep('VDR',length(VDR_hrs)),rep('BDR',length(BDR_hrs)))) %>% arrange(hr)
record_data(regimes_seq)

# bacteria / phage diversity ----------------------------------------------
notify('Generating abundance profile plots...')
dom_strains_num <- 100
cols <- c('gray50',gg_color_hue(dom_strains_num))

dominant_strains <- bacteria_abundance %>% 
  group_by(label) %>% 
  summarise(mean_abund=mean(Bdensity)) %>% 
  top_n(n=dom_strains_num, wt = mean_abund) %>% 
  arrange(label) %>% 
  mutate(dominant=T)
dom <- bacteria_abundance %>%
  filter(label%in%dominant_strains$label) %>% 
  select(timesOfRecord, Bdensity, label)
non_dom <- bacteria_abundance %>%
  filter(!label%in%dominant_strains$label) %>% 
  group_by(timesOfRecord) %>% summarise(Bdensity=sum(Bdensity)) %>% 
  mutate(label='0')
to_plot <- rbind(dom,non_dom)
to_plot$label <- factor(to_plot$label, levels=sort(as.numeric(unique(to_plot$label))))
plt_bact_abund <- 
  standard_plot(
    ggplot(to_plot, aes(x=timesOfRecord,y=Bdensity/10^5, fill=label))+
    geom_area(stat="identity", color='black', size=0.2)+
    scale_fill_manual(values = cols)+
    labs(y='Bact. abund. (*10^5)')+
    theme(legend.position = 'none')+
    geom_hline(yintercept = 10^5.5/10^5, linetype='dashed', size=0.4)
  )
  
# bacteria_abundance %>% ggplot(aes(x=Bdensity))+geom_density()
dominant_strains <- phage_abundance %>% 
  group_by(label) %>% 
  summarise(mean_abund=mean(Pdensity)) %>% 
  top_n(n=dom_strains_num, wt = mean_abund) %>% 
  arrange(label) %>% 
  mutate(dominant=T)
dom <- phage_abundance %>% 
  filter(label%in%dominant_strains$label) %>% 
  select(timesOfRecord, Pdensity, label)
non_dom <- phage_abundance %>%
  filter(!label%in%dominant_strains$label) %>% 
  group_by(timesOfRecord) %>% 
  summarise(Pdensity=sum(Pdensity)) %>% 
  mutate(label='0')
to_plot <- rbind(dom,non_dom)
to_plot$label <- factor(to_plot$label, levels=sort(as.numeric(unique(to_plot$label))))
plt_virus_abund <- 
  standard_plot(
    ggplot(to_plot, aes(x=timesOfRecord,y=Pdensity/10^7, fill=label))+
    geom_area(stat="identity", color='black',size=0.2)+
    scale_fill_manual(values = cols)+
    labs(y='Virus abund. (*10^7)')+
    theme(legend.position = 'none')
  )

plt_abundance_profiles <- plot_grid(plt_bact_abund, plt_virus_abund, nrow=2, align = 'vh')

make_png(plt_abundance_profiles)
make_svg(plt_abundance_profiles)

# Generate networks -------------------------------------------------------
print('Generating networks...')
all_networks <- vector(mode = 'list', length = length(hr_seq))
for (hr in hr_seq){
  nets <- create_networks_hr(virus_data, bacteria_data, hr)
  all_networks[[which(hr_seq==hr)]] <- nets
  notify(paste('Generated networks for time ',hr,'/',stop_time,sep=''))
}



# Measures of diversity ---------------------------------------------------

# Data frame for virus density (abundance)
virus_density <- NULL
for (hr in hr_seq){
  tmp <- all_networks[[which(hr_seq==hr)]]$virus_abund_hr
  tmp$hr <- hr
  virus_density <- rbind(virus_density, tmp)
}
total_density <- virus_density %>% group_by(hr) %>% summarise(total_density_hr=sum(density))
virus_density %<>% left_join(total_density) %>% select(hr,V_ID,density,total_density_hr)
record_data(virus_density)

# Data frame for bacteria density (abundance)
bacteria_density <- NULL
for (hr in hr_seq){
  tmp <- all_networks[[which(hr_seq==hr)]]$bacteria_abund_hr
  tmp$hr <- hr
  bacteria_density <- rbind(bacteria_density, tmp)
}
total_density <- bacteria_density %>% group_by(hr) %>% summarise(total_density_hr=sum(density))
bacteria_density %<>% left_join(total_density) %>% select(hr,B_ID,density,total_density_hr)
record_data(bacteria_density)

# Immunity network density and size
imm_density_size <- NULL
for (hr in hr_seq){
  x <- all_networks[[which(hr_seq==hr)]]$immunity_matrix
  d <- sum(x!=0)/(nrow(x)*ncol(x))
  imm_density_size <- rbind(imm_density_size, tibble(hr=hr, Density=d, Links=sum(x!=0), B=nrow(x), V=ncol(x), Size=nrow(x)+ncol(x)))
}
record_data(imm_density_size)
plt_immunity_network_size <- 
  standard_plot(
    imm_density_size %>% 
      select(hr,Density,Links,Size) %>% 
      gather(key = 'variable', value = 'value', -hr) %>% 
      ggplot(aes(hr, value))+
      geom_line()+
      facet_grid(variable~., scales = 'free_y')+
      labs(title='Immunity network information')
  )

make_png(plt_immunity_network_size)
make_svg(plt_immunity_network_size)

# Richness
richness <- NULL
for (hr in hr_seq){
  x <- all_networks[[which(hr_seq==hr)]]
  x <- x$nodes %>% group_by(type) %>% count()
  x$hr <- hr
  richness <- rbind(richness, x)
}
richness %<>% mutate(node_type=case_when(type==1~'viruses',
                                         type==2~'bacteria',
                                         type==3~'spacers'))
record_data(richness)

plt_richness <- 
  standard_plot(
    ggplot(richness, aes(hr, n, color=node_type))+
      geom_line(size=1.5)+
      scale_color_manual(values = c('blue','brown', 'red'))+
      theme(legend.position = 'none')+
      labs(y='Richness')
  )
make_png(plt_richness)

plt_bacteria_richness_abundace <- 
  plot_grid(plt_bact_abund,
            standard_plot(
              ggplot(richness %>% filter(node_type=='bacteria'), aes(hr, n))+
                geom_line(color='blue')+
                theme(legend.position = 'none')+
                labs(y='Bacteria richness')
            ), nrow=2, align = 'vh')
make_png(plt_bacteria_richness_abundace)
make_svg(plt_bacteria_richness_abundace)


# Phage and bacteria diversification ----------------------------------
print(' ------- Creating virus and bacteria persistence data frames -------')

# These data frames contain the extinction/mutation events and persistence of viruses/bacteria
virus_dynamics_list <- virus_data %>% 
  filter(timesOfRecord<=max(hr_seq)) %>% 
  select(timesOfRecord, label) %>% 
  rename(V_ID=label) %>% 
  mutate(V_ID=paste('V_',str_pad(V_ID, 4, 'left', '0'),sep=''), w=1) %>% 
  arrange(timesOfRecord, V_ID)
virus_dynamics_list %<>% group_by(V_ID) %>% summarise(birth=first(timesOfRecord), death=last(timesOfRecord)) %>% 
  mutate(persistence=death-birth+1)
virus_dynamics_list %<>% filter(birth>=min(hr_seq), death<=max(hr_seq))

record_data(virus_dynamics_list)

plt_virus_persistence <- 
  standard_plot(
      ggplot(virus_dynamics_list)+
      geom_rect(aes(ymin=parse_number(virus_dynamics_list$V_ID),
                    ymax=parse_number(virus_dynamics_list$V_ID), 
                    xmin=birth, 
                    xmax=death), color='red',
                size=0.7)+
      labs(y='Virus ID')+
      theme(legend.position = 'none')
  )
make_png(plt_virus_persistence)
make_svg(plt_virus_persistence)


bacteria_dynamics_list <-
  bacteria_data %>% 
  filter(timesOfRecord<=max(hr_seq)) %>% 
  select(timesOfRecord, label) %>% 
  rename(B_ID=label) %>% 
  mutate(B_ID=paste('B_',str_pad(B_ID, 4, 'left', '0'),sep=''), w=1) %>% 
  arrange(timesOfRecord, B_ID)
bacteria_dynamics_list %<>% group_by(B_ID) %>% summarise(birth=first(timesOfRecord), death=last(timesOfRecord)) %>% 
  mutate(persistence=death-birth+1)
bacteria_dynamics_list %<>% filter(birth>=min(hr_seq), death<=max(hr_seq))

record_data(bacteria_dynamics_list)

plt_bacteria_persistence <- 
  standard_plot(
    ggplot(bacteria_dynamics_list)+
      geom_rect(aes(ymin=parse_number(B_ID),
                    ymax=parse_number(B_ID), 
                    xmin=birth, 
                    xmax=death), color='blue',
                size=0.7)+
      labs(y='Bacteria ID')+
      theme(legend.position = 'none')
  )
make_png(plt_bacteria_persistence)
make_svg(plt_bacteria_persistence)

# Trees -------------------------------------------------------------------
tree_data <- read_delim(paste(base_name,'_Phage-TREE.txt',sep=''), delim = '\t',col_names = c("Recordtime","id","parent_id","creation_time"))
tree_data$id <- paste('V_',str_pad(tree_data$id, 4, 'left', '0'),sep='')
tree_data$parent_id <- paste('V_',str_pad(tree_data$parent_id, 4, 'left', '0'),sep='')
tree_data %<>% left_join(virus_dynamics_list %>% select(V_ID,death), by=c('id'='V_ID'))
tree_data[1,3] <- NA
tree_data$death[is.na(tree_data$death)]<-tree_data$creation_time[is.na(tree_data$death)]+1
tree <- nodes_dataframe_to_newick(tree_data)
writeLines(tree, 'tree.nwk')
tree_viruses <- treeio::read.tree('tree.nwk')
plt_viruses_tree <- 
  standard_plot(
    ggtree::ggtree(tree_viruses, ladderize = T) +
      ggtree::theme_tree2()
  )

make_png(plt_viruses_tree)
make_svg(plt_viruses_tree)


# Bacteria tree
if (file.exists(paste(base_name,'_Bacteria-TREE.txt',sep=''))){
  tree_data <- read_delim(paste(base_name,'_Bacteria-TREE.txt',sep=''), delim = '\t',col_names = c("Recordtime","id","parent_id","creation_time"))
  tree_data$id <- paste('B_',str_pad(tree_data$id, 4, 'left', '0'),sep='')
  tree_data$parent_id <- paste('B_',str_pad(tree_data$parent_id, 4, 'left', '0'),sep='')
  tree_data %<>% left_join(bacteria_dynamics_list %>% select(B_ID,death), by=c('id'='B_ID'))
  tree_data[1,3] <- NA
  tree_data$death[is.na(tree_data$death)]<-tree_data$creation_time[is.na(tree_data$death)]+1
  tree <- nodes_dataframe_to_newick(tree_data)
  writeLines(tree, 'tree_bacteria.nwk')
  tree_bacteria <- treeio::read.tree('tree_bacteria.nwk')
  plt_bacteria_tree <- 
    standard_plot(
      ggtree::ggtree(tree_bacteria, ladderize = F) +
        ggtree::theme_tree2()
    )
    
  make_png(plt_bacteria_tree)
  make_svg(plt_bacteria_tree)
}

# Modularity of infection networks ----------------------------------------
modules_df_infection <- NULL
for (hr in hr_seq){
  notify(paste('Infection networks (modularity) ',hr,'/',stop_time,sep=''))
  edges <- all_networks[[which(hr_seq==hr)]]$infection_list
  if(is.null(edges)){next}
  nodes <- all_networks[[which(hr_seq==hr)]]$nodes
  x <- create_monolayer_object(edges, bipartite = T)
  modules <- run_infomap_monolayer(x, infomap_executable = 'Infomap', flow_model = 'undirected', silent = T, 
                                   trials = 100, two_level = T, seed = 123, 
                                   signif = F)
  if (!is.null(modules)){
    x <- modules$modules %>% select(node_id, node_name, m=module_level1) %>% mutate(hr=hr)
    modules_df_infection <- rbind(modules_df_infection, x)
    if (make_plots){
      png(paste('plots/infection_modules_',str_pad(hr,4,'left','0'),'.png',sep=''),1920,1080,res=150)
      print(
        plot_modular_matrix(modules)+
          theme(axis.text = element_text(size=6),
                axis.text.x = element_text(angle=-90))+
          labs(title=hr)
      )
      dev.off()
    }
  }
}
record_data(modules_df_infection)

plt_modules_infection <- 
  standard_plot(
    modules_df %>% group_by(hr) %>% summarise(num_mod=length(unique(m))) %>% 
      ggplot(aes(hr,num_mod))+
      geom_line(size=1.2)+
      labs(y='Number of modules')+
      theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.position = 'none'
      )
  )

# Significance of modularity of infection networks --------------------------------------
infection_modularity_signif <- NULL
# Test at the end of each VDR
for (i in 1:nrow(VDRs)){
  hr <- VDRs$end[i]
  edges <- all_networks[[which(hr_seq==hr)]]$infection_list
  x <- create_monolayer_object(edges, directed = F, bipartite = T)
  test <- run_infomap_monolayer(x, infomap_executable = 'Infomap', flow_model = 'undirected', silent = T, trials = 100, two_level = T, seed = 123, signif = T, shuff_method = 'r00', nsim = 1000)
  infection_modularity_signif %<>% bind_rows(tibble(hr=hr, 
                                                    pvalue=test$pvalue,
                                                    n_hosts=ncol(x$mat),
                                                    n_spacer=nrow(x$mat),
                                                    n_interactions=sum(x$mat>0),
                                                    n_modules=max(test$modules$module_level1)))
  
  plot_modular_matrix(test)
}
record_data(infection_modularity_signif)



# Phylogenetic signal in infection networks --------------------------------------
phylogenetic_signal_infection <- NULL
# Test at the end of each VDR
for (i in 1:nrow(VDRs)){
  hr <- VDRs$end[i]
  edges <- all_networks[[which(hr_seq==hr)]]$infection_list
  x <- create_monolayer_object(edges, directed = F, bipartite = T)
  infection_modularity <- run_infomap_monolayer(x, infomap_executable = 'Infomap', flow_model = 'undirected', silent = T, trials = 100, two_level = T, seed = 123, signif = F)
  pd_results <- test_PD_modules(tree = tree_bacteria, module_object = infection_modularity, node_start_letter = 'B')
  out <- pd_results$result_within_moduels %>% left_join(pd_results$D_obs)
  out$VDR <- i
  phylogenetic_signal_infection <- rbind(phylogenetic_signal_infection, out)
}
record_data(phylogenetic_signal_infection)



# Significance of modularity of host-spacer networks --------------------------------------

# Agregate host-spacer networks within each VDR and test for significance of modularity
host_sp_modularity <- NULL
for (i in 1:nrow(VDRs)){
  metaweb <- NULL
  for (hr in VDRs$start[i]:VDRs$end[i]){
    edges <- all_networks[[which(hr_seq==hr)]]$bacteria_sp_list
    if(is.null(edges)){next}
    metaweb %<>% bind_rows(edges) %>% distinct(B_ID, SP)
  }
  metaweb %<>% mutate(w=1) 
  x <- create_monolayer_object(metaweb, directed = F, bipartite = T)
  test <- run_infomap_monolayer(x, infomap_executable = 'Infomap', flow_model = 'undirected', silent = T, trials = 100, two_level = T, seed = 123, signif = T, shuff_method = 'r00', nsim = 100)
  host_sp_modularity %<>% bind_rows(tibble(hr=hr, 
                                             pvalue=test$pvalue,
                                             n_hosts=ncol(x$mat),
                                             n_spacer=nrow(x$mat),
                                             n_interactions=sum(x$mat>0),
                                             n_modules=max(test$modules$module_level1)))
  
}
record_data(host_sp_modularity)


# Phylogenetic signal in host-spacer modules ------------------------------
# Agregate host-spacer networks within each VDR and test for significance of modularity
phylogenetic_signal_hs <- NULL
for (i in 1:nrow(VDRs)){
  metaweb <- NULL
  for (hr in VDRs$start[i]:VDRs$end[i]){
    edges <- all_networks[[which(hr_seq==hr)]]$bacteria_sp_list
    if(is.null(edges)){next}
    metaweb %<>% bind_rows(edges) %>% distinct(B_ID, SP)
  }
  metaweb %<>% mutate(w=1) 
  x <- create_monolayer_object(metaweb, directed = F, bipartite = T)
  host_sp_modularity <- run_infomap_monolayer(x, infomap_executable = 'Infomap', flow_model = 'undirected', silent = T, 
                                              trials = 100, two_level = T, seed = 123, 
                                              signif = F)
  pd_results <- test_PD_modules(tree = tree_bacteria, module_object = host_sp_modularity, node_start_letter = 'B')
  out <- pd_results$result_within_moduels %>% left_join(pd_results$D_obs)
  out$VDR <- i
  phylogenetic_signal_hs <- rbind(phylogenetic_signal_hs, out)
}
record_data(phylogenetic_signal_hs)





# WNODF -------------------------------------------------------------------

WNODF_df <- NULL
for (hr in hr_seq){
  x <- all_networks[[which(hr_seq==hr)]]$immunity_matrix
  notify(paste('Immunity networks (WNODF nestedness) ',hr,'/',stop_time,' | dimensions: ',nrow(x),' by ',ncol(x),sep=''))
 
  if (nrow(x)<2 || ncol(x)<2 || all(x==0)){
    WNODF <- NA
  } else {
    WNODF <- bipartite::networklevel(x, index = 'weighted NODF')
  }
  WNODF_df <- rbind(WNODF_df, tibble(hr=hr, WNODF))
}
record_data(WNODF_df)
plt_WNODF <- 
  standard_plot(
    ggplot(WNODF_df, aes(hr, WNODF))+
      geom_line(color='#1E5984', size=1.1)+
      labs(y='WNODF')
  )
make_png(plt_WNODF)
make_svg(plt_WNODF)

# Spacer matches ------------------------------------
spacer_matches <- NULL
match_cols <- tibble(value=0:6, col=c('gray','#f9ed69','#f08a5d','#b83b5e','#6a2c70','#61210F','#3C3438'))
# ggplot(match_cols, aes(value,value,color=col))+geom_point(size=10)+scale_color_identity()
for (hr in hr_seq){
  x <- all_networks[[which(hr_seq==hr)]]$immunity_matrix
  notify(paste('Spacer matches ',hr,'/',stop_time,' | dimensions: ',nrow(x),' by ',ncol(x),sep=''))
  y <- tibble(hr=hr, num_matches=as.integer(names(table(x))), cnt=as.vector(table(x)))
  spacer_matches <- rbind(spacer_matches, y)

  # Plot with fixed colors for matches
  if(on_Midway() && make_plots){
    M <- t(x)
    M_orig <- M
    if (nrow(M)==1 & ncol(M)==1){
      M <- tibble(Var1=colnames(M), Var2=rownames(M), value=M[1,1])
    } else {
      M <- M[order(rowSums(M), decreasing = T), order(colSums(M),decreasing = F)]
      M <- reshape2::melt(M)
      if (ncol(M)==1){ # If M has only one dimension after ordering M for the nestedness then reshaping does not work well.
        M$Var2 <- rownames(M)
        M$Var1 <- dimnames(M_orig)[[which(dim(M_orig)==1)]]
      }
    }
    suppressMessages(suppressWarnings(M <- as_tibble(M) %>% left_join(match_cols)))
    M %<>% mutate(col=ifelse(value>6,'black',col))
    png(paste('plots/immunity_nested_',str_pad(hr,4,'left','0'),'.png',sep=''),1920,1080,res=150)
    print(
      M %>% 
        ggplot()+
        geom_tile(aes(Var1,Var2,fill=col))+
        scale_fill_identity('# matches',drop=FALSE, guide = 'legend', labels=match_cols$value, breaks=match_cols$col)+
        coord_fixed()+
        theme(
          # legend.position = 'none',
          axis.text.x=element_text(angle=-90),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size=6))+
        labs(title=hr) 
    )
    dev.off()
  }
}

total_spacer_matches <- spacer_matches %>% group_by(hr) %>% summarise(tot=sum(cnt))
spacer_matches %<>% 
  left_join(total_spacer_matches) %>% 
  mutate(prop=cnt/tot)
record_data(spacer_matches)

# Calculate the proportion of spacer matches
plt_spacer_matches <-
  standard_plot(
    spacer_matches %>% 
      left_join(match_cols, by=c('num_matches'='value')) %>% 
      mutate(col=ifelse(num_matches>6,'black',col)) %>% # more than 6 matches are colored black
      ggplot(aes(x=hr, prop, color=col))+geom_line(size=1.5)+
      scale_color_identity('# matches',drop=FALSE, guide = 'legend', labels=match_cols$value, breaks=match_cols$col)+
      labs(y='% matches')+
      theme(
        legend.position = c(0.9,0.9),
            legend.direction='horizontal',
            legend.justification = c("right", "top"),
            legend.background = element_rect(fill='white'),
            legend.title=element_blank()
      )
  )
make_png(plt_spacer_matches)
make_svg(plt_spacer_matches)


# Extinctions -------------------------------------------------------------

# Get the order in which extinctions happen: Is it from the most protected to the least protected?
notify('Analyzing extinction ranks')
ext_immunity_rank <- NULL
for (i in 1:nrow(virus_dynamics_list)){
  hr <- virus_dynamics_list[i,3]$death
  v_ext <- virus_dynamics_list[i,1]$V_ID
  # Pull out the immunity network at the time. 
  x <- all_networks[[which(hr_seq==hr)]]$immunity_matrix
  # Get the immunity ranks. Duplicates eliminated to avoid ties
  virus_immunity_total <- sort(unique(colSums(x)), decreasing = F)
  virus_immunity_rank <- which(virus_immunity_total%in%colSums(x)[v_ext])
  # Calculate the relative rank
  ext_immunity_rank <- rbind(ext_immunity_rank, tibble(V_ID=v_ext, death=hr, ext_immunity_rank=virus_immunity_rank/length(virus_immunity_total)))
  # ext_immunity_rank <- c(ext_immunity_rank, virus_immunity_rank/length(virus_immunity_total))
}
record_data(ext_immunity_rank)

plt_extinction_rank <- ggplot(ext_immunity_rank, aes(ext_immunity_rank))+
  geom_histogram(fill=NA, color='tomato')+
  geom_vline(xintercept = median(virus_dynamics_list$ext_immunity_rank, na.rm = T), linetype='dashed')+
  labs(x='Immunity rank of extinct viruses')

make_png(plt_extinction_rank)
make_svg(plt_extinction_rank)


virus_dynamics_list %<>%
  left_join(ext_immunity_rank) %>%
  mutate(high_rank=ifelse(ext_immunity_rank>median(ext_immunity_rank),T,F))



# R0 for 0 and 1 matches --------------------------------------------------

# R0 (whether for 0 or 1 matches) is calculated per bacteria-virus pair because
# it is the growth rate of a virus strain in a bacteria strain. To obtain the R0
# of a virus, R0_j, R0_ij needs to be summed across bacteria. R_0m is calcualted
# for interactions with 0-matches (where viruses can grow). R_1m is calcualted
# for interactions with 1-matches, and represents the potential addition to R0
# if a virus escaped a bacteria.

notify('--> Calculate R0 for 0 and 1 matches...')
R_0m_df <- NULL
for (hr in hr_seq){
  x <- all_networks[[which(hr_seq==hr)]]
  imm <- x$immunity_matrix
  if (!0%in%imm){print('No 0-matches in immunity matrix, skipping');next}
  
  # Prepare a matrix for R_0m
  R_0m <- imm
  R_0m[,] <- 0
  bact <- rownames(imm)
  bact_abund <- x$bacteria_abund_hr %>% filter(B_ID%in%bact)
  N_T <- sum(x$bacteria_abund_hr$density) # Total bacteria abundance
  Z <- (beta*phi*(1-q))/(phi*N_T+m)
  
  # Calculate the R_0m for interactions with 0-matches (where viruses can grow)
  idx <- which(imm==0, arr.ind = T) # Get only the 0-matches
  for (r in 1:nrow(idx)){
    i <- idx[r,1]
    j <- idx[r,2]
    bact <- rownames(imm)[i]
    virus <- colnames(imm)[j]
    N_i <- subset(bact_abund, B_ID==bact)$density
    R_0m[i,j] <- Z*N_i
  }
  
  # Sum across bacteria to get the R_0m of a virus
  tmp <- data.frame(hr=hr, V_ID=names(colSums(R_0m)), R_0m=colSums(R_0m), stringsAsFactors = F)
  
  # Weigh R_0m by relative virus abundance
  virus_abund <- x$virus_abund_hr %>% filter(V_ID%in%tmp$V_ID)
  virus_abund$rel_density <- virus_abund$density/sum(virus_abund$density)
  suppressMessages(tmp %<>% 
                     left_join(virus_abund) %>%
                     mutate(R_0m_w_j=R_0m*rel_density))
  
  R_0m_df <- rbind(R_0m_df, as_tibble(tmp))

  notify(paste('Calculated R_0m for time',hr))
}

record_data(R_0m_df)

# plot population measures of R_0m
plt_R0m <- 
standard_plot(
  R_0m_df %>% group_by(hr) %>%
    summarise(R_0m_mean=mean(R_0m),
              R_0m_mean_w=sum(R_0m_w_j),# Need to SUM the R_0m_w across viruses because this gives the weighted average for the population
              R_0m_max=max(R_0m)) %>% 
    gather(key='variable', value='value', -hr) %>% 
    ggplot(aes(hr, value))+
    geom_line(color='#8ACAFE')+
    facet_grid(variable~., scales = 'free_y')+
    geom_hline(yintercept = 1, linetype='dashed', color='gray50')+
    labs(title = 'Virus population measures for R with 0 matches')
)
make_png(plt_R0m)
make_svg(plt_R0m)

R_1m_df <- NULL
for (hr in hr_seq){
  x <- all_networks[[which(hr_seq==hr)]]
  imm <- x$immunity_matrix
  if (!1%in%imm){print('No 1-matches in immunity matrix, skipping');next}
  
  # Prepare a matrix for R_1m
  R_1m <- imm
  R_1m[,] <- 0
  bact <- rownames(imm)
  bact_abund <- x$bacteria_abund_hr %>% filter(B_ID%in%bact)
  
  N_T <- sum(x$bacteria_abund_hr$density) # Total bacteria abundance
  Z <- (beta*phi*(1-q))/(phi*N_T+m)
  
  # Get only the 1-matches
  idx <- which(imm==1,arr.ind = T) 
  for (r in 1:nrow(idx)){
    i <- idx[r,1]
    j <- idx[r,2]
    bact <- rownames(imm)[i]
    virus <- colnames(imm)[j]
    N_i <- subset(bact_abund, B_ID==bact)$density
    R_1m[i,j] <- Z*N_i/protospacer_len # Need to divide by the length of the protospacer cassette because an escape mutation must hit the right protospacer
  }

  # Sum across bacteria to get the R_1m of a virus
  tmp <- data.frame(hr=hr, V_ID=names(colSums(R_1m)), R_1m=colSums(R_1m), stringsAsFactors = F)
  
  # Weigh R_1m by relative virus abundance
  virus_abund <- x$virus_abund_hr %>% filter(V_ID%in%tmp$V_ID)
  virus_abund$rel_density <- virus_abund$density/sum(virus_abund$density)
  suppressMessages(tmp %<>%
                     left_join(virus_abund) %>%
                     mutate(R_1m_w_j=R_1m*rel_density))
  
  R_1m_df <- rbind(R_1m_df, as_tibble(tmp))

  notify(paste('Calculated R_1m for time',hr))
}

record_data(R_1m_df)

# plot population measures of R_1m
plt_R1m <- 
  standard_plot(
    R_1m_df %>% group_by(hr) %>%
      summarise(R_1m_mean=mean(R_1m),
                R_1m_mean_w=sum(R_1m_w_j),# Need to SUM the R_0m_w across viruses because this gives the weighted average for the population
                R_1m_max=max(R_1m)) %>% 
      gather(key='variable', value='value', -hr) %>% 
      ggplot(aes(hr, value))+
      geom_line(color='#F066D8')+
      facet_grid(variable~., scales = 'free_y')+
      geom_hline(yintercept = 1, linetype='dashed', color='gray50')+
      labs(title = 'Virus population measures for R with 1 matches')
  )

make_png(plt_R1m)
make_svg(plt_R1m)


# Calculate the R_potential
R_pot_df <-
  # Need to join 0-matches with 1-matches to sum them
  full_join(R_0m_df,R_1m_df,by=c('hr','V_ID')) %>% 
  # Convert NAs resulting from the join to 0s
  mutate(R_0m=ifelse(is.na(R_0m),0,R_0m),
         R_1m=ifelse(is.na(R_1m),0,R_1m),
         R_0m_w_j=ifelse(is.na(R_0m_w_j),0,R_0m_w_j),
         R_1m_w_j=ifelse(is.na(R_1m_w_j),0,R_1m_w_j)) %>% 
  # Remove variables that are now replicated
  mutate(density=ifelse(is.na(density.x), density.y, density.x)) %>% 
  mutate(rel_density=ifelse(is.na(rel_density.x), rel_density.y, rel_density.x)) %>% 
  # Select and organize the columns
  select(-density.x, -density.y, -rel_density.x, -rel_density.y) %>% 
  select("hr","V_ID","density","rel_density","R_0m","R_0m_w_j","R_1m","R_1m_w_j") %>% 
  # Caclulate R_pot and weighted R_pot
  mutate(Rpot=R_0m+R_1m,
         Rpot_w_j=R_0m_w_j+R_1m_w_j)

record_data(R_pot_df)

plt_R_pot <- 
standard_plot(
  R_pot_df %>% group_by(hr) %>%
    summarise(R_pot_mean=mean(Rpot),
              R_pot_mean_w=sum(Rpot_w_j),# Need to SUM the R_0m_w across viruses because this gives the weighted average for the population
              R_pot_max=max(Rpot)) %>% 
    gather(key='variable', value='value', -hr) %>% 
    ggplot(aes(hr, value))+
    geom_line()+
    facet_grid(variable~., scales = 'free_y')+
    geom_hline(yintercept = 1, linetype='dashed', color='gray50')+
    labs(title = 'Virus population measures for R_pot')
)
make_png(plt_R_pot)
make_svg(plt_R_pot)


# Plot areas where an escape mutation can lead to Rpot>1
R_pot_summary <- 
  R_pot_df %>% group_by(hr) %>% 
  summarise(R_0m_mean_w=sum(R_0m_w_j),
            R_1m_mean_w=sum(R_1m_w_j),
            R_pot_mean_w=sum(Rpot_w_j)) 
plt_R_pot_effect <- 
standard_plot(
  ggplot()+
    geom_hline(yintercept = 1, linetype='dashed', color='gray50')+
    geom_line(data=R_pot_summary, aes(hr,R_0m_mean_w),color='#8ACAFE')+
    geom_line(data=R_pot_summary, aes(hr,R_1m_mean_w),color='#F066D8')+
    geom_point(data=R_pot_summary %>% group_by(hr) %>% filter(R_0m_mean_w<=1, R_pot_mean_w>1, hr%in%BDR_hrs),
               aes(hr,R_pot_mean_w),color='red')+
    labs(y='Value of R0 (weighted means)', title='Times when R_pot>1 in BDRs, potentially causing an escape')
)
make_png(plt_R_pot_effect)
make_svg(plt_R_pot_effect)


# Plot -----------------------------------------------------------


notify('--> Generate final plots...')

 # A PDF with all the main plots
pdf(paste(base_name,'_main_figures.pdf',sep=''), 16, 10, onefile = T)
  title <- ggdraw() + draw_label("Viral and bacterial strain abundance", fontface='bold')
  grid.arrange(plot_grid(title, plt_abundance_profiles, ncol=1, rel_heights=c(0.1, 1)))
  grid.arrange(plt_richness+labs(title = 'Richness of viruses, hosts and spacers'))
  title <- ggdraw() + draw_label("Bacteria abundance and richness", fontface='bold')
  grid.arrange(plot_grid(title, plt_bacteria_richness_abundace, ncol=1, rel_heights=c(0.1, 1)))
  grid.arrange(plt_spacer_matches+labs(title = 'Proportion of spacer matches'))
  grid.arrange(plt_WNODF+labs(title = 'Nestedness of immunity network'))
  grid.arrange(plt_extinction_rank+labs(title = 'Distribution of extinction ranks'))
  grid.arrange(plt_R0m)
  grid.arrange(plt_R1m)
  grid.arrange(plt_R_pot)
  grid.arrange(plt_R_pot_effect)
  grid.arrange(plt_virus_persistence+labs(title = 'Diversification and persistence of viruses'))
  grid.arrange(plt_bacteria_persistence+labs(title = 'Diversification and persistence of bacteria'))
  grid.arrange(plt_viruses_tree+labs(title = 'Viruses phylogenetic tree'))
  if ('plt_bacteria_tree' %in% ls(pattern = 'plt_bacteria_tree')){
    grid.arrange(plt_bacteria_tree+labs(title = 'Bacteria phylogenetic tree'))
  }
dev.off()

writeLines('success!', paste(base_name,'_success.txt',sep=''))
notify('--------- DONE! ---------')