# !diagnostics off

# parameters of the model; used by the functions, so need to initialize here
beta=50;phi=10^-7;p=10^-5;q=10^-5;m=0.1;spacer_len=10

on_Midway <- function(){ifelse(Sys.getenv('USER')=='pilosofs',T,F)}

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
complete_analysis <- as.logical(args[5])

if(on_Midway()){system('module load gcc/6.1')}
if(!on_Midway()){setwd(paste('/Users/Shai/Dropbox (BGU)/Projects/CRISPR/simulation_data/',base_name,sep=''))}

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

empty <- function (web, count = FALSE) 
{
  web[is.na(web)] <- 0
  if (NCOL(web) == 1 | NROW(web) == 1) {
    if (NCOL(web) == 1 & NROW(web) != 1) {
      nr <- sum(web > 0)
      nc <- 1
    }
    if (NROW(web) == 1 & NCOL(web) != 1) {
      nc <- sum(web > 0)
      nr <- 1
    }
    if (NROW(web) == 1 & NCOL(web) == 1) {
      nr <- 1
      nc <- 1
    }
    out <- web[1:nr, 1:nc, drop = FALSE]
    if (count) 
      attr(out, "empty") <- c(`empty rows` = NROW(web) - 
                                nr, `empty columns` = NCOL(web) - nc)
    return(out)
  }
  cempty <- which(colSums(web) == 0)
  rempty <- which(rowSums(web) == 0)
  cind <- if (length(cempty) == 0) 
    1:NCOL(web)
  else (1:NCOL(web))[-cempty]
  rind <- if (length(rempty) == 0) 
    1:NROW(web)
  else (1:NROW(web))[-rempty]
  out <- web[rind, cind, drop = FALSE]
  if (count) 
    attr(out, "empty") <- c(`empty rows` = length(rempty), 
                            `empty columns` = length(cempty))
  return(out)
}

# Matrix trace
tr <- function(x){
  sum(diag(x))
}

calculate_ev_nestedness <- function(B){
  # It is faster to calculate the ev for smaller matrices. Because the leading
  # ev of BB^T and B^TB is the same, we first check how to produce A.
  if (nrow(B)<ncol(B)){
    A <- B%*%t(B)
  } else {
    A <- t(B)%*%B
  }
  ev_max <- max(eigen(A, symmetric = T, only.values = T)$values) # Not calculating eigenvectors speeds calculatoins remarkably.
  return(ev_max)
}

# The upper bound is taken from Stanić Z. Inequalities for Graph Eigenvalues.
# Cambridge University Press; 2015, page 51. It is calculated as:
# \lambda_1^2\leq \frac{tr(BB^T)}{n_1}+\sqrt{\frac{n_1-1}{n_1}\left( tr((BB^T)^2)-\frac{tr(BB^T)^2}{n_1} \right)}
# where n1 is the shorter dimension of B (also see function calculate_ev_nestedness)
calculate_upper_bound <- function(B){
  if (nrow(B)<ncol(B)){
    A <- B%*%t(B)
  } else {
    A <- t(B)%*%B
  }
  n1 <- nrow(A)
  x <- tr(A)/n1
  y <- (n1-1)/n1
  z <- tr(A%*%A)-tr(A)^2/n1
  upper_bound <- x+sqrt(y*z)
  return(upper_bound)
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

plot_matrix <- function(M, layout='random', method = 'ggplot', binary_cols=c('gray','red'), title='', x_title='', y_title='', legend_title=''){
  M_orig <- M
  # M=M_orig
  if (layout=='random'){
    M <- M[sample(1:nrow(M), size = nrow(M), replace = F), sample(1:ncol(M), size = ncol(M), replace = F)]
  }
  
  if (layout == "diagonal") {
    ca <- cca(M)
    M <- M[order(summary(ca)$sites[, 1], decreasing = TRUE), 
           order(summary(ca)$species[, 1], decreasing = TRUE)]
  }
  
  if (layout=='nested' & method == 'ggplot'){
    tmp <- as.matrix(M[order(rowSums(M), decreasing = T), order(colSums(M),decreasing = F)])
    rnames <- rownames(M)[order(rowSums(M), decreasing = T)]
    cnames <- colnames(M)[order(colSums(M), decreasing = F)]
    rownames(tmp) <- rnames
    colnames(tmp) <- cnames
    M <- tmp
  }
  if (layout=='nested' & method == 'heatmap'){
    M <- M[order(rowSums(M), decreasing = T), order(colSums(M),decreasing = T)]
  }
  
  if (identical(unique(as.vector(M)), c(0,1))){ # If binary matrix
    colors=binary_cols
  } else {
    if (all(M==0)){ # If empty matrix
      colors='gray'
    } else {
        colors=c('gray',gg_color_hue(n=max(M))) # If neither binary nor empty
      }
  }
    
  
  if (method=='heatmap'){
    heatmap(M, Rowv = NA, Colv = NA, 
            symm = F, scale = 'none', revC = T, margins=c(7,7), 
            # labRow = F, labCol = F, 
            col=colors, main=title)
    return(invisible())
  }
  
  if (method=='ggplot'){
    M <- reshape2::melt(M)
    if (ncol(M)==1){ # If M has only one dimension after ordering M for the nestedness then reshaping does not work well.
      M$Var1 <- rownames(M)
      M$Var2 <- dimnames(M_orig)[[which(dim(M_orig)==1)]]
    }
    plt=as_tibble(M) %>% 
      ggplot()+
      geom_tile(aes(Var1,Var2,fill=value))+
      scale_fill_gradientn(colours=colors, name=legend_title)+
      theme(
        axis.text.x=element_text(angle=-90),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 18))+
      labs(title=title, x=x_title,y=y_title)+
      coord_fixed()
    return(list(M=M,plt=plt))
  }
}

gg_color_hue <- function(n, hue_min = 10, hue_max = 280, tune1 = 62, tune2 = 100) {
  hues = seq(hue_min, hue_max, length=n+1)
  hcl(h=hues, l=tune1, c=tune2)[1:n]
}

matrix_to_list_bipartite <- function(x){
  g <- graph.incidence(t(x), weighted = T)
  l_bip <- as_tibble(igraph::as_data_frame(g, "edges"))
  return(l_bip)
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

# A function to run Infomap on a unipartite/bipartite network.
run_infomap <- function(edges, remove_zero_edges = T, return_L=F){
  if (nrow(edges)==0 || is.null(edges)){
    print('No edges!! Returning NULL')
    return (NULL)
  } else {
    if (remove_zero_edges){
      edges %<>% filter(w!=0)
    }
    edges <- as.data.frame(edges[,1:3])
    # Create internal index for nodes
    Set1 <- sort(unique(edges[,1]))
    Set2 <- sort(unique(edges[,2]))
    nodes <- data.frame(nodeID=1:(length(Set1)+length(Set2)), nodeName=c(Set1,Set2))
    
    # Transform node names to IDs
    edges[,1] <- nodes$nodeID[match(edges[,1],nodes$nodeName)]
    edges[,2] <- nodes$nodeID[match(edges[,2],nodes$nodeName)]
    
    write_delim(edges, 'infomap.txt', delim = ' ', col_names = F)
    # system('./Infomap_v01926 infomap.txt . -N 20 --tree -2 --undirected --seed 123 -i link-list --silent --out-name infomap_out')
    system('./Infomap_mac_01925 infomap.txt . -N 20 --tree -2 --undirected --seed 123 -i link-list --silent --out-name infomap_out')
    modules_hr <- suppressMessages(read_delim('infomap_out.tree', delim = ' ', skip = 2, col_names = c('path', 'flow', 'name', 'node')))
    modules_hr %<>% dplyr::select(name, path, nodeID=node) %>%
      rowwise %>%
      mutate(module=str_split(path,':')[[1]][1]) %>% 
      dplyr::select(nodeID, module)
    suppressMessages(
      nodes %<>% inner_join(modules_hr) %>% mutate(module=as.integer(module), hr) %>% 
        mutate(mod_hr=paste(hr,module,sep='_')) %>%  # This is a unique identifier of the module because the number of module will repeat in other hours.
        dplyr::select(-nodeID)   # Remove that because it was an internal index for that function, so to not confuse with values outside the function
    )
    
    if (return_L){
      L <- read_lines('infomap_out.tree', n_max = 1)
      L <- parse_number(str_sub(L, str_locate(L, 'to codelength')[1], str_length(L)))
      return(list(L=L, modules=nodes))
    } else {
      return(as.tibble(nodes))
    }
  }
}



plot_modules <- function(network_list, m_df, remove_zero_edges=T){
  if (remove_zero_edges){
    network_list %<>% filter(w!=0)
  }
  M_set1 <- M_set2 <- network_list <- network_list[,1:3]
  names(M_set1) <- names(M_set2) <- names(network_list) <- c('Set1','Set2','w')
  suppressMessages(suppressWarnings(
    M_set1 %<>% left_join(m_df, by=c('Set1'='nodeName')) %>% rename(module1=module,mod_hr1=mod_hr)
  ))
  suppressMessages(suppressWarnings(
    M_set2 %<>% left_join(m_df, by=c('Set2'='nodeName')) %>% rename(module2=module,mod_hr2=mod_hr)
  ))
  suppressMessages(suppressWarnings(
    M <- full_join(M_set1, M_set2, by = c("Set1", "Set2", "w")) %>% 
      dplyr::select(hr=hr.x,Set1,Set2,w,module1,module2,mod_hr1,mod_hr2)
  ))
  Set1_modules <- unique(M_set1[,c('Set1','module1')])
  Set1_modules <- with(Set1_modules, Set1_modules[order(module1,Set1),])
  Set2_modules <- unique(M_set2[,c('Set2','module2')])
  Set2_modules <- with(Set2_modules, Set2_modules[order(module2,Set2),])
  
  # If there are no interactions outside the module then do not need the gray
  # color. Otherwise, it will plot the first module in gray.
  M %<>% mutate(edge_in_out=ifelse(module1==module2,'in','out')) %>% 
    mutate(value_mod=ifelse(edge_in_out=='in',module1,0)) %>% 
    mutate(Set1=factor(Set1, levels=Set1_modules$Set1), Set2=factor(Set2, levels=Set2_modules$Set2))
  
  module_colors <- tibble(module1=unique(M$module1), col=gg_color_hue(n=length(unique(M$module1))))
  
  suppressMessages(suppressWarnings(
  M %<>% left_join(module_colors) %>% mutate(col=ifelse(edge_in_out=='in',col,'gray'))
  ))
  
  plt <- 
    ggplot()+
    geom_tile(data=M %>% filter(w!=0), aes(Set1, Set2,fill=col), colour='black') + # Interactions within modules
     geom_tile(data=M %>% filter(w==0), aes(Set1, Set2), fill='white') + # Add nodes with no interactions, if they exist
    scale_fill_identity()+
    theme(legend.position='none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle=-90),
          axis.title = element_blank(),
          axis.ticks = element_blank())+
    coord_fixed()

  # geom_text(aes(V_ID,B_ID,label = value_mod))
  return(plt)
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
dir.create('figures')

library(ggtree)
library(ape)
library(treeio)
library(igraph)
library(tidyverse)
library(magrittr)
library(bipartite)
library(cowplot)
library(gridExtra)
library(grid)
library(egg)
library(infomapecology)

virus_data <- read_delim(paste(base_name,'_data-phage.txt',sep=''), delim=' ', col_names = T)
bacteria_data <- read_delim(paste(base_name,'_data-bact.txt',sep=''), delim=' ', col_names = T)
bacteria_abundance <- read_delim(paste(base_name,'_Bacteria-abundance.txt',sep=''), delim = ' ')
phage_abundance <- read_delim(paste(base_name,'_Phage-abundance.txt',sep=''), delim = ' ')

stop_time <- min(max(virus_data$timesOfRecord), max(phage_abundance$timesOfRecord))
hr_seq <- seq(1, stop_time, 1)

print(paste('-------- Working simulation:',base_name,' | stop time: ',stop_time,' | mu: ',mu,' | protospacers: ', protospacer_len, ' | complete analysis: ',complete_analysis,' | make plots: ',make_plots,'---------'))

if (complete_analysis){
  print('Complete analysis')
} else {
  print('REDUCED analysis')
}


virus_data %<>% filter(timesOfRecord<=stop_time)
bacteria_data %<>% filter(timesOfRecord<=stop_time)
bacteria_abundance %<>% filter(timesOfRecord<=stop_time)
phage_abundance %<>% filter(timesOfRecord<=stop_time)
 
regimes_df <- get_regimes(phage_time_series = phage_abundance, do_smoothing = T)$regimes_df
if (regimes_df[nrow(regimes_df),]$duration==1){regimes_df <- regimes_df[-nrow(regimes_df),]} # remove the end of simulation spurious effect
if(on_Midway()){record_data(regimes_df)}
BDRs <- subset(regimes_df, regime_type=='BDR')
VDRs <- subset(regimes_df, regime_type=='VDR')

# Vectors with the time points of VDRs and BDRs
VDR_hrs <- unlist(apply(VDRs, MARGIN = 1, FUN = function(x) seq(x[1],x[2])))
BDR_hrs <- unlist(apply(BDRs, MARGIN = 1, FUN = function(x) seq(x[1],x[2])))

# This is for the x axis labels when plotting
label_seq <- pretty(hr_seq, n=10)
label_seq <- subset(label_seq, label_seq<stop_time)
regimes_seq <- tibble(hr=c(VDR_hrs,BDR_hrs), regime_type=c(rep('VDR',length(VDR_hrs)),rep('BDR',length(BDR_hrs)))) %>% arrange(hr)
if(on_Midway()){record_data(regimes_seq)}

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
  # if(on_Midway()){
  #   write_csv(nets$virus_ps_list, make_name('virus_ps_list',hr), col_names = T)
  #   if(!is.null(nets$bacteria_sp_list)){
  #     write_csv(nets$bacteria_sp_list, make_name('bacteria_sp_list',hr), col_names = T)
  #   }
  # if(!is.null(nets$immunity_list)){
  #   write_csv(nets$immunity_list, make_name('immunity_list',hr), col_names = T)
  # }
  #   write.csv(as.data.frame(nets$immunity_matrix), make_name('immunity_matrix',hr))
  #   if(!is.null(nets$infection_list)){
  #     write_csv(nets$infection_list, make_name('infection_list',hr), col_names = T)
  #   }
  #   write.csv(as.data.frame(nets$infection_matrix), make_name('infection_matrix',hr))
  #   write.csv(as.data.frame(nets$mutation_matrix), make_name('mutation_matrix',hr))
  #   write_csv(nets$bacteria_no_sp, make_name('bacteria_no_sp',hr), col_names = T)
  #   write_csv(nets$nodes, make_name('nodes',hr), col_names = T)
  # }
  # if(!is.null(nets$immunity_matrix)){
  #   write.csv(nets$immunity_matrix, make_name('immunity_matrix',hr))
  # }
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
if(on_Midway()){record_data(virus_density)}

# Data frame for bacteria density (abundance)
bacteria_density <- NULL
for (hr in hr_seq){
  tmp <- all_networks[[which(hr_seq==hr)]]$bacteria_abund_hr
  tmp$hr <- hr
  bacteria_density <- rbind(bacteria_density, tmp)
}
total_density <- bacteria_density %>% group_by(hr) %>% summarise(total_density_hr=sum(density))
bacteria_density %<>% left_join(total_density) %>% select(hr,B_ID,density,total_density_hr)
if(on_Midway()){record_data(bacteria_density)}
  
# Immunity network density and size
imm_density_size <- NULL
for (hr in hr_seq){
  x <- all_networks[[which(hr_seq==hr)]]$immunity_matrix
  d <- sum(x!=0)/(nrow(x)*ncol(x))
  imm_density_size <- rbind(imm_density_size, tibble(hr=hr, Density=d, Links=sum(x!=0), B=nrow(x), V=ncol(x), Size=nrow(x)+ncol(x)))
}
if(on_Midway()){record_data(imm_density_size)}
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
if(on_Midway()){record_data(richness)}

plt_richness_no_spacers <- 
  standard_plot(
    ggplot(richness %>% filter(node_type!='spacers'), aes(hr, n, color=node_type))+
    geom_line(size=1.5)+
    scale_color_manual(values = c('#27AE60','purple'))+
    theme(legend.position = 'none')+
    labs(y='Richness')
  )
make_png(plt_richness_no_spacers)
make_svg(plt_richness)

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
            
diversity_measures <- NULL
for (hr in hr_seq){
  x <- all_networks[[which(hr_seq==hr)]]
  H_bact <- vegan::diversity(x$bacteria_abund_hr$density)
  J_bact <- vegan::diversity(x$bacteria_abund_hr$density)/log(nrow(x$bacteria_abund_hr))
  H_virus <- vegan::diversity(x$virus_abund_hr$density)
  J_virus <- vegan::diversity(x$virus_abund_hr$density)/log(nrow(x$virus_abund_hr))
  
  # Hill_1_bact <- vegan::renyi(x$bacteria_abund_hr$density, scales = 1, hill = T)[1]
  # Hill_1_virus <- vegan::renyi(x$virus_abund_hr$density, scales = 1, hill = T)[1]
  # Hill_2_bact <- vegan::renyi(x$bacteria_abund_hr$density, scales = 2, hill = T)[1]
  # Hill_2_virus <- vegan::renyi(x$virus_abund_hr$density, scales = 2, hill = T)[1]
  # Hill_3_bact <- vegan::renyi(x$bacteria_abund_hr$density, scales = 3, hill = T)[1]
  # Hill_3_virus <- vegan::renyi(x$virus_abund_hr$density, scales = 3, hill = T)[1]
  
  out <- tibble(hr=hr, H_bact=H_bact, J_bact=J_bact, H_virus=H_virus,J_virus=J_virus)
                # Hill_1_bact=Hill_1_bact,
                # Hill_2_bact=Hill_2_bact,
                # Hill_3_bact=Hill_3_bact,
                # Hill_1_virus=Hill_1_virus,
                # Hill_2_virus=Hill_2_virus,
                # Hill_3_virus=Hill_3_virus)
  diversity_measures <- rbind(diversity_measures, out)
}
diversity_measures %<>% gather(key = 'measure', value = 'value', -hr) %>% 
  mutate(kind=case_when(str_detect(measure,'H_')~'Diversity',
                        str_detect(measure,'J_')~'Evenness',
                        str_detect(measure,'Hill_1')~'Hill 1',
                        str_detect(measure,'Hill_2')~'Hill 2',
                        str_detect(measure,'Hill_3')~'Hill 3'),
         entity=ifelse(str_detect(measure,'bact'),'Bacteria','Viruses'))

if(on_Midway()){record_data(diversity_measures)}

plt_diversity <- 
  standard_plot(
    ggplot(diversity_measures, aes(hr, value, color=entity))+geom_line()+
    scale_color_manual(values = c('blue','red'))+
    facet_grid(kind~., scales = 'free_y')+
    theme(legend.position = 'none')
  )
make_png(plt_diversity)
make_svg(plt_diversity)

if (complete_analysis){
  
  # Modularity of infection networks ----------------------------------------
  modules_df <- NULL
  for (hr in hr_seq){
    notify(paste('Infection networks (modularity) ',hr,'/',stop_time,sep=''))
    edges <- all_networks[[which(hr_seq==hr)]]$infection_list
    nodes <- all_networks[[which(hr_seq==hr)]]$nodes
    modules <- run_infomap(edges, remove_zero_edges = T)
    if (!is.null(modules)){
      suppressMessages(suppressWarnings(modules %<>% left_join(nodes)))
      modules_df <- rbind(modules_df, modules)
      if (make_plots){
        png(paste('plots/infection_modules_',str_pad(hr,4,'left','0'),'.png',sep=''),1920,1080,res=150)
        print(
          plot_modules(network_list = edges, m_df = modules, remove_zero_edges = T)+
            theme(axis.text = element_text(size=6),
                  axis.text.x = element_text(angle=-90))+
            labs(title=hr)
        )
        dev.off()
      }
    }
  }
  
  if(on_Midway()){write_csv(modules_df, 'modules_df_infection.csv', col_names = T)}
  
  plt_modules_infection <- 
    standard_plot(
      modules_df %>% group_by(hr) %>% summarise(num_mod=length(unique(module))) %>% 
        ggplot(aes(hr,num_mod))+
        geom_line(size=1.2)+
        labs(y='Number of modules')+
        theme(
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.position = 'none'
        )
    )
  
  make_png(plt_modules_infection)
  make_svg(plt_modules_infection)
  
  
  # Plot infection networks without modularity
  if (make_plots){
    for (hr in hr_seq){
      notify(paste('Infection networks (plots) ',hr,'/',stop_time,sep=''))
      x <- all_networks[[which(hr_seq==hr)]]$infection_matrix
      
      if (!all(x==0)){
        plt <- plot_matrix(t(log(x+1)), layout = 'nested', method = 'ggplot', binary_cols = c('gray',gg_color_hue(1)), legend_title='log(w)')$plt
        png(paste('plots/infection_nested_',str_pad(hr,4,'left','0'),'.png',sep=''),1920,1080,res=150)
        print(
          plt+theme(axis.text = element_text(size=6),
                    axis.text.x = element_text(angle=-90))+
            labs(title=hr)
        )
        dev.off()
      }
    }
  }

  # Modularity of bacteria-spacer networks ----------------------------------
  modules_df <- NULL
  for (hr in hr_seq){
    notify(paste('Bacteria-spacer networks (modularity) ',hr,'/',stop_time,sep=''))
    edges <- all_networks[[which(hr_seq==hr)]]$bacteria_sp_list
    nodes <- all_networks[[which(hr_seq==hr)]]$nodes
    modules <- run_infomap(edges, remove_zero_edges = F)
    if (!is.null(modules)){
      suppressMessages(suppressWarnings(modules %<>% left_join(nodes)))
      modules_df <- rbind(modules_df, modules)
      if (make_plots){
        png(paste('plots/bacteria_sp_',str_pad(hr,4,'left','0'),'.png',sep=''),1920,1080,res=150)
        print(
          plot_modules(network_list = edges, m_df = modules, remove_zero_edges = F)+
            labs(title=hr)+
            theme(axis.text = element_text(size=5))
        )
        dev.off()
      }
    }
  }
  
  if(on_Midway()){write_csv(modules_df, 'modules_df_bacteria_ps.csv', col_names = T)}
  
  plt_modules_BS <- 
    standard_plot(
      modules_df %>% group_by(hr) %>% summarise(num_mod=length(unique(module))) %>% 
        ggplot(aes(hr,num_mod))+
        geom_line()+
        labs(y='Number of modules',title = 'Number of modules in the bacteria-spacer networks')
    )
  make_png(plt_modules_BS)
  make_svg(plt_modules_BS)
  
  plt_modules_bs_infection <- plot_grid(plt_modules_BS,plt_modules_infection,nrow = 2, align = 'vh')
  make_png(plt_modules_bs_infection)
  make_svg(plt_modules_bs_infection)


    
  # Modularity of virus-protospacer networks --------------------------------
  modules_df <- NULL
  for (hr in hr_seq){
    notify(paste('Virus-protospacer networks (modularity) ',hr,'/',stop_time,sep=''))
    edges <- all_networks[[which(hr_seq==hr)]]$virus_ps_list
    nodes <- all_networks[[which(hr_seq==hr)]]$nodes
    modules <- run_infomap(edges)
    suppressMessages(suppressWarnings(modules %<>% left_join(nodes)))
    modules_df <- rbind(modules_df, modules)
    if (make_plots){
    png(paste('plots/virus_ps_',str_pad(hr,4,'left','0'),'.png',sep=''),1920,1080,res=150)
    print(
      plot_modules(edges, modules)+
      labs(title=hr)+
      theme(axis.text = element_text(size=5))
    )
    dev.off()
    }
  }
  if(on_Midway()){write_csv(modules_df, 'modules_df_virus_ps.csv', col_names = T)}
  plt_modules_VP <- 
    standard_plot(
      modules_df %>% group_by(hr) %>% summarise(num_mod=length(unique(module))) %>% 
      ggplot(aes(hr,num_mod))+
        geom_line()+
        labs(y='Number of modules',title = 'Number of modules in the virus-protospacer networks')
    )
  make_png(plt_modules_VP)
  make_svg(plt_modules_VP)
  
  
  
  
  
  # Modularity of immunity networks -----------------------------------------
  modules_df <- NULL
  for (hr in hr_seq){
    notify(paste('Immunity networks (modularity) ',hr,'/',stop_time,sep=''))
    edges <- all_networks[[which(hr_seq==hr)]]$immunity_list
    nodes <- all_networks[[which(hr_seq==hr)]]$nodes
    modules <- run_infomap(edges, remove_zero_edges = T)
    if (!is.null(modules)){
      suppressMessages(suppressWarnings(modules %<>% left_join(nodes)))
      modules_df <- rbind(modules_df, modules)
      if (make_plots){
        png(paste('plots/immunity_modules_',str_pad(hr,4,'left','0'),'.png',sep=''),1920,1080,res=150)
        print(
          plot_modules(edges, modules, remove_zero_edges = T)+
            labs(title=hr)+
            theme(axis.text = element_text(size=5))
        )
        dev.off()
      }
    }
  }
  if(on_Midway()){write_csv(modules_df, 'modules_df_immunity.csv', col_names = T)}
  plt_modules_immunity <- 
    standard_plot(
      modules_df %>% group_by(hr) %>% summarise(num_mod=length(unique(module))) %>% 
        ggplot(aes(hr,num_mod))+
        geom_line()+
        labs(y='Number of modules',title = 'Number of modules in the immunity networks')
    )
  make_png(plt_modules_immunity)
  make_svg(plt_modules_immunity)
  
  
  
} # End complete_analysis

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

if(on_Midway()){record_data(virus_dynamics_list)}

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

if(on_Midway()){record_data(bacteria_dynamics_list)}

plt_bacteria_persistence <- 
  standard_plot(
    ggplot(bacteria_dynamics_list)+
      geom_rect(aes(ymin=parse_number(bacteria_dynamics_list$B_ID),
                    ymax=parse_number(bacteria_dynamics_list$B_ID), 
                    xmin=birth, 
                    xmax=death), color='blue',
                size=0.7)+
      labs(y='Bacteria ID')+
      theme(legend.position = 'none')
  )
make_png(plt_bacteria_persistence)
make_svg(plt_bacteria_persistence)

# Trees -------------------------------------------------------------------
if (complete_analysis){
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
}


# Significance of modularity of host-spacer networks --------------------------------------
# hs_modularity_pvalues <- NULL
# # for (hr in hr_seq){
# # Test for significance within VDRs
# for (hr in VDR_hrs){
#   edges <- all_networks[[which(hr_seq==hr)]]$bacteria_sp_list
#   if(is.null(edges)){next}
#   if(nrow(edges)<50){next}
#   x <- create_monolayer_object(edges, directed = F, bipartite = T)
#   test <- run_infomap_monolayer_nonrandom(x, infomap_executable = 'Infomap', flow_model = 'undirected', silent = T, trials = 100, two_level = T, seed = 123, signif = T, shuff_method = 'r00', nsim = 100)
#   notify(paste('Significance of host-spacer networks (modularity) ',hr,'/',stop_time,max(test$modules$module_level1),' modules | pvalue:',test$pvalue,sep=''))
#   hs_modularity_pvalues <- bind_rows(hs_modularity_pvalues, tibble(hr=hr,
#                                                                    n_hosts=ncol(x$mat),
#                                                                    n_spacer=nrow(x$mat),
#                                                                    n_interactions=sum(x$mat>0),
#                                                                    n_modules=max(test$modules$module_level1),
#                                                                    pvalue=ifelse(class(test$pvalue)!='numeric',0,test$pvalue)))
# }
# record_data(hs_modularity_pvalues)
# # hs_modularity_pvalues %>% 
# mutate(signif=ifelse(pvalue<0.05,T,F)) %>% 
#   ggplot(aes(x=hr, y=pvalue))+geom_line()+geom_point(aes(color=signif))+geom_hline(yintercept = 0.05)


# Agregate host-spacer networks within each VDR and test for significance
hs_modularity_aggregated <- NULL
library(infomapecology)
# install_infomap()
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
                                signif = T, shuff_method = 'r00', nsim = 100)
  
  print(host_sp_modularity$pvalue) # Is modularity significant?
  
  plot_modular_matrix(host_sp_modularity)
  
  tibble(L_sim=test$L_sim) %>%
    ggplot(aes(L_sim))+
    geom_histogram(fill='plum')+
    geom_vline(xintercept = test$L, linetype='dashed')+
    theme_bw()+
    labs(x='Map equation L', y='Count')+
    theme(legend.position='none', axis.text = element_text(size=20), axis.title = element_text(size=20))
  
  
  # Phylogenetic signal analysis
  tree <- tree_bacteria
  D <- ape::cophenetic.phylo(tree) # Phyloegentic distance
  D <- matrix_to_list_unipartite(D, directed = T) # Use directed to make sure that the from column has all the nodes (need it for joining later)
  D <- D$edge_list
  
  ## @knitr strains_tree_network
  # Difference between tree and matrix
  nodes_in_modules <- host_sp_modularity$modules %>%
    filter(str_starts(node_name, 'B')) %>%
    distinct(node_name) %>% 
    mutate(node_name=str_replace_all(node_name, pattern = '\\.', ''))
  nodes_in_modules <- nodes_in_modules$node_name
  nodes_in_tree <- tree$tip.label
  print(setdiff(nodes_in_modules, nodes_in_tree)) # In modules but not in tree
  print(setdiff(nodes_in_tree, nodes_in_modules)) # In tree but not in modules
  # Overlapping nodes:
  (overlapping <- intersect(nodes_in_tree, nodes_in_modules))
  ## @knitr END
  
  
  ## @knitr PD_within_modules
  # Observed modules
  M_obs <- host_sp_modularity$modules %>%
    filter(str_starts(node_name, 'B')) %>% 
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
  ## @knitr END
  
  ## @knitr permute_modules
  #Shuffle to create permuted modules of the same size,
  #and recalculate the meand PD within modules. The shuffling permutes the ID of the strains.
  D_perm <- NULL
  nperm <- 500
  for (i in 1:nperm){
    print(i)
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
  print(D_perm)
  ## @knitr END
  
  # Null hypothesis is that the permuted distance is smaller than the observed for
  # each module (i.e., no signal). If we reject this hypothesis then there is
  # phylogenetic signal because the observed PD beteween hosts within each module
  # would be smaller than expected by chance (closely related hosts share a module).
  
  ## @knitr test_significance_PD_across_modules
  # Plot the means 
  D_perm %>% group_by(run) %>% 
    summarise(D_perm_mean=mean(d_mean)) %>% 
    ggplot(aes(x=D_perm_mean))+geom_histogram()+geom_vline(xintercept = D_obs_mean)
  
  D_perm %>% group_by(run) %>% 
    summarise(D_perm=mean(d_mean)) %>% 
    mutate(test=D_perm<D_obs_mean) %>%
    summarise(pvalue=sum(test)/nperm) %>% 
    mutate(res=ifelse(pvalue<0.05,'Signal','No signal'))
  ## @knitr END
  
  ## @knitr test_significance_PD_within_modules
  # This can also be tested per module
  observations <- nperm*nrow(D_obs)
  D_perm %>% 
    full_join(D_obs, by='m') %>% 
    rename(d_perm=d_mean.x, d_obs=d_mean.y) %>% 
    ggplot(aes(x=d_perm))+
    geom_histogram()+
    facet_wrap(~m)+
    geom_vline(data = D_obs, aes(xintercept = d_mean))
  
  D_perm %>% 
    full_join(D_obs, by='m') %>% 
    rename(d_perm=d_mean.x, d_obs=d_mean.y) %>% 
    mutate(test=d_perm<d_obs) %>%  
    group_by(m) %>% 
    summarise(pvalue=sum(test)/nperm) %>% 
    mutate(Signif=ifelse(pvalue<0.05,'Signal','No signal'),
           Signif_Bonferroni=ifelse(pvalue<0.05/nrow(D_obs),'Signal','No signal')) # Need to divide by number of modules for Bonferroni correction
  ## @knitr END
}
record_data(hs_modularity_aggregated)




### Correlate phylogenetic distance to cooccurrence (similarity in module association)
# D <- ape::cophenetic.phylo(tree_bacteria) # Phyloegentic distance
# # Build data for cooccurrence
# c_dat <- modules_df %>% 
#   # filter(hr>200, hr<300) %>% 
#   filter(str_starts(nodeName, 'B')) %>% 
#   dplyr::select(nodeName, mod_hr) %>% 
#   mutate(in_module=1) %>% 
#   spread(key = mod_hr, value = in_module, fill = 0)
# rnames <- c_dat$nodeName
# c_dat <- data.matrix(c_dat[,-1])
# rownames(c_dat) <- rnames
# dim(c_dat)
# # Jaccard index for coocurrence (similarity in module affiliation)
# C <- as.matrix(vegan::vegdist(c_dat, method = 'jaccard'))
# # Order phylogeny and distance
# dim(D)
# dim(C)
# D_nodes <- rownames(D)
# C_nodes <- rownames(C)
# setdiff(D_nodes,C_nodes)
# setdiff(C_nodes,D_nodes)
# overlapping_nodes <- intersect(D_nodes,C_nodes)
# C <- C[overlapping_nodes, overlapping_nodes]
# D <- D[overlapping_nodes, overlapping_nodes]
# dim(C);dim(D)
# all(rownames(D)==rownames(C))
# 
# cor.test(D[lower.tri(D)], C[lower.tri(C)], method = 'pearson')
# ggplot(data = tibble(phylo=D[lower.tri(D)], cooccur=C[lower.tri(C)]), aes(x=phylo, y=cooccur))+geom_point()





# Significance of modularity of infection networks --------------------------------------

# Test at the end of each VDR
infection_modularity <- NULL
for (i in 1:nrow(VDRs)){
  hr <- VDRs$end[i]
  edges <- all_networks[[which(hr_seq==hr)]]$infection_list
  x <- create_monolayer_object(edges, directed = F, bipartite = T)
  test <- run_infomap_monolayer(x, infomap_executable = 'Infomap', flow_model = 'undirected', silent = T, trials = 100, two_level = T, seed = 123, signif = T, shuff_method = 'r00', nsim = 100)
  infection_modularity %<>% bind_rows(tibble(hr=hr, 
                                             pvalue=test$pvalue,
                                             n_hosts=ncol(x$mat),
                                             n_spacer=nrow(x$mat),
                                             n_interactions=sum(x$mat>0),
                                             n_modules=max(test$modules$module_level1)))
  
  plot_modular_matrix(test)
  tibble(L_sim=test$L_sim) %>%
    ggplot(aes(L_sim))+
    geom_histogram(fill='plum')+
    geom_vline(xintercept = test$L, linetype='dashed')+
    theme_bw()+
    labs(x='Map equation L', y='Count')+
    theme(legend.position='none', axis.text = element_text(size=20), axis.title = element_text(size=20))
  
  pd_results <- test_PD_modules(tree = tree_bacteria, module_object = test, node_start_letter = 'B')
  pd_results$D_obs
  pd_results$result_across_moduels
  pd_results$result_within_moduels
}
record_data(infection_modularity)





# Nestedness with EV ------------------------------------

# We use the leading eigenvalue (rho) as a measurement for nestedness, following:
# Staniczenko PPA, Kopp JC, Allesina S. The ghost of nestedness in ecological
# networks. Nat Commun. 2013;4: 1391. Because this measure is not comparable for
# networks of different sizes we also calculate the upper bound for the
# eigenvalue (see Stanić Z. Inequalities for Graph Eigenvalues. Cambridge
# University Press; 2015, page 51). If rho is equal to the upper
# bound that means that the network is organized in the most nested way possible
# and so the full "potential" of the nestedness is realized, giving a ratio of
# rho/upper_bound=1. This method avoids the necessity to randomize
# networks test for significance because if a network is organized in the most
# possible nested way, it will naturally have higher nestedness (larger leading ev) than any
# randomized version.

ev_nestedness <- NULL

for (hr in hr_seq){
  x <- all_networks[[which(hr_seq==hr)]]$immunity_matrix
  notify(paste('Immunity networks (EV nestedness) ',hr,'/',stop_time,' | dimensions: ',nrow(x),' by ',ncol(x),sep=''))
  if (nrow(x)<2 || ncol(x)<2 || all(x==0)){
    rho <- NA
    ev_upper <- NA
  } else {
    rho <- calculate_ev_nestedness(x)
    ev_upper <- calculate_upper_bound(x)
  }
  ev_nestedness <- rbind(ev_nestedness, tibble(hr=hr, 
                                               rho=rho, 
                                               ev_upper=ev_upper,
                                               nestedness_potential=rho/ev_upper
                                               ))
  
}

if(on_Midway()){record_data(ev_nestedness)}

plt_ev_nestedness <- 
  standard_plot(
    ggplot(ev_nestedness, aes(hr, rho))+
    geom_line(color='#1E5984', size=1.1)+
    labs(y=expression(lambda[max]))
  )
make_png(plt_ev_nestedness)
make_svg(plt_ev_nestedness)

plt_potential_nestedness <- ggplot(ev_nestedness, aes(x=nestedness_potential))+
  geom_density(fill='gray')+
  labs(y='Density', x='Potential nestedness')
make_png(plt_potential_nestedness)
make_svg(plt_potential_nestedness)

plt_ev_nestedness_w_inset <- 
plt_ev_nestedness + 
  annotation_custom(
    ggplotGrob(plt_potential_nestedness+
                 theme(axis.text = element_text(size=8),
                 axis.title = element_text(size=8))
               ), 
    xmin = 0.75*max(ev_nestedness$hr), xmax = 0.98*max(ev_nestedness$hr), ymin = 0.75*max(ev_nestedness$rho, na.rm = T), ymax = 0.98*max(ev_nestedness$rho, na.rm = T)
  )
make_png(plt_ev_nestedness_w_inset)
make_svg(plt_ev_nestedness_w_inset)


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
if(on_Midway()){record_data(WNODF_df)}
plt_WNODF <- 
  standard_plot(
    ggplot(WNODF_df, aes(hr, WNODF))+
      geom_line(color='#1E5984', size=1.1)+
      labs(y='WNODF')
  )
make_png(plt_WNODF)
make_svg(plt_WNODF)

plt_nestedness_comparison <- 
  standard_plot(
    WNODF_df %>% 
      inner_join(ev_nestedness) %>% 
      select(hr, rho, WNODF) %>% 
      gather(key = 'nestedness', value = 'value', -hr) %>% 
      ggplot(aes(hr, value, color=nestedness))+
        geom_line()+
        facet_grid(nestedness~.,scales = 'free_y')
)
make_png(plt_nestedness_comparison)
make_svg(plt_nestedness_comparison)
# 
# nestedness_comparison_VDR <- 
#   WNODF_df %>%
#   inner_join(ev_nestedness) %>% 
#   select(hr, rho, WNODF) %>% 
#   filter(!is.na(rho)) %>% 
#   filter(hr%in%VDR_hrs)
# ccf(nestedness_comparison_VDR$rho, nestedness_comparison_VDR$WNODF)
# 
# nestedness_comparison_BDR <- 
#   WNODF_df %>% inner_join(ev_nestedness) %>% 
#   select(hr, rho, WNODF) %>% 
#   filter(!is.na(rho)) %>% 
#   filter(hr%in%BDR_hrs)
# ccf(nestedness_comparison_BDR$rho, nestedness_comparison_BDR$WNODF)

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
if(on_Midway()){record_data(spacer_matches)}

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
# plt_spacer_matches_01 <-
#   standard_plot(
#     spacer_matches %>% 
#       filter(num_matches%in%c(0,1)) %>% 
#       left_join(match_cols, by=c('num_matches'='value')) %>% 
#       mutate(col=ifelse(num_matches>6,'black',col)) %>% # more than 6 matches are colored black
#       ggplot(aes(x=hr, prop, color=col))+geom_line(size=1.5)+
#       scale_color_identity('# matches',drop=FALSE, guide = 'legend', labels=match_cols$value, breaks=match_cols$col)+
#       labs(y='% matches')+
#       theme(
#             axis.text = element_text(size = 16),
#             axis.title = element_text(size = 16),
#             legend.position = 'none'
#             )
#   )
make_png(plt_spacer_matches)
make_svg(plt_spacer_matches)


# Filled spacers/protospacers ---------------------------------------------
spacers_per_host <- protosp_per_virus <- NULL
for (hr in hr_seq){
  notify(paste('Spacers/protospacers in hosts/viruses ',hr,'/',stop_time,sep=''))
  
  # x <- all_networks[[which(hr_seq==hr)]]$virus_ps_list
  # tmp <- x %>% group_by(V_ID) %>% count()
  # tmp$prop <- tmp$n/protospacer_len
  # tmp$hr <- hr
  # protosp_per_virus <- rbind(protosp_per_virus, tmp)
  
  x <- all_networks[[which(hr_seq==hr)]]$bacteria_sp_list
  if (is.null(x)){next}
  x  <- subset(x, w!=0)
  tmp <- x %>% group_by(B_ID) %>% count()
  tmp$prop <- tmp$n/spacer_len
  tmp$hr <- hr
  spacers_per_host <- rbind(spacers_per_host, tmp)
}

plt_spacers_per_host <- standard_plot(
spacers_per_host %>%
  group_by(hr) %>% summarise(mean_n=mean(n),
                                   min_n=min(n),
                                   max_n=max(n),
                                   n=n()) %>% 
  ggplot+
  geom_line(aes(x=hr, y=mean_n))+
  geom_line(aes(x=hr, y=min_n), linetype='dashed')+
  geom_line(aes(x=hr, y=max_n), linetype='dashed')+
  labs(x='Time', y='Mean number of spacers filled')
)
make_png(plt_spacers_per_host)
make_svg(plt_spacers_per_host)

# spacers_per_host$type <- 'host'
# protosp_per_virus$type <- 'virus'
# names(protosp_per_virus)[1] <- names(spacers_per_host)[1] <- 'ID'
# d <- ungroup(spacers_per_host) %>% bind_rows(ungroup(protosp_per_virus))
# d %<>% group_by(type, hr) %>% summarise(mean_prop=mean(prop),
#                                        min_prop=min(prop),
#                                        max_prop=max(prop),
#                                        n=n())
# d %>% ggplot+
#   geom_line(aes(x=hr, y=mean_prop, color=type))+
#   geom_line(aes(x=hr, y=min_prop, color=type, linetype='dashed'))+
#   geom_line(aes(x=hr, y=max_prop, color=type, linetype='dashed'))+
#   facet_grid(~type)
  


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
if(on_Midway()){record_data(ext_immunity_rank)}

plt_extinction_rank <- ggplot(ext_immunity_rank, aes(ext_immunity_rank))+
  geom_histogram(fill=NA, color='tomato')+
  geom_vline(xintercept = median(virus_dynamics_list$ext_immunity_rank, na.rm = T), linetype='dashed')+
  labs(x='Immunity rank of extinct viruses')

make_png(plt_extinction_rank)
make_svg(plt_extinction_rank)


virus_dynamics_list %<>%
  left_join(ext_immunity_rank) %>%
  mutate(high_rank=ifelse(ext_immunity_rank>median(ext_immunity_rank),T,F))


# Plots of extinction events
if(on_Midway() && make_plots){
  for (hr in hr_seq){
    notify(paste('Plot extinction events ',hr,'/',stop_time,sep=''))
    x <- all_networks[[which(hr_seq==hr)]]$immunity_matrix
    v_ext <- subset(virus_dynamics_list, death==hr)$V_ID
    if (length(v_ext)>0){
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
      png(paste('plots/immunity_extinction_event_',str_pad(hr,4,'left','0'),'.png',sep=''),1920,1080,res=150)
      print(
        M %>% 
          ggplot()+
          geom_tile(aes(Var1,Var2,fill=col))+
          geom_tile(data=M %>% filter(Var1%in%v_ext), aes(Var1,Var2), fill=NA, color='red', size=0.6)+
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
    } else {print('no extinction event')}
  }
}


# Effective mutation --------------------------------------------
mu_eff_df <- NULL
notify('--> Calculate effective mutation...')
for (hr in hr_seq){
  notify(paste('Effective mutation at time ',hr,'/',stop_time,sep=''))
  x <- all_networks[[which(hr_seq==hr)]]$mutation_matrix
  mu_eff_df <- rbind(mu_eff_df, tibble(hr=hr, V_ID=colnames(x), mu_eff_j=colSums(x)))
}
mu_eff_total <- mu_eff_df %>% group_by(hr) %>% summarise(mu_total_hr=sum(mu_eff_j)) # Total effective mutation rate across viruses for any hour
mu_eff_df %<>% left_join(mu_eff_total)

if(on_Midway()){record_data(mu_eff_df)}

plt_mu_eff <- 
  standard_plot(
    ggplot(mu_eff_df, aes(hr, mu_total_hr))+
    geom_line(color='brown')+
    labs(y=expression(mu[eff]))+
    labs(title='Total effective mutation rate (across viruses)')
  )
make_png(plt_mu_eff)
make_svg(plt_mu_eff)


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
  
  # R_0m <- empty(R_0m)
  
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

if(on_Midway()){record_data(R_0m_df)}

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
  # R_1m <- empty(R_1m)
  
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

if(on_Midway()){record_data(R_1m_df)}

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

if(on_Midway()){record_data(R_pot_df)}

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


# Distribution of Rpot of all the viruses
Q90 <- quantile(R_pot_df$Rpot, probs = 0.9)
plt_Rpot_PDF <- ggplot(R_pot_df, aes(Rpot))+geom_density(fill='red',alpha=0.7)+
  labs(y='Density', x=expression(R[pot]))+
  geom_vline(xintercept = Q90, linetype='dashed')
# Calculate what proportion of the virus density (abundance) is concentrated in
# viruses with Rpot above the 90% percentile.
R_pot_rel_density <- R_pot_df %>% filter(Rpot>Q90) %>% group_by(hr) %>% summarise(D=sum(rel_density))
plt_Rpot_density <- 
standard_plot(
  R_pot_rel_density %>% 
    ggplot(aes(hr,D))+geom_line()+
    labs(y='Relative density',title = 'Relative density of viruses concentrated in viruses with Rpot higher than the 90% percentile')
)+scale_x_continuous(limits = range(hr_seq), breaks=label_seq)
make_png(plt_Rpot_density)
make_svg(plt_Rpot_density)



# Probability of escape ---------------------------------------------------
# This part is not used, but I kept the code, just in case.
# 
# P_E_df <- NULL
# for (t in hr_seq){
#   notify(paste('Calculating P_E for time',t))
#   imm <- all_networks[[which(hr_seq==t)]]$immunity_matrix
#   if (!1%in%imm){next}
#   imm[imm!=1] <- 0
#   
#   virus_abund_hr <- all_networks[[which(hr_seq==t)]]$virus_abund_hr
#   virus_abund_hr <- subset(virus_abund_hr, V_ID%in%colnames(imm))
#   V_T <- sum(virus_abund_hr$density)
#   virus_abund_hr$virus_rel_abund <- virus_abund_hr$density/V_T
#   P_E <- (imm%*%diag(virus_abund_hr$virus_rel_abund))/protospacer_len # This is to multiply columns by the vector
#   colnames(P_E) <- colnames(imm)
#   g <- graph.incidence(t(P_E), directed = F, weighted = T) # Need to transpose so the "from" would be viruses and the "to" would be bacteria
#   tmp <- as.tibble(igraph::as_data_frame(g, what = 'edges'))
#   names(tmp) <- c('V_ID', 'B_ID', 'PE_t')
#   tmp %<>% distinct(V_ID,PE_t)
#   tmp$hr <- t
#   
#   P_E_df <- rbind(P_E_df, tmp)
#   
#   # all_networks[[which(hr_seq==t)]]$PE_matrix <- P_E
# }
# if(on_Midway()){record_data(P_E_df)}

# PE <- P_E_df %>% group_by(hr) %>% summarise(P_E=sum(PE_t)) %>% 
#   left_join(escapes) %>% 
#   left_join(no_escapes)
# plt_PE <- ggplot()+
#   geom_line(data=PE, aes(hr, P_E),color='dark green')+
#   geom_point(data=PE %>% filter(!is.na(num_escaped_bact)) %>% filter(num_escaped_bact>0), aes(x=hr,y=P_E,color=num_escaped_bact), size=2)+
#   geom_point(data=PE %>% filter(escaped_bacteria=='none'), aes(x=hr,y=P_E), color='gray50', size=2)+
#   # scale_color_gradient2(low = "#DA4453", mid='blue',high = "#89216B")+
#   scale_color_gradientn(colours = gg_color_hue(max(PE$num_escaped_bact, na.rm = T)))+
#   labs(x='Time', y='Prob. escape')+
#   theme(legend.position = c(0.9,1),
#         legend.direction='horizontal',
#         legend.justification = c("right", "top"),
#         legend.title=element_blank())

if (complete_analysis){
  
  # Escape process --------------------------------
  # This part uses the virus phylogenetic tree to follow mutations.
  notify('--> Follow virus escape...')
  tree <- read_delim(paste(base_name,'_Phage-TREE.txt',sep=''), delim = '\t', 
                   col_names = c("Recordtime","V_ID","parent_id","creation_time"))
  tree %<>% filter(Recordtime<=max(hr_seq))
  tree %<>% mutate(V_ID=paste('V_',str_pad(V_ID, 4, 'left', '0'),sep=''),
                 parent_id=paste('V_',str_pad(parent_id, 4, 'left', '0'),sep=''))
  tree$parent_id[1] <- NA
  
  escape_process <- vector(mode = 'list', length = nrow(tree))
  escape_process[[1]] <- NA # Root has no escape
  escape_process_df <- NULL
  
  for (i in 2:nrow(tree)){
    hr <- tree$Recordtime[i]
    V <- tree$V_ID[i]
    parent_id <- tree$parent_id[i]
    x <- all_networks[[which(hr_seq==hr)]]$mutation_matrix
    parent_mu_eff <- colSums(x)[parent_id]
    parent_mu_eff_prop <- parent_mu_eff/sum(x)
    
    parent_spacers <- all_networks[[which(hr_seq==hr)]]$virus_ps_list %>% filter(V_ID==parent_id)
    parent_spacers <- parent_spacers$PS
    V_spacers <- all_networks[[which(hr_seq==hr)]]$virus_ps_list %>% filter(V_ID==V)
    if (nrow(V_spacers)==0){print('Newly mutated virus could not be found, skipping.');next}
    V_spacers <- V_spacers$PS
    mutated_spacers <- setdiff(parent_spacers,V_spacers)
    new_spacers <- setdiff(V_spacers, parent_spacers)
    
    # Which bacteria have the mutated spacer
    bacteria_spacers <- all_networks[[which(hr_seq==hr)]]$bacteria_sp_list
    if (is.null(bacteria_spacers)) {
      escaped_bacteria <- NULL
    } else {
      bacteria_spacers %<>% filter(SP%in%mutated_spacers)
      if (nrow(bacteria_spacers)==0) {
        escaped_bacteria <- NULL
        bacteria_spacers <- NULL
      } else {
        bacteria_spacers %<>% filter(w==1)
        bacteria_spacers <- bacteria_spacers$B_ID
        escaped_bacteria <- all_networks[[which(hr_seq==hr)]]$immunity_list %>% filter(B_ID%in%bacteria_spacers, V_ID==parent_id, w==1)
        escaped_bacteria <- escaped_bacteria$B_ID
      }
    }
    
    cat(i);cat(' (');cat(hr);cat(' hr)');cat(': \n')
    cat(paste('Mutated virus:','\t',parent_id,' --> ',V,sep=''));cat(' | Proportion of mu_eff: ');cat(round(parent_mu_eff_prop,3));cat('\n')
    cat(paste('Mutated spacer:','\t',mutated_spacers,' --> ',new_spacers,sep=''));cat('\n')
    cat('Bacteria with matching spacers: ');cat(bacteria_spacers);cat('\n')
    cat('Escaped from bacteria: ');cat(escaped_bacteria);cat('\n')
    cat('------------------');cat('\n')
    
    escape_process_df <- rbind(escape_process_df, 
        expand.grid(mutation_event=i,
                    hr=hr,
                    parent_id=parent_id,
                    child=V,
                    mutated_spacer=mutated_spacers,
                    new_spacers=new_spacers,
                    bacteria_matching_spacers=returnnull(bacteria_spacers),
                    escaped_bacteria=returnnull(escaped_bacteria))
    )
    
    if(on_Midway() && make_plots){
      # Plot with fixed colors for matches
      x <- all_networks[[which(hr_seq==hr)]]$immunity_matrix
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
      plt <- 
        M %>% 
        ggplot()+
        geom_tile(aes(Var1,Var2,fill=col))+
        # geom_tile(data=M %>% filter(Var1%in%parent_id, Var2%in%escaped_bacteria), aes(Var1,Var2), fill=NA, color='purple', size=0.8)+
        geom_tile(data=M %>% filter(Var1%in%parent_id, Var2%in%bacteria_spacers), aes(Var1,Var2), fill=NA, color='purple', size=0.8)+
        geom_tile(data=M %>% filter(Var1%in%V, Var2%in%escaped_bacteria), aes(Var1,Var2), fill=NA, color='black', size=0.8)+
        scale_fill_identity('# matches',drop=FALSE, guide = 'legend', labels=match_cols$value, breaks=match_cols$col)+
        coord_fixed()+
        theme(
          legend.position = 'none',
          axis.text.x=element_text(angle=-90),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size=6))+
        labs(title=paste('Mutation event ',i,' (hr=',hr,')',sep=''))
    
      png(paste('plots/mutation_',str_pad(i,4,'left','0'),'_',str_pad(hr,4,'left','0'),'hr.png',sep=''),1920,1080,res=150)
      print(plt)
      dev.off()
    }
     # Store results
    escape_process[[i]]$hr <- hr
    escape_process[[i]]$parent_id <- parent_id
    escape_process[[i]]$mutatnt_id <- V
    escape_process[[i]]$parent_mu_eff <- parent_mu_eff
    escape_process[[i]]$parent_mu_eff_prop <- parent_mu_eff_prop
    escape_process[[i]]$mutated_spacers <- mutated_spacers
    escape_process[[i]]$new_spacers <- new_spacers
    escape_process[[i]]$bacteria_spacers <- bacteria_spacers
    escape_process[[i]]$escaped_bacteria <- escaped_bacteria
  }
  escape_process_df <- as_tibble(escape_process_df)
  if(on_Midway()){record_data(escape_process_df)}
  
  # Add the information on escapes to the virus dynamics data frame
  virus_dynamics_list$num_escaped_bact <- NA
  virus_dynamics_list$mu_eff_prop <- NA
  virus_dynamics_list$mu_eff <- NA
  for (i in 2:length(escape_process)){
    x <- escape_process[[i]]
    hr <- x$hr
    mutatnt_id <- x$mutatnt_id
    num_escaped_bact <- length(unique(x$escaped_bacteria))
    mu_eff_prop <- x$parent_mu_eff_prop
    mu_eff <- x$parent_mu_eff
    virus_dynamics_list[virus_dynamics_list$V_ID==mutatnt_id,'num_escaped_bact'] <- num_escaped_bact
    virus_dynamics_list[virus_dynamics_list$V_ID==mutatnt_id,'mu_eff_prop'] <- mu_eff_prop
    virus_dynamics_list[virus_dynamics_list$V_ID==mutatnt_id,'mu_eff'] <- mu_eff
  }
  
  if(on_Midway()){record_data(virus_dynamics_list)}
  
  # The points in the beginning of each virus depicts how many bacteria THIS
  # PARTICULAR MUTATION escaped
  start_colors <- tibble(num_escaped_bact=sort(unique(virus_dynamics_list$num_escaped_bact)))
  start_colors$start_col <- c('gray50',gg_color_hue(nrow(start_colors)-1))
  
  plt_virus_persistence_with_mutations <- 
    standard_plot(
      virus_dynamics_list %>% 
        left_join(start_colors, by = "num_escaped_bact") %>%
        ggplot()+
        geom_rect(aes(ymin=parse_number(virus_dynamics_list$V_ID),
                      ymax=parse_number(virus_dynamics_list$V_ID), 
                      xmin=birth, 
                      xmax=death), color='red',
                  size=0.7)+
        geom_point(aes(x=birth, y=parse_number(virus_dynamics_list$V_ID), color=start_col), size=1.5)+
        scale_color_identity(drop=FALSE, guide = 'legend', labels=start_colors$num_escaped_bact, breaks=start_colors$start_col)+
        labs(y='Virus ID')+
        theme(legend.position = c(0.1,0.9), legend.direction = 'horizontal', legend.title = element_blank(), legend.background = element_rect(fill = 'white'))
    )
  
  make_png(plt_virus_persistence_with_mutations)
  make_svg(plt_virus_persistence_with_mutations)
  
  plt_summary_virus_dynamics <- 
    plot_grid(
    # Plot the distribution of viruses persistence
    ggplot(virus_dynamics_list, aes(persistence))+
    geom_histogram(color='red', fill=NA)+
    geom_vline(xintercept = median(virus_dynamics_list$persistence, na.rm = T), linetype='dashed')+
    labs(x='Viruses persistence'),
  
    # Plot the distribution of bacteria persistence
    ggplot(bacteria_dynamics_list, aes(persistence))+
      geom_histogram(color='blue', fill=NA)+
      geom_vline(xintercept = median(virus_dynamics_list$persistence, na.rm = T), linetype='dashed')+
      labs(x='Bacteria persistence'),
    
    # Plot the distribution of ranks of viruses that went extinct
    ggplot(virus_dynamics_list, aes(ext_immunity_rank))+
    geom_histogram(fill=NA, color='tomato')+
    geom_vline(xintercept = median(virus_dynamics_list$ext_immunity_rank, na.rm = T), linetype='dashed')+
    labs(x='Immunity rank of extinct viruses'),
  
    # Plot the distribution of number of bacteria escaped
    ggplot(virus_dynamics_list, aes(num_escaped_bact))+
    geom_histogram(color='#B9770E', fill=NA)+
    geom_vline(xintercept = median(virus_dynamics_list$num_escaped_bact, na.rm = T), linetype='dashed')+
    labs(x='Number of escaped bacteria'),
  
    nrow=2,ncol=2)
  
  make_png(plt_summary_virus_dynamics)
  make_svg(plt_summary_virus_dynamics)
  
}

# Plot -----------------------------------------------------------


if (complete_analysis){
notify('--> Generate final plots...')

 # A PDF with all the main plots
pdf(paste(base_name,'_main_figures.pdf',sep=''), 16, 10, onefile = T)
  title <- ggdraw() + draw_label("Viral and bacterial strain abundance", fontface='bold')
  grid.arrange(plot_grid(title, plt_abundance_profiles, ncol=1, rel_heights=c(0.1, 1)))
  grid.arrange(plt_richness+labs(title = 'Richness of viruses bacteria and spacers'))
  title <- ggdraw() + draw_label("Bacteria abundance and richness", fontface='bold')
  grid.arrange(plot_grid(title, plt_bacteria_richness_abundace, ncol=1, rel_heights=c(0.1, 1)))
  grid.arrange(plt_diversity+labs(title = 'Diversity of bacteria and viruses'))
  grid.arrange(plt_modules_bs_infection)
  grid.arrange(plt_spacer_matches+labs(title = 'Proportion of spacer matches'))
  grid.arrange(plt_WNODF+labs(title = 'Nestedness of immunity network'))
  grid.arrange(plt_extinction_rank+labs(title = 'Distribution of extinction ranks'))
  grid.arrange(plt_mu_eff)
  grid.arrange(plt_R0m)
  grid.arrange(plt_R1m)
  grid.arrange(plt_R_pot)
  grid.arrange(plt_R_pot_effect)
  grid.arrange(plt_Rpot_density)
  grid.arrange(plt_virus_persistence+labs(title = 'Diversification and persistence of viruses'))
  grid.arrange(plt_bacteria_persistence+labs(title = 'Diversification and persistence of bacteria'))
  grid.arrange(plt_viruses_tree+labs(title = 'Viruses phylogenetic tree'))
  if ('plt_bacteria_tree' %in% ls(pattern = 'plt_bacteria_tree')){
    grid.arrange(plt_bacteria_tree+labs(title = 'Bacteria phylogenetic tree'))
  }
  grid.arrange(plt_summary_virus_dynamics)
dev.off()
  
pdf(paste(base_name,'_figures.pdf',sep=''), 16, 10, onefile = T)
  grid.arrange(plt_abundance_profiles)
  grid.arrange(plt_richness)
  grid.arrange(plt_bacteria_richness_abundace)
  grid.arrange(plt_diversity)
  grid.arrange(plt_modules_bs_infection)
  grid.arrange(plt_immunity_network_size)
  grid.arrange(plt_spacer_matches)
  grid.arrange(plt_WNODF)
  grid.arrange(plt_ev_nestedness)
  grid.arrange(plt_ev_nestedness_w_inset)
  grid.arrange(plt_potential_nestedness)
  grid.arrange(plt_mu_eff)
  grid.arrange(plt_R0m)
  grid.arrange(plt_R1m)
  grid.arrange(plt_R_pot)
  grid.arrange(plt_R_pot_effect)
  grid.arrange(plt_Rpot_density)
  grid.arrange(plt_virus_persistence)
  grid.arrange(plt_virus_persistence_with_mutations)
  grid.arrange(plt_bacteria_persistence)
  grid.arrange(plt_viruses_tree)
  if ('plt_bacteria_tree' %in% ls(pattern = 'plt_bacteria_tree')){
    grid.arrange(plt_bacteria_tree)
  }
  grid.arrange(plt_summary_virus_dynamics)
dev.off()

writeLines('success!', paste(base_name,'_success.txt',sep=''))
notify('--------- DONE! ---------')

} else {


pdf(paste(base_name,'_figures.pdf',sep=''), 16, 10, onefile = T)
  grid.arrange(plt_abundance_profiles)
  grid.arrange(plt_richness)
  grid.arrange(plt_bacteria_richness_abundace)
  grid.arrange(plt_diversity)
  grid.arrange(plt_immunity_network_size)
  grid.arrange(plt_spacer_matches)
  grid.arrange(plt_WNODF)
  grid.arrange(plt_ev_nestedness)
  grid.arrange(plt_ev_nestedness_w_inset)
  grid.arrange(plt_potential_nestedness)
  grid.arrange(plt_mu_eff)
  grid.arrange(plt_R0m)
  grid.arrange(plt_R1m)
  grid.arrange(plt_R_pot)
  grid.arrange(plt_R_pot_effect)
  grid.arrange(plt_Rpot_density)
  grid.arrange(plt_virus_persistence)
  grid.arrange(plt_bacteria_persistence)
dev.off()

writeLines('success!', paste(base_name,'_success.txt',sep=''))
notify('--------- DONE! ---------')
  
}
