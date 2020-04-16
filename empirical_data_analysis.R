library(tidyverse)
library(magrittr)
library(sqldf)
library(cowplot)
library(igraph)
library(bipartite)
library(infomapecology)

# setwd('/home/shai/Documents/CRISPR/Data/')
setwd('/Users/Shai/GitHub/ecomplab/CRISPR_networks/')
db <- dbConnect(SQLite(), dbname = '/Users/Shai/GitHub/ecomplab/CRISPR_networks/CRISPR_database_V2.sqlite')

# install_infomap(target_folder = getwd())

## @knitr load_functions

# Functions ---------------------------------------------------------------


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

bipartite_density <- function(M){
  sum(M!=0)/(nrow(M)*ncol(M))
}

plot_nested_matrix <- function(M){
  tmp <- as.matrix(M[order(rowSums(M), decreasing = T), order(colSums(M),decreasing = F)])
  rnames <- rownames(M)[order(rowSums(M), decreasing = T)]
  cnames <- colnames(M)[order(colSums(M), decreasing = F)]
  rownames(tmp) <- rnames
  colnames(tmp) <- cnames
  M <- tmp
  M <- as_tibble(reshape2::melt(M))
  num_matches <- sort(unique(M$value))
  match_colors <- tibble(value=num_matches, col=c('white',metafolio::gg_color_hue(length(num_matches)-1)))
  plt <- M %>% left_join(match_colors) %>% 
    ggplot()+
    geom_tile(aes(Var1,Var2,fill=col))+
    labs(x='Virus strain',y='Host strain')+
    scale_fill_identity('# matches',drop=FALSE, guide = 'legend', labels=match_colors$value, breaks=match_colors$col)+
    theme(
      # legend.position = 'none',
      axis.text.x = element_text(angle=-90),
      # axis.ticks = element_blank(),
      axis.title = element_text(size = 18))+
    coord_fixed()
  return(plt)
}


# A function to get the spacer set of a virus or a host from the host-spacer/virus-protospacer matrix
get_strain_spacers <- function(z, id, method) {
  stopifnot(method%in%c('dataframe','matrix'))
  if (method=='matrix'){
    x <- z[,id]
    return(names(which(x!=0)))
  } else {
    # From a dataframe
    if ('strain_id'%in%names(z)){ # If a host data frame
      x <- z %>% filter(strain_id==id)
    } else {  # If a virus data frame
      x <- z %>% filter(virus_strain==id)
    }
    return(x$spacer)
  }
}

# Match spacers and protospacers and produce an immunity network
get_matches <- function(host_data,virus_data, method){
  stopifnot(method%in%c('dataframe','matrix'))
  matching_spacers <- NULL
  if (method=='dataframe'){
    host_virus_pairs <- expand.grid(strain_id=unique(host_data$strain_id), virus_strain=unique(virus_data$virus_strain), matches=NA)
    for (p in 1:nrow(host_virus_pairs)){
      spacer_set <- get_strain_spacers(host_data, host_virus_pairs[p,1], method = 'dataframe')
      protospacer_set <- get_strain_spacers(virus_data, host_virus_pairs[p,2], method = 'dataframe')
      host_virus_pairs[p,3] <- length(intersect(protospacer_set,spacer_set))
      matching_spacers <- unique(c(matching_spacers,intersect(protospacer_set,spacer_set)))
    }
  } else {
    host_virus_pairs <- expand.grid(strain_id=colnames(host_data), virus_strain=colnames(virus_data), matches=NA)
    for (p in 1:nrow(host_virus_pairs)){
      spacer_set <- get_strain_spacers(host_data, host_virus_pairs[p,1], method = 'matrix')
      protospacer_set <- get_strain_spacers(virus_data, host_virus_pairs[p,2], method = 'matrix')
      host_virus_pairs[p,3] <- length(intersect(protospacer_set,spacer_set))
      matching_spacers <- unique(c(matching_spacers,intersect(protospacer_set,spacer_set)))
    }
  }
  return(list(immunity_network=host_virus_pairs, matching_spacers=as.character(sort(as.numeric(matching_spacers)))))
}

# Function to calculate significance of a weighted-nested networks

weighted_nestedness_significance <- function(M, weighted=F, nsim=10^3, shuff_method = 'r00_samp', make_plots=T){
  
  print('Observed nestedness...')
  ev_obs <- calculate_ev_nestedness(M)
  
  if (weighted){
    WNODF_obs <- bipartite::networklevel(M, index = 'weighted NODF')
  } else {
    WNODF_obs <- bipartite::networklevel(M, index = 'NODF')
  }
  
  # Make shuffled versions by shuffling the virus-protospacer network
  print('Shuffling...')
  null <- vegan::nullmodel(M, shuff_method)
  shuffled <- simulate(null, nsim = nsim)
  
  # If randomizing creates rows or cols that sum to 0 ignore the matrix because this seems to halt the ev calculation
  # any(apply(shuffled, MARGIN = 3, FUN = rowSums)==0) || any(apply(shuffled, MARGIN = 3, FUN = colSums)==0)
  print('Shuffled nestedness...')
  ev_shuffled <- tibble(ev_shuffled=apply(shuffled, MARGIN = 3, FUN = calculate_ev_nestedness))
  
  if (weighted){
    WNODF_shuffled <- tibble(WNODF_shuffled=apply(shuffled, MARGIN = 3, FUN = function (z) bipartite::networklevel(z, index = 'weighted NODF')))
  } else {
    WNODF_shuffled <- tibble(WNODF_shuffled=apply(shuffled, MARGIN = 3, FUN = function (z) bipartite::networklevel(z, index = 'NODF')))
  }
  
  # Calculate observed matrices
  p_value_ev <- sum(ev_shuffled>=ev_obs)/nsim
  p_value_WNODF <- sum(WNODF_shuffled>=WNODF_obs)/nsim
  print(paste('P-value (rho): ',p_value_ev,sep=''))
  print(paste('P-value (NODF or WNODF): ',p_value_WNODF,sep=''))
  
  out <- list(ev_obs=ev_obs,
              ev_shuffled=ev_shuffled$ev_shuffled,
              p_value_ev=p_value_ev,
              WNODF_obs=WNODF_obs,
              WNODF_shuffled=WNODF_shuffled$WNODF_shuffled,
              p_value_WNODF=p_value_WNODF
  )
  
  if (make_plots){
    print('Plotting...')
    
    p_shuffled_ev <- ggplot(ev_shuffled, aes(ev_shuffled))+geom_histogram(fill='steelblue')+
      geom_vline(xintercept = ev_obs, linetype='dashed')+labs(x='Shuffled rho')
    p_shuffled_WNODF <- ggplot(WNODF_shuffled, aes(WNODF_shuffled))+geom_histogram(fill='navy')+
      geom_vline(xintercept = WNODF_obs, linetype='dashed')+labs(x='Shuffled WNODF')
    out$plot_ev <- p_shuffled_ev
    out$plot_WNODF <- p_shuffled_WNODF
  }
  
  return(out)  
}


# Function for main analysis ----------------------------------------------

main <- function(dataset_id, nsim=10, font_size=20){
  # Get data
  
  # The criteria for matches matters A LOT !!!
  # criterium <- '0mm/-PAM'
  # criterium <- '4mm/-PAM'
  criterium <- '4mm/-PAM'
  
  as_tibble(dbGetQuery(db, paste('SELECT name FROM data_sets WHERE id=',dataset_id,sep='')))
  data_virus <- as_tibble(dbGetQuery(db, paste('SELECT virus_strain, spacer, criteria FROM virus WHERE dataset_id=',dataset_id,sep='')))
  data_virus %<>% filter(criteria==criterium) %>% select(-criteria) %>% mutate(w=1)
  data_host <- as_tibble(dbGetQuery(db, paste('SELECT strain_id, spacer FROM bacteria WHERE dataset_id=',dataset_id,sep='')))
  data_host %<>% mutate(w=1)
  
  # Create monolayer objects for infomapecology
  data_virus <- create_monolayer_object(data_virus, directed = F, bipartite = T, group_names = names(data_virus)[1:2])
  data_host <- create_monolayer_object(data_host, directed = F, bipartite = T, group_names = names(data_host)[1:2])
  
  # Is host-spacer network modular?
  host_sp_shuffled <- shuffle_infomap(data_host, shuff_method = 'r00', nsim = nsim, burnin=1000)
  host_sp_modularity <- run_infomap_monolayer(x = data_host, infomap_executable = 'Infomap', flow_model = 'undirected', silent = T, trials = 100, two_level = T, seed = 109743, signif = T, shuff_method = host_sp_shuffled, nsim = NULL)
  print(host_sp_modularity$pvalue)
  p1 <- plot_modular_matrix(host_sp_modularity, fix_coordinates=F, axes_titles=c('Host', 'Spacer'), transpose = T)
  p1 <- p1+theme_bw()+theme(legend.position='none', 
                            panel.grid = element_blank(),
                            axis.ticks = element_blank(),
                            axis.text = element_blank(), 
                            axis.title = element_text(size=font_size))
  
  p2 <- tibble(L_sim=host_sp_modularity$L_sim) %>% 
    ggplot(aes(L_sim))+
    geom_histogram(fill='plum')+
    geom_vline(xintercept = host_sp_modularity$L, linetype='dashed')+
    theme_bw()+
    labs(x='Map equation L', y='Count')+
    theme(legend.position='none', axis.text = element_text(size=font_size), axis.title = element_text(size=font_size))
  
  
  
  # Get immunity network
  matches_df <- get_matches(host_data = data_host$edge_list, virus_data = data_virus$edge_list, method = 'dataframe')
  immunity_network <- create_monolayer_object(matches_df$immunity_network %>% filter(matches!=0), directed = F, bipartite = T, group_names = c('host','virus'))
  
  # Is the immunity network nested?
  immunity_nestedness <- weighted_nestedness_significance(immunity_network$mat, nsim = nsim, weighted = T, shuff_method = 'r00_samp')
  print(immunity_nestedness[1:6])
  p3 <- immunity_nestedness$plot_WNODF+
    theme_bw()+
    theme(legend.position='none', axis.text = element_text(size=font_size), axis.title = element_text(size=font_size))
  
  p4 <- immunity_nestedness$plot_ev+
    theme_bw()+
    theme(legend.position='none', axis.text = element_text(size=font_size), axis.title = element_text(size=font_size))
  
  p5 <- plot_nested_matrix(immunity_network$mat)
  return(list(host_sp_modularity=host_sp_modularity,immunity_nestedness=immunity_nestedness,p1=p1,p2=p2,p3=p3,p4=p4,p5=p5))
}

## @knitr END


# Analysis ------------------------------------------------------------------

## @knitr analysis
Yellowstone <- main(dataset_id = 1, nsim = 250, font_size = 10)
Pseudomonas <- main(dataset_id = 3, nsim = 250, font_size = 10)
Russia2010 <- main(dataset_id = 6, nsim = 250, font_size = 10)
## @knitr END

pdf('/Users/Shai/Dropbox (BGU)/Apps/Overleaf/CRISPR-Networks-NEE/figures/SI_host_spacer_empirical.pdf',12,8)
plot_grid(Yellowstone$p1, Yellowstone$p2, 
          Pseudomonas$p1, Pseudomonas$p2,
          Russia2010$p1, Russia2010$p2, 
          ncol=2, nrow=3, align = 'vh', labels = letters[1:6], label_size = 16, scale = 0.95)
dev.off()

pdf('/Users/Shai/Dropbox (BGU)/Apps/Overleaf/CRISPR-Networks-NEE/figures/SI_immunity_nestedness_signif.pdf',12,8)
plot_grid(Yellowstone$p4, NULL, 
          Pseudomonas$p4, Pseudomonas$p3,
          Russia2010$p4, Russia2010$p3, 
          ncol=2, nrow=3, align = 'vh', labels = c('a','',letters[2:5]), label_size = 16, scale = 0.95)
dev.off()

pdf('/Users/Shai/Dropbox (BGU)/Apps/Overleaf/CRISPR-Networks-NEE/figures/empirical_data.pdf',12,8)
plot_grid(Yellowstone$p5,Pseudomonas$p5,Russia2010$p5, 
          ncol=3, nrow=1, labels = letters[1:3], label_size = 16, scale = 1)
dev.off()

write_res <- function(x){
  write_lines(x,'empirical_data_results.txt',append = T)
}
write_lines(Sys.time(),'empirical_data_results.txt',append = F)
write_res('YELLOWSTONE')
write_res(Yellowstone$host_sp_modularity$L)
write_res(Yellowstone$host_sp_modularity$L_sim[1:5])
write_res(Yellowstone$host_sp_modularity$pvalue)
write_res(Yellowstone$immunity_nestedness$ev_obs)
write_res(Yellowstone$immunity_nestedness$ev_shuffled[1:5])
write_res(Yellowstone$immunity_nestedness$p_value_ev)
write_res('\n')
write_res('PSEUDOMONAS')
write_res(Pseudomonas$host_sp_modularity$L)
write_res(Pseudomonas$host_sp_modularity$L_sim[1:5])
write_res(Pseudomonas$host_sp_modularity$pvalue)
write_res(Pseudomonas$immunity_nestedness$ev_obs)
write_res(Pseudomonas$immunity_nestedness$ev_shuffled[1:5])
write_res(Pseudomonas$immunity_nestedness$p_value_ev)
write_res(Pseudomonas$immunity_nestedness$WNODF_obs)
write_res(Pseudomonas$immunity_nestedness$WNODF_shuffled[1:5])
write_res(Pseudomonas$immunity_nestedness$p_value_WNODF)
write_res('\n')
write_res('Russia2010')
write_res(Russia2010$host_sp_modularity$L)
write_res(Russia2010$host_sp_modularity$L_sim[1:5])
write_res(Russia2010$host_sp_modularity$pvalue)
write_res(Russia2010$immunity_nestedness$ev_obs)
write_res(Russia2010$immunity_nestedness$ev_shuffled[1:5])
write_res(Russia2010$immunity_nestedness$p_value_ev)
write_res(Russia2010$immunity_nestedness$WNODF_obs)
write_res(Russia2010$immunity_nestedness$WNODF_shuffled[1:5])
write_res(Russia2010$immunity_nestedness$p_value_WNODF)


# Phylogenetic signal for Russia2010 data set -----------------------------------------

## @knitr load_data_phylogenetic_analysis
criterium <- '4mm/-PAM'
dataset_id <- 6
as_tibble(dbGetQuery(db, paste('SELECT name FROM data_sets WHERE id=',dataset_id,sep='')))
data_virus <- as_tibble(dbGetQuery(db, paste('SELECT virus_strain, spacer, criteria FROM virus WHERE dataset_id=',dataset_id,sep='')))
data_virus %<>% filter(criteria==criterium) %>% select(-criteria) %>% mutate(w=1)
data_host <- as_tibble(dbGetQuery(db, paste('SELECT strain_id, spacer FROM bacteria WHERE dataset_id=',dataset_id,sep='')))
data_host %<>% mutate(w=1)

# Create monolayer objects for infomapecology
data_virus <- create_monolayer_object(data_virus, directed = F, bipartite = T, group_names = names(data_virus)[1:2])
data_host <- create_monolayer_object(data_host, directed = F, bipartite = T, group_names = names(data_host)[1:2])
host_sp_modularity <- run_infomap_monolayer(x = data_host, infomap_executable = 'Infomap', flow_model = 'undirected', silent = T, trials = 100, two_level = T, seed = 109743, signif = F)

# Get the tree
tree <- treeio::read.tree('~/GitHub/ecomplab/CRISPR_networks/M2010_uzonroot_gtrgamma_raxml_1000boot.nwk')
D <- ape::cophenetic.phylo(tree) # Phyloegentic distance
D <- matrix_to_list_unipartite(D, directed = T)
D <- D$edge_list
## @knitr END

plot_modular_matrix(host_sp_modularity, fix_coordinates = F)
plot(tree)

## @knitr strains_tree_network
# Difference between tree and matrix
nodes_in_modules <- host_sp_modularity$modules %>%
  filter(str_starts(node_name, 'M')) %>%
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
  filter(str_starts(node_name, 'M')) %>% 
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


# 
# # Try fitDiscrete
# tree$tip.label
# 
# S <- host_sp_modularity$modules %>%
#   filter(str_starts(node_name, 'M')) %>% 
#   select(node_name, module_level1) %>% 
#   drop_na() %>% 
#   dplyr::select(node_name, module_level1) %>% 
#   mutate(node_name=str_replace_all(node_name, pattern = '\\.', '')) %>% 
#   filter(node_name %in% tree$tip.label) %>% 
#   filter(node_name!='M03042')
# traits <- as.character(as.vector(S$module_level1))
# names(traits) <- S$node_name
# traits <- traits[tree$tip.label]
# 
# tree <- drop.tip(tree, setdiff(tree$tip.label, names(traits)))
# plot(tree)
# traits <- as.vector(na.omit(traits))
# names(traits) <- tree$tip.label
# m_obs <- fitDiscrete(phy = tree, dat=traits, model = 'ER')
# obs <- m_obs$opt$aicc
# 
# nperm=100
# res <- obs
# 
# 
# for (i in 1:99){
#   traits_perm <- sample(traits, replace = F)
#   names(traits_perm) <- tree$tip.label
#   m_perm <- fitDiscrete(phy = tree, dat=traits_perm, model = 'ER')
#   res <- c(res, m_perm$opt$aicc)
#   print(m_perm$opt$aicc)
# }
# 
# ltail <- res[nperm*0.025] # lower tail
# htail <- res[nperm*0.975] # higher tail
# which(res==obs)
# qplot(res)+
#   geom_vline(xintercept =  obs)+
#   geom_vline(xintercept =  c(ltail, htail), color='red')
# signif <- obs<ltail | obs>htail # Is significant?
# 
# 
# 
# # Get module sharing
# S <- host_sp_modularity$modules %>%
#   filter(str_starts(node_name, 'M')) %>% 
#   select(node_name, module_level1) %>% 
#   drop_na() %>% 
#   dplyr::select(node_name, module_level1) %>% 
#   mutate(in_module=1, node_name=str_replace_all(node_name, pattern = '\\.', '')) %>% 
#   spread(key = module_level1, value = in_module, fill = 0)
# rnames <- S$node_name
# S <- data.matrix(S[,-1])
# rownames(S) <- rnames
# dim(S)
# S <- tcrossprod(S)
# 
# # Only work with overlapping strains
# overlapping <- intersect(rownames(S), rownames(D))
# D <- D[overlapping, overlapping]
# S <- S[overlapping, overlapping]
# isSymmetric(D)
# isSymmetric(S)
# 
# # Logistic regression
# logreg_dat <- data.frame(x= D[lower.tri(D)], y=S[lower.tri(S)])
# ggplot(logreg_dat, aes(x,y))+geom_point()+
#   geom_smooth(method = "glm",
#               method.args = list(family = "binomial"),
#               se = FALSE)
# m <- glm(y~x, family=binomial(link='logit'), data=logreg_dat)
# 
# summary(m)
# coef(m)
# exp(coef(m))
# 
# # Compare the coefficient of the model to a random expectation.
# nperm <- 999
# coeff_perm <- coef(m)[2] # Add the observed coefficient
# for(p in 1:nperm) {
#   newSample <- sample(logreg_dat$y,size = nrow(logreg_dat), replace = F)
#   newdata <- data.frame(Y=newSample, X=logreg_dat$x)
#   newfit <- glm(Y ~ X, data=newdata, family=binomial(link = "logit"))
#   coeff_perm <- c(coeff_perm, coef(newfit)[2])
# }
# coeff_perm <- sort(coeff_perm)
# ltail <- coeff_perm[nperm*0.025] # lower tail
# htail <- coeff_perm[nperm*0.975] # higher tail
# which(coeff_perm==coef(m)[2])
# qplot(coeff_perm)+
#   geom_vline(xintercept =  coef(m)[2])+
#   geom_vline(xintercept =  c(ltail, htail), color='red')
# signif <- coef(m)[2]<ltail |  coef(m)[2]>htail # Is significant?
# 

# Other analyses ----------------------------------------------------------

# Any other analysis of modularity/nestedness can follow the examples in the
# main function. Example for virus-protospacer network:

# Is virus-protospacer modular?
virus_ps_shuffled <- shuffle_infomap(data_virus, shuff_method = 'r00', nsim = 50, burnin=1000)
virus_ps_modularity <- run_infomap_monolayer(x = data_virus, infomap_executable = 'Infomap', flow_model = 'undirected', silent = T, trials = 100, two_level = T, seed = 109743, signif = T, shuff_method = virus_ps_shuffled, nsim = NULL)
virus_ps_modularity$call
virus_ps_modularity$pvalue
tibble(L_sim=virus_ps_modularity$L_sim) %>% ggplot(aes(L_sim))+geom_histogram()+geom_vline(xintercept = virus_ps_modularity$L)
class(virus_ps_modularity)
plot_modular_matrix(x, fix_coordinates=T, axes_titles=c('Virus', 'Protospacer'), transpose = T)

# Is it nested?
virus_ps_nestedness <- weighted_nestedness_significance(data_virus$mat, nsim = 100, weighted = F, shuff_method = 'r00')
virus_ps_nestedness[1:6]
virus_ps_nestedness$plot_WNODF
virus_ps_nestedness$plot_ev

# Obtain infection / immunity networks ------------------------------------
matches_df <- get_matches(host_data = data_host$edge_list, virus_data = data_virus$edge_list, method = 'dataframe')
# Create infection edge list and matrices
infection_df <- as_tibble(matches_df$immunity_network %>% 
                            filter(matches==0) %>% 
                            rename(w=matches) %>% 
                            select(strain_id, virus_strain, w) %>% 
                            mutate_at(c(1,2), as.character) %>% 
                            mutate(w=1))

immunity_network <- create_monolayer_object(matches_df$immunity_network %>% filter(matches!=0), directed = F, bipartite = T, group_names = c('host','virus'))
infection_network <- create_monolayer_object(infection_df, directed = F, bipartite = T, group_names = c('host','virus'))

# Hosts that appear both in infection and immunity networks
intersect(rownames(immunity_network$mat),rownames(infection_network$mat))
# Viruses that appear both in infection and immunity networks
intersect(colnames(immunity_network$mat),colnames(infection_network$mat))