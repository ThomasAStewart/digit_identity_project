# initialize packages ####
library(ggplot2); library(plyr)

# load data ####
correlation_table <- read.table("output_files/corr_table_orth_vs_random_pentadactyl_copy3.txt", header = T) # both stages, but remove self-comparisons
#correlation_table <- read.table("output_files/corr_table_orth_vs_random_pentadactyl.txt", header = T) # both alligator early and late
#correlation_table <- read.table("output_files/corr_table_orth_vs_random_pentadactyl_copy.txt", header = T) # alligator early only
#correlation_table <- read.table("output_files/corr_table_orth_vs_random_pentadactyl_copy2.txt", header = T) # alligator late only

#correlation_table[correlation_table==1.00000000] <- NA # exclude the correlations of one vector to itself (those with corr = 1.0)

position <- c("d1_d2", "d2_d3", "d3_d4", "d4_d5")
t_test_stats_for_figure <- vector("list",length(position))
u_test_stats_for_figure <- vector("list",length(position))

for(k in 1:length(position)){
  
  corr_at_one_position <- correlation_table[correlation_table$position==position[k],]
  real_dist <-subset(corr_at_one_position, comparison=="orthologs")
  random_dist <- subset(corr_at_one_position , comparison=="random")
  
  real_dist <- as.vector(rbind(real_dist[,4], real_dist[,5], real_dist[,6], real_dist[,7]))
  random_dist <- as.vector(rbind(random_dist[,4], random_dist[,5], random_dist[,6], random_dist[,7]))
  
  dataforplotting <- rbind(data.frame("data_type"=rep("random", length(random_dist)), "corval"= random_dist),
                           data.frame("data_type"=rep("orthologs", length(real_dist)), "corval"= real_dist))
  dataforplotting <- na.omit(dataforplotting)
  
  q <- ggplot(dataforplotting, aes(corval, fill = data_type)) + 
    geom_histogram(alpha = 0.5, position = 'identity', binwidth = 0.075) +
    theme_bw() + 
    scale_colour_grey() +
    ggtitle(position[k]) +
    geom_vline(xintercept=mean(real_dist, na.rm = T), color="red") +
    geom_vline(xintercept=mean(random_dist, na.rm = T), color="red") +
    xlim(-1,1)
  plot(q)
  
  t_test_stats_for_figure[[k]] <- t.test(real_dist,random_dist)
  u_test_stats_for_figure[[k]] <-   wilcox.test(real_dist,random_dist, alternative = "two.sided")
}
