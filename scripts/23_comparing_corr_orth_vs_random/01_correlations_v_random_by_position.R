library(ggplot2)
library(reshape)

correlation_table <- read.delim("output_files/cor_of_DE_vs_random_NAs_for_hist.txt")
correlation_table$am_e_real=as.numeric(levels(correlation_table$am_e_real))[correlation_table$am_e_real]
correlation_table$am_l_real=as.numeric(levels(correlation_table$am_l_real))[correlation_table$am_l_real]
correlation_table$gg_e_real=as.numeric(levels(correlation_table$gg_e_real))[correlation_table$gg_e_real]
correlation_table$gg_l_real=as.numeric(levels(correlation_table$gg_l_real))[correlation_table$gg_l_real]
correlation_table$mm_real=as.numeric(levels(correlation_table$mm_real))[correlation_table$mm_real]


#test <- correlation_table[1:10,]
#test <- correlation_table[11:20,]
#test <- correlation_table[21:30,]
test <- correlation_table[31:40,]
correlation_table <- test

real_dist <- subset(correlation_table , gene_list_type=="real")
real_dist <- rbind(real_dist$mm_real, real_dist$am_e_real, real_dist$am_l_real, real_dist$gg_e_real, real_dist$gg_l_real)
#hist(real_dist, breaks=20, xlim=c(-1,1))

random_dist <- subset(correlation_table , gene_list_type=="random")
random_dist <- rbind(random_dist$mm_real, random_dist$am_e_real, random_dist$am_l_real, random_dist$gg_e_real, random_dist$gg_l_real)
#hist(random_dist, breaks=20, xlim=c(-1,1))

# combine these to be a single histogram
real_dist <- as.vector(real_dist)
random_dist <- as.vector(random_dist)

dataforplotting <- rbind(data.frame("data_type"=rep("random", length(random_dist)), "corval"= random_dist),
                         data.frame("data_type"=rep("real", length(real_dist)), "corval"= real_dist))
test <- na.omit(dataforplotting)

#ggplot(test, aes(corval, fill = data_type)) + geom_density(alpha = 0.4)
#ggplot(test, aes(corval, fill = data_type)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity', binwidth = 0.075)
ggplot(test, aes(corval, fill = data_type)) + 
  geom_histogram(alpha = 0.5, position = 'identity', binwidth = 0.075) +
  theme_bw() + 
  geom_vline(xintercept=mean(real_dist, na.rm = T), color="red") +
  geom_vline(xintercept=mean(random_dist, na.rm = T), color="red") +
  xlim(-1,1)



t.test(real_dist,random_dist)
1.96*sd(random_dist, na.rm = T)


wilcox.test(real_dist,random_dist, alternative = "two.sided")
