library(data.table)
library(ggplot2)
library(dplyr)
library(hrbrthemes)

df_sml <- fread("df_sml.csv") %>% as.data.frame()
head(df_sml[,c(1:5,893:897)])

scores_avg <- fread('Ensemb_ScoresAveraged.csv') %>% as.data.frame()
scores_bin <-  fread('Ensemb_ScoresPerCell.csv') %>% as.data.frame()
scores_avg
scores_bin

boxplot(scores_bin)
boxplot(scores_avg)

scores_eln <- fread('elnet_f1_scores.csv') %>% 
  as.data.frame() %>%
  select(-1)
colnames(scores_eln) <- paste0("eln_",colnames(scores_eln))

boxplot(scores_eln)

colnameses <- c()
for(col in 1:5){
assign(paste0('w_',col), wilcox.test(scores_bin[,col+15] - scores_eln[,col], exact = FALSE))
  assign(paste0('bind_',col), cbind(scores_bin[,col+15],scores_eln[,col]))
  colnameses <- c(colnameses, colnames(scores_bin)[col+15], colnames(scores_eln)[col])
}

total_mat <- cbind(bind_1, bind_2, bind_3, bind_4, bind_5)
colnames(total_mat) <- colnameses
total_mat
boxplot(total_mat)

nameses <- gsub("f1_","",melt(total_mat)[,2])
nameses <- gsub("eln_","",nameses)
nameses <- gsub("[_]", " ", nameses)
nameses2 <- paste0(nameses,c(rep(", Ensemble",10),rep(", Elastic_Net",10)))

plot_mat <- as.data.frame(nameses2)
plot_mat$values <-  melt(total_mat)[,3]
plot_mat$group <- nameses
plot_mat$ticks <- rep(c(rep("XGB-EN",10),rep("EN",10)),5)

plot_mat

plot_mat %>%
  ggplot(aes(x=as.factor(nameses2), y = values, fill=nameses2)) + 
  geom_boxplot() +
  theme_ipsum() +
  theme(legend.position="none",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.text=element_text(size=9),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c(rep(c("EN","XGB-EN"),5))) +
  facet_wrap(.~group, ncol=5, scales = 'free_x') +
  ylab('F1 Scores') + xlab("") + 
  ggtitle('F1 Scores By Model and Cell Type') +
  scale_fill_brewer(name = 'Legend\nCell_Type.Model', palette = "Paired") -> cell_boxes

cell_boxes 

ggsave(
  "box_per_cell.png",
  plot = cell_boxes,
  device = 'png',
  scale = 1,
  width = 12.5,
  height = 6.5,
  units = 'in',
  dpi = 300
)


statistic_value <- c(w_1$statistic, w_2$statistic, w_3$statistic, w_4$statistic, w_5$statistic)
p_values <- c(w_1$p.value, w_2$p.value, w_3$p.value, w_4$p.value, w_5$p.value)

statistic_value
p_values

mean_difference <- c(mean(scores_bin[,16] - scores_eln[,1]),
  mean(scores_bin[,17] - scores_eln[,2]),
  mean(scores_bin[,18] - scores_eln[,3]),
  mean(scores_bin[,19] - scores_eln[,4]),
  mean(scores_bin[,20] - scores_eln[,5]))

mean_values <- c(mean(mean_difference), mean(statistic_value), mean(p_values))

tab1 <- rbind(mean_difference, statistic_value, p_values)
tab2 <- cbind(tab1, mean_values) %>% as.data.frame()
colnames(tab2) <- c(unique(nameses),"averages")
fwrite(tab2, paste0("p_vals_table.csv"))



wilcox.test(scores_avg[,4] - apply(scores_eln,1,mean),exact=FALSE)


cbind(apply(scores_eln,1,mean),scores_avg[,4]) %>%
  as.data.frame %>% 
  melt() %>%
  ggplot(aes(x=variable, y = value, fill=variable)) + 
  geom_boxplot() +
  theme_ipsum() +
  theme(legend.position="none",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.text=element_text(size=9),
        plot.title = element_text(hjust = 0.5)) +
  ylab('F1 Scores') + xlab("") + 
  scale_x_discrete(labels=c("EN","XGB-EN")) +
  ggtitle('Average F1 Scores Across Cell Types') +
  scale_fill_brewer(name = 'Legend\nCell_Type.Model', palette = "Paired") -> cell_boxes_2

cell_boxes_2

ggsave(
  "box_avg.png",
  plot = cell_boxes_2,
  device = 'png',
  scale = 1,
  width = 8,
  height = 6,
  units = 'in',
  dpi = 300
)
