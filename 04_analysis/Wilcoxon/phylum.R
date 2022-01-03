
rm(list = ls())
library(tidyverse)
library(stringr)
data <- read.table("./data/table_phylum.tsv", header = TRUE, row.names = 1, sep = "\t")
group <- read.table("./data/metadata.txt", header = TRUE, sep = "\t")

# 精简物种注释
rownames(data)<-sapply(str_split(rownames(data), "__"),  "[", 3)
# data <- data[, -ncol(data)]
dim(data)
# 把绝对丰度变为相对丰度
# colSums(data / colSums(data)) # 不能直接这样求相对丰度，见https://www.jianshu.com/p/26a926ad1b76
data<-scale(data, center = FALSE, scale = colSums(data))
colSums(data) # 看来相对丰度就是每个ASV在每个样本的比例？
data <- data* 100
# 过滤了平均丰度低于1%的功能分类。
# data <- data %>% filter(apply(data, 1, mean) > 1)

data <- t(data)
data1 <- data.frame(data, group$Group)
colnames(data1) <- c(colnames(data), "Group")
data1$Group <- as.factor(data1$Group)

diff <- data1 %>%
    select_if(is.numeric) %>%
    map_df(~ broom::tidy(wilcox.test(. ~ Group, data = data1)), .id = "var")


diff$p.value <- p.adjust(diff$p.value, "bonferroni")
# 查看p值小于0.05
diff <- diff %>% filter(p.value < 0.05)
# Actinobacteria

abun.bar <- data1[, c(diff$var, "Group")] %>%
    gather(variable, value, -Group) %>%
    group_by(variable, Group) %>%
    summarise(Mean = mean(value),Std= sd(value))
abun.bar$variable <- paste0("*", abun.bar$variable)


# 只提取AD列的
AD_data <- abun.bar[abun.bar$Group == "AD", ]
# 根据AD的菌种进行排序
AD_order<-AD_data[order(AD_data$Mean, decreasing = TRUE), ]$variable
# 给因子排序
abun.bar$variable <- factor(abun.bar$variable, levels = AD_order)
library(ggplot2)
cbbPalette <- c("#f8756b", "#00bec6")

p <- ggplot(abun.bar, aes(variable, Mean, fill = Group)) +
    coord_flip() +
    xlab("") +
    ylab("Mean abundance (%)") +
    theme_set(theme_bw()) +
    scale_y_continuous(expand = c(0, 0))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_bar(stat = "identity",position = "dodge",  colour = "black") +
    geom_errorbar(aes(ymin=Mean-0.1*Std, ymax=Mean+0.1*Std), size=.2,position=position_dodge(.9),width=0.2) +
    scale_fill_manual(values=cbbPalette)

p
ggsave("plot/phylum.svg", plot=p,width = 20, height = 4, units = "cm")
