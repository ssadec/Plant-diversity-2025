library(ggplot2)
####extended fig2c####

data_fig2c <- read.csv('/Users/yusha/Fig.2c.csv',
                          header = TRUE, sep = ",")
View(data_fig2c)
# 转换列类型
data_fig2c$V2 <- as.numeric(data_fig2c$V2)
data_fig2c$V3 <- as.numeric(data_fig2c$V3)
data_fig2c$V4 <- as.numeric(data_fig2c$V4)
data_fig2c$V5 <- as.factor(data_fig2c$V5)


data_fig2c$V1_unique <- make.unique(as.character(data_fig2c$V1))


y_order <- rev(data_fig2c$V1_unique)
y_labels <- rev(data_fig2c$V1)

# 绘图
fig2c <- ggplot(data = data_fig2c, aes(x = V2, y = V1_unique, color = V5)) +  
  geom_point(size = 4) +  
  geom_errorbar(aes(xmin = V3, xmax = V4), linewidth = 1, width = 0.4) +  
  scale_y_discrete(limits = y_order, labels = y_labels) +  
  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1)) +  
  xlab("") +  
  ylab("") +  
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +  
  theme_bw() +  
  theme(
    legend.position = "none",  
    plot.margin = unit(c(1, 1, 1, 0), 'cm'),  
    panel.grid = element_blank(),  
    panel.border = element_rect(colour = "black", fill = NA, size = 3),  
    panel.spacing = unit(0.1, "lines"),  
    strip.background = element_blank(),  
    strip.text = element_text(size = 48),  
    axis.title.x = element_text(hjust = 0.5),  
    text = element_text(color = "black", family = "Arial", size = 24),  
    axis.text.x = element_text(color = "black"),  
    axis.text.y = element_text(color = "black",
                               margin = margin(r = 1, unit = "pt"),
                               lineheight = 50, size = 18),  
    axis.ticks = element_line(color = "black")
  ) +  
  scale_color_manual(values = c("#BC3A24", "#EFB882", "#06798F", "#7A9A01", "#384D73")) +
  coord_fixed(ratio = 0.295)


fig2c



ggsave("extended_fig2c.png", width = 18, height = 18, units = "cm", dpi = 300)
