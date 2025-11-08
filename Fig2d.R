rm(list = ls())

data_fig2d<-read.csv('/Users/yusha/Fig.2d.csv',header = TRUE,sep = ",")
View(data_fig2d)

data_fig2d$V2<-as.numeric(data_fig2d$V2)
data_fig2d$V3<-as.numeric(data_fig2d$V3)
data_fig2d$V4<-as.numeric(data_fig2d$V4)
data_fig2d$V5<-as.factor(data_fig2d$V5)

b<-unlist(lapply(data_fig2d,is.numeric))
b


fig2d <- ggplot(data = data_fig2d, family = font_family, aes(x = V2, y = as.factor(V1))) +  
  geom_point(aes(x = V2, y = as.factor(V1), color = V5), size = 3.5) +  
  geom_errorbar(aes(x = V2, y = as.factor(V1), xmin = V3, xmax = V4, color = V5), linewidth = 1, width = 0.3) +  
  scale_y_discrete(limits = c("Total crop response (84/22)","Total herbivore response (84/25)","Total NEs response (9/3)","Total parasitoid response (35/11)","Total predator response (76/17)",
                              "","Total crop response (76/16)","Total herbivore response (148/32)","Total NEs response (0/0)","Total parasitoid response (14/1)","Total predator response (121/22)",
                              "","Total crop response (232/69)","Total herbivore response (454/78)","Total NEs response (10/4)","Total parasitoid response (56/21)","Total predator response (125/40)")) +  
  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1)) +  
  xlab("") +  
  ylab("") +  
  #geom_hline(yintercept = 4, color="red", linetype="dashed", size=1) + 
  geom_vline(aes(xintercept = 0), color = "black", linetype = "dashed", linewidth = 2) +  
  theme_bw() +  
  theme(legend.position = "none",  
        plot.margin = unit(c(1, 1, 1, 0), 'cm'),  
        panel.grid = element_blank(),  
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 3),
        panel.spacing = unit(0.1, "lines"), 
        strip.background = element_blank(),  
        strip.text = element_text(size = 48),  
        axis.title.x = element_text(hjust = 0.5),  
        text = element_text(color = "black", family = "Arial", size = 24),  
        axis.text.x = element_text(color = "black"),  
        axis.text.y = element_text(color = "black",margin = margin(r = 1, unit = "pt"), 
                                   lineheight = 50, size = 18),  
        axis.ticks = element_line(color = "black")) +  
  scale_color_manual(values = c("#BC3A24","#EFB882","#06798F","#7A9A01","#384D73")) 

coord_cartesian(clip = 'off')
#annotate("text", x = c(0.3786, -0.0719, 0.8212, 0.6097), y = c(1, 2, 3, 4), 
# label = c("**", "***", "***", "**"), fontface = "bold", size = 6) 

fig2d<- fig2d + coord_fixed(ratio = 0.4)
fig2d


ggsave("fig2d.png", fig2d, width = 18, height = 18)

