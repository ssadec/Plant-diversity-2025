rm(list = ls())

data_fig2e<-read.csv('/Users/yusha/Fig.2e.csv',header = TRUE,sep = ",")
View(data_fig2e)

data_fig2e$V2<-as.numeric(data_fig2e$V2)
data_fig2e$V3<-as.numeric(data_fig2e$V3)
data_fig2e$V4<-as.numeric(data_fig2e$V4)
data_fig2e$V5<-as.factor(data_fig2e$V5)

b<-unlist(lapply(data_fig2e,is.numeric))
b


fig2e <- ggplot(data = data_fig2e, family = font_family, aes(x = V2, y = as.factor(V1))) +  
  geom_point(aes(x = V2, y = as.factor(V1), color = V5), size = 3.5) +  
  geom_errorbar(aes(x = V2, y = as.factor(V1), xmin = V3, xmax = V4, color = V5), linewidth = 1, width = 0.2) +  
  scale_y_discrete(limits = c("Total crop response (0/0)","Total herbivore response (53/5)","Total NEs response (0/0)","Total parasitoid response (12/2)","Total predator response (15/3)",
                              "","Total crop response (1/1)","Total herbivore response (14/4)","Total NEs response (0/0)","Total parasitoid response (0/0)","Total predator response (1/1)",
                              "","Total crop response (391/104)","Total herbivore response (617/124)","Total NEs response (19/7)","Total parasitoid response (93/31)","Total predator response (304/73)")) +  
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

fig2e<- fig2e + coord_fixed(ratio = 0.415)
fig2e


ggsave("fig2e.png", fig2, width = 18, height = 18)

