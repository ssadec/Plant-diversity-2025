rm(list = ls())

data_figcrop<-read.csv('/Users/yusha/Fig.2b.csv',header = TRUE,sep = ",")

data_figcrop$V2<-as.numeric(data_figcrop$V2)
data_figcrop$V3<-as.numeric(data_figcrop$V3)
data_figcrop$V4<-as.numeric(data_figcrop$V4)
data_figcrop$V5<-as.factor(data_figcrop$V5)

b<-unlist(lapply(data_figcrop,is.numeric))
b


figcrop<- ggplot(data = data_figcrop, family = font_family, aes(x = V2, y = as.factor(V1))) +  
  geom_point(aes(x = V2, y = as.factor(V1), color = V5), size = 3.5) +  
  geom_errorbar(aes(x = V2, y = as.factor(V1), xmin = V3, xmax = V4, color = V5), linewidth = 1, width = 0.3) +  
  scale_y_discrete(limits = c("Total crop response (82/26)","Total herbivore response (170/137)","Total NEs response (9/3)","Total parasitoid response (34/8)","Total predator response (108/23)",
                              "","Total crop response (310/81)","Total herbivore response (516/97)","Total NEs response (10/5)","Total parasitoid response (71/25)","Total predator response (214/55)",
                              "","Total crop response (25/6)","Total herbivore response (115/25)","Total NEs response (8/2)","Total parasitoid response (29/7)","Total predator response (111/19)",
                              "","Total crop response (367/99)","Total herbivore response (571/108)","Total NEs response (11/5)","Total parasitoid response (76/25)","Total predator response (211/58)"
                             )) +  
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

figcrop<- figcrop + coord_fixed(ratio = 0.3)
figcrop


ggsave("figcrop.png", figcrop, width = 18, height = 18)

