rm(list = ls())

data_fig2a<-read.csv('/Users/yusha/Fig.2a.csv',header = TRUE,sep = ",")

data_fig2a$V2<-as.numeric(data_fig2a$V2)
data_fig2a$V3<-as.numeric(data_fig2a$V3)
data_fig2a$V4<-as.numeric(data_fig2a$V4)
data_fig2a$V5<-as.factor(data_fig2a$V5)

b<-unlist(lapply(data_fig2a,is.numeric))
b


fig2a <- ggplot(data = data_fig2a, family = font_family, aes(x = V2, y = as.factor(V1))) +  
  geom_point(aes(x = V2, y = as.factor(V1), color = V5), size = 3) +  
  geom_errorbar(aes(x = V2, y = as.factor(V1), xmin = V3, xmax = V4, color = V5), linewidth = 1, width = 0.4) +  
  scale_y_discrete(limits = c("Crop quality (11/3)","Crop reproduction (281/80)","Crop growth (100/32)",
                              "","Herbivore damage (239/48)","Herbivore reproduction (445/109)","Herbivore growth (1/1)",
                              "","NEs diversity (6/3)","NEs reproduction (13/5)",
                              "","Parasitoid parasitism (36/18)","Parasitoid diversity (14/2)","Parasitoid reproduction (52/21)","Parasitoid growth (3/1)",
                              "","Predator predation (39/12)","Predator diversity (62/16)","Predator reproduction (221/67)",
                              "","Total crop response (392/105)","Total herbivore response (686/132)","Total NEs response (19/7)","Total parasitoid response (105/32)","Total predator response (322/77)")) +  
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
#scale_color_manual(values = c("#A8D08D","#F1F1F1","#ADD8E6","#FFEB8D")) 
coord_cartesian(clip = 'off')
#annotate("text", x = c(0.3786, -0.0719, 0.8212, 0.6097), y = c(1, 2, 3, 4), 
# label = c("**", "***", "***", "**"), fontface = "bold", size = 6) 

fig2a<- fig2a + coord_fixed(ratio = 0.265)
fig2a


ggsave("fig2a.png", fig1b, width = 18, height = 18)

