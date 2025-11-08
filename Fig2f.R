rm(list = ls())

data_fig2f<-read.csv('/Users/yusha/Fig.2f.csv',header = TRUE,sep = ",")
View(data_fig2f)

data_fig2f$V2<-as.numeric(data_fig2f$V2)
data_fig2f$V3<-as.numeric(data_fig2f$V3)
data_fig2f$V4<-as.numeric(data_fig2f$V4)
data_fig2f$V5<-as.factor(data_fig2f$V5)

b<-unlist(lapply(data_fig2f,is.numeric))
b


fig2f <- ggplot(data = data_fig2f, family = font_family, aes(x = V2, y = as.factor(V1))) +  
  geom_point(aes(x = V2, y = as.factor(V1), color = V5), size = 3.5) +  
  geom_errorbar(aes(x = V2, y = as.factor(V1), xmin = V3, xmax = V4, color = V5), linewidth = 1, width = 0.2) +  
  scale_y_discrete(limits = c("Total crop response (151/35)","Total herbivore response (256/36)","Total NEs response (0/0)","Total parasitoid response (15/6)","Total predator response (15/7)",
                              "","Total crop response (239/68)","Total herbivore response (355/89)","Total NEs response (19/7)","Total parasitoid response (78/25)","Total predator response (290/65)")) +  
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

fig2f<- fig2f + coord_fixed(ratio = 0.61)
fig2f


ggsave("fig2f.png", fig2f, width = 18, height = 18)

