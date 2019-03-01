lightgrey <- brewer.pal(9, "Greys")[2]
plottheme <- theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(face = "bold", size = 12),
        panel.grid.major.y = element_line(colour = lightgrey),
        panel.grid.major.x = element_line(colour = lightgrey),
        panel.border = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 13, face = "bold"),
        axis.title.x = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"))