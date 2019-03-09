lightgrey <- brewer.pal(9, "Greys")[2]
plottheme <- theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(face = "bold", size = 8),
        panel.grid.major.y = element_line(colour = lightgrey),
        panel.grid.major.x = element_line(colour = lightgrey),
        panel.border = element_blank(),
        panel.spacing = unit(0.6, "lines"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
        legend.position = "top",
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.6, "line"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 7, face = "bold"),
        axis.text.x = element_text(size = 7, face = "bold"))