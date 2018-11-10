 ggplot()+
 geom_segment(aes(x = 1, y = 1, xend = 1000, yend = 1, colour = "chromosome 1")) +
 geom_segment(aes(x = 200, y = 0.75, xend = 200, yend = 1.25, colour = "Lethal")) +
 geom_segment(aes(x = 300, y = 0.75, xend = 300, yend = 1.25, colour = "Viable")) +
 scale_color_manual(values=c("black","green"))

 geom_segment(aes(x = 1, y = 1, xend = 1000, yend = 1, colour = "chromosome 2")) +
 geom_segment(aes(x = 200, y = 0.75, xend = 200, yend = 1.25, colour = "Lethal")) +
 geom_segment(aes(x = 300, y = 0.75, xend = 300, yend = 1.25, colour = "Viable")) +


