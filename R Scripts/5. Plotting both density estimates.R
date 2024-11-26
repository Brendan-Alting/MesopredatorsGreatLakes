library(ggplot2)

#Density estimate plot results

Density <- c(0.12,0.38,0.13,0.28,0.30,1.58)
LCI <- c(0.07,0.28,0.07,0.18,0.2,1.09)
UCI <- c(0.18,0.51,0.69,0.73,0.67,2.32)

Species <- c("Dingo","Quoll","Dingo","Quoll","Fox","Lace Monitor")
ModelType <- c("SECR","SECR","SUN","SUN","SUN","SUN")

Densityestimates <- data.frame(cbind(Density,LCI,UCI,Species,ModelType))

Densityestimates$Density <- as.numeric(Densityestimates$Density)
Densityestimates$LCI <- as.numeric(Densityestimates$LCI)
Densityestimates$UCI <- as.numeric(Densityestimates$UCI)
Densityestimates$Species <- factor(Densityestimates$Species, levels = c("Dingo", "Quoll","Fox","Lace Monitor"))
Densityestimates$ModelType <- factor(Densityestimates$ModelType, levels = c("SECR","SUN"))


densityplotests <- ggplot(Densityestimates, aes(x = Species, y = Density, color = ModelType)) +
  geom_point(position = position_dodge(width = 0.5),size = 3)+ 
  geom_errorbar(aes(ymin = LCI, ymax = UCI),position = position_dodge(width = 0.5), width = 0.2) +
  scale_color_manual(values = c("SECR" = "blue", "SUN" = "red")) +  # Set colors
  theme_minimal() +  # Use a minimal theme
  labs(x = "Species", y = expression(Density~(Individuals/km^2)), color = "Model") +  # Labels
  theme(axis.title = element_text(size = 20),       
              axis.title.y = element_text(vjust = 0.5),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 16),
              axis.text.x = element_text(size = 16,angle=45,hjust=1),
              axis.text.y = element_text(size = 16),
              panel.background = element_blank(),      # Remove background panel
              axis.ticks.x = element_line(size = 0.5),
              axis.ticks.y = element_line(size = 0.5),
              axis.line.x = element_line(color = "black", size = 0.5), # Add x-axis line
              axis.line.y = element_line(color = "black", size = 0.5)) # Add y-axis line))


#plot

##
png("Figures/Density Estimates.png", width = 10, height = 7, res= 300, units = "in")

densityplotests
dev.off()

