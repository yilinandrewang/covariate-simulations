setwd("/Users/andrewang/Documents/iLibrary/UC DAVIS/My Research/Covariate Simulation/Results and Outputs")
full <- read.csv("Graphing Data.csv", TRUE, ",")
View(full)
full[,6] <- as.factor(full[,6])
str(full)
library(ggplot2);library(grid);library(gridExtra)

#### rho = 0.75 ####
p1.0 <- ggplot(full[1:18,1:6], aes(x=n, color=effect.size))
p1.1 <- p1.0 + geom_ribbon(aes(ymin=s1, ymax=s2, fill=effect.size), alpha=0.3, colour=NA)
p1.2 <- p1.1 + geom_line(aes(x=n, y=s1), size = 0.5)
p1.3 <- p1.2 + geom_hline(yintercept = c(0.8,0.1)) + ggtitle("Power Boost vs. Type I Error Inflation at rho(y,c) = 0.75")
p1.4 <- p1.3 + scale_x_continuous(breaks = c(20,50,100), labels = c("n = 20", "n = 50", "n = 100"))
p1.5 <- p1.4 + scale_y_continuous(breaks = c(0,0.05,0.10,0.20,0.40,0.60,0.80,1),
                                  labels = c("0","0.05","0.10","0.20","0.40","0.60","0.80","1"))
p1.6 <- p1.5 + theme(axis.text=element_text(size=9,color="black"),
             axis.title=element_blank(),
             title=element_text(size=10,color="black",vjust=1),
             legend.text=element_text(size=8,color="black"),
             plot.margin=grid::unit(c(1, 0, 1, 1), "cm"))
p1.6
#### rho = 0.5 ####
p2.0 <- ggplot(full[19:36,1:6], aes(x=n, color=effect.size))
p2.1 <- p2.0 + geom_ribbon(aes(ymin=s1, ymax=s2, fill=effect.size), alpha=0.3, colour=NA)
p2.2 <- p2.1 + geom_line(aes(x=n, y=s1), size = 0.5)
p2.3 <- p2.2 + geom_hline(yintercept = c(0.8,0.1)) + ggtitle("Power Boost vs. Type I Error Inflation at rho(y,c) = 0.50")
p2.4 <- p2.3 + scale_x_continuous(breaks = c(20,50,100), labels = c("n = 20", "n = 50", "n = 100"))
p2.5 <- p2.4 + scale_y_continuous(breaks = c(0,0.05,0.10,0.20,0.40,0.60,0.80,1),
                                  labels = c("0","0.05","0.10","0.20","0.40","0.60","0.80","1"))
p2.6 <- p2.5 + theme(axis.text=element_text(size=9,color="black"),
                     axis.title=element_blank(),
                     title=element_text(size=10,color="black",vjust=1),
                     legend.text=element_text(size=8,color="black"),
                     plot.margin=grid::unit(c(1, 0, 1, 1), "cm"))

#### rho = 0.25 ####
p3.0 <- ggplot(full[37:54,1:6], aes(x=n, color=effect.size))
p3.1 <- p3.0 + geom_ribbon(aes(ymin=s1, ymax=s2, fill=effect.size), alpha=0.3, colour=NA)
p3.2 <- p3.1 + geom_line(aes(x=n, y=s1), size = 0.5)
p3.3 <- p3.2 + geom_hline(yintercept = c(0.8,0.1)) + ggtitle("Power Boost vs. Type I Error Inflation at rho(y,c) = 0.25")
p3.4 <- p3.3 + scale_x_continuous(breaks = c(20,50,100), labels = c("n = 20", "n = 50", "n = 100"))
p3.5 <- p3.4 + scale_y_continuous(breaks = c(0,0.05,0.10,0.20,0.40,0.60,0.80,1),
                                  labels = c("0","0.05","0.10","0.20","0.40","0.60","0.80","1"))
p3.6 <- p3.5 + theme(axis.text=element_text(size=9,color="black"),
                     axis.title=element_blank(),
                     title=element_text(size=10,color="black",vjust=1),
                     legend.text=element_text(size=8,color="black"),
                     plot.margin=grid::unit(c(1, 0, 1, 1), "cm"))

#### rho = 0 ####
p4.0 <- ggplot(full[55:72,1:6], aes(x=n, color=effect.size))
p4.1 <- p4.0 + geom_ribbon(aes(ymin=s1, ymax=s2, fill=effect.size), alpha=0.3, colour=NA)
p4.2 <- p4.1 + geom_line(aes(x=n, y=s1), size = 0.5)
p4.3 <- p4.2 + geom_hline(yintercept = c(0.8,0.1)) + ggtitle("Power Boost vs. Type I Error Inflation at rho(y,c) = 0")
p4.4 <- p4.3 + scale_x_continuous(breaks = c(20,50,100), labels = c("n = 20", "n = 50", "n = 100"))
p4.5 <- p4.4 + scale_y_continuous(breaks = c(0,0.05,0.10,0.20,0.40,0.60,0.80,1),
                                  labels = c("0","0.05","0.10","0.20","0.40","0.60","0.80","1"))
p4.6 <- p4.5 + theme(axis.text=element_text(size=9,color="black"),
                     axis.title=element_blank(),
                     title=element_text(size=10,color="black",vjust=1),
                     legend.text=element_text(size=8,color="black"),
                     plot.margin=grid::unit(c(1, 0, 1, 1), "cm"))

#### Individual Graphs ####
p1.6
p2.6
p3.6
p4.6

#### Four Graphs ####
grid.arrange(p4.6, p3.6, p2.6, p1.6, ncol = 2, 
             main=textGrob("Strategy 1 (y ~ x) vs. Strategy 2 (y ~ x OR y ~ x + c) over 10000 Simulations",
                           gp=gpar(fontsize=16, fontface="bold", col="navy")))