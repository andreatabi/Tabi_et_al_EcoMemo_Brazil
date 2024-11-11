rm(list=ls()) 
 
source("~/Toolbox.R")     
source("~/IC.R")     
library(pcalg)
library(igraph)
library(dplyr)
library(ppcor)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(corrplot)
library(psychometric)
library(cowplot)
library(grid)
library(gridExtra)
library(colorspace)

load("~/Dropbox/Macro_BR/docs/df.RData")

####--------------------------------------------------------------------------------
####--------------------------------------------------------------------------------
#### FIGURE 2 
####--------------------------------------------------------------------------------
####--------------------------------------------------------------------------------

pp.cols <- c("#98E3D5", "#47D7AC", "#41b6e6", "#B5D7DD", 
             "#8986CA", "#D9BFDC", "#F57EB6", "#F8C4CC", 
             "#EDDD87")

df.fig2a1 <- dd %>% dplyr::select(Year, ID, Shannon_entropy) %>% distinct() %>% 
              group_by(Year) %>% summarise(m=mean(Shannon_entropy), s=sd(Shannon_entropy)/sqrt(length(Year)) )

Fig2A.1 <- ggplot(df.fig2a1) + 
  geom_errorbar(aes(x=Year, ymin = m - 1.96*s, ymax = m + 1.96*s), color="#B5D7DD", width=0)+
  geom_point(aes(x=Year, y=m), fill="#98E3D5", pch=21)+
  theme_bw()+ylab("LCS")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position="none", axis.title = element_text(size=14), 
        axis.text = element_text(size=14), strip.text = element_text(size = 14))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black"))+
  scale_color_distiller(palette="Pastel1")
Fig2A.1 

df.fig2a2 <- dd %>% dplyr::select(Year, ID, Wasserstein_entropy) %>% distinct() %>% 
  group_by(Year) %>% summarise(m=mean(Wasserstein_entropy), s=sd(Wasserstein_entropy)/sqrt(length(Year)) )

Fig2A.2 <- ggplot(df.fig2a2) + 
  geom_errorbar(aes(x=Year, ymin = m - 1.96*s, ymax = m + 1.96*s), color="#B5D7DD", width=0)+
  geom_point(aes(x=Year, y=m), fill="#98E3D5", pch=21)+
  theme_bw()+ylab("LCW")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position="none", axis.title = element_text(size=14), 
        axis.text = element_text(size=14), strip.text = element_text(size = 14))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black"))+
  scale_color_distiller(palette="Pastel1")
Fig2A.2

df.fig2a3 <- dd %>% dplyr::select(Year, ID, Cushman_entropy_norm) %>% distinct() %>% 
  group_by(Year) %>% summarise(m=mean(Cushman_entropy_norm), s=sd(Cushman_entropy_norm)/sqrt(length(Year)) )

Fig2A.3 <- ggplot(df.fig2a3) + 
  geom_errorbar(aes(x=Year, ymin = m - 1.96*s, ymax = m + 1.96*s), color="#B5D7DD", width=0)+
  geom_point(aes(x=Year, y=m), fill="#98E3D5", pch=21)+
  theme_bw()+ylab("LCB")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position="none", axis.title = element_text(size=14), 
        axis.text = element_text(size=14), strip.text = element_text(size = 14))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black"))+
  scale_color_distiller(palette="Pastel1")
Fig2A.3


names(dd)
cor_df <- dd %>% dplyr::select(Shannon_entropy, Wasserstein_entropy, Cushman_entropy_norm,
                             Forest_cover, Forest_cooccur, Agri_cover, Agri_cooccur) 
colnames(cor_df) <- c( "LCS", "LCW", "LCB", "FC", "FCO" , "AC", "ACO")
col1 <- palette(hcl.colors(3, "Pastel1"))
pp.cols <- c("#98E3D5", "#47D7AC", "#41b6e6", "#B5D7DD", 
             "#8986CA", "#D9BFDC", "#F57EB6", "#F8C4CC", 
             "#EDDD87")

corrmatrix <- cor(cor_df, method = "spearman", use="pairwise.complete.obs")
library(ggplotify)
library(ggcorrplot)
pdf("~/Dropbox/Macro_BR/figs/Fig2F.pdf", width=7, height=4)
corrplot.mixed(corrmatrix, upper = "circle", lower = "number",
               tl.cex = 1.1, number.cex = 1, 
               tl.col = c("black","black","black","forestgreen","forestgreen","#F57EB6", "#F57EB6"), 
               font=2,
               lower.col = "black",
               number.font=1,
               upper.col = colorRampPalette(col1)(100) )
dev.off()

setwd("~/Dropbox/Macro_BR/figs/")
p <- ggplot()
grid2 <- plot_grid( Fig2A.1, Fig2A.2, Fig2A.3, p,
                     labels = c("A", "B", "C", ""), 
                     nrow=2, align="hv", label_size = 20, #vjust = 0,
                     label_fontface = "plain")
grid2

save_plot("Fig2.pdf", grid2,
          ncol = 2, 
          nrow = 2, 
          base_aspect_ratio = 1.1)



