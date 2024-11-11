rm(list=ls()) 
load("~/Dropbox/Macro_BR/docs/df.RData")
library(cowplot)
library(grid)
library(gridExtra)
library(ggplot2)
setwd("~/Dropbox/Macro_BR/figs/")

names(d)
names(dd)
data <- dd %>% dplyr::select("ID", "Year", "sr", "Stream",
                            "cwm_pred", "cwm_size", 
                            "Shannon_entropy", "Wasserstein_entropy", "Boltzmann_entropy",
                            "Cushman_entropy_norm", 
                            "mean_Shannon", "mean_Wasserstein", "mean_Cushman_norm", 
                            "Natural_cover", "Forest_cover", "Forest_cooccur",
                            "Agri_cover", "Agri_cooccur",
                            "mean_Nat_cover", "mean_Forest_cover", "mean_Forest_cooccur",
                            "patch_richness", "area",
                            "shading", "temp", "turb", "litter", "condut","PT", "NT",
                            "flow_veloc","OD", "width") %>% distinct()
data$Area <- round(data$area)

## Landscape composition -------------------------------------------------------
load("~/Dropbox/Macro_BR/scripts/github/Memo_weights_SE.RData") 

res
n <- 500
cors <- lapply(seq(n), function(x) res[[x]]$rho)
cors_rand <- lapply(seq(n), function(x) -res[[x]]$rho_random)

## calc CI for weights

CI <- lapply(seq(31), function(x){
  weights <- list()
  for(j in seq(n)){
    # j <- x <- 31
    weights[[j]] <- res[[j]]$w[x]  
  }
  w <- unlist(weights)
  margin <- 1.96*sd(w)/sqrt(length(w))
  mean <- mean(w)
  return(data.frame(mean,margin))
})
CI

ci_df <- do.call( "rbind",CI)

df0 <- data %>% dplyr::select(ID, Stream, sr, Year, Shannon_entropy) %>% filter(Year==2015) 
df <- reshape(df0, idvar = c("sr", "ID", "Stream"), timevar = "Year", direction = "wide")
df_var <- df[,-c(1:2)]
cor_cc <- cor.test(df0$sr, df0$Shannon_entropy, method = "spearman" )$est

pp.cols <- c("#98E3D5", "#47D7AC", "#41b6e6", "#B5D7DD", 
             "#8986CA", "#D9BFDC", "#F57EB6", "#F8C4CC", 
             "#EDDD87", "#a4b6dd")
dat <- data.frame(c=unlist(cors), c2=unlist(cors_rand))

SE1 <- ggplot(dat)+
  geom_histogram(aes(x=c), color="white", fill="#98E3D5")+
  geom_histogram(aes(x=c2), color="white", fill="#F8C4CC")+
  theme_classic()+
  geom_vline(aes(xintercept=cor_cc), size=1, linetype="dotted", color=pp.cols[7])+
  ylab("Count")+xlab("")
SE1

dat_w <- data.frame(time=c(1985:2015), mean_w=ci_df$mean, margin=ci_df$margin )
SE2 <- ggplot(dat_w)+
  geom_line(aes(x=time, y=mean_w), color=pp.cols[6])+
  geom_errorbar( aes(x=time, ymin=mean_w-margin, ymax=mean_w+margin), width=0, colour="grey20", alpha=0.9, size=0.5)+
  geom_point(aes(x=time, y=mean_w), pch=21, color="black", fill="grey", size=2)+
  theme_classic()+xlab(" ")+ ylab("Memory weights (LCS)")
SE2

df2 <- data %>% dplyr::select(ID, Stream,  Year, Shannon_entropy) %>% distinct()
df2 <- reshape(df2, idvar = c( "ID", "Stream"), timevar = "Year", direction = "wide")
se_opt <- c(dat_w$mean_w %*% t(df2[,c(-1, -2)]))
cor.test( df0$sr, se_opt, method="spearman")

SE3 <- ggplot(df0)+
  geom_point(aes(x=sr, y=se_opt), pch=21, color="grey", fill=pp.cols[1], size=2)+
  geom_smooth(aes(x=sr, y=se_opt), method="lm",color="grey50", size=0.5, se=T, 
              linetype="dashed")+
  theme_classic()+xlab(" ") +   ylab(expression(paste("LCS"[opt] )))

SE3


grid1 <- plot_grid( SE1,SE2,SE3,
                    labels = c("A", "", ""), 
                    nrow=1, align="hv", label_size = 20, #vjust = 0,
                    label_fontface = "plain")
grid1

##  Landscape configuration (W) -------------------------------------------------------
load("~/Dropbox/Macro_BR/scripts/github/Memo_weights_WE.RData") 

n <- 500
cors <- lapply(seq(n), function(x) res[[x]]$rho)
cors_rand <- lapply(seq(n), function(x) res[[x]]$rho_random)

## calc CI for weights

CI <- lapply(seq(31), function(x){
  weights <- list()
  for(j in seq(n)){
    # j <- x <- 31
    weights[[j]] <- res[[j]]$w[x]  
  }
  w <- unlist(weights)
  margin <- 1.96*sd(w)/sqrt(length(w))
  mean <- mean(w)
  return(data.frame(mean,margin))
})
CI

ci_df <- do.call( "rbind",CI)

df0 <- data %>% dplyr::select(ID, sr, Stream, Year, Wasserstein_entropy) %>% filter(Year==2015)                                                                                                                       
cor_cc <- cor.test(df0$sr, df0$Wasserstein_entropy, method = "spearman" )$est

dat <- data.frame(c=unlist(cors), c2=unlist(cors_rand))
WE1 <- ggplot(dat)+
  geom_histogram(aes(x=c), color="white", fill="#98E3D5")+
  geom_histogram(aes(x=c2), color="white", fill="#F8C4CC")+
  theme_classic()+xlab(" ")+ylab("Count")+
  geom_vline(aes(xintercept=cor_cc), size=1, linetype="dotted", color=pp.cols[7])
WE1

dat_w <- data.frame(time=c(1985:2015), mean_w=ci_df$mean, margin=ci_df$margin )
WE2 <- ggplot(dat_w)+
  geom_line(aes(x=time, y=mean_w), color=pp.cols[6])+
  geom_errorbar( aes(x=time, ymin=mean_w-margin, ymax=mean_w+margin), width=0, colour="grey20", alpha=0.9, size=0.5)+
  geom_point(aes(x=time, y=mean_w), pch=21, color="black", fill="grey", size=2)+
  theme_classic()+xlab(" ")+ ylab("Memory weights (LCW)")
WE2

df2 <- data %>% dplyr::select(ID, Stream,  Year, Wasserstein_entropy) %>% distinct()
df2 <- reshape(df2, idvar = c( "ID", "Stream"), timevar = "Year", direction = "wide")
se_opt <- c(dat_w$mean_w %*% t(df2[,c(-1, -2)]))
cor.test( df0$sr, se_opt, method="spearman")

WE3 <- ggplot(df0)+
  geom_point(aes(x=sr, y=se_opt), pch=21, color="grey", fill=pp.cols[1], size=2)+
  geom_smooth(aes(x=sr, y=se_opt), method="lm",color="grey50", size=0.5, se=T, linetype="dashed")+
  theme_classic()+xlab(" ") +   ylab(expression(paste("LCW"[opt] )))
WE3

grid2 <- plot_grid( WE1,WE2,WE3,
                    labels = c("B", "", ""), 
                    nrow=1, align="hv", label_size = 20, #vjust = 0,
                    label_fontface = "plain")
grid2

## Landscape configuration (C) -------------------------------------------------------
load("~/Dropbox/Macro_BR/scripts/github/Memo_weights_CE.RData") 

n <- 500
cors <- lapply(seq(n), function(x) res[[x]]$rho)
cors_rand <- lapply(seq(n), function(x) res[[x]]$rho_random)

## calc CI for weights

CI <- lapply(seq(31), function(x){
  weights <- list()
  for(j in seq(n)){
    # j <- x <- 31
    weights[[j]] <- res[[j]]$w[x]  
  }
  w <- unlist(weights)
  margin <- 1.96*sd(w)/sqrt(length(w))
  mean <- mean(w)
  return(data.frame(mean,margin))
})
CI

ci_df <- do.call( "rbind",CI)

df0 <- data %>% dplyr::select(Stream, sr, Year, Cushman_entropy_norm) %>% filter(Year==2015)
cor_cc <- cor.test(df0$sr, df0$Cushman_entropy_norm, method = "spearman" )$est

dat <- data.frame(c=unlist(cors), c2=unlist(cors_rand))
CE1 <- ggplot(dat)+
  geom_histogram(aes(x=c), color="white", fill="#98E3D5")+
  geom_histogram(aes(x=c2), color="white", fill="#F8C4CC")+
  theme_classic()+xlab(expression(paste("Spearman's ",rho )))+
  geom_vline(aes(xintercept=cor_cc), size=1, linetype="dotted", color=pp.cols[7])+
  ylab("Count")
CE1

dat_w <- data.frame(time=c(1985:2015), mean_w=ci_df$mean, margin=ci_df$margin )
CE2 <- ggplot(dat_w)+
  geom_line(aes(x=time, y=mean_w), color=pp.cols[6], )+
  geom_errorbar( aes(x=time, ymin=mean_w-margin, ymax=mean_w+margin), width=0, colour="grey20", alpha=0.9, size=0.5)+
  geom_point(aes(x=time, y=mean_w), pch=21, color="black", fill="grey", size=2)+
  theme_classic()+xlab(" ")+ ylab("Memory weights (LCB)")+
  xlab("Time")
CE2

df2 <- data %>% dplyr::select(ID, Stream, Year, Cushman_entropy_norm) %>% distinct()
df2 <- reshape(df2, idvar = c( "ID", "Stream"), timevar = "Year", direction = "wide")
se_opt <- c(dat_w$mean_w %*% t(df2[,c(-1, -2)]))
cor.test( df0$sr, se_opt, method="spearman")

CE3 <- ggplot(df0)+
  geom_point(aes(x=sr, y=se_opt), pch=21, color="grey", fill=pp.cols[1], size=2)+
  geom_smooth(aes(x=sr, y=se_opt), method="lm",color="grey50", size=0.5, se=T, linetype="dashed")+
  theme_classic()+xlab(" ") + 
  ylab(expression(paste("LCB"[opt] )))+
  xlab("Species richness")
CE3

grid3 <- plot_grid( CE1,CE2,CE3,
                    labels = c("C", "", ""), 
                    nrow=1, align="hv", label_size = 20, #vjust = 0,
                    label_fontface = "plain")
grid3

## Agri cover  -------------------------------------------------------
load("~/Dropbox/Macro_BR/scripts/github/Memo_weights_AC.RData") 

n <- 500
cors <- lapply(seq(n), function(x) res[[x]]$rho)
cors_rand <- lapply(seq(n), function(x) res[[x]]$rho_random)

## calc CI for weights

CI <- lapply(seq(31), function(x){
  weights <- list()
  for(j in seq(n)){
    # j <- x <- 31
    weights[[j]] <- res[[j]]$w[x]  
  }
  w <- unlist(weights)
  margin <- 1.96*sd(w)/sqrt(length(w))
  mean <- mean(w)
  return(data.frame(mean,margin))
})
CI

ci_df <- do.call( "rbind",CI)

df0 <- data %>% dplyr::select(Stream, sr, Year, Agri_cover) %>% filter(Year==2015)
df <- reshape(df0, idvar = c("sr", "Stream"), timevar = "Year", direction = "wide")
cor_cc <- cor.test(df0$sr, df0$Agri_cover, method = "spearman" )$est

dat <- data.frame(c=unlist(cors), c2=unlist(cors_rand))
AC1 <- ggplot(dat)+
  geom_histogram(aes(x=c), color="white", fill="#98E3D5")+
  geom_histogram(aes(x=c2), color="white", fill="#F8C4CC")+
  theme_classic()+xlab(expression(paste("Spearman's ",rho)))+
  geom_vline(aes(xintercept=cor_cc), size=1, linetype="dotted", color=pp.cols[7])+
  ylab("Count")
AC1

dat_w <- data.frame(time=c(1985:2015), mean_w=ci_df$mean, margin=ci_df$margin )
AC2 <- ggplot(dat_w)+
  geom_line(aes(x=time, y=mean_w), color=pp.cols[6], )+
  geom_errorbar( aes(x=time, ymin=mean_w-margin, ymax=mean_w+margin), width=0, colour="grey20", alpha=0.9, size=0.5)+
  geom_point(aes(x=time, y=mean_w), pch=21, color="black", fill="grey", size=2)+
  theme_classic()+xlab(" ")+ ylab("Memory weights")+
  xlab("Time")
AC2

df2 <- data %>% dplyr::select(ID, Stream, Year, Agri_cover) %>% distinct()
df2 <- reshape(df2, idvar = c( "ID", "Stream"), timevar = "Year", direction = "wide")
se_opt <- c(dat_w$mean_w %*% t(df2[,c(-1, -2)]))
cor.test( df0$sr, se_opt, method="spearman")

AC3 <- ggplot(df0)+
  geom_point(aes(x=sr, y=se_opt), pch=21, color="grey", fill=pp.cols[1], size=2)+
  geom_smooth(aes(x=sr, y=se_opt), method="lm",color="grey50", size=0.5, se=T, linetype="dashed")+
  theme_classic()+xlab(" ") + 
  ylab(expression(paste("AC"[opt] )))+
  xlab("Species richness")
AC3


grid <- plot_grid( grid1, grid2, grid3,
                   labels = c("", "", "", ""), 
                   nrow=3, align="hv", label_size = 20, #vjust = 0,
                   label_fontface = "plain" )
grid

save_plot("Fig3.pdf", grid,
          ncol = 3, 
          nrow = 3, 
          scale = 0.8,
          base_aspect_ratio = 1.1)


###


grid4 <- plot_grid( AC1,AC2,AC3,
                    labels = c("A", "B", "C"), 
                    nrow=1, align="hv", label_size = 20, #vjust = 0,
                    label_fontface = "plain")
grid4
save_plot("FigS8.pdf", grid4,
          ncol = 3, 
          nrow = 1, 
          scale = 0.8,
          base_aspect_ratio = 1.1)

