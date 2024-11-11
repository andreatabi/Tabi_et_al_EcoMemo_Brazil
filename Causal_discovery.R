rm(list=ls()) 

source("~/Dropbox/Macro_BR/scripts/Toolbox.R")     
source("~/Dropbox/Macro_BR/scripts/github/IC.R") 
load("/Users/Andrea/Dropbox/Macro_BR/docs/df.RData")
load("/Users/Andrea/Dropbox/Macro_BR/docs/df_buffer.RData")
library(pcalg)
library(igraph)
library(dplyr)
library(ppcor)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(corrplot)

names(d)
names(dd)
dd1 <- dd %>% dplyr::select("ID", "Year", "sr", "Stream",
                     "cwm_pred", "cwm_size", 
                     "Shannon_entropy", "Wasserstein_entropy", "Boltzmann_entropy_rel",
                     "Cushman_entropy_norm", "Cushman_entropy", 
                     "mean_Shannon", "mean_Wasserstein", "mean_Cushman_norm", "mean_Cushman", "mean_Boltzmann_rel",
                     "Natural_cover", "Forest_cover",
                     "mean_Nat_cover", "mean_Forest_cover",
                     mean_Forest_cooccur, Forest_cooccur,
                     "mean_Agri_cover", "Agri_cover",
                     Agri_cooccur, mean_Agri_cooccur,
                     "patch_richness", "area",
                     "shading", "temp", "turb", "litter", "condut","PT", "NT",
                     "flow_veloc","OD", "width") %>% distinct()
d_buff <- merge(d, dd1 %>% dplyr::select(ID, Stream, sr) %>% distinct(), by="Stream")
d1 <- d_buff %>% dplyr::select("ID", "Year", "sr", Stream, area_buffer, 
                               patch_richness_buffer, Natural_cover_buffer, mean_Nat_cover_buffer,
                               Forest_cover_buffer, mean_Forest_cover_buffer, 
                               Forest_cooccur_buffer, mean_Forest_cooccur_buffer) %>% distinct()
data <- merge(dd1, d1, by=c("Year", "ID", "Stream", "sr") )
data$Area <- round(data$area)


names(data)
col1 <- colorRampPalette(brewer.pal(9,"BrBG"))
#col1 <- palette(hcl.colors(10, "Cyan-Magenta"))
d <- data %>%  filter(Year==2015) %>% dplyr::select("ID", "sr", 
                                                   "cwm_pred", "cwm_size", 
                                                   "mean_Shannon", "mean_Wasserstein", "mean_Cushman_norm", 
                                                   "mean_Forest_cover", 
                                                   "mean_Agri_cover",
                                                   "patch_richness", "Area",
                                                   "shading", "temp", "turb", "litter", "condut","PT", "NT",
                                                   "flow_veloc","OD", "width") %>% distinct()
#d <- aggregate(d0[,-1], by=list(d0$ID), function(x) mean(x, na.rm=TRUE))
d <- d[,-1]
colnames(d) <- c( "Species richness", 
                  "Predation CWM", "Body size CWM", 
                  "LCS", "LCW","LCB", 
                  "FC", 
                  "AC",
                  "Patch richness", "Catchment area",
                  "Shading", "Temperature", "Turbidity", "Organic detritus", "Conductivity", "Phosphorus", "Nitrates",
                  "Flow velocity", "Dissolved oxygen", "Stream width")
corrmatrix <- cor(d, method = "spearman")
pdf("~/Dropbox/Macro_BR/figs/FigS4.pdf", width=8, height=6)
corrplot.mixed(corrmatrix, upper = "square", lower = "number",
               tl.cex = 0.75, lower.col = "black", number.cex = .5, 
               tl.pos = "lt", tl.col = "black", tl.offset=1, diag="u", #font=1,
               number.font=1,
               upper.col = col1(100) )
dev.off()


names(dd)
alpha <- seq(0.01, 0.1, by=0.01)
variables <- c("ID", "sr", 
               "Shannon_entropy", "Wasserstein_entropy", "Cushman_entropy_norm",
               "cwm_size", "cwm_pred", "litter","shading",
               "temp",  "PT", "NT"
               )  

for(j in 1:length(alpha) ){
          #j <- 1
          print(j)
              names(data)
              dat0 <- data %>% dplyr::filter(Year==2015) %>% dplyr::select(variables) %>% distinct() 
              dat0[is.na(dat0)] <- 0
              dat0 <- aggregate(dat0[,-1], by=list(dat0$ID), function(x) mean(x, na.rm=TRUE))
              cor.test(dat0$sr, dat0$Shannon_entropy)
              cor.test(dat0$sr, dat0$cwm_pred)
              pcor.test(dat0$sr, dat0$cwm_pred, dat0$cwm_size, method = "spearman")
              
              df <- list()
              dg <- IC(dat0[,-1], alpha = alpha[j])
              plot(dg$G, layout=layout.auto, main=paste0("alpha = ", alpha[j] ),
                     vertex.size = 1,      
                     vertex.label = V(dg$G)$name,
                     vertex.label.cex = 1,
                     vertex.label.color = "black",
                     vertex.color="white" )
              out <- do.call("rbind", strsplit(dg$Skeleton, "--"))
              links <- paste0(out[,1], "_", out[,2])
              targ <- ifelse( any(links=="1_2"), 1, 0)
              df <- data.frame(link1=out[,1], link2=out[,2], alpha=alpha[j], all_links=nrow(out), target_link=targ)
         
         if(j==1) dat <- df
         if(j>1)  dat <- rbind(dat, df)
}

cg <- bind_rows(dat)
head(cg)

#save(cg, file ="/Users/Andrea/Dropbox/Macro_BR/scripts/github/cg.RData") 

dat0 <- data %>% dplyr::filter(Year==2015) %>% dplyr::select(variables) %>% distinct()
lab <- sample(colnames(dat0))
dat1 <- dat0[, lab]
names(dat1)
suffStat <- list(C = cor(dat1, method = "spearman"), n = nrow(dat1))
pc.B <- pc(suffStat,
            indepTest = gaussCItest, alpha = 0.06, labels = lab, verbose = T,
            skel.method = "stable",
            #  u2pd = "retry"
            maj.rule = T
            # conservative =T
 )
 
require(Rgraphviz)
#dev.off()  
#plot(pc.B, main = "Estimated DAG", cex.text=3)
iplotPC(pc.B)

### FIGURE 4 --------------------------------------------------------------------------------------------------

pp.cols <- c("#98E3D5", "#47D7AC", "#41b6e6", "#B5D7DD", 
             "#8986CA", "#D9BFDC", "#F57EB6", "#F8C4CC", 
             "#EDDD87")

cg
num_of_links <- lapply(alpha, function(x){
              #x <- 0.01
                temp <- cg %>% filter(alpha==x)
                l1 <- length(which(temp$link1 ==1)) 
                l2 <- length(which(temp[ which(temp$link1 ==1), "link2"] %in% temp[-which(temp$link1 ==1),"link1"]))
                l3 <- length(which(temp[ which(temp$link1 ==1), "link2"] %in% temp[-which(temp$link1 ==1),"link2"]))
                l <- l1+l2+l3
})

df <- data.frame(alpha=alpha, nof=unlist(num_of_links))

ggplot(df,aes(x=factor(alpha), y=nof, group=1 )) + 
  annotate("rect", xmin = 9, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#98E3D5", alpha=0.2 ) +
  geom_point( pch=21, size=2)+
  geom_line(linewidth=0.5, color="black")+
  theme_classic()+  
  theme(legend.position="none", axis.title = element_text(size=14), 
                          axis.text = element_text(size=10), strip.text = element_text(size = 14))+
  geom_vline(xintercept=9,  colour="#47D7AC", linetype = "longdash")+
  ylab("Number of causal links")+ 
  xlab(bquote(alpha))

ggsave("~/Dropbox/Macro_BR/figs/Fig4A.pdf", width=4, height = 3.5)

### Causal discovery with memory -------------------------------------------------------------
load("~/Dropbox/Macro_BR/scripts/github/Memo_weights_SE.RData") 
se <- res
mean_se <- lapply(seq(31), function(x){
  weights <- list()
  for(j in seq(n)){
    # j <- x <- 31
    weights[[j]] <- se[[j]]$w[x]  
  }
  w <- unlist(weights)
  mean <- mean(w)
  return(data.frame(mean))
})
mean_se
df0 <- data %>% dplyr::select(ID, Stream, sr, Year, Shannon_entropy)
df_se0 <- reshape(df0, idvar = c( "ID", "sr", "Stream"), timevar = "Year", direction = "wide")
df_se <- df_se0[,-c(1:3)]
SEmemo <- c(unlist(mean_se) %*% t(df_se))
df.se <- data.frame(ID=df_se0$ID, Stream=df_se0$Stream, sr=df_se0$sr,  SEmemo=SEmemo)

data_memo <- merge(data, df.se, by=c( "ID", "Stream", "sr") )

load("~/Dropbox/Macro_BR/scripts/github/Memo_weights_WE.RData") 
we <- res
mean_we <- lapply(seq(31), function(x){
  weights <- list()
  for(j in seq(n)){
    # j <- x <- 31
    weights[[j]] <- we[[j]]$w[x]  
  }
  w <- unlist(weights)
  mean <- mean(w)
  return(data.frame(mean))
})
mean_we
df0 <- data %>% dplyr::select(ID, Stream, sr, Year, Wasserstein_entropy)
df_we0 <- reshape(df0, idvar = c( "ID", "sr", "Stream"), timevar = "Year", direction = "wide")
df_we <- df_we0[,-c(1:3)]
WEmemo <- c(unlist(mean_we) %*% t(df_we))
df.we <- data.frame(ID=df_we0$ID, Stream=df_we0$Stream, sr=df_we0$sr,  WEmemo=WEmemo)

data_memo <- merge(data_memo, df.we, by=c( "ID", "Stream", "sr") )

load("~/Dropbox/Macro_BR/scripts/github/Memo_weights_CE.RData") 
ce <- res
mean_ce <- lapply(seq(31), function(x){
  weights <- list()
  for(j in seq(n)){
    # j <- x <- 31
    weights[[j]] <- ce[[j]]$w[x]  
  }
  w <- unlist(weights)
  mean <- mean(w)
  return(data.frame(mean))
})
mean_ce
df0 <- data %>% dplyr::select(ID, Stream, sr, Year,Cushman_entropy_norm)
df_ce0 <- reshape(df0, idvar = c( "ID", "sr", "Stream"), timevar = "Year", direction = "wide")
df_ce <- df_ce0[,-c(1:3)]
CEmemo <- c(unlist(mean_ce) %*% t(df_ce))
df.ce <- data.frame(ID=df_ce0$ID, Stream=df_ce0$Stream, sr=df_ce0$sr,  CEmemo=CEmemo)

data_memo <- merge(data_memo, df.ce, by=c( "ID", "Stream", "sr") )
names(data_memo)

alpha <- seq(0.01, 0.1, by=0.01)

variables <- c("sr", 
               "SEmemo", "WEmemo","CEmemo",
               "cwm_size", "cwm_pred", "litter","shading", "temp",  "PT")  

for(j in 1:length(alpha)){
  #j <- 7
  print(j)
  names(data)
  dat0 <- data_memo %>% dplyr::filter(Year==2015) %>% dplyr::select(variables) %>% distinct() 
  dat0[is.na(dat0)] <- 0
  df <- list()
  dg <- IC(dat0, alpha = alpha[j])
  plot(dg$G, layout=layout.auto, main=paste0("alpha = ", alpha[j] ),
       vertex.size = 1,      
       vertex.label = V(dg$G)$name,
       vertex.label.cex = 1,
       vertex.label.color = "black",
       vertex.color="white" )
  out <- do.call("rbind", strsplit(dg$Skeleton, "--"))
  links <- paste0(out[,1], "_", out[,2])
  targ <- ifelse( any(links=="1_2"), 1, 0)
  df <- data.frame(link1=out[,1], link2=out[,2], alpha=alpha[j], all_links=nrow(out), target_link=targ)
  
  if(j==1) dat <- df
  if(j>1)  dat <- rbind(dat, df)
}

cg.memo <- bind_rows(dat)
head(cg.memo)


num_of_links <- lapply(alpha, function(x){
  #x <- 0.01
  temp <- cg.memo %>% filter(alpha==x)
  l1 <- length(which(temp$link1 ==1)) 
  l2 <- length(which(temp[ which(temp$link1 ==1), "link2"] %in% temp[-which(temp$link1 ==1),"link1"]))
  l3 <- length(which(temp[ which(temp$link1 ==1), "link2"] %in% temp[-which(temp$link1 ==1),"link2"]))
  l <- l1+l2+l3
})

df.memo <- data.frame(alpha=alpha, nof=unlist(num_of_links))

ggplot(df.memo,aes(x=factor(alpha), y=nof, group=1 )) + 
  annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#98E3D5", alpha=0.2 ) +
  geom_point( pch=21, size=2)+
  geom_line(linewidth=0.5, color="black")+
  theme_classic()+  
  theme(legend.position="none", axis.title = element_text(size=14), 
        axis.text = element_text(size=10), strip.text = element_text(size = 14))+
  geom_vline(xintercept=1,  colour="#47D7AC", linetype = "longdash")+
  ylab("Number of causal links")+ 
  xlab(bquote(alpha))

ggsave("~/Dropbox/Macro_BR/figs/Fig4B.pdf", width=4, height = 3.5)

