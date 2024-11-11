rm(list=ls()) 

library(neighbours)  
library(NMOF)        
library(cowplot)
library(grid)
library(gridExtra)
library(dplyr)
  

load("~/Dropbox/Macro_BR/docs/df.RData")

names(d)
names(dd)
data <- dd %>% dplyr::select("ID", "Year", "sr", "Stream","size","even",
                            "cwm_pred", "cwm_size", 
                            "Shannon_entropy", "Wasserstein_entropy", "Boltzmann_entropy",
                            "Cushman_entropy_norm", 
                            "mean_Shannon", "mean_Wasserstein", "mean_Cushman_norm", 
                            "all_Forest_cover",
                            "Natural_cover", "Forest_cover", "Forest_cooccur",
                            "Agri_cover", "Agri_cooccur",
                            "mean_Nat_cover", "mean_Forest_cover", "mean_Forest_cooccur",
                            "patch_richness", "area",
                            "shading", "temp", "turb", "litter", "condut","PT", "NT",
                            "flow_veloc","OD", "width") %>% distinct()
data$Area <- round(data$area)

names(data)
var <- "Shannon_entropy"

df0 <- data %>% dplyr::select(ID, Stream, sr, Year, var) %>% arrange(Year)
df0 <- reshape(df0, idvar = c("sr", "ID", "Stream"), timevar = "Year", direction = "wide")
df_var <- df0[,-c(1:3)]

Data <- list( R = t(df_var),
              na = dim(df_var)[2],
              ns = dim(df_var)[1],
              eps= 0.1/100,
              wmin=0,
              wmax=1, 
              resample = function(x, ...) x[sample.int(length(x), ...)] )

OF <- function(w, Data) {
  -c(cor.test( df0$sr, c(w %*% Data$R), method="spearman")$est)
}

neighbor <- function(w, Data){
  toSell <- w > Data$wmin
  toBuy <- w < Data$wmax
  i <- Data$resample(which(toSell), size=1)
  j <- Data$resample(which(toBuy), size=1)
  eps <- runif(1) * Data$eps
  eps <- min(w[i] - Data$wmin, Data$wmax - w[j], eps)
  w[i] <- w[i] - eps
  w[j] <- w[j] + eps
  return(w)
}

num <- 500
iter <- 500000
res <- lapply(seq(num), function(i){
          print(i)
          w0 <- runif(31)
          w0 <- w0/sum(w0)
          algo <- list(x0 = w0,
                       neighbour = neighbor,
                       nI=iter,
                       printBar = FALSE
                       )
          out <- TAopt(OF, algo, Data)
          init <- -OF(w0, Data)
          df <- list(n=i, w=out$xbest, rho=-out$OFvalue, rho_random= init , iteration=iter)
})

res

save(res, file ="~/Memo_weights_SE.RData") 


CI <- lapply(seq(31), function(x){
  weights <- list()
  for(j in seq(num)){
    weights[[j]] <- res[[j]]$w[x]  
  }
  w <- unlist(weights)
  margin <- 1.96*sd(w)/sqrt(length(w))
  mean <- mean(w)
  return(data.frame(mean,margin))
})
CI
ci_df <- do.call( "rbind",CI)

plot(c(1985:2015), ci_df$mean, type="b")
abline(lm(ci_df$mean ~ c(1985:2015) ))

