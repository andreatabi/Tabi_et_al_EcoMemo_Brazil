library(pcalg)
library(igraph)
# V <- dat0; type="con"; alpha=0.05
IC <- function(V, alpha=0.05, type="con"){
  
  l <- ncol(V)
  co <- combn(l, 2)
  E <- v_struc <- list()
  
  # Step 1: finding the skeleton
  
  for(i in 1:ncol(co)){
    #i <- 1
    vars <- co[,i]
    a <- vars[1]
    b <- vars[2]
    c <- seq(l)[-c(a,b)]
    
    if(length(c)==1){ 
      S <- c
    } else{
      S <- do.call("c", lapply(seq(length(c)), function(i) combn(c, i, FUN = list)))
    }
    
    if(type=="bin"){
      t1 <- gSquareBin(a,b,c(),V)   # null hypothesis -> independence
      t2 <- lapply(S, function(x) gSquareBin(a,b,c(x),V) )
    }
    
    if(type=="con"){
      t1 <- cor.test(V[,a], V[,b], method = c("spearman"))$p.value 
      if(t1<alpha){
          t2 <- lapply(S, function(x) pcor.test(V[,a],V[,b],c(V[,x]), method = c("spearman"))$p.value )
          E[[i]] <- ifelse(all(unlist(t2) < alpha), 1, 0)
      }else{
          E[[i]] <- 0
      }
      

      if(t1 > alpha){
        t3 <- lapply(c, function(x) pcor.test(V[,a],V[,b],c(V[,x]),  method = c("spearman"))$p.value )
        if(any(unlist(t3) < alpha)) v_struc[[i]] <- paste0(a, "-->", c[which(unlist(t3) < alpha)], "<--", b)
      } 
    }
  }
  
  Edges <- as.matrix(co[, which( unlist(E)==1 )], bycol=T)
  if( all(unlist(E)==0) ) return("no edges")
  
  Sk <- lapply(seq(ncol(Edges)), function(y) paste0(Edges[1,y], "--", Edges[2,y]  ))
  Sk <- unlist(Sk)
  
  # Step 2: finding colliders
  
  v_struc <- unlist(v_struc)
  if( length(v_struc) > 0 ){
    check_neighbor_V <- lapply(seq(length(v_struc)), function(z){ 
      vec <- as.numeric(strsplit(v_struc[z], "-->|<--" )[[1]])
      a <- vec[1]
      c <- vec[2]
      b <- vec[3]
      neigh1 <- lapply(seq(ncol(Edges)), function(u) ifelse( all(c(a,c) %in% Edges[,u]),1,0))
      neigh2 <- lapply(seq(ncol(Edges)), function(u) ifelse( all(c(b,c) %in% Edges[,u]),1,0))
      N <- ifelse( sum(unlist(neigh1))==1 & sum(unlist(neigh2))==1, 1, 0)
    }) 
    V_struc <- v_struc[which(check_neighbor_V==1)]
  } else{
    V_struc <- "none"
  } 
  
  
  # Step 3: allocating directions to remaining undirected edges
  
  if(  !(length(V_struc)==0) &  all(!(V_struc == "none")) ){
    v0 <- lapply(seq(length(V_struc)), function(g){
      n <- as.numeric(strsplit(V_struc[g], "-->|<--" )[[1]]) 
      m <- matrix( c(n[1], n[2], n[3], n[2]), ncol=2)  })
    mat_E <- unique(do.call("cbind", v0), MARGIN=2)
    rem_E <- lapply(seq(ncol(Edges)), function(u) 
      ifelse( any(unlist(lapply(seq(ncol(mat_E)), function(k) all(mat_E[,k] %in% Edges[,u]) ))) ,0,1))
    remaining_edges <- matrix(Edges[,which(unlist(rem_E)==1)], nrow=2)
    
    n <- 0
    counter <- 0
    new_edge <- rep(NA, ncol(remaining_edges) )
    while( n < ncol(remaining_edges) && counter < ncol(remaining_edges)^2 ){
      for(z in seq(ncol(remaining_edges)) ){
        #z <- 1  
        if( is.na(new_edge[[z]]) ){  
          vec <- remaining_edges[,z]
          b <- vec[1]
          c <- vec[2]
          
          # 1. Orient b − c into b → c whenever there is a → b such that a and c are non-adjacents
          neigh_b <- lapply(seq(ncol(mat_E)), function(u){ 
            x <- match(b , mat_E[,u])
            adj <- lapply(seq(ncol(Edges)), function(m) ifelse( all(c(mat_E[1,u],c) %in% Edges[,m]),1,0))
            ifelse(is.na(x) | any(unlist(adj)==1), 0, 1)  
          })
          neigh_c <- lapply(seq(ncol(mat_E)), function(u){ 
            x <- match(c , mat_E[,u])
            adj <- lapply(seq(ncol(Edges)), function(m) ifelse( all(c(mat_E[1,u],b) %in% Edges[,m]),1,0))
            ifelse(is.na(x) | any(unlist(adj)==1), 0, 1)  
          })
          
          if( any(unlist(neigh_b)==1) ){
            new_edge[z] <- paste0(b, "-->", c)
            mat_E <- cbind(mat_E, c(b,c))
          } 
          
          if( any(unlist(neigh_c)==1) ){  
            new_edge[z] <- paste0(c, "-->", b)
            mat_E <- cbind(mat_E, c(c,b))
          } 
          
          # 2. Orient b − c into b → c whenever there is a chain b → a → c
          
          chain_b <- lapply(seq(ncol(mat_E)), function(u){ 
            x1 <- match(b , mat_E[1,u])
            x2 <- match(c , mat_E[2,u])
            ifelse(is.na(x1) | is.na(x2), 0, 1)  
          })
          
          chain_c <- lapply(seq(ncol(mat_E)), function(u){ 
            x1 <- match(c , mat_E[1,u])
            x2 <- match(b , mat_E[2,u])
            ifelse(is.na(x1) | is.na(x2), 0, 1)  
          })
          
          if( any(unlist(chain_b)==1) ){
            new_edge[z] <- paste0(b, "-->", c)
            mat_E <- cbind(mat_E, c(b,c))
          } 
          
          if( any(unlist(chain_c)==1) ){  
            new_edge[z] <- paste0(c, "-->", b)
            mat_E <- cbind(mat_E, c(c,b))
          } 
        }    
      }
      
      # 3. Orient a — b into a -> b whenever there are two chains a — c -> b and 
      #    a — d -> b such that c and d are nonadjacent.
      
      chain_acb <- lapply(seq(ncol(mat_E)), function(u){ 
        x1 <- match(b , mat_E[1,u])
        x2 <- match(c , mat_E[2,u])
        ifelse(is.na(x1) | is.na(x2), 0, 1)  
      })
      
      chain_adb <- lapply(seq(ncol(mat_E)), function(u){ 
        x1 <- match(b , mat_E[1,u])
        x2 <- match(c , mat_E[2,u])
        ifelse(is.na(x1) | is.na(x2), 0, 1)  
      })
      
      
      # 4. Orient a — b into a -> b whenever there are two chains a — c -> d and 
      #    c —> d -> b such that c and b are nonadjacent and a and d are adjacent.
      
      n <- length(na.omit(new_edge))
      #  print(counter)
      counter <- sum(counter, 1)
    }  
    
    ## Directed edges
    dir_E <- unique(mat_E, MARGIN=2)
    dir_edges <- c(paste0(t(dir_E)[,1], "-->", t(dir_E)[,2]))
    
  }else{
    dir_E <- NA
    dir_edges <- 0
  }    
  
  
  which_undir_E <- lapply(seq(ncol(Edges)), function(u) 
    ifelse( any(unlist(lapply(seq(ncol(dir_E)), function(k) all(dir_E[,k] %in% Edges[,u]) ))) ,0,1))
  undir_E1 <- matrix(Edges[,which(unlist(which_undir_E)==1)], nrow=2)
  undir_E2 <- rbind(undir_E1[2,], undir_E1[1,])
  undir_E <- cbind(undir_E1, undir_E2)
  
  all_edges <- cbind(dir_E, undir_E)
  sep_ver <- which(!(1:ncol(V) %in% c(all_edges)))
  sep_E <- rbind(sep_ver, sep_ver)
  all_edges <- cbind(all_edges, sep_E)
  
  # plotting
  all_edges <- na.omit(t(all_edges))
  g <- igraph::graph.data.frame(all_edges, directed=T  )
  V(g)$name <- colnames(V)[as.numeric(V(g)$name)]
  
  return(list(Skeleton=Sk, V_Struc=V_struc, Dir_Edges=dir_edges, G=g))
  
}
