genBaseProbs <- function(n, base, similarity, digits = 8) {
  
  n_levels <- length(base)
  x <- rdirichlet(n, similarity * base) 
  
  #--- ensure that each vector of probabilities sums exactly to 1
  
  x <- round(floor(x*1e8)/1e8, digits)   # round the generated probabilities
  xpart <- x[, 1:(n_levels-1)]           # delete the base prob of the final level
  partsum <- apply(xpart, 1, sum)        # add the values of levels 1 to K-1
  x[, n_levels] <- 1 - partsum           # the base prob of the level K = 1 - sum(1:[K-1])
  
  return(x)
}

make_data <- function(nreef,nyear){
  
  d_study <- genData(nreef, id = "reef")
  d_ind <- genCluster(d_study, cLevelVar = "reef", numIndsVar = 500, level1ID = "id") # 500 observations per reef
  d_ind[, z := 0] # set to zero to be used after 
  
  
  # Get the probability of occurrence for each year 
  basestudy <- genBaseProbs(
    n = nyear,
    base =  c(0.4, 0.1, 0.5), # prob of occurrence 
    similarity = 50 # variation of prob across years (highest number = more similarities between years)
  )
  
  list_ind_y <- list() 
  abundance_matrix <- list()
  
  for(j in 1:nyear) {
    list_ind_y[[j]] <- lapply(
      X = 1:nreef, 
      function(i) {
        b <- basestudy[j,]
        d_x <- d_ind[reef == i]
        genOrdCat(d_x, adjVar = "z", b, catVar = "ordY")
      })
    
    list_ind_y <- list_ind_y[[j]]
    abundance_matrix[[j]] <- t(sapply(list_ind_y, function(x) x[, prop.table(table(ordY))]))
  }
  
  # transform into proportion matrix 
  abundance_table <<- do.call(rbind,abundance_matrix)%>% data.frame()%>%
    rename("HC" = "X1",
           "SC" = "X2",
           "Algae" = "X3") %>%
    dplyr::mutate(year = factor(rep(1:nyear, each=nreef)))%>%
    dplyr::mutate(reef = factor(rep(paste('reef_', 1:`nreef`),nyear))
    )
  return(abundance_table)
}

# Inverse matrix
gettVinv <- function(V){
  return(ginv(t(V)))
}

