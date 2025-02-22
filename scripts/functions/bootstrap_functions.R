
block_bootstrap_pca <- function(data, block_length = 50, R = 1000, reference_pca = NULL) {
  n <- nrow(data)
  results <- list()
  
  if (is.null(reference_pca)) {
    reference_pca <- prcomp(data, scale. = FALSE)
  }
  
  for(i in 1:R) {
    blocks <- ceiling(n/block_length)
    start_points <- sample(1:(n - block_length + 1), blocks, replace = TRUE)
    
    boot_indices <- unlist(lapply(start_points, function(x) x:(x + block_length - 1)))
    boot_indices <- boot_indices[boot_indices <= n][1:n]
    
    boot_data <- data[boot_indices, ]
    boot_pca <- prcomp(boot_data, scale. = FALSE)
    
    # Align signs with reference PCA
    for(j in 1:ncol(boot_pca$rotation)) {
      cor_sign <- sign(sum(boot_pca$rotation[,j] * reference_pca$rotation[,j] ))
      if(cor_sign < 0) {
        boot_pca$rotation[,j] <- -boot_pca$rotation[,j]
        boot_pca$x[,j] <- -boot_pca$x[,j]
      }
    }
    
    results[[i]] <- list(
      rotation = boot_pca$rotation,
      sdev = boot_pca$sdev
    )
  }
  
  return(results)
}

bbootstrap_pca_ts_par <- function(data,
                             block_length = 50,
                             replicates = 1000,
                             reference_pca = NULL,
                             scale = 1,
                             eps = 0.01,
                             workers = 6) {
    start_time <- Sys.time()
    n <- nrow(data)
    data <- data * scale

    D <- ncol(data)
    basis_vectors <- lapply(1:(D - 1), generate_orthonormal_basis, D)

    if (is.null(reference_pca)) {
        reference_pca <- prcomp(clr(data), scale. = FALSE)
    }
    
    plan(multisession, workers = workers)

    results <- future_lapply(1:replicates, function(i) {
      cat("Starts Replicate:", i, "\n")
      blocks <- ceiling(n/block_length)
      start_points <- sample(1:(n - block_length + 1), blocks, replace = TRUE)
      
      boot_indices <- unlist(lapply(start_points, function(x) x:(x + block_length - 1)))
      boot_indices <- boot_indices[boot_indices <= n][1:n]
      
      boot_data <- data[boot_indices, ]
      boot_results <- co_pca_mcem_bfgs(boot_data,
                                      max_iter = 50,
                                      r = 10,
                                      lambda = 1,
                                      eps = eps)
      boot_pca <- boot_results$pca
      
      for(j in seq_len(ncol(boot_pca$rotation))) {
          cor_sign <- sign(sum(boot_pca$rotation[, j] * reference_pca$rotation[, j]))
          if (cor_sign < 0) {
              boot_pca$rotation[, j] <- -boot_pca$rotation[, j]
          }
      }
      cat("Ends Replication:", i, "\n")      
      list(
          rotation = boot_pca$rotation,
          sdev = boot_pca$sdev
      )
    }, future.seed = TRUE)
    
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    print(paste("The bootstrap finished after:", elapsed_time, attr(elapsed_time, "units")))
    return(results)
}
