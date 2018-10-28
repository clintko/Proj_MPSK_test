### load library
suppressWarnings(suppressMessages(library("tidyverse")))
suppressWarnings(suppressMessages(library("MPSK")))
suppressWarnings(suppressMessages(library("mvtnorm")))

cat("Goal: run MPSK on small multivariate normal dataset\n")
cat("in order to see if MPSK is able to handle large dataset\n")

### function that generate a sample with clusters size = n and centered at mu
generate_sample <- function(mu, n){
    Y = lapply(
        mu, 
        function(x){
            rmvnorm(n, x + c(rnorm(1), rnorm(1)), diag(c(1, 1)))
    }) # end lapply
    
    Y = do.call(rbind, Y)
    Z = rep(1:length(mu), each = n)
    mat = cbind(Y, Z)
    colnames(mat) = c("x", "y", "z")
    
    return(mat)
} # end func

### set parameters
K  = 12   # number of samples
mu = list(
    c( 0,  0),
    c( 0,  5),
    c( 5,  0),
    c( 5,  5))

### generate data
set.seed(0)
N   = 100
dat = replicate(K, generate_sample(mu, N))
dat = lapply(1:K, function(idx){dat[,,idx]})
dat = do.call(rbind, dat)

Y_small   = dat[, c(1, 2)]
C_small   = rep(1:K, each = N * length(mu))

### standardization
cat("\nStandardization\n")
Y_small_scaled <- scale(Y_small, center = T, scale = T)

### store the data
write.table(Y_small,        "Y_small_raw.txt",    sep = "\t", row.names = FALSE)
write.table(Y_small_scaled, "Y_small_scaled.txt", sep = "\t", row.names = FALSE)
write.table(C_small,        "C_small_raw.txt",    sep = "\t", row.names = FALSE)

### run MPSK model
cat("\nApply MPSK model\n")

cat("===== run small dataset =====\n")
tt_start = Sys.time()
res      = mpsk(Y_small_scaled, C_small)
tt_stop  = Sys.time()

out = list(
    time_start = tt_start,
    time_stop  = tt_stop,
    res        = res)

saveRDS(out, paste0("test_mpsk_blob_small", ".RDS"))
cat(difftime(tt_stop, tt_start), "\n")
