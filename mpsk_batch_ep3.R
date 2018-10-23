#!/usr/bin/env Rscript

########## SET ENVIRONMENT #########

### library tools
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(reticulate)))

### plots
#suppressWarnings(suppressMessages(library(gridExtra)))
#suppressWarnings(suppressMessages(library(tsne)))

### Main model
suppressWarnings(suppressMessages(library(MPSK)))

library("optparse")

### directory path
datadir = "/data/flow/EP3/files"

########## PASSING ARGUMENTS #########
# https://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/

### set options
option_list = list(
    
    make_option(c("-o", "--out"), 
                type    = "character", 
                default = "out", 
                help    = "header of output file name [default= %default]", 
                metavar = "character"),
    
    make_option(c("-t", "--trt"), 
                type    = "character", 
                default = "01", 
                help    = "which treatment (01 = 02 = 03 | 05 = 06 = 07 | 09 = 10 = 11); [default= %default]"),
    
    make_option(c("-p", "--patient"), 
                type    = "character", 
                default = "F6901PRY_C", 
                help    = "one of patients (F6901PRY_C, G6904VJT_E, or K6902C85_A); [default= %default]"),
    
    make_option(c("-n", "--number"), 
                type    = "integer", 
                default = 100, 
                help    = "number of events [default= %default]", 
                metavar = "integer"),
    
    make_option(c("-s", "--seed"), 
                type    = "integer", 
                default = 0, 
                help    = "number of seeds [default= %default]", 
                metavar = "integer")
); 
 
### parse
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);
cat("\n================================================\n")
cat("Parsing Arguments\n\n")
print(opt)

########## IMPORT DATA #########
cat("\n================================================\n")
cat("Import Data\n\n")

### filename criteria
trt     = paste0("[ACE]", opt$trt) #"[ACE]06"
patient = opt$patient              #"F6901PRY_C"

### select files
fnames = dir(datadir) %>% grep(trt, ., value = TRUE) %>% grep(patient, ., value = TRUE)
fnames = sort(fnames)

cat("File names\n")
print(fnames)

### import data --- marker names
markers = c('FSC_A', 'SSC_A', 'CD3', 'CD4', 'CD8', 'IFN+IL2')

### import data --- EP3 flow data
# initiate numpy downloader
np = import('numpy')

# read each npy file
dat_raw = sapply(fnames, function(fn, col_names){
    mat = np$load(file.path(datadir, fn))
    colnames(mat) = col_names
    return(mat)
}, col_names = markers)


########## SELECT SUBSET OF EACH MATRIX AND COMBINE #########
cat("\n================================================\n")
cat("Select subset of each matrix and combine\n\n")

# intialization
set.seed(opt$seed)
N = opt$number

# select a proportion of cells for each flow experiment at random
dat_sub = lapply(dat_raw, function(mat){
    idx = sample(1:nrow(mat), N)
    return(mat[idx, ])
})

# combine data points by row
Y_ref = do.call(rbind, dat_sub[fnames])

# set labels for each sample
C_ref = c()
for (idx in 1:length(fnames)){
    C_ref <- c(C_ref, rep(idx, N))
}

######### APPLY MPSK ##########
cat("\n================================================\n")
cat("Apply MPSK\n\n")

### Check data before proceed
cat("Check data before proceed\n")

if (nrow(Y_ref) != length(C_ref)){
    print(nrow(Y_ref))
    print(length(C_ref))
    stop("Number of observations does not match length of labels", call.=FALSE)
}
cat("   Pass!\n")

if (nrow(Y_ref) != (N * length(fnames))){
    print(nrow(Y_ref))
    print(N * length(fnames))
    stop("Number of observations does not match (number of events * number of files)", call.=FALSE)
}
cat("   Pass!\n")

### Summarize data 
cat("\n")
cat("Summarize data before applying MPSK model:\n")
cat("- Patient:                         ", opt$patient,    "\n")
cat("- Treatment:                       ", opt$trt,        "\n")
cat("- Number of samples:               ", length(fnames), "\n")
cat("- Number of cells from each sample:", N,              "\n")
cat("- Dimension of Data Matrix:        ", dim(Y_ref),     "\n")
cat("- Legnth of labels:                ", length(C_ref),  "\n")


### Apply MPSK
cat("\n")
cat("Standardization\n")
Y   <- scale(Y_ref, center = T, scale = T)

cat("Apply MPSK model\n")

tt_start = Sys.time()
res      = mpsk(Y, C_ref)
tt_stop  = Sys.time()

out = list(time_start = tt_start,
           time_stop  = tt_stop,
           res        = res)

saveRDS(out, paste0(opt$out, ".RDS"))
print(difftime(tt_stop, tt_start))
#res           = mpsk(Y, C_ref)
#resRelab      = relabelChain(res)
#resCalibrated = calibrate(resRelab)
#chainSummary  = summarizeChain(resRelab)
