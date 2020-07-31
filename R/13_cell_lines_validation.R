# pre processing of validation cell lines

## read gated files and build a labeled DF


# Set working directory
setwd("~/Documents/Massion_lab/manuscripts/prelim_cytof_adc/data/cell_lines/processed_files/mix")

# read cell lines samples and build a labeled df
files_list <- list.files(pattern='.fcs|.FCS')

# Extract cell line ID
cl_ID <- sapply(strsplit(files_list, "_"), "[[", 2)

# Read first file to correct column order later
smp <- flowCore::read.FCS(files_list[1], transformation = FALSE)
descrp <- smp@parameters@data$desc
smp <- data.frame(smp@exprs)
col_nms <- colnames(smp)

big_df <- data.frame(matrix(ncol = ncol(smp)+1, nrow=0))
colnames(big_df) <- c(col_nms, 'cl_ID')

for (i in 1:length(files_list)){
    # Read exprs from FCS file
    dt <- data.frame(flowCore::read.FCS(files_list[i], transformation = FALSE)@exprs)

    # Order columns
    dt <- dt[col_nms]

    dt['cl_ID'] <- rep(cl_ID[i], nrow(dt))

    big_df <- rbind(big_df, dt)
}

colnames(big_df)<- c(descrp, 'cl_ID')

# Save RData file
save(big_df, file = paste0('mix_celllines_labeled', '.RData'))