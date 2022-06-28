###########################         Import, remove duplicates and compress NJR        ###########################

pacman::p_load(pacman, data.table, rio, tidyverse)

directory <- "A:\DATA"
setwd(directory)

## Import
# fread was quicker
# However, unfortunately, I couldn't resolve the problem of some unquoted lines with the csv files, and - as such - some records were lost which was not acceptable

rk <- import("A:\DATA\njr\KneeAllRevisions.txt", quote="")
rk <- as.data.table(rk)

pk <- import("A:\DATA\njr\KneePrimaryOutcomes.txt", quote="")
pk <- as.data.table(pk)

## Change all names to lower case
rk <- rk %>% janitor::clean_names()
pk <- pk %>% janitor::clean_names()

## Sort out dates
rk[, op_date := as.Date(fasttime::fastPOSIXct(op_date))]
pk[, primary_op_date := as.Date(fasttime::fastPOSIXct(primary_op_date))]
pk[, revision_date := as.Date(fasttime::fastPOSIXct(revision_date))]

## Drop 'cat_no', 'details' and 'text' columns (by setting them to NULL)
rk[, grep("details", names(rk)) := NULL]
rk[, grep("text", names(rk)) := NULL]
pk[, grep("details", names(pk)) := NULL]
pk[, grep("cat_no", names(pk)) := NULL]
pk[, grep("text", names(pk)) := NULL]

## Replace "NULL", "" or "na" with NAs 
# Need to exclude date columns from this

# Revision dataset
p_load(car)
.cols <- setdiff(colnames(rk), "op_date")
rk[,(.cols):=lapply(.SD, recode, '"NULL"=NA'), .SDcols = .cols]
rk[,(.cols):=lapply(.SD, recode, '""=NA'), .SDcols = .cols]
rk[,(.cols):=lapply(.SD, recode, '"na"=NA'), .SDcols = .cols]

# Primary dataset
.cols <- setdiff(colnames(pk), c("primary_op_date", "revision_date"))
pk[,(.cols):=lapply(.SD, recode, '"NULL"=NA'), .SDcols = .cols]
pk[,(.cols):=lapply(.SD, recode, '""=NA'), .SDcols = .cols]
pk[,(.cols):=lapply(.SD, recode, '"na"=NA'), .SDcols = .cols]

## Convert all character fields to factor

# Revision dataset
changeCols <- colnames(rk)[which(as.vector(rk[,lapply(.SD, class)]) == "character")]
rk[,(changeCols):= lapply(.SD, as.factor), .SDcols = changeCols]

# Primary dataset
changeCols <- colnames(pk)[which(as.vector(pk[,lapply(.SD, class)]) == "character")]
pk[,(changeCols):= lapply(.SD, as.factor), .SDcols = changeCols]

# Convert numeric fields
# Except for those which need to remain numeric
# Revision dataset
changeCols <- colnames(rk)[which(as.vector(rk[,lapply(.SD, class)]) == "numeric")]
changeCols <- setdiff(changeCols, "bmi")
rk[,(changeCols):= lapply(.SD, as.factor), .SDcols = changeCols]

# Primary dataset
changeCols <- colnames(pk)[which(as.vector(pk[,lapply(.SD, class)]) == "numeric")]
changeCols <- setdiff(changeCols, c("bmi", "primary_to_outcome_years"))
pk[,(changeCols):= lapply(.SD, as.factor), .SDcols = changeCols]

# Remove duplicates based on ALL variables
pkslim <- unique(pk)
rkslim <- unique(rk)

# Attrition data for flowchart
no_pkr_supp <- pk[,.N]
no_rkr_supp <- rk[,.N]
no_uniq_pkr_supp <- pkslim[,.N]
no_uniq_rkr_supp <- rkslim[,.N]

# Save the environment for quickly re-loading the compressed NJR with all fields
rm(list=setdiff(ls(), c("pkslim","rkslim", "no_pkr_supp", "no_rkr_supp", "no_uniq_pkr_supp", "no_uniq_rkr_supp")))

save.image("njr.RData")
save.image("njr backup.RData")