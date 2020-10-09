#!/usr/bin/env R

# Low variance in autosomal DNAm analysis in 7 noncancer tissues

# This script describes low variance analysis on the ANOVA-filtered
# probe data from 7 tissues. First, run the script `nct_anova_bc.R`
# from a directory containing the `cgall` and `anova_results` folders
# with tissue-level ANOVA results and cg Beta-values summaries, which 
# generates `data.table.filt` files for each tissue. Next, run this 
# script to identify probes with cross-tissue low variances and 
# within-tissue high variances.

#--------------------
# summaries by tissue
#--------------------
ptid <- get(load())
dat <- matrix(nrow = 0, ncol = 9)
for(t in names(lfilt)){
  lft <- lfilt[[t]]
  lft <- lft[lft[,1] %in% ptid,]
  dt <- c(t, mean(lft$var), var(lft$var), mean(lft$mean), var(lft$mean),
          min(lft$var), max(lft$var), min(lft$mean), max(lft$mean))
  dat <- rbind(dat, dt)
}
colnames(dat) <- c("tissue", "mean_of_var", "var_of_var", 
                   "mean_of_mean", "var_of_mean", "min_var", "max_var",
                   "min_mean", "max_mean")
dat <- as.data.frame(dat, stringsAsFactors = FALSE)
for(c in 2:ncol(dat)){dat[,c] <- as.numeric(dat[,c])}

summary(dat$var_of_mean)
summary(dat$mean_of_mean)

min(dat$min_var)
max(dat$min_var)

min(dat$min_mean)
max(dat$max_mean)

#----------------------
# summaries by probe id
#----------------------
ldat.mean <- ldat.var <- list()
for(id in ptid){
  nr.mean <- nr.var <- c()
  for(t in names(lfilt)){
    lft <- lfilt[[t]]
    lft <- lft[lft$cgid == id,]
    nr.mean <- c(nr.mean, lft$mean)
    nr.var <- c(nr.var, lft$var)
  }
  ldat.mean[[id]] <- nr.mean
  ldat.var[[id]] <- nr.var
  message(id)
}

did.mean <- do.call(rbind, ldat.mean)
did.var <- do.call(rbind, ldat.var)
colnames(did.mean) <- colnames(did.var) <- names(lfilt)
for(c in seq(ncol(did.mean))){did.mean[,c] <- as.numeric(did.mean[,c])}
for(c in seq(ncol(did.var))){did.var[,c] <- as.numeric(did.var[,c])}

mean.min <- as.numeric(apply(did.mean, 1, min))
mean.max <- as.numeric(apply(did.mean, 1, max))
mean.range <- mean.max - mean.min

which(mean.range > 0.1)

#---------------------------------
# mean-filtered hypovar probe list
#---------------------------------
ptid.fn <- "cgidfilt_txmean-rangefilt01.rda"
meanrange.max <- 0.01
cgid.filt <- which(mean.range < meanrange.max)
ptidf <- ptid[cgid.filt]
length(ptidf) # 4577 # 4577/5188 = 0.8822282, 88%
save(ptid.fn, file = ptid.fn)

#-----------------
# mapping analysis
#-----------------
sta$is.prom <- ifelse(grepl("TSS|5'", sta$UCSC_RefGene_Group), TRUE, FALSE)
sta$is.isl <- ifelse(sta$Relation_to_Island == "Island", TRUE, FALSE)
table(sta$is.prom, sta$is.isl)
3107/nrow(sta) # 0.6788289
