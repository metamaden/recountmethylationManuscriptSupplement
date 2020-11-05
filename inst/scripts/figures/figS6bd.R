#!/usr/bin/env R

# Get supplemental plots of DNAm variances between PCA comparator 
# groups (cancers vs. non-cancers).

library(ggplot2)

#----------
# load data
#----------
pkg.name <- "recountmethylationManuscriptSupplement"
pca.dir <- system.file("extdata", "pcadata", package = pkgname) 
mv1 <- get(load(file.path(pca.dir, "m-var_blood-vs-leukemia.rda")))
mv2 <- get(load(file.path(pca.dir, "m-var_brain-vs-braintumor.rda")))

#------------------
# make violin plots
#------------------
# blood versus leukemias
vp <- data.frame(var = c(mv1[,1], mv1[,2]), stringsAsFactors = FALSE)
vp$label <- rep(c("blood", "leukemia"), each = nrow(mv1))

figS6b <- ggplot(vp, aes(x = label, y = var, fill = label)) + 
  geom_violin(show.legend = FALSE) +
  scale_fill_manual(values = c("red", "purple")) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90),
                     text = element_text(size=20)) +
  xlab("Tissue") + ylab("Variance")

# brain versus brain tumors
vp <- data.frame(var = c(mv2[,1], mv2[,2]), stringsAsFactors = FALSE)
vp$label <- rep(c("brain", "braintumor"), each = nrow(mv2))
figS6d <- ggplot(vp, aes(x = label, y = var, fill = label)) + 
  geom_violin(show.legend = FALSE) +
  scale_fill_manual(values = c("blue", "cyan4")) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90),
                     text = element_text(size=20)) +
  xlab("Tissue") + ylab("Variance")

#-----------
# save plots
#-----------
#pdf("sfig_dnam-var_blood-vs-leukemia.pdf", 3, 4)
#figS6b
#dev.off()

#pdf("sfig_dnam-var_brain-vs-braintumor.pdf", 3, 4)
#figS6d
#dev.off()