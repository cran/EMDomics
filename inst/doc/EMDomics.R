## ------------------------------------------------------------------------
exp_data <- rnorm(100)
names(exp_data) <- paste("sample", 1:100)

groupA <- names(exp_data)[1:50]
groupB <- names(exp_data)[51:100]

## ------------------------------------------------------------------------
library(ggplot2)
df <- as.data.frame(exp_data)
df$group[1:50] <- "A"
df$group[51:100] <- "B"
ggplot(df, aes(exp_data, fill=group)) + geom_density(alpha=0.5)

## ------------------------------------------------------------------------
library(EMDomics)
calculate_emd_gene(exp_data, groupA, groupB)

## ------------------------------------------------------------------------
exp_data2 <- exp_data
mod_vec <- sample(c(2,-2), 50, replace=TRUE)
exp_data2[1:50] <- exp_data2[1:50] + mod_vec

## ------------------------------------------------------------------------
df <- as.data.frame(exp_data2)
df$group[1:50] <- "A"
df$group[51:100] <- "B"
ggplot(df, aes(exp_data2, fill=group)) + geom_density(alpha=0.5)

calculate_emd_gene(exp_data2, groupA, groupB)

## ------------------------------------------------------------------------
data <- matrix(rnorm(10000), nrow=100, ncol=100)
rownames(data) <- paste("gene", 1:100, sep="")
colnames(data) <- paste("sample", 1:100, sep="")

groupA <- colnames(data)[1:50]
groupB <- colnames(data)[51:100]

## ----, message=FALSE-----------------------------------------------------
results <- calculate_emd(data, groupA, groupB, nperm=10)

## ------------------------------------------------------------------------
emd <- results$emd
head(emd)

## ------------------------------------------------------------------------
emd2 <- emd[(order(emd[,"q-value"])),]
head(emd2)

## ------------------------------------------------------------------------
emd3 <- emd[(order(emd[,"emd"])),]
smallest_gene <- rownames(emd3)[1]
biggest_gene <- rownames(emd3)[nrow(emd3)]

plot_density(results, smallest_gene)
plot_density(results, biggest_gene)

## ----, message=FALSE-----------------------------------------------------
plot_perms(results)

## ------------------------------------------------------------------------
plot_emdnull(results)

