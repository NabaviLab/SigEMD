# SigEMD needs log2(TPM+1) as input data.
# Example: Prepare data, row represents gene, column represents sample
data <- matrix(rnorm(10000), nrow=100, ncol=100)
data[data<0] <- 0
rownames(data) <- paste("gene", 1:100, sep="")
colnames(data) <- paste("sample", 1:100, sep="")

# Prepare the groups of the samples
condition<- c(rep("A",50),rep("B",50))
names(condition)<-colnames(data)


# Call calculate_single. 
#Here we only use 5 permutations as an example, but in actual experiments using at least 100 permutations is advised. 
results<- calculate_single(data, condition, binSize=0.2,nperm=5)
# The function returns a matrix with EMD score and the adjust p-value of each gene.
emd<- results$emdall
head(emd)

# The function plot_emd_density_sig will display the density distributions of each of the groups for a given gene.
plot_emd_density_sig(results,"gene1")

# We can plot the gene with the largest EMD score:
plot_emd_density_sig(results,rownames(emd[order(-emd[,"emd"]),])[1])
