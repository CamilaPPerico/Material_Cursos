library(rSWeeP)
setwd('~/UFPR/EVENTOS/EPB2024/minicursoSWeeP/')
# vetorize as sequências de SARS-CoV-2
sw = SWeePlite('SARSCoV2data/',seqtype='AA',psz=600,mask=c(2,1,2),bin=FALSE,verbose=T)

# Veja os dados
sw$info

# carregue os metadados
mt = read.csv('metadata_SARSCoV2.csv')

# visaualização com PCA
pca_output = prcomp (sw$proj , scale = FALSE)
par(mfrow=c(1,2))
plot(pca_output$x[,1] , pca_output$x [,2] , xlab = 'PC-1' , ylab = 'PC-2' , pch =20 , col = mt$id)
legend("topleft" , unique(mt$variant) , col = as.character(c(1:max( mt$id))) , pch =20)
plot(pca_output$x[,3] , pca_output$x [,4] , xlab = 'PC-3' , ylab = 'PC-4' , pch =20 , col =  mt$id)

# visaualização com tSNE
library(Rtsne)
# caso haja duplicações nos dados
idx=which(duplicated(sw$proj))
tsne_output <- Rtsne::Rtsne (sw$proj[-idx,] , dims =2 , pca = FALSE , perplexity =5 , dist = euclidean)
plot(tsne_output$Y[,1] , tsne_output$Y [,2] , xlab = 'tSNE-1' , ylab = 'tSNE-2' , pch =20,col=mt$id[-idx])
legend("topleft" , unique(mt$variant) , col = as.character(c(1:max(mt$id))) , pch =20)

# visaualização com UMAP
library(umap)
umap_output <- umap(sw$proj)
plot(umap_output$layout[,1],umap_output$layout[,2] , xlab = 'UMAP-1' , ylab = 'UMAP-2' , pch =20, col=mt$id)
legend("topleft" , unique(mt$variant) , col = as.character(c(1:max(mt$id))) , pch =20)

# Árvore filogenética
library(ape)
library(ggtree)
mdist = dist(sw$proj,method='euclidean')
tr = nj(mdist)
p<- ggtree(tr)+ geom_tiplab()
p



# Heatmap
data = pca_output$x[,1:5] 
rownames(data) = mt$variant
mdist = dist(data,method='euclidean')
heatmap(as.matrix(mdist),symm=T,col=colorRampPalette(c('#000000','#FFFFFF'))(10))


# Clusterização
# Escolho o númeoro clusters como 5 (valor arbitrário)
library(stats)
kclass = kmeans(sw$proj,5,iter.max=10,nstart = 1)

par(mfrow=c(1,2))
plot(pca_output$x[,1] , pca_output$x [,2] , xlab = 'PC-1' , ylab = 'PC-2' , pch =20 , col = mt$id,main="variantes")
legend("topleft" , unique(mt$variant) , col = as.character(c(1:max( mt$id))) , pch =20)
plot(pca_output$x[,1] , pca_output$x [,2] , xlab = 'PC-1' , ylab = 'PC-2' , pch =20 , col = kclass$cluster,main="clusters")
legend("topleft" , as.character(unique(kclass$cluster)) , col = as.character(c(1:max( kclass$cluster))) , pch =20)

