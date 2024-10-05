# 1. SweeP dos proteomas mitocondriais
# Considere um conjunto de 13 proteomas mitocondriais (CDSs traduzidos) 
# depositados na pasta no endereço `path` em formato FASTA. 
# Estas sequências podem ser vectorizadas
# **arquivos inclusos no pacote rSWeeP instalado

library(rSWeeP)
path = paste (system.file("examples/aaMitochondrial/",package = "rSWeeP"),'/', sep = '')
sw = SWeePlite(path,seqtype='AA',mask=c(4),psz=1000)

# veja as informações
sw$info
# veja os vetores
sw$proj

# 2. Visualização gráfica
# Para visualizar os dados vectorizados, fornecemos os metadados com as classes 
# ao nível taxonómico da família.
# metadados também inclusos no pacote
pathmetadata <- system.file(package = "rSWeeP" , "examples" , "metadata_mitochondrial.csv")
mt = read.csv(pathmetadata,header=TRUE)


vunq = unique(mt$family) # lista de famíliar únicas
ncl = length(vunq) # número de  familias
cl = mt$id         # classes numericas


## 2.1. Principal Component Analysis (PCA)
pca_output <- prcomp (sw$proj , scale = FALSE)
par(mfrow=c(1,2))
plot(pca_output$x[,1] , pca_output$x [,2] , xlab = 'PC-1' , ylab = 'PC-2' , pch =20 , col = cl)
legend("bottomright" , vunq , col = as.character(c(1:ncl)) , pch =20)
plot(pca_output$x[,3] , pca_output$x [,4] , xlab = 'PC-3' , ylab = 'PC-4' , pch =20 , col = cl)



## 2.2. tSNE
library(Rtsne)
tsne_output <- Rtsne::Rtsne (sw$proj , dims =2 , pca = FALSE , perplexity =2 , dist = euclidean)
plot(tsne_output$Y[,1] , tsne_output$Y [,2] , xlab = 'tSNE-1' , ylab = 'tSNE-2' , pch =20,col=cl)
legend("topleft" , vunq , col = as.character(c(1:ncl)) , pch =20)


## 2.3. Heatmap
#O heatmap torna visíveis os grupos relacionados no conjunto de dados, 
# com base na matriz de distância euclidiana. 
# Estes grupos coincidem com as famílias taxonómicas. 
# A cor preta indica maior semelhança e a cor branca nenhuma semelhança.
data = sw$proj
rownames(data) = mt$family
mdist = dist(data,method='euclidean')
heatmap(as.matrix(mdist),symm=T,col=colorRampPalette(c('#000000','#FFFFFF'))(10))

## 2.4. Relações filogenéticas
# Também podemos obter a relação filogenética entre os organismos vectorizados
#  utilizando o método Neighbour Joining (NJ).
library(ape)
library(ggtree)
# matriz de distâncias
mdist = dist(sw$proj,method='euclidean')
# roda o agortimo de NJ
tr = nj(mdist) 
# enraiza a árvore no grupo externo
tr = root(tr,outgroup='14_Rhazya_stricta') 
# plot da árvore
p<- ggtree(tr)+ geom_tiplab()
p

# nomeando de acordo com a família
data = sw$proj
rownames(data) = mt$family # renomeio os vetores
mdist = dist(data,method='euclidean')
tr = nj(mdist)
tr = root(tr,outgroup='Apocynaceae')
p<- ggtree(tr)+ geom_tiplab()
p
