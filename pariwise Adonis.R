#########################################################################
Load the libraries and convert the matrix to a data frame
#########################################################################
library(EcolUtils)
library(vegan)
library(pairwiseAdonis)
setwd("C:/Users/homar/OneDrive/Escritorio/Microbioma/paco")
peces= read.delim('peces.txt',sep='\t',row.names=1,header=T) 
pecesm=as.data.frame(peces)
#########################################################################
Estimate the pairwise comparision matrix among species from Bray Curtis and Bonferroni correction
#########################################################################
Especieb= pairwise.adonis(x=pecesm[,3:4898],factors=pecesm$Especie,sim.function='vegdist',
                sim.method='bray',p.adjust.m='bonferroni')
#########################################################################
Same procedure as above but among trophic guilds
#########################################################################
Grupob=pairwise.adonis(x=pecesm[,3:4898],factors=pecesm$Grupo,sim.function='vegdist',
                      sim.method='bray',p.adjust.m='bonferroni')
#########################################################################
Export the results to tables
#########################################################################
write.table(Especieb, file="Especieb.txt", row.names=TRUE, col.names=TRUE)
write.table(Grupob, file="Grupob.txt", row.names=TRUE, col.names=TRUE)
#########################################################################
Estimate pairwise comparison using Unifrac distance
#########################################################################
Read the distance matrix of weighted Unifrac
#########################################################################
dist_uni=read.delim("weiadonis mod.txt",sep='\t',row.names=1,header=T)
df_dist_uni=as.data.frame(dist_uni)
#########################################################################
Read the table with the labels
#########################################################################
factores=read.delim("etiquetas.txt",sep='\t',row.names=1,header=T)
facdt=as.data.frame(factores)
especie=as.vector(unlist(facdt$Especie))
grupo=as.vector(unlist(facdt$Grupo))
#########################################################################
Perform the pairwise comparison among species and trophic guilds
#########################################################################
especie=pairwise.adonis(df_dist_uni,factors=especie,p.adjust.m='bonferroni')
grupo=pairwise.adonis(df_dist_uni,factors=grupo,p.adjust.m='bonferroni')
#########################################################################
Export the results to tables
#########################################################################
write.table(especie, file="especie weighted.txt", row.names=TRUE, col.names=TRUE)
write.table(grupo, file="grupo weighted.txt", row.names=TRUE, col.names=TRUE)
#########################################################################
Now the same procedure wit the Unifrac Unweighted distance
#########################################################################
dist_uni=read.delim("unwadonis mod.txt",sep='\t',row.names=1,header=T)
df_dist_uni=as.data.frame(dist_uni)
factores=read.delim("etiquetas.txt",sep='\t',row.names=1,header=T)
facdt=as.data.frame(factores)
especie=as.vector(unlist(facdt$Especie))
grupo=as.vector(unlist(facdt$Grupo))
especie=pairwise.adonis(df_dist_uni,factors=especie,p.adjust.m='bonferroni')
grupo=pairwise.adonis(df_dist_uni,factors=grupo,p.adjust.m='bonferroni')
write.table(especie, file="especie unweighted.txt", row.names=TRUE, col.names=TRUE)
write.table(grupo, file="grupo unweighted.txt", row.names=TRUE, col.names=TRUE)
#########################################################################


