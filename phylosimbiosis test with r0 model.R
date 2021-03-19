#########################################################################
Load the libraries
#########################################################################
library(ape)
library(vegan)
library(paco)
library(plyr)
library(phangorn)
#########################################################################
Read the Unifrac weighted distance table among species, convert to matrix 
anf build the NJ tree
#########################################################################
unw_phylum=read.table('weipacotodas.txt', sep = "\t", head = T, row=1)
betaunphy=as.matrix(unw_phylum)
treeNJunphy=nj(betaunphy)
#########################################################################
Read the association matrix and the complete phylogeny of Herichthys
HP=as.matrix(read.table("herlinks.txt", header=TRUE))
NLinks = sum(HP)
TreeH <- read.tree("peces.phy")
#########################################################################
Estimate the cophenetic matrices
#########################################################################
gdist=cophenetic(TreeH)
bdist=cophenetic(treeNJunphy)
#########################################################################
Prepare the data for PACo test
#########################################################################
preparados=prepare_paco_data(gdist, bdist,HP )
componentes= add_pcoord(preparados, correction='cailliez')
#########################################################################
Perform the PAcO test and print the results
#########################################################################
ajuste=PACo(componentes,nperm = 1000,seed = 25,method = "r0",symmetric = FALSE,proc.warnings = TRUE,shuffled = FALSE)
print(ajuste$gof)
#########################################################################
Perform the Parafit test of Legendre
#########################################################################
Legendre=parafit(gdist, bdist, HP, nperm = 999, test.links = TRUE,seed = NULL, correction = "cailliez", silent = FALSE)
#########################################################################
Now the same procedure for the Hericthys bartoni species group
#########################################################################
unw_phylum1=read.table('weipacobartoni.txt', sep = "\t", head = T, row=1)
betaunphy1=as.matrix(unw_phylum1)
treeNJunphy1=nj(betaunphy1)
HP1=as.matrix(read.table("pame.txt", header=TRUE))
NLinks = sum(HP1)
TreeH1 <- read.tree("bartoni.phy")
gdist1=cophenetic(TreeH1)
bdist1=cophenetic(treeNJunphy1)
preparados1=prepare_paco_data(gdist1, bdist1,HP1)
componentes1= add_pcoord(preparados1, correction='cailliez')
ajuste1=PACo(componentes1,nperm = 1000,seed = 25,method = "r0",symmetric = FALSE,proc.warnings = TRUE,shuffled = FALSE)
print(ajuste1$gof)
Legendre2=parafit(gdist1, bdist1, HP1, nperm = 999, test.links = TRUE,seed = NULL, correction = "cailliez", silent = FALSE)
#########################################################################
And finally for the group of species included in the H. cyanoguttatus species group
#########################################################################
unw_phylum2=read.table('weipacocyano.txt', sep = "\t", head = T, row=1)
betaunphy2=as.matrix(unw_phylum2)
treeNJunphy2=nj(betaunphy2)
HP2=as.matrix(read.table("tepehua.txt", header=TRUE))
NLinks = sum(HP2)
TreeH2 <- read.tree("cyano.phy")
gdist2=cophenetic(TreeH2)
bdist2=cophenetic(treeNJunphy2)
preparados2=prepare_paco_data(gdist2, bdist2,HP2 )
componentes2= add_pcoord(preparados2, correction='cailliez')
ajuste2=PACo(componentes2,nperm = 1000,seed = 25,method = "r0",symmetric = FALSE,proc.warnings = TRUE,shuffled = FALSE)
print(ajuste2$gof)
Legendre3=parafit(gdist2, bdist2, HP2, nperm = 999, test.links = TRUE,seed = NULL, correction = "cailliez", silent = FALSE)
#########################################################################





