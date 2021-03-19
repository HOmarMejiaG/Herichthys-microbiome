#########################################################################
Read the libraries
########################################################################
library(ape)
library(vegan)
library(paco)
library(plyr)
library(phangorn)
#########################################################################
Read the association matrix (HP) and the phylogenetic trees of fishes (TreeH) and microbiome (TreeP)

HP=as.matrix(read.table("generosp.txt", header=TRUE, check.names = FALSE))
NLinks = sum(HP)
TreeH <- read.tree("peces.phy")
TreeP <- read.tree("generon.nwk")
##########################################################################
Generate the cophenetic matrices
###########################################################################
gdist=cophenetic(TreeH)
ldist =cophenetic(TreeP)
#########################################################################
Prepare the data for PACo test 
#######################################################################
preparados=prepare_paco_data(gdist, ldist,HP)HP2 )
componentes= add_pcoord(preparados, correction='cailliez')
###################################################################
Perform the PACo test, print and save the results
#############################################################################
ajuste=PACo(componentes,nperm = 1000,seed = 25,method = "quasiswap",symmetric = FALSE,proc.warnings = TRUE,shuffled = FALSE)
print(ajuste$gof)
restodos=write.csv(ajuste$gof, file="todospaco.txt")
###############################################################
Obtain the residuals of the analysis and save
###############################################################
residuales=residuals_paco(ajuste$proc)
print(ajuste$proc)
###############################################################
salida=write.table(residuales, "salida.txt", row.names = TRUE, col.names = TRUE)
###############################################################
Elaborate the residuals graphic in ggplot2
###############################################################
load the libraries and theme
library(ggplot2)
theme_set(
  theme_classic() 
)
###############################################################
Read the residuals table for species and trophic guilds
###############################################################
datos=read.table("residualesespecie.txt", header=TRUE, sep= "\t",dec = ".")
###############################################################
Convert to adataframe
###############################################################
datos1=as.data.frame(datos)
###############################################################
Convert the species variable to a factor and order the species
datos1$Species <- as.factor(datos1$Species)
datos1$Species <- factor(datos1$Species, levels = c(" bartoni", " labridens", " pame", " steindachneri", " pantostictus",  " carpintis"," cyanoguttatus"," deppii", " minckleyi"," tamasopoensis"," tepehua"))
###############################################################
Make the graph and plot
###############################################################
residuales=ggplot(datos1, aes(x=Species,y=residuals, color=Species)) + geom_boxplot()
plot(residuales)
###############################################################
Now the graph for the ten most abundant phyla
###############################################################
datos2=read.table("residualesphylum.txt", header=TRUE, sep= "\t",dec = ".")
datos3=as.data.frame(datos2)
datos2$Phylum <- as.factor(datos2$Phylum)
atos2$Phylum=as.factor(datos2$Phylum)
residuales1=ggplot(datos2, aes(x=Phylum,y=residuals, color=Phylum)) + geom_boxplot()
plot(residuales1)
###############################################################





