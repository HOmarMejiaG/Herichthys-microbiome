#########################################################################
Load the files
#########################################################################
matriz <- read.delim('peces_table.from_biom_sinpame8.txt.csv',sep='\t',row.names=1,header=T)
OTU = otu_table(matriz, taxa_are_rows=TRUE)
suelos = phyloseq(OTU)
read.tree ("rooted_peces_tree.nwk")->soil_tree
metadata=read.table("peces_sample_data_bueno3.3_sin_pame8_2.csv", header=T)
rownames(metadata)<-sample_names(suelos)
sampledata=sample_data(metadata)
taximat = as.matrix(read.table("Taxonomy_table__peces_editada_16s.tsv.csv", header=T, row.names=1, sep="\t", na.strings="NA" ))
taxi=tax_table(taximat)
suelos=phyloseq(OTU,taxi,soil_tree, sampledata)
suelos
#########################################################################
Estimate Alpha diversity
#########################################################################
alpha_meas = c("Observed", "Shannon", "Simpson", "InvSimpson")
(p <- plot_richness(suelos, "Especie", measures=alpha_meas))
p + geom_boxplot(data=p$data, aes(x=Especie, y=value, color=NULL), alpha=0.1)
#########################################################################
Estimate distance matrices
#########################################################################
Unweighted unifrac
unifrac_dist = phyloseq::distance(suelos, method="unifrac", weighted=F)
ordination = ordinate(suelos, method="PCoA", distance=unifrac_dist)
plot_ordination(suelos, ordination, color="Especie") + theme(aspect.ratio=1)
adonis(unifrac_dist ~ sample_data(suelos)$Especie)
#########################################################################
Bray-Curtis
#########################################################################
bray_dist = phyloseq::distance(suelos, method="bray")
ordination = ordinate(suelos, method="PCoA", distance=bray_dist)
plot_ordination(suelos, ordination, color="Especie") + theme(aspect.ratio=1)
adonis(bray_dist ~ sample_data(suelos)$Especie)
#########################################################################
Weighted unifrac
wunifrac_dist = phyloseq::distance(suelos, method="unifrac", weighted=T)
ordination = ordinate(suelos, method="PCoA", distance=wunifrac_dist)
plot_ordination(suelos, ordination, color="Especie") + theme(aspect.ratio=1)
adonis(wunifrac_dist ~ sample_data(suelos)$Especie)
#########################################################################
Perform the PERMANOVA analysis
#########################################################################
meta = metadata[order(row.names(metadata)),]
#BC.dist<- phyloseq::distance(suelos,method = "unifrac", weighted=TRUE, normalized=TRUE)
BC.dist<- phyloseq::distance(suelos,method = "unifrac", weighted=FALSE, normalized=TRUE)
BC.dist2<- phyloseq::distance(suelos,method = "bray", weighted=TRUE, normalized=TRUE)
adonis(BC.dist ~ Especie, data = meta, permutations = 9999)
adonis(BC.dist ~ Grupo_trofico, data = meta, permutations = 9999)
adonis(BC.dist2 ~ Especie, data = meta, permutations = 9999)
adonis(BC.dist2 ~ Grupo_trofico, data = meta, permutations = 9999)
#########################################################################
Make the comparison between sympatric species in Rio Verde
#########################################################################
Unweighted UniFrac
#########################################################################
unifrac_dist = phyloseq::distance(suelos_barvslabri, method="unifrac", weighted=F)
ordination = ordinate(suelos_barvslabri, method="PCoA", distance=unifrac_dist)
plot_ordination(suelos_barvslabri, ordination, color="Especie") + theme(aspect.ratio=1) + theme_bw()
plot_ordination(suelos_barvslabri, ordination, color="Grupo_trofico") + theme(aspect.ratio=1) + theme_bw()
adonis(unifrac_dist ~ sample_data(suelos_barvslabri)$Especie)
adonis(unifrac_dist ~ sample_data(suelos_barvslabri)$Grupo_trofico)
#########################################################################
Weighted UniFrac
#########################################################################
unifrac_dist = phyloseq::distance(suelos_barvslabri, method="unifrac", weighted=T)
ordination = ordinate(suelos_barvslabri, method="PCoA", distance=unifrac_dist)
plot_ordination(suelos_barvslabri, ordination, color="Especie") + theme(aspect.ratio=1) + theme_bw()
adonis(unifrac_dist ~ sample_data(suelos_barvslabri)$Especie)
adonis(unifrac_dist ~ sample_data(suelos_barvslabri)$Grupo_trofico)
#########################################################################
Bray-Curtis
#########################################################################
bray_dist = phyloseq::distance(suelos_barvslabri, method="bray")
ordination = ordinate(suelos_barvslabri, method="PCoA", distance=bray_dist)
plot_ordination(suelos_barvslabri, ordination, color="Especie") + theme(aspect.ratio=1) + theme_bw()
adonis(bray_dist ~ sample_data(suelos_barvslabri)$Especie)
adonis(bray_dist ~ sample_data(suelos_barvslabri)$Grupo_trofico)
#########################################################################
Make the comparison among sympatric species in Rio Gallinas
#########################################################################
Unweighted
#########################################################################
unifrac_dist = phyloseq::distance(suelos_tam_pam_stein, method="unifrac", weighted=F)
ordination = ordinate(suelos_tam_pam_stein, method="PCoA", distance=unifrac_dist)
plot_ordination(suelos_tam_pam_stein, ordination, color="Especie") + theme(aspect.ratio=1) + theme_bw()
adonis(unifrac_dist ~ sample_data(suelos_tam_pam_stein)$Especie)
adonis(unifrac_dist ~ sample_data(suelos_tam_pam_stein)$Grupo_trofico)
#########################################################################
Weighted
#########################################################################
wunifrac_dist = phyloseq::distance(suelos_tam_pam_stein, method="unifrac", weighted=T)
ordination = ordinate(suelos_tam_pam_stein, method="PCoA", distance=wunifrac_dist)
plot_ordination(suelos_tam_pam_stein, ordination, color="Especie") + theme(aspect.ratio=1) + theme_bw()
adonis(wunifrac_dist ~ sample_data(suelos_tam_pam_stein)$Especie)
adonis(wunifrac_dist ~ sample_data(suelos_tam_pam_stein)$Grupo_Trofico)
#########################################################################
Bray-Curtis
#########################################################################
bray_dist = phyloseq::distance(suelos_tam_pam_stein, method="bray")
ordination = ordinate(suelos_tam_pam_stein, method="PCoA", distance=bray_dist)
plot_ordination(suelos_tam_pam_stein, ordination, color="Especie") + theme(aspect.ratio=1) + theme_bw()
adonis(bray_dist ~ sample_data(suelos_tam_pam_stein)$Especie)
adonis(bray_dist ~ sample_data(suelos_tam_pam_stein)$Grupo_trofico)
#########################################################################
Estimate the Negative Binomial
#########################################################################
bartoni vs labridens
#########################################################################
bd1 <- read.delim('Para_estadistica_peces.txt', sep='\t', row.names=1, header=T)
# Model con bd as ref
m2<-glm.nb(bd2$OTUS~bd2$Fish)
summary(m2)
is.factor(bd2$Fish)
boxplot(bd2$OTUS~bd2$Fish)
# Changing Reference Value
bd2$Sitio <- relevel(bd2$Fish, ref = "bartoni")
str(bd2)
m3<-glm.nb(bd2$OTUS~bd2$Fish)
summary(m3)
# Changing Reference Value
bd2$Fish <- relevel(bd2$Fish, ref = "labridens")
str(bd2)
m4<-glm.nb(bd2$OTUS~bd2$Fish)
summary(m4)
##Binomial negative pamvstamvsstein
bd2 <- read.delim('Para_estadistica_peces_pamvstamvsstein.txt', sep='\t', row.names=1, header=T)
bd2
m2<-glm.nb(bd2$OTUS~bd2$Fish)
summary(m2)
is.factor(bd2$Fish)
boxplot(bd2$OTUS~bd2$Fish)
#########################################################################
pame vs steindachneri
#########################################################################
# Changing reference
bd2$Sitio <- relevel(bd2$Fish, ref = "pame")
str(bd2)
m3<-glm.nb(bd2$OTUS~bd2$Fish)
summary(m3)
# Changing Reference
bd2$Fish <- relevel(bd2$Fish, ref = "steindachneri")
str(bd2)
m4<-glm.nb(bd2$OTUS~bd2$Fish)
summary(m4)
#########################################################################




