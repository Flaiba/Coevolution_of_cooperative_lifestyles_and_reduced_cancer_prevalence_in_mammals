
setwd("")

library(glmmTMB)
library(ggplot2)
library(car)
library(phytools)
library(nlme)
library(phylolm) 
library(emmeans)
library(phyr)
library(DHARMa)
library(dplyr)
library(plyr)
library(GGally)

#data base db1 
#The dataset initially included 191 species, but D. byrnei was removed due to its extremely high CMR, which was considered an outlier. 
DB1=read.csv2("db1.csv", header = TRUE,stringsAsFactors = T, sep=",", dec=".")
str(DB1)

# Charge the  phylogenetic tree
phy <- read.nexus("consensus_phylogeny.tre") 
phy <- bind.tip(phy, "Cervus_canadensis", where = which(phy$tip.label=="Cervus_elaphus"),
                edge.length=0.5, position = 0.5)
phy <- bind.tip(phy, "Gazella_marica", where = which(phy$tip.label=="Gazella_subgutturosa"),
                edge.length=0.5, position = 0.5)

#Previous investigations with the species included in these datasets showed there is a phylogenetic signal for CMR and neoplasia among mammal species (Vince et al. 2022, Boddy et al. 2020)

# O. Vincze, F. Colchero, J.-F. Lema?tre, D. A. Conde, S. Pavard, M. Bieuville, A. O. Urrutia, B. Ujvari, A. M. Boddy, C. C. Maley, Cancer risk across mammals. Nature 601, 263-267 (2022).
# A. M. Boddy, L. M. Abegglen, A. P. Pessier, A. Aktipis, J. D. Schiffman, C. C. Maley, C. Witte, Lifetime cancer prevalence and life history traits in mammals. Evolution, medicine, and public health 2020, 187-195 (2020)

#Detection of philogenetyc signal for CMR was probed in the paper Vincze, et al. Cancer risk across mammals. Nature

#Transform variables to detect lineal effect apply log 
DB1$gestation_length_log=log(DB1$gestation_length_d)
DB1$total_litters_log=log(DB1$total_litters)
DB1$litter_size_n_log=log(DB1$litter_size_n)
DB1$body_mass_kg_log=log(DB1$body_mass_kg)
DB1$A_Lexp_log=log(DB1$life_expectancy_d)
DB1$Met_rate_log=log(DB1$metabolic_rate)

DB1$animal=as.factor(DB1$animal) #Defining animal as categorical variable

poly=subset(DB1,litters=="Polytocous")
mono=subset(DB1,litters=="Monotocous")

#Model S1_a phylGLMM for continuous variables
S1_a=pglmm(cbind(n_neoplasia, n_no_neoplasia) ~ log(body_mass_kg)+litter_size_n_log+log(gestation_length_d)+log(life_expectancy_d )+(1|species), data = DB1, family = "binomial",cov_ranef = list(sp = phy),add.obs.re = FALSE)
summary(S1_a)
simulationOut<- simulateResiduals(fittedModel = S1_a, n = 1000,refit = FALSE)
testZeroInflation(simulationOut)
Pre_S1_a_model=pglmm_predicted_values(S1_a, type="response")

#Model S1_b Metabolic rate
S1_b=pglmm(cbind(n_neoplasia, n_no_neoplasia) ~ log(metabolic_rate)+(1|species), data = DB1, family = "binomial",cov_ranef = list(sp = phy),add.obs.re = FALSE)
summary(S1_b)
simulationOut<- simulateResiduals(fittedModel = S1_b, n = 1000,refit = FALSE)
testZeroInflation(simulationOut)
PreMet_model=pglmm_predicted_values(S1_b, type="response")

#Modelto correlate Total litters vs Litters size
S2=pglmm(litter_size_n_log~ total_litters_log+(1|Species), data = Naturesout,cov_ranef = list(sp = phy),add.obs.re = FALSE)
summary(S2)
simulationOut<- simulateResiduals(fittedModel = S2, n = 1000,refit = FALSE)
testZeroInflation(simulationOut)


#Model S2_a (Litters as categorical factor)
S2_a=pglmm(cbind(n_neoplasia, n_no_neoplasia) ~ litters+(1|species), data = DB1, family = "binomial",cov_ranef = list(sp = phy),add.obs.re = FALSE)
summary(S2_a)
simulationOut<- simulateResiduals(fittedModel = S2_a, n = 1000,refit = FALSE)
testZeroInflation(simulationOut)
Pre_S2_a_model=pglmm_predicted_values(S2_a, type="response")

#Model S2 b (Litters as categorical factor and continous variables)
S2_b=pglmm(cbind(n_neoplasia, n_no_neoplasia) ~ litters* log(body_mass_kg)+log(life_expectancy_d)+(1|species), data = DB1, family = "binomial",cov_ranef = list(sp = phy),add.obs.re = FALSE)
summary(S2_b)
simulationOut<- simulateResiduals(fittedModel = S2_b, n = 1000,refit = FALSE)
testZeroInflation(simulationOut)
Pre_S2_a_model=pglmm_predicted_values(S2_b, type="response")


#Model S3 a Group living
S3_a=pglmm(cbind(n_neoplasia, n_no_neoplasia) ~group_living +(1|species), data = DB1, family = "binomial",cov_ranef = list(sp = phy),add.obs.re = FALSE)
summary(S3_a)
simulationOut<- simulateResiduals(fittedModel = S3_a, n = 1000,refit = FALSE)
testZeroInflation(simulationOut)
Pre_S3_a_model=pglmm_predicted_values(S3_a, type="response")

S3_b=pglmm(cbind(n_neoplasia, n_no_neoplasia) ~litter_size_n_log+A_Lexp_log+group_living *body_mass_kg_log+(1|species), data = DB1, family = "binomial",cov_ranef = list(sp = phy),add.obs.re = FALSE)
summary(S3_b)
simulationOut<- simulateResiduals(fittedModel = S3_b, n = 1000,refit = FALSE)
testZeroInflation(simulationOut)
Pre_S3_a_model=pglmm_predicted_values(S3_b, type="response")

# Model S4 a Breeding System
S4_a=pglmm(cbind(n_neoplasia, n_no_neoplasia) ~breeding_system  +(1|species), data = DB1, family = "binomial",cov_ranef = list(sp = phy),add.obs.re = FALSE)
summary(S4_a) 
simulationOut<- simulateResiduals(fittedModel = S4_a, n = 1000,refit = FALSE)
testZeroInflation(simulationOut)
Pre_S4_a_model=pglmm_predicted_values(S4_a, type="response")

S4_b=pglmm(cbind(n_neoplasia, n_no_neoplasia) ~litter_size_n_log+A_Lexp_log+breeding_system  *body_mass_kg_log+(1|species), data = DB1, family = "binomial",cov_ranef = list(sp = phy),add.obs.re = FALSE)
summary(S4_b) 
simulationOut<- simulateResiduals(fittedModel = S4_a, n = 1000,refit = FALSE)
testZeroInflation(simulationOut)
Pre_S4_b_model=pglmm_predicted_values(S4_b, type="response")

#Paternal Care 
S5_a=pglmm(cbind(n_neoplasia, n_no_neoplasia) ~paternal_care +(1|species), data = DB1, family = "binomial",cov_ranef = list(sp = phy),add.obs.re = FALSE)
summary(S5_a) 
simulationOut<- simulateResiduals(fittedModel = S5_a, n = 1000,refit = FALSE)
testZeroInflation(simulationOut)
Pre_S5_a_model=pglmm_predicted_values(S5_a, type="response")

#Animal diet (yes/no)
S7_Diet_animal=pglmm(cbind(n_neoplasia, n_no_neoplasia) ~ animal+ (1|species), data = DB1, family = "binomial",cov_ranef = list(sp = phy),add.obs.re = FALSE)
summary(S7_Diet_animal)
simulationOut<- simulateResiduals(fittedModel = S7_Diet_animal, n = 1000,refit = FALSE)
testZeroInflation(simulationOut)
Pre_S7_a_model=pglmm_predicted_values(S7_Diet_animal, type="response")

### Separate the dataset according of Group Living (Solitary sp NO, and  Gregarious sp.YES) 
GL_yes=subset(DB1,group_living=="yes")
GL_no=subset(DB1,group_living=="no")

S7_Diet_animal_Gregarious=pglmm(cbind(n_neoplasia, n_no_neoplasia) ~ animal+(1|species), data = GL_yes, family = "binomial",cov_ranef = list(sp = phy),add.obs.re = FALSE)
summary(S7_Diet_animal_Gregarious)
simulationOut<- simulateResiduals(fittedModel = S7_Diet_animal_Gregarious, n = 1000,refit = FALSE)
testZeroInflation(simulationOut)
Pre_S7_GLyes_model=pglmm_predicted_values(S7_Diet_animal_Gregarious, type="response")

S7_Diet_animal_Solitary=pglmm(cbind(n_neoplasia, n_no_neoplasia) ~ animal+(1|species), data = GL_no, family = "binomial",cov_ranef = list(sp = phy),add.obs.re = FALSE)
summary(S7_Diet_animal_Solitary)
simulationOut<- simulateResiduals(fittedModel = S7_Diet_animal_Solitary, n = 1000,refit = FALSE)
testZeroInflation(simulationOut)
Pre_S7_GLno_model=pglmm_predicted_values(S7_Diet_animal_Solitary, type="response")

### Separate the dataset according of Breeding System (Singular and Plural)

BS_Pl=subset(DB1,breeding_system=="PluralBreeder")
BS_Si=subset(DB1,breeding_system=="SingularBreeder")

S7_Diet_animal_BS_Pl=pglmm(cbind(n_neoplasia, n_no_neoplasia) ~ animal+(1|species), data = BS_Pl, family = "binomial",cov_ranef = list(sp = phy),add.obs.re = FALSE)
summary(S7_Diet_animal_BS_Pl)
simulationOut<- simulateResiduals(fittedModel = S7_Diet_animal_BS_Pl, n = 1000,refit = FALSE)
testZeroInflation(simulationOut)
Pre_S7_Diet_animal_BS_Pl=pglmm_predicted_values(S7_Diet_animal_BS_Pl, type="response")

S7_Diet_animal_BS_S=pglmm(cbind(n_neoplasia, n_no_neoplasia) ~ animal+(1|species), data = BS_Si, family = "binomial",cov_ranef = list(sp = phy),add.obs.re = FALSE)
summary(S7_Diet_animal_BS_S)
simulationOut<- simulateResiduals(fittedModel = S7_Diet_animal_BS_S, n = 1000,refit = FALSE)
testZeroInflation(simulationOut)
Pre_S7_Diet_animal_BS_S=pglmm_predicted_values(S7_Diet_animal_BS_S, type="response")

# Analysis of indexes for traits : (a) Litters Index: ratio between monotocous and polytocous species within each order, (b) Group Living Index: ratio between the species with and without Group Living within each order, and (c) Breeding System Index: ratio between plural and singular breeding species within each order

Index_orders=read.csv2("table_index_orders.csv", header = TRUE,stringsAsFactors = T, sep=";", dec=".")
str(Index_orders)

I_GL<- glmmTMB(cbind(n_cancer, n_no_cancer)~  index_gl, data=Index_orders,family = "binomial")
summary(I_GL)

I_BS<- glmmTMB(cbind(n_cancer, n_no_cancer)~  index_BS, data=Index_orders,family = "binomial")
summary(I_BS)

I_Litt<- glmmTMB(cbind(n_cancer, n_no_cancer)~  index_litters, data=Index_orders,family = "binomial")
summary(I_Litt)

#data base db2 
DB2=read.csv2("db2.csv", header = TRUE,stringsAsFactors = T, sep=",", dec=".")
str(DB2)

# Charge the  phylogenetic tree
tree2 <- read.tree("tree_for_database2.txt")

#Previous investigations with the species included in these datasets showed there is a phylogenetic signal for CMR and neoplasia among mammal species (Vince et al. 2022, Boddy et al. 2020)

DB2$M_rate_log=log(DB2$metabolic_rate)
DB2$gestation_length_d_log=log(DB2$gestation_length_d)
DB2$litter_size_n_log=log(DB2$litter_size_n)
DB2$adult_mass_log=log(DB2$adult_mass_kg)
DB2$max_lifespan_yr_log=log(DB2$max_lifespan_yr)
DB2$total_litters_log=log(DB2$total_litters)


#Gestation length
m1=pglmm(cbind(any_neoplasia , any_no_neoplasia ) ~gestation_length_d_log+(1|species), data = DB2, family = "binomial",cov_ranef = list(sp = tree2),add.obs.re = FALSE)
summary(m1)

#Adult mass
m2=pglmm(cbind(any_neoplasia , any_no_neoplasia ) ~adult_mass_log+(1|species), data = DB2, family = "binomial",cov_ranef = list(sp = tree2),add.obs.re = FALSE)
summary(m2)

#Life span
m3=pglmm(cbind(any_neoplasia , any_no_neoplasia ) ~max_lifespan_yr_log+(1|species), data = DB2, family = "binomial",cov_ranef = list(sp = tree2),add.obs.re = FALSE)
summary(m3) 

#Litter Size
m4=pglmm(cbind(any_neoplasia , any_no_neoplasia ) ~litter_size_n_log+(1|species), data = DB2, family = "binomial",cov_ranef = list(sp = tree2),add.obs.re = FALSE)
summary(m4) 

#Metabolic rate
m5=pglmm(cbind(any_neoplasia , any_no_neoplasia ) ~M_rate_log+(1|species), data = DB2, family = "binomial",cov_ranef = list(sp = tree2),add.obs.re = FALSE)
summary(m5) 

# Litters 
m6=pglmm(cbind(any_neoplasia , any_no_neoplasia ) ~litters+(1|species), data = DB2, family = "binomial",cov_ranef = list(sp = tree2),add.obs.re = FALSE)
summary(m6)

# Paternal care
m7=pglmm(cbind(any_neoplasia , any_no_neoplasia ) ~paternal_care+(1|species), data = DB2, family = "binomial",cov_ranef = list(sp = tree2),add.obs.re = FALSE)
summary(m7)

# to correct the p-values, we used the false discovery method 

p.adjust(pvalues_list,method="fdr")
#Due the Value_List object is a list of all P values corresponds to all models made for each data set
