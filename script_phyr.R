# Set working directory
setwd("")

# Load libraries
library(ape)
library(glmmTMB)
library(ggplot2)
library(geiger)
library(car)
library(phytools)
library(nlme)
library(phylolm)  # for phyloglm
library(emmeans)
library(phyr)
library(DHARMa)
library(dplyr)
library(plyr)
library(sjPlot)
library(GGally)
library(ggeffects)
library(Rmisc)
library(performance)

# Load data from VinczeEtal2021Nature
# I combined all the variables from both original files into a single dataset to allow the development of more complex models for this example.

data <- read.csv2("data.csv", header = TRUE, stringsAsFactors = TRUE, sep = ";", dec = ".")
str(data)

# Convert variables to factors
data$Animal <- as.factor(data$Animal)
data$Vertebrate <- as.factor(data$Vertebrate)
data$Invertebrate <- as.factor(data$Invertebrate)
data$Mammal <- as.factor(data$Mammal)

# Create frequency table
freq_table <- table(data$order)

# Identify orders with more than one observation
valid_orders <- names(freq_table[freq_table > 1])

# Filter the dataset to keep only those orders
filtered_data <- subset(data, order %in% valid_orders)
filtered_data <- droplevels(filtered_data)

# Load phylogeny
phy <- read.nexus("consensus_phylogeny.tre")  # consensus VertLife tree
phy <- bind.tip(phy, "Cervus_canadensis", where = which(phy$tip.label == "Cervus_elaphus"),
                edge.length = 0.5, position = 0.5)
phy <- bind.tip(phy, "Gazella_marica", where = which(phy$tip.label == "Gazella_subgutturosa"),
                edge.length = 0.5, position = 0.5)

# For the CMR model, it is necessary to use a binomial distribution with a logit link function.
# In the phyr package, the binomial response must be specified as the number of successes (Neoplasia)
# and the number of failures (i.e., animals that died from causes other than Neoplasia).

# One important point to consider is that the phyr package only allows mixed models, so it is mandatory to include a random effect.
# "Fortunately", overdispersion is a common issue in binomial models.
# Thus, as a first approximation, we assume that overdispersion is present and include a random effect at the species level.
# This is a common method for accounting for overdispersion, known as OLRE (observation-level random effects).

# Before fitting Model 1 (M1), which includes phylogeny, we can fit a model that includes "Species" as a random effect, but does not account for phylogeny. This preliminary model allows us to test for overdispersion.
# However, it’s important to consider that any detected overdispersion might stem from not including phylogenetic structure.
# Unfortunately, there is currently no reliable way to assess overdispersion in a model that includes phylogeny.

summary(filtered_data)

# Remove NAs from Mammal variable
filtered_data2 <- filtered_data %>% filter(!is.na(Mammal))

# Fit Model 1 (includes phylogeny)

M1 <- pglmm(cbind(Neoplasia, knownDeaths - Neoplasia) ~ Adult_life_expectancy + order + body.mass + Vertebrate + (1 | Species),
            data = filtered_data2, family = "binomial", cov_ranef = list(sp = phy), add.obs.re = FALSE)
summary(M1)

# Check model residuals and zero-inflation
simulationOut_M1 <- simulateResiduals(fittedModel = M1, n = 1000)
plot(simulationOut_M1)
testZeroInflation(simulationOut_M1)

## Model diagnostics and interpretation
# First, we observe that the model does not fit the binomial distribution well, and the residuals are not randomly distributed.
# Nevertheless, this model provides the best possible option within a frequentist framework.
# An alternative would be to use Bayesian methods.
# Here, we model both the overdispersion and the phylogenetic non-independence between taxa.
# While not perfect, this model is consistent with the simpler model that excludes phylogeny (Model 0), and it provides a better fit overall.
# I prefer a suboptimal model whose limitations are clearly acknowledged over not modeling the data at all.

# Second, note that we included all taxa in the analysis, even those with a CMR of 0.
# This allowed us to test for zero-inflation, and no evidence of excess zeros was found.
# Using the complete dataset results in a more comprehensive analysis.

# Fit Model 0 (no phylogeny)
M0 <- glmmTMB(cbind(Neoplasia, knownDeaths - Neoplasia) ~ Adult_life_expectancy + order + body.mass + Vertebrate + (1 | Species),
              data = filtered_data2, family = "binomial")
summary(M0)

# Residual diagnostics for Model 0
simulationOut_M0 <- simulateResiduals(fittedModel = M0, n = 1000)
plot(simulationOut_M0)
testZeroInflation(simulationOut_M0)

# This model appears to perform adequately, but it does not account for phylogeny.
# Therefore, its inferences are unreliable due to the known phylogenetic signal in the data,
# which violates the assumption of independence between taxa.
                  



