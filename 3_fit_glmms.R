#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# R Code supporting MS
#
# Vocal mimicry in corvids
# 
# Claudia A.F. Wascher*, Gemini Waterhouse & Bret A. Beheim
# corresponding author: claudia.wascher@gmail.com
#
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

options(digits=10)

library(dplyr)

# load data frame with info on species vocal behaviour, morphology, ecology, social factors etc

mimicry <- read.csv("data/corvid_mimicry_covariates.csv", header = TRUE, sep = ",", dec = ".")


################################################################################################
#Corvid Phylogenetic Tree
################################################################################################

# Package to get phylogeny from Open Tree Taxonomy
# OpenTreeOfLife, Benjamin Redelings, Luna Luisa Sanchez Reyes, Karen A. Cranston, Jim Allman, Mark T. Holder, & Emily Jane McTavish (2019) Open tree of life synthetic tree, Zenodo (12.3). https://doi.org/10.5281/zenodo.3937742
library(rotl)
library(ape)
# Create the data frame; match taxonomic names to the Open Tree Taxonomy
taxon_search <- tnrs_match_names(names = mimicry$scientific_name, context_name = "All life")
taxa <- mimicry$scientific_name
# Warning message: Cyanocorax colliei, Pica serica are not matched 
# Get the taxonomy ids for the taxa
resolved_names <- tnrs_match_names(taxa)
# Filter out taxa without a match
matched_names <- resolved_names[!is.na(resolved_names$ott_id), ]
# Extract the ott ids from the matched names
ott_ids <- matched_names$ott_id
# Fetch the phylogeny for the given ott ids
phylo <- tol_induced_subtree(ott_ids = ott_ids)
# Print the phylogeny
print(phylo)
# Save the tree to a file in Newick format
write.tree(phylo, file = "data/phylogeny_tree.newick")
# Define the path to Newick file
file_path <- 'data/phylogeny_tree.newick'
# Read the tree from the Newick file
treeCC <- read.tree(file = file_path)
# Apply the function to remove "OTT" from the tip labels
treeCC$tip.label <- gsub("_ott[0-9]+$", "", treeCC$tip.label, ignore.case = TRUE)

#Create a column animal, which is required for MCMCglmm
mimicry$animal <- as.factor(mimicry$scientific_name) # Pedigree option in MCMCglmm requires column called 'animal' 
#In column animal replace space with _ 
mimicry$animal <- gsub(" ", "_", mimicry$animal)


#rename species in dataset to match scientific names from phyologeny 

# Define the renaming list (modify with actual names)
rename_list <- c(
  "Coloeus_dauuricus" = "Corvus_dauuricus",
  "Coloeus_monedula" = "Corvus_monedula",
  "Cyanocorax_formosus" = "Gymnothorax_formosus"
)

# Apply renaming to the 'species' column
mimicry$animal <- ifelse(mimicry$animal %in% names(rename_list),
                          rename_list[mimicry$animal],
                          mimicry$animal)

  
################################################################################################
#Model investigating socio-ecological factors of specialist and generalist caching
################################################################################################
library(MCMCglmm)
# removing levels from data which are not in ginverse/phylogeny 
levels_in_data <- unique(mimicry$animal)
# Generate the ginverse matrix if not already done
ginverse <- inverseA(treeCC)$Ainv
# Extract the row names (animal labels)
levels_in_ginverse <- rownames(ginverse)
missing_levels <- setdiff(levels_in_data, levels_in_ginverse)
if (length(missing_levels) > 0) {
  cat("The following levels are missing in the ginverse matrix:\n", missing_levels, "\n")
} else {
  cat("All levels in the dataset are present in the ginverse matrix.\n")
}
mimicry <- mimicry[mimicry$animal %in% levels_in_ginverse, ]


# remove 'not reported' cases and change factors from character to integer
mimicry<-mimicry[!(mimicry$bodymass=="not reported"),]
mimicry$bodymass <- as.integer(mimicry$bodymass)
mimicry<-mimicry[!(mimicry$repertoire_combined=="not reported"),]
mimicry$repertoire_combined <- as.integer(mimicry$repertoire_combined)
mimicry<-mimicry[!(mimicry$breeding_system=="not reported"),]
mimicry<-mimicry[!(mimicry$habitat_breath=="not reported"),]
mimicry<-mimicry[!(mimicry$trophic_niche=="not reported"),]



# Count fixed effects: 1 for intercept + continuous predictors + (levels - 1) for each factor
num_fixed_effects <- 1 +  # Intercept
  9 

### Prior
Sigma <- diag(10) * 0.5  # Default variance 0.5
Sigma[1,1] <- 1  # Stronger prior on intercept

prior <- list(
  R = list(V = 1, nu = 0.002),
  G = list(G1 = list(V = 1, nu = 10)),
  B = list(mu = rep(0, 10), V = Sigma)
)

# MCMCglmm investigating whether different socio-ecological affect occurrence of vocal mimicry using real data; 
glmm_data <- MCMCglmm(mimicry~bodymass+repertoire_combined+habitat_breath+breeding_system+trophic_niche, random= ~genus, data=mimicry,
                 nitt=3000000, thin=2000, burnin=2000000, pr=TRUE, prior=prior, pedigree=treeCC, family = "categorical")
summary(glmm_data)

# MCMCglmm investigating whether different socio-ecological affect occurrence of vocal mimicry using imputed data 1
glmm_imputed1 <- MCMCglmm(mimicry_imputed1~bodymass+repertoire_combined+habitat_breath+breeding_system+trophic_niche, random= ~genus, data=mimicry,
                      nitt=3000000, thin=2000, burnin=2000000, pr=TRUE, prior=prior, pedigree=treeCC, family = "categorical")
# MCMCglmm investigating whether different socio-ecological affect occurrence of vocal mimicry using imputed data 2
glmm_imputed2 <- MCMCglmm(mimicry_imputed2~bodymass+repertoire_combined+habitat_breath+breeding_system+trophic_niche, random= ~genus, data=mimicry,
                          nitt=3000000, thin=2000, burnin=2000000, pr=TRUE, prior=prior, pedigree=treeCC, family = "categorical")
# MCMCglmm investigating whether different socio-ecological affect occurrence of vocal mimicry using imputed data 3
glmm_imputed3 <- MCMCglmm(mimicry_imputed3~bodymass+repertoire_combined+habitat_breath+breeding_system+trophic_niche, random= ~genus, data=mimicry,
                          nitt=3000000, thin=2000, burnin=2000000, pr=TRUE, prior=prior, pedigree=treeCC, family = "categorical")
# MCMCglmm investigating whether different socio-ecological affect occurrence of vocal mimicry using imputed data 4
glmm_imputed4 <- MCMCglmm(mimicry_imputed4~bodymass+repertoire_combined+habitat_breath+breeding_system+trophic_niche, random= ~genus, data=mimicry,
                          nitt=3000000, thin=2000, burnin=2000000, pr=TRUE, prior=prior, pedigree=treeCC, family = "categorical")
# MCMCglmm investigating whether different socio-ecological affect occurrence of vocal mimicry using imputed data 5
glmm_imputed5 <- MCMCglmm(mimicry_imputed5~bodymass+repertoire_combined+habitat_breath+breeding_system+trophic_niche, random= ~genus, data=mimicry,
                          nitt=3000000, thin=2000, burnin=2000000, pr=TRUE, prior=prior, pedigree=treeCC, family = "categorical")
# MCMCglmm investigating whether different socio-ecological affect occurrence of vocal mimicry using imputed data 6
glmm_imputed6 <- MCMCglmm(mimicry_imputed6~bodymass+repertoire_combined+habitat_breath+breeding_system+trophic_niche, random= ~genus, data=mimicry,
                          nitt=3000000, thin=2000, burnin=2000000, pr=TRUE, prior=prior, pedigree=treeCC, family = "categorical")
# MCMCglmm investigating whether different socio-ecological affect occurrence of vocal mimicry using imputed data 7
glmm_imputed7 <- MCMCglmm(mimicry_imputed7~bodymass+repertoire_combined+habitat_breath+breeding_system+trophic_niche, random= ~genus, data=mimicry,
                          nitt=3000000, thin=2000, burnin=2000000, pr=TRUE, prior=prior, pedigree=treeCC, family = "categorical")
# MCMCglmm investigating whether different socio-ecological affect occurrence of vocal mimicry using imputed data 8
glmm_imputed8 <- MCMCglmm(mimicry_imputed8~bodymass+repertoire_combined+habitat_breath+breeding_system+trophic_niche, random= ~genus, data=mimicry,
                          nitt=3000000, thin=2000, burnin=2000000, pr=TRUE, prior=prior, pedigree=treeCC, family = "categorical")
# MCMCglmm investigating whether different socio-ecological affect occurrence of vocal mimicry using imputed data 9
glmm_imputed9 <- MCMCglmm(mimicry_imputed9~bodymass+repertoire_combined+habitat_breath+breeding_system+trophic_niche, random= ~genus, data=mimicry,
                          nitt=3000000, thin=2000, burnin=2000000, pr=TRUE, prior=prior, pedigree=treeCC, family = "categorical")
# MCMCglmm investigating whether different socio-ecological affect occurrence of vocal mimicry using imputed data 10
glmm_imputed10 <- MCMCglmm(mimicry_imputed10~bodymass+repertoire_combined+habitat_breath+breeding_system+trophic_niche, random= ~genus, data=mimicry,
                          nitt=3000000, thin=2000, burnin=2000000, pr=TRUE, prior=prior, pedigree=treeCC, family = "categorical")

library(ggplot2)
# Mimicry and body mass
fig1<-ggplot(mimicry, aes(x = factor(mimicry), y = bodymass)) + geom_boxplot()+
  geom_jitter() +
  geom_smooth(
    method = "glm", 
    method.args = list(family = poisson)) + scale_y_continuous(labels = scales::comma)
fig1 + theme_classic() + labs(x="mimicry", y = "bodymass")


################################################################################################
#Model type (predator versus non-predator) and seasonal effects
################################################################################################
library(dplyr)
# load data frame with info on species vocal behaviour, morphology, ecology, social factors etc
mimicry <- read.csv("data/corvid_mimicry_covariates.csv", header = TRUE, sep = ",", dec = ".")
#only mimics
mimicry_filtered <- mimicry %>% 
  filter(mimicry != 0)


# Wilcoxon signed-rank test to test whether predators are mimicked more compared to non-predators
predator_outcome <- wilcox.test(mimicry_filtered$nr_predator, 
            mimicry_filtered$nr_nonpredator, 
            paired = TRUE, 
            exact = TRUE)  # Use exact=FALSE for larger datasets
# Ensure nrow() correctly calculates sample size
n <- as.numeric(nrow(mimicry_filtered))

# Print results with sample size
cat("Wilcoxon Signed-Rank Test Results:\n")
cat("V =", predator_outcome$statistic, "\n")  # V-statistic
cat("p-value =", predator_outcome$p.value, "\n")  # p-value
cat("Sample Size (n) =", n, "\n")  # Sample size


# Remove NAs for data of xeno-canto recordings to check for breeding versus non-breeding
mimicry_filtered_season <- mimicry %>% 
  filter(!is.na(nr_breedingseason))

# Wilcoxon signed-rank test to test whether mimicked calls are recorded more during the breeding season (xeno-canto)
breedingseason_outcome <- wilcox.test(mimicry_filtered_season$nr_breedingseason, 
                                mimicry_filtered_season$nr_outsidebreedingseason, 
                                paired = TRUE, 
                                exact = TRUE)  # Use exact=FALSE for larger datasets
# Ensure nrow() correctly calculates sample size
n <- as.numeric(nrow(mimicry_filtered_season))

# Print results with sample size
cat("Wilcoxon Signed-Rank Test Results:\n")
cat("V =", breedingseason_outcome$statistic, "\n")  # V-statistic
cat("p-value =", breedingseason_outcome$p.value, "\n")  # p-value
cat("Sample Size (n) =", n, "\n")  # Sample size

# Boxplot for visualization
boxplot(mimicry_filtered_season$nr_breedingseason, mimicry_filtered_season$nr_outsidebreedingseason,
        names = c("breeding season", "non-breeding season"),
        col = c("red", "blue"))





