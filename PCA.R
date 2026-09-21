# ======================================================================
# PCA Lab Meeting Script
# ======================================================================
#
# Main dataset:
#   Arabidopsis thaliana accessions × 13 climate variables
#   948 accessions
#
# Main question:
#   How does PCA turn many correlated variables into a small number of
#   synthetic axes, and how can variable choice change the ecological story?
#
# Sequence:
#   1. Prepare the variables
#   2. Fit a correlation-based PCA
#   3. Inspect scores and loadings
#   4. Combine them in a biplot
#   5. Inspect eigenvalues / scree plot and scaling
#   6. Change the variable set and watch the interpretation change
# ======================================================================


# ======================================================================
# Section 0. Packages
# ======================================================================

# Install these packages the first time you run the script:
# install.packages(c("FactoMineR", "factoextra", "ggplot2", "ggExtra", "ggrepel"))

library(FactoMineR)
library(factoextra)
library(ggplot2)


# ======================================================================
# Section 1. Data: Arabidopsis accessions × climate
# ======================================================================

clim <- read.csv("data/arabidopsis_climate.csv", row.names = 1, check.names = FALSE)

# One row = one Arabidopsis accession.
# LAT / LONG = geographic origin of the accession.
# The remaining 13 columns = climate variables at the accession's origin.

dim(clim)
colnames(clim)
head(clim)
summary(clim)


# ----------------------------------------------------------------------
# 1.1 Select the 13 climate variables used in the PCA
# ----------------------------------------------------------------------

clim_vars <- c("seasonal", "tempWarmest", "tempColdest", "preciWettest", "preciDriest", "preciCV", "PAR_SPRING", "growingL", "conseqCold", "conseqFrFree", "RelHumidSp", "dayLSp", "aridity")

clim_pca <- clim[, clim_vars]

colnames(clim_pca)


# ======================================================================
# Section 2. Preprocessing: transformation and standardisation
# ======================================================================

# Transformation and standardisation solve different problems.
#
# A transformation changes the shape of a variable and can reduce the
# influence of strong right skew or very large positive values.
#
# Standardisation changes the scale: each variable is expressed relative
# to its own mean and standard deviation.
#
# PCA does not require every variable to be normally distributed.
# The purpose of the log transformation here is not to "make PCA normal";
# it is to make strongly skewed positive climate variables less dominated
# by their largest values before examining their linear correlation structure.


# ----------------------------------------------------------------------
# 2.1 Inspect one example before transformation
# ----------------------------------------------------------------------

hist(clim_pca$preciCV, main = "Precipitation CV before transformation", xlab = "preciCV")
qqnorm(clim_pca$preciCV)
qqline(clim_pca$preciCV)


# ----------------------------------------------------------------------
# 2.2 Log-transform selected non-negative, right-skewed variables
# ----------------------------------------------------------------------

log_vars <- c("preciWettest", "preciDriest", "preciCV", "aridity")

clim_pca[log_vars] <- lapply(clim_pca[log_vars], function(x) log10(x + 1))

hist(clim_pca$preciCV, main = "Precipitation CV after log10(x + 1)", xlab = "log10(preciCV + 1)")
qqnorm(clim_pca$preciCV)
qqline(clim_pca$preciCV)


# ======================================================================
# Section 3. Fit the main PCA
# ======================================================================

# scale.unit = TRUE standardises every variable internally before PCA.
# Therefore this is a correlation-based PCA: each variable starts with
# variance = 1 and variables measured in different units become comparable.

pca_clim <- PCA(clim_pca, scale.unit = TRUE, graph = FALSE)
pca_clim


# ======================================================================
# Section 4. Scores, loadings, biplot
# ======================================================================

# scores   = coordinates of each observation in the new PC space
# loadings = relationships between the original variables and the PCs
#
# Once scores and loadings are understood separately, a biplot simply
# places both kinds of information in the same figure.


# ----------------------------------------------------------------------
# 4.1 Scores: where are the accessions in PC space?
# ----------------------------------------------------------------------

scores_pc12 <- pca_clim$ind$coord[, 1:2]

head(scores_pc12)

# Each row is one accession.
# An accession with a large positive PC1 score lies toward the positive
# end of the multivariate gradient represented by PC1.
#
# PCA axis signs are arbitrary. Multiplying every PC1 score and every PC1
# loading by -1 gives exactly the same PCA geometry and interpretation.

fviz_pca_ind(pca_clim, geom = "point")


# ----------------------------------------------------------------------
# 4.2 Loadings: which climate variables define each PC?
# ----------------------------------------------------------------------

loadings_pc12 <- pca_clim$var$cor[, 1:2]

round(loadings_pc12, 3)

# FactoMineR stores these as correlations between each original variable
# and each PC axis. They are therefore correlation loadings and range
# from -1 to +1.
#
# Large absolute value = strong alignment with that PC.
# Positive versus negative sign = opposite ends of the same axis.
# The absolute size matters more than the arbitrary orientation of the axis.


# ----------------------------------------------------------------------
# 4.3 Contributions: which variables contribute most to each PC?
# ----------------------------------------------------------------------

contrib_pc12 <- pca_clim$var$contrib[, 1:2]

round(contrib_pc12, 1)

fviz_contrib(pca_clim, choice = "var", axes = 1, top = length(clim_vars))

# Loadings and contributions are related but answer different questions.
#
# Loading:
#   How strongly is a variable associated with this PC?
#
# Contribution:
#   How much does this variable help define this PC relative to the other
#   variables included in this PCA?


# ----------------------------------------------------------------------
# 4.4 Biplot: scores + loadings in one view
# ----------------------------------------------------------------------

# First show the standard factoextra biplot.
fviz_pca_biplot(pca_clim, geom.ind = "point", repel = TRUE)

# The biplot combines:
#   points  = accession scores
#   arrows  = original climate variables
#
# Arrow direction:
#   shows the direction in PC space in which that variable increases.
#
# Arrow length:
#   shows how well that variable is represented in the displayed PC1-PC2
#   plane; a short arrow can mean that much of its information lies in
#   later PCs rather than that the variable is biologically unimportant.
#
# Arrow angles:
#   small angle     -> positive association in this projection
#   opposite arrows -> negative association in this projection
#   near 90 degrees -> weak association in this projection
#
# Important:
# these angles describe the PC1-PC2 projection, not necessarily the exact
# full-data correlation between every pair of variables.

# A contour-density version can make the dense centre of the score cloud easier

# Extract accession scores for PC1 and PC2.
scores_density <- data.frame(PC1 = pca_clim$ind$coord[, 1], PC2 = pca_clim$ind$coord[, 2])

# Store the PC1 and PC2 loadings for each climate variable.
loadings_density <- data.frame(variable = rownames(loadings_pc12), PC1 = loadings_pc12[, 1], PC2 = loadings_pc12[, 2])

# Rescale loading arrows so they fit naturally inside the score plot.
arrow_scale <- 0.70 * min(max(abs(scores_density$PC1)) / max(abs(loadings_density$PC1)), max(abs(scores_density$PC2)) / max(abs(loadings_density$PC2)))

# Convert loadings to plotting coordinates for the arrows and labels.
loadings_density$PC1_plot <- loadings_density$PC1 * arrow_scale
loadings_density$PC2_plot <- loadings_density$PC2 * arrow_scale

# Extract the percentage of variance explained by PC1 and PC2.
pc1_pct <- round(pca_clim$eig[1, 2], 1)
pc2_pct <- round(pca_clim$eig[2, 2], 1)

p_density <- ggplot(scores_density, aes(x = PC1, y = PC2)) + 
  geom_point(alpha = 0.18, size = 0.55, colour = "grey20") + 
  stat_density_2d(colour = "dodgerblue3", linewidth = 0.5) + 
  geom_segment(data = loadings_density, aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot), 
               inherit.aes = FALSE, arrow = grid::arrow(length = grid::unit(0.15, "cm")), colour = "black") + 
  ggrepel::geom_text_repel(data = loadings_density, aes(x = PC1_plot, y = PC2_plot, label = variable), 
                           inherit.aes = FALSE, size = 3, seed = 22092026, max.overlaps = Inf, colour = "black") + 
  labs(x = paste0("PC1 (", pc1_pct, "%)"), y = paste0("PC2 (", pc2_pct, "%)"), title = "") + 
  theme_classic()

p_density

# Add the one-dimensional score distributions to the top and right.
# Histogram = frequency view; density = smoothed distribution view.
ggExtra::ggMarginal(p_density, type = "histogram", bins = 50)
ggExtra::ggMarginal(p_density, type = "density")

# The quality of representation of
# observations in the displayed PC1-PC2 plane. High cos2 means that PC1 and PC2
# capture more of that accession's multivariate position.
fviz_pca_ind(pca_clim, geom = "point", col.ind = "cos2")


# ----------------------------------------------------------------------
# 4.5 Eigenvalues, explained variance and scree plot
# ----------------------------------------------------------------------

pca_clim$eig

fviz_eig(pca_clim, addlabels = TRUE)

# Eigenvalue:
#   amount of standardised variance captured by a PC.
#
# Percent variance:
#   eigenvalue divided by the total variance across all standardised variables.


# ======================================================================
# Section 5. Can adding only one or two variables change the story?
# ======================================================================

# PC1, PC2 and later PCs are not fixed biological components.
# They are recalculated from the complete set of variables supplied to PCA.
#
# This section deliberately changes the input variable set to demonstrate
# that PCA is descriptive: the axes summarise the data matrix we choose.


# ----------------------------------------------------------------------
# 5.1 Start with a five-variable climate PCA
# ----------------------------------------------------------------------

core_vars <- c("tempWarmest", "tempColdest", "preciWettest", "preciDriest", "aridity")

pca_core <- PCA(clim_pca[, core_vars], scale.unit = TRUE, graph = FALSE)

pca_core$eig[1:3, , drop = FALSE]
round(pca_core$var$cor[, 1:2], 3)
round(pca_core$var$contrib[, 1:2], 1)

fviz_pca_biplot(pca_core, geom.ind = "point", repel = TRUE)

# In this five-variable PCA, the leading interpretation is approximately:
#
# PC1 = water availability
# PC2 = temperature
#
# With the current data, PC1 explains about 44.6% and PC2 about 27.5%.
# The exact sign of each axis may be flipped without changing the result.


# ----------------------------------------------------------------------
# 5.2 Add two real variables: temperature seasonality and cold exposure
# ----------------------------------------------------------------------

# seasonal:
#   temperature seasonality
#
# conseqCold:
#   number of cold days
#
# Both represent a cold-seasonality
# dimension that was only indirectly represented in the five-variable PCA.

core_plus_cold <- c(core_vars, "seasonal", "conseqCold")

pca_core_plus_cold <- PCA(clim_pca[, core_plus_cold], scale.unit = TRUE, graph = FALSE)

pca_core_plus_cold$eig[1:3, , drop = FALSE]
round(pca_core_plus_cold$var$cor[, 1:2], 3)
round(pca_core_plus_cold$var$contrib[, 1:2], 1)

fviz_pca_biplot(pca_core_plus_cold, geom.ind = "point", repel = TRUE)


# ----------------------------------------------------------------------
# 5.3 The ecological story has changed, not just the axis labels
# ----------------------------------------------------------------------

# Compare the accession scores from the old and new PCAs.
# Because PCA signs are arbitrary, inspect both signed and absolute
# correlations before deciding whether an old axis corresponds to a new one.

score_axis_cor <- cor(pca_core$ind$coord[, 1:3], pca_core_plus_cold$ind$coord[, 1:3])

round(score_axis_cor, 2)
round(abs(score_axis_cor), 2)

# Before adding seasonal + conseqCold:
#
# PC1 ≈ 44.6%  -> water / aridity
# PC2 ≈ 27.5%  -> temperature
#
# After adding seasonal + conseqCold:
#
# PC1 ≈ 45.2%  -> temperature seasonality + cold exposure
# PC2 ≈ 28.8%  -> warm-dry / aridity gradient
#
# Approximate leading contributions after adding the two variables:
#
# New PC1:
#   seasonal      ≈ 27%
#   conseqCold    ≈ 23%
#   tempColdest   ≈ 21%
#
# New PC2:
#   tempWarmest   ≈ 30%
#   aridity       ≈ 26%
#   preciDriest   ≈ 15%
#
# The old axes are mixed across the new axes rather than showing a simple
# PC2 <-> PC3 relabelling. This is the key point:
#
# adding variables can change the ecological interpretation of the leading
# dimensions because PCA optimises variance for the current input matrix.


# ----------------------------------------------------------------------
# 5.4 Why did the story change?
# ----------------------------------------------------------------------

cold_cor <- cor(clim_pca[, c("tempColdest", "seasonal", "conseqCold")])

round(cold_cor, 2)

# In this dataset, approximately:
#
# seasonal vs tempColdest    r ≈ -0.88
# conseqCold vs tempColdest  r ≈ -0.87
# seasonal vs conseqCold     r ≈  0.84
#
# We therefore added two variables that strongly reinforce the same
# cold-seasonality dimension.
#
# PCA treats every input column as a variable. If one ecological process is
# represented by several strongly correlated columns while another process
# is represented by only one or two, the first process can receive more
# weight in the covariance/correlation structure simply because it has been
# represented repeatedly.
#
# This is a useful distinction:
#
# More variables != more independent ecological information.


# ----------------------------------------------------------------------
# 5.5 Variable selection should follow the ecological question
# ----------------------------------------------------------------------

# Example water variables:
#   precipitation amount
#   precipitation CV
#   soil moisture
#   
#   aridity
#   dry-spell duration
#
# Example temperature variables:
#   mean temperature
#   Tmax
#   Tmin
#   temperature seasonality
#   cold days
#   growing degree days
#   soil temp
#
# All may be biologically meaningful, but including many redundant variables
# from one process can make that process dominate the PCA.
#
# Before adding variables, ask:
#
# 1. What ecological process does this variable represent?
# 2. Does it add genuinely new information, or mostly duplicate another variable?
# 3. Am I giving one process more weight because I measured it in many ways?
# 4. Does my biological interpretation survive sensible changes in variable choice?


# ======================================================================
# Section 6. Correlation PCA versus covariance PCA
# ======================================================================

# The main analysis above uses correlation PCA because the climate variables
# have different units and numerical scales.
#
# To see why this choice matters, fit the same transformed data without
# standardisation. This is effectively covariance-based PCA.

pca_unscaled <- PCA(clim_pca, scale.unit = FALSE, graph = FALSE)

fviz_pca_biplot(pca_unscaled, geom.ind = "point", repel = TRUE)
fviz_pca_biplot(pca_clim, geom.ind = "point", repel = TRUE)

# Without scaling, variables with large raw variances can dominate.
# With scaling, every variable starts with variance = 1.
#
# Neither choice is universally correct: they answer different questions.
# For variables in very different units, correlation PCA is usually easier


# ======================================================================
# Section 7. Final take-home messages
# ======================================================================

# 1. PCA is descriptive dimension reduction for correlated multivariate data.
#
# 2. Scores are coordinates of observations on the new PC axes.
#
# 3. Loadings describe how strongly original variables align with each PC.
#
# 4. A biplot combines scores and loadings in the same geometric space.
#
# 5. Eigenvalues tell us how much variance each PC captures; scree and Kaiser
#    criteria are useful summaries, not substitutes for biological judgement.
#
# 6. Correlation PCA standardises variables; covariance PCA does not.
#    The choice changes the question.
#
# 7. PCA axis signs are arbitrary.
#
# 8. PC1 and PC2 are not fixed biological entities. They depend on the variables
#    included, their transformations, scaling, and correlation structure.
#
# 9. Adding correlated variables can overweight one ecological process and
#    genuinely change the apparent ecological story.
#
# 10. PCA is primarily exploratory. It does not by itself provide causal
#     inference or hypothesis-test P values.
