#
# Plot order parameter moments for different simulated conditions
#

# Load libraries
library(ggplot2)
library(data.table)

# Set working directory
setwd("/Users/Igor/Personal/Github/microtubule_simulation/")

# Load files
decay_ordered.dt = fread("decay_or_ordered_state_phi_moments.txt")
