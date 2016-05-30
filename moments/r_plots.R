#
# Plot order parameter moments for different simulated conditions
#

# Load libraries
library(ggplot2)
library(data.table)

# Set working directory
setwd("/Users/Igor/Personal/Github/microtubule_simulation/")

# Load files, annotate, melt & tidyfy.
decay_ordered.dt = fread("decay_or_ordered_state_phi_moments.txt")
names(decay_ordered.dt) = c("step", "moment_1", "moment_2", "moment_3")
decay_tidy.dt = melt(decay_ordered.dt, id.vars=c("step"), measure.vars=names(decay_ordered.dt)[2:4], value.name="moment")

# PLot
decay_ordered.plot = ggplot(decay_tidy.dt, aes(x=step, y=moment, group=variable, color=variable)) +
  geom_point(size=0.7) +
  scale_x_continuous() +
  scale_y_continuous("Distribution moment") +
  scale_color_discrete("Distribution\nmoment", labels=c(1,2,3))
decay_ordered.plot

