# Load in all colour schemes for consistency across the entire analysis

library(here)

# Datasets
ALL_DATASETS_COLOURS <- readRDS(here("data", "colours", "ALL_DATASETS_COLOURS.rds"))

# Samples
ALL_SAMPLES_COLOURS <- readRDS(here("data", "colours", "ALL_SAMPLES_COLOURS.rds"))
ALL_SAMPLES_COLOURS.PAPER <- readRDS(here("data", "colours", "ALL_SAMPLES_COLOURS.PAPER.rds"))

# Donors
ALL_DONORS_COLOURS <- readRDS(here("data", "colours", "ALL_DONORS_COLOURS.rds"))
ALL_DONORS_COLOURS.PAPER <- readRDS(here("data", "colours", "ALL_DONORS_COLOURS.PAPER.rds"))

# Foetal sex
ALL_FOETAL_SEX_COLOURS <- readRDS(here("data", "colours", "ALL_FOETAL_SEX_COLOURS.rds"))

# Gestational Age Groups
ALL_GEST_COLOURS <- readRDS(here("data", "colours", "ALL_GEST_COLOURS.rds"))

# Clusters
ALL_CLUSTER_COLOURS <- readRDS(here("data", "colours", "ALL_CLUSTER_COLOURS.rds"))

REFINED_CLUSTER_COLOURS <- readRDS(here("data", "colours", "REFINED_CLUSTER_COLOURS.rds"))
CLUSTER_NAME_COLOURS <- readRDS(here("data", "colours", "CLUSTER_NAME_COLOURS.rds"))

REFINED_CLUSTER_COLOURS.GERLI <- readRDS(here("data", "colours", "REFINED_CLUSTER_COLOURS.GERLI.rds"))
CLUSTER_NAME_COLOURS.GERLI <- readRDS(here("data", "colours", "CLUSTER_NAME_COLOURS.GERLI.rds"))

# Cell type labels
ALL_CELL_TYPE_COLOURS <- readRDS(here("data", "colours", "ALL_CELL_TYPE_COLOURS.rds"))