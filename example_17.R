#
# Difference from standard: 
# - Uses copy_tree for twin tree
# - Uses unlinked node-subst
#
library(pirouette)
suppressMessages(library(ggplot2))
library(beautier)

################################################################################
# Constants
################################################################################
is_testing <- is_on_travis()
example_no <- 17
rng_seed <- 314
folder_name <- paste0("example_", example_no, "_", rng_seed)

################################################################################
# Create phylogeny
################################################################################
phylogeny  <- ape::read.tree(
  text = "(((A:8, B:8):1, C:9):1, ((D:8, E:8):1, F:9):1);"
)

################################################################################
# Setup pirouette
################################################################################
pir_params <- create_std_pir_params(
  folder_name = folder_name
)
pir_params$alignment_params$sim_tral_fun <-
  get_sim_tral_with_uns_nsm_fun(
    branch_mutation_rate = 0.1,
    node_mutation_rate = 0.1
  )
pir_params$twinning_params$sim_twin_tree_fun <- create_copy_twtr_from_true_fun()

if (is_testing) {
  pir_params <- shorten_pir_params(pir_params)
}

################################################################################
# Run pirouette
################################################################################
pir_out <- pir_run(
  phylogeny,
  pir_params = pir_params
)

pir_save(
  phylogeny = phylogeny,
  pir_params = pir_params,
  pir_out = pir_out,
  folder_name = folder_name
)

