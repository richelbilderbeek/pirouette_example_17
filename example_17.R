# Code of example 17
#
# Works under Windows
#
#
#

# Set the RNG seed
rng_seed <- 314

library(pirouette)
suppressMessages(library(ggplot2))
library(beautier)

root_folder <- getwd()
example_no <- 17
example_folder <- file.path(root_folder, paste0("example_", example_no, "_", rng_seed))
dir.create(example_folder, showWarnings = FALSE, recursive = TRUE)
setwd(example_folder)
set.seed(rng_seed)
testit::assert(is_beast2_installed())
phylogeny  <- ape::read.tree(
  text = "(((A:8, B:8):1, C:9):1, ((D:8, E:8):1, F:9):1);"
)
ape::write.tree(phylogeny, file = "tree_true.fas")

alignment_params <- create_alignment_params(
  sim_true_alignment_fun =
    get_sim_true_alignment_with_unlinked_node_sub_site_model_fun(
      branch_mutation_rate = 0.1,
      node_mutation_rate = 0.1
    ),
  fasta_filename = "alignment_gen.fas",
  root_sequence = create_blocked_dna(length = 1000),
  rng_seed = rng_seed
)

experiment <- create_gen_experiment()
experiments <- list(experiment)

# Set the RNG seed
for (i in seq_along(experiments)) {
  experiments[[i]]$beast2_options$rng_seed <- rng_seed
}

# Testing
if (1 == 1) {
  for (i in seq_along(experiments)) {
    experiments[[i]]$inference_model$mcmc <- create_mcmc(chain_length = 20000, store_every = 1000)
  }
}

pir_params <- create_pir_params(
  alignment_params = alignment_params,
  experiments = experiments,
  twinning_params = create_twinning_params(
    rng_seed_twin_tree = rng_seed,
    sim_twin_tree_fun = create_copy_twin_tree_from_true_fun(),
    rng_seed_twin_alignment = rng_seed,
    sim_twin_alignment_fun = get_sim_twin_alignment_with_std_site_model_fun(),
    twin_tree_filename = "tree_twin.newick",
    twin_alignment_filename = "alignment_twin.fas"
  )
)

rm_pir_param_files(pir_params)

errors <- pir_run(
  phylogeny,
  pir_params = pir_params
)

utils::write.csv(
  x = errors,
  file = file.path(example_folder, "errors.csv"),
  row.names = FALSE
)

pir_plot(errors) +
  ggsave(file.path(example_folder, "errors.png"))

pir_to_pics(
  phylogeny = phylogeny,
  pir_params = pir_params,
  folder = example_folder
)

pir_to_tables(
  pir_params = pir_params,
  folder = example_folder
)
