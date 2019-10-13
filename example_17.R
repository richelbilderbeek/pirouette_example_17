# Code of example 17
#
# Works under Windows
#
#
#

# Set the RNG seed
rng_seed <- 314
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  arg <- suppressWarnings(as.numeric(args[1]))
  if (is.na(arg)) {
    stop(
      "Please supply a numerical value for the RNG seed. \n",
      "Actual value: ", args[1]
    )
  }
  rng_seed <- arg
  if (rng_seed < 1) {
    stop("Please supply an RNG seed with a positive non-zero value")
  }
}
if (length(args) > 1) {
  stop(
    "Please supply only 1 argument for the RNG seed. \n",
    "Number of arguments given: ", length(args) - 1
  )
}

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
ape::write.tree(phylogeny, file = "tree_true.fasta")

alignment_params <- create_alignment_params(
  site_model = "unlinked_node_sub",
  fasta_filename = "alignment_gen.fasta",
  root_sequence = create_blocked_dna(length = 1000),
  mutation_rate = 0.01,
  rng_seed = rng_seed
)

experiment <- create_gen_experiment(
  beast2_options = beastier::create_beast2_options(
    input_filename = "input_gen.xml",
    output_log_filename = "output_gen.log",
    output_trees_filenames = "output_gen.trees",
    output_state_filename = "output_state.xml.state"
  )
)
experiments <- list(experiment)

# Set the RNG seed
for (i in seq_along(experiments)) {
  experiments[[i]]$beast2_options$rng_seed <- rng_seed
}

# Testing
if (1 == 2) {
  for (i in seq_along(experiments)) {
    experiments[[i]]$inference_model$mcmc <- create_mcmc(chain_length = 20000, store_every = 1000)
  }
}

pir_params <- create_pir_params(
  alignment_params = alignment_params,
  experiments = experiments,
  twinning_params = create_twinning_params(
    twin_model = "copy_true",
    rng_seed_twin_tree = rng_seed,
    rng_seed_twin_alignment = rng_seed,
    twin_tree_filename = "tree_twin.newick",
    twin_alignment_filename = "alignment_twin.fasta"
  )
)

# Make Peregrine friendly
pir_params <- peregrine::to_pff_pir_params(pir_params)
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