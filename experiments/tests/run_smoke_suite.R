source(file.path("experiments", "script_utils.R"))

project.root <- sbf_project_root()
scripts <- c(
  "check_benchmark_sbf_vs_origin.R",
  "check_additive_binning_experiment.R",
  "check_multiplicative_binning_experiment.R",
  "check_prediction_methods_method_check.R",
  "check_multiplicative_identification_method_check.R"
)

for (script in scripts) {
  path <- file.path(project.root, "experiments", "tests", script)
  cat(sprintf("Running %s\n", path))
  source(path)
}

cat("Experiment smoke suite passed.\n")
