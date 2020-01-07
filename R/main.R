library(stringr)

# Get Functions
source(file = "R/functions.R")

### CHECKING VERSIONS ####
correctVersions <- TRUE
ifelse(packageVersion("dplyr") >= "0.8.3", correctVersions <- TRUE, correctVersions <- FALSE)
ifelse(packageVersion("psych") >= "1.8.1", correctVersions <- TRUE, correctVersions <- FALSE)
ifelse(packageVersion("kableExtra") >= "1.1.0", correctVersions <- TRUE, correctVersions <- FALSE)
ifelse(packageVersion("knitr") >= "1.23", correctVersions <- TRUE, correctVersions <- FALSE)
ifelse(packageVersion("stringr") >= "1.4.0", correctVersions <- TRUE, correctVersions <- FALSE)
ifelse(packageVersion("purrr") >= "0.3.2", correctVersions <- TRUE, correctVersions <- FALSE)
ifelse(packageVersion("tidyr") >= "1.0.0", correctVersions <- TRUE, correctVersions <- FALSE)
ifelse(packageVersion("reshape2") >= "1.4.3", correctVersions <- TRUE, correctVersions <- FALSE)
ifelse(packageVersion("ggplot2") >= "3.2.0", correctVersions <- TRUE, correctVersions <- FALSE)
ifelse(packageVersion("readr") >= "1.3.1", correctVersions <- TRUE, correctVersions <- FALSE)

# SET RESULTS PATH
results_path <- "experiment_results/"

# Saved combained file for speed of analysis
mutationanalysis <- read_csv(file.path(results_path, dir(path = results_path, pattern = "mutationanalysis.dat")))
mergesFile <- read_csv(file.path(results_path, dir(path = results_path, pattern = "numberOfMerges.dat")))


# Change additional greedy -> greedy to match the function
mutationanalysis <- mutationanalysis %>% mutate(reducedwith = ifelse(reducedwith == "additionalGreedy", "greedy", reducedwith))
mutationanalysis <- mutationanalysis %>% mutate(reducedwith = if_else(fullreduce == FALSE, "original", reducedwith))

# Change the case study name to match
# Changes IsoFlav_R2Repaired to IsoFlav_R2
mutationanalysis <- mutationanalysis %>%
  mutate(casestudy = replace(casestudy,
                             casestudy == "IsoFlav_R2Repaired",
                             "IsoFlav_R2"))

# Get all the reduction time
if (!("reductionTime" %in% colnames(mutationanalysis))) {
  reduction_time <- reduction_time_analysis(mutationanalysis, func = "median")
  mutationanalysis <- left_join(mutationanalysis, reduction_time)
  rm(reduction_time)
} else {
  print("reductionTime is already in the mutationanalysis data frame")
}

# Get all the Mutation Score
if (!("mutationScore" %in% colnames(mutationanalysis))) {
  mutanttiming <- read_mutanttiming_files(results_path)
  # Change additional greedy -> greedy to match the functions
  mutanttiming <- mutanttiming %>% mutate(reducedwith = ifelse(reducedwith == "additionalGreedy", "greedy", reducedwith))
  mutationScores <- mutanttiming %>%
    filter(type == "NORMAL") %>%
    group_by(casestudy, datagenerator, dbms, randomseed, reducedwith) %>%
    summarise(killed_mutants = sum(killed == TRUE),
              total_mutants = (sum(killed == TRUE) + sum(killed == FALSE))) %>%
    mutate(mutationScore = round((killed_mutants/total_mutants) * 100, 2))

  # remove mutant timing
  rm(mutanttiming)

  mutationanalysis <- left_join(mutationanalysis, mutationScores)
  rm(mutationScores)
} else {
  print("mutationScore is already in the mutationanalysis data frame")
}

if (correctVersions == FALSE) {
  print("One of the packages has a version lower than required. Anyways, the files are executed but note there might be issues.")
} else {
  print("The packages versions are equal or greater than the required")
}