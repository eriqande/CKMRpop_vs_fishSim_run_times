# This is just like 001, except that for pair-finding in
# fishSim, quickin is used rather than

#### Get the command line arg to Rscript to set the cores ####
args = commandArgs(trailingOnly=TRUE)
useCores = as.integer(args[1])

#### Get the necessary packages ####

if(!("remotes" %in% rownames(installed.packages()))) {
  install.packages("remotes")
}
# get version of fishSim with four-generation relationship finding
remotes::install_github(
  "eriqande/fishSim",
  ref = "12c85524b0c8"
)
remotes::install_github(
  "eriqande/CKMRpop",
  ref = "8c18687c71e9"
)
# install the spip binary
if(!CKMRpop::spip_exists())
  CKMRpop::install_spip(Dir = system.file(package = "CKMRpop"))
# get tidyverse if it aint there
if(!("tidyverse" %in% rownames(installed.packages()))) {
  install.packages("tidyverse")
}



#### Load the packages ####
# Want to try out Shane Baylis' package (that Paul Conn just told me about)
# a little bit more.
library(fishSim)
library(CKMRpop)
library(tidyverse)



#### Parameters for simulations ####

# just adapting from the fishSim vignette
p_maxAge <- 20  # max age of founders to produce
p_ageMort <- matrix(c(0.47, 0.37, 0.27, rep(0.23, 97)), ncol = 1) # age specific mortality rate
p_maleCurve <- c(0,0,0,0.1,0.5,0.8,0.85,0.9,0.95,rep(1, 91)) # age specific male prob of reproducing
p_femaleCurve <- c(0,0,0.5,0.9,0.95,rep(1,95)) # age specific female prob of reproducing

#### Function to translate the same scenario into spip ####
translate_to_CKMRpop <- function(
  N = 1000, # desired pop size
  maxAge = p_maxAge,  # max age of founders to produce
  ageMort = p_ageMort, # age specific mortality rate
  maleCurve = p_maleCurve, # age specific male prob of reproducing
  femaleCurve = p_femaleCurve
) {
  S <- list()
  S$`max-age` <- maxAge
  S$`fem-surv-probs` <- 1 - ageMort[,1][1:maxAge]
  S$`male-surv-probs` <- 1 - ageMort[,1][1:maxAge]
  S$`fem-prob-repro` <- femaleCurve[1:maxAge]
  S$`male-prob-repro` <- maleCurve[1:maxAge]
  S$`fem-asrf` <- rep(1, length.out = 20)
  S$`male-asrp` <- rep(1, length.out = 20)
  S$`offsp-dsn` <- "pois"
  S$`mate-fidelity` <- -1
  S$`sex-ratio` <- 0.5
  S$`number-of-years` <- 100  # hard wire this to be the same as the fishSim

  # now, get the stable age distribution.  It is the same for males and females
  # in this scenario.  Here it is for a cohort size of N:
  SAD <- leslie_from_spip(S, N)$stable_age_distro_fem

  # so, the desired cohort size would be
  dcs <- round(N * (N/sum(SAD)))

  # so put that in there, and set the initial number of males and females
  init <- round(SAD / sum(SAD) * N)
  S$`initial-males` <- init
  S$`initial-females` <- init

  # tell spip to use the cohort size
  S$`cohort-size` <- paste("const", dcs, collapse = " ")

  # now, do the sampling.  We are running this simulation to 120 years,
  # so the equivalent of years 94 to 100 will be
  # years 114 to 120.  We will sample all age classes at that same rate...
  samp_frac <- 0.04
  samp_start_year <- 114
  samp_stop_year <- 120
  S$`discard-all` <- 0
  S$`gtyp-ppn-fem-pre` <- paste(
    samp_start_year, "-", samp_stop_year, " ",
    paste(rep(samp_frac, maxAge), collapse = " "),
    sep = "")
  S$`gtyp-ppn-male-pre` <- S$`gtyp-ppn-fem-pre`

  S
}

#### Take a moment to compute the stable age ratio from CKMRpop ####

# It turns out that population size regulation is done differently
# In spip than in fishSim.  In spip, an appropriate number of
# offspring each year is chosen, whereas in fishSim, a number of offspring
# per female is set.  This makes fishSim far more sensitive to
# small random variations in population size, and therefore it
# is a little less stable.

# So, I can't get each simulation to have exactly the same numbers
# of individuals or samples, but I will record those values
# and note them in the comparison.

# And also find the number of offspring per reproductive female
spip_pars <- translate_to_CKMRpop(N = 1000)
SAD <- leslie_from_spip(spip_pars, 1000)$stable_age_distro_fem
p_SAD <- SAD/sum(SAD)

CZ <- as.numeric(strsplit(spip_pars$`cohort-size`,split = " +")[[1]][2])

offspring_per_reproductive_female <- CZ / sum(spip_pars$`initial-females` * spip_pars$`fem-prob-repro`)

#### Function to implement fishSim scenario ####

# the mortality specification in fishSim seems to be required
# to carry out well
# past the max age of the founders. spip, on the other hand,
# just kills everyone off at a max age.

# altMate doesn't seem set up for fecundity differences in different age groups.
# That is OK.  I just want a speed comparison here.

run_fishSim <- function(
  N = 1000, # desired pop size
  maxAge = p_maxAge,  # max age of founders to produce
  ageMort = p_ageMort, # age specific mortality rate
  maleCurve = p_maleCurve, # age specific male prob of reproducing
  femaleCurve = p_femaleCurve, # age specific female prob of reproducing
  SAD = p_SAD,  # the stable age distribution that we want to start from
  batch_Size = offspring_per_reproductive_female
) {

  # make the founders
  indiv <- makeFounders(
    pop = N,
    stocks = c(1),
    maxAge = maxAge,
    survCurv = SAD
  )

  # set up the archive
  archive <- make_archive()

  # check the growth rate
  #check_growthrate(
  #  mateType = "ageSex",
  #  mortType = "ageStock",
  #  batchSize = 1.67436, # I fiddled this to get it as close as possible to 1
  #  femaleCurve = femaleCurve,
  #  ageStockMort = ageMort
  #)

  # now, run it forward for 100 years.
  for (k in 1:100) {

    ## mate animals using the age-specific, sex-specific curves we set up before
    indiv <- altMate(
      indiv = indiv,
      batchSize = batch_Size,
      type = "ageSex",
      maleCurve = maleCurve,
      femaleCurve = femaleCurve,
      singlePaternity = FALSE, # make it full on random mating
      year = k
    )

    ## kill animals on the basis of their age and stock, as set up before
    indiv <- mort(
      indiv = indiv,
      year = k,
      type = "ageStock",
      ageStockMort = ageMort
    )

    # archive every 10 years to improve speed
    if(k %% 10 == 0) {
      archive <- archive_dead(indiv = indiv, archive = archive)
      indiv <- remove_dead(indiv = indiv)
      cat( sprintf( '\rIteration %i complete', k))
    }

    # in the last 6 years of simulation, lethally sample 4% of the population
    if(k %in% c(94:100)) {
      indiv <- capture(
        indiv,
        n = ceiling(0.04 * N),
        fatal = TRUE,
        year = k)
    }

    # everyone has a birthday!
    indiv <- birthdays(indiv = indiv)
  }
  archive <- rbind(archive, indiv)
  indiv <- archive  ## merge 'indiv' and 'archive', since they were only separated for speed.

  indiv

}



#### Make a function to run both fishSim and CKMRpop and record times ####

# We do this in two different steps which are separate:
#  1. Simulate the population (this includes slurping it up in spip)
#  2. Find the kin pairs

#' @param n the size of the initial population.
run_and_time <- function(n, nCores = 1) {

  #### do the spip sim and get the times, etc. ####
  spars <- translate_to_CKMRpop(N=n)

  # time the populations sim:
  spip_sim_time <- system.time(spip_dir <- run_spip(pars = spars))

  # time the pair finding:
  spip_pair_time <- system.time(
    {
      slurped <- slurp_spip(spip_dir, num_generations = 2)
      crel <- compile_related_pairs(slurped$samples)
    }
  )
  # summarise the simulation results conditions
  spip_size <- nrow(slurped$pedigree)
  spip_sample_size <- nrow(slurped$samples)


  #### do the fishSim sims and get the times, etc ####
  # time the population simulation
  fishSim_sim_time <- system.time(indivs <- run_fishSim(N = n))

  # get the simulation results
  fishSim_size <- nrow(indivs)
  fishSim_samples <- sum(!is.na(indivs$SampY))

  # time the pair-finding
  fishSim_pair_time <- system.time(pairs <- quickin(indivs, max_gen = 2))

  # put it all in a tibble
  summary <- tibble(
    N = n,
    spip_size = spip_size,
    fish_sim_size = fishSim_size,

    spip_sim_time = list(spip_sim_time),
    fishSim_sim_time = list(fishSim_sim_time),

    spip_sample_size = spip_sample_size,
    fish_sim_sample_size = fishSim_samples,

    spip_pair_time = list(spip_pair_time),
    fishSim_pair_time = list(fishSim_pair_time)
  )

  summary

}

#### Now, set up a grid of parameters to run and do it ####

sim_conditions <- expand_grid(
  n = c(500, 1000, 2500, 5000, 10000, 25000, 50000, 100000),
  nCores = useCores
)

results <- sim_conditions %>%
  mutate(
    res = map2(
      .x = n,
      .y = nCores,
      .f = ~ run_and_time(n = .x, nCores = .y)
    )
  )

write_rds(results, file = paste0("sim_results_cores_", useCores, ".rds"))

