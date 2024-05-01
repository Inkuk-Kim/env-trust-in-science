
## April 7 2024

#This is the script for the study INKUK is leading on place of science in science and environmmental attitudes
# for saving models. # INKUK set path fo your computer
push_mods <-
  fs::path_expand(
    "/Users/kimin/Documents/2024-env-trust-science-place-in-society/files"
  )


# libraries for jb (when internet is not accessible)
# read libraries
# source("/Users/joseph/GIT/templates/functions/libs2.R")

# read functions
#source("/Users/joseph/GIT/templates/functions/funs.R")


library(margot)
library(tidyverse)
library(lmtp)
library(kableExtra)
# experimental functions (more functions)
source(
  "https://raw.githubusercontent.com/go-bayes/templates/main/functions/libs2.R"
)

# experimental functions (more functions)
# source(
#   "https://raw.githubusercontent.com/go-bayes/templates/main/functions/experimental_funs.R"
# )


# read data/ set to path in your computer
pull_path <-
  fs::path_expand(
    "nzavs_data_qs"
  )


# read data: note that you need use the arrow package in R
dat <- qs::qread(pull_path)



# check path:is this correct?  check so you know you are not overwriting other directors
push_mods


# for information on LMPT go to
#https://github.com/nt-williams/lmtp
#devtools::install_github("nt-williams/lmtp@devel")

# for modified treatment policies
library("lmtp")
# push_mods


# set number of folds for ML here. use a minimum of 5 and a max of 10
SL_folds = 10

#this will allow you to track progress
progressr::handlers(global = TRUE)

# set seed for reproducing results
set.seed(0112358)

# set cores for estimation
library(future)
plan(multisession)
n_cores = 10  #<- parallel::detectCores() - 2 # save two cores for other work while these models run

# min of 10
n_cores

# super learner libraries
# these are useful for high-dimensional data
sl_lib <- c("SL.glmnet",
            "SL.ranger", # forests
            "SL.xgboost") # grandient boost

# libraries
library(SuperLearner)
library(ranger)
library(xgboost)
library(glmnet)


# check options
listWrappers()



# data preparation --------------------------------------------------------

# ensure there are no 'haven' labels
dat <- as.data.frame(dat)
dat <- haven::zap_formats(dat)
dat <- haven::zap_label(dat)
dat <- haven::zap_widths(dat)



# save total n ------------------------------------------------------------

n_total <- skimr::n_unique( dat$id )
n_total
# margot::here_save(n_total, "n_total")

nzavs_exposure <- "trust_science_our_society_places_too_much_emphasis_reversed"


# filter the original dataset for these IDs three waves
# get ids
ids_2019 <- dat |>
  filter( year_measured == 1 & wave == 2019 & !is.na(!!sym(nzavs_exposure))) %>% 
           pull(id)
ids_2020 <- dat |>
  filter( year_measured == 1 & wave == 2020  & !is.na(!!sym(nzavs_exposure))) %>% 
  pull(id)

ids_2019_2020 <- intersect(ids_2019, ids_2020)

# ids_2019 <- dat %>%
#   filter(year_measured == 1, wave %in% c("2019", "2020") &
#            !is.na(!!sym(nzavs_exposure))) |> # criteria, no missing
#   pull(id)
# 



dat$trust_science_high_confidence_scientific_community
# criteria for inclusion:  enrolled in 2018, might have been lost to follow up at any point after
# outcomes measured one year + exposure year
dat_long <- dat |>
  dplyr::filter(id %in% ids_2019_2020 &
                  wave %in% c(2019, 2020, 2021)) |>
  arrange(id, wave) |>
  # ensure this is a factor
  dplyr::rename(sample_weights = w_gend_age_euro) |>
  mutate(
    covid19_timeline = as.factor(covid19_timeline),
    # church attendancy as binary
    # religion_church_binary = ifelse(religion_church > 0, 1, religion_church),  # not in 2019
    # make numeric for easier imputation
    male = as.numeric(male),
    #asier imputation
    education_level_coarsen = as.integer(education_level_coarsen),
    # someone gave neg number
    household_inc_log = log(household_inc + 1),
    hours_children_log = log(hours_children + 1),
    hours_work_log = log(hours_work + 1),
    hours_housework_log = log(hours_housework + 1),
    hours_exercise_log = log(hours_exercise + 1),
    rural_gch_2018_l = as.numeric(as.character(rural_gch_2018_l)),
    #has_siblings = as.numeric(as.character(has_siblings)),
    parent = as.numeric(as.character(parent)),
    partner = as.numeric(as.character(partner)),
    born_nz = as.numeric(as.character(born_nz)),
    employed = as.numeric(as.character(employed)),
    hlth_disability = as.numeric(as.character(hlth_disability)
    )
  ) |>
  dplyr::mutate(sample_origin = sample_origin_names_combined) |>  #shorter name
  mutate(
    #initialize 'censored'
    censored = ifelse(lead(year_measured) == 1, 1, 0),
    
    # modify 'censored' based on the condition; no need to check for NA here as 'censored' is already defined in the previous step
    censored =  ifelse(is.na(censored) &
                         year_measured == 1, 1, censored)
    
    # # Apply the case_when condition for setting 'censored' based on 'wave' and the dynamic column specified by 'nzavs_exposure'
    # censored = case_when(
    #   # Add this condition to keep previous modifications unless the specific condition is met!is.na(censored) ~ censored,
    #
    #   # Then check if 'wave' is 2019 and the specified exposure is NA, adjusting the condition to reflect the accurate logic
    #   wave == 2019 & !is.na(!!sym(nzavs_exposure)) ~ 1,
    #
    #   # Default case if none of the above apply; might not be necessary if all possibilities are covered
    #   TRUE ~ 0
    # )
  )|>
  select(
    "id",
    "wave",
    "male",
    "age",
    "sample_origin",
    "sample_frame_opt_in", # opted in
    "censored",
    "alert_level_combined",
    # "edu",
    "education_level_coarsen",
    "born_nz",
    "rural_gch_2018_l", # rural urban units
    "hlth_disability",
    "kessler_latent_anxiety",
    "kessler_latent_depression",
    "eth_cat",
    #factor(EthCat, labels = c("Euro", "Maori", "Pacific", "Asian")),
    "employed",
    # Are you currently employed? (this includes self-employment or casual work)
    "household_inc_log",
    # Please estimate your total household income (before tax) for the last year.
    "nz_dep2018",
    # see nzavs materials
    "nzsei_13_l", # longtigudinal status (less imissingness)
    # see nzavs materials
    "partner",
    # 0 = no, 1 = yes
    "parent",
    # 0 = no, 1 = yes
    "political_conservative",
    #Please rate how politically liberal versus conservative you see yourself as being.
    "pol_wing",
    # Please rate how politically left-wing versus right-wing you see yourself as being.
    # see NZAVS,
    #  "has_siblings",
    #Do you have siblings?
    # sum siblings
    "hours_children_log",
    #Hours - Looking after children
    "hours_work_log",
    #Hours - Working in paid employment
    "hours_housework_log",
    # Hours - Housework/cooking
    "hours_exercise_log",
    "agreeableness",
    # Mini-IPIP6 Agreeableness (also modelled as empathy facet)
    # Sympathize with others' feelings.
    # Am not interested in other people's problems.
    # Feel others' emotions.
    # Am not really interested in others.
    "conscientiousness",
    # see mini ipip6
    # Get chores done right away.
    # Like order.
    # Make a mess of things.
    # Often forget to put things back in their proper place.
    "extraversion",
    # Mini-IPIP6 Extraversion
    # Am the life of the party.
    # Don't talk a lot.
    # Keep in the background.
    # Talk to a lot of different people at parties.
    "honesty_humility",
    # see mini ipip6
    # Would like to be seen driving around in a very expensive car.
    # Would get a lot of pleasure from owning expensive luxury goods.
    # Feel entitled to more of everything.
    # Deserve more things in life.
    "openness",
    # see mini ipip6
    # Have a vivid imagination.
    # Have difficulty understanding abstract ideas.
    # Do not have a good imagination.
    # Am not interested in abstract ideas.
    "neuroticism",
    # see mini ipip6
    # Have frequent mood swings.
    # Am relaxed most of the time.
    # Get upset easily.
    # Seldom feel blue.
    "modesty",
    # for unmeasured confounders
    "hlth_fatigue",
    "hlth_sleep_hours",
    "support",
    #   "support_help",
    # 'There are people I can depend on to help me if I really need it.
    #  "support_turnto",
    # There is no one I can turn to for guidance in times of stress.
    #  "support_noguidance_reverseed",
    #There is no one I can turn to for guidance in times of stress.
    "belong",
    #    "belong_accept",
    #Know that people in my life accept and value me.
    #    "belong_routside_reversed",
    # Feel like an outsider.
    #  "belong_beliefs",
    # Know that people around me share my attitudes and beliefs.
    "religion_identification_level",
    #How important is your religion to how you see yourself?"
    # "religion_church_binary",  not in 2019
    # How many times did you attend a church or place of worship in the last month?
    # "religion_spiritual_identification", not in 2019
    #w8,w10,w12-13 "I identify as a spiritual person."
    "env_climate_chg_real",
    "env_climate_chg_cause",
    # "Climate change is caused by humans"
    "env_climate_chg_concern",
    "env_sat_nz_environment",
    "trust_science_high_confidence_scientific_community",
    "trust_science_our_society_places_too_much_emphasis_reversed",
    "sample_weights"
  ) |>
  droplevels() |>
  data.frame()

#naniar::vis_miss(dat_long, warn_large_data = FALSE)
#dev.off()
# get n -------------------------------------------------------------------


n_participants<- n_unique(dat_long$id) # 32737
n_participants

table1::table1(~ trust_science_our_society_places_too_much_emphasis_reversed|wave, data = dat_long)
# save N for manus
# install margot package if needed

# devtools::install_github("go-bayes/margot")
here_save(n_participants, "n_participants")



# set baseline exposure outcome -------------------------------------------

# initial names without id, wave, etc
dat_long_names <- sort( colnames(dat_long) )

# check
dat_long_names


# exposure
exposure_var = c("trust_science_our_society_places_too_much_emphasis_reversed",
                 "censored")


outcome_vars <- c("env_climate_chg_real", "env_climate_chg_cause","env_climate_chg_concern", "env_sat_nz_environment")


exposure_var

# just core baseline variables
baseline_vars <-
  setdiff(dat_long_names, c("id", "wave", exposure_var, outcome_vars))


baseline_vars

# for tables
base_var <-
  setdiff(baseline_vars, c(outcome_vars, "sample_weights"))
base_var <- sort(base_var)
base_var



# positivity --------------------------------------------------------------

dt_positivity_full <- dat_long|>
  filter(wave == 2018| wave == 2019) |>
  select(wave, id, trust_science_our_society_places_too_much_emphasis_reversed)

dt_positivity_full


out <- msm::statetable.msm(trust_science_our_society_places_too_much_emphasis_reversed, id, data = dat_long)


# set lables if you like
# tab_labels <- c("< weekly", ">= weekly")

# transition table
transition_table  <- margot::transition_table(out)


# view
transition_table

# for import later
here_save(transition_table, "transition_table")



#  sd values --------------------------------------------------------------



# sd values ---------------------------------------------------------------

dt_outcome <-
  dat_long |>
  filter(wave == 2021)


sd_env_climate_chg_real <-
  sd(dt_outcome$env_climate_chg_real, na.rm = TRUE)
sd_env_climate_chg_cause <-
  sd(dt_outcome$env_climate_chg_cause, na.rm = TRUE)

sd_env_climate_chg_concern <-
  sd(dt_outcome$env_climate_chg_cause, na.rm = TRUE)

sd_env_sat_nz_environment <-
  sd(dt_outcome$env_sat_nz_environment, na.rm = TRUE)


# save for manuscript
here_save(sd_env_climate_chg_real, "sd_env_climate_chg_real")
here_save(sd_env_climate_chg_cause, "sd_env_climate_chg_cause")
here_save(sd_env_climate_chg_concern, "sd_env_climate_chg_concern")
here_save(sd_env_climate_chg_concern, "sd_env_climate_chg_concern")


# ordinary regressions ----------------------------------------------------
outcome_vars
dt_19 <-
  dat_long |>
  select(-censored) |> # factor has only one level in 2018
  filter(wave == 2019)


# check
base_var
# trust_

str(dt_19)
table(dt_19$sample_frame_opt_in)

naniar::vis_miss(dt_19, warn_large_data = FALSE)
base_var

#regress_var <- setdiff(base_var, "trust_science_high_confidence_scientific_community")
#regress_var
# base_vars set above
fit_env_climate_chg_real <-
  margot::regress_with_covariates(
    dt_19,
    outcome = "env_climate_chg_real",
    exposure = "trust_science_our_society_places_too_much_emphasis_reversed",
    baseline_vars = base_var
  )
parameters::model_parameters(fit_env_climate_chg_real,  ci_method="wald")[2, ]



fit_env_climate_chg_cause <-
  margot::regress_with_covariates(
    dt_19,
    outcome = "env_climate_chg_cause",
    exposure = "trust_science_our_society_places_too_much_emphasis_reversed",
    baseline_vars = base_var
  )
parameters::model_parameters(fit_env_climate_chg_cause,  ci_method="wald")[2, ]


fit_env_climate_chg_concern <-
  margot::regress_with_covariates(
    dt_19,
    outcome = "env_climate_chg_concern",
    exposure = "trust_science_our_society_places_too_much_emphasis_reversed",
    baseline_vars = base_var
  )
parameters::model_parameters(fit_env_climate_chg_concern,  ci_method="wald")[2, ]


fit_env_sat_nz_environment <-
  margot::regress_with_covariates(
    dt_19,
    outcome = "env_sat_nz_environment",
    exposure = "trust_science_our_society_places_too_much_emphasis_reversed",
    baseline_vars = base_var
  )
parameters::model_parameters(fit_env_sat_nz_environment,  ci_method="wald")[2, ]





# histograms --------------------------------------------------------------
# generate bar plot
# graph_density_of_exposure <-
#   coloured_histogram(dt_20,
#                      col_name = "trust_science_our_society_places_too_much_emphasis_reversed",
#                      scale_min = 1,
#                      scale_max = 7)


# histogram exposure ------------------------------------------------------

dt_19 <- dat_long |>
  filter(wave == 2019)

library(ggplot2)
library(dplyr)
#
mean_exposure <-mean(dt_19$trust_science_our_society_places_too_much_emphasis_reversed, na.rm=TRUE)

median_exposure <-median(dt_19$trust_science_our_society_places_too_much_emphasis_reversed, na.rm=TRUE)
median_exposure
mean_exposure

# # generate bar plot
graph_density_of_exposure_up <- margot::coloured_histogram_shift(
  dt_19,
  col_name = "trust_science_our_society_places_too_much_emphasis_reversed",
  binwidth = 1,
  range_highlight = c(0,mean_exposure)
)
graph_density_of_exposure_up


here_save(graph_density_of_exposure_up, "graph_density_of_exposure_up")




# tables ------------------------------------------------------------------

library(gtsummary)



# table baseline ----------------------------------------------------------

# get names
base_var


# prepare df
selected_base_cols <-
  dt_19 |> select(all_of(base_var)) #


#check
colnames(selected_base_cols)

# tabls
library(gtsummary)

table_baseline <- selected_base_cols |>
  janitor::clean_names(case = "title") |>
  tbl_summary(
    missing = "ifany",
    percent = "column",
    statistic = list(
      all_continuous() ~ c(
        "{mean} ({sd})", # Mean and SD
        "{min}, {max}", # Range (Min, Max)
        "{p25}, {p75}" # IQR (25th percentile, 75th percentile)
      )
    ),
    type = all_continuous() ~ "continuous2"
  ) |>
  modify_header(label = "**Exposure + Demographic Variables**") |> # update the column header
  bold_labels()



table_baseline
# save baseline
here_save(table_baseline, "table_baseline")
#

# table exposure ----------------------------------------------------------

# get first and second wave
dt_exposure<- dat_long|>
  dplyr::filter(wave == 2019 | wave == 2020) |>
  droplevels()

exposure_var
# get vars.
selected_exposure_cols <-
  dt_exposure %>% select(
    c(
      "trust_science_our_society_places_too_much_emphasis_reversed",
      "wave"
    )
  )

# check
str(selected_exposure_cols)


library(gtsummary)

table_exposures <- selected_exposure_cols %>%
  janitor::clean_names(case = "title") %>%
  labelled::to_factor() %>%  # ensure consistent use of pipe operator
  tbl_summary(
    by = "Wave",  #specify the grouping variable. Adjust "Wave" to match the cleaned column name
    missing = "always",
    percent = "column",
    # statistic = list(all_continuous() ~ "{mean} ({sd})")  # Uncomment and adjust if needed for continuous variables
  ) %>%
  #  add_n() %>%  # Add column with total number of non-missing observations
  modify_header(label = "**Exposure Variables by Wave**") %>%  # Update the column header
  bold_labels()

table_exposures


# save baseline
here_save(table_exposures, "table_exposures")

table_exposures


# outcome table -----------------------------------------------------------
dt_outcomes <- dat_long|>
  dplyr::filter(wave == 2019 | wave == 2021) |>
  droplevels()

names_outcomes_tab <- setdiff(outcome_vars, dt_outcomes)
names_outcomes_sorted <- sort(names_outcomes_tab)
names_outcomes_final <- names_outcomes_sorted # consistent workflow

names_outcomes_final

names_outcomes_final

base_var

# names_outcomes_final
# better names
selected_outcome_cols <-
  dt_outcomes %>% select(all_of(names_outcomes_final),
                         wave) #|>
#mutate(Volunteers_binary = factor(ifelse(hours_charity > 0, "yes", "no"),
#levels = c("no", "yes"))) #%>% rename vars if desired
#   rename(
#     Social_belonging = belong,
#     Annual_charity = charity_donate,
#     Volunteering_hours = hours_charity,
#     Community_gives_money = community_money_binary,
#     Community_gives_time = community_time_binary,
#     Family_gives_money = family_money_binary,
#     Family_gives_time = family_time_binary,
#     Friends_give_money = friends_money_binary,
#     Friends_give_time = friends_time_binary,
#     Social_support = support,
#     Sense_neighbourhood_community = neighbourhood_community
#   )
#
# # order names correctly
selected_outcome_cols <- selected_outcome_cols %>%
  select(sort(names(selected_outcome_cols)))

# checks
str(selected_outcome_cols)
colnames(selected_outcome_cols)

table_outcomes <- selected_outcome_cols %>%
  janitor::clean_names(case = "title") %>%
  labelled::to_factor() %>%  # ensure consistent use of pipe operator
  tbl_summary(
    by = "Wave",  #specify the grouping variable. Adjust "Wave" to match the cleaned column name
    missing = "always",
    percent = "column",
    # statistic = list(all_continuous() ~ "{mean} ({sd})")  # Uncomment and adjust if needed for continuous variables
  ) %>%
  #  add_n() %>%  # Add column with total number of non-missing observations
  modify_header(label = "**Outcome Variables by Wave**") %>%  # Update the column header
  bold_labels()

table_outcomes


# save
here_save(table_outcomes, "table_outcomes")



# impute ------------------------------------------------------------------

naniar::vis_miss(dat_long, warn_large_data = FALSE)


# missing values data
prep_coop_all <- margot::margot_wide_impute_baseline(
  dat_long,
  baseline_vars = baseline_vars,
  exposure_var = exposure_var,
  outcome_vars = outcome_vars
)

# check
push_mods

# save function -- will save to your "push_mod" directory
here_save(prep_coop_all, "prep_coop_all")

# read function
prep_coop_all <- here_read("prep_coop_all")

head(prep_coop_all)

#check
naniar::vis_miss(prep_coop_all, warn_large_data = FALSE)

# arrange data for analysis -----------------------------------------------
# spit and shine

df_wide_censored <- prep_coop_all |>
  select(-id) |> 
  mutate(
    t0_eth_cat = as.factor(t0_eth_cat),
    t0_education_level_coarsen = as.factor(t0_education_level_coarsen),
    t0_sample_frame_opt_in = as.factor(t0_sample_frame_opt_in),
  ) |>
  relocate("t0_censored", .before = starts_with("t1_")) |>
  relocate("t1_censored", .before = starts_with("t2_")) |>
  relocate(starts_with("t0_"), .before = starts_with("t1_")) |>
  relocate("t0_censored", .before = starts_with("t1_"))  |>
  relocate("t1_censored", .before = starts_with("t2_"))

# check
naniar::vis_miss(df_wide_censored, warn_large_data = FALSE)

str(df_wide_censored)
# Assuming df_wide_censored is your dataframe

# Calculate the conditions before the mutate steps
t0_na_condition <-  # Won't be needed,
  rowSums(is.na(select(df_wide_censored, starts_with("t1_")))) > 0
t1_na_condition <-
  rowSums(is.na(select(df_wide_censored, starts_with("t2_")))) > 0

df_clean <- df_wide_censored %>%
  mutate(t0_censored = ifelse(t0_na_condition, 0, t0_censored)) %>%
  mutate(t1_censored = ifelse(t1_na_condition, 0, t1_censored)) %>%
  mutate(across(starts_with("t1_"), ~ ifelse(t0_censored == 0, NA_real_, .)),
         across(starts_with("t2_"), ~ ifelse(t0_censored == 0, NA_real_, .))) %>%
  mutate(across(starts_with("t2_"), ~ ifelse(t1_censored == 0, NA_real_, .))) |>
  # select variables
  dplyr::mutate(
    across(
      .cols = where(is.numeric) &
        !ends_with("_censored") &
        !t0_sample_weights &
        !t0_alert_level_combined &
        !t0_nzsei_13_l &
        !t0_rural_gch_2018_l &
        !t0_trust_science_our_society_places_too_much_emphasis_reversed &
        !starts_with("t1_"),
      #  !t2_hours_charity,
      .fns = ~ scale(.),
      .names = "{.col}_z"
    )
  ) |>
  # select(-t0_charity_donate,
  #        -t0_hours_charity) |>
  select(
    where(is.factor),
    ends_with("_censored"),
    t0_sample_weights,
    t0_nzsei_13_l,
    t0_rural_gch_2018_l,
    t0_alert_level_combined,
    t0_trust_science_our_society_places_too_much_emphasis_reversed,
    starts_with("t1_"),
    ends_with("_z")
  ) |>
  dplyr::select(-t0_alert_level_combined) %>%  # not needed
  relocate(starts_with("t0_"), .before = starts_with("t1_")) |>
  relocate("t0_censored", .before = starts_with("t1_"))  |>
  relocate("t1_censored", .before = starts_with("t2_"))

# save
here_save(df_clean, "df_clean")

# read
df_clean<-here_read("df_clean")


# checks
head(df_clean)
str(df_clean)
# vis missing
naniar::vis_miss(df_clean, warn_large_data = FALSE)

table(df_clean$t1_trust_science_our_society_places_too_much_emphasis_reversed)

# check path
push_mods



# start analysis here -----------------------------------------------------
# read if needed
df_clean <- here_read("df_clean")

# check
df_clean$t0_sample_weights
# get names
names_base <-
  df_clean |> select(starts_with("t0"),
                     -t0_sample_weights,
                    # -t0_alert_level_combined,
                     -t0_censored
                    ) |> colnames()

names_base

# save
here_save(names_base, "names_base")

# if needed
names_base <- here_read("names_base")


#### SET VARIABLE NAMES: Customise for each outcomewide model
#  model

# define exposure (not we have the baseline exposure)
A <- c( "t1_trust_science_our_society_places_too_much_emphasis_reversed")

A
C <- c("t1_censored")

#naniar::vis_miss(prep_coop_all, warn_large_data = F)
# L <- list(c(NULL), c("L2"))
# W <- c(paste(names_base, collapse = ", "))


# define shift function (if any one is average make them average, otherwise leave alone)
shift_up <- function(data, trt) {
  ifelse(data[[trt]] < 6, 6,  data[[trt]])
}


median(df_clean$t1_trust_science_our_society_places_too_much_emphasis_reversed, na.rm=TRUE)
# test model --------------------------------------------------------------

df_clean_slice <- df_clean |>
  slice_head(n = 500) |>
  as.data.frame()

A

df_clean$t1_trust_science_our_society_places_too_much_emphasis_reversed
t2_env_climate_chg_real_z_test <- lmtp_tmle(
  data = df_clean_slice,
  trt = A,
  baseline = names_base,
  outcome = "t2_env_climate_chg_real_z",
  cens = C,
  shift = shift_up,
  outcome_type = "continuous",
  mtp = TRUE,
  folds = 5,
  weights = df_clean_slice$t0_sample_weights,
  learners_trt = sl_lib,
  learners_outcome = sl_lib,
  parallel = n_cores
)

# print
t2_env_climate_chg_real_z_test

A
# models start here -------------------------------------------------------

A
# outcome wide models
t2_env_climate_chg_real_z <- lmtp_tmle(
  data = df_clean,
  trt = A,
  baseline = names_base,
  outcome = "t2_env_climate_chg_real_z",
  cens = C,
  shift = shift_up,
  outcome_type = "continuous",
  mtp = TRUE,
  folds = 5,
  weights = df_clean$t0_sample_weights,
  learners_trt= sl_lib,
  learners_outcome= sl_lib,
  parallel = n_cores
)
here_save(t2_env_climate_chg_real_z, "t2_env_climate_chg_real_z")


t2_env_climate_chg_real_z_null <- lmtp_tmle(
  data = df_clean,
  trt = A,
  baseline = names_base,
  outcome = "t2_env_climate_chg_real_z",
  cens = C,
  shift = NULL,
  outcome_type = "continuous",
  mtp = TRUE,
  folds = 5,
  weights = df_clean$t0_sample_weights,
  learners_trt= sl_lib,
  learners_outcome= sl_lib,
  parallel = n_cores
)
here_save(t2_env_climate_chg_real_z_null, "t2_env_climate_chg_real_z_null")


# model 2 t2_env_climate_chg_cause_z
t2_env_climate_chg_cause_z <- lmtp_tmle(
  data = df_clean,
  trt = A,
  baseline = names_base,
  outcome = "t2_env_climate_chg_cause_z",
  cens = C,
  shift = shift_up,
  outcome_type = "continuous",
  mtp = TRUE,
  folds = 5,
  weights = df_clean$t0_sample_weights,
  learners_trt= sl_lib,
  learners_outcome= sl_lib,
  parallel = n_cores
)

here_save(t2_env_climate_chg_cause_z, "t2_env_climate_chg_cause_z")


t2_env_climate_chg_cause_z_null  <- lmtp_tmle(
  data = df_clean,
  trt = A,
  baseline = names_base,
  outcome = "t2_env_climate_chg_cause_z",
  cens = C,
  shift = NULL,
  outcome_type = "continuous",
  mtp = TRUE,
  folds = 5,
  weights = df_clean$t0_sample_weights,
  learners_trt= sl_lib,
  learners_outcome= sl_lib,
  parallel = n_cores
)

here_save(t2_env_climate_chg_cause_z_null, "t2_env_climate_chg_cause_z_null")


# model 3  t2_env_climate_chg_concern_z

t2_env_climate_chg_concern_z <- lmtp_tmle(
  data = df_clean,
  trt = A,
  baseline = names_base,
  outcome = "t2_env_climate_chg_concern_z",
  cens = C,
  shift = shift_up,
  outcome_type = "continuous",
  mtp = TRUE,
  folds = 5,
  weights = df_clean$t0_sample_weights,
  learners_trt= sl_lib,
  learners_outcome= sl_lib,
  parallel = n_cores
)
here_save(t2_env_climate_chg_concern_z, "t2_env_climate_chg_concern_z")


t2_env_climate_chg_concern_z_null  <- lmtp_tmle(
  data = df_clean,
  trt = A,
  baseline = names_base,
  outcome = "t2_env_climate_chg_concern_z",
  cens = C,
  shift = NULL,
  outcome_type = "continuous",
  mtp = TRUE,
  folds = 5,
  weights = df_clean$t0_sample_weights,
  learners_trt= sl_lib,
  learners_outcome= sl_lib,
  parallel = n_cores
)

t2_env_climate_chg_concern_z_null
here_save(t2_env_climate_chg_concern_z_null, "t2_env_climate_chg_concern_z_null")


# model 4 t2_env_sat_nz_environment_z

t2_env_sat_nz_environment_z <- lmtp_tmle(
  data = df_clean,
  trt = A,
  baseline = names_base,
  outcome = "t2_env_sat_nz_environment_z",
  cens = C,
  shift = shift_up,
  outcome_type = "continuous",
  mtp = TRUE,
  folds = 5,
  weights = df_clean$t0_sample_weights,
  learners_trt= sl_lib,
  learners_outcome= sl_lib,
  parallel = n_cores
)
here_save(t2_env_sat_nz_environment_z, "t2_env_sat_nz_environment_z")

t2_env_sat_nz_environment_z_null  <- lmtp_tmle(
  data = df_clean,
  trt = A,
  baseline = names_base,
  outcome = "t2_env_sat_nz_environment_z",
  cens = C,
  shift = NULL,
  outcome_type = "continuous",
  mtp = TRUE,
  folds = 5,
  weights = df_clean$t0_sample_weights,
  learners_trt= sl_lib,
  learners_outcome= sl_lib,
  parallel = n_cores
)

t2_env_sat_nz_environment_z_null
here_save(t2_env_sat_nz_environment_z_null, "t2_env_sat_nz_environment_z_null")

# contrasts

# contrasts 1 climate change is real
t2_env_climate_chg_real_z <- here_read("t2_env_climate_chg_real_z")
t2_env_climate_chg_real_z_null <-
  here_read("t2_env_climate_chg_real_z_null")

# first contrast
contrast_t2_env_climate_chg_real_z <-
  lmtp_contrast(t2_env_climate_chg_real_z,
                ref = t2_env_climate_chg_real_z_null,
                type = "additive")


tab_contrast_t2_env_climate_chg_real_z <-
  margot_lmtp_tab(contrast_t2_env_climate_chg_real_z,
                  scale = "RD",
                  new_name = "Climate change is real.")

tab_contrast_t2_env_climate_chg_real_z


out_tab_contrast_t2_env_climate_chg_real_z<-
  lmtp_evalue_tab(tab_contrast_t2_env_climate_chg_real_z,
                  scale = c("RD"))

out_tab_contrast_t2_env_climate_chg_real_z



# contrasts: model 2 t2_env_climate_chg_cause_z
# climate change is real
t2_env_climate_chg_cause_z <- here_read("t2_env_climate_chg_cause_z")
t2_env_climate_chg_cause_z_null <-
  here_read("t2_env_climate_chg_cause_z_null")

# first contrast
contrast_t2_env_climate_chg_cause_z <-
  lmtp_contrast(t2_env_climate_chg_cause_z,
                ref = t2_env_climate_chg_cause_z_null,
                type = "additive")


tab_contrast_t2_env_climate_chg_cause_z <-
  margot_lmtp_tab(contrast_t2_env_climate_chg_cause_z,
                  scale = "RD",
                  new_name = "Climate change caused by humans")

tab_contrast_t2_env_climate_chg_cause_z


out_tab_contrast_t2_env_climate_chg_cause_z<-
  lmtp_evalue_tab(tab_contrast_t2_env_climate_chg_cause_z,
                  scale = c("RD"))

out_tab_contrast_t2_env_climate_chg_cause_z

#  contrasts: 3: t2_env_climate_chg_concern_z
t2_env_climate_chg_concern_z <- here_read("t2_env_climate_chg_concern_z")
t2_env_climate_chg_concern_z_null <-
  here_read("t2_env_climate_chg_concern_z_null")

# first contrast
contrast_t2_env_climate_chg_concern_z <-
  lmtp_contrast(t2_env_climate_chg_concern_z,
                ref = t2_env_climate_chg_concern_z_null,
                type = "additive")


tab_contrast_t2_env_climate_chg_concern_z <-
  margot_lmtp_tab(contrast_t2_env_climate_chg_concern_z,
                  scale = "RD",
                  new_name = "Deeply concerned about climate change")

tab_contrast_t2_env_climate_chg_concern_z


out_tab_contrast_t2_env_climate_chg_concern_z<-
  lmtp_evalue_tab(tab_contrast_t2_env_climate_chg_concern_z,
                  scale = c("RD"))

out_tab_contrast_t2_env_climate_chg_concern_z

# contrasts: 4 t2_env_sat_nz_environment_z

t2_env_sat_nz_environment_z <- here_read("t2_env_sat_nz_environment_z")
t2_env_sat_nz_environment_z_null <-
  here_read("t2_env_sat_nz_environment_z_null")


# first contrast
contrast_t2_env_sat_nz_environment_z <-
  lmtp_contrast(t2_env_sat_nz_environment_z,
                ref = t2_env_sat_nz_environment_z_null,
                type = "additive")


tab_contrast_t2_env_sat_nz_environment_z <-
  margot_lmtp_tab(contrast_t2_env_sat_nz_environment_z,
                  scale = "RD",
                  new_name = "Satisfied with NZ natural environment")

tab_contrast_t2_env_sat_nz_environment_z

out_tab_contrast_t2_env_sat_nz_environment_z<-
  lmtp_evalue_tab(tab_contrast_t2_env_sat_nz_environment_z,
                  scale = c("RD"))

out_tab_contrast_t2_env_sat_nz_environment_z


# bind individual tables
tab_envir <- rbind(
  out_tab_contrast_t2_env_climate_chg_real_z,
  out_tab_contrast_t2_env_climate_chg_cause_z,
  out_tab_contrast_t2_env_climate_chg_concern_z,
  out_tab_contrast_t2_env_sat_nz_environment_z
)
t2_env_climate_chg_real_z
out_tab_contrast_t2_env_climate_chg_cause_z
# make group table
group_tab_envir<- group_tab(tab_envir  , type = "RD")

# save
here_save(group_tab_envir, "group_tab_envir")

# read
group_tab_envir <- here_read("group_tab_envir")


# create plots -------------------------------------------------------------
# check N
n_participants

sub_title = "Trust in science (place in society): shift those below median to median, otherwise do not shift: N = 32,737"


# graph health
plot_group_tab_envir <- margot_plot(
  group_tab_envir,
  type = "RD",
  title = "Environmental attitudes",
  subtitle = sub_title,
  estimate_scale = 1,
  base_size = 12,
  text_size = 3.0,
  point_size = .5,
  title_size = 15,
  subtitle_size = 11,
  legend_text_size = 8,
  legend_title_size = 10,
  x_offset = -1,
  x_lim_lo = -1,
  x_lim_hi =  .5
)
plot_group_tab_envir
dev.off()


# save graph
ggsave(
  plot_group_tab_envir,
  path = here::here(here::here(push_mods, "figs")),
  width = 12,
  height = 8,
  units = "in",
  filename = "plot_group_tab_envir.png",
  device = 'png',
  limitsize = FALSE,
  dpi = 600
)

# second analysis Graph

