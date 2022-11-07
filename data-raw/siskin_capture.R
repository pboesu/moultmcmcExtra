## code to prepare `siskin_capture` dataset goes here
siskin_capture <- readRDS("../../2022_moult_methods_paper/data/empirical_capture_prob_siskins.rds")
usethis::use_data(siskin_capture, overwrite = TRUE)
