print("Loading data...")
ges <- readRDS("/Users/AlexanderWang/Desktop/Califano_Lab_Fall_2021/LNCaPF21_project/data/lncap_clones_R6/AJ081_R6/gExpr/AJ081_STEP02_gExpr_PISCES02_preProcess_modular/AJ081_gExpr_filt.rds")
net <- readRDS("/Users/AlexanderWang/Desktop/Califano_Lab_Fall_2021/Pyther_project/pyther_main/test/unit_tests/test_3/test_3_inputs/AJ081_pruned.rds")
print("Sampling...")
# set.seed(0)
# ran_rows <- sample(1:nrow(ges), 10000, replace = FALSE)
# ran_rows <- sample(1:nrow(ges), 2500, replace = FALSE)
j <- 1
ran_rows <- rep(NA, 25000)
for(i in 1:500){
  ran_rows[j:(j+49)] <- names(net[[i]]$tfmode)
  # print(length(names(net[[1]]$tfmode)))
  # print(j:j+49)
  j <- j + 50
}
ran_rows <- unique(ran_rows)[1:2000]
set.seed(0)
# ran_cols <- sample(1:ncol(ges), 250, replace = FALSE)
# ran_cols <- sample(1:ncol(ges), 100, replace = FALSE)
ran_cols <- sample(1:ncol(ges), 25, replace = FALSE)
ges_small <- ges[ran_rows, ran_cols]
print(dim(ges_small))
print("Saving data...")
saveRDS(ges_small, "/Users/AlexanderWang/Desktop/Califano_Lab_Fall_2021/Pyther_project/pyther_main/test/unit_tests/test_3/test_3_inputs/AJ081_gExpr_GES.rds")
print("Done.")

# setwd("/Users/AlexanderWang/Desktop/Califano_Lab_Fall_2021/LNCaPF21_project/LNCaP_project/job_outputs/")
# opt <- list()
# opt$file="/Users/AlexanderWang/Desktop/Califano_Lab_Fall_2021/Pyther_project/pyther_main/test/unit_tests/test_3/test_3_inputs/AJ081_gExpr_GES.rds"
# opt$ext="tsv"