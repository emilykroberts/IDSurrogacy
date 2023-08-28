library(here); library(xtable)
setwd("./R/results")

############################ STEP 2: Aggregate Results ############################
# post-processing to aggregate .txt files and recreate tables in manuscript
# to be executed after desired simulation replications are run
# (pre-run files can be found in Table2 folder) for reproducibility check

############################ Table 2 ############################
b = NULL
fil = "Table2" # indicates files are saved in for this particular code replication
for(i in 1:8){ # loop over 8 data generating scenarios to create Table 2
  # load in files for a given scenario 
  params = do.call(rbind, sapply(list.files(path=fil, pattern = paste0("params,", i), full.names=TRUE), read.table, head = TRUE, simplify = F))
  # create tabular form of results (mean and standard deviation of point estimates, mean of standard error estimates)
  a = (rbind(rbind(round(colMeans(params[c(
    "int", "slope", "shape12_0", "shape13_0", "shape23_0",
    "shape12_1", "shape13_1", "shape23_1",
    "scale12_0", "scale13_0", "scale23_0",
    "scale12_1", "scale13_1", "scale23_1",
    "theta23_0", "theta23_1"
  )], na.rm = T), 3)),
  rbind(round(colMeans(params[c(
    "intse", "slopese", "shape12_0SE", "shape13_0SE", "shape23_0SE",
    "shape12_1SE", "shape13_1SE", "shape23_1SE",
    "scale12_0SE", "scale13_0SE", "scale23_0SE",
    "scale12_1SE", "scale13_1SE", "scale23_1SE",
    "theta23_0SE", "theta23_1SE"
  )], na.rm = T), 3)),
  rbind(round(apply(params[c(
    "int", "slope", "shape12_0", "shape13_0", "shape23_0",
    "shape12_1", "shape13_1", "shape23_1",
    "scale12_0", "scale13_0", "scale23_0",
    "scale12_1", "scale13_1", "scale23_1",
    "theta23_0", "theta23_1"
  )], 2, FUN = sd, na.rm = T), 3))))
  b = rbind(b, a)
}

# save table results table2results.csv
write.csv(b, file = "table2results.csv") 

############################ Table S2 ############################

setwd("./supplement") # results about supplemental simulations in supplement folder

# read in results from many simulations for rows 1-3 (scenario 2 in this case as indicated by the pattern argument below)
params = do.call(rbind, sapply(list.files(path = "TableS2_1_3", pattern = "params,2", full.names = TRUE), read.table, head = TRUE, simplify = F))

# load in files for a given scenario and create tabular form of results (mean and standard deviation of point estimates, mean of standard error estimates)
# produces Latex code for manuscript tables that corresponds to Supplemental Table 2 rows 1-3
xtable(rbind(rbind(round(colMeans(params[c(
  "int", "slope", 
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1"
)], na.rm = T), 3)), 
rbind(round(colMeans(params[c(
  "intse", "slopese", 
  "scale12_0SE", "scale13_0SE", "scale23_0SE", 
  "scale12_1SE", "scale13_1SE", "scale23_1SE"
)], na.rm = T), 3)), 
rbind(round(apply(params[c(
  "int", "slope",
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1"
)], 2, FUN = sd, na.rm = T), 3))), digits = 3)

tabS2a = (rbind(rbind(round(colMeans(params[c(
  "int", "slope", "shape12_0", "shape13_0", "shape23_0", 
  "shape12_1", "shape13_1", "shape23_1", 
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1", 
  "theta23_0", "theta23_1"
)], na.rm = T), 3)), 
rbind(round(colMeans(params[c(
  "intse", "slopese", "shape12_0SE", "shape13_0SE", "shape23_0SE", 
  "shape12_1SE", "shape13_1SE", "shape23_1SE", 
  "scale12_0SE", "scale13_0SE", "scale23_0SE", 
  "scale12_1SE", "scale13_1SE", "scale23_1SE", 
  "theta23_0SE", "theta23_1SE"
)], na.rm = T), 3)), 
rbind(round(apply(params[c(
  "int", "slope", "shape12_0", "shape13_0", "shape23_0", 
  "shape12_1", "shape13_1", "shape23_1", 
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1", 
  "theta23_0", "theta23_1"
)], 2, FUN = sd, na.rm = T), 3))))

# read in results from many simulations for rows 4-6 (scenario 2 in this case as indicated by the pattern argument below)
params = do.call(rbind, sapply(list.files(path = "TableS2_4_6", pattern = "params,2", full.names = TRUE), read.table, head = TRUE, simplify = F))

# produces Latex code for manuscript tables, corresponds to Supplemental Table 2 rows 4-6 of results
xtable(rbind(rbind(round(colMeans(params[c(
  "int", "slope", 
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1"
)], na.rm = T), 3)), 
rbind(round(colMeans(params[c(
  "intse", "slopese", 
  "scale12_0SE", "scale13_0SE", "scale23_0SE", 
  "scale12_1SE", "scale13_1SE", "scale23_1SE"
)], na.rm = T), 3)), 
rbind(round(apply(params[c(
  "int", "slope",
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1"
)], 2, FUN = sd, na.rm = T), 3))), digits = 3)

tabS2b = (rbind(rbind(round(colMeans(params[c(
  "int", "slope", "shape12_0", "shape13_0", "shape23_0", 
  "shape12_1", "shape13_1", "shape23_1", 
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1", 
  "theta23_0", "theta23_1"
)], na.rm = T), 3)), 
rbind(round(colMeans(params[c(
  "intse", "slopese", "shape12_0SE", "shape13_0SE", "shape23_0SE", 
  "shape12_1SE", "shape13_1SE", "shape23_1SE", 
  "scale12_0SE", "scale13_0SE", "scale23_0SE", 
  "scale12_1SE", "scale13_1SE", "scale23_1SE", 
  "theta23_0SE", "theta23_1SE"
)], na.rm = T), 3)), 
rbind(round(apply(params[c(
  "int", "slope", "shape12_0", "shape13_0", "shape23_0", 
  "shape12_1", "shape13_1", "shape23_1", 
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1", 
  "theta23_0", "theta23_1"
)], 2, FUN = sd, na.rm = T), 3))))

tabS2 = rbind(tabS2a, tabS2b)
# save csv file for table S2 (both rows 1-3 and 4-6 combined)
write.csv(tabS2, file = "tableS2results.csv")

############################ Table S3 ############################

# read in results from many simulations (scenario 2 in this case as indicated by the pattern argument below)
params = do.call(rbind, sapply(list.files(path = "TableS3", pattern = "params,2", full.names = TRUE), read.table, head = TRUE, simplify = F))

# load in files for a given scenario and create tabular form of results (mean and standard deviation of point estimates, mean of standard error estimates)
# produces Latex code for manuscript tables, corresponds to Supplemental Table 3
xtable(rbind(rbind(round(colMeans(params[c(
  "int", "slope", "shape12_0", "shape13_0", "shape23_0", 
  "shape12_1", "shape13_1", "shape23_1", 
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1", 
  "theta23_0", "theta23_1"
)], na.rm = T), 3)), 
rbind(round(colMeans(params[c(
  "intse", "slopese", "shape12_0SE", "shape13_0SE", "shape23_0SE", 
  "shape12_1SE", "shape13_1SE", "shape23_1SE", 
  "scale12_0SE", "scale13_0SE", "scale23_0SE", 
  "scale12_1SE", "scale13_1SE", "scale23_1SE", 
  "theta23_0SE", "theta23_1SE"
)], na.rm = T), 3)), 
rbind(round(apply(params[c(
  "int", "slope", "shape12_0", "shape13_0", "shape23_0", 
  "shape12_1", "shape13_1", "shape23_1", 
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1", 
  "theta23_0", "theta23_1"
)], 2, FUN = sd, na.rm = T), 3))), digits = 3)

tabS3 = (rbind(rbind(round(colMeans(params[c(
  "int", "slope", "shape12_0", "shape13_0", "shape23_0", 
  "shape12_1", "shape13_1", "shape23_1", 
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1", 
  "theta23_0", "theta23_1"
)], na.rm = T), 3)), 
rbind(round(colMeans(params[c(
  "intse", "slopese", "shape12_0SE", "shape13_0SE", "shape23_0SE", 
  "shape12_1SE", "shape13_1SE", "shape23_1SE", 
  "scale12_0SE", "scale13_0SE", "scale23_0SE", 
  "scale12_1SE", "scale13_1SE", "scale23_1SE", 
  "theta23_0SE", "theta23_1SE"
)], na.rm = T), 3)), 
rbind(round(apply(params[c(
  "int", "slope", "shape12_0", "shape13_0", "shape23_0", 
  "shape12_1", "shape13_1", "shape23_1", 
  "scale12_0", "scale13_0", "scale23_0", 
  "scale12_1", "scale13_1", "scale23_1", 
  "theta23_0", "theta23_1"
)], 2, FUN = sd, na.rm = T), 3)))
)
# save as .csv for Table S3
write.csv(tabS3, file = "tableS3results.csv")

############################ Table S4 ############################

b = NULL
for(i in 1:8){ # for loop over 8 data generating scenarios
  
  params = do.call(rbind, sapply(list.files(path = "TableS4", pattern = paste0("params,", i), full.names=TRUE), read.table, head = TRUE, simplify = F))
  
  # load in files for a given scenario and create tabular form of results (mean and standard deviation of point estimates, mean of standard error estimates)
  a = (rbind(rbind(round(colMeans(params[c(
    "int", "slope", "shape12_0", "shape13_0", "shape23_0",
    "shape12_1", "shape13_1", "shape23_1",
    "scale12_0", "scale13_0", "scale23_0",
    "scale12_1", "scale13_1", "scale23_1",
    "theta23_0", "theta23_1"
  )], na.rm = T), 3)),
  rbind(round(colMeans(params[c(
    "intse", "slopese", "shape12_0SE", "shape13_0SE", "shape23_0SE",
    "shape12_1SE", "shape13_1SE", "shape23_1SE",
    "scale12_0SE", "scale13_0SE", "scale23_0SE",
    "scale12_1SE", "scale13_1SE", "scale23_1SE",
    "theta23_0SE", "theta23_1SE"
  )], na.rm = T), 3)),
  rbind(round(apply(params[c(
    "int", "slope", "shape12_0", "shape13_0", "shape23_0",
    "shape12_1", "shape13_1", "shape23_1",
    "scale12_0", "scale13_0", "scale23_0",
    "scale12_1", "scale13_1", "scale23_1",
    "theta23_0", "theta23_1"
  )], 2, FUN = sd, na.rm = T), 3))))
  
  b = rbind(b, a)
  
}
# save as .csv for Table S4
write.csv(b, file = "tableS4results.csv")


