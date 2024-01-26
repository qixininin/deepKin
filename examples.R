install.packages("/Users/zqx/ZJU-PhD/2-Work/Rpackage/deepKin_0.1.0.tar.gz", repos = NULL, type = "source")
library(deepKin)


#### Simple examples -----------------------------------------------------------
## UKB Chinese
n = 1435
me = 38948.2

## Maize
n = 207
me = 19995.37

## Wheat
n = 811
me = 281.86

## Cattle
n = 1490
me = 242.73


npairs = n*(n-1)/2
minMe = minMe(theta = (1/2)^(0:4), alpha = 0.05/npairs, beta = 0.1)
deepDegree = log(deepTheta(me = me, alpha = 0.05/npairs, beta = 0.1), base = 1/2)
thrd = -logpThreshold(me = me, (1/2)^seq(0,floor(deepDegree)), beta = 0.1)



#### Example One: one cohort ---------------------------------------------------
#### ++ Prepare ----------------------------------------------------------------
plink_path = "/Users/zqx/Downloads/plink2"
# gear_path = "/Users/zqx/Downloads/gear"
bfileprefix = "./inst/1KG-EUR.example1"

#### ++ Step 1 -----------------------------------------------------------------
## deepKin Principle I:
## Evaluate the minimal number of me required for detecting target relatedness
n = as.numeric(system(paste0("awk 'END{print NR}' ", bfileprefix, ".fam"), intern = T))
m = as.numeric(system(paste0("awk 'END{print NR}' ", bfileprefix, ".bim"), intern = T))
id = read.table(file = paste0(bfileprefix,".fam"), header = F)[,1]
npairs = n*(n-1)/2
degree = 0:5
theta = (1/2)^degree
minMe = minMe(theta = theta, alpha = 0.05/npairs, beta = 0.1)

#### ++ Step 2 -----------------------------------------------------------------
## Perform your data QC and calculate me
## GRM method or LB method
me = calculateMe(bfileprefix = bfileprefix, method = "GRM", plink_path = plink_path, pop_size = n)  # not suggested for biobank data
# me = calculateMe(bfileprefix = bfileprefix, method = "LB", gear_path = gear_path, pop_size = n)     # suggested for biobank data

#### ++ Step 3 -----------------------------------------------------------------
## deepKin Principle II:
deeptheta = deepTheta(me = me, alpha = 0.05/npairs, beta = 0.1)
deepDegree = log(deeptheta, base = 1/2)
## deepKin Principle III:
thrd = -logpThreshold(me = me, (1/2)^seq(0,floor(deepDegree)), beta = 0.1)
## Set target degree
targetDegree = floor(deepDegree)
targetTheta = (1/2)^targetDegree

#### ++ Step 4  -------------------------------------------------------------------
## Perform plink GRM
if(!file.exists(paste0(bfileprefix,".rel.bin"))){
  system(paste0(plink_path, " --silent --bfile ", bfileprefix, " --make-rel triangle bin4 --out ", bfileprefix))
}
## Extract plink GRM
grm.rst = extractPlinkGRM(bfileprefix, xcohort = F, pop_size = n)
grm.diag = grm.rst$diag
grm.tri  = grm.rst$tri

## Perform deepKin
deepkin = deepKin_estimation(grm.diag = grm.diag, grm.tri = grm.tri, xcohort = F, me = me)

#### ++ Step 5 --------------------------------------------------------------------
## deepkin inference
deepkin.qc = deepkin[which(deepkin$`-logp` > (-logpThreshold(me = me, theta = targetTheta, beta = 0.1))), ]
deepkin.qc$tag = cut(deepkin.qc$`-logp`,
                     breaks = c(0,thrd,Inf),
                     labels = c("Unrelated", paste0("Degree ", (length(thrd)-1):0)))
index = as.numeric(rownames(deepkin.qc))
deepkin.qc.id = extract_individual_id(id = id, index = index, xcohort = F)
deepkin.qc.rst = cbind(deepkin.qc.id, deepkin.qc)


#### Example Two: cross cohorts ------------------------------------------------
#### ++ Prepare ----------------------------------------------------------------
plink_path = "/Users/zqx/Downloads/plink2"
# gear_path = "/Users/zqx/Downloads/gear"
bfileprefix = "./inst/1KG-EUR.example2"

#### ++ Step 1 -----------------------------------------------------------------
## deepKin Principle I:
## Evaluate the minimal number of me required for detecting target relatedness
n = as.numeric(system(paste0("awk 'END{print NR}' ", bfileprefix, ".fam"), intern = T))
m = as.numeric(system(paste0("awk 'END{print NR}' ", bfileprefix, ".bim"), intern = T))
id = read.table(file = paste0(bfileprefix,".fam"), header = F)[,1]
n1 = 99
n2 = n-n1
npairs = n1*n2
degree = 1:5
theta = (1/2)^degree
minMe = minMe(theta = theta, alpha = 0.05/npairs, beta = 0.1)

#### ++ Step 2 -----------------------------------------------------------------
## Perform your data QC and calculate me
## GRM method or LB method
me = calculateMe(bfileprefix = bfileprefix, method = "GRM", plink_path = plink_path, pop_size = n)  # not suggested for biobank data
# me = calculateMe(bfileprefix = bfileprefix, method = "LB", gear_path = gear_path, pop_size = n)     # suggested for biobank data

#### ++ Step 3 -----------------------------------------------------------------
## deepKin Principle II:
deeptheta = deepTheta(me = me, alpha = 0.05/npairs, beta = 0.1)
deepDegree = log(deeptheta, base = 1/2)
## deepKin Principle III:
thrd = -logpThreshold(me = me, (1/2)^seq(0,floor(deepDegree)), beta = 0.1)
## Set target degree
targetDegree = floor(deepDegree)
targetTheta = (1/2)^targetDegree

#### ++ Step 4  -------------------------------------------------------------------
## Perform plink GRM
if(!file.exists(paste0(bfileprefix,".rel.bin"))){
  system(paste0(plink_path, " --silent --bfile ", bfileprefix, " --make-rel triangle bin4 --out ", bfileprefix))
}
## Extract plink GRM
grm.rst = extractPlinkGRM(bfileprefix, xcohort = T, pop_size = n1, pop_size2 = n2)
grm.diag = grm.rst$diag
grm.tri  = grm.rst$tri

## Perform deepKin
deepkin = deepKin_estimation(grm.diag = grm.diag, grm.tri = grm.tri, xcohort = T, me = me, pop_size1 = n1, pop_size2 = n2)

#### ++ Step 5 --------------------------------------------------------------------
## deepkin inference
deepkin.qc = deepkin[which(deepkin$`-logp` > (-logpThreshold(me = me, theta = targetTheta, beta = 0.1))), ]
deepkin.qc$tag = cut(deepkin.qc$`-logp`,
                     breaks = c(0,thrd,Inf),
                     labels = c("Unrelated", paste0("Degree ", (length(thrd)-1):0)))
index = as.numeric(rownames(deepkin.qc))
deepkin.qc.id = extract_individual_id(id = id, index = index, xcohort = T, pop_size1 = n1, pop_size2 = n2)
deepkin.qc.rst = cbind(deepkin.qc.id, deepkin.qc)
