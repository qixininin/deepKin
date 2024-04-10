#### Install package -----------------------------------------------------------
library(devtools)
install_github("qixininin/deepKin")

install.packages("/public3/zqx/kingless/deepKin_0.1.0.tar.gz", repos = NULL, type = "source")

#### Load package --------------------------------------------------------------
library(deepKin)


#### Simple examples -----------------------------------------------------------
## UKB british white
n = 427287
me = 56945

## UKB Chinese
n = 1435
me = 38948

## Middle East
n = 137
me = 22049

## Report summary
dK = deepKin(n = n, me = me, alpha = 0.05/(n*(n-1)/2), beta = 0.1, max.degree = 5)
cat(deepKin.summary(dK))

## Visualization two guidelines
plot(dK$me.min$degree, dK$me.min$Me.min, ylab = "me", xlab = "Degree")
plot(dK$power.max$degree, dK$power.max$Power.max, ylab = "Power", xlab = "Degree")
abline(h=0.9, col = "blue")
abline(v=dK$delta, col = "red")


#### Example One: one cohort ---------------------------------------------------
#### ++ Prepare ----------------------------------------------------------------
plink_path = "/usr/bin/plink2"
gear_path = "/usr/bin/gear"
# bfileprefix = "/public3/zqx/kingless/oxford3K/oxford3K.qc.maf005"
# bfileprefix = "/public3/zqx/kingless/oxford3K/oxford3K.qc.maf02"
# bfileprefix = "/public3/zqx/kingless/oxford3K/oxford3K.qc.maf02.noaim005"
bfileprefix = "/public3/zqx/kingless/oxford3K/oxford3K.qc.maf02.prune01"

#### ++ Step 1 -----------------------------------------------------------------
## deepKin Guideline I:
## Evaluate the minimal number of me required for detecting target relatedness
n = as.numeric(system(paste0("awk 'END{print NR}' ", bfileprefix, ".fam"), intern = T))
npairs = n*(n-1)/2
me.min = me.min(theta = (1/2)^(0:4), alpha = 0.05/npairs, beta = 0.1)

#### ++ Step 2 -----------------------------------------------------------------
## Perform your data QC and calculate me
## GRM method or RDM method
# me = calculate.me(bfileprefix = bfileprefix, method = "GRM", plink_path = plink_path, pop_size = n)  # not suggested for biobank data
me = calculate.me(bfileprefix = bfileprefix, method = "RDM", gear_path = gear_path, pop_size = n)     # suggested for biobank data

## Report summary
dK = deepKin(n = n, me = me, alpha = 0.05/(n*(n-1)/2), beta = 0.1, max.degree = 5)
cat(deepKin.summary(dK))

## deepKin Critical value
theta.min = dK$theta.min
degree.deep = dK$delta

#### ++ Step 3  -------------------------------------------------------------------
## Perform deepKin
deepkin = deepKin.estimation(bfileprefix, plink_path, me, xcohort = F, pop_size1 = n)

## Output estimation and p-values (This file can be large according to the number of pairs)
write.table(deepkin, file = paste0(bfileprefix, ".deepkin"), quote = F, col.names = T, row.names = F)

#### ++ Step 4 --------------------------------------------------------------------
## deepkin inference
deepkin.qc = deepkin[which(deepkin$theta > theta.min), ]
## individual ID
id = read.table(file = paste0(bfileprefix,".fam"), header = F)[,2]
index = as.numeric(rownames(deepkin.qc))
deepkin.qc.id = extract.indi.id(id = id, index = index, xcohort = F)
## classification
deepkin.qc.class = deepKin.classification(deepkin.qc$theta, me = me, alpha = 0.05/npairs)
## output results
deepkin.qc.rst = cbind(deepkin.qc.id, deepkin.qc.class)
write.table(deepkin.qc.rst, file = paste0(bfileprefix, ".related"), quote = F, col.names = T, row.names = F)

#### Example Two: cross cohorts ------------------------------------------------
#### ++ Prepare ----------------------------------------------------------------
plink_path = "/usr/bin/plink2"
gear_path = "/usr/bin/gear"
bfileprefix = "/public3/zqx/kingless/oxford3K/oxford3K.qc.maf02.prune01"

#### ++ Step 1 -----------------------------------------------------------------
## deepKin Principle I:
## Evaluate the minimal number of me required for detecting target relatedness
n = as.numeric(system(paste0("awk 'END{print NR}' ", bfileprefix, ".fam"), intern = T))
m = as.numeric(system(paste0("awk 'END{print NR}' ", bfileprefix, ".bim"), intern = T))
id = read.table(file = paste0(bfileprefix,".fam"), header = F)[,2]
n1 = 100
n2 = n-n1
npairs = n1*n2
me.min = me.min(theta = (1/2)^(0:4), alpha = 0.05/npairs, beta = 0.1)

#### ++ Step 2 -----------------------------------------------------------------
## Perform your data QC and calculate me
## GRM method or RDM method
# me = calculate.me(bfileprefix = bfileprefix, method = "GRM", plink_path = plink_path, pop_size = n)  # not suggested for biobank data
me = calculate.me(bfileprefix = bfileprefix, method = "RDM", gear_path = gear_path, pop_size = n)     # suggested for biobank data

#### ++ Step 3 -----------------------------------------------------------------
## deepKin Principle II:
theta.min = theta.min(me = me, alpha = 0.05/npairs)
degree.deep = log(theta.min, base = 1/2)

## deepKin Principle III:
power.max = power.max(me, theta = (1/2)^(0:6), alpha = 0.05/npairs)

#### ++ Step 4  -------------------------------------------------------------------
## Perform plink GRM
if(!file.exists(paste0(bfileprefix,".rel.bin"))){
  system(paste0(plink_path, " --silent --bfile ", bfileprefix, " --make-rel triangle bin4 --out ", bfileprefix))
}
## Extract plink GRM
grm.rst = extract.plink.grm(bfileprefix, xcohort = T, pop_size1 = n1, pop_size2 = n2)
grm.diag = grm.rst$diag
grm.tri  = grm.rst$tri

## Perform deepKin
deepkin = deepKin.estimation(grm.diag = grm.diag, grm.tri = grm.tri, xcohort = T, me = me, pop_size1 = n1, pop_size2 = n2)

#### ++ Step 5 --------------------------------------------------------------------
## deepkin inference
deepkin.qc = deepkin[which(deepkin$king > theta.min), ]
## individual ID
index = as.numeric(rownames(deepkin.qc))
deepkin.qc.id = extract.indi.id(id = id, index = index, xcohort = T, pop_size1 = n1, pop_size2 = n2)
## classification
deepkin.qc.class = deepKin.classification(deepkin.qc$king, me = me, alpha = 0.05/npairs)
## output results
deepkin.qc.rst = cbind(deepkin.qc.id, deepkin.qc.class)
write.table(deepkin.qc.rst, file = paste0(bfileprefix, ".related"), quote = F, col.names = T, row.names = F)
