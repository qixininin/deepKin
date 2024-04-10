
# deepKin

## Package installation
The current GitHub version of **deepKin** can be installed via:
```
library(devtools)
install_github("qixininin/deepKin")
```

## Load the library
```
library(deepKin)
```

## Simple examples - to illustrate the critical value of significant estimation and two guidelines

- Function **theta.min()**:

The critical value or the deepest relatedness that the data would support to be detected significantly related pairs at a given significant level. If we assume a significant level of $\alpha$, the critical value ($\theta_D^\delta$) is    

$$\theta_D^\delta=z_{1-\alpha}\ \sqrt{2/m_e}$$

- Function **me.min()**:

Guideline I: The minimum number of effective markers is required for detecting the target degree of relatives from unrelated pairs, and it is of economic utility to minimize a budget or industrial application.    

$$m_e \geq 2 \left[ \frac{{z_{1-\alpha} + z_{1-\beta} (1-\theta_D^t)}}{{\theta_D^t}} \right]^2$$

- Function **power.max()**:    

Guideline II: Given the target degree of relatedness, how much of the power ($\pi$) could be compromised or improved.  

$$\pi = 1-\beta = z^{-1} \left( \frac{{\sqrt{\frac{{m_e}}{2}} \theta_D^t - z_{1-\alpha}}}{{1-\theta_D^t}} \right)$$

These three functions are also in-built in the following summary function.    

- Function **deepKin()** allows calculation for the critical value and two guidelines, and returns a result list.
- Function **deepKin.summary()** provides summary of the results.
```
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
dK = deepKin(n = n, me = me, alpha = 0.05, beta = 0.1, max.degree = 5)
cat(deepKin.summary(dK))

## Visualization two guidelines
plot(dK$me.min$degree, dK$me.min$Me.min, ylab = "me", xlab = "Degree")
plot(dK$power.max$degree, dK$power.max$Power.max, ylab = "Power", xlab = "Degree")
abline(h=0.9, col = "blue")
abline(v=dK$delta, col = "red")
```

In a real-world analysis, we provided more functions on performing both relatedness estimation and relatedness inference.    

The relatedness estimation is based on a moment estimator and required the installation of **plink2**.    

## Example one
Example one is an example to conduct deepKin on a single dataset, to calculate the n*(n-1)/2 pairs of relatedness scores and perform relatedness inference on all these results.

Please prepare plink bfiles - **bfileprefix.bed**, **bfileprefix.bim**, and **bfileprefix.fam**.   
Please download [GEAR](https://github.com/gc5k/GEAR) if you want to calculate me using biobank-scale data.   

```
plink_path = "/usr/bin/plink2"
gear_path = "/usr/bin/gear"
bfileprefix = "./inst/1KG-EUR.example1"
```

### Step 1 Preparation
Evaluate the minimum number of me required for detecting target relatedness.    
This will give you an outline of how to design your genotyping strategy based on your target relatedness.     

```
## deepKin Guideline I:
## Evaluate the minimal number of me required for detecting target relatedness
n = as.numeric(system(paste0("awk 'END{print NR}' ", bfileprefix, ".fam"), intern = T))
npairs = n*(n-1)/2
me.min = me.min(theta = (1/2)^(0:4), alpha = 0.05/npairs, beta = 0.1)
```

### Step 2 The number of effective markers
If you have already got both the sample and their genotypes. Perform your data QC and calculate me using GRM method or RDM (randomization) method.    

```
## GRM method or RDM method
me = calculate.me(bfileprefix = bfileprefix, method = "GRM", plink_path = plink_path, pop_size = n)  # not suggested for biobank data
me = calculate.me(bfileprefix = bfileprefix, method = "RDM", gear_path = gear_path, pop_size = n)     # suggested for biobank data
```
Now it is able to give a summary report with $m_e$ provided.

```
## Report summary
dK = deepKin(n = n, me = me, alpha = 0.05/(n*(n-1)/2), beta = 0.1, max.degree = 5)
cat(deepKin.summary(dK))

## deepKin Critical value
theta.min = dK$theta.min
degree.deep = dK$delta
```

### Step 3 deepKin Estimation
Once you are ready for performing relatedness estimation, we employed plink2 **--make-rel** command with argument **triangle** and **bin4** to calculate GRM, which will save us a lot of storage, because the result is stored in binary format.    

- Function **deepKin.estimation()**: perform deepKin relatedness estimation based on GRM elements.    

```
## Perform deepKin
deepkin = deepKin.estimation(bfileprefix, plink_path, me, xcohort = F, pop_size1 = n)

## Output estimation and p-values (This file can be large according to the number of pairs)
write.table(deepkin, file = paste0(bfileprefix, ".deepkin"), quote = F, col.names = T, row.names = F)
```

### Step 4 deepKin Classification
Perform deepkin inference and classification.   
Based on the critical value, we are able to perform relatedness classification.    

- Function **deepKin.classification()**: performs relatedness classification, also returns p values for individual pairs belonging to t or t+1 degrees    

- Function **extract.indi.id()**: retrieves individual ID from the original plink.fam file.    

```
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
```
