
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

## Simple examples - to illustrate three principles
As long as we get the knowledge of the sample size and the number of effective markers, we offered three functions to calculate the three principles, including:    
I) The minimum number of effective markers is required for detecting the target degree of relatives from unrelated pairs, and it is of economic utility to minimize a budget or industrial application.    
II) The deepest relatedness ($\delta$) that the data would support to be detected from unrelated pairs, it characterizes the full potential of a given dataset.    
III) Given the target degree of relatedness, how much of the power ($\pi$) could be compromised or improved.


- Function **deepKin()** allows calculation for all three principles and returns a result list.
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


dK = deepKin(n = n, me = me, alpha = 0.05, beta = 0.1, max.degree = 5)
cat(deepKin.summary(dK))

plot(dK$power.max$degree, dK$power.max$Power.max, ylab = "Power", xlab = "Degree")
abline(h=0.9, col = "blue")
abline(v=dK$delta, col = "red")
```

In a real-world analysis, we provided both relatedness estimation and related inference.    
The relatedness estimation is based on a moment estimator and required the installation of **plink2**.    

## Example one
Example one is an example to conduct deepKin on single dataset, to calculate the n*(n-1)/2 pairs of relatedness scores and perform relatedness inference on all these results.

### Prepare
Please prepare plink bfiles - **bfileprefix.bed**, **bfileprefix.bim**, and **bfileprefix.fam**.   
Please download [GEAR](https://github.com/gc5k/GEAR) if you want to calculate me using biobank-scale data.   
```
plink_path = "/usr/bin/plink2"
gear_path = "/usr/bin/gear"
bfileprefix = "./inst/1KG-EUR.example1"
```

### Step 1
Evaluate the minimum number of me required for detecting target relatedness.    
This will give you an outline of how to design your genotyping strategy based on your target relatedness.     
- Function **me.min()**:

$$m_e \geq 2 \left[ \frac{{z_{1-\alpha} + z_{1-\beta} (1-\theta_D^t)}}{{\theta_D^t}} \right]^2$$


```
## deepKin Principle I:  
n = as.numeric(system(paste0("awk 'END{print NR}' ", bfileprefix, ".fam"), intern = T))
m = as.numeric(system(paste0("awk 'END{print NR}' ", bfileprefix, ".bim"), intern = T))
id = read.table(file = paste0(bfileprefix,".fam"), header = F)[,1]
npairs = n*(n-1)/2
me.min = me.min(theta = (1/2)^(0:4), alpha = 0.05/npairs, beta = 0.1)
```

### Step 2
If you have already got both the sample and their genotypes. Perform your data QC and calculate me using GRM method or RDM (randomization) method.    

```
## GRM method or RDM method
me = calculate.me(bfileprefix = bfileprefix, method = "GRM", plink_path = plink_path, pop_size = n)  # not suggested for biobank data
me = calculate.me(bfileprefix = bfileprefix, method = "RDM", gear_path = gear_path, pop_size = n)     # suggested for biobank data
```

### Step 3
After calculating me for the data, principle II and III will give you the basic idea on how to make relatedness inference confidentially or perform individual QC with reliable cutoffs.    

- Function **theta.min()**:    
$$\theta_D^\delta \geq \frac{{z_{1-\alpha} + z_{1-\beta}}}{{\sqrt{\frac{{m_e}}{2}} + z_{1-\beta}}}$$  

- Function **power.max()**:    
$$\pi = 1-\beta = z^{-1} \left( \frac{{\sqrt{\frac{{m_e}}{2}} \theta_D^t - z_{1-\alpha}}}{{1-\theta_D^t}} \right)$$

```
## deepKin Principle II:
theta.min = theta.min(me = me, alpha = 0.05/npairs, beta = 0.1)
degree.deep = log(theta.min, base = 1/2)

## deepKin Principle III:
power.max = power.max(me, theta = (1/2)^(0:6), alpha = 0.05/npairs)
```

### Step 4
Once you are ready for performing relatedness estimation, make sure you use plink2 **--make-rel** command with argument **triangle** and **bin4** to calculate GRM, which will save you a lot of storage, because the result is stored in binary format.    

- Function **extractPlinkGRM()**: extract the diagonal and upper triangular content of GRM for further calculation.    

- Function **deepKin_estimation()**: perform relatedness estimation, based on GRM elements.    

Input me to calculate sampling variance and p-values.    


```
## Perform plink GRM
if(!file.exists(paste0(bfileprefix,".rel.bin"))){
  system(paste0(plink_path, " --silent --bfile ", bfileprefix, " --make-rel triangle bin4 --out ", bfileprefix))
}
## Extract plink GRM
grm.rst = extract.plink.grm(bfileprefix, xcohort = F, pop_size1 = n)
grm.diag = grm.rst$diag
grm.tri  = grm.rst$tri

## Perform deepKin
deepkin = deepKin.estimation(grm.diag = grm.diag, grm.tri = grm.tri, xcohort = F, me = me)

## Output estimation and p-values (This file can be large according to the number of pairs)
write.table(deepkin, file = paste0(bfileprefix, ".deepkin"), quote = F, col.names = T, row.names = F)
```

### Step 5
Perform deepkin inference.    
```
thrd = (1/2)^seq(0.5, 10.5, 1)
thrd = thrd[which(thrd>theta.min)]
deepkin.qc = deepkin[which(deepkin$king > theta.min), ]
deepkin.qc$degree = cut(deepkin.qc$king,
                        breaks = c(theta.min,thrd,1),
                        labels = c(paste0("Deepest-",length(thrd)-1), (length(thrd)-1):0))
index = as.numeric(rownames(deepkin.qc))
deepkin.qc.id = extract.indi.id(id = id, index = index, xcohort = F)
deepkin.qc.rst = cbind(deepkin.qc.id, deepkin.qc)

write.table(deepkin.qc.rst, file = paste0(bfileprefix, ".related"), quote = F, col.names = T, row.names = F)
```
