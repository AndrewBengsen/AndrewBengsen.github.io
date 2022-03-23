Andrew Bengsen, NSW Department of Primary Industries  
<andrew.bengsen@dpi.nsw.gov.au>  
updated 2022-03-23  
This document is compiled from model output saved on 2022-01-20.

## Background

This file is adapted from the telemetry-augmented spatial mark-resight
model outlined in chapter 19 of Royle *et al.* (2013). Here, we use
camera trap data augmented with feral pig location data provided by GPS
tracking collars to estimate pig density around the Gwydir wetlands. The
camera survey ran from August 2019 to May 2020, but pig collars were
recycled onto different pigs in November 2019, so we have at least two
distinct survey periods: August 2019 to 11 November 2019 and 27 November
2019 to May 2020. This file deals with the first of these periods.

## Setup

First, load the dependencies, including a source file `SMR_functions.R`
that contains processing and modelling functions.

Then load the detection and spatial data. The three detection data
objects were all created by extracting exif metadata from tagged images
using the file `3-DataWrangling_Marked_gwydir.R`.

Spatial coordinates for the state space and detections are centred on
the centre of the state space to try and reduce autocorrelation.  
The main outputs are:

-   cam.region1 = state space polygon in GDA94 NSW Lambert CC
-   cam.regionS = state space polygon, centred, in km scale
-   **cam.region** = state space polygon, as.owin(), in km scale
-   camlocs = camera locations, centred with cam.regionS, in km scale
-   **locs** = camera location coordinates, centred with cam.regionS, in
    km scale
-   **Yk** = 3D array of known individual detections x camera x day (r,
    c, slice)
-   **Yu** = 2D matrix showing counts of unknown individual detections
    camera x day (r, c)

Items in **bold** are model inputs.

``` r
# Start day of first period, expressed as number of days since the start of the survey
startDay1  <- 3 # Day 3 = 2019-08-23
# End day of subset
endDay1   <- startDay1+81 
days1 <- c(startDay1 : endDay1)
dates1 <- days1 + max(as.Date(cams$Start, format="%d/%m/%Y")) -1
startDate1 <- min(dates1)
endDate1 <- max(dates1)
  
Yk <- Yk[,, days1]
Nk <- rowSums(Yk)
Yu <- Yu[, days1]

Nu <- rowSums(Yu[,2:length(Yu)], na.rm=T) 
nnDist <- round(mean(nndist(camlocs))*1000)

# Add empty rows to Yk for four collared pigs never detected
Yk_empty <- array(dim=dim(Yk))
 row.names(Yk_empty) <- c("C_42229", "C_42230", "C_42243", "C_42244")
 Yk_empty[is.na(Yk_empty)] <- 0

Yk <- abind::abind(Yk, Yk_empty, along=1)

# Evenness of detections
sums2 <- rowSums(Yu[which(rowSums(Yu, na.rm=T) > 0), ], na.rm=T)
sums2cv <- sd(sums2)/mean(sums2)

# Number of individual captures
r_yk <- sum(Yk, na.rm=T)

# How many individual captures by pig
n_indivs <- rowSums(Yk)[which(rowSums(Yk)>0)]

# How many individual captures by camera
n_indiv_cams <- rowSums(colSums(Yk))[which(rowSums(colSums(Yk))>0)]
```

There were 41 camera stations spaced about 733 m from their nearest
neighbour using a hexagonal grid. Camera stations operated for up to 82
days. There were 1065 detections of unidentifiable pigs across 23 camera
stations. The state space was defined using a 3km buffer around the
outer cameras.

``` r
# Load, thin to daily fixes (for autocorrelation) and trim to survey period
teldat_df <- read.csv("Geographic/Moree/moree_pigs_trimmed_period1_4326.csv", 
                      header=T) %>%
  mutate(date = as.Date(dtg_fix, format="%d/%m/%Y")) %>%
  group_by(collar_serial_number, date) %>%
  summarise(lat = latitude[1],
            lon = longitude[1]) %>%
  filter(date >= startDate1,
         date <= endDate1) %>%
  mutate(collar_serial_number = factor(paste0("C_", collar_serial_number)))

# Convert to sf spatial object
teldat_sf <- teldat_df %>%
  st_as_sf(coords=c("lon", "lat"), crs=4326) 

# Transform to epsg 3308, centre, and scale to km 
teldatCRS <- WGS84toCRS(teldat_df$lon, teldat_df$lat)
 teldat_df$X <- teldatCRS$X
 teldat_df$Y <- teldatCRS$Y

teldat_df <- teldat_df %>%
  mutate(X = X-centroidS$x, Y = Y-centroidS$y) %>%
  mutate(X = X/1000, Y = Y/1000)
 head(teldat_df)
```

    ## # A tibble: 6 x 6
    ## # Groups:   collar_serial_number [1]
    ##   collar_serial_number date         lat   lon      X     Y
    ##   <fct>                <date>     <dbl> <dbl>  <dbl> <dbl>
    ## 1 C_42135              2019-08-25 -29.4  149. -0.584 -1.86
    ## 2 C_42135              2019-08-26 -29.4  149. -0.549 -1.62
    ## 3 C_42135              2019-08-27 -29.4  149. -0.529 -1.83
    ## 4 C_42135              2019-08-28 -29.4  149. -0.479 -1.87
    ## 5 C_42135              2019-08-29 -29.4  149. -0.666 -1.66
    ## 6 C_42135              2019-08-30 -29.4  149. -0.478 -1.87

``` r
# Create list of location data matrices for scrPID.tel
pig_num <- unique(teldat_df$collar_serial_number)
tel_locs <- vector(mode = "list", length = length(pig_num))

for(i in 1:length(pig_num)){
  tel_locs[[i]] <- teldat_df %>%
    filter(collar_serial_number==pig_num[i]) %>%
    as.data.frame() %>%
    select(X, Y) %>%
    as.matrix()
}
# Create a vector of telemetry references for Yk for scrPID.tel
telID <- c(1, 2, 5, 6, 3, 7, 8, 4)
```

![**Camera array and pig detections at Gwydir. Point size reflects the
relative number of pig detections, increasing on a log
scale.**](Attachment1_Gwydir_SMRtel_period1_2019_files/figure-markdown_github/Plot%20pig%20detections-1.png)

There were 8 pigs collared, but only 1 identifiable individual detected
at 2 sites. This pig was recaptured 5 times.

## Prepare and run the model

All of the necessary data are now loaded, formatted and ready to be
used. We now have to specifiy the conditions for the model. This
includes:

-   setting a data augmentation parameter `M` to provide a ceiling for
    the maximum number of unidentifiable pigs
-   setting the precision of the candidate distrutions for detection
    parameters (delta\[1\] = *σ*, delta\[2\] = $lambda$0) and activity
    range centres for collared and uncollared animals (delta\[3\],
    delta\[4\])
-   setting prior distributions for *σ* and $$0
-   defining the limits of the state space (`xlim`, `ylim`)
-   setting the number of MCMC draws `ni`, burn-in draws to discard
    `nb`, and parallel chains `nc` to use
-   setting a seed value to ensure reproducibility

``` r
M <- 300     # data augmentation for unmarked portion
m <- nrow(Yk) # known number of marked animals

y <- apply(Yk, 1:2, sum) # marked detections
X <- locs

delta <- c(0.14, 0.01, 2.5, 2)
sigma.prior <- list("uniform",0, 5)
lam0.prior <- list("uniform",0, 5)

# Define limits of state space
xlim <- cam.regionS@bbox[1,]
ylim <- cam.regionS@bbox[2,]

ni <- 60000
nb <- 10000
nc <- 6

set.seed(1080)
```

The following chunk will run the model fresh if `modTypeM` == `"new"`.
This may take several days to complete depending on the number of MCMC
draws required to achieve adequate convergence and coverage.

``` r
if(modTypeM == "new"){

library(parallel)
  startDate <- Sys.Date()
  toc <- Sys.time()
  logfile <- paste0("Logfiles/GW_pig_tel_log_", startDate, ".txt")
  cl<- makePSOCKcluster(nc, outfile=logfile)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl, {library(scrbook)})
  clusterExport(cl, c("scrPID.tel"))

  ## Basic model (no telemetry), known number of marks
  inits_pid.tel <- function(){
    list(psi=runif(1), sigma=0.5, lam0=0.5, 
         S=cbind(runif(M1+m, xlim[1], xlim[2]),
                 runif(M1+m, ylim[1], ylim[2])))
    }

  mod <- clusterCall(cl, scrPID.tel, n=Yu, X=X, y=Yk, M=M1+m, locs = tel_locs, 
                     telID = telID, obsmod="pois", niters=ni, xlims=xlim, ylims=ylim, 
                     inits=inits_pid.tel(), delta=delta)
  

  stopCluster(cl)
  fullmod <- mod
  
  tic <- Sys.time()
  (runTime <- tic-toc)
  beepr::beep()
   
  for(i in 1:nc) {
    # Add density and convert to mcmc object
    D <- fullmod[[i]][,5] / area_s
    fullmod[[i]] <- cbind(fullmod[[i]], D)
    fullmod[[i]] <- as.mcmc(fullmod[[i]])
    }
  fullmod <- mcmc.list(fullmod)
  saveRDS(mod, paste0("Output/SCRpid_GW_PIG_tel", startDate))
  for(i in 1:nc) {
    mod[[i]] <- mod[[i]][-c(1:nb),] # drop burnin
    mod[[i]] <- as.mcmc(mod[[i]])
  }
  mod <- mcmc.list(mod)
} # modTypeM
```

Alternatively, we can examine the output of a model that has already
been run by specifiying `modTypeM` as `"old"` and providing the date on
which the model was started to `modDateM`.

``` r
modDateM<-"2022-01-20"
if(modTypeM == "old"){
  fullmod <- readRDS(paste0("Output/SCRpid_GW_PIG_", modDateM))
  for(i in 1:nc) {
    # Add density and convert to mcmc object
    D <- fullmod[[i]][,5] / area_s
    fullmod[[i]] <- cbind(fullmod[[i]], D)
    fullmod[[i]] <- as.mcmc(fullmod[[i]])
    }
  fullmod <- mcmc.list(fullmod)
  mod <- fullmod
  for(i in 1:nc) {
    mod[[i]] <- mod[[i]][-c(1:nb),] # drop burnin
    mod[[i]] <- as.mcmc(mod[[i]])
  }
  mod<- mcmc.list(mod)
}
```

## Results

We can now plot the MCMC traces to check whether burn in, convergence
and coverage are adequate. In this case, they are. We can also check
that the traces for psi are not bumping up against 1.0, which would
indicate that the data augmentation parameter `M` might be too small. We
also check the effective sample size to make sure that we have
sufficient independent data to be confident of point and credible
interval estimates and the acceptance rates for *σ* and $$0 to see
whether they are in the favourable window. Finally, we check the
Gelman-Rubin statistic (potential scale reduction factor) *R̂* for the
key parameters.

![**Trace plots for main parameters. Vertical red lines indicate the
burn in
threshold.**](Attachment1_Gwydir_SMRtel_period1_2019_files/figure-markdown_github/Plot%20traces-1.png)

    ## Effective Sample Size =

    ##    sigma     lam0        c      psi        N        D 
    ## 2538.059 1627.797    0.000 3852.293 2802.131 2802.131

    ## Acceptance rate =

    ## sigma  lam0     c   psi     N     D 
    ## 0.186 0.198 0.000 1.000 0.964 0.964

|       | Point est. | Upper C.I. |
|:------|-----------:|-----------:|
| sigma |          1 |        1.0 |
| lam0  |          1 |        1.1 |
| psi   |          1 |        1.1 |
| N     |          1 |        1.1 |
| D     |          1 |        1.1 |

Potential scale reduction factors

We now extract, collate and examine the parameter estimates. We can also
check that the ratio of buffer distance to estimated *σ* is \> 2, which
is a good rule of thumb for minimising the chance of animals with range
centres outside the state space being detected on the survey grid (Royle
et al 2013).

``` r
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
modmat <- as.matrix(mod, chains=T)
 modmat <- modmat[,c(1:3, 5:7)]

modes <- numeric(ncol(modmat)-1)
 for(i in c(1:5)){
   modes[i] <- getmode(modmat[,i+1])
  }
medians <- numeric(ncol(modmat)-1)
 for(i in c(1:5)){
   medians[i] <- median(modmat[,i+1])
  }
sumTab <- as.data.frame(cbind(Summ$statistics, Summ$quantiles)) 

# Add row for density
Dtab <- sumTab["N",] %>%
  mutate(Mean = mean(D),
         "2.5%" = quantile(D, probs = 0.025),
         "25%" = quantile(D, probs = 0.25),
         "50%" = quantile(D, probs = 0.25),
         "75%" = quantile(D, probs = 0.25),
         "97.5%"= quantile(D, probs = 0.975),
         SD = sd(D), 
         "Naive SE" = NA, "Time-series SE" = NA)

sumTab <- sumTab %>%
  mutate(Par = row.names(Summ$statistics), 
         ESS = ES) %>%
  mutate(MCSE = SD/sqrt(ESS),
         Acc = as.numeric(acceptance),
         Rhat = as.numeric(Rhat$psrf[,1])) %>%
  filter(Par != "c") %>%
  mutate(Mode = modes,
         Median = round(medians, 3),
         CV = SD/Mean) %>%
  select(Par, Mean,  Median, Mode, "2.5%", "97.5%", SD, CV, ESS, MCSE, Rhat, Acc)
row.names(sumTab) <- NULL

buffer2sigma <- (maskBuffer/1000) / sumTab$Median[which(sumTab$Par == "sigma")] 
```

In the table below, *m* represents the estimate for the number of marked
individuals and *N* is the estimated total population of unmarked and
marked individuals (*D* is density). *Psi* is the data augmentation
parameter for the unmarked population.

The ratio of state space buffer distance to sigma is 3.05 so we can be
confident that the state space is large enough to produce reliable
results.

Finally, we plot the detection function and the posterior distribution
of feral pig density.

``` r
dHN <- function(lam0, sigma, d){
  lam0 * exp(-d^2/(2*sigma^2))
}

lam0 <- sumTab$Mode[sumTab$Par == "lam0"]
sig  <- sumTab$Mode[sumTab$Par == "sigma"]
d<-seq(0,3,0.01)

# NB: the dHN function takes lam0 as argument to describe the intercept,
# but we can convert it to g0 by dividing by ndays * 100
hn <- data.frame(d = d, pij = dHN(lam0 = lam0, sigma = sig, d = d))
g0 <- lam0 / dim(Yu)[2] * 100
```

``` r
Mat <- as.data.frame(as.matrix(mod)) %>%
  mutate(D=N/area_s)
Mat$ni <- c(1:nrow(Mat))
modeN <- sumTab$Mode[sumTab$Par == "N"]
modeD <- sumTab$Mode[sumTab$Par == "D"]
medianN <- sumTab$Median[sumTab$Par == "N"]
medianD <- sumTab$Median[sumTab$Par == "D"]
```

![**(a) Half-normal detection function and (b) posterior distribution of
population density from pigs detections at
Gwydir**](Attachment1_Gwydir_SMRtel_period1_2019_files/figure-markdown_github/unnamed-chunk-2-1.png)

## References

Royle, J. A., Chandler, R. B., Sollmann, R. and Gardner, B. (2013).
*‘Spatial capture-recapture.’* (Academic Press: Waltham.)
