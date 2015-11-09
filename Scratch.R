## Scratch



## Model Brainstorming

Let's think through this process a bit. At it's most basic level, you want something like

$$b_{delta} = \alpha + \beta(MPA)$$

  where $b_{delta}$ is the change in biomass over some time period as a function of whether or not there is an MPA (1,0). We would then randomly assign MPAs of varying size all over the world and measure the effect. But, no one let's us do that

### What is the dependent variable?

* mean density of fish of species *f* in time *t*
* `r length(unique(mlpa.dat$CommonName)) * length(unique(mlpa.dat$year))`
samples (`r length(unique(mlpa.dat$CommonName))` groups)
* pros:
1. Species level allows for control by "fishing"
2. Non-spatial nature means you don't have to deal as explicitly
with spatial dynamics: you're already looking at the whole area
3. Decent sample size
4. Control over all kinds of species things
* cons:
1. How to define MPA treatment: no longer binary. % of species range currently in an MPA?
2. How to control for spatial covariates?
* mean density of fish of species *f* in time *t* at site *s*
* `r length(unique(mlpa.dat$CommonName)) * length(unique(mlpa.dat$year)) * length(unique(mlpa.dat$site))`
samples (`r length(unique(mlpa.dat$CommonName)) * length(unique(mlpa.dat$site))` groups)
* pros
1. Lots of samples
2. You now get binary treatment of MPA or not by site
3. You can do stuff like "distance from MPA"
* cons
1. You now have to deal with the spatial issues: you have correlation in time and space to deal with
2. Do you need another?
* mean density of fish in time *t* at site *s*
* `r length(unique(mlpa.dat$site)) * length(unique(mlpa.dat$year))`
samples (`r length(unique(mlpa.dat$site))` groups)
*  pros
1. simpler, gets more at the direct question (did MLPA provide a net boost in biomass)
2. Control by time?
* cons
1. Can't control by fishing anymore..


### What are the independent variables?

You've got a lot to work with, especially once you start messing with interactions...

* Fished
* MPA
* MPA start date
* Campus surveying
* All the life history mess (growth, weight, habitat,trophic group, feeding type, taxonomy)
* Temperature
* Surge
* visibility
* time
* site
* survey type
* kelp canopy
* juv/adult

Things you could create

* Distance from an MPA center
* Abundance of predators/prey
* #/area of MPAs

### Identification

This is the even harder part. A few options

__Propensity score matching of treated-untreated__ : "control" is inside and outside of MPAs. Matching by propensity score for being an MPA. Estimate the effect of a site being an MPA on the density inside (potentially broken out into fished and unfished species)

__Propensity score matching of fished-unfished__ : "control" is whether a given species is fished or unfished. Find comparable species that are fished and unfished, and then analyze density over time, weighting by propensity to be fished. The coefficients you're interested in are the effect of fished*MPA (defined somehow) within each "block"

__Difference - in -  difference between fished and unfished__ : Basically the same as above, but now using a DiD estimator, with the term of interest being fished*MPA


One potential model, treat each species as a group, then look for a difference in difference estimator i the increase in biomass following treatment of an MPA? (Lecture 13 slide "Equivalent...")
* The other idea: Can we say the causal effect of the MLPA of the biomass of fished species in the whole system? The "control" will be the biomass of unfished species (some issues with that)
* Another idea. Match MPAs with outside areas via. propensity scores, then controlling for that, regress the effect of the MLPA on biomass inside the reserve
* You could potentially include factors such as the time-varying abundance of different trophic groups?

##Data

```{r,corr plots, cache=F}

# pairs(processed.dat[,c('biomass','temp','vis','month','depth','max_total_length','pctcnpy','fish_tl','year')], upper.panel=panel.smooth2,
#       lower.panel=panel.cor,cex.labels=1.5)

```


## Prelim Regression



Let's look at number of species by site by year

```{r}

# year.site.targeted <- processed.dat %>%
#   subset(is.na(Targeted) == F & Targeted != '') %>%
#   group_by(year,site,Targeted) %>%
#   summarise(mpa.period = unique(period), total.biomass = sum(biomass, na.rm = T),
#             site.type = unique(MPA_STATUS),years.mpa = max(0,(year - (YEAR_MPA-1)) * as.numeric(YEAR_MPA>0)), region = unique(REGION)) %>%
# mutate(log.biomass = log(pmax(total.biomass, .0001)),
#        fishable = as.numeric(Targeted == 'Targeted'),
#        MPA = as.numeric(mpa.period == 'AFTER'))
#
# datatable(year.site.targeted)
#
#  varnames <- colnames(year.site.targeted)
#
# year.site.targeted$MLPA <- ordered((year.site.targeted$mpa.period), levels = c('BEFORE','AFTER'))
#
#  total.biomass.pre.post.mlpa.plot <- ggplot(year.site.targeted, aes(MLPA,total.biomass, fill = Targeted)) +
#    geom_boxplot() +
#    ggtitle('Total biomass before and after MLPA by fished or unfished')
#
#
#  total.biomass.pre.post.mlpa.plot
#
#   total.biomass.by.time.in.mpa <-
#     ggplot(year.site.targeted, aes(factor(years.mpa),total.biomass)) +
#     geom_boxplot(aes(fill = factor(Targeted))) +
#    ggtitle('Total biomass vs years since MPA by fished or unfished')
#
#   total.biomass.by.time.in.mpa
#
#
#   year.site.targeted.regression <- lm(log.biomass ~  factor(year) + fishable + fishable*MPA +
#                                            years.mpa + region, data = year.site.targeted)

```

```{r species-site-year,cache=F}

mlpa.dat$RESERVE[mlpa.dat$RESERVE == 'YES'] <- 'IN'

species.site.year <- mlpa.dat %>%
ungroup() %>%
subset(is.na(RESERVE) == F & is.na(Targeted) == F) %>%
group_by(CommonName,site,year) %>%
summarize(fished = unique(Targeted), mean.numbers = mean(count, na.rm = T), mean.length = mean(fish_tl, na.rm = T),
mean.temp = mean(temp, na.rm = T), site.size = mean(MPAAREANM2, na.rm = T),
years.mpa = max(0,(year - (YEAR_MPA-1)) * as.numeric(YEAR_MPA>0)), inside.outside.mpa = (RESERVE)[1], year.mpa = unique(YEAR_MPA)) %>%
mutate(fished.mpa = as.numeric(fished == 'Targeted') * as.numeric(years.mpa >0))

site.year <- mlpa.dat %>%
ungroup() %>%
subset(is.na(RESERVE) == F & is.na(Targeted) == F) %>%
group_by(site,year) %>%
summarize(fished = mean(Targeted == 'Targeted'), total.numbers = sum(count, na.rm = T), mean.length = mean(fish_tl, na.rm = T),
mean.temp = mean(temp, na.rm = T), site.size = mean(MPAAREANM2, na.rm = T),
years.mpa = max(0,(year - (YEAR_MPA-1)) * as.numeric(YEAR_MPA>0)), inside.outside.mpa = (RESERVE)[1], year.mpa = unique(YEAR_MPA)) %>%
mutate(fished.mpa = as.numeric(fished == 'Targeted') * as.numeric(years.mpa >0))

fished.site.year <- mlpa.dat %>%
ungroup() %>%
subset(is.na(RESERVE) == F & is.na(Targeted) == F) %>%
group_by(Targeted,site,year) %>%
summarize(total.numbers = sum(count, na.rm = T), mean.length = mean(fish_tl, na.rm = T), mean.temp = mean(temp, na.rm = T), site.size = mean(MPAAREANM2, na.rm = T),
years.mpa = max(0,(year - (YEAR_MPA-1)) * as.numeric(YEAR_MPA>0)), inside.outside.mpa = (RESERVE)[1], year.mpa = unique(YEAR_MPA)) %>%
mutate(fished.mpa = as.numeric(Targeted == 'Targeted') * as.numeric(years.mpa >0))
```

```{r species site year regress}
varnames <- colnames(species.site.year)

# datatable(species.site.year)

fmla <- as.formula(paste("mean.numbers ~ ", paste( varnames[!varnames %in% c('mean.numbers','CommonName')], collapse= "+")))

species.site.year.regression <- lm(fmla,data=species.site.year)
```

```{r site year regress}

varnames <- colnames(site.year)

# datatable(site.year)

fmla <- as.formula(paste("total.numbers ~ ", paste( varnames[!varnames %in% c('mean.numbers','CommonName')], collapse= "+")))

site.year.regression <- lm(fmla,data=site.year)

```

```{r fished site year regress}
varnames <- colnames(fished.site.year)

# datatable(fished.site.year)

fmla <- as.formula(paste("total.numbers ~ ", paste( varnames[!varnames %in% c('total.numbers','CommonName','mean.length','year.mpa','inside.outside.mpa')], collapse= "+")))

fished.site.year.regression <- lm(fmla,data=fished.site.year)

```


```{r, results='asis', echo=F}

# stargazer(year.site.targeted.regression, title = 'Prelim regression on processed biomass by site targetting over time', type = 'html')

#   stargazer(species.site.year.regression, title = 'Prelim regression on numbers of fish by species by site over time', type = 'html')
#
#     stargazer(site.year.regression, title = 'Prelim regression on numbers of fish by site over time', type = 'html')
#
#         stargazer(site.year.regression, title = 'Prelim regression on numbers of fish by fished by site over time', type = 'html')

```


$$ D_{f,s,t} = \beta_{0} + \alpha{R_{s,t}} + \beta_{1}Targeted + \sum_{s}{\beta_{2:s+1} S_{s}} +
\sum_{t}{\beta_{s+2:(s+2+t)} Y_{t}} + \beta_{3+s+t}TrophicGroup_{f} + \beta_{4+s+t}MeanTemp_{s,t} +
\sum_{g}{\beta_{5+s+t:5+s+t+g} Region_{s}} + ... + u_{f,s,t}$$
