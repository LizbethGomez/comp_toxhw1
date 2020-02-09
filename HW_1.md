HW 1 - Curve fitting to cytotoxicity data
================
Lizbeth Gomez & Nastasia Gernat
2/9/2020

``` r
library(drc) ## for dose-response logistic regression
```

    ## Loading required package: MASS

    ## 
    ## 'drc' has been loaded.

    ## Please cite R and 'drc' if used for a publication,

    ## for references type 'citation()' and 'citation('drc')'.

    ## 
    ## Attaching package: 'drc'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     gaussian, getInitial

``` r
library(tidyverse) #needs to be listed after drc otherwise you run into issues
```

    ## ── Attaching packages ───────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.4
    ## ✓ tidyr   1.0.2     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()
    ## x dplyr::select() masks MASS::select()

``` r
library(dplyr)
library(janitor) ## to make column names have the same format
```

    ## 
    ## Attaching package: 'janitor'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     chisq.test, fisher.test

``` r
library(ggplot2)
library(readxl)
library(raster) ## to compute the coefficient of variation
```

    ## Loading required package: sp

    ## 
    ## Attaching package: 'raster'

    ## The following object is masked from 'package:janitor':
    ## 
    ##     crosstab

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     extract

    ## The following objects are masked from 'package:MASS':
    ## 
    ##     area, select

``` r
cyto_death = read_xlsx("data/cytotox_data.xlsx") %>%
  clean_names() 

summary(cyto_death)
```

    ##  concentration        rep_1            rep_2            rep_3       
    ##  Min.   :0.0000   Min.   : 1.457   Min.   : 4.959   Min.   : 6.819  
    ##  1st Qu.:0.0200   1st Qu.: 6.322   1st Qu.: 9.437   1st Qu.:16.693  
    ##  Median :0.1750   Median :22.739   Median :24.584   Median :26.319  
    ##  Mean   :0.9411   Mean   :17.850   Mean   :21.174   Mean   :22.690  
    ##  3rd Qu.:0.8750   3rd Qu.:25.049   3rd Qu.:30.043   3rd Qu.:28.091  
    ##  Max.   :5.0000   Max.   :32.960   Max.   :39.766   Max.   :35.285  
    ##      rep_4            rep_5            rep_6        avg_cell_death  
    ##  Min.   : 4.182   Min.   : 2.779   Min.   : 2.729   Min.   : 6.204  
    ##  1st Qu.: 6.522   1st Qu.:12.731   1st Qu.: 8.121   1st Qu.: 8.733  
    ##  Median :26.014   Median :21.094   Median :17.799   Median :24.481  
    ##  Mean   :20.353   Mean   :20.273   Mean   :19.731   Mean   :20.345  
    ##  3rd Qu.:29.160   3rd Qu.:28.153   3rd Qu.:27.656   3rd Qu.:28.268  
    ##  Max.   :40.276   Max.   :38.021   Max.   :44.000   Max.   :34.972

``` r
cyto_death_plot <- cyto_death %>% 
  ggplot(aes(x = concentration, y = avg_cell_death)) +
  ylim(0, 45) + xlim(0,6) +
  geom_point(aes(color = concentration)) 
cyto_death_plot
```

![](HW_1_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
cyto_tidy_death = 
  pivot_longer(
    cyto_death, 
    rep_1:rep_6,
    names_to = "replicates", 
    values_to = "rep")

summary(cyto_tidy_death)
```

    ##  concentration    avg_cell_death    replicates             rep        
    ##  Min.   :0.0000   Min.   : 6.204   Length:60          Min.   : 1.457  
    ##  1st Qu.:0.0100   1st Qu.: 8.219   Class :character   1st Qu.: 9.117  
    ##  Median :0.1750   Median :24.481   Mode  :character   Median :22.648  
    ##  Mean   :0.9411   Mean   :20.345                      Mean   :20.345  
    ##  3rd Qu.:1.0000   3rd Qu.:28.884                      3rd Qu.:28.666  
    ##  Max.   :5.0000   Max.   :34.972                      Max.   :44.000

``` r
cyto_death_plot <- cyto_tidy_death %>% 
  ggplot(aes(x = concentration, y = rep)) +
  ylim(0, 45) + xlim(0,6) +
  geom_point(aes(color = replicates)) 
cyto_death_plot
```

![](HW_1_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
cv_mean <- cyto_tidy_death %>%
 group_by(concentration)%>%
  summarise(mean = mean(rep), #plate to plate variation
            cv = cv(rep)) %>%
 knitr::kable()
  cv_mean
```

| concentration |      mean |       cv |
| ------------: | --------: | -------: |
|         0.000 |  8.218574 | 72.75949 |
|         0.001 |  7.537044 | 77.25855 |
|         0.010 |  6.204400 | 82.65797 |
|         0.050 | 10.275707 | 36.46630 |
|         0.100 | 24.226595 | 24.12076 |
|         0.250 | 31.977673 | 15.89939 |
|         0.500 | 34.971719 | 21.06594 |
|         1.000 | 26.418933 | 16.30543 |
|         2.500 | 28.884271 | 21.96060 |
|         5.000 | 24.735796 | 14.80658 |

``` r
regress_cyto<- drm(rep ~concentration, data = cyto_tidy_death, fct = LL.4())


regress_cyto <- drm(rep ~concentration, data = cyto_tidy_death, fct = LL.3(names = c("Slope", "Upper Limit", "LC50" )))
summary(regress_cyto) 
```

    ## 
    ## Model fitted: Log-logistic (ED50 as parameter) with lower limit at 0 (3 parms)
    ## 
    ## Parameter estimates:
    ## 
    ##                           Estimate Std. Error t-value   p-value    
    ## Slope:(Intercept)       -3.2234353  1.1921681 -2.7038  0.009019 ** 
    ## Upper Limit:(Intercept) 29.4058373  1.3221236 22.2414 < 2.2e-16 ***
    ## LC50:(Intercept)         0.0605162  0.0073582  8.2243 2.918e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error:
    ## 
    ##  7.104517 (57 degrees of freedom)
