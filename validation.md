Validation of Noncompartmental Analysis Performed by NonCompart R
package
================
Sungpil Han <shan@acp.kr>
2018-05-11



-----

# Introduction

NonCompart R package (Bae 2018; Kim et al. 2018) can conduct a
noncompartmental analysis as similar as possible to the most widely used
commercial software for pharmacokinetic analysis, i.e.
[Phoenix<sup>®</sup>
WinNonlin<sup>®</sup>](https://www.certara.com/software/pkpd-modeling-and-simulation/phoenix-winnonlin/)
(Certara USA 2018). This document provides validation of
noncompartmental analysis performed by NonCompart R package version
0.4.2 as compared to the results from the commercial software,
WinNonlin<sup>®</sup> version 6.3 and 7.0.

# Results

A function, `Equal()` will return `TRUE` if there is no difference
between results from NonCompart and WinNonlin.

``` r
# install.packages("NonCompart", repos="http://r.acr.kr")
library(NonCompart)
RptCfg = read.csv("RptCfg.csv", as.is=TRUE)

Equal = function(Wres, Rres, Tol=0.001)
{
  Wres[,"ID"] = as.character(Wres[,"Subject"])
  ColName0 = colnames(Rres)
  rownames(RptCfg) = RptCfg[,"PPTESTCD"]
  colnames(Rres) = c(ColName0[1], RptCfg[ColName0[-1],"WNL"])
  Inter = intersect(colnames(Wres), colnames(Rres))
  
  IsSame = TRUE
  for (i in 1:nrow(Wres)) {
    for (j in Inter) {
      R = as.numeric(Rres[i,j])
      W = as.numeric(Wres[i,j])
      if (W != 0) {
        if(abs((R - W)/W) > Tol) {
          print(Wres[i,j])
          print(Rres[i,j])
          IsSame = FALSE
        }
      }
    }
  }
  return(IsSame)
}
```

Eight comparison tests were performed using `Theoph` and `Indometh`
default datasets. (Table 1) Detailed side-by-side comparison is in
Appendix
A.

| No. | Dataset        | Down   | Route                | Hyperlink                                                                                                                                 |
| --: | :------------- | :----- | :------------------- | :---------------------------------------------------------------------------------------------------------------------------------------- |
|   1 | Theoph (n=12)  | Linear | Extravascular        | [CSV](https://raw.githubusercontent.com/asancpt/NonCompart-tests/master/Final_Parameters_Pivoted_Theoph_Linear.csv)                       |
|   2 | Theoph (n=12)  | Log    | Extravascular        | [CSV](https://raw.githubusercontent.com/asancpt/NonCompart-tests/master/Final_Parameters_Pivoted_Theoph_Log.csv)                          |
|   3 | Indometh (n=6) | Linear | IV Bolus             | [CSV](https://raw.githubusercontent.com/asancpt/NonCompart-tests/master/Final_Parameters_Pivoted_Indometh_Linear.csv)                     |
|   4 | Indometh (n=6) | Log    | IV Bolus             | [CSV](https://raw.githubusercontent.com/asancpt/NonCompart-tests/master/Final_Parameters_Pivoted_Indometh_Log.csv)                        |
|   5 | Indometh (n=6) | Linear | IV Infusion (0.25hr) | [CSV](https://raw.githubusercontent.com/asancpt/NonCompart-tests/master/Final_Parameters_Pivoted_Indometh_Linear_Infusion.csv)            |
|   6 | Indometh (n=6) | Log    | IV Infusion (0.25hr) | [CSV](https://raw.githubusercontent.com/asancpt/NonCompart-tests/master/Final_Parameters_Pivoted_Indometh_Log_Infusion.csv)               |
|   7 | Indometh (n=6) | Linear | Extravascular        | [CSV](https://raw.githubusercontent.com/asancpt/NonCompart-tests/master/Final_Parameters_Pivoted_Indometh_Linear_Wrong_Extravascular.csv) |
|   8 | Indometh (n=6) | Log    | Extravascular        | [CSV](https://raw.githubusercontent.com/asancpt/NonCompart-tests/master/Final_Parameters_Pivoted_Indometh_Log_Wrong_Extravascular.csv)    |

Description of settings for the noncompartmental analysis performed in
WinNonlin and links to the raw data

``` r
Theoph[,"Subject"] = as.numeric(as.character(Theoph[,"Subject"]))
Indometh[,"Subject"] = as.numeric(as.character(Indometh[,"Subject"]))

Wres1 = read.csv("Final_Parameters_Pivoted_Theoph_Linear.csv")
Rres1 = tblNCA(Theoph, "Subject", "Time", "conc", dose=320, concUnit="mg/L")
Equal(Wres1, Rres1)
```

    ## [1] TRUE

``` r
Wres2 = read.csv("Final_Parameters_Pivoted_Theoph_Log.csv")
Rres2 = tblNCA(Theoph, "Subject", "Time", "conc", dose=320, down="Log", 
               concUnit="mg/L")
Equal(Wres2, Rres2)
```

    ## [1] TRUE

``` r
Wres3 = read.csv("Final_Parameters_Pivoted_Indometh_Linear.csv")
Rres3 = tblNCA(Indometh, "Subject", "time", "conc", dose=25, adm="Bolus", 
               concUnit="mg/L", R2ADJ=0.8)
Equal(Wres3, Rres3)
```

    ## [1] TRUE

``` r
Wres4 = read.csv("Final_Parameters_Pivoted_Indometh_Log.csv")
Rres4 = tblNCA(Indometh, "Subject", "time", "conc", dose=25, adm="Bolus", 
               down="Log", concUnit="mg/L", R2ADJ=0.8)
Equal(Wres4, Rres4)
```

    ## [1] TRUE

``` r
Wres5 = read.csv("Final_Parameters_Pivoted_Indometh_Linear_Infusion.csv")
Rres5 = tblNCA(Indometh, "Subject", "time", "conc", dose=25, adm="Infusion", 
               dur=0.25, concUnit="mg/L", R2ADJ=0.8)
Equal(Wres5, Rres5)
```

    ## [1] TRUE

``` r
Wres6 = read.csv("Final_Parameters_Pivoted_Indometh_Log_Infusion.csv")
Rres6 = tblNCA(Indometh, "Subject", "time", "conc", dose=25, adm="Infusion", 
               dur=0.25, down="Log", concUnit="mg/L", R2ADJ=0.8)
Equal(Wres6, Rres6)
```

    ## [1] TRUE

``` r
Wres7 = read.csv("Final_Parameters_Pivoted_Indometh_Linear_Wrong_Extravascular.csv")
Rres7 = tblNCA(Indometh, "Subject", "time", "conc", dose=25, concUnit="mg/L", 
               R2ADJ=0.8)
Equal(Wres7, Rres7)
```

    ## [1] TRUE

``` r
Wres8 = read.csv("Final_Parameters_Pivoted_Indometh_Log_Wrong_Extravascular.csv")
Rres8 = tblNCA(Indometh, "Subject", "time", "conc", dose=25, down="Log", 
               concUnit="mg/L", R2ADJ=0.8)
Equal(Wres8, Rres8)
```

    ## [1] TRUE

\pagebreak

# Conclusion

*There is no discrepancy* between results from NonCompart and WinNonlin.
We also performed multiple analyses with the real clinical trial
datasets and have found no differences (data not shown: confidential).
Noncompartmental analysis performed by the open-source R package,
NonCompart can be **qualified and validated** enough to acquire the
identical results of the commercial software, WinNonlin.

*Please report issues regarding validation of the R package to
<https://github.com/asancpt/NonCompart-tests/issues>.*

-----

**Affiliation**:  
Sungpil Han  
M.D/Ph.D, Resident  
Department of Clinical Pharmacology and Therapeutics,  
Asan Medical Center, University of Ulsan,  
Seoul 05505, Republic of Korea  
E-mail: <shan@acp.kr>  
URL: www.github.com/shanmdphd

\pagebreak

# (APPENDIX) Appendix

# Side-by-side comparison of results

``` r
library(dplyr)
library(tidyr)

table_wres_rres <- function(wres, rres, Caption){
  wres %>% 
    gather(WNL, WinNonlin, -Subject) %>% 
    left_join(RptCfg %>% select(PPTESTCD, WNL), by = "WNL") %>% 
    left_join(rres %>% as.data.frame() %>% gather(PPTESTCD, NonCompart, -Subject),
              by = c("Subject", "PPTESTCD")) %>% 
    select(Subject, PPTESTCD, WNL, NonCompart, WinNonlin) %>% 
    mutate(NonCompart = as.numeric(NonCompart),
           WinNonlin = as.numeric(WinNonlin)) %>% 
    mutate(Difference = ifelse((NonCompart - WinNonlin) / WinNonlin < 0.001, 
                               yes = 0, no = 'Larger than tol.')) %>% 
    mutate(Difference = ifelse(is.na(Difference), 0, Difference)) %>% 
    filter(!is.na(WinNonlin) & !is.na(NonCompart)) %>% 
    knitr::kable(longtable = TRUE, booktabs = TRUE, # format = "latex", 
                 caption = Caption) %>% 
    add_header_above(c(" ", "Pharmacokinetic Parameters" = 2, "Values" = 2, " ")) %>% 
    kable_styling(latex_options = c("repeat_header"))
}
```

## Test 1: Theoph (n=12), Linear, Extravascular

``` r
table_wres_rres(Wres1, Rres1, 
                Caption = 'Theoph (n=12), Linear, Extravascular')
```

| Subject | PPTESTCD | WNL                   |   NonCompart |    WinNonlin | Difference |
| ------: | :------- | :-------------------- | -----------: | -----------: | ---------: |
|       1 | R2       | Rsq                   |    0.9999997 |    0.9999997 |          0 |
|       2 | R2       | Rsq                   |    0.9971954 |    0.9971954 |          0 |
|       3 | R2       | Rsq                   |    0.9993250 |    0.9993250 |          0 |
|       4 | R2       | Rsq                   |    0.9989241 |    0.9989241 |          0 |
|       5 | R2       | Rsq                   |    0.9986472 |    0.9986472 |          0 |
|       6 | R2       | Rsq                   |    0.9982413 |    0.9982413 |          0 |
|       7 | R2       | Rsq                   |    0.9986702 |    0.9986702 |          0 |
|       8 | R2       | Rsq                   |    0.9910124 |    0.9910124 |          0 |
|       9 | R2       | Rsq                   |    0.9994437 |    0.9994437 |          0 |
|      10 | R2       | Rsq                   |    0.9995087 |    0.9995087 |          0 |
|      11 | R2       | Rsq                   |    0.9999983 |    0.9999983 |          0 |
|      12 | R2       | Rsq                   |    0.9993968 |    0.9993968 |          0 |
|       1 | R2ADJ    | Rsq\_adjusted         |    0.9999995 |    0.9999995 |          0 |
|       2 | R2ADJ    | Rsq\_adjusted         |    0.9957931 |    0.9957931 |          0 |
|       3 | R2ADJ    | Rsq\_adjusted         |    0.9986499 |    0.9986499 |          0 |
|       4 | R2ADJ    | Rsq\_adjusted         |    0.9978483 |    0.9978483 |          0 |
|       5 | R2ADJ    | Rsq\_adjusted         |    0.9979708 |    0.9979708 |          0 |
|       6 | R2ADJ    | Rsq\_adjusted         |    0.9978896 |    0.9978896 |          0 |
|       7 | R2ADJ    | Rsq\_adjusted         |    0.9980053 |    0.9980053 |          0 |
|       8 | R2ADJ    | Rsq\_adjusted         |    0.9887655 |    0.9887655 |          0 |
|       9 | R2ADJ    | Rsq\_adjusted         |    0.9988873 |    0.9988873 |          0 |
|      10 | R2ADJ    | Rsq\_adjusted         |    0.9990174 |    0.9990174 |          0 |
|      11 | R2ADJ    | Rsq\_adjusted         |    0.9999965 |    0.9999965 |          0 |
|      12 | R2ADJ    | Rsq\_adjusted         |    0.9987936 |    0.9987936 |          0 |
|       1 | CORRXY   | Corr\_XY              |  \-0.9999999 |  \-0.9999999 |          0 |
|       2 | CORRXY   | Corr\_XY              |  \-0.9985967 |  \-0.9985967 |          0 |
|       3 | CORRXY   | Corr\_XY              |  \-0.9996624 |  \-0.9996624 |          0 |
|       4 | CORRXY   | Corr\_XY              |  \-0.9994619 |  \-0.9994619 |          0 |
|       5 | CORRXY   | Corr\_XY              |  \-0.9993234 |  \-0.9993234 |          0 |
|       6 | CORRXY   | Corr\_XY              |  \-0.9991203 |  \-0.9991203 |          0 |
|       7 | CORRXY   | Corr\_XY              |  \-0.9993349 |  \-0.9993349 |          0 |
|       8 | CORRXY   | Corr\_XY              |  \-0.9954961 |  \-0.9954961 |          0 |
|       9 | CORRXY   | Corr\_XY              |  \-0.9997218 |  \-0.9997218 |          0 |
|      10 | CORRXY   | Corr\_XY              |  \-0.9997543 |  \-0.9997543 |          0 |
|      11 | CORRXY   | Corr\_XY              |  \-0.9999991 |  \-0.9999991 |          0 |
|      12 | CORRXY   | Corr\_XY              |  \-0.9996984 |  \-0.9996984 |          0 |
|       1 | LAMZNPT  | No\_points\_lambda\_z |    3.0000000 |    3.0000000 |          0 |
|       2 | LAMZNPT  | No\_points\_lambda\_z |    4.0000000 |    4.0000000 |          0 |
|       3 | LAMZNPT  | No\_points\_lambda\_z |    3.0000000 |    3.0000000 |          0 |
|       4 | LAMZNPT  | No\_points\_lambda\_z |    3.0000000 |    3.0000000 |          0 |
|       5 | LAMZNPT  | No\_points\_lambda\_z |    4.0000000 |    4.0000000 |          0 |
|       6 | LAMZNPT  | No\_points\_lambda\_z |    7.0000000 |    7.0000000 |          0 |
|       7 | LAMZNPT  | No\_points\_lambda\_z |    4.0000000 |    4.0000000 |          0 |
|       8 | LAMZNPT  | No\_points\_lambda\_z |    6.0000000 |    6.0000000 |          0 |
|       9 | LAMZNPT  | No\_points\_lambda\_z |    3.0000000 |    3.0000000 |          0 |
|      10 | LAMZNPT  | No\_points\_lambda\_z |    3.0000000 |    3.0000000 |          0 |
|      11 | LAMZNPT  | No\_points\_lambda\_z |    3.0000000 |    3.0000000 |          0 |
|      12 | LAMZNPT  | No\_points\_lambda\_z |    3.0000000 |    3.0000000 |          0 |
|       1 | LAMZ     | Lambda\_z             |    0.0484570 |    0.0484570 |          0 |
|       2 | LAMZ     | Lambda\_z             |    0.1040864 |    0.1040864 |          0 |
|       3 | LAMZ     | Lambda\_z             |    0.1024443 |    0.1024443 |          0 |
|       4 | LAMZ     | Lambda\_z             |    0.0992870 |    0.0992870 |          0 |
|       5 | LAMZ     | Lambda\_z             |    0.0866189 |    0.0866189 |          0 |
|       6 | LAMZ     | Lambda\_z             |    0.0877957 |    0.0877957 |          0 |
|       7 | LAMZ     | Lambda\_z             |    0.0883365 |    0.0883365 |          0 |
|       8 | LAMZ     | Lambda\_z             |    0.0814505 |    0.0814505 |          0 |
|       9 | LAMZ     | Lambda\_z             |    0.0824586 |    0.0824586 |          0 |
|      10 | LAMZ     | Lambda\_z             |    0.0749598 |    0.0749598 |          0 |
|      11 | LAMZ     | Lambda\_z             |    0.0954586 |    0.0954586 |          0 |
|      12 | LAMZ     | Lambda\_z             |    0.1102595 |    0.1102595 |          0 |
|       1 | LAMZLL   | Lambda\_z\_lower      |    9.0500000 |    9.0500000 |          0 |
|       2 | LAMZLL   | Lambda\_z\_lower      |    7.0300000 |    7.0300000 |          0 |
|       3 | LAMZLL   | Lambda\_z\_lower      |    9.0000000 |    9.0000000 |          0 |
|       4 | LAMZLL   | Lambda\_z\_lower      |    9.0200000 |    9.0200000 |          0 |
|       5 | LAMZLL   | Lambda\_z\_lower      |    7.0200000 |    7.0200000 |          0 |
|       6 | LAMZLL   | Lambda\_z\_lower      |    2.0300000 |    2.0300000 |          0 |
|       7 | LAMZLL   | Lambda\_z\_lower      |    6.9800000 |    6.9800000 |          0 |
|       8 | LAMZLL   | Lambda\_z\_lower      |    3.5300000 |    3.5300000 |          0 |
|       9 | LAMZLL   | Lambda\_z\_lower      |    8.8000000 |    8.8000000 |          0 |
|      10 | LAMZLL   | Lambda\_z\_lower      |    9.3800000 |    9.3800000 |          0 |
|      11 | LAMZLL   | Lambda\_z\_lower      |    9.0300000 |    9.0300000 |          0 |
|      12 | LAMZLL   | Lambda\_z\_lower      |    9.0300000 |    9.0300000 |          0 |
|       1 | LAMZUL   | Lambda\_z\_upper      |   24.3700000 |   24.3700000 |          0 |
|       2 | LAMZUL   | Lambda\_z\_upper      |   24.3000000 |   24.3000000 |          0 |
|       3 | LAMZUL   | Lambda\_z\_upper      |   24.1700000 |   24.1700000 |          0 |
|       4 | LAMZUL   | Lambda\_z\_upper      |   24.6500000 |   24.6500000 |          0 |
|       5 | LAMZUL   | Lambda\_z\_upper      |   24.3500000 |   24.3500000 |          0 |
|       6 | LAMZUL   | Lambda\_z\_upper      |   23.8500000 |   23.8500000 |          0 |
|       7 | LAMZUL   | Lambda\_z\_upper      |   24.2200000 |   24.2200000 |          0 |
|       8 | LAMZUL   | Lambda\_z\_upper      |   24.1200000 |   24.1200000 |          0 |
|       9 | LAMZUL   | Lambda\_z\_upper      |   24.4300000 |   24.4300000 |          0 |
|      10 | LAMZUL   | Lambda\_z\_upper      |   23.7000000 |   23.7000000 |          0 |
|      11 | LAMZUL   | Lambda\_z\_upper      |   24.0800000 |   24.0800000 |          0 |
|      12 | LAMZUL   | Lambda\_z\_upper      |   24.1500000 |   24.1500000 |          0 |
|       1 | LAMZHL   | HL\_Lambda\_z         |   14.3043776 |   14.3043776 |          0 |
|       2 | LAMZHL   | HL\_Lambda\_z         |    6.6593416 |    6.6593416 |          0 |
|       3 | LAMZHL   | HL\_Lambda\_z         |    6.7660874 |    6.7660874 |          0 |
|       4 | LAMZHL   | HL\_Lambda\_z         |    6.9812467 |    6.9812467 |          0 |
|       5 | LAMZHL   | HL\_Lambda\_z         |    8.0022640 |    8.0022640 |          0 |
|       6 | LAMZHL   | HL\_Lambda\_z         |    7.8949979 |    7.8949979 |          0 |
|       7 | LAMZHL   | HL\_Lambda\_z         |    7.8466683 |    7.8466683 |          0 |
|       8 | LAMZHL   | HL\_Lambda\_z         |    8.5100379 |    8.5100379 |          0 |
|       9 | LAMZHL   | HL\_Lambda\_z         |    8.4059988 |    8.4059988 |          0 |
|      10 | LAMZHL   | HL\_Lambda\_z         |    9.2469158 |    9.2469158 |          0 |
|      11 | LAMZHL   | HL\_Lambda\_z         |    7.2612365 |    7.2612365 |          0 |
|      12 | LAMZHL   | HL\_Lambda\_z         |    6.2865082 |    6.2865082 |          0 |
|       1 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       2 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       3 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       4 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       5 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       6 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       7 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       8 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       9 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|      10 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|      11 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|      12 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       1 | TMAX     | Tmax                  |    1.1200000 |    1.1200000 |          0 |
|       2 | TMAX     | Tmax                  |    1.9200000 |    1.9200000 |          0 |
|       3 | TMAX     | Tmax                  |    1.0200000 |    1.0200000 |          0 |
|       4 | TMAX     | Tmax                  |    1.0700000 |    1.0700000 |          0 |
|       5 | TMAX     | Tmax                  |    1.0000000 |    1.0000000 |          0 |
|       6 | TMAX     | Tmax                  |    1.1500000 |    1.1500000 |          0 |
|       7 | TMAX     | Tmax                  |    3.4800000 |    3.4800000 |          0 |
|       8 | TMAX     | Tmax                  |    2.0200000 |    2.0200000 |          0 |
|       9 | TMAX     | Tmax                  |    0.6300000 |    0.6300000 |          0 |
|      10 | TMAX     | Tmax                  |    3.5500000 |    3.5500000 |          0 |
|      11 | TMAX     | Tmax                  |    0.9800000 |    0.9800000 |          0 |
|      12 | TMAX     | Tmax                  |    3.5200000 |    3.5200000 |          0 |
|       1 | CMAX     | Cmax                  |   10.5000000 |   10.5000000 |          0 |
|       2 | CMAX     | Cmax                  |    8.3300000 |    8.3300000 |          0 |
|       3 | CMAX     | Cmax                  |    8.2000000 |    8.2000000 |          0 |
|       4 | CMAX     | Cmax                  |    8.6000000 |    8.6000000 |          0 |
|       5 | CMAX     | Cmax                  |   11.4000000 |   11.4000000 |          0 |
|       6 | CMAX     | Cmax                  |    6.4400000 |    6.4400000 |          0 |
|       7 | CMAX     | Cmax                  |    7.0900000 |    7.0900000 |          0 |
|       8 | CMAX     | Cmax                  |    7.5600000 |    7.5600000 |          0 |
|       9 | CMAX     | Cmax                  |    9.0300000 |    9.0300000 |          0 |
|      10 | CMAX     | Cmax                  |   10.2100000 |   10.2100000 |          0 |
|      11 | CMAX     | Cmax                  |    8.0000000 |    8.0000000 |          0 |
|      12 | CMAX     | Cmax                  |    9.7500000 |    9.7500000 |          0 |
|       1 | CMAXD    | Cmax\_D               |    0.0328125 |    0.0328125 |          0 |
|       2 | CMAXD    | Cmax\_D               |    0.0260312 |    0.0260312 |          0 |
|       3 | CMAXD    | Cmax\_D               |    0.0256250 |    0.0256250 |          0 |
|       4 | CMAXD    | Cmax\_D               |    0.0268750 |    0.0268750 |          0 |
|       5 | CMAXD    | Cmax\_D               |    0.0356250 |    0.0356250 |          0 |
|       6 | CMAXD    | Cmax\_D               |    0.0201250 |    0.0201250 |          0 |
|       7 | CMAXD    | Cmax\_D               |    0.0221562 |    0.0221562 |          0 |
|       8 | CMAXD    | Cmax\_D               |    0.0236250 |    0.0236250 |          0 |
|       9 | CMAXD    | Cmax\_D               |    0.0282188 |    0.0282188 |          0 |
|      10 | CMAXD    | Cmax\_D               |    0.0319063 |    0.0319062 |          0 |
|      11 | CMAXD    | Cmax\_D               |    0.0250000 |    0.0250000 |          0 |
|      12 | CMAXD    | Cmax\_D               |    0.0304688 |    0.0304688 |          0 |
|       1 | TLST     | Tlast                 |   24.3700000 |   24.3700000 |          0 |
|       2 | TLST     | Tlast                 |   24.3000000 |   24.3000000 |          0 |
|       3 | TLST     | Tlast                 |   24.1700000 |   24.1700000 |          0 |
|       4 | TLST     | Tlast                 |   24.6500000 |   24.6500000 |          0 |
|       5 | TLST     | Tlast                 |   24.3500000 |   24.3500000 |          0 |
|       6 | TLST     | Tlast                 |   23.8500000 |   23.8500000 |          0 |
|       7 | TLST     | Tlast                 |   24.2200000 |   24.2200000 |          0 |
|       8 | TLST     | Tlast                 |   24.1200000 |   24.1200000 |          0 |
|       9 | TLST     | Tlast                 |   24.4300000 |   24.4300000 |          0 |
|      10 | TLST     | Tlast                 |   23.7000000 |   23.7000000 |          0 |
|      11 | TLST     | Tlast                 |   24.0800000 |   24.0800000 |          0 |
|      12 | TLST     | Tlast                 |   24.1500000 |   24.1500000 |          0 |
|       1 | CLST     | Clast                 |    3.2800000 |    3.2800000 |          0 |
|       2 | CLST     | Clast                 |    0.9000000 |    0.9000000 |          0 |
|       3 | CLST     | Clast                 |    1.0500000 |    1.0500000 |          0 |
|       4 | CLST     | Clast                 |    1.1500000 |    1.1500000 |          0 |
|       5 | CLST     | Clast                 |    1.5700000 |    1.5700000 |          0 |
|       6 | CLST     | Clast                 |    0.9200000 |    0.9200000 |          0 |
|       7 | CLST     | Clast                 |    1.1500000 |    1.1500000 |          0 |
|       8 | CLST     | Clast                 |    1.2500000 |    1.2500000 |          0 |
|       9 | CLST     | Clast                 |    1.1200000 |    1.1200000 |          0 |
|      10 | CLST     | Clast                 |    2.4200000 |    2.4200000 |          0 |
|      11 | CLST     | Clast                 |    0.8600000 |    0.8600000 |          0 |
|      12 | CLST     | Clast                 |    1.1700000 |    1.1700000 |          0 |
|       1 | AUCLST   | AUClast               |  148.9230500 |  148.9230500 |          0 |
|       2 | AUCLST   | AUClast               |   91.5268000 |   91.5268000 |          0 |
|       3 | AUCLST   | AUClast               |   99.2865000 |   99.2865000 |          0 |
|       4 | AUCLST   | AUClast               |  106.7963000 |  106.7963000 |          0 |
|       5 | AUCLST   | AUClast               |  121.2944000 |  121.2944000 |          0 |
|       6 | AUCLST   | AUClast               |   73.7755500 |   73.7755500 |          0 |
|       7 | AUCLST   | AUClast               |   90.7534000 |   90.7534000 |          0 |
|       8 | AUCLST   | AUClast               |   88.5599500 |   88.5599500 |          0 |
|       9 | AUCLST   | AUClast               |   86.3261500 |   86.3261500 |          0 |
|      10 | AUCLST   | AUClast               |  138.3681000 |  138.3681000 |          0 |
|      11 | AUCLST   | AUClast               |   80.0936000 |   80.0936000 |          0 |
|      12 | AUCLST   | AUClast               |  119.9775000 |  119.9775000 |          0 |
|       1 | AUCALL   | AUCall                |  148.9230500 |  148.9230500 |          0 |
|       2 | AUCALL   | AUCall                |   91.5268000 |   91.5268000 |          0 |
|       3 | AUCALL   | AUCall                |   99.2865000 |   99.2865000 |          0 |
|       4 | AUCALL   | AUCall                |  106.7963000 |  106.7963000 |          0 |
|       5 | AUCALL   | AUCall                |  121.2944000 |  121.2944000 |          0 |
|       6 | AUCALL   | AUCall                |   73.7755500 |   73.7755500 |          0 |
|       7 | AUCALL   | AUCall                |   90.7534000 |   90.7534000 |          0 |
|       8 | AUCALL   | AUCall                |   88.5599500 |   88.5599500 |          0 |
|       9 | AUCALL   | AUCall                |   86.3261500 |   86.3261500 |          0 |
|      10 | AUCALL   | AUCall                |  138.3681000 |  138.3681000 |          0 |
|      11 | AUCALL   | AUCall                |   80.0936000 |   80.0936000 |          0 |
|      12 | AUCALL   | AUCall                |  119.9775000 |  119.9775000 |          0 |
|       1 | AUCIFO   | AUCINF\_obs           |  216.6119330 |  216.6119330 |          0 |
|       2 | AUCIFO   | AUCINF\_obs           |  100.1734591 |  100.1734591 |          0 |
|       3 | AUCIFO   | AUCINF\_obs           |  109.5359707 |  109.5359707 |          0 |
|       4 | AUCIFO   | AUCINF\_obs           |  118.3788814 |  118.3788814 |          0 |
|       5 | AUCIFO   | AUCINF\_obs           |  139.4197778 |  139.4197778 |          0 |
|       6 | AUCIFO   | AUCINF\_obs           |   84.2544183 |   84.2544183 |          0 |
|       7 | AUCIFO   | AUCINF\_obs           |  103.7718018 |  103.7718018 |          0 |
|       8 | AUCIFO   | AUCINF\_obs           |  103.9066868 |  103.9066868 |          0 |
|       9 | AUCIFO   | AUCINF\_obs           |   99.9087179 |   99.9087179 |          0 |
|      10 | AUCIFO   | AUCINF\_obs           |  170.6520606 |  170.6520606 |          0 |
|      11 | AUCIFO   | AUCINF\_obs           |   89.1027449 |   89.1027449 |          0 |
|      12 | AUCIFO   | AUCINF\_obs           |  130.5888316 |  130.5888316 |          0 |
|       1 | AUCIFOD  | AUCINF\_D\_obs        |    0.6769123 |    0.6769123 |          0 |
|       2 | AUCIFOD  | AUCINF\_D\_obs        |    0.3130421 |    0.3130421 |          0 |
|       3 | AUCIFOD  | AUCINF\_D\_obs        |    0.3422999 |    0.3422999 |          0 |
|       4 | AUCIFOD  | AUCINF\_D\_obs        |    0.3699340 |    0.3699340 |          0 |
|       5 | AUCIFOD  | AUCINF\_D\_obs        |    0.4356868 |    0.4356868 |          0 |
|       6 | AUCIFOD  | AUCINF\_D\_obs        |    0.2632951 |    0.2632951 |          0 |
|       7 | AUCIFOD  | AUCINF\_D\_obs        |    0.3242869 |    0.3242869 |          0 |
|       8 | AUCIFOD  | AUCINF\_D\_obs        |    0.3247084 |    0.3247084 |          0 |
|       9 | AUCIFOD  | AUCINF\_D\_obs        |    0.3122147 |    0.3122147 |          0 |
|      10 | AUCIFOD  | AUCINF\_D\_obs        |    0.5332877 |    0.5332877 |          0 |
|      11 | AUCIFOD  | AUCINF\_D\_obs        |    0.2784461 |    0.2784461 |          0 |
|      12 | AUCIFOD  | AUCINF\_D\_obs        |    0.4080901 |    0.4080901 |          0 |
|       1 | AUCPEO   | AUC\_.Extrap\_obs     |   31.2489169 |   31.2489169 |          0 |
|       2 | AUCPEO   | AUC\_.Extrap\_obs     |    8.6316867 |    8.6316867 |          0 |
|       3 | AUCPEO   | AUC\_.Extrap\_obs     |    9.3571734 |    9.3571734 |          0 |
|       4 | AUCPEO   | AUC\_.Extrap\_obs     |    9.7843309 |    9.7843309 |          0 |
|       5 | AUCPEO   | AUC\_.Extrap\_obs     |   13.0005786 |   13.0005786 |          0 |
|       6 | AUCPEO   | AUC\_.Extrap\_obs     |   12.4371737 |   12.4371737 |          0 |
|       7 | AUCPEO   | AUC\_.Extrap\_obs     |   12.5452209 |   12.5452209 |          0 |
|       8 | AUCPEO   | AUC\_.Extrap\_obs     |   14.7697297 |   14.7697297 |          0 |
|       9 | AUCPEO   | AUC\_.Extrap\_obs     |   13.5949777 |   13.5949777 |          0 |
|      10 | AUCPEO   | AUC\_.Extrap\_obs     |   18.9180022 |   18.9180022 |          0 |
|      11 | AUCPEO   | AUC\_.Extrap\_obs     |   10.1109623 |   10.1109623 |          0 |
|      12 | AUCPEO   | AUC\_.Extrap\_obs     |    8.1257573 |    8.1257573 |          0 |
|       1 | VZFO     | Vz\_F\_obs            |   30.4867482 |   30.4867482 |          0 |
|       2 | VZFO     | Vz\_F\_obs            |   30.6904416 |   30.6904416 |          0 |
|       3 | VZFO     | Vz\_F\_obs            |   28.5170999 |   28.5170999 |          0 |
|       4 | VZFO     | Vz\_F\_obs            |   27.2259641 |   27.2259641 |          0 |
|       5 | VZFO     | Vz\_F\_obs            |   26.4979947 |   26.4979946 |          0 |
|       6 | VZFO     | Vz\_F\_obs            |   43.2597345 |   43.2597345 |          0 |
|       7 | VZFO     | Vz\_F\_obs            |   34.9084408 |   34.9084408 |          0 |
|       8 | VZFO     | Vz\_F\_obs            |   37.8105081 |   37.8105081 |          0 |
|       9 | VZFO     | Vz\_F\_obs            |   38.8427934 |   38.8427934 |          0 |
|      10 | VZFO     | Vz\_F\_obs            |   25.0155401 |   25.0155401 |          0 |
|      11 | VZFO     | Vz\_F\_obs            |   37.6221852 |   37.6221852 |          0 |
|      12 | VZFO     | Vz\_F\_obs            |   22.2242936 |   22.2242936 |          0 |
|       1 | CLFO     | Cl\_F\_obs            |    1.4772963 |    1.4772963 |          0 |
|       2 | CLFO     | Cl\_F\_obs            |    3.1944589 |    3.1944589 |          0 |
|       3 | CLFO     | Cl\_F\_obs            |    2.9214147 |    2.9214147 |          0 |
|       4 | CLFO     | Cl\_F\_obs            |    2.7031849 |    2.7031849 |          0 |
|       5 | CLFO     | Cl\_F\_obs            |    2.2952267 |    2.2952267 |          0 |
|       6 | CLFO     | Cl\_F\_obs            |    3.7980204 |    3.7980204 |          0 |
|       7 | CLFO     | Cl\_F\_obs            |    3.0836893 |    3.0836894 |          0 |
|       8 | CLFO     | Cl\_F\_obs            |    3.0796863 |    3.0796863 |          0 |
|       9 | CLFO     | Cl\_F\_obs            |    3.2029237 |    3.2029237 |          0 |
|      10 | CLFO     | Cl\_F\_obs            |    1.8751605 |    1.8751605 |          0 |
|      11 | CLFO     | Cl\_F\_obs            |    3.5913596 |    3.5913596 |          0 |
|      12 | CLFO     | Cl\_F\_obs            |    2.4504393 |    2.4504393 |          0 |
|       1 | AUCIFP   | AUCINF\_pred          |  216.6149558 |  216.6149558 |          0 |
|       2 | AUCIFP   | AUCINF\_pred          |  100.0643176 |  100.0643176 |          0 |
|       3 | AUCIFP   | AUCINF\_pred          |  109.5857218 |  109.5857218 |          0 |
|       4 | AUCIFP   | AUCINF\_pred          |  118.4435586 |  118.4435586 |          0 |
|       5 | AUCIFP   | AUCINF\_pred          |  139.2546304 |  139.2546304 |          0 |
|       6 | AUCIFP   | AUCINF\_pred          |   84.4966986 |   84.4966986 |          0 |
|       7 | AUCIFP   | AUCINF\_pred          |  103.8931470 |  103.8931470 |          0 |
|       8 | AUCIFP   | AUCINF\_pred          |  103.6430515 |  103.6430515 |          0 |
|       9 | AUCIFP   | AUCINF\_pred          |   99.8660677 |   99.8660677 |          0 |
|      10 | AUCIFP   | AUCINF\_pred          |  170.5679125 |  170.5679125 |          0 |
|      11 | AUCIFP   | AUCINF\_pred          |   89.1007190 |   89.1007190 |          0 |
|      12 | AUCIFP   | AUCINF\_pred          |  130.6390680 |  130.6390680 |          0 |
|       1 | AUCIFPD  | AUCINF\_D\_pred       |    0.6769217 |    0.6769217 |          0 |
|       2 | AUCIFPD  | AUCINF\_D\_pred       |    0.3127010 |    0.3127010 |          0 |
|       3 | AUCIFPD  | AUCINF\_D\_pred       |    0.3424554 |    0.3424554 |          0 |
|       4 | AUCIFPD  | AUCINF\_D\_pred       |    0.3701361 |    0.3701361 |          0 |
|       5 | AUCIFPD  | AUCINF\_D\_pred       |    0.4351707 |    0.4351707 |          0 |
|       6 | AUCIFPD  | AUCINF\_D\_pred       |    0.2640522 |    0.2640522 |          0 |
|       7 | AUCIFPD  | AUCINF\_D\_pred       |    0.3246661 |    0.3246661 |          0 |
|       8 | AUCIFPD  | AUCINF\_D\_pred       |    0.3238845 |    0.3238845 |          0 |
|       9 | AUCIFPD  | AUCINF\_D\_pred       |    0.3120815 |    0.3120815 |          0 |
|      10 | AUCIFPD  | AUCINF\_D\_pred       |    0.5330247 |    0.5330247 |          0 |
|      11 | AUCIFPD  | AUCINF\_D\_pred       |    0.2784397 |    0.2784397 |          0 |
|      12 | AUCIFPD  | AUCINF\_D\_pred       |    0.4082471 |    0.4082471 |          0 |
|       1 | AUCPEP   | AUC\_.Extrap\_pred    |   31.2498763 |   31.2498763 |          0 |
|       2 | AUCPEP   | AUC\_.Extrap\_pred    |    8.5320300 |    8.5320300 |          0 |
|       3 | AUCPEP   | AUC\_.Extrap\_pred    |    9.3983245 |    9.3983245 |          0 |
|       4 | AUCPEP   | AUC\_.Extrap\_pred    |    9.8335939 |    9.8335939 |          0 |
|       5 | AUCPEP   | AUC\_.Extrap\_pred    |   12.8974027 |   12.8974027 |          0 |
|       6 | AUCPEP   | AUC\_.Extrap\_pred    |   12.6882455 |   12.6882455 |          0 |
|       7 | AUCPEP   | AUC\_.Extrap\_pred    |   12.6473665 |   12.6473664 |          0 |
|       8 | AUCPEP   | AUC\_.Extrap\_pred    |   14.5529307 |   14.5529307 |          0 |
|       9 | AUCPEP   | AUC\_.Extrap\_pred    |   13.5580763 |   13.5580763 |          0 |
|      10 | AUCPEP   | AUC\_.Extrap\_pred    |   18.8780012 |   18.8780012 |          0 |
|      11 | AUCPEP   | AUC\_.Extrap\_pred    |   10.1089184 |   10.1089184 |          0 |
|      12 | AUCPEP   | AUC\_.Extrap\_pred    |    8.1610870 |    8.1610870 |          0 |
|       1 | VZFP     | Vz\_F\_pred           |   30.4863228 |   30.4863228 |          0 |
|       2 | VZFP     | Vz\_F\_pred           |   30.7239161 |   30.7239161 |          0 |
|       3 | VZFP     | Vz\_F\_pred           |   28.5041534 |   28.5041534 |          0 |
|       4 | VZFP     | Vz\_F\_pred           |   27.2110972 |   27.2110972 |          0 |
|       5 | VZFP     | Vz\_F\_pred           |   26.5294196 |   26.5294196 |          0 |
|       6 | VZFP     | Vz\_F\_pred           |   43.1356944 |   43.1356944 |          0 |
|       7 | VZFP     | Vz\_F\_pred           |   34.8676684 |   34.8676684 |          0 |
|       8 | VZFP     | Vz\_F\_pred           |   37.9066862 |   37.9066862 |          0 |
|       9 | VZFP     | Vz\_F\_pred           |   38.8593822 |   38.8593822 |          0 |
|      10 | VZFP     | Vz\_F\_pred           |   25.0278813 |   25.0278813 |          0 |
|      11 | VZFP     | Vz\_F\_pred           |   37.6230406 |   37.6230406 |          0 |
|      12 | VZFP     | Vz\_F\_pred           |   22.2157473 |   22.2157473 |          0 |
|       1 | CLFP     | Cl\_F\_pred           |    1.4772757 |    1.4772757 |          0 |
|       2 | CLFP     | Cl\_F\_pred           |    3.1979432 |    3.1979432 |          0 |
|       3 | CLFP     | Cl\_F\_pred           |    2.9200884 |    2.9200884 |          0 |
|       4 | CLFP     | Cl\_F\_pred           |    2.7017088 |    2.7017088 |          0 |
|       5 | CLFP     | Cl\_F\_pred           |    2.2979487 |    2.2979487 |          0 |
|       6 | CLFP     | Cl\_F\_pred           |    3.7871302 |    3.7871302 |          0 |
|       7 | CLFP     | Cl\_F\_pred           |    3.0800877 |    3.0800877 |          0 |
|       8 | CLFP     | Cl\_F\_pred           |    3.0875201 |    3.0875201 |          0 |
|       9 | CLFP     | Cl\_F\_pred           |    3.2042916 |    3.2042916 |          0 |
|      10 | CLFP     | Cl\_F\_pred           |    1.8760856 |    1.8760856 |          0 |
|      11 | CLFP     | Cl\_F\_pred           |    3.5914413 |    3.5914413 |          0 |
|      12 | CLFP     | Cl\_F\_pred           |    2.4494970 |    2.4494970 |          0 |
|       1 | AUMCLST  | AUMClast              | 1459.0711035 | 1459.0711040 |          0 |
|       2 | AUMCLST  | AUMClast              |  706.5865660 |  706.5865660 |          0 |
|       3 | AUMCLST  | AUMClast              |  803.1858700 |  803.1858700 |          0 |
|       4 | AUMCLST  | AUMClast              |  901.0842105 |  901.0842105 |          0 |
|       5 | AUMCLST  | AUMClast              | 1017.1143165 | 1017.1143170 |          0 |
|       6 | AUMCLST  | AUMClast              |  609.1523875 |  609.1523875 |          0 |
|       7 | AUMCLST  | AUMClast              |  782.4198600 |  782.4198600 |          0 |
|       8 | AUMCLST  | AUMClast              |  739.5345980 |  739.5345980 |          0 |
|       9 | AUMCLST  | AUMClast              |  705.2296255 |  705.2296255 |          0 |
|      10 | AUMCLST  | AUMClast              | 1278.1800420 | 1278.1800420 |          0 |
|      11 | AUMCLST  | AUMClast              |  617.2422125 |  617.2422125 |          0 |
|      12 | AUMCLST  | AUMClast              |  977.8807235 |  977.8807235 |          0 |
|       1 | AUMCIFO  | AUMCINF\_obs          | 4505.5348194 | 4505.5348190 |          0 |
|       2 | AUMCIFO  | AUMCINF\_obs          |  999.7722880 |  999.7722880 |          0 |
|       3 | AUMCIFO  | AUMCINF\_obs          | 1150.9647687 | 1150.9647690 |          0 |
|       4 | AUMCIFO  | AUMCINF\_obs          | 1303.2524014 | 1303.2524010 |          0 |
|       5 | AUMCIFO  | AUMCINF\_obs          | 1667.7216119 | 1667.7216120 |          0 |
|       6 | AUMCIFO  | AUMCINF\_obs          |  978.4284857 |  978.4284857 |          0 |
|       7 | AUMCIFO  | AUMCINF\_obs          | 1245.0984083 | 1245.0984080 |          0 |
|       8 | AUMCIFO  | AUMCINF\_obs          | 1298.1157547 | 1298.1157550 |          0 |
|       9 | AUMCIFO  | AUMCINF\_obs          | 1201.7715381 | 1201.7715380 |          0 |
|      10 | AUMCIFO  | AUMCINF\_obs          | 2473.9934274 | 2473.9934270 |          0 |
|      11 | AUMCIFO  | AUMCINF\_obs          |  928.5599714 |  928.5599714 |          0 |
|      12 | AUMCIFO  | AUMCINF\_obs          | 1330.3840024 | 1330.3840020 |          0 |
|       1 | AUMCPEO  | AUMC\_.Extrap\_obs    |   67.6160287 |   67.6160287 |          0 |
|       2 | AUMCPEO  | AUMC\_.Extrap\_obs    |   29.3252499 |   29.3252499 |          0 |
|       3 | AUMCPEO  | AUMC\_.Extrap\_obs    |   30.2162940 |   30.2162940 |          0 |
|       4 | AUMCPEO  | AUMC\_.Extrap\_obs    |   30.8588107 |   30.8588107 |          0 |
|       5 | AUMCPEO  | AUMC\_.Extrap\_obs    |   39.0117446 |   39.0117446 |          0 |
|       6 | AUMCPEO  | AUMC\_.Extrap\_obs    |   37.7417567 |   37.7417567 |          0 |
|       7 | AUMCPEO  | AUMC\_.Extrap\_obs    |   37.1599984 |   37.1599984 |          0 |
|       8 | AUMCPEO  | AUMC\_.Extrap\_obs    |   43.0301500 |   43.0301500 |          0 |
|       9 | AUMCPEO  | AUMC\_.Extrap\_obs    |   41.3174965 |   41.3174965 |          0 |
|      10 | AUMCPEO  | AUMC\_.Extrap\_obs    |   48.3353501 |   48.3353501 |          0 |
|      11 | AUMCPEO  | AUMC\_.Extrap\_obs    |   33.5269416 |   33.5269415 |          0 |
|      12 | AUMCPEO  | AUMC\_.Extrap\_obs    |   26.4963558 |   26.4963558 |          0 |
|       1 | AUMCIFP  | AUMCINF\_pred         | 4505.6708646 | 4505.6708650 |          0 |
|       2 | AUMCIFP  | AUMCINF\_pred         |  996.0715835 |  996.0715835 |          0 |
|       3 | AUMCIFP  | AUMCINF\_pred         | 1152.6528903 | 1152.6528900 |          0 |
|       4 | AUMCIFP  | AUMCINF\_pred         | 1305.4981092 | 1305.4981090 |          0 |
|       5 | AUMCIFP  | AUMCINF\_pred         | 1661.7936744 | 1661.7936740 |          0 |
|       6 | AUMCIFP  | AUMCINF\_pred         |  986.9664597 |  986.9664597 |          0 |
|       7 | AUMCIFP  | AUMCINF\_pred         | 1249.4110601 | 1249.4110600 |          0 |
|       8 | AUMCIFP  | AUMCINF\_pred         | 1288.5201162 | 1288.5201160 |          0 |
|       9 | AUMCIFP  | AUMCINF\_pred         | 1200.2123597 | 1200.2123600 |          0 |
|      10 | AUMCIFP  | AUMCINF\_pred         | 2470.8765418 | 2470.8765420 |          0 |
|      11 | AUMCIFP  | AUMCINF\_pred         |  928.4899636 |  928.4899636 |          0 |
|      12 | AUMCIFP  | AUMCINF\_pred         | 1332.0528341 | 1332.0528340 |          0 |
|       1 | AUMCPEP  | AUMC\_.Extrap\_pred   |   67.6170065 |   67.6170065 |          0 |
|       2 | AUMCPEP  | AUMC\_.Extrap\_pred   |   29.0626720 |   29.0626720 |          0 |
|       3 | AUMCPEP  | AUMC\_.Extrap\_pred   |   30.3184960 |   30.3184960 |          0 |
|       4 | AUMCPEP  | AUMC\_.Extrap\_pred   |   30.9777468 |   30.9777468 |          0 |
|       5 | AUMCPEP  | AUMC\_.Extrap\_pred   |   38.7941877 |   38.7941877 |          0 |
|       6 | AUMCPEP  | AUMC\_.Extrap\_pred   |   38.2803355 |   38.2803355 |          0 |
|       7 | AUMCPEP  | AUMC\_.Extrap\_pred   |   37.3769062 |   37.3769062 |          0 |
|       8 | AUMCPEP  | AUMC\_.Extrap\_pred   |   42.6058943 |   42.6058943 |          0 |
|       9 | AUMCPEP  | AUMC\_.Extrap\_pred   |   41.2412629 |   41.2412629 |          0 |
|      10 | AUMCPEP  | AUMC\_.Extrap\_pred   |   48.2701778 |   48.2701778 |          0 |
|      11 | AUMCPEP  | AUMC\_.Extrap\_pred   |   33.5219295 |   33.5219295 |          0 |
|      12 | AUMCPEP  | AUMC\_.Extrap\_pred   |   26.5884432 |   26.5884432 |          0 |
|       1 | MRTEVLST | MRTlast               |    9.7974834 |    9.7974834 |          0 |
|       2 | MRTEVLST | MRTlast               |    7.7199964 |    7.7199964 |          0 |
|       3 | MRTEVLST | MRTlast               |    8.0895778 |    8.0895778 |          0 |
|       4 | MRTEVLST | MRTlast               |    8.4374104 |    8.4374104 |          0 |
|       5 | MRTEVLST | MRTlast               |    8.3855010 |    8.3855010 |          0 |
|       6 | MRTEVLST | MRTlast               |    8.2568329 |    8.2568329 |          0 |
|       7 | MRTEVLST | MRTlast               |    8.6213834 |    8.6213834 |          0 |
|       8 | MRTEVLST | MRTlast               |    8.3506664 |    8.3506664 |          0 |
|       9 | MRTEVLST | MRTlast               |    8.1693626 |    8.1693627 |          0 |
|      10 | MRTEVLST | MRTlast               |    9.2375341 |    9.2375341 |          0 |
|      11 | MRTEVLST | MRTlast               |    7.7065110 |    7.7065110 |          0 |
|      12 | MRTEVLST | MRTlast               |    8.1505343 |    8.1505343 |          0 |
|       1 | MRTEVIFO | MRTINF\_obs           |   20.8000305 |   20.8000305 |          0 |
|       2 | MRTEVIFO | MRTINF\_obs           |    9.9804109 |    9.9804109 |          0 |
|       3 | MRTEVIFO | MRTINF\_obs           |   10.5076420 |   10.5076420 |          0 |
|       4 | MRTEVIFO | MRTINF\_obs           |   11.0091630 |   11.0091630 |          0 |
|       5 | MRTEVIFO | MRTINF\_obs           |   11.9618725 |   11.9618725 |          0 |
|       6 | MRTEVIFO | MRTINF\_obs           |   11.6127855 |   11.6127855 |          0 |
|       7 | MRTEVIFO | MRTINF\_obs           |   11.9984272 |   11.9984272 |          0 |
|       8 | MRTEVIFO | MRTINF\_obs           |   12.4930916 |   12.4930916 |          0 |
|       9 | MRTEVIFO | MRTINF\_obs           |   12.0286954 |   12.0286954 |          0 |
|      10 | MRTEVIFO | MRTINF\_obs           |   14.4972959 |   14.4972959 |          0 |
|      11 | MRTEVIFO | MRTINF\_obs           |   10.4212275 |   10.4212274 |          0 |
|      12 | MRTEVIFO | MRTINF\_obs           |   10.1875787 |   10.1875787 |          0 |
|       1 | MRTEVIFP | MRTINF\_pred          |   20.8003683 |   20.8003683 |          0 |
|       2 | MRTEVIFP | MRTINF\_pred          |    9.9543135 |    9.9543135 |          0 |
|       3 | MRTEVIFP | MRTINF\_pred          |   10.5182762 |   10.5182762 |          0 |
|       4 | MRTEVIFP | MRTINF\_pred          |   11.0221115 |   11.0221115 |          0 |
|       5 | MRTEVIFP | MRTINF\_pred          |   11.9334895 |   11.9334895 |          0 |
|       6 | MRTEVIFP | MRTINF\_pred          |   11.6805328 |   11.6805328 |          0 |
|       7 | MRTEVIFP | MRTINF\_pred          |   12.0259237 |   12.0259237 |          0 |
|       8 | MRTEVIFP | MRTINF\_pred          |   12.4322866 |   12.4322866 |          0 |
|       9 | MRTEVIFP | MRTINF\_pred          |   12.0182199 |   12.0182199 |          0 |
|      10 | MRTEVIFP | MRTINF\_pred          |   14.4861745 |   14.4861745 |          0 |
|      11 | MRTEVIFP | MRTINF\_pred          |   10.4206787 |   10.4206787 |          0 |
|      12 | MRTEVIFP | MRTINF\_pred          |   10.1964355 |   10.1964355 |          0 |

Theoph (n=12), Linear, Extravascular

## Test 2: Theoph (n=12), Log, Extravascular

``` r
table_wres_rres(Wres2, Rres2,
                Caption = 'Theoph (n=12), Log, Extravascular')
```

| Subject | PPTESTCD | WNL                   |   NonCompart |    WinNonlin | Difference |
| ------: | :------- | :-------------------- | -----------: | -----------: | ---------: |
|       1 | R2       | Rsq                   |    0.9999997 |    0.9999997 |          0 |
|       2 | R2       | Rsq                   |    0.9971954 |    0.9971954 |          0 |
|       3 | R2       | Rsq                   |    0.9993250 |    0.9993250 |          0 |
|       4 | R2       | Rsq                   |    0.9989241 |    0.9989241 |          0 |
|       5 | R2       | Rsq                   |    0.9986472 |    0.9986472 |          0 |
|       6 | R2       | Rsq                   |    0.9982413 |    0.9982413 |          0 |
|       7 | R2       | Rsq                   |    0.9986702 |    0.9986702 |          0 |
|       8 | R2       | Rsq                   |    0.9910124 |    0.9910124 |          0 |
|       9 | R2       | Rsq                   |    0.9994437 |    0.9994437 |          0 |
|      10 | R2       | Rsq                   |    0.9995087 |    0.9995087 |          0 |
|      11 | R2       | Rsq                   |    0.9999983 |    0.9999983 |          0 |
|      12 | R2       | Rsq                   |    0.9993968 |    0.9993968 |          0 |
|       1 | R2ADJ    | Rsq\_adjusted         |    0.9999995 |    0.9999995 |          0 |
|       2 | R2ADJ    | Rsq\_adjusted         |    0.9957931 |    0.9957931 |          0 |
|       3 | R2ADJ    | Rsq\_adjusted         |    0.9986499 |    0.9986499 |          0 |
|       4 | R2ADJ    | Rsq\_adjusted         |    0.9978483 |    0.9978483 |          0 |
|       5 | R2ADJ    | Rsq\_adjusted         |    0.9979708 |    0.9979708 |          0 |
|       6 | R2ADJ    | Rsq\_adjusted         |    0.9978896 |    0.9978896 |          0 |
|       7 | R2ADJ    | Rsq\_adjusted         |    0.9980053 |    0.9980053 |          0 |
|       8 | R2ADJ    | Rsq\_adjusted         |    0.9887655 |    0.9887655 |          0 |
|       9 | R2ADJ    | Rsq\_adjusted         |    0.9988873 |    0.9988873 |          0 |
|      10 | R2ADJ    | Rsq\_adjusted         |    0.9990174 |    0.9990174 |          0 |
|      11 | R2ADJ    | Rsq\_adjusted         |    0.9999965 |    0.9999965 |          0 |
|      12 | R2ADJ    | Rsq\_adjusted         |    0.9987936 |    0.9987936 |          0 |
|       1 | CORRXY   | Corr\_XY              |  \-0.9999999 |  \-0.9999999 |          0 |
|       2 | CORRXY   | Corr\_XY              |  \-0.9985967 |  \-0.9985967 |          0 |
|       3 | CORRXY   | Corr\_XY              |  \-0.9996624 |  \-0.9996624 |          0 |
|       4 | CORRXY   | Corr\_XY              |  \-0.9994619 |  \-0.9994619 |          0 |
|       5 | CORRXY   | Corr\_XY              |  \-0.9993234 |  \-0.9993234 |          0 |
|       6 | CORRXY   | Corr\_XY              |  \-0.9991203 |  \-0.9991203 |          0 |
|       7 | CORRXY   | Corr\_XY              |  \-0.9993349 |  \-0.9993349 |          0 |
|       8 | CORRXY   | Corr\_XY              |  \-0.9954961 |  \-0.9954961 |          0 |
|       9 | CORRXY   | Corr\_XY              |  \-0.9997218 |  \-0.9997218 |          0 |
|      10 | CORRXY   | Corr\_XY              |  \-0.9997543 |  \-0.9997543 |          0 |
|      11 | CORRXY   | Corr\_XY              |  \-0.9999991 |  \-0.9999991 |          0 |
|      12 | CORRXY   | Corr\_XY              |  \-0.9996984 |  \-0.9996984 |          0 |
|       1 | LAMZNPT  | No\_points\_lambda\_z |    3.0000000 |    3.0000000 |          0 |
|       2 | LAMZNPT  | No\_points\_lambda\_z |    4.0000000 |    4.0000000 |          0 |
|       3 | LAMZNPT  | No\_points\_lambda\_z |    3.0000000 |    3.0000000 |          0 |
|       4 | LAMZNPT  | No\_points\_lambda\_z |    3.0000000 |    3.0000000 |          0 |
|       5 | LAMZNPT  | No\_points\_lambda\_z |    4.0000000 |    4.0000000 |          0 |
|       6 | LAMZNPT  | No\_points\_lambda\_z |    7.0000000 |    7.0000000 |          0 |
|       7 | LAMZNPT  | No\_points\_lambda\_z |    4.0000000 |    4.0000000 |          0 |
|       8 | LAMZNPT  | No\_points\_lambda\_z |    6.0000000 |    6.0000000 |          0 |
|       9 | LAMZNPT  | No\_points\_lambda\_z |    3.0000000 |    3.0000000 |          0 |
|      10 | LAMZNPT  | No\_points\_lambda\_z |    3.0000000 |    3.0000000 |          0 |
|      11 | LAMZNPT  | No\_points\_lambda\_z |    3.0000000 |    3.0000000 |          0 |
|      12 | LAMZNPT  | No\_points\_lambda\_z |    3.0000000 |    3.0000000 |          0 |
|       1 | LAMZ     | Lambda\_z             |    0.0484570 |    0.0484570 |          0 |
|       2 | LAMZ     | Lambda\_z             |    0.1040864 |    0.1040864 |          0 |
|       3 | LAMZ     | Lambda\_z             |    0.1024443 |    0.1024443 |          0 |
|       4 | LAMZ     | Lambda\_z             |    0.0992870 |    0.0992870 |          0 |
|       5 | LAMZ     | Lambda\_z             |    0.0866189 |    0.0866189 |          0 |
|       6 | LAMZ     | Lambda\_z             |    0.0877957 |    0.0877957 |          0 |
|       7 | LAMZ     | Lambda\_z             |    0.0883365 |    0.0883365 |          0 |
|       8 | LAMZ     | Lambda\_z             |    0.0814505 |    0.0814505 |          0 |
|       9 | LAMZ     | Lambda\_z             |    0.0824586 |    0.0824586 |          0 |
|      10 | LAMZ     | Lambda\_z             |    0.0749598 |    0.0749598 |          0 |
|      11 | LAMZ     | Lambda\_z             |    0.0954586 |    0.0954586 |          0 |
|      12 | LAMZ     | Lambda\_z             |    0.1102595 |    0.1102595 |          0 |
|       1 | LAMZLL   | Lambda\_z\_lower      |    9.0500000 |    9.0500000 |          0 |
|       2 | LAMZLL   | Lambda\_z\_lower      |    7.0300000 |    7.0300000 |          0 |
|       3 | LAMZLL   | Lambda\_z\_lower      |    9.0000000 |    9.0000000 |          0 |
|       4 | LAMZLL   | Lambda\_z\_lower      |    9.0200000 |    9.0200000 |          0 |
|       5 | LAMZLL   | Lambda\_z\_lower      |    7.0200000 |    7.0200000 |          0 |
|       6 | LAMZLL   | Lambda\_z\_lower      |    2.0300000 |    2.0300000 |          0 |
|       7 | LAMZLL   | Lambda\_z\_lower      |    6.9800000 |    6.9800000 |          0 |
|       8 | LAMZLL   | Lambda\_z\_lower      |    3.5300000 |    3.5300000 |          0 |
|       9 | LAMZLL   | Lambda\_z\_lower      |    8.8000000 |    8.8000000 |          0 |
|      10 | LAMZLL   | Lambda\_z\_lower      |    9.3800000 |    9.3800000 |          0 |
|      11 | LAMZLL   | Lambda\_z\_lower      |    9.0300000 |    9.0300000 |          0 |
|      12 | LAMZLL   | Lambda\_z\_lower      |    9.0300000 |    9.0300000 |          0 |
|       1 | LAMZUL   | Lambda\_z\_upper      |   24.3700000 |   24.3700000 |          0 |
|       2 | LAMZUL   | Lambda\_z\_upper      |   24.3000000 |   24.3000000 |          0 |
|       3 | LAMZUL   | Lambda\_z\_upper      |   24.1700000 |   24.1700000 |          0 |
|       4 | LAMZUL   | Lambda\_z\_upper      |   24.6500000 |   24.6500000 |          0 |
|       5 | LAMZUL   | Lambda\_z\_upper      |   24.3500000 |   24.3500000 |          0 |
|       6 | LAMZUL   | Lambda\_z\_upper      |   23.8500000 |   23.8500000 |          0 |
|       7 | LAMZUL   | Lambda\_z\_upper      |   24.2200000 |   24.2200000 |          0 |
|       8 | LAMZUL   | Lambda\_z\_upper      |   24.1200000 |   24.1200000 |          0 |
|       9 | LAMZUL   | Lambda\_z\_upper      |   24.4300000 |   24.4300000 |          0 |
|      10 | LAMZUL   | Lambda\_z\_upper      |   23.7000000 |   23.7000000 |          0 |
|      11 | LAMZUL   | Lambda\_z\_upper      |   24.0800000 |   24.0800000 |          0 |
|      12 | LAMZUL   | Lambda\_z\_upper      |   24.1500000 |   24.1500000 |          0 |
|       1 | LAMZHL   | HL\_Lambda\_z         |   14.3043776 |   14.3043776 |          0 |
|       2 | LAMZHL   | HL\_Lambda\_z         |    6.6593416 |    6.6593416 |          0 |
|       3 | LAMZHL   | HL\_Lambda\_z         |    6.7660874 |    6.7660874 |          0 |
|       4 | LAMZHL   | HL\_Lambda\_z         |    6.9812467 |    6.9812467 |          0 |
|       5 | LAMZHL   | HL\_Lambda\_z         |    8.0022640 |    8.0022640 |          0 |
|       6 | LAMZHL   | HL\_Lambda\_z         |    7.8949979 |    7.8949979 |          0 |
|       7 | LAMZHL   | HL\_Lambda\_z         |    7.8466683 |    7.8466683 |          0 |
|       8 | LAMZHL   | HL\_Lambda\_z         |    8.5100379 |    8.5100379 |          0 |
|       9 | LAMZHL   | HL\_Lambda\_z         |    8.4059988 |    8.4059988 |          0 |
|      10 | LAMZHL   | HL\_Lambda\_z         |    9.2469158 |    9.2469158 |          0 |
|      11 | LAMZHL   | HL\_Lambda\_z         |    7.2612365 |    7.2612365 |          0 |
|      12 | LAMZHL   | HL\_Lambda\_z         |    6.2865082 |    6.2865082 |          0 |
|       1 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       2 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       3 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       4 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       5 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       6 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       7 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       8 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       9 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|      10 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|      11 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|      12 | TLAG     | Tlag                  |    0.0000000 |    0.0000000 |          0 |
|       1 | TMAX     | Tmax                  |    1.1200000 |    1.1200000 |          0 |
|       2 | TMAX     | Tmax                  |    1.9200000 |    1.9200000 |          0 |
|       3 | TMAX     | Tmax                  |    1.0200000 |    1.0200000 |          0 |
|       4 | TMAX     | Tmax                  |    1.0700000 |    1.0700000 |          0 |
|       5 | TMAX     | Tmax                  |    1.0000000 |    1.0000000 |          0 |
|       6 | TMAX     | Tmax                  |    1.1500000 |    1.1500000 |          0 |
|       7 | TMAX     | Tmax                  |    3.4800000 |    3.4800000 |          0 |
|       8 | TMAX     | Tmax                  |    2.0200000 |    2.0200000 |          0 |
|       9 | TMAX     | Tmax                  |    0.6300000 |    0.6300000 |          0 |
|      10 | TMAX     | Tmax                  |    3.5500000 |    3.5500000 |          0 |
|      11 | TMAX     | Tmax                  |    0.9800000 |    0.9800000 |          0 |
|      12 | TMAX     | Tmax                  |    3.5200000 |    3.5200000 |          0 |
|       1 | CMAX     | Cmax                  |   10.5000000 |   10.5000000 |          0 |
|       2 | CMAX     | Cmax                  |    8.3300000 |    8.3300000 |          0 |
|       3 | CMAX     | Cmax                  |    8.2000000 |    8.2000000 |          0 |
|       4 | CMAX     | Cmax                  |    8.6000000 |    8.6000000 |          0 |
|       5 | CMAX     | Cmax                  |   11.4000000 |   11.4000000 |          0 |
|       6 | CMAX     | Cmax                  |    6.4400000 |    6.4400000 |          0 |
|       7 | CMAX     | Cmax                  |    7.0900000 |    7.0900000 |          0 |
|       8 | CMAX     | Cmax                  |    7.5600000 |    7.5600000 |          0 |
|       9 | CMAX     | Cmax                  |    9.0300000 |    9.0300000 |          0 |
|      10 | CMAX     | Cmax                  |   10.2100000 |   10.2100000 |          0 |
|      11 | CMAX     | Cmax                  |    8.0000000 |    8.0000000 |          0 |
|      12 | CMAX     | Cmax                  |    9.7500000 |    9.7500000 |          0 |
|       1 | CMAXD    | Cmax\_D               |    0.0328125 |    0.0328125 |          0 |
|       2 | CMAXD    | Cmax\_D               |    0.0260312 |    0.0260312 |          0 |
|       3 | CMAXD    | Cmax\_D               |    0.0256250 |    0.0256250 |          0 |
|       4 | CMAXD    | Cmax\_D               |    0.0268750 |    0.0268750 |          0 |
|       5 | CMAXD    | Cmax\_D               |    0.0356250 |    0.0356250 |          0 |
|       6 | CMAXD    | Cmax\_D               |    0.0201250 |    0.0201250 |          0 |
|       7 | CMAXD    | Cmax\_D               |    0.0221562 |    0.0221562 |          0 |
|       8 | CMAXD    | Cmax\_D               |    0.0236250 |    0.0236250 |          0 |
|       9 | CMAXD    | Cmax\_D               |    0.0282188 |    0.0282188 |          0 |
|      10 | CMAXD    | Cmax\_D               |    0.0319063 |    0.0319062 |          0 |
|      11 | CMAXD    | Cmax\_D               |    0.0250000 |    0.0250000 |          0 |
|      12 | CMAXD    | Cmax\_D               |    0.0304688 |    0.0304688 |          0 |
|       1 | TLST     | Tlast                 |   24.3700000 |   24.3700000 |          0 |
|       2 | TLST     | Tlast                 |   24.3000000 |   24.3000000 |          0 |
|       3 | TLST     | Tlast                 |   24.1700000 |   24.1700000 |          0 |
|       4 | TLST     | Tlast                 |   24.6500000 |   24.6500000 |          0 |
|       5 | TLST     | Tlast                 |   24.3500000 |   24.3500000 |          0 |
|       6 | TLST     | Tlast                 |   23.8500000 |   23.8500000 |          0 |
|       7 | TLST     | Tlast                 |   24.2200000 |   24.2200000 |          0 |
|       8 | TLST     | Tlast                 |   24.1200000 |   24.1200000 |          0 |
|       9 | TLST     | Tlast                 |   24.4300000 |   24.4300000 |          0 |
|      10 | TLST     | Tlast                 |   23.7000000 |   23.7000000 |          0 |
|      11 | TLST     | Tlast                 |   24.0800000 |   24.0800000 |          0 |
|      12 | TLST     | Tlast                 |   24.1500000 |   24.1500000 |          0 |
|       1 | CLST     | Clast                 |    3.2800000 |    3.2800000 |          0 |
|       2 | CLST     | Clast                 |    0.9000000 |    0.9000000 |          0 |
|       3 | CLST     | Clast                 |    1.0500000 |    1.0500000 |          0 |
|       4 | CLST     | Clast                 |    1.1500000 |    1.1500000 |          0 |
|       5 | CLST     | Clast                 |    1.5700000 |    1.5700000 |          0 |
|       6 | CLST     | Clast                 |    0.9200000 |    0.9200000 |          0 |
|       7 | CLST     | Clast                 |    1.1500000 |    1.1500000 |          0 |
|       8 | CLST     | Clast                 |    1.2500000 |    1.2500000 |          0 |
|       9 | CLST     | Clast                 |    1.1200000 |    1.1200000 |          0 |
|      10 | CLST     | Clast                 |    2.4200000 |    2.4200000 |          0 |
|      11 | CLST     | Clast                 |    0.8600000 |    0.8600000 |          0 |
|      12 | CLST     | Clast                 |    1.1700000 |    1.1700000 |          0 |
|       1 | AUCLST   | AUClast               |  147.2347485 |  147.2347485 |          0 |
|       2 | AUCLST   | AUClast               |   88.7312755 |   88.7312755 |          0 |
|       3 | AUCLST   | AUClast               |   95.8781978 |   95.8781978 |          0 |
|       4 | AUCLST   | AUClast               |  102.6336232 |  102.6336232 |          0 |
|       5 | AUCLST   | AUClast               |  118.1793538 |  118.1793538 |          0 |
|       6 | AUCLST   | AUClast               |   71.6970150 |   71.6970150 |          0 |
|       7 | AUCLST   | AUClast               |   87.9692274 |   87.9692274 |          0 |
|       8 | AUCLST   | AUClast               |   86.8065635 |   86.8065635 |          0 |
|       9 | AUCLST   | AUClast               |   83.9374360 |   83.9374360 |          0 |
|      10 | AUCLST   | AUClast               |  135.5760701 |  135.5760701 |          0 |
|      11 | AUCLST   | AUClast               |   77.8934723 |   77.8934723 |          0 |
|      12 | AUCLST   | AUClast               |  115.2202082 |  115.2202082 |          0 |
|       1 | AUCALL   | AUCall                |  147.2347485 |  147.2347485 |          0 |
|       2 | AUCALL   | AUCall                |   88.7312755 |   88.7312755 |          0 |
|       3 | AUCALL   | AUCall                |   95.8781978 |   95.8781978 |          0 |
|       4 | AUCALL   | AUCall                |  102.6336232 |  102.6336232 |          0 |
|       5 | AUCALL   | AUCall                |  118.1793538 |  118.1793538 |          0 |
|       6 | AUCALL   | AUCall                |   71.6970150 |   71.6970150 |          0 |
|       7 | AUCALL   | AUCall                |   87.9692274 |   87.9692274 |          0 |
|       8 | AUCALL   | AUCall                |   86.8065635 |   86.8065635 |          0 |
|       9 | AUCALL   | AUCall                |   83.9374360 |   83.9374360 |          0 |
|      10 | AUCALL   | AUCall                |  135.5760701 |  135.5760701 |          0 |
|      11 | AUCALL   | AUCall                |   77.8934723 |   77.8934723 |          0 |
|      12 | AUCALL   | AUCall                |  115.2202082 |  115.2202082 |          0 |
|       1 | AUCIFO   | AUCINF\_obs           |  214.9236316 |  214.9236316 |          0 |
|       2 | AUCIFO   | AUCINF\_obs           |   97.3779346 |   97.3779346 |          0 |
|       3 | AUCIFO   | AUCINF\_obs           |  106.1276685 |  106.1276685 |          0 |
|       4 | AUCIFO   | AUCINF\_obs           |  114.2162046 |  114.2162046 |          0 |
|       5 | AUCIFO   | AUCINF\_obs           |  136.3047316 |  136.3047316 |          0 |
|       6 | AUCIFO   | AUCINF\_obs           |   82.1758833 |   82.1758833 |          0 |
|       7 | AUCIFO   | AUCINF\_obs           |  100.9876292 |  100.9876292 |          0 |
|       8 | AUCIFO   | AUCINF\_obs           |  102.1533003 |  102.1533003 |          0 |
|       9 | AUCIFO   | AUCINF\_obs           |   97.5200039 |   97.5200039 |          0 |
|      10 | AUCIFO   | AUCINF\_obs           |  167.8600307 |  167.8600307 |          0 |
|      11 | AUCIFO   | AUCINF\_obs           |   86.9026173 |   86.9026173 |          0 |
|      12 | AUCIFO   | AUCINF\_obs           |  125.8315397 |  125.8315397 |          0 |
|       1 | AUCIFOD  | AUCINF\_D\_obs        |    0.6716363 |    0.6716363 |          0 |
|       2 | AUCIFOD  | AUCINF\_D\_obs        |    0.3043060 |    0.3043060 |          0 |
|       3 | AUCIFOD  | AUCINF\_D\_obs        |    0.3316490 |    0.3316490 |          0 |
|       4 | AUCIFOD  | AUCINF\_D\_obs        |    0.3569256 |    0.3569256 |          0 |
|       5 | AUCIFOD  | AUCINF\_D\_obs        |    0.4259523 |    0.4259523 |          0 |
|       6 | AUCIFOD  | AUCINF\_D\_obs        |    0.2567996 |    0.2567996 |          0 |
|       7 | AUCIFOD  | AUCINF\_D\_obs        |    0.3155863 |    0.3155863 |          0 |
|       8 | AUCIFOD  | AUCINF\_D\_obs        |    0.3192291 |    0.3192291 |          0 |
|       9 | AUCIFOD  | AUCINF\_D\_obs        |    0.3047500 |    0.3047500 |          0 |
|      10 | AUCIFOD  | AUCINF\_D\_obs        |    0.5245626 |    0.5245626 |          0 |
|      11 | AUCIFOD  | AUCINF\_D\_obs        |    0.2715707 |    0.2715707 |          0 |
|      12 | AUCIFOD  | AUCINF\_D\_obs        |    0.3932236 |    0.3932236 |          0 |
|       1 | AUCPEO   | AUC\_.Extrap\_obs     |   31.4943883 |   31.4943883 |          0 |
|       2 | AUCPEO   | AUC\_.Extrap\_obs     |    8.8794850 |    8.8794850 |          0 |
|       3 | AUCPEO   | AUC\_.Extrap\_obs     |    9.6576801 |    9.6576801 |          0 |
|       4 | AUCPEO   | AUC\_.Extrap\_obs     |   10.1409266 |   10.1409266 |          0 |
|       5 | AUCPEO   | AUC\_.Extrap\_obs     |   13.2976879 |   13.2976879 |          0 |
|       6 | AUCPEO   | AUC\_.Extrap\_obs     |   12.7517562 |   12.7517562 |          0 |
|       7 | AUCPEO   | AUC\_.Extrap\_obs     |   12.8910857 |   12.8910857 |          0 |
|       8 | AUCPEO   | AUC\_.Extrap\_obs     |   15.0232413 |   15.0232413 |          0 |
|       9 | AUCPEO   | AUC\_.Extrap\_obs     |   13.9279813 |   13.9279813 |          0 |
|      10 | AUCPEO   | AUC\_.Extrap\_obs     |   19.2326669 |   19.2326669 |          0 |
|      11 | AUCPEO   | AUC\_.Extrap\_obs     |   10.3669431 |   10.3669432 |          0 |
|      12 | AUCPEO   | AUC\_.Extrap\_obs     |    8.4329665 |    8.4329665 |          0 |
|       1 | VZFO     | Vz\_F\_obs            |   30.7262325 |   30.7262325 |          0 |
|       2 | VZFO     | Vz\_F\_obs            |   31.5715024 |   31.5715024 |          0 |
|       3 | VZFO     | Vz\_F\_obs            |   29.4329299 |   29.4329299 |          0 |
|       4 | VZFO     | Vz\_F\_obs            |   28.2182304 |   28.2182304 |          0 |
|       5 | VZFO     | Vz\_F\_obs            |   27.1035678 |   27.1035677 |          0 |
|       6 | VZFO     | Vz\_F\_obs            |   44.3539348 |   44.3539348 |          0 |
|       7 | VZFO     | Vz\_F\_obs            |   35.8708471 |   35.8708471 |          0 |
|       8 | VZFO     | Vz\_F\_obs            |   38.4594978 |   38.4594978 |          0 |
|       9 | VZFO     | Vz\_F\_obs            |   39.7942323 |   39.7942323 |          0 |
|      10 | VZFO     | Vz\_F\_obs            |   25.4316257 |   25.4316257 |          0 |
|      11 | VZFO     | Vz\_F\_obs            |   38.5746722 |   38.5746722 |          0 |
|      12 | VZFO     | Vz\_F\_obs            |   23.0645237 |   23.0645237 |          0 |
|       1 | CLFO     | Cl\_F\_obs            |    1.4889010 |    1.4889010 |          0 |
|       2 | CLFO     | Cl\_F\_obs            |    3.2861654 |    3.2861654 |          0 |
|       3 | CLFO     | Cl\_F\_obs            |    3.0152363 |    3.0152363 |          0 |
|       4 | CLFO     | Cl\_F\_obs            |    2.8017040 |    2.8017040 |          0 |
|       5 | CLFO     | Cl\_F\_obs            |    2.3476808 |    2.3476808 |          0 |
|       6 | CLFO     | Cl\_F\_obs            |    3.8940865 |    3.8940865 |          0 |
|       7 | CLFO     | Cl\_F\_obs            |    3.1687049 |    3.1687049 |          0 |
|       8 | CLFO     | Cl\_F\_obs            |    3.1325469 |    3.1325469 |          0 |
|       9 | CLFO     | Cl\_F\_obs            |    3.2813780 |    3.2813780 |          0 |
|      10 | CLFO     | Cl\_F\_obs            |    1.9063502 |    1.9063502 |          0 |
|      11 | CLFO     | Cl\_F\_obs            |    3.6822827 |    3.6822827 |          0 |
|      12 | CLFO     | Cl\_F\_obs            |    2.5430826 |    2.5430826 |          0 |
|       1 | AUCIFP   | AUCINF\_pred          |  214.9266543 |  214.9266543 |          0 |
|       2 | AUCIFP   | AUCINF\_pred          |   97.2687931 |   97.2687931 |          0 |
|       3 | AUCIFP   | AUCINF\_pred          |  106.1774195 |  106.1774195 |          0 |
|       4 | AUCIFP   | AUCINF\_pred          |  114.2808818 |  114.2808818 |          0 |
|       5 | AUCIFP   | AUCINF\_pred          |  136.1395842 |  136.1395842 |          0 |
|       6 | AUCIFP   | AUCINF\_pred          |   82.4181636 |   82.4181636 |          0 |
|       7 | AUCIFP   | AUCINF\_pred          |  101.1089745 |  101.1089745 |          0 |
|       8 | AUCIFP   | AUCINF\_pred          |  101.8896649 |  101.8896649 |          0 |
|       9 | AUCIFP   | AUCINF\_pred          |   97.4773537 |   97.4773537 |          0 |
|      10 | AUCIFP   | AUCINF\_pred          |  167.7758826 |  167.7758826 |          0 |
|      11 | AUCIFP   | AUCINF\_pred          |   86.9005913 |   86.9005913 |          0 |
|      12 | AUCIFP   | AUCINF\_pred          |  125.8817762 |  125.8817762 |          0 |
|       1 | AUCIFPD  | AUCINF\_D\_pred       |    0.6716458 |    0.6716458 |          0 |
|       2 | AUCIFPD  | AUCINF\_D\_pred       |    0.3039650 |    0.3039650 |          0 |
|       3 | AUCIFPD  | AUCINF\_D\_pred       |    0.3318044 |    0.3318044 |          0 |
|       4 | AUCIFPD  | AUCINF\_D\_pred       |    0.3571278 |    0.3571278 |          0 |
|       5 | AUCIFPD  | AUCINF\_D\_pred       |    0.4254362 |    0.4254362 |          0 |
|       6 | AUCIFPD  | AUCINF\_D\_pred       |    0.2575568 |    0.2575568 |          0 |
|       7 | AUCIFPD  | AUCINF\_D\_pred       |    0.3159655 |    0.3159655 |          0 |
|       8 | AUCIFPD  | AUCINF\_D\_pred       |    0.3184052 |    0.3184052 |          0 |
|       9 | AUCIFPD  | AUCINF\_D\_pred       |    0.3046167 |    0.3046167 |          0 |
|      10 | AUCIFPD  | AUCINF\_D\_pred       |    0.5242996 |    0.5242996 |          0 |
|      11 | AUCIFPD  | AUCINF\_D\_pred       |    0.2715643 |    0.2715643 |          0 |
|      12 | AUCIFPD  | AUCINF\_D\_pred       |    0.3933806 |    0.3933806 |          0 |
|       1 | AUCPEP   | AUC\_.Extrap\_pred    |   31.4953518 |   31.4953518 |          0 |
|       2 | AUCPEP   | AUC\_.Extrap\_pred    |    8.7772423 |    8.7772423 |          0 |
|       3 | AUCPEP   | AUC\_.Extrap\_pred    |    9.7000114 |    9.7000114 |          0 |
|       4 | AUCPEP   | AUC\_.Extrap\_pred    |   10.1917822 |   10.1917822 |          0 |
|       5 | AUCPEP   | AUC\_.Extrap\_pred    |   13.1925116 |   13.1925116 |          0 |
|       6 | AUCPEP   | AUC\_.Extrap\_pred    |   13.0082352 |   13.0082352 |          0 |
|       7 | AUCPEP   | AUC\_.Extrap\_pred    |   12.9956288 |   12.9956288 |          0 |
|       8 | AUCPEP   | AUC\_.Extrap\_pred    |   14.8033674 |   14.8033674 |          0 |
|       9 | AUCPEP   | AUC\_.Extrap\_pred    |   13.8903213 |   13.8903213 |          0 |
|      10 | AUCPEP   | AUC\_.Extrap\_pred    |   19.1921580 |   19.1921580 |          0 |
|      11 | AUCPEP   | AUC\_.Extrap\_pred    |   10.3648535 |   10.3648535 |          0 |
|      12 | AUCPEP   | AUC\_.Extrap\_pred    |    8.4695087 |    8.4695087 |          0 |
|       1 | VZFP     | Vz\_F\_pred           |   30.7258003 |   30.7258003 |          0 |
|       2 | VZFP     | Vz\_F\_pred           |   31.6069275 |   31.6069275 |          0 |
|       3 | VZFP     | Vz\_F\_pred           |   29.4191386 |   29.4191386 |          0 |
|       4 | VZFP     | Vz\_F\_pred           |   28.2022603 |   28.2022603 |          0 |
|       5 | VZFP     | Vz\_F\_pred           |   27.1364464 |   27.1364464 |          0 |
|       6 | VZFP     | Vz\_F\_pred           |   44.2235499 |   44.2235499 |          0 |
|       7 | VZFP     | Vz\_F\_pred           |   35.8277969 |   35.8277969 |          0 |
|       8 | VZFP     | Vz\_F\_pred           |   38.5590101 |   38.5590101 |          0 |
|       9 | VZFP     | Vz\_F\_pred           |   39.8116439 |   39.8116439 |          0 |
|      10 | VZFP     | Vz\_F\_pred           |   25.4443810 |   25.4443809 |          0 |
|      11 | VZFP     | Vz\_F\_pred           |   38.5755715 |   38.5755715 |          0 |
|      12 | VZFP     | Vz\_F\_pred           |   23.0553192 |   23.0553192 |          0 |
|       1 | CLFP     | Cl\_F\_pred           |    1.4888800 |    1.4888800 |          0 |
|       2 | CLFP     | Cl\_F\_pred           |    3.2898527 |    3.2898527 |          0 |
|       3 | CLFP     | Cl\_F\_pred           |    3.0138235 |    3.0138235 |          0 |
|       4 | CLFP     | Cl\_F\_pred           |    2.8001184 |    2.8001184 |          0 |
|       5 | CLFP     | Cl\_F\_pred           |    2.3505287 |    2.3505287 |          0 |
|       6 | CLFP     | Cl\_F\_pred           |    3.8826393 |    3.8826393 |          0 |
|       7 | CLFP     | Cl\_F\_pred           |    3.1649020 |    3.1649020 |          0 |
|       8 | CLFP     | Cl\_F\_pred           |    3.1406522 |    3.1406522 |          0 |
|       9 | CLFP     | Cl\_F\_pred           |    3.2828138 |    3.2828138 |          0 |
|      10 | CLFP     | Cl\_F\_pred           |    1.9073063 |    1.9073063 |          0 |
|      11 | CLFP     | Cl\_F\_pred           |    3.6823685 |    3.6823685 |          0 |
|      12 | CLFP     | Cl\_F\_pred           |    2.5420677 |    2.5420677 |          0 |
|       1 | AUMCLST  | AUMClast              | 1499.1290852 | 1499.1290850 |          0 |
|       2 | AUMCLST  | AUMClast              |  716.2787279 |  716.2787279 |          0 |
|       3 | AUMCLST  | AUMClast              |  810.8726830 |  810.8726830 |          0 |
|       4 | AUMCLST  | AUMClast              |  911.7828093 |  911.7828093 |          0 |
|       5 | AUMCLST  | AUMClast              | 1038.8799844 | 1038.8799840 |          0 |
|       6 | AUMCLST  | AUMClast              |  618.6659191 |  618.6659191 |          0 |
|       7 | AUMCLST  | AUMClast              |  795.6267785 |  795.6267785 |          0 |
|       8 | AUMCLST  | AUMClast              |  756.3619816 |  756.3619816 |          0 |
|       9 | AUMCLST  | AUMClast              |  723.3794155 |  723.3794155 |          0 |
|      10 | AUMCLST  | AUMClast              | 1306.7406149 | 1306.7406150 |          0 |
|      11 | AUMCLST  | AUMClast              |  626.6357849 |  626.6357849 |          0 |
|      12 | AUMCLST  | AUMClast              |  982.6343023 |  982.6343023 |          0 |
|       1 | AUMCIFO  | AUMCINF\_obs          | 4545.5928011 | 4545.5928010 |          0 |
|       2 | AUMCIFO  | AUMCINF\_obs          | 1009.4644499 | 1009.4644500 |          0 |
|       3 | AUMCIFO  | AUMCINF\_obs          | 1158.6515817 | 1158.6515820 |          0 |
|       4 | AUMCIFO  | AUMCINF\_obs          | 1313.9510002 | 1313.9510000 |          0 |
|       5 | AUMCIFO  | AUMCINF\_obs          | 1689.4872798 | 1689.4872800 |          0 |
|       6 | AUMCIFO  | AUMCINF\_obs          |  987.9420173 |  987.9420173 |          0 |
|       7 | AUMCIFO  | AUMCINF\_obs          | 1258.3053268 | 1258.3053270 |          0 |
|       8 | AUMCIFO  | AUMCINF\_obs          | 1314.9431383 | 1314.9431380 |          0 |
|       9 | AUMCIFO  | AUMCINF\_obs          | 1219.9213281 | 1219.9213280 |          0 |
|      10 | AUMCIFO  | AUMCINF\_obs          | 2502.5540002 | 2502.5540000 |          0 |
|      11 | AUMCIFO  | AUMCINF\_obs          |  937.9535438 |  937.9535438 |          0 |
|      12 | AUMCIFO  | AUMCINF\_obs          | 1335.1375811 | 1335.1375810 |          0 |
|       1 | AUMCPEO  | AUMC\_.Extrap\_obs    |   67.0201632 |   67.0201632 |          0 |
|       2 | AUMCPEO  | AUMC\_.Extrap\_obs    |   29.0436897 |   29.0436897 |          0 |
|       3 | AUMCPEO  | AUMC\_.Extrap\_obs    |   30.0158308 |   30.0158308 |          0 |
|       4 | AUMCPEO  | AUMC\_.Extrap\_obs    |   30.6075486 |   30.6075486 |          0 |
|       5 | AUMCPEO  | AUMC\_.Extrap\_obs    |   38.5091562 |   38.5091562 |          0 |
|       6 | AUMCPEO  | AUMC\_.Extrap\_obs    |   37.3783169 |   37.3783169 |          0 |
|       7 | AUMCPEO  | AUMC\_.Extrap\_obs    |   36.7699745 |   36.7699745 |          0 |
|       8 | AUMCPEO  | AUMC\_.Extrap\_obs    |   42.4794913 |   42.4794914 |          0 |
|       9 | AUMCPEO  | AUMC\_.Extrap\_obs    |   40.7027815 |   40.7027815 |          0 |
|      10 | AUMCPEO  | AUMC\_.Extrap\_obs    |   47.7837196 |   47.7837196 |          0 |
|      11 | AUMCPEO  | AUMC\_.Extrap\_obs    |   33.1911704 |   33.1911704 |          0 |
|      12 | AUMCPEO  | AUMC\_.Extrap\_obs    |   26.4020191 |   26.4020191 |          0 |
|       1 | AUMCIFP  | AUMCINF\_pred         | 4545.7288462 | 4545.7288460 |          0 |
|       2 | AUMCIFP  | AUMCINF\_pred         | 1005.7637454 | 1005.7637450 |          0 |
|       3 | AUMCIFP  | AUMCINF\_pred         | 1160.3397033 | 1160.3397030 |          0 |
|       4 | AUMCIFP  | AUMCINF\_pred         | 1316.1967080 | 1316.1967080 |          0 |
|       5 | AUMCIFP  | AUMCINF\_pred         | 1683.5593423 | 1683.5593420 |          0 |
|       6 | AUMCIFP  | AUMCINF\_pred         |  996.4799913 |  996.4799913 |          0 |
|       7 | AUMCIFP  | AUMCINF\_pred         | 1262.6179786 | 1262.6179790 |          0 |
|       8 | AUMCIFP  | AUMCINF\_pred         | 1305.3474998 | 1305.3475000 |          0 |
|       9 | AUMCIFP  | AUMCINF\_pred         | 1218.3621498 | 1218.3621500 |          0 |
|      10 | AUMCIFP  | AUMCINF\_pred         | 2499.4371146 | 2499.4371150 |          0 |
|      11 | AUMCIFP  | AUMCINF\_pred         |  937.8835360 |  937.8835360 |          0 |
|      12 | AUMCIFP  | AUMCINF\_pred         | 1336.8064129 | 1336.8064130 |          0 |
|       1 | AUMCPEP  | AUMC\_.Extrap\_pred   |   67.0211503 |   67.0211503 |          0 |
|       2 | AUMCPEP  | AUMC\_.Extrap\_pred   |   28.7826061 |   28.7826061 |          0 |
|       3 | AUMCPEP  | AUMC\_.Extrap\_pred   |   30.1176474 |   30.1176474 |          0 |
|       4 | AUMCPEP  | AUMC\_.Extrap\_pred   |   30.7259467 |   30.7259467 |          0 |
|       5 | AUMCPEP  | AUMC\_.Extrap\_pred   |   38.2926424 |   38.2926424 |          0 |
|       6 | AUMCPEP  | AUMC\_.Extrap\_pred   |   37.9148679 |   37.9148679 |          0 |
|       7 | AUMCPEP  | AUMC\_.Extrap\_pred   |   36.9859457 |   36.9859457 |          0 |
|       8 | AUMCPEP  | AUMC\_.Extrap\_pred   |   42.0566568 |   42.0566568 |          0 |
|       9 | AUMCPEP  | AUMC\_.Extrap\_pred   |   40.6268969 |   40.6268969 |          0 |
|      10 | AUMCPEP  | AUMC\_.Extrap\_pred   |   47.7186040 |   47.7186040 |          0 |
|      11 | AUMCPEP  | AUMC\_.Extrap\_pred   |   33.1861835 |   33.1861835 |          0 |
|      12 | AUMCPEP  | AUMC\_.Extrap\_pred   |   26.4938967 |   26.4938967 |          0 |
|       1 | MRTEVLST | MRTlast               |   10.1818973 |   10.1818973 |          0 |
|       2 | MRTEVLST | MRTlast               |    8.0724494 |    8.0724494 |          0 |
|       3 | MRTEVLST | MRTlast               |    8.4573209 |    8.4573209 |          0 |
|       4 | MRTEVLST | MRTlast               |    8.8838607 |    8.8838607 |          0 |
|       5 | MRTEVLST | MRTlast               |    8.7907063 |    8.7907063 |          0 |
|       6 | MRTEVLST | MRTlast               |    8.6288937 |    8.6288937 |          0 |
|       7 | MRTEVLST | MRTlast               |    9.0443761 |    9.0443761 |          0 |
|       8 | MRTEVLST | MRTlast               |    8.7131889 |    8.7131889 |          0 |
|       9 | MRTEVLST | MRTlast               |    8.6180785 |    8.6180785 |          0 |
|      10 | MRTEVLST | MRTlast               |    9.6384311 |    9.6384311 |          0 |
|      11 | MRTEVLST | MRTlast               |    8.0447792 |    8.0447792 |          0 |
|      12 | MRTEVLST | MRTlast               |    8.5283156 |    8.5283156 |          0 |
|       1 | MRTEVIFO | MRTINF\_obs           |   21.1498046 |   21.1498045 |          0 |
|       2 | MRTEVIFO | MRTINF\_obs           |   10.3664599 |   10.3664599 |          0 |
|       3 | MRTEVIFO | MRTINF\_obs           |   10.9175260 |   10.9175260 |          0 |
|       4 | MRTEVIFO | MRTINF\_obs           |   11.5040681 |   11.5040681 |          0 |
|       5 | MRTEVIFO | MRTINF\_obs           |   12.3949276 |   12.3949276 |          0 |
|       6 | MRTEVIFO | MRTINF\_obs           |   12.0222866 |   12.0222866 |          0 |
|       7 | MRTEVIFO | MRTINF\_obs           |   12.4599947 |   12.4599947 |          0 |
|       8 | MRTEVIFO | MRTINF\_obs           |   12.8722531 |   12.8722531 |          0 |
|       9 | MRTEVIFO | MRTINF\_obs           |   12.5094471 |   12.5094471 |          0 |
|      10 | MRTEVIFO | MRTINF\_obs           |   14.9085758 |   14.9085758 |          0 |
|      11 | MRTEVIFO | MRTINF\_obs           |   10.7931564 |   10.7931564 |          0 |
|      12 | MRTEVIFO | MRTINF\_obs           |   10.6105161 |   10.6105161 |          0 |
|       1 | MRTEVIFP | MRTINF\_pred          |   21.1501401 |   21.1501401 |          0 |
|       2 | MRTEVIFP | MRTINF\_pred          |   10.3400455 |   10.3400455 |          0 |
|       3 | MRTEVIFP | MRTINF\_pred          |   10.9283095 |   10.9283095 |          0 |
|       4 | MRTEVIFP | MRTINF\_pred          |   11.5172082 |   11.5172082 |          0 |
|       5 | MRTEVIFP | MRTINF\_pred          |   12.3664205 |   12.3664205 |          0 |
|       6 | MRTEVIFP | MRTINF\_pred          |   12.0905386 |   12.0905386 |          0 |
|       7 | MRTEVIFP | MRTINF\_pred          |   12.4876944 |   12.4876944 |          0 |
|       8 | MRTEVIFP | MRTINF\_pred          |   12.8113828 |   12.8113828 |          0 |
|       9 | MRTEVIFP | MRTINF\_pred          |   12.4989252 |   12.4989252 |          0 |
|      10 | MRTEVIFP | MRTINF\_pred          |   14.8974756 |   14.8974756 |          0 |
|      11 | MRTEVIFP | MRTINF\_pred          |   10.7926025 |   10.7926025 |          0 |
|      12 | MRTEVIFP | MRTINF\_pred          |   10.6195388 |   10.6195388 |          0 |

Theoph (n=12), Log, Extravascular

## Test 3: Indometh (n=6), Linear, IV Bolus

``` r
table_wres_rres(Wres3, Rres3,
                Caption = 'Indometh (n=6), Linear, IV Bolus')
```

| Subject | PPTESTCD | WNL                   |  NonCompart |   WinNonlin | Difference |
| ------: | :------- | :-------------------- | ----------: | ----------: | ---------: |
|       1 | R2       | Rsq                   |   0.9970667 |   0.9970667 |          0 |
|       2 | R2       | Rsq                   |   0.9476691 |   0.9476691 |          0 |
|       3 | R2       | Rsq                   |   0.8758261 |   0.8758261 |          0 |
|       4 | R2       | Rsq                   |   0.8728249 |   0.8728249 |          0 |
|       5 | R2       | Rsq                   |   0.8752442 |   0.8752442 |          0 |
|       6 | R2       | Rsq                   |   0.9039538 |   0.9039538 |          0 |
|       1 | R2ADJ    | Rsq\_adjusted         |   0.9941335 |   0.9941335 |          0 |
|       2 | R2ADJ    | Rsq\_adjusted         |   0.9401933 |   0.9401933 |          0 |
|       3 | R2ADJ    | Rsq\_adjusted         |   0.8603043 |   0.8603043 |          0 |
|       4 | R2ADJ    | Rsq\_adjusted         |   0.8586943 |   0.8586943 |          0 |
|       5 | R2ADJ    | Rsq\_adjusted         |   0.8544516 |   0.8544516 |          0 |
|       6 | R2ADJ    | Rsq\_adjusted         |   0.8902329 |   0.8902329 |          0 |
|       1 | CORRXY   | Corr\_XY              | \-0.9985323 | \-0.9985323 |          0 |
|       2 | CORRXY   | Corr\_XY              | \-0.9734830 | \-0.9734830 |          0 |
|       3 | CORRXY   | Corr\_XY              | \-0.9358558 | \-0.9358558 |          0 |
|       4 | CORRXY   | Corr\_XY              | \-0.9342510 | \-0.9342510 |          0 |
|       5 | CORRXY   | Corr\_XY              | \-0.9355449 | \-0.9355449 |          0 |
|       6 | CORRXY   | Corr\_XY              | \-0.9507649 | \-0.9507649 |          0 |
|       1 | LAMZNPT  | No\_points\_lambda\_z |   3.0000000 |   3.0000000 |          0 |
|       2 | LAMZNPT  | No\_points\_lambda\_z |   9.0000000 |   9.0000000 |          0 |
|       3 | LAMZNPT  | No\_points\_lambda\_z |  10.0000000 |  10.0000000 |          0 |
|       4 | LAMZNPT  | No\_points\_lambda\_z |  11.0000000 |  11.0000000 |          0 |
|       5 | LAMZNPT  | No\_points\_lambda\_z |   8.0000000 |   8.0000000 |          0 |
|       6 | LAMZNPT  | No\_points\_lambda\_z |   9.0000000 |   9.0000000 |          0 |
|       1 | LAMZ     | Lambda\_z             |   0.1583205 |   0.1583205 |          0 |
|       2 | LAMZ     | Lambda\_z             |   0.3022800 |   0.3022800 |          0 |
|       3 | LAMZ     | Lambda\_z             |   0.4218926 |   0.4218926 |          0 |
|       4 | LAMZ     | Lambda\_z             |   0.4554455 |   0.4554455 |          0 |
|       5 | LAMZ     | Lambda\_z             |   0.2527478 |   0.2527478 |          0 |
|       6 | LAMZ     | Lambda\_z             |   0.3535205 |   0.3535205 |          0 |
|       1 | LAMZLL   | Lambda\_z\_lower      |   5.0000000 |   5.0000000 |          0 |
|       2 | LAMZLL   | Lambda\_z\_lower      |   0.7500000 |   0.7500000 |          0 |
|       3 | LAMZLL   | Lambda\_z\_lower      |   0.5000000 |   0.5000000 |          0 |
|       4 | LAMZLL   | Lambda\_z\_lower      |   0.2500000 |   0.2500000 |          0 |
|       5 | LAMZLL   | Lambda\_z\_lower      |   1.0000000 |   1.0000000 |          0 |
|       6 | LAMZLL   | Lambda\_z\_lower      |   0.7500000 |   0.7500000 |          0 |
|       1 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       2 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       3 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       4 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       5 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       6 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       1 | LAMZHL   | HL\_Lambda\_z         |   4.3781270 |   4.3781270 |          0 |
|       2 | LAMZHL   | HL\_Lambda\_z         |   2.2930632 |   2.2930632 |          0 |
|       3 | LAMZHL   | HL\_Lambda\_z         |   1.6429468 |   1.6429468 |          0 |
|       4 | LAMZHL   | HL\_Lambda\_z         |   1.5219104 |   1.5219104 |          0 |
|       5 | LAMZHL   | HL\_Lambda\_z         |   2.7424461 |   2.7424461 |          0 |
|       6 | LAMZHL   | HL\_Lambda\_z         |   1.9606986 |   1.9606986 |          0 |
|       1 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       2 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       3 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       4 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       5 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       6 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       1 | CMAX     | Cmax                  |   1.5000000 |   1.5000000 |          0 |
|       2 | CMAX     | Cmax                  |   2.0300000 |   2.0300000 |          0 |
|       3 | CMAX     | Cmax                  |   2.7200000 |   2.7200000 |          0 |
|       4 | CMAX     | Cmax                  |   1.8500000 |   1.8500000 |          0 |
|       5 | CMAX     | Cmax                  |   2.0500000 |   2.0500000 |          0 |
|       6 | CMAX     | Cmax                  |   2.3100000 |   2.3100000 |          0 |
|       1 | CMAXD    | Cmax\_D               |   0.0600000 |   0.0600000 |          0 |
|       2 | CMAXD    | Cmax\_D               |   0.0812000 |   0.0812000 |          0 |
|       3 | CMAXD    | Cmax\_D               |   0.1088000 |   0.1088000 |          0 |
|       4 | CMAXD    | Cmax\_D               |   0.0740000 |   0.0740000 |          0 |
|       5 | CMAXD    | Cmax\_D               |   0.0820000 |   0.0820000 |          0 |
|       6 | CMAXD    | Cmax\_D               |   0.0924000 |   0.0924000 |          0 |
|       1 | C0       | C0                    |   2.3936170 |   2.3936170 |          0 |
|       2 | C0       | C0                    |   2.5281595 |   2.5281595 |          0 |
|       3 | C0       | C0                    |   4.9653691 |   4.9653691 |          0 |
|       4 | C0       | C0                    |   2.4622302 |   2.4622302 |          0 |
|       5 | C0       | C0                    |   4.0408654 |   4.0408654 |          0 |
|       6 | C0       | C0                    |   3.7056250 |   3.7056250 |          0 |
|       1 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       2 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       3 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       4 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       5 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       6 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       1 | CLST     | Clast                 |   0.0500000 |   0.0500000 |          0 |
|       2 | CLST     | Clast                 |   0.0800000 |   0.0800000 |          0 |
|       3 | CLST     | Clast                 |   0.0800000 |   0.0800000 |          0 |
|       4 | CLST     | Clast                 |   0.0700000 |   0.0700000 |          0 |
|       5 | CLST     | Clast                 |   0.0600000 |   0.0600000 |          0 |
|       6 | CLST     | Clast                 |   0.0900000 |   0.0900000 |          0 |
|       1 | AUCLST   | AUClast               |   2.0404521 |   2.0404521 |          0 |
|       2 | AUCLST   | AUClast               |   3.2485199 |   3.2485199 |          0 |
|       3 | AUCLST   | AUClast               |   3.5544211 |   3.5544211 |          0 |
|       4 | AUCLST   | AUClast               |   2.7852788 |   2.7852788 |          0 |
|       5 | AUCLST   | AUClast               |   2.4588582 |   2.4588582 |          0 |
|       6 | AUCLST   | AUClast               |   3.3357031 |   3.3357031 |          0 |
|       1 | AUCALL   | AUCall                |   2.0404521 |   2.0404521 |          0 |
|       2 | AUCALL   | AUCall                |   3.2485199 |   3.2485199 |          0 |
|       3 | AUCALL   | AUCall                |   3.5544211 |   3.5544211 |          0 |
|       4 | AUCALL   | AUCall                |   2.7852788 |   2.7852788 |          0 |
|       5 | AUCALL   | AUCall                |   2.4588582 |   2.4588582 |          0 |
|       6 | AUCALL   | AUCall                |   3.3357031 |   3.3357031 |          0 |
|       1 | AUCIFO   | AUCINF\_obs           |   2.3562672 |   2.3562672 |          0 |
|       2 | AUCIFO   | AUCINF\_obs           |   3.5131752 |   3.5131752 |          0 |
|       3 | AUCIFO   | AUCINF\_obs           |   3.7440428 |   3.7440428 |          0 |
|       4 | AUCIFO   | AUCINF\_obs           |   2.9389745 |   2.9389745 |          0 |
|       5 | AUCIFO   | AUCINF\_obs           |   2.6962490 |   2.6962490 |          0 |
|       6 | AUCIFO   | AUCINF\_obs           |   3.5902852 |   3.5902852 |          0 |
|       1 | AUCIFOD  | AUCINF\_D\_obs        |   0.0942507 |   0.0942507 |          0 |
|       2 | AUCIFOD  | AUCINF\_D\_obs        |   0.1405270 |   0.1405270 |          0 |
|       3 | AUCIFOD  | AUCINF\_D\_obs        |   0.1497617 |   0.1497617 |          0 |
|       4 | AUCIFOD  | AUCINF\_D\_obs        |   0.1175590 |   0.1175590 |          0 |
|       5 | AUCIFOD  | AUCINF\_D\_obs        |   0.1078500 |   0.1078500 |          0 |
|       6 | AUCIFOD  | AUCINF\_D\_obs        |   0.1436114 |   0.1436114 |          0 |
|       1 | AUCPEO   | AUC\_.Extrap\_obs     |  13.4031956 |  13.4031956 |          0 |
|       2 | AUCPEO   | AUC\_.Extrap\_obs     |   7.5332215 |   7.5332215 |          0 |
|       3 | AUCPEO   | AUC\_.Extrap\_obs     |   5.0646241 |   5.0646241 |          0 |
|       4 | AUCPEO   | AUC\_.Extrap\_obs     |   5.2295685 |   5.2295685 |          0 |
|       5 | AUCPEO   | AUC\_.Extrap\_obs     |   8.8044838 |   8.8044838 |          0 |
|       6 | AUCPEO   | AUC\_.Extrap\_obs     |   7.0908603 |   7.0908603 |          0 |
|       1 | AUCPBEO  | AUC\_.Back\_Ext\_obs  |  20.6556421 |  20.6556421 |          0 |
|       2 | AUCPBEO  | AUC\_.Back\_Ext\_obs  |  16.2180906 |  16.2180906 |          0 |
|       3 | AUCPBEO  | AUC\_.Back\_Ext\_obs  |  25.6586578 |  25.6586578 |          0 |
|       4 | AUCPBEO  | AUC\_.Back\_Ext\_obs  |  18.3407098 |  18.3407098 |          0 |
|       5 | AUCPBEO  | AUC\_.Back\_Ext\_obs  |  28.2376805 |  28.2376805 |          0 |
|       6 | AUCPBEO  | AUC\_.Back\_Ext\_obs  |  20.9441054 |  20.9441054 |          0 |
|       1 | VZO      | Vz\_obs               |  67.0159780 |  67.0159780 |          0 |
|       2 | VZO      | Vz\_obs               |  23.5413171 |  23.5413171 |          0 |
|       3 | VZO      | Vz\_obs               |  15.8269504 |  15.8269504 |          0 |
|       4 | VZO      | Vz\_obs               |  18.6770303 |  18.6770303 |          0 |
|       5 | VZO      | Vz\_obs               |  36.6853493 |  36.6853493 |          0 |
|       6 | VZO      | Vz\_obs               |  19.6968341 |  19.6968341 |          0 |
|       1 | CLO      | Cl\_obs               |  10.6100020 |  10.6100020 |          0 |
|       2 | CLO      | Cl\_obs               |   7.1160698 |   7.1160698 |          0 |
|       3 | CLO      | Cl\_obs               |   6.6772740 |   6.6772740 |          0 |
|       4 | CLO      | Cl\_obs               |   8.5063686 |   8.5063686 |          0 |
|       5 | CLO      | Cl\_obs               |   9.2721407 |   9.2721407 |          0 |
|       6 | CLO      | Cl\_obs               |   6.9632351 |   6.9632351 |          0 |
|       1 | AUCIFP   | AUCINF\_pred          |   2.3578369 |   2.3578369 |          0 |
|       2 | AUCIFP   | AUCINF\_pred          |   3.4958268 |   3.4958268 |          0 |
|       3 | AUCIFP   | AUCINF\_pred          |   3.6491670 |   3.6491670 |          0 |
|       4 | AUCIFP   | AUCINF\_pred          |   2.8554521 |   2.8554521 |          0 |
|       5 | AUCIFP   | AUCINF\_pred          |   2.6549884 |   2.6549884 |          0 |
|       6 | AUCIFP   | AUCINF\_pred          |   3.4947956 |   3.4947956 |          0 |
|       1 | AUCIFPD  | AUCINF\_D\_pred       |   0.0943135 |   0.0943135 |          0 |
|       2 | AUCIFPD  | AUCINF\_D\_pred       |   0.1398331 |   0.1398331 |          0 |
|       3 | AUCIFPD  | AUCINF\_D\_pred       |   0.1459667 |   0.1459667 |          0 |
|       4 | AUCIFPD  | AUCINF\_D\_pred       |   0.1142181 |   0.1142181 |          0 |
|       5 | AUCIFPD  | AUCINF\_D\_pred       |   0.1061995 |   0.1061995 |          0 |
|       6 | AUCIFPD  | AUCINF\_D\_pred       |   0.1397918 |   0.1397918 |          0 |
|       1 | AUCPEP   | AUC\_.Extrap\_pred    |  13.4608442 |  13.4608442 |          0 |
|       2 | AUCPEP   | AUC\_.Extrap\_pred    |   7.0743442 |   7.0743442 |          0 |
|       3 | AUCPEP   | AUC\_.Extrap\_pred    |   2.5963692 |   2.5963692 |          0 |
|       4 | AUCPEP   | AUC\_.Extrap\_pred    |   2.4575198 |   2.4575198 |          0 |
|       5 | AUCPEP   | AUC\_.Extrap\_pred    |   7.3872362 |   7.3872362 |          0 |
|       6 | AUCPEP   | AUC\_.Extrap\_pred    |   4.5522694 |   4.5522694 |          0 |
|       1 | AUCPBEP  | AUC\_.Back\_Ext\_pred |  20.6418914 |  20.6418914 |          0 |
|       2 | AUCPBEP  | AUC\_.Back\_Ext\_pred |  16.2985748 |  16.2985748 |          0 |
|       3 | AUCPBEP  | AUC\_.Back\_Ext\_pred |  26.3257654 |  26.3257654 |          0 |
|       4 | AUCPBEP  | AUC\_.Back\_Ext\_pred |  18.8771782 |  18.8771782 |          0 |
|       5 | AUCPBEP  | AUC\_.Back\_Ext\_pred |  28.6765156 |  28.6765156 |          0 |
|       6 | AUCPBEP  | AUC\_.Back\_Ext\_pred |  21.5163690 |  21.5163690 |          0 |
|       1 | VZP      | Vz\_pred              |  66.9713647 |  66.9713647 |          0 |
|       2 | VZP      | Vz\_pred              |  23.6581437 |  23.6581437 |          0 |
|       3 | VZP      | Vz\_pred              |  16.2384403 |  16.2384403 |          0 |
|       4 | VZP      | Vz\_pred              |  19.2233361 |  19.2233361 |          0 |
|       5 | VZP      | Vz\_pred              |  37.2554675 |  37.2554675 |          0 |
|       6 | VZP      | Vz\_pred              |  20.2350180 |  20.2350180 |          0 |
|       1 | CLP      | Cl\_pred              |  10.6029388 |  10.6029388 |          0 |
|       2 | CLP      | Cl\_pred              |   7.1513841 |   7.1513841 |          0 |
|       3 | CLP      | Cl\_pred              |   6.8508786 |   6.8508786 |          0 |
|       4 | CLP      | Cl\_pred              |   8.7551811 |   8.7551811 |          0 |
|       5 | CLP      | Cl\_pred              |   9.4162369 |   9.4162369 |          0 |
|       6 | CLP      | Cl\_pred              |   7.1534941 |   7.1534941 |          0 |
|       1 | AUMCLST  | AUMClast              |   3.2712500 |   3.2712500 |          0 |
|       2 | AUMCLST  | AUMClast              |   6.3987500 |   6.3987500 |          0 |
|       3 | AUMCLST  | AUMClast              |   5.0062500 |   5.0062500 |          0 |
|       4 | AUMCLST  | AUMClast              |   4.3818750 |   4.3818750 |          0 |
|       5 | AUMCLST  | AUMClast              |   3.7075000 |   3.7075000 |          0 |
|       6 | AUMCLST  | AUMClast              |   5.5325000 |   5.5325000 |          0 |
|       1 | AUMCIFO  | AUMCINF\_obs          |   7.7925545 |   7.7925545 |          0 |
|       2 | AUMCIFO  | AUMCINF\_obs          |   9.3915223 |   9.3915223 |          0 |
|       3 | AUMCIFO  | AUMCINF\_obs          |   6.9726784 |   6.9726784 |          0 |
|       4 | AUMCIFO  | AUMCINF\_obs          |   5.9489028 |   5.9489028 |          0 |
|       5 | AUMCIFO  | AUMCINF\_obs          |   6.5458663 |   6.5458663 |          0 |
|       6 | AUMCIFO  | AUMCINF\_obs          |   8.2892908 |   8.2892908 |          0 |
|       1 | AUMCPEO  | AUMC\_.Extrap\_obs    |  58.0208261 |  58.0208261 |          0 |
|       2 | AUMCPEO  | AUMC\_.Extrap\_obs    |  31.8667432 |  31.8667432 |          0 |
|       3 | AUMCPEO  | AUMC\_.Extrap\_obs    |  28.2019090 |  28.2019090 |          0 |
|       4 | AUMCPEO  | AUMC\_.Extrap\_obs    |  26.3414589 |  26.3414589 |          0 |
|       5 | AUMCPEO  | AUMC\_.Extrap\_obs    |  43.3612023 |  43.3612023 |          0 |
|       6 | AUMCPEO  | AUMC\_.Extrap\_obs    |  33.2572574 |  33.2572574 |          0 |
|       1 | AUMCIFP  | AUMCINF\_pred         |   7.8150259 |   7.8150259 |          0 |
|       2 | AUMCIFP  | AUMCINF\_pred         |   9.1953427 |   9.1953427 |          0 |
|       3 | AUMCIFP  | AUMCINF\_pred         |   5.9887901 |   5.9887901 |          0 |
|       4 | AUMCIFP  | AUMCINF\_pred         |   5.0973376 |   5.0973376 |          0 |
|       5 | AUMCIFP  | AUMCINF\_pred         |   6.0525342 |   6.0525342 |          0 |
|       6 | AUMCIFP  | AUMCINF\_pred         |   7.2552635 |   7.2552635 |          0 |
|       1 | AUMCPEP  | AUMC\_.Extrap\_pred   |  58.1415337 |  58.1415337 |          0 |
|       2 | AUMCPEP  | AUMC\_.Extrap\_pred   |  30.4131425 |  30.4131425 |          0 |
|       3 | AUMCPEP  | AUMC\_.Extrap\_pred   |  16.4063210 |  16.4063210 |          0 |
|       4 | AUMCPEP  | AUMC\_.Extrap\_pred   |  14.0360055 |  14.0360055 |          0 |
|       5 | AUMCPEP  | AUMC\_.Extrap\_pred   |  38.7446663 |  38.7446663 |          0 |
|       6 | AUMCPEP  | AUMC\_.Extrap\_pred   |  23.7450164 |  23.7450164 |          0 |
|       1 | MRTIVLST | MRTlast               |   1.6031986 |   1.6031986 |          0 |
|       2 | MRTIVLST | MRTlast               |   1.9697432 |   1.9697432 |          0 |
|       3 | MRTIVLST | MRTlast               |   1.4084572 |   1.4084572 |          0 |
|       4 | MRTIVLST | MRTlast               |   1.5732267 |   1.5732267 |          0 |
|       5 | MRTIVLST | MRTlast               |   1.5078137 |   1.5078137 |          0 |
|       6 | MRTIVLST | MRTlast               |   1.6585709 |   1.6585709 |          0 |
|       1 | MRTIVIFO | MRTINF\_obs           |   3.3071607 |   3.3071607 |          0 |
|       2 | MRTIVIFO | MRTINF\_obs           |   2.6732291 |   2.6732291 |          0 |
|       3 | MRTIVIFO | MRTINF\_obs           |   1.8623394 |   1.8623394 |          0 |
|       4 | MRTIVIFO | MRTINF\_obs           |   2.0241424 |   2.0241424 |          0 |
|       5 | MRTIVIFO | MRTINF\_obs           |   2.4277678 |   2.4277678 |          0 |
|       6 | MRTIVIFO | MRTINF\_obs           |   2.3088112 |   2.3088112 |          0 |
|       1 | MRTIVIFP | MRTINF\_pred          |   3.3144897 |   3.3144897 |          0 |
|       2 | MRTIVIFP | MRTINF\_pred          |   2.6303771 |   2.6303771 |          0 |
|       3 | MRTIVIFP | MRTINF\_pred          |   1.6411390 |   1.6411390 |          0 |
|       4 | MRTIVIFP | MRTINF\_pred          |   1.7851245 |   1.7851245 |          0 |
|       5 | MRTIVIFP | MRTINF\_pred          |   2.2796838 |   2.2796838 |          0 |
|       6 | MRTIVIFP | MRTINF\_pred          |   2.0760194 |   2.0760194 |          0 |
|       1 | VSSO     | Vss\_obs              |  35.0889819 |  35.0889819 |          0 |
|       2 | VSSO     | Vss\_obs              |  19.0228851 |  19.0228851 |          0 |
|       3 | VSSO     | Vss\_obs              |  12.4353504 |  12.4353504 |          0 |
|       4 | VSSO     | Vss\_obs              |  17.2181012 |  17.2181012 |          0 |
|       5 | VSSO     | Vss\_obs              |  22.5106044 |  22.5106044 |          0 |
|       6 | VSSO     | Vss\_obs              |  16.0767951 |  16.0767951 |          0 |
|       1 | VSSP     | Vss\_pred             |  35.1433309 |  35.1433309 |          0 |
|       2 | VSSP     | Vss\_pred             |  18.8108371 |  18.8108371 |          0 |
|       3 | VSSP     | Vss\_pred             |  11.2432438 |  11.2432438 |          0 |
|       4 | VSSP     | Vss\_pred             |  15.6290886 |  15.6290886 |          0 |
|       5 | VSSP     | Vss\_pred             |  21.4660427 |  21.4660427 |          0 |
|       6 | VSSP     | Vss\_pred             |  14.8507925 |  14.8507925 |          0 |

Indometh (n=6), Linear, IV Bolus

## Test 4: Indometh (n=6), Log, IV Bolus

``` r
table_wres_rres(Wres4, Rres4,
                Caption = 'Indometh (n=6), Log, IV Bolus')
```

| Subject | PPTESTCD | WNL                   |  NonCompart |   WinNonlin | Difference |
| ------: | :------- | :-------------------- | ----------: | ----------: | ---------: |
|       1 | R2       | Rsq                   |   0.9970667 |   0.9970667 |          0 |
|       2 | R2       | Rsq                   |   0.9476691 |   0.9476691 |          0 |
|       3 | R2       | Rsq                   |   0.8758261 |   0.8758261 |          0 |
|       4 | R2       | Rsq                   |   0.8728249 |   0.8728249 |          0 |
|       5 | R2       | Rsq                   |   0.8752442 |   0.8752442 |          0 |
|       6 | R2       | Rsq                   |   0.9039538 |   0.9039538 |          0 |
|       1 | R2ADJ    | Rsq\_adjusted         |   0.9941335 |   0.9941335 |          0 |
|       2 | R2ADJ    | Rsq\_adjusted         |   0.9401933 |   0.9401933 |          0 |
|       3 | R2ADJ    | Rsq\_adjusted         |   0.8603043 |   0.8603043 |          0 |
|       4 | R2ADJ    | Rsq\_adjusted         |   0.8586943 |   0.8586943 |          0 |
|       5 | R2ADJ    | Rsq\_adjusted         |   0.8544516 |   0.8544516 |          0 |
|       6 | R2ADJ    | Rsq\_adjusted         |   0.8902329 |   0.8902329 |          0 |
|       1 | CORRXY   | Corr\_XY              | \-0.9985323 | \-0.9985323 |          0 |
|       2 | CORRXY   | Corr\_XY              | \-0.9734830 | \-0.9734830 |          0 |
|       3 | CORRXY   | Corr\_XY              | \-0.9358558 | \-0.9358558 |          0 |
|       4 | CORRXY   | Corr\_XY              | \-0.9342510 | \-0.9342510 |          0 |
|       5 | CORRXY   | Corr\_XY              | \-0.9355449 | \-0.9355449 |          0 |
|       6 | CORRXY   | Corr\_XY              | \-0.9507649 | \-0.9507649 |          0 |
|       1 | LAMZNPT  | No\_points\_lambda\_z |   3.0000000 |   3.0000000 |          0 |
|       2 | LAMZNPT  | No\_points\_lambda\_z |   9.0000000 |   9.0000000 |          0 |
|       3 | LAMZNPT  | No\_points\_lambda\_z |  10.0000000 |  10.0000000 |          0 |
|       4 | LAMZNPT  | No\_points\_lambda\_z |  11.0000000 |  11.0000000 |          0 |
|       5 | LAMZNPT  | No\_points\_lambda\_z |   8.0000000 |   8.0000000 |          0 |
|       6 | LAMZNPT  | No\_points\_lambda\_z |   9.0000000 |   9.0000000 |          0 |
|       1 | LAMZ     | Lambda\_z             |   0.1583205 |   0.1583205 |          0 |
|       2 | LAMZ     | Lambda\_z             |   0.3022800 |   0.3022800 |          0 |
|       3 | LAMZ     | Lambda\_z             |   0.4218926 |   0.4218926 |          0 |
|       4 | LAMZ     | Lambda\_z             |   0.4554455 |   0.4554455 |          0 |
|       5 | LAMZ     | Lambda\_z             |   0.2527478 |   0.2527478 |          0 |
|       6 | LAMZ     | Lambda\_z             |   0.3535205 |   0.3535205 |          0 |
|       1 | LAMZLL   | Lambda\_z\_lower      |   5.0000000 |   5.0000000 |          0 |
|       2 | LAMZLL   | Lambda\_z\_lower      |   0.7500000 |   0.7500000 |          0 |
|       3 | LAMZLL   | Lambda\_z\_lower      |   0.5000000 |   0.5000000 |          0 |
|       4 | LAMZLL   | Lambda\_z\_lower      |   0.2500000 |   0.2500000 |          0 |
|       5 | LAMZLL   | Lambda\_z\_lower      |   1.0000000 |   1.0000000 |          0 |
|       6 | LAMZLL   | Lambda\_z\_lower      |   0.7500000 |   0.7500000 |          0 |
|       1 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       2 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       3 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       4 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       5 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       6 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       1 | LAMZHL   | HL\_Lambda\_z         |   4.3781270 |   4.3781270 |          0 |
|       2 | LAMZHL   | HL\_Lambda\_z         |   2.2930632 |   2.2930632 |          0 |
|       3 | LAMZHL   | HL\_Lambda\_z         |   1.6429468 |   1.6429468 |          0 |
|       4 | LAMZHL   | HL\_Lambda\_z         |   1.5219104 |   1.5219104 |          0 |
|       5 | LAMZHL   | HL\_Lambda\_z         |   2.7424461 |   2.7424461 |          0 |
|       6 | LAMZHL   | HL\_Lambda\_z         |   1.9606986 |   1.9606986 |          0 |
|       1 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       2 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       3 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       4 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       5 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       6 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       1 | CMAX     | Cmax                  |   1.5000000 |   1.5000000 |          0 |
|       2 | CMAX     | Cmax                  |   2.0300000 |   2.0300000 |          0 |
|       3 | CMAX     | Cmax                  |   2.7200000 |   2.7200000 |          0 |
|       4 | CMAX     | Cmax                  |   1.8500000 |   1.8500000 |          0 |
|       5 | CMAX     | Cmax                  |   2.0500000 |   2.0500000 |          0 |
|       6 | CMAX     | Cmax                  |   2.3100000 |   2.3100000 |          0 |
|       1 | CMAXD    | Cmax\_D               |   0.0600000 |   0.0600000 |          0 |
|       2 | CMAXD    | Cmax\_D               |   0.0812000 |   0.0812000 |          0 |
|       3 | CMAXD    | Cmax\_D               |   0.1088000 |   0.1088000 |          0 |
|       4 | CMAXD    | Cmax\_D               |   0.0740000 |   0.0740000 |          0 |
|       5 | CMAXD    | Cmax\_D               |   0.0820000 |   0.0820000 |          0 |
|       6 | CMAXD    | Cmax\_D               |   0.0924000 |   0.0924000 |          0 |
|       1 | C0       | C0                    |   2.3936170 |   2.3936170 |          0 |
|       2 | C0       | C0                    |   2.5281595 |   2.5281595 |          0 |
|       3 | C0       | C0                    |   4.9653691 |   4.9653691 |          0 |
|       4 | C0       | C0                    |   2.4622302 |   2.4622302 |          0 |
|       5 | C0       | C0                    |   4.0408654 |   4.0408654 |          0 |
|       6 | C0       | C0                    |   3.7056250 |   3.7056250 |          0 |
|       1 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       2 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       3 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       4 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       5 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       6 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       1 | CLST     | Clast                 |   0.0500000 |   0.0500000 |          0 |
|       2 | CLST     | Clast                 |   0.0800000 |   0.0800000 |          0 |
|       3 | CLST     | Clast                 |   0.0800000 |   0.0800000 |          0 |
|       4 | CLST     | Clast                 |   0.0700000 |   0.0700000 |          0 |
|       5 | CLST     | Clast                 |   0.0600000 |   0.0600000 |          0 |
|       6 | CLST     | Clast                 |   0.0900000 |   0.0900000 |          0 |
|       1 | AUCLST   | AUClast               |   2.0098984 |   2.0098984 |          0 |
|       2 | AUCLST   | AUClast               |   3.2028878 |   3.2028878 |          0 |
|       3 | AUCLST   | AUClast               |   3.4743971 |   3.4743971 |          0 |
|       4 | AUCLST   | AUClast               |   2.7483832 |   2.7483832 |          0 |
|       5 | AUCLST   | AUClast               |   2.3983736 |   2.3983736 |          0 |
|       6 | AUCLST   | AUClast               |   3.2908266 |   3.2908266 |          0 |
|       1 | AUCALL   | AUCall                |   2.0098984 |   2.0098984 |          0 |
|       2 | AUCALL   | AUCall                |   3.2028878 |   3.2028878 |          0 |
|       3 | AUCALL   | AUCall                |   3.4743971 |   3.4743971 |          0 |
|       4 | AUCALL   | AUCall                |   2.7483832 |   2.7483832 |          0 |
|       5 | AUCALL   | AUCall                |   2.3983736 |   2.3983736 |          0 |
|       6 | AUCALL   | AUCall                |   3.2908266 |   3.2908266 |          0 |
|       1 | AUCIFO   | AUCINF\_obs           |   2.3257135 |   2.3257135 |          0 |
|       2 | AUCIFO   | AUCINF\_obs           |   3.4675431 |   3.4675431 |          0 |
|       3 | AUCIFO   | AUCINF\_obs           |   3.6640188 |   3.6640188 |          0 |
|       4 | AUCIFO   | AUCINF\_obs           |   2.9020789 |   2.9020789 |          0 |
|       5 | AUCIFO   | AUCINF\_obs           |   2.6357645 |   2.6357645 |          0 |
|       6 | AUCIFO   | AUCINF\_obs           |   3.5454087 |   3.5454087 |          0 |
|       1 | AUCIFOD  | AUCINF\_D\_obs        |   0.0930285 |   0.0930285 |          0 |
|       2 | AUCIFOD  | AUCINF\_D\_obs        |   0.1387017 |   0.1387017 |          0 |
|       3 | AUCIFOD  | AUCINF\_D\_obs        |   0.1465608 |   0.1465608 |          0 |
|       4 | AUCIFOD  | AUCINF\_D\_obs        |   0.1160832 |   0.1160832 |          0 |
|       5 | AUCIFOD  | AUCINF\_D\_obs        |   0.1054306 |   0.1054306 |          0 |
|       6 | AUCIFOD  | AUCINF\_D\_obs        |   0.1418163 |   0.1418163 |          0 |
|       1 | AUCPEO   | AUC\_.Extrap\_obs     |  13.5792780 |  13.5792780 |          0 |
|       2 | AUCPEO   | AUC\_.Extrap\_obs     |   7.6323571 |   7.6323571 |          0 |
|       3 | AUCPEO   | AUC\_.Extrap\_obs     |   5.1752381 |   5.1752381 |          0 |
|       4 | AUCPEO   | AUC\_.Extrap\_obs     |   5.2960545 |   5.2960545 |          0 |
|       5 | AUCPEO   | AUC\_.Extrap\_obs     |   9.0065258 |   9.0065258 |          0 |
|       6 | AUCPEO   | AUC\_.Extrap\_obs     |   7.1806138 |   7.1806138 |          0 |
|       1 | AUCPBEO  | AUC\_.Back\_Ext\_obs  |  20.5542573 |  20.5542573 |          0 |
|       2 | AUCPBEO  | AUC\_.Back\_Ext\_obs  |  16.3658871 |  16.3658871 |          0 |
|       3 | AUCPBEO  | AUC\_.Back\_Ext\_obs  |  25.4552663 |  25.4552663 |          0 |
|       4 | AUCPBEO  | AUC\_.Back\_Ext\_obs  |  18.4484084 |  18.4484084 |          0 |
|       5 | AUCPBEO  | AUC\_.Back\_Ext\_obs  |  27.8259014 |  27.8259014 |          0 |
|       6 | AUCPBEO  | AUC\_.Back\_Ext\_obs  |  20.8230657 |  20.8230657 |          0 |
|       1 | VZO      | Vz\_obs               |  67.8963898 |  67.8963898 |          0 |
|       2 | VZO      | Vz\_obs               |  23.8511160 |  23.8511160 |          0 |
|       3 | VZO      | Vz\_obs               |  16.1726192 |  16.1726192 |          0 |
|       4 | VZO      | Vz\_obs               |  18.9144805 |  18.9144805 |          0 |
|       5 | VZO      | Vz\_obs               |  37.5271908 |  37.5271908 |          0 |
|       6 | VZO      | Vz\_obs               |  19.9461495 |  19.9461495 |          0 |
|       1 | CLO      | Cl\_obs               |  10.7493892 |  10.7493892 |          0 |
|       2 | CLO      | Cl\_obs               |   7.2097158 |   7.2097158 |          0 |
|       3 | CLO      | Cl\_obs               |   6.8231092 |   6.8231092 |          0 |
|       4 | CLO      | Cl\_obs               |   8.6145142 |   8.6145142 |          0 |
|       5 | CLO      | Cl\_obs               |   9.4849143 |   9.4849143 |          0 |
|       6 | CLO      | Cl\_obs               |   7.0513732 |   7.0513732 |          0 |
|       1 | AUCIFP   | AUCINF\_pred          |   2.3272832 |   2.3272832 |          0 |
|       2 | AUCIFP   | AUCINF\_pred          |   3.4501946 |   3.4501946 |          0 |
|       3 | AUCIFP   | AUCINF\_pred          |   3.5691429 |   3.5691429 |          0 |
|       4 | AUCIFP   | AUCINF\_pred          |   2.8185565 |   2.8185565 |          0 |
|       5 | AUCIFP   | AUCINF\_pred          |   2.5945039 |   2.5945039 |          0 |
|       6 | AUCIFP   | AUCINF\_pred          |   3.4499191 |   3.4499191 |          0 |
|       1 | AUCIFPD  | AUCINF\_D\_pred       |   0.0930913 |   0.0930913 |          0 |
|       2 | AUCIFPD  | AUCINF\_D\_pred       |   0.1380078 |   0.1380078 |          0 |
|       3 | AUCIFPD  | AUCINF\_D\_pred       |   0.1427657 |   0.1427657 |          0 |
|       4 | AUCIFPD  | AUCINF\_D\_pred       |   0.1127423 |   0.1127423 |          0 |
|       5 | AUCIFPD  | AUCINF\_D\_pred       |   0.1037802 |   0.1037802 |          0 |
|       6 | AUCIFPD  | AUCINF\_D\_pred       |   0.1379968 |   0.1379968 |          0 |
|       1 | AUCPEP   | AUC\_.Extrap\_pred    |  13.6375646 |  13.6375646 |          0 |
|       2 | AUCPEP   | AUC\_.Extrap\_pred    |   7.1679092 |   7.1679092 |          0 |
|       3 | AUCPEP   | AUC\_.Extrap\_pred    |   2.6545826 |   2.6545826 |          0 |
|       4 | AUCPEP   | AUC\_.Extrap\_pred    |   2.4896893 |   2.4896893 |          0 |
|       5 | AUCPEP   | AUC\_.Extrap\_pred    |   7.5594516 |   7.5594516 |          0 |
|       6 | AUCPEP   | AUC\_.Extrap\_pred    |   4.6114853 |   4.6114853 |          0 |
|       1 | AUCPBEP  | AUC\_.Back\_Ext\_pred |  20.5403945 |  20.5403945 |          0 |
|       2 | AUCPBEP  | AUC\_.Back\_Ext\_pred |  16.4481790 |  16.4481790 |          0 |
|       3 | AUCPBEP  | AUC\_.Back\_Ext\_pred |  26.1319245 |  26.1319245 |          0 |
|       4 | AUCPBEP  | AUC\_.Back\_Ext\_pred |  18.9950907 |  18.9950907 |          0 |
|       5 | AUCPBEP  | AUC\_.Back\_Ext\_pred |  28.2684182 |  28.2684182 |          0 |
|       6 | AUCPBEP  | AUC\_.Back\_Ext\_pred |  21.3994230 |  21.3994230 |          0 |
|       1 | VZP      | Vz\_pred              |  67.8505969 |  67.8505969 |          0 |
|       2 | VZP      | Vz\_pred              |  23.9710455 |  23.9710455 |          0 |
|       3 | VZP      | Vz\_pred              |  16.6025238 |  16.6025238 |          0 |
|       4 | VZP      | Vz\_pred              |  19.4749739 |  19.4749739 |          0 |
|       5 | VZP      | Vz\_pred              |  38.1239878 |  38.1239878 |          0 |
|       6 | VZP      | Vz\_pred              |  20.4982349 |  20.4982349 |          0 |
|       1 | CLP      | Cl\_pred              |  10.7421392 |  10.7421392 |          0 |
|       2 | CLP      | Cl\_pred              |   7.2459681 |   7.2459681 |          0 |
|       3 | CLP      | Cl\_pred              |   7.0044827 |   7.0044827 |          0 |
|       4 | CLP      | Cl\_pred              |   8.8697884 |   8.8697884 |          0 |
|       5 | CLP      | Cl\_pred              |   9.6357534 |   9.6357534 |          0 |
|       6 | CLP      | Cl\_pred              |   7.2465467 |   7.2465467 |          0 |
|       1 | AUMCLST  | AUMClast              |   3.3047961 |   3.3047961 |          0 |
|       2 | AUMCLST  | AUMClast              |   6.4131687 |   6.4131687 |          0 |
|       3 | AUMCLST  | AUMClast              |   5.0552993 |   5.0552993 |          0 |
|       4 | AUMCLST  | AUMClast              |   4.4049718 |   4.4049718 |          0 |
|       5 | AUMCLST  | AUMClast              |   3.7472994 |   3.7472994 |          0 |
|       6 | AUMCLST  | AUMClast              |   5.5904206 |   5.5904206 |          0 |
|       1 | AUMCIFO  | AUMCINF\_obs          |   7.8261005 |   7.8261005 |          0 |
|       2 | AUMCIFO  | AUMCINF\_obs          |   9.4059410 |   9.4059410 |          0 |
|       3 | AUMCIFO  | AUMCINF\_obs          |   7.0217278 |   7.0217278 |          0 |
|       4 | AUMCIFO  | AUMCINF\_obs          |   5.9719996 |   5.9719996 |          0 |
|       5 | AUMCIFO  | AUMCINF\_obs          |   6.5856658 |   6.5856658 |          0 |
|       6 | AUMCIFO  | AUMCINF\_obs          |   8.3472113 |   8.3472113 |          0 |
|       1 | AUMCPEO  | AUMC\_.Extrap\_obs    |  57.7721236 |  57.7721236 |          0 |
|       2 | AUMCPEO  | AUMC\_.Extrap\_obs    |  31.8178935 |  31.8178935 |          0 |
|       3 | AUMCPEO  | AUMC\_.Extrap\_obs    |  28.0049084 |  28.0049084 |          0 |
|       4 | AUMCPEO  | AUMC\_.Extrap\_obs    |  26.2395827 |  26.2395827 |          0 |
|       5 | AUMCPEO  | AUMC\_.Extrap\_obs    |  43.0991557 |  43.0991557 |          0 |
|       6 | AUMCPEO  | AUMC\_.Extrap\_obs    |  33.0264883 |  33.0264883 |          0 |
|       1 | AUMCIFP  | AUMCINF\_pred         |   7.8485720 |   7.8485720 |          0 |
|       2 | AUMCIFP  | AUMCINF\_pred         |   9.2097614 |   9.2097614 |          0 |
|       3 | AUMCIFP  | AUMCINF\_pred         |   6.0378395 |   6.0378395 |          0 |
|       4 | AUMCIFP  | AUMCINF\_pred         |   5.1204344 |   5.1204344 |          0 |
|       5 | AUMCIFP  | AUMCINF\_pred         |   6.0923336 |   6.0923336 |          0 |
|       6 | AUMCIFP  | AUMCINF\_pred         |   7.3131841 |   7.3131841 |          0 |
|       1 | AUMCPEP  | AUMC\_.Extrap\_pred   |  57.8930274 |  57.8930274 |          0 |
|       2 | AUMCPEP  | AUMC\_.Extrap\_pred   |  30.3655279 |  30.3655279 |          0 |
|       3 | AUMCPEP  | AUMC\_.Extrap\_pred   |  16.2730417 |  16.2730417 |          0 |
|       4 | AUMCPEP  | AUMC\_.Extrap\_pred   |  13.9726930 |  13.9726930 |          0 |
|       5 | AUMCPEP  | AUMC\_.Extrap\_pred   |  38.4915588 |  38.4915588 |          0 |
|       6 | AUMCPEP  | AUMC\_.Extrap\_pred   |  23.5569554 |  23.5569554 |          0 |
|       1 | MRTIVLST | MRTlast               |   1.6442602 |   1.6442602 |          0 |
|       2 | MRTIVLST | MRTlast               |   2.0023083 |   2.0023083 |          0 |
|       3 | MRTIVLST | MRTlast               |   1.4550148 |   1.4550148 |          0 |
|       4 | MRTIVLST | MRTlast               |   1.6027502 |   1.6027502 |          0 |
|       5 | MRTIVLST | MRTlast               |   1.5624335 |   1.5624335 |          0 |
|       6 | MRTIVLST | MRTlast               |   1.6987892 |   1.6987892 |          0 |
|       1 | MRTIVIFO | MRTINF\_obs           |   3.3650320 |   3.3650320 |          0 |
|       2 | MRTIVIFO | MRTINF\_obs           |   2.7125665 |   2.7125665 |          0 |
|       3 | MRTIVIFO | MRTINF\_obs           |   1.9164006 |   1.9164006 |          0 |
|       4 | MRTIVIFO | MRTINF\_obs           |   2.0578350 |   2.0578350 |          0 |
|       5 | MRTIVIFO | MRTINF\_obs           |   2.4985790 |   2.4985790 |          0 |
|       6 | MRTIVIFO | MRTINF\_obs           |   2.3543721 |   2.3543721 |          0 |
|       1 | MRTIVIFP | MRTINF\_pred          |   3.3724181 |   3.3724181 |          0 |
|       2 | MRTIVIFP | MRTINF\_pred          |   2.6693455 |   2.6693455 |          0 |
|       3 | MRTIVIFP | MRTINF\_pred          |   1.6916777 |   1.6916777 |          0 |
|       4 | MRTIVIFP | MRTINF\_pred          |   1.8166868 |   1.8166868 |          0 |
|       5 | MRTIVIFP | MRTINF\_pred          |   2.3481690 |   2.3481690 |          0 |
|       6 | MRTIVIFP | MRTINF\_pred          |   2.1198132 |   2.1198132 |          0 |
|       1 | VSSO     | Vss\_obs              |  36.1720388 |  36.1720388 |          0 |
|       2 | VSSO     | Vss\_obs              |  19.5568334 |  19.5568334 |          0 |
|       3 | VSSO     | Vss\_obs              |  13.0758105 |  13.0758105 |          0 |
|       4 | VSSO     | Vss\_obs              |  17.7272490 |  17.7272490 |          0 |
|       5 | VSSO     | Vss\_obs              |  23.6988080 |  23.6988080 |          0 |
|       6 | VSSO     | Vss\_obs              |  16.6015562 |  16.6015562 |          0 |
|       1 | VSSP     | Vss\_pred             |  36.2269851 |  36.2269851 |          0 |
|       2 | VSSP     | Vss\_pred             |  19.3419923 |  19.3419923 |          0 |
|       3 | VSSP     | Vss\_pred             |  11.8493272 |  11.8493272 |          0 |
|       4 | VSSP     | Vss\_pred             |  16.1136274 |  16.1136274 |          0 |
|       5 | VSSP     | Vss\_pred             |  22.6263772 |  22.6263772 |          0 |
|       6 | VSSP     | Vss\_pred             |  15.3613252 |  15.3613252 |          0 |

Indometh (n=6), Log, IV Bolus

## Test 5: Indometh (n=6), Linear, IV Infusion (0.25hr)

``` r
table_wres_rres(Wres5, Rres5,
                Caption = 'Indometh (n=6), Linear, IV Infusion (0.25hr)')
```

| Subject | PPTESTCD | WNL                   |  NonCompart |   WinNonlin | Difference |
| ------: | :------- | :-------------------- | ----------: | ----------: | ---------: |
|       1 | R2       | Rsq                   |   0.9970667 |   0.9970667 |          0 |
|       2 | R2       | Rsq                   |   0.9476691 |   0.9476691 |          0 |
|       3 | R2       | Rsq                   |   0.8758261 |   0.8758261 |          0 |
|       4 | R2       | Rsq                   |   0.8671179 |   0.8671179 |          0 |
|       5 | R2       | Rsq                   |   0.8752442 |   0.8752442 |          0 |
|       6 | R2       | Rsq                   |   0.9039538 |   0.9039538 |          0 |
|       1 | R2ADJ    | Rsq\_adjusted         |   0.9941335 |   0.9941335 |          0 |
|       2 | R2ADJ    | Rsq\_adjusted         |   0.9401933 |   0.9401933 |          0 |
|       3 | R2ADJ    | Rsq\_adjusted         |   0.8603043 |   0.8603043 |          0 |
|       4 | R2ADJ    | Rsq\_adjusted         |   0.8505077 |   0.8505077 |          0 |
|       5 | R2ADJ    | Rsq\_adjusted         |   0.8544516 |   0.8544516 |          0 |
|       6 | R2ADJ    | Rsq\_adjusted         |   0.8902329 |   0.8902329 |          0 |
|       1 | CORRXY   | Corr\_XY              | \-0.9985323 | \-0.9985323 |          0 |
|       2 | CORRXY   | Corr\_XY              | \-0.9734830 | \-0.9734830 |          0 |
|       3 | CORRXY   | Corr\_XY              | \-0.9358558 | \-0.9358558 |          0 |
|       4 | CORRXY   | Corr\_XY              | \-0.9311917 | \-0.9311917 |          0 |
|       5 | CORRXY   | Corr\_XY              | \-0.9355449 | \-0.9355449 |          0 |
|       6 | CORRXY   | Corr\_XY              | \-0.9507649 | \-0.9507649 |          0 |
|       1 | LAMZNPT  | No\_points\_lambda\_z |   3.0000000 |   3.0000000 |          0 |
|       2 | LAMZNPT  | No\_points\_lambda\_z |   9.0000000 |   9.0000000 |          0 |
|       3 | LAMZNPT  | No\_points\_lambda\_z |  10.0000000 |  10.0000000 |          0 |
|       4 | LAMZNPT  | No\_points\_lambda\_z |  10.0000000 |  10.0000000 |          0 |
|       5 | LAMZNPT  | No\_points\_lambda\_z |   8.0000000 |   8.0000000 |          0 |
|       6 | LAMZNPT  | No\_points\_lambda\_z |   9.0000000 |   9.0000000 |          0 |
|       1 | LAMZ     | Lambda\_z             |   0.1583205 |   0.1583205 |          0 |
|       2 | LAMZ     | Lambda\_z             |   0.3022800 |   0.3022800 |          0 |
|       3 | LAMZ     | Lambda\_z             |   0.4218926 |   0.4218926 |          0 |
|       4 | LAMZ     | Lambda\_z             |   0.4290762 |   0.4290762 |          0 |
|       5 | LAMZ     | Lambda\_z             |   0.2527478 |   0.2527478 |          0 |
|       6 | LAMZ     | Lambda\_z             |   0.3535205 |   0.3535205 |          0 |
|       1 | LAMZLL   | Lambda\_z\_lower      |   5.0000000 |   5.0000000 |          0 |
|       2 | LAMZLL   | Lambda\_z\_lower      |   0.7500000 |   0.7500000 |          0 |
|       3 | LAMZLL   | Lambda\_z\_lower      |   0.5000000 |   0.5000000 |          0 |
|       4 | LAMZLL   | Lambda\_z\_lower      |   0.5000000 |   0.5000000 |          0 |
|       5 | LAMZLL   | Lambda\_z\_lower      |   1.0000000 |   1.0000000 |          0 |
|       6 | LAMZLL   | Lambda\_z\_lower      |   0.7500000 |   0.7500000 |          0 |
|       1 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       2 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       3 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       4 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       5 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       6 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       1 | LAMZHL   | HL\_Lambda\_z         |   4.3781270 |   4.3781270 |          0 |
|       2 | LAMZHL   | HL\_Lambda\_z         |   2.2930632 |   2.2930632 |          0 |
|       3 | LAMZHL   | HL\_Lambda\_z         |   1.6429468 |   1.6429468 |          0 |
|       4 | LAMZHL   | HL\_Lambda\_z         |   1.6154409 |   1.6154409 |          0 |
|       5 | LAMZHL   | HL\_Lambda\_z         |   2.7424461 |   2.7424461 |          0 |
|       6 | LAMZHL   | HL\_Lambda\_z         |   1.9606986 |   1.9606986 |          0 |
|       1 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       2 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       3 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       4 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       5 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       6 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       1 | CMAX     | Cmax                  |   1.5000000 |   1.5000000 |          0 |
|       2 | CMAX     | Cmax                  |   2.0300000 |   2.0300000 |          0 |
|       3 | CMAX     | Cmax                  |   2.7200000 |   2.7200000 |          0 |
|       4 | CMAX     | Cmax                  |   1.8500000 |   1.8500000 |          0 |
|       5 | CMAX     | Cmax                  |   2.0500000 |   2.0500000 |          0 |
|       6 | CMAX     | Cmax                  |   2.3100000 |   2.3100000 |          0 |
|       1 | CMAXD    | Cmax\_D               |   0.0600000 |   0.0600000 |          0 |
|       2 | CMAXD    | Cmax\_D               |   0.0812000 |   0.0812000 |          0 |
|       3 | CMAXD    | Cmax\_D               |   0.1088000 |   0.1088000 |          0 |
|       4 | CMAXD    | Cmax\_D               |   0.0740000 |   0.0740000 |          0 |
|       5 | CMAXD    | Cmax\_D               |   0.0820000 |   0.0820000 |          0 |
|       6 | CMAXD    | Cmax\_D               |   0.0924000 |   0.0924000 |          0 |
|       1 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       2 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       3 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       4 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       5 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       6 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       1 | CLST     | Clast                 |   0.0500000 |   0.0500000 |          0 |
|       2 | CLST     | Clast                 |   0.0800000 |   0.0800000 |          0 |
|       3 | CLST     | Clast                 |   0.0800000 |   0.0800000 |          0 |
|       4 | CLST     | Clast                 |   0.0700000 |   0.0700000 |          0 |
|       5 | CLST     | Clast                 |   0.0600000 |   0.0600000 |          0 |
|       6 | CLST     | Clast                 |   0.0900000 |   0.0900000 |          0 |
|       1 | AUCLST   | AUClast               |   1.7412500 |   1.7412500 |          0 |
|       2 | AUCLST   | AUClast               |   2.9325000 |   2.9325000 |          0 |
|       3 | AUCLST   | AUClast               |   2.9337500 |   2.9337500 |          0 |
|       4 | AUCLST   | AUClast               |   2.4775000 |   2.4775000 |          0 |
|       5 | AUCLST   | AUClast               |   1.9537500 |   1.9537500 |          0 |
|       6 | AUCLST   | AUClast               |   2.8725000 |   2.8725000 |          0 |
|       1 | AUCALL   | AUCall                |   1.7412500 |   1.7412500 |          0 |
|       2 | AUCALL   | AUCall                |   2.9325000 |   2.9325000 |          0 |
|       3 | AUCALL   | AUCall                |   2.9337500 |   2.9337500 |          0 |
|       4 | AUCALL   | AUCall                |   2.4775000 |   2.4775000 |          0 |
|       5 | AUCALL   | AUCall                |   1.9537500 |   1.9537500 |          0 |
|       6 | AUCALL   | AUCall                |   2.8725000 |   2.8725000 |          0 |
|       1 | AUCIFO   | AUCINF\_obs           |   2.0570651 |   2.0570651 |          0 |
|       2 | AUCIFO   | AUCINF\_obs           |   3.1971553 |   3.1971553 |          0 |
|       3 | AUCIFO   | AUCINF\_obs           |   3.1233717 |   3.1233717 |          0 |
|       4 | AUCIFO   | AUCINF\_obs           |   2.6406412 |   2.6406412 |          0 |
|       5 | AUCIFO   | AUCINF\_obs           |   2.1911408 |   2.1911408 |          0 |
|       6 | AUCIFO   | AUCINF\_obs           |   3.1270821 |   3.1270821 |          0 |
|       1 | AUCIFOD  | AUCINF\_D\_obs        |   0.0822826 |   0.0822826 |          0 |
|       2 | AUCIFOD  | AUCINF\_D\_obs        |   0.1278862 |   0.1278862 |          0 |
|       3 | AUCIFOD  | AUCINF\_D\_obs        |   0.1249349 |   0.1249349 |          0 |
|       4 | AUCIFOD  | AUCINF\_D\_obs        |   0.1056256 |   0.1056256 |          0 |
|       5 | AUCIFOD  | AUCINF\_D\_obs        |   0.0876456 |   0.0876456 |          0 |
|       6 | AUCIFOD  | AUCINF\_D\_obs        |   0.1250833 |   0.1250833 |          0 |
|       1 | AUCPEO   | AUC\_.Extrap\_obs     |  15.3527035 |  15.3527035 |          0 |
|       2 | AUCPEO   | AUC\_.Extrap\_obs     |   8.2778360 |   8.2778360 |          0 |
|       3 | AUCPEO   | AUC\_.Extrap\_obs     |   6.0710577 |   6.0710577 |          0 |
|       4 | AUCPEO   | AUC\_.Extrap\_obs     |   6.1780905 |   6.1780905 |          0 |
|       5 | AUCPEO   | AUC\_.Extrap\_obs     |  10.8341191 |  10.8341191 |          0 |
|       6 | AUCPEO   | AUC\_.Extrap\_obs     |   8.1412032 |   8.1412032 |          0 |
|       1 | VZO      | Vz\_obs               |  76.7635175 |  76.7635175 |          0 |
|       2 | VZO      | Vz\_obs               |  25.8682374 |  25.8682374 |          0 |
|       3 | VZO      | Vz\_obs               |  18.9720552 |  18.9720552 |          0 |
|       4 | VZO      | Vz\_obs               |  22.0646091 |  22.0646091 |          0 |
|       5 | VZO      | Vz\_obs               |  45.1421631 |  45.1421631 |          0 |
|       6 | VZO      | Vz\_obs               |  22.6144534 |  22.6144534 |          0 |
|       1 | CLO      | Cl\_obs               |  12.1532371 |  12.1532371 |          0 |
|       2 | CLO      | Cl\_obs               |   7.8194513 |   7.8194513 |          0 |
|       3 | CLO      | Cl\_obs               |   8.0041706 |   8.0041706 |          0 |
|       4 | CLO      | Cl\_obs               |   9.4673975 |   9.4673975 |          0 |
|       5 | CLO      | Cl\_obs               |  11.4095817 |  11.4095817 |          0 |
|       6 | CLO      | Cl\_obs               |   7.9946733 |   7.9946733 |          0 |
|       1 | AUCIFP   | AUCINF\_pred          |   2.0586347 |   2.0586347 |          0 |
|       2 | AUCIFP   | AUCINF\_pred          |   3.1798068 |   3.1798068 |          0 |
|       3 | AUCIFP   | AUCINF\_pred          |   3.0284958 |   3.0284958 |          0 |
|       4 | AUCIFP   | AUCINF\_pred          |   2.5577884 |   2.5577884 |          0 |
|       5 | AUCIFP   | AUCINF\_pred          |   2.1498803 |   2.1498803 |          0 |
|       6 | AUCIFP   | AUCINF\_pred          |   3.0315925 |   3.0315925 |          0 |
|       1 | AUCIFPD  | AUCINF\_D\_pred       |   0.0823454 |   0.0823454 |          0 |
|       2 | AUCIFPD  | AUCINF\_D\_pred       |   0.1271923 |   0.1271923 |          0 |
|       3 | AUCIFPD  | AUCINF\_D\_pred       |   0.1211398 |   0.1211398 |          0 |
|       4 | AUCIFPD  | AUCINF\_D\_pred       |   0.1023115 |   0.1023115 |          0 |
|       5 | AUCIFPD  | AUCINF\_D\_pred       |   0.0859952 |   0.0859952 |          0 |
|       6 | AUCIFPD  | AUCINF\_D\_pred       |   0.1212637 |   0.1212637 |          0 |
|       1 | AUCPEP   | AUC\_.Extrap\_pred    |  15.4172443 |  15.4172443 |          0 |
|       2 | AUCPEP   | AUC\_.Extrap\_pred    |   7.7774164 |   7.7774164 |          0 |
|       3 | AUCPEP   | AUC\_.Extrap\_pred    |   3.1284787 |   3.1284787 |          0 |
|       4 | AUCPEP   | AUC\_.Extrap\_pred    |   3.1389787 |   3.1389787 |          0 |
|       5 | AUCPEP   | AUC\_.Extrap\_pred    |   9.1228460 |   9.1228460 |          0 |
|       6 | AUCPEP   | AUC\_.Extrap\_pred    |   5.2478198 |   5.2478198 |          0 |
|       1 | VZP      | Vz\_pred              |  76.7049878 |  76.7049878 |          0 |
|       2 | VZP      | Vz\_pred              |  26.0093699 |  26.0093699 |          0 |
|       3 | VZP      | Vz\_pred              |  19.5664063 |  19.5664063 |          0 |
|       4 | VZP      | Vz\_pred              |  22.7793336 |  22.7793336 |          0 |
|       5 | VZP      | Vz\_pred              |  46.0085322 |  46.0085322 |          0 |
|       6 | VZP      | Vz\_pred              |  23.3267671 |  23.3267671 |          0 |
|       1 | CLP      | Cl\_pred              |  12.1439707 |  12.1439707 |          0 |
|       2 | CLP      | Cl\_pred              |   7.8621128 |   7.8621128 |          0 |
|       3 | CLP      | Cl\_pred              |   8.2549230 |   8.2549230 |          0 |
|       4 | CLP      | Cl\_pred              |   9.7740688 |   9.7740688 |          0 |
|       5 | CLP      | Cl\_pred              |  11.6285546 |  11.6285546 |          0 |
|       6 | CLP      | Cl\_pred              |   8.2464909 |   8.2464909 |          0 |
|       1 | AUMCLST  | AUMClast              |   3.2712500 |   3.2712500 |          0 |
|       2 | AUMCLST  | AUMClast              |   6.3987500 |   6.3987500 |          0 |
|       3 | AUMCLST  | AUMClast              |   5.0062500 |   5.0062500 |          0 |
|       4 | AUMCLST  | AUMClast              |   4.3818750 |   4.3818750 |          0 |
|       5 | AUMCLST  | AUMClast              |   3.7075000 |   3.7075000 |          0 |
|       6 | AUMCLST  | AUMClast              |   5.5325000 |   5.5325000 |          0 |
|       1 | AUMCIFO  | AUMCINF\_obs          |   7.7925545 |   7.7925545 |          0 |
|       2 | AUMCIFO  | AUMCINF\_obs          |   9.3915223 |   9.3915223 |          0 |
|       3 | AUMCIFO  | AUMCINF\_obs          |   6.9726784 |   6.9726784 |          0 |
|       4 | AUMCIFO  | AUMCINF\_obs          |   6.0672197 |   6.0672197 |          0 |
|       5 | AUMCIFO  | AUMCINF\_obs          |   6.5458663 |   6.5458663 |          0 |
|       6 | AUMCIFO  | AUMCINF\_obs          |   8.2892908 |   8.2892908 |          0 |
|       1 | AUMCPEO  | AUMC\_.Extrap\_obs    |  58.0208261 |  58.0208261 |          0 |
|       2 | AUMCPEO  | AUMC\_.Extrap\_obs    |  31.8667432 |  31.8667432 |          0 |
|       3 | AUMCPEO  | AUMC\_.Extrap\_obs    |  28.2019090 |  28.2019090 |          0 |
|       4 | AUMCPEO  | AUMC\_.Extrap\_obs    |  27.7778746 |  27.7778746 |          0 |
|       5 | AUMCPEO  | AUMC\_.Extrap\_obs    |  43.3612023 |  43.3612023 |          0 |
|       6 | AUMCPEO  | AUMC\_.Extrap\_obs    |  33.2572574 |  33.2572574 |          0 |
|       1 | AUMCIFP  | AUMCINF\_pred         |   7.8150259 |   7.8150259 |          0 |
|       2 | AUMCIFP  | AUMCINF\_pred         |   9.1953427 |   9.1953427 |          0 |
|       3 | AUMCIFP  | AUMCINF\_pred         |   5.9887901 |   5.9887901 |          0 |
|       4 | AUMCIFP  | AUMCINF\_pred         |   5.2113018 |   5.2113018 |          0 |
|       5 | AUMCIFP  | AUMCINF\_pred         |   6.0525342 |   6.0525342 |          0 |
|       6 | AUMCIFP  | AUMCINF\_pred         |   7.2552635 |   7.2552635 |          0 |
|       1 | AUMCPEP  | AUMC\_.Extrap\_pred   |  58.1415337 |  58.1415337 |          0 |
|       2 | AUMCPEP  | AUMC\_.Extrap\_pred   |  30.4131425 |  30.4131425 |          0 |
|       3 | AUMCPEP  | AUMC\_.Extrap\_pred   |  16.4063210 |  16.4063210 |          0 |
|       4 | AUMCPEP  | AUMC\_.Extrap\_pred   |  15.9159231 |  15.9159231 |          0 |
|       5 | AUMCPEP  | AUMC\_.Extrap\_pred   |  38.7446663 |  38.7446663 |          0 |
|       6 | AUMCPEP  | AUMC\_.Extrap\_pred   |  23.7450164 |  23.7450164 |          0 |
|       1 | MRTIVLST | MRTlast               |   1.7536791 |   1.7536791 |          0 |
|       2 | MRTIVLST | MRTlast               |   2.0570119 |   2.0570119 |          0 |
|       3 | MRTIVLST | MRTlast               |   1.5814337 |   1.5814337 |          0 |
|       4 | MRTIVLST | MRTlast               |   1.6436680 |   1.6436680 |          0 |
|       5 | MRTIVLST | MRTlast               |   1.7726328 |   1.7726328 |          0 |
|       6 | MRTIVLST | MRTlast               |   1.8010226 |   1.8010226 |          0 |
|       1 | MRTIVIFO | MRTINF\_obs           |   3.6631905 |   3.6631905 |          0 |
|       2 | MRTIVIFO | MRTINF\_obs           |   2.8124621 |   2.8124621 |          0 |
|       3 | MRTIVIFO | MRTINF\_obs           |   2.1074203 |   2.1074203 |          0 |
|       4 | MRTIVIFO | MRTINF\_obs           |   2.1726312 |   2.1726312 |          0 |
|       5 | MRTIVIFO | MRTINF\_obs           |   2.8624239 |   2.8624239 |          0 |
|       6 | MRTIVIFO | MRTINF\_obs           |   2.5258069 |   2.5258069 |          0 |
|       1 | MRTIVIFP | MRTINF\_pred          |   3.6712178 |   3.6712178 |          0 |
|       2 | MRTIVIFP | MRTINF\_pred          |   2.7667929 |   2.7667929 |          0 |
|       3 | MRTIVIFP | MRTINF\_pred          |   1.8524801 |   1.8524801 |          0 |
|       4 | MRTIVIFP | MRTINF\_pred          |   1.9124249 |   1.9124249 |          0 |
|       5 | MRTIVIFP | MRTINF\_pred          |   2.6902890 |   2.6902890 |          0 |
|       6 | MRTIVIFP | MRTINF\_pred          |   2.2682186 |   2.2682186 |          0 |
|       1 | VSSO     | Vss\_obs              |  44.5196227 |  44.5196227 |          0 |
|       2 | VSSO     | Vss\_obs              |  21.9919102 |  21.9919102 |          0 |
|       3 | VSSO     | Vss\_obs              |  16.8681518 |  16.8681518 |          0 |
|       4 | VSSO     | Vss\_obs              |  20.5691634 |  20.5691634 |          0 |
|       5 | VSSO     | Vss\_obs              |  32.6590590 |  32.6590590 |          0 |
|       6 | VSSO     | Vss\_obs              |  20.1930009 |  20.1930009 |          0 |
|       1 | VSSP     | Vss\_pred             |  44.5831617 |  44.5831617 |          0 |
|       2 | VSSP     | Vss\_pred             |  21.7528377 |  21.7528377 |          0 |
|       3 | VSSP     | Vss\_pred             |  15.2920802 |  15.2920802 |          0 |
|       4 | VSSP     | Vss\_pred             |  18.6921722 |  18.6921722 |          0 |
|       5 | VSSP     | Vss\_pred             |  31.2841719 |  31.2841719 |          0 |
|       6 | VSSP     | Vss\_pred             |  18.7048438 |  18.7048438 |          0 |

Indometh (n=6), Linear, IV Infusion (0.25hr)

## Test 6: Indometh (n=6), Log, IV Infusion (0.25hr)

``` r
table_wres_rres(Wres6, Rres6,
                Caption = 'Indometh (n=6), Log, IV Infusion (0.25hr)')
```

| Subject | PPTESTCD | WNL                   |  NonCompart |   WinNonlin | Difference |
| ------: | :------- | :-------------------- | ----------: | ----------: | ---------: |
|       1 | R2       | Rsq                   |   0.9970667 |   0.9970667 |          0 |
|       2 | R2       | Rsq                   |   0.9476691 |   0.9476691 |          0 |
|       3 | R2       | Rsq                   |   0.8758261 |   0.8758261 |          0 |
|       4 | R2       | Rsq                   |   0.8671179 |   0.8671179 |          0 |
|       5 | R2       | Rsq                   |   0.8752442 |   0.8752442 |          0 |
|       6 | R2       | Rsq                   |   0.9039538 |   0.9039538 |          0 |
|       1 | R2ADJ    | Rsq\_adjusted         |   0.9941335 |   0.9941335 |          0 |
|       2 | R2ADJ    | Rsq\_adjusted         |   0.9401933 |   0.9401933 |          0 |
|       3 | R2ADJ    | Rsq\_adjusted         |   0.8603043 |   0.8603043 |          0 |
|       4 | R2ADJ    | Rsq\_adjusted         |   0.8505077 |   0.8505077 |          0 |
|       5 | R2ADJ    | Rsq\_adjusted         |   0.8544516 |   0.8544516 |          0 |
|       6 | R2ADJ    | Rsq\_adjusted         |   0.8902329 |   0.8902329 |          0 |
|       1 | CORRXY   | Corr\_XY              | \-0.9985323 | \-0.9985323 |          0 |
|       2 | CORRXY   | Corr\_XY              | \-0.9734830 | \-0.9734830 |          0 |
|       3 | CORRXY   | Corr\_XY              | \-0.9358558 | \-0.9358558 |          0 |
|       4 | CORRXY   | Corr\_XY              | \-0.9311917 | \-0.9311917 |          0 |
|       5 | CORRXY   | Corr\_XY              | \-0.9355449 | \-0.9355449 |          0 |
|       6 | CORRXY   | Corr\_XY              | \-0.9507649 | \-0.9507649 |          0 |
|       1 | LAMZNPT  | No\_points\_lambda\_z |   3.0000000 |   3.0000000 |          0 |
|       2 | LAMZNPT  | No\_points\_lambda\_z |   9.0000000 |   9.0000000 |          0 |
|       3 | LAMZNPT  | No\_points\_lambda\_z |  10.0000000 |  10.0000000 |          0 |
|       4 | LAMZNPT  | No\_points\_lambda\_z |  10.0000000 |  10.0000000 |          0 |
|       5 | LAMZNPT  | No\_points\_lambda\_z |   8.0000000 |   8.0000000 |          0 |
|       6 | LAMZNPT  | No\_points\_lambda\_z |   9.0000000 |   9.0000000 |          0 |
|       1 | LAMZ     | Lambda\_z             |   0.1583205 |   0.1583205 |          0 |
|       2 | LAMZ     | Lambda\_z             |   0.3022800 |   0.3022800 |          0 |
|       3 | LAMZ     | Lambda\_z             |   0.4218926 |   0.4218926 |          0 |
|       4 | LAMZ     | Lambda\_z             |   0.4290762 |   0.4290762 |          0 |
|       5 | LAMZ     | Lambda\_z             |   0.2527478 |   0.2527478 |          0 |
|       6 | LAMZ     | Lambda\_z             |   0.3535205 |   0.3535205 |          0 |
|       1 | LAMZLL   | Lambda\_z\_lower      |   5.0000000 |   5.0000000 |          0 |
|       2 | LAMZLL   | Lambda\_z\_lower      |   0.7500000 |   0.7500000 |          0 |
|       3 | LAMZLL   | Lambda\_z\_lower      |   0.5000000 |   0.5000000 |          0 |
|       4 | LAMZLL   | Lambda\_z\_lower      |   0.5000000 |   0.5000000 |          0 |
|       5 | LAMZLL   | Lambda\_z\_lower      |   1.0000000 |   1.0000000 |          0 |
|       6 | LAMZLL   | Lambda\_z\_lower      |   0.7500000 |   0.7500000 |          0 |
|       1 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       2 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       3 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       4 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       5 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       6 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       1 | LAMZHL   | HL\_Lambda\_z         |   4.3781270 |   4.3781270 |          0 |
|       2 | LAMZHL   | HL\_Lambda\_z         |   2.2930632 |   2.2930632 |          0 |
|       3 | LAMZHL   | HL\_Lambda\_z         |   1.6429468 |   1.6429468 |          0 |
|       4 | LAMZHL   | HL\_Lambda\_z         |   1.6154409 |   1.6154409 |          0 |
|       5 | LAMZHL   | HL\_Lambda\_z         |   2.7424461 |   2.7424461 |          0 |
|       6 | LAMZHL   | HL\_Lambda\_z         |   1.9606986 |   1.9606986 |          0 |
|       1 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       2 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       3 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       4 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       5 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       6 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       1 | CMAX     | Cmax                  |   1.5000000 |   1.5000000 |          0 |
|       2 | CMAX     | Cmax                  |   2.0300000 |   2.0300000 |          0 |
|       3 | CMAX     | Cmax                  |   2.7200000 |   2.7200000 |          0 |
|       4 | CMAX     | Cmax                  |   1.8500000 |   1.8500000 |          0 |
|       5 | CMAX     | Cmax                  |   2.0500000 |   2.0500000 |          0 |
|       6 | CMAX     | Cmax                  |   2.3100000 |   2.3100000 |          0 |
|       1 | CMAXD    | Cmax\_D               |   0.0600000 |   0.0600000 |          0 |
|       2 | CMAXD    | Cmax\_D               |   0.0812000 |   0.0812000 |          0 |
|       3 | CMAXD    | Cmax\_D               |   0.1088000 |   0.1088000 |          0 |
|       4 | CMAXD    | Cmax\_D               |   0.0740000 |   0.0740000 |          0 |
|       5 | CMAXD    | Cmax\_D               |   0.0820000 |   0.0820000 |          0 |
|       6 | CMAXD    | Cmax\_D               |   0.0924000 |   0.0924000 |          0 |
|       1 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       2 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       3 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       4 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       5 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       6 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       1 | CLST     | Clast                 |   0.0500000 |   0.0500000 |          0 |
|       2 | CLST     | Clast                 |   0.0800000 |   0.0800000 |          0 |
|       3 | CLST     | Clast                 |   0.0800000 |   0.0800000 |          0 |
|       4 | CLST     | Clast                 |   0.0700000 |   0.0700000 |          0 |
|       5 | CLST     | Clast                 |   0.0600000 |   0.0600000 |          0 |
|       6 | CLST     | Clast                 |   0.0900000 |   0.0900000 |          0 |
|       1 | AUCLST   | AUClast               |   1.7193653 |   1.7193653 |          0 |
|       2 | AUCLST   | AUClast               |   2.8891436 |   2.8891436 |          0 |
|       3 | AUCLST   | AUClast               |   2.8817113 |   2.8817113 |          0 |
|       4 | AUCLST   | AUClast               |   2.4442459 |   2.4442459 |          0 |
|       5 | AUCLST   | AUClast               |   1.9211984 |   1.9211984 |          0 |
|       6 | AUCLST   | AUClast               |   2.8413138 |   2.8413138 |          0 |
|       1 | AUCALL   | AUCall                |   1.7193653 |   1.7193653 |          0 |
|       2 | AUCALL   | AUCall                |   2.8891436 |   2.8891436 |          0 |
|       3 | AUCALL   | AUCall                |   2.8817113 |   2.8817113 |          0 |
|       4 | AUCALL   | AUCall                |   2.4442459 |   2.4442459 |          0 |
|       5 | AUCALL   | AUCall                |   1.9211984 |   1.9211984 |          0 |
|       6 | AUCALL   | AUCall                |   2.8413138 |   2.8413138 |          0 |
|       1 | AUCIFO   | AUCINF\_obs           |   2.0351804 |   2.0351804 |          0 |
|       2 | AUCIFO   | AUCINF\_obs           |   3.1537989 |   3.1537989 |          0 |
|       3 | AUCIFO   | AUCINF\_obs           |   3.0713330 |   3.0713330 |          0 |
|       4 | AUCIFO   | AUCINF\_obs           |   2.6073871 |   2.6073871 |          0 |
|       5 | AUCIFO   | AUCINF\_obs           |   2.1585892 |   2.1585892 |          0 |
|       6 | AUCIFO   | AUCINF\_obs           |   3.0958959 |   3.0958959 |          0 |
|       1 | AUCIFOD  | AUCINF\_D\_obs        |   0.0814072 |   0.0814072 |          0 |
|       2 | AUCIFOD  | AUCINF\_D\_obs        |   0.1261520 |   0.1261520 |          0 |
|       3 | AUCIFOD  | AUCINF\_D\_obs        |   0.1228533 |   0.1228533 |          0 |
|       4 | AUCIFOD  | AUCINF\_D\_obs        |   0.1042955 |   0.1042955 |          0 |
|       5 | AUCIFOD  | AUCINF\_D\_obs        |   0.0863436 |   0.0863436 |          0 |
|       6 | AUCIFOD  | AUCINF\_D\_obs        |   0.1238358 |   0.1238358 |          0 |
|       1 | AUCPEO   | AUC\_.Extrap\_obs     |  15.5177942 |  15.5177942 |          0 |
|       2 | AUCPEO   | AUC\_.Extrap\_obs     |   8.3916343 |   8.3916343 |          0 |
|       3 | AUCPEO   | AUC\_.Extrap\_obs     |   6.1739217 |   6.1739217 |          0 |
|       4 | AUCPEO   | AUC\_.Extrap\_obs     |   6.2568848 |   6.2568848 |          0 |
|       5 | AUCPEO   | AUC\_.Extrap\_obs     |  10.9974979 |  10.9974979 |          0 |
|       6 | AUCPEO   | AUC\_.Extrap\_obs     |   8.2232127 |   8.2232127 |          0 |
|       1 | VZO      | Vz\_obs               |  77.5889712 |  77.5889712 |          0 |
|       2 | VZO      | Vz\_obs               |  26.2238573 |  26.2238573 |          0 |
|       3 | VZO      | Vz\_obs               |  19.2935053 |  19.2935053 |          0 |
|       4 | VZO      | Vz\_obs               |  22.3460171 |  22.3460171 |          0 |
|       5 | VZO      | Vz\_obs               |  45.8229078 |  45.8229078 |          0 |
|       6 | VZO      | Vz\_obs               |  22.8422576 |  22.8422576 |          0 |
|       1 | CLO      | Cl\_obs               |  12.2839234 |  12.2839234 |          0 |
|       2 | CLO      | Cl\_obs               |   7.9269481 |   7.9269481 |          0 |
|       3 | CLO      | Cl\_obs               |   8.1397881 |   8.1397881 |          0 |
|       4 | CLO      | Cl\_obs               |   9.5881430 |   9.5881430 |          0 |
|       5 | CLO      | Cl\_obs               |  11.5816384 |  11.5816384 |          0 |
|       6 | CLO      | Cl\_obs               |   8.0752068 |   8.0752068 |          0 |
|       1 | AUCIFP   | AUCINF\_pred          |   2.0367500 |   2.0367500 |          0 |
|       2 | AUCIFP   | AUCINF\_pred          |   3.1364504 |   3.1364504 |          0 |
|       3 | AUCIFP   | AUCINF\_pred          |   2.9764572 |   2.9764572 |          0 |
|       4 | AUCIFP   | AUCINF\_pred          |   2.5245343 |   2.5245343 |          0 |
|       5 | AUCIFP   | AUCINF\_pred          |   2.1173287 |   2.1173287 |          0 |
|       6 | AUCIFP   | AUCINF\_pred          |   3.0004063 |   3.0004063 |          0 |
|       1 | AUCIFPD  | AUCINF\_D\_pred       |   0.0814700 |   0.0814700 |          0 |
|       2 | AUCIFPD  | AUCINF\_D\_pred       |   0.1254580 |   0.1254580 |          0 |
|       3 | AUCIFPD  | AUCINF\_D\_pred       |   0.1190583 |   0.1190583 |          0 |
|       4 | AUCIFPD  | AUCINF\_D\_pred       |   0.1009814 |   0.1009814 |          0 |
|       5 | AUCIFPD  | AUCINF\_D\_pred       |   0.0846931 |   0.0846931 |          0 |
|       6 | AUCIFPD  | AUCINF\_D\_pred       |   0.1200163 |   0.1200163 |          0 |
|       1 | AUCPEP   | AUC\_.Extrap\_pred    |  15.5829013 |  15.5829013 |          0 |
|       2 | AUCPEP   | AUC\_.Extrap\_pred    |   7.8849267 |   7.8849267 |          0 |
|       3 | AUCPEP   | AUC\_.Extrap\_pred    |   3.1831752 |   3.1831752 |          0 |
|       4 | AUCPEP   | AUC\_.Extrap\_pred    |   3.1803265 |   3.1803265 |          0 |
|       5 | AUCPEP   | AUC\_.Extrap\_pred    |   9.2630996 |   9.2630996 |          0 |
|       6 | AUCPEP   | AUC\_.Extrap\_pred    |   5.3023656 |   5.3023656 |          0 |
|       1 | VZP      | Vz\_pred              |  77.5291765 |  77.5291765 |          0 |
|       2 | VZP      | Vz\_pred              |  26.3689077 |  26.3689077 |          0 |
|       3 | VZP      | Vz\_pred              |  19.9084941 |  19.9084941 |          0 |
|       4 | VZP      | Vz\_pred              |  23.0793917 |  23.0793917 |          0 |
|       5 | VZP      | Vz\_pred              |  46.7158621 |  46.7158621 |          0 |
|       6 | VZP      | Vz\_pred              |  23.5692251 |  23.5692251 |          0 |
|       1 | CLP      | Cl\_pred              |  12.2744566 |  12.2744566 |          0 |
|       2 | CLP      | Cl\_pred              |   7.9707940 |   7.9707940 |          0 |
|       3 | CLP      | Cl\_pred              |   8.3992473 |   8.3992473 |          0 |
|       4 | CLP      | Cl\_pred              |   9.9028165 |   9.9028165 |          0 |
|       5 | CLP      | Cl\_pred              |  11.8073306 |  11.8073306 |          0 |
|       6 | CLP      | Cl\_pred              |   8.3322048 |   8.3322048 |          0 |
|       1 | AUMCLST  | AUMClast              |   3.2965543 |   3.2965543 |          0 |
|       2 | AUMCLST  | AUMClast              |   6.4082620 |   6.4082620 |          0 |
|       3 | AUMCLST  | AUMClast              |   5.0353383 |   5.0353383 |          0 |
|       4 | AUMCLST  | AUMClast              |   4.3990453 |   4.3990453 |          0 |
|       5 | AUMCLST  | AUMClast              |   3.7299741 |   3.7299741 |          0 |
|       6 | AUMCLST  | AUMClast              |   5.5775672 |   5.5775672 |          0 |
|       1 | AUMCIFO  | AUMCINF\_obs          |   7.8178588 |   7.8178588 |          0 |
|       2 | AUMCIFO  | AUMCINF\_obs          |   9.4010343 |   9.4010343 |          0 |
|       3 | AUMCIFO  | AUMCINF\_obs          |   7.0017667 |   7.0017667 |          0 |
|       4 | AUMCIFO  | AUMCINF\_obs          |   6.0843899 |   6.0843899 |          0 |
|       5 | AUMCIFO  | AUMCINF\_obs          |   6.5683405 |   6.5683405 |          0 |
|       6 | AUMCIFO  | AUMCINF\_obs          |   8.3343579 |   8.3343579 |          0 |
|       1 | AUMCPEO  | AUMC\_.Extrap\_obs    |  57.8330281 |  57.8330281 |          0 |
|       2 | AUMCPEO  | AUMC\_.Extrap\_obs    |  31.8345005 |  31.8345005 |          0 |
|       3 | AUMCPEO  | AUMC\_.Extrap\_obs    |  28.0847466 |  28.0847466 |          0 |
|       4 | AUMCPEO  | AUMC\_.Extrap\_obs    |  27.6994849 |  27.6994849 |          0 |
|       5 | AUMCPEO  | AUMC\_.Extrap\_obs    |  43.2128382 |  43.2128382 |          0 |
|       6 | AUMCPEO  | AUMC\_.Extrap\_obs    |  33.0774222 |  33.0774222 |          0 |
|       1 | AUMCIFP  | AUMCINF\_pred         |   7.8403303 |   7.8403303 |          0 |
|       2 | AUMCIFP  | AUMCINF\_pred         |   9.2048546 |   9.2048546 |          0 |
|       3 | AUMCIFP  | AUMCINF\_pred         |   6.0178784 |   6.0178784 |          0 |
|       4 | AUMCIFP  | AUMCINF\_pred         |   5.2284721 |   5.2284721 |          0 |
|       5 | AUMCIFP  | AUMCINF\_pred         |   6.0750083 |   6.0750083 |          0 |
|       6 | AUMCIFP  | AUMCINF\_pred         |   7.3003307 |   7.3003307 |          0 |
|       1 | AUMCPEP  | AUMC\_.Extrap\_pred   |  57.9538845 |  57.9538845 |          0 |
|       2 | AUMCPEP  | AUMC\_.Extrap\_pred   |  30.3817147 |  30.3817147 |          0 |
|       3 | AUMCPEP  | AUMC\_.Extrap\_pred   |  16.3270188 |  16.3270188 |          0 |
|       4 | AUMCPEP  | AUMC\_.Extrap\_pred   |  15.8636553 |  15.8636553 |          0 |
|       5 | AUMCPEP  | AUMC\_.Extrap\_pred   |  38.6013327 |  38.6013327 |          0 |
|       6 | AUMCPEP  | AUMC\_.Extrap\_pred   |  23.5984312 |  23.5984312 |          0 |
|       1 | MRTIVLST | MRTlast               |   1.7923089 |   1.7923089 |          0 |
|       2 | MRTIVLST | MRTlast               |   2.0930490 |   2.0930490 |          0 |
|       3 | MRTIVLST | MRTlast               |   1.6223430 |   1.6223430 |          0 |
|       4 | MRTIVLST | MRTlast               |   1.6747556 |   1.6747556 |          0 |
|       5 | MRTIVLST | MRTlast               |   1.8164830 |   1.8164830 |          0 |
|       6 | MRTIVLST | MRTlast               |   1.8380240 |   1.8380240 |          0 |
|       1 | MRTIVIFO | MRTINF\_obs           |   3.7163591 |   3.7163591 |          0 |
|       2 | MRTIVIFO | MRTINF\_obs           |   2.8558604 |   2.8558604 |          0 |
|       3 | MRTIVIFO | MRTINF\_obs           |   2.1547159 |   2.1547159 |          0 |
|       4 | MRTIVIFO | MRTINF\_obs           |   2.2085200 |   2.2085200 |          0 |
|       5 | MRTIVIFO | MRTINF\_obs           |   2.9178858 |   2.9178858 |          0 |
|       6 | MRTIVIFO | MRTINF\_obs           |   2.5670666 |   2.5670666 |          0 |
|       1 | MRTIVIFP | MRTINF\_pred          |   3.7244318 |   3.7244318 |          0 |
|       2 | MRTIVIFP | MRTINF\_pred          |   2.8098000 |   2.8098000 |          0 |
|       3 | MRTIVIFP | MRTINF\_pred          |   1.8968260 |   1.8968260 |          0 |
|       4 | MRTIVIFP | MRTINF\_pred          |   1.9460640 |   1.9460640 |          0 |
|       5 | MRTIVIFP | MRTINF\_pred          |   2.7441853 |   2.7441853 |          0 |
|       6 | MRTIVIFP | MRTINF\_pred          |   2.3081140 |   2.3081140 |          0 |
|       1 | VSSO     | Vss\_obs              |  45.6514707 |  45.6514707 |          0 |
|       2 | VSSO     | Vss\_obs              |  22.6382575 |  22.6382575 |          0 |
|       3 | VSSO     | Vss\_obs              |  17.5389306 |  17.5389306 |          0 |
|       4 | VSSO     | Vss\_obs              |  21.1756058 |  21.1756058 |          0 |
|       5 | VSSO     | Vss\_obs              |  33.7938980 |  33.7938980 |          0 |
|       6 | VSSO     | Vss\_obs              |  20.7295934 |  20.7295934 |          0 |
|       1 | VSSP     | Vss\_pred             |  45.7153760 |  45.7153760 |          0 |
|       2 | VSSP     | Vss\_pred             |  22.3963368 |  22.3963368 |          0 |
|       3 | VSSP     | Vss\_pred             |  15.9319103 |  15.9319103 |          0 |
|       4 | VSSP     | Vss\_pred             |  19.2715146 |  19.2715146 |          0 |
|       5 | VSSP     | Vss\_pred             |  32.4015028 |  32.4015028 |          0 |
|       6 | VSSP     | Vss\_pred             |  19.2316785 |  19.2316785 |          0 |

Indometh (n=6), Log, IV Infusion (0.25hr)

## Test 7: Indometh (n=6), Linear, Extravascular

``` r
table_wres_rres(Wres7, Rres7,
                Caption = 'Indometh (n=6), Linear, Extravascular')
```

| Subject | PPTESTCD | WNL                   |  NonCompart |   WinNonlin | Difference |
| ------: | :------- | :-------------------- | ----------: | ----------: | ---------: |
|       1 | R2       | Rsq                   |   0.9970667 |   0.9970667 |          0 |
|       2 | R2       | Rsq                   |   0.9476691 |   0.9476691 |          0 |
|       3 | R2       | Rsq                   |   0.8758261 |   0.8758261 |          0 |
|       4 | R2       | Rsq                   |   0.8671179 |   0.8671179 |          0 |
|       5 | R2       | Rsq                   |   0.8752442 |   0.8752442 |          0 |
|       6 | R2       | Rsq                   |   0.9039538 |   0.9039538 |          0 |
|       1 | R2ADJ    | Rsq\_adjusted         |   0.9941335 |   0.9941335 |          0 |
|       2 | R2ADJ    | Rsq\_adjusted         |   0.9401933 |   0.9401933 |          0 |
|       3 | R2ADJ    | Rsq\_adjusted         |   0.8603043 |   0.8603043 |          0 |
|       4 | R2ADJ    | Rsq\_adjusted         |   0.8505077 |   0.8505077 |          0 |
|       5 | R2ADJ    | Rsq\_adjusted         |   0.8544516 |   0.8544516 |          0 |
|       6 | R2ADJ    | Rsq\_adjusted         |   0.8902329 |   0.8902329 |          0 |
|       1 | CORRXY   | Corr\_XY              | \-0.9985323 | \-0.9985323 |          0 |
|       2 | CORRXY   | Corr\_XY              | \-0.9734830 | \-0.9734830 |          0 |
|       3 | CORRXY   | Corr\_XY              | \-0.9358558 | \-0.9358558 |          0 |
|       4 | CORRXY   | Corr\_XY              | \-0.9311917 | \-0.9311917 |          0 |
|       5 | CORRXY   | Corr\_XY              | \-0.9355449 | \-0.9355449 |          0 |
|       6 | CORRXY   | Corr\_XY              | \-0.9507649 | \-0.9507649 |          0 |
|       1 | LAMZNPT  | No\_points\_lambda\_z |   3.0000000 |   3.0000000 |          0 |
|       2 | LAMZNPT  | No\_points\_lambda\_z |   9.0000000 |   9.0000000 |          0 |
|       3 | LAMZNPT  | No\_points\_lambda\_z |  10.0000000 |  10.0000000 |          0 |
|       4 | LAMZNPT  | No\_points\_lambda\_z |  10.0000000 |  10.0000000 |          0 |
|       5 | LAMZNPT  | No\_points\_lambda\_z |   8.0000000 |   8.0000000 |          0 |
|       6 | LAMZNPT  | No\_points\_lambda\_z |   9.0000000 |   9.0000000 |          0 |
|       1 | LAMZ     | Lambda\_z             |   0.1583205 |   0.1583205 |          0 |
|       2 | LAMZ     | Lambda\_z             |   0.3022800 |   0.3022800 |          0 |
|       3 | LAMZ     | Lambda\_z             |   0.4218926 |   0.4218926 |          0 |
|       4 | LAMZ     | Lambda\_z             |   0.4290762 |   0.4290762 |          0 |
|       5 | LAMZ     | Lambda\_z             |   0.2527478 |   0.2527478 |          0 |
|       6 | LAMZ     | Lambda\_z             |   0.3535205 |   0.3535205 |          0 |
|       1 | LAMZLL   | Lambda\_z\_lower      |   5.0000000 |   5.0000000 |          0 |
|       2 | LAMZLL   | Lambda\_z\_lower      |   0.7500000 |   0.7500000 |          0 |
|       3 | LAMZLL   | Lambda\_z\_lower      |   0.5000000 |   0.5000000 |          0 |
|       4 | LAMZLL   | Lambda\_z\_lower      |   0.5000000 |   0.5000000 |          0 |
|       5 | LAMZLL   | Lambda\_z\_lower      |   1.0000000 |   1.0000000 |          0 |
|       6 | LAMZLL   | Lambda\_z\_lower      |   0.7500000 |   0.7500000 |          0 |
|       1 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       2 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       3 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       4 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       5 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       6 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       1 | LAMZHL   | HL\_Lambda\_z         |   4.3781270 |   4.3781270 |          0 |
|       2 | LAMZHL   | HL\_Lambda\_z         |   2.2930632 |   2.2930632 |          0 |
|       3 | LAMZHL   | HL\_Lambda\_z         |   1.6429468 |   1.6429468 |          0 |
|       4 | LAMZHL   | HL\_Lambda\_z         |   1.6154409 |   1.6154409 |          0 |
|       5 | LAMZHL   | HL\_Lambda\_z         |   2.7424461 |   2.7424461 |          0 |
|       6 | LAMZHL   | HL\_Lambda\_z         |   1.9606986 |   1.9606986 |          0 |
|       1 | TLAG     | Tlag                  |   0.0000000 |   0.0000000 |          0 |
|       2 | TLAG     | Tlag                  |   0.0000000 |   0.0000000 |          0 |
|       3 | TLAG     | Tlag                  |   0.0000000 |   0.0000000 |          0 |
|       4 | TLAG     | Tlag                  |   0.0000000 |   0.0000000 |          0 |
|       5 | TLAG     | Tlag                  |   0.0000000 |   0.0000000 |          0 |
|       6 | TLAG     | Tlag                  |   0.0000000 |   0.0000000 |          0 |
|       1 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       2 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       3 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       4 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       5 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       6 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       1 | CMAX     | Cmax                  |   1.5000000 |   1.5000000 |          0 |
|       2 | CMAX     | Cmax                  |   2.0300000 |   2.0300000 |          0 |
|       3 | CMAX     | Cmax                  |   2.7200000 |   2.7200000 |          0 |
|       4 | CMAX     | Cmax                  |   1.8500000 |   1.8500000 |          0 |
|       5 | CMAX     | Cmax                  |   2.0500000 |   2.0500000 |          0 |
|       6 | CMAX     | Cmax                  |   2.3100000 |   2.3100000 |          0 |
|       1 | CMAXD    | Cmax\_D               |   0.0600000 |   0.0600000 |          0 |
|       2 | CMAXD    | Cmax\_D               |   0.0812000 |   0.0812000 |          0 |
|       3 | CMAXD    | Cmax\_D               |   0.1088000 |   0.1088000 |          0 |
|       4 | CMAXD    | Cmax\_D               |   0.0740000 |   0.0740000 |          0 |
|       5 | CMAXD    | Cmax\_D               |   0.0820000 |   0.0820000 |          0 |
|       6 | CMAXD    | Cmax\_D               |   0.0924000 |   0.0924000 |          0 |
|       1 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       2 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       3 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       4 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       5 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       6 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       1 | CLST     | Clast                 |   0.0500000 |   0.0500000 |          0 |
|       2 | CLST     | Clast                 |   0.0800000 |   0.0800000 |          0 |
|       3 | CLST     | Clast                 |   0.0800000 |   0.0800000 |          0 |
|       4 | CLST     | Clast                 |   0.0700000 |   0.0700000 |          0 |
|       5 | CLST     | Clast                 |   0.0600000 |   0.0600000 |          0 |
|       6 | CLST     | Clast                 |   0.0900000 |   0.0900000 |          0 |
|       1 | AUCLST   | AUClast               |   1.7412500 |   1.7412500 |          0 |
|       2 | AUCLST   | AUClast               |   2.9325000 |   2.9325000 |          0 |
|       3 | AUCLST   | AUClast               |   2.9337500 |   2.9337500 |          0 |
|       4 | AUCLST   | AUClast               |   2.4775000 |   2.4775000 |          0 |
|       5 | AUCLST   | AUClast               |   1.9537500 |   1.9537500 |          0 |
|       6 | AUCLST   | AUClast               |   2.8725000 |   2.8725000 |          0 |
|       1 | AUCALL   | AUCall                |   1.7412500 |   1.7412500 |          0 |
|       2 | AUCALL   | AUCall                |   2.9325000 |   2.9325000 |          0 |
|       3 | AUCALL   | AUCall                |   2.9337500 |   2.9337500 |          0 |
|       4 | AUCALL   | AUCall                |   2.4775000 |   2.4775000 |          0 |
|       5 | AUCALL   | AUCall                |   1.9537500 |   1.9537500 |          0 |
|       6 | AUCALL   | AUCall                |   2.8725000 |   2.8725000 |          0 |
|       1 | AUCIFO   | AUCINF\_obs           |   2.0570651 |   2.0570651 |          0 |
|       2 | AUCIFO   | AUCINF\_obs           |   3.1971553 |   3.1971553 |          0 |
|       3 | AUCIFO   | AUCINF\_obs           |   3.1233717 |   3.1233717 |          0 |
|       4 | AUCIFO   | AUCINF\_obs           |   2.6406412 |   2.6406412 |          0 |
|       5 | AUCIFO   | AUCINF\_obs           |   2.1911408 |   2.1911408 |          0 |
|       6 | AUCIFO   | AUCINF\_obs           |   3.1270821 |   3.1270821 |          0 |
|       1 | AUCIFOD  | AUCINF\_D\_obs        |   0.0822826 |   0.0822826 |          0 |
|       2 | AUCIFOD  | AUCINF\_D\_obs        |   0.1278862 |   0.1278862 |          0 |
|       3 | AUCIFOD  | AUCINF\_D\_obs        |   0.1249349 |   0.1249349 |          0 |
|       4 | AUCIFOD  | AUCINF\_D\_obs        |   0.1056256 |   0.1056256 |          0 |
|       5 | AUCIFOD  | AUCINF\_D\_obs        |   0.0876456 |   0.0876456 |          0 |
|       6 | AUCIFOD  | AUCINF\_D\_obs        |   0.1250833 |   0.1250833 |          0 |
|       1 | AUCPEO   | AUC\_.Extrap\_obs     |  15.3527035 |  15.3527035 |          0 |
|       2 | AUCPEO   | AUC\_.Extrap\_obs     |   8.2778360 |   8.2778360 |          0 |
|       3 | AUCPEO   | AUC\_.Extrap\_obs     |   6.0710577 |   6.0710577 |          0 |
|       4 | AUCPEO   | AUC\_.Extrap\_obs     |   6.1780905 |   6.1780905 |          0 |
|       5 | AUCPEO   | AUC\_.Extrap\_obs     |  10.8341191 |  10.8341191 |          0 |
|       6 | AUCPEO   | AUC\_.Extrap\_obs     |   8.1412032 |   8.1412032 |          0 |
|       1 | VZFO     | Vz\_F\_obs            |  76.7635175 |  76.7635175 |          0 |
|       2 | VZFO     | Vz\_F\_obs            |  25.8682374 |  25.8682374 |          0 |
|       3 | VZFO     | Vz\_F\_obs            |  18.9720552 |  18.9720552 |          0 |
|       4 | VZFO     | Vz\_F\_obs            |  22.0646091 |  22.0646091 |          0 |
|       5 | VZFO     | Vz\_F\_obs            |  45.1421631 |  45.1421631 |          0 |
|       6 | VZFO     | Vz\_F\_obs            |  22.6144534 |  22.6144534 |          0 |
|       1 | CLFO     | Cl\_F\_obs            |  12.1532371 |  12.1532371 |          0 |
|       2 | CLFO     | Cl\_F\_obs            |   7.8194513 |   7.8194513 |          0 |
|       3 | CLFO     | Cl\_F\_obs            |   8.0041706 |   8.0041706 |          0 |
|       4 | CLFO     | Cl\_F\_obs            |   9.4673975 |   9.4673975 |          0 |
|       5 | CLFO     | Cl\_F\_obs            |  11.4095817 |  11.4095817 |          0 |
|       6 | CLFO     | Cl\_F\_obs            |   7.9946733 |   7.9946733 |          0 |
|       1 | AUCIFP   | AUCINF\_pred          |   2.0586347 |   2.0586347 |          0 |
|       2 | AUCIFP   | AUCINF\_pred          |   3.1798068 |   3.1798068 |          0 |
|       3 | AUCIFP   | AUCINF\_pred          |   3.0284958 |   3.0284958 |          0 |
|       4 | AUCIFP   | AUCINF\_pred          |   2.5577884 |   2.5577884 |          0 |
|       5 | AUCIFP   | AUCINF\_pred          |   2.1498803 |   2.1498803 |          0 |
|       6 | AUCIFP   | AUCINF\_pred          |   3.0315925 |   3.0315925 |          0 |
|       1 | AUCIFPD  | AUCINF\_D\_pred       |   0.0823454 |   0.0823454 |          0 |
|       2 | AUCIFPD  | AUCINF\_D\_pred       |   0.1271923 |   0.1271923 |          0 |
|       3 | AUCIFPD  | AUCINF\_D\_pred       |   0.1211398 |   0.1211398 |          0 |
|       4 | AUCIFPD  | AUCINF\_D\_pred       |   0.1023115 |   0.1023115 |          0 |
|       5 | AUCIFPD  | AUCINF\_D\_pred       |   0.0859952 |   0.0859952 |          0 |
|       6 | AUCIFPD  | AUCINF\_D\_pred       |   0.1212637 |   0.1212637 |          0 |
|       1 | AUCPEP   | AUC\_.Extrap\_pred    |  15.4172443 |  15.4172443 |          0 |
|       2 | AUCPEP   | AUC\_.Extrap\_pred    |   7.7774164 |   7.7774164 |          0 |
|       3 | AUCPEP   | AUC\_.Extrap\_pred    |   3.1284787 |   3.1284787 |          0 |
|       4 | AUCPEP   | AUC\_.Extrap\_pred    |   3.1389787 |   3.1389787 |          0 |
|       5 | AUCPEP   | AUC\_.Extrap\_pred    |   9.1228460 |   9.1228460 |          0 |
|       6 | AUCPEP   | AUC\_.Extrap\_pred    |   5.2478198 |   5.2478198 |          0 |
|       1 | VZFP     | Vz\_F\_pred           |  76.7049878 |  76.7049878 |          0 |
|       2 | VZFP     | Vz\_F\_pred           |  26.0093699 |  26.0093699 |          0 |
|       3 | VZFP     | Vz\_F\_pred           |  19.5664063 |  19.5664063 |          0 |
|       4 | VZFP     | Vz\_F\_pred           |  22.7793336 |  22.7793336 |          0 |
|       5 | VZFP     | Vz\_F\_pred           |  46.0085322 |  46.0085322 |          0 |
|       6 | VZFP     | Vz\_F\_pred           |  23.3267671 |  23.3267671 |          0 |
|       1 | CLFP     | Cl\_F\_pred           |  12.1439707 |  12.1439707 |          0 |
|       2 | CLFP     | Cl\_F\_pred           |   7.8621128 |   7.8621128 |          0 |
|       3 | CLFP     | Cl\_F\_pred           |   8.2549230 |   8.2549230 |          0 |
|       4 | CLFP     | Cl\_F\_pred           |   9.7740688 |   9.7740687 |          0 |
|       5 | CLFP     | Cl\_F\_pred           |  11.6285546 |  11.6285546 |          0 |
|       6 | CLFP     | Cl\_F\_pred           |   8.2464909 |   8.2464909 |          0 |
|       1 | AUMCLST  | AUMClast              |   3.2712500 |   3.2712500 |          0 |
|       2 | AUMCLST  | AUMClast              |   6.3987500 |   6.3987500 |          0 |
|       3 | AUMCLST  | AUMClast              |   5.0062500 |   5.0062500 |          0 |
|       4 | AUMCLST  | AUMClast              |   4.3818750 |   4.3818750 |          0 |
|       5 | AUMCLST  | AUMClast              |   3.7075000 |   3.7075000 |          0 |
|       6 | AUMCLST  | AUMClast              |   5.5325000 |   5.5325000 |          0 |
|       1 | AUMCIFO  | AUMCINF\_obs          |   7.7925545 |   7.7925545 |          0 |
|       2 | AUMCIFO  | AUMCINF\_obs          |   9.3915223 |   9.3915223 |          0 |
|       3 | AUMCIFO  | AUMCINF\_obs          |   6.9726784 |   6.9726784 |          0 |
|       4 | AUMCIFO  | AUMCINF\_obs          |   6.0672197 |   6.0672197 |          0 |
|       5 | AUMCIFO  | AUMCINF\_obs          |   6.5458663 |   6.5458663 |          0 |
|       6 | AUMCIFO  | AUMCINF\_obs          |   8.2892908 |   8.2892908 |          0 |
|       1 | AUMCPEO  | AUMC\_.Extrap\_obs    |  58.0208261 |  58.0208261 |          0 |
|       2 | AUMCPEO  | AUMC\_.Extrap\_obs    |  31.8667432 |  31.8667432 |          0 |
|       3 | AUMCPEO  | AUMC\_.Extrap\_obs    |  28.2019090 |  28.2019090 |          0 |
|       4 | AUMCPEO  | AUMC\_.Extrap\_obs    |  27.7778746 |  27.7778746 |          0 |
|       5 | AUMCPEO  | AUMC\_.Extrap\_obs    |  43.3612023 |  43.3612023 |          0 |
|       6 | AUMCPEO  | AUMC\_.Extrap\_obs    |  33.2572574 |  33.2572574 |          0 |
|       1 | AUMCIFP  | AUMCINF\_pred         |   7.8150259 |   7.8150259 |          0 |
|       2 | AUMCIFP  | AUMCINF\_pred         |   9.1953427 |   9.1953427 |          0 |
|       3 | AUMCIFP  | AUMCINF\_pred         |   5.9887901 |   5.9887901 |          0 |
|       4 | AUMCIFP  | AUMCINF\_pred         |   5.2113018 |   5.2113018 |          0 |
|       5 | AUMCIFP  | AUMCINF\_pred         |   6.0525342 |   6.0525342 |          0 |
|       6 | AUMCIFP  | AUMCINF\_pred         |   7.2552635 |   7.2552635 |          0 |
|       1 | AUMCPEP  | AUMC\_.Extrap\_pred   |  58.1415337 |  58.1415337 |          0 |
|       2 | AUMCPEP  | AUMC\_.Extrap\_pred   |  30.4131425 |  30.4131425 |          0 |
|       3 | AUMCPEP  | AUMC\_.Extrap\_pred   |  16.4063210 |  16.4063210 |          0 |
|       4 | AUMCPEP  | AUMC\_.Extrap\_pred   |  15.9159231 |  15.9159230 |          0 |
|       5 | AUMCPEP  | AUMC\_.Extrap\_pred   |  38.7446663 |  38.7446663 |          0 |
|       6 | AUMCPEP  | AUMC\_.Extrap\_pred   |  23.7450164 |  23.7450164 |          0 |
|       1 | MRTEVLST | MRTlast               |   1.8786791 |   1.8786791 |          0 |
|       2 | MRTEVLST | MRTlast               |   2.1820119 |   2.1820119 |          0 |
|       3 | MRTEVLST | MRTlast               |   1.7064337 |   1.7064337 |          0 |
|       4 | MRTEVLST | MRTlast               |   1.7686680 |   1.7686680 |          0 |
|       5 | MRTEVLST | MRTlast               |   1.8976328 |   1.8976328 |          0 |
|       6 | MRTEVLST | MRTlast               |   1.9260226 |   1.9260226 |          0 |
|       1 | MRTEVIFO | MRTINF\_obs           |   3.7881905 |   3.7881905 |          0 |
|       2 | MRTEVIFO | MRTINF\_obs           |   2.9374621 |   2.9374621 |          0 |
|       3 | MRTEVIFO | MRTINF\_obs           |   2.2324203 |   2.2324203 |          0 |
|       4 | MRTEVIFO | MRTINF\_obs           |   2.2976312 |   2.2976312 |          0 |
|       5 | MRTEVIFO | MRTINF\_obs           |   2.9874239 |   2.9874239 |          0 |
|       6 | MRTEVIFO | MRTINF\_obs           |   2.6508069 |   2.6508069 |          0 |
|       1 | MRTEVIFP | MRTINF\_pred          |   3.7962178 |   3.7962178 |          0 |
|       2 | MRTEVIFP | MRTINF\_pred          |   2.8917929 |   2.8917929 |          0 |
|       3 | MRTEVIFP | MRTINF\_pred          |   1.9774801 |   1.9774801 |          0 |
|       4 | MRTEVIFP | MRTINF\_pred          |   2.0374249 |   2.0374249 |          0 |
|       5 | MRTEVIFP | MRTINF\_pred          |   2.8152890 |   2.8152890 |          0 |
|       6 | MRTEVIFP | MRTINF\_pred          |   2.3932186 |   2.3932186 |          0 |

Indometh (n=6), Linear, Extravascular

## Test 8: Indometh (n=6), Log, Extravascular

``` r
table_wres_rres(Wres8, Rres8,
                Caption = 'Indometh (n=6), Log, Extravascular')
```

| Subject | PPTESTCD | WNL                   |  NonCompart |   WinNonlin | Difference |
| ------: | :------- | :-------------------- | ----------: | ----------: | ---------: |
|       1 | R2       | Rsq                   |   0.9970667 |   0.9970667 |          0 |
|       2 | R2       | Rsq                   |   0.9476691 |   0.9476691 |          0 |
|       3 | R2       | Rsq                   |   0.8758261 |   0.8758261 |          0 |
|       4 | R2       | Rsq                   |   0.8671179 |   0.8671179 |          0 |
|       5 | R2       | Rsq                   |   0.8752442 |   0.8752442 |          0 |
|       6 | R2       | Rsq                   |   0.9039538 |   0.9039538 |          0 |
|       1 | R2ADJ    | Rsq\_adjusted         |   0.9941335 |   0.9941335 |          0 |
|       2 | R2ADJ    | Rsq\_adjusted         |   0.9401933 |   0.9401933 |          0 |
|       3 | R2ADJ    | Rsq\_adjusted         |   0.8603043 |   0.8603043 |          0 |
|       4 | R2ADJ    | Rsq\_adjusted         |   0.8505077 |   0.8505077 |          0 |
|       5 | R2ADJ    | Rsq\_adjusted         |   0.8544516 |   0.8544516 |          0 |
|       6 | R2ADJ    | Rsq\_adjusted         |   0.8902329 |   0.8902329 |          0 |
|       1 | CORRXY   | Corr\_XY              | \-0.9985323 | \-0.9985323 |          0 |
|       2 | CORRXY   | Corr\_XY              | \-0.9734830 | \-0.9734830 |          0 |
|       3 | CORRXY   | Corr\_XY              | \-0.9358558 | \-0.9358558 |          0 |
|       4 | CORRXY   | Corr\_XY              | \-0.9311917 | \-0.9311917 |          0 |
|       5 | CORRXY   | Corr\_XY              | \-0.9355449 | \-0.9355449 |          0 |
|       6 | CORRXY   | Corr\_XY              | \-0.9507649 | \-0.9507649 |          0 |
|       1 | LAMZNPT  | No\_points\_lambda\_z |   3.0000000 |   3.0000000 |          0 |
|       2 | LAMZNPT  | No\_points\_lambda\_z |   9.0000000 |   9.0000000 |          0 |
|       3 | LAMZNPT  | No\_points\_lambda\_z |  10.0000000 |  10.0000000 |          0 |
|       4 | LAMZNPT  | No\_points\_lambda\_z |  10.0000000 |  10.0000000 |          0 |
|       5 | LAMZNPT  | No\_points\_lambda\_z |   8.0000000 |   8.0000000 |          0 |
|       6 | LAMZNPT  | No\_points\_lambda\_z |   9.0000000 |   9.0000000 |          0 |
|       1 | LAMZ     | Lambda\_z             |   0.1583205 |   0.1583205 |          0 |
|       2 | LAMZ     | Lambda\_z             |   0.3022800 |   0.3022800 |          0 |
|       3 | LAMZ     | Lambda\_z             |   0.4218926 |   0.4218926 |          0 |
|       4 | LAMZ     | Lambda\_z             |   0.4290762 |   0.4290762 |          0 |
|       5 | LAMZ     | Lambda\_z             |   0.2527478 |   0.2527478 |          0 |
|       6 | LAMZ     | Lambda\_z             |   0.3535205 |   0.3535205 |          0 |
|       1 | LAMZLL   | Lambda\_z\_lower      |   5.0000000 |   5.0000000 |          0 |
|       2 | LAMZLL   | Lambda\_z\_lower      |   0.7500000 |   0.7500000 |          0 |
|       3 | LAMZLL   | Lambda\_z\_lower      |   0.5000000 |   0.5000000 |          0 |
|       4 | LAMZLL   | Lambda\_z\_lower      |   0.5000000 |   0.5000000 |          0 |
|       5 | LAMZLL   | Lambda\_z\_lower      |   1.0000000 |   1.0000000 |          0 |
|       6 | LAMZLL   | Lambda\_z\_lower      |   0.7500000 |   0.7500000 |          0 |
|       1 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       2 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       3 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       4 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       5 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       6 | LAMZUL   | Lambda\_z\_upper      |   8.0000000 |   8.0000000 |          0 |
|       1 | LAMZHL   | HL\_Lambda\_z         |   4.3781270 |   4.3781270 |          0 |
|       2 | LAMZHL   | HL\_Lambda\_z         |   2.2930632 |   2.2930632 |          0 |
|       3 | LAMZHL   | HL\_Lambda\_z         |   1.6429468 |   1.6429468 |          0 |
|       4 | LAMZHL   | HL\_Lambda\_z         |   1.6154409 |   1.6154409 |          0 |
|       5 | LAMZHL   | HL\_Lambda\_z         |   2.7424461 |   2.7424461 |          0 |
|       6 | LAMZHL   | HL\_Lambda\_z         |   1.9606986 |   1.9606986 |          0 |
|       1 | TLAG     | Tlag                  |   0.0000000 |   0.0000000 |          0 |
|       2 | TLAG     | Tlag                  |   0.0000000 |   0.0000000 |          0 |
|       3 | TLAG     | Tlag                  |   0.0000000 |   0.0000000 |          0 |
|       4 | TLAG     | Tlag                  |   0.0000000 |   0.0000000 |          0 |
|       5 | TLAG     | Tlag                  |   0.0000000 |   0.0000000 |          0 |
|       6 | TLAG     | Tlag                  |   0.0000000 |   0.0000000 |          0 |
|       1 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       2 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       3 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       4 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       5 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       6 | TMAX     | Tmax                  |   0.2500000 |   0.2500000 |          0 |
|       1 | CMAX     | Cmax                  |   1.5000000 |   1.5000000 |          0 |
|       2 | CMAX     | Cmax                  |   2.0300000 |   2.0300000 |          0 |
|       3 | CMAX     | Cmax                  |   2.7200000 |   2.7200000 |          0 |
|       4 | CMAX     | Cmax                  |   1.8500000 |   1.8500000 |          0 |
|       5 | CMAX     | Cmax                  |   2.0500000 |   2.0500000 |          0 |
|       6 | CMAX     | Cmax                  |   2.3100000 |   2.3100000 |          0 |
|       1 | CMAXD    | Cmax\_D               |   0.0600000 |   0.0600000 |          0 |
|       2 | CMAXD    | Cmax\_D               |   0.0812000 |   0.0812000 |          0 |
|       3 | CMAXD    | Cmax\_D               |   0.1088000 |   0.1088000 |          0 |
|       4 | CMAXD    | Cmax\_D               |   0.0740000 |   0.0740000 |          0 |
|       5 | CMAXD    | Cmax\_D               |   0.0820000 |   0.0820000 |          0 |
|       6 | CMAXD    | Cmax\_D               |   0.0924000 |   0.0924000 |          0 |
|       1 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       2 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       3 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       4 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       5 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       6 | TLST     | Tlast                 |   8.0000000 |   8.0000000 |          0 |
|       1 | CLST     | Clast                 |   0.0500000 |   0.0500000 |          0 |
|       2 | CLST     | Clast                 |   0.0800000 |   0.0800000 |          0 |
|       3 | CLST     | Clast                 |   0.0800000 |   0.0800000 |          0 |
|       4 | CLST     | Clast                 |   0.0700000 |   0.0700000 |          0 |
|       5 | CLST     | Clast                 |   0.0600000 |   0.0600000 |          0 |
|       6 | CLST     | Clast                 |   0.0900000 |   0.0900000 |          0 |
|       1 | AUCLST   | AUClast               |   1.7193653 |   1.7193653 |          0 |
|       2 | AUCLST   | AUClast               |   2.8891436 |   2.8891436 |          0 |
|       3 | AUCLST   | AUClast               |   2.8817113 |   2.8817113 |          0 |
|       4 | AUCLST   | AUClast               |   2.4442459 |   2.4442459 |          0 |
|       5 | AUCLST   | AUClast               |   1.9211984 |   1.9211984 |          0 |
|       6 | AUCLST   | AUClast               |   2.8413138 |   2.8413138 |          0 |
|       1 | AUCALL   | AUCall                |   1.7193653 |   1.7193653 |          0 |
|       2 | AUCALL   | AUCall                |   2.8891436 |   2.8891436 |          0 |
|       3 | AUCALL   | AUCall                |   2.8817113 |   2.8817113 |          0 |
|       4 | AUCALL   | AUCall                |   2.4442459 |   2.4442459 |          0 |
|       5 | AUCALL   | AUCall                |   1.9211984 |   1.9211984 |          0 |
|       6 | AUCALL   | AUCall                |   2.8413138 |   2.8413138 |          0 |
|       1 | AUCIFO   | AUCINF\_obs           |   2.0351804 |   2.0351804 |          0 |
|       2 | AUCIFO   | AUCINF\_obs           |   3.1537989 |   3.1537989 |          0 |
|       3 | AUCIFO   | AUCINF\_obs           |   3.0713330 |   3.0713330 |          0 |
|       4 | AUCIFO   | AUCINF\_obs           |   2.6073871 |   2.6073871 |          0 |
|       5 | AUCIFO   | AUCINF\_obs           |   2.1585892 |   2.1585892 |          0 |
|       6 | AUCIFO   | AUCINF\_obs           |   3.0958959 |   3.0958959 |          0 |
|       1 | AUCIFOD  | AUCINF\_D\_obs        |   0.0814072 |   0.0814072 |          0 |
|       2 | AUCIFOD  | AUCINF\_D\_obs        |   0.1261520 |   0.1261520 |          0 |
|       3 | AUCIFOD  | AUCINF\_D\_obs        |   0.1228533 |   0.1228533 |          0 |
|       4 | AUCIFOD  | AUCINF\_D\_obs        |   0.1042955 |   0.1042955 |          0 |
|       5 | AUCIFOD  | AUCINF\_D\_obs        |   0.0863436 |   0.0863436 |          0 |
|       6 | AUCIFOD  | AUCINF\_D\_obs        |   0.1238358 |   0.1238358 |          0 |
|       1 | AUCPEO   | AUC\_.Extrap\_obs     |  15.5177942 |  15.5177942 |          0 |
|       2 | AUCPEO   | AUC\_.Extrap\_obs     |   8.3916343 |   8.3916343 |          0 |
|       3 | AUCPEO   | AUC\_.Extrap\_obs     |   6.1739217 |   6.1739217 |          0 |
|       4 | AUCPEO   | AUC\_.Extrap\_obs     |   6.2568848 |   6.2568848 |          0 |
|       5 | AUCPEO   | AUC\_.Extrap\_obs     |  10.9974979 |  10.9974979 |          0 |
|       6 | AUCPEO   | AUC\_.Extrap\_obs     |   8.2232127 |   8.2232127 |          0 |
|       1 | VZFO     | Vz\_F\_obs            |  77.5889712 |  77.5889712 |          0 |
|       2 | VZFO     | Vz\_F\_obs            |  26.2238573 |  26.2238573 |          0 |
|       3 | VZFO     | Vz\_F\_obs            |  19.2935053 |  19.2935053 |          0 |
|       4 | VZFO     | Vz\_F\_obs            |  22.3460171 |  22.3460171 |          0 |
|       5 | VZFO     | Vz\_F\_obs            |  45.8229078 |  45.8229078 |          0 |
|       6 | VZFO     | Vz\_F\_obs            |  22.8422576 |  22.8422576 |          0 |
|       1 | CLFO     | Cl\_F\_obs            |  12.2839234 |  12.2839234 |          0 |
|       2 | CLFO     | Cl\_F\_obs            |   7.9269481 |   7.9269481 |          0 |
|       3 | CLFO     | Cl\_F\_obs            |   8.1397881 |   8.1397881 |          0 |
|       4 | CLFO     | Cl\_F\_obs            |   9.5881430 |   9.5881430 |          0 |
|       5 | CLFO     | Cl\_F\_obs            |  11.5816384 |  11.5816384 |          0 |
|       6 | CLFO     | Cl\_F\_obs            |   8.0752068 |   8.0752068 |          0 |
|       1 | AUCIFP   | AUCINF\_pred          |   2.0367500 |   2.0367500 |          0 |
|       2 | AUCIFP   | AUCINF\_pred          |   3.1364504 |   3.1364504 |          0 |
|       3 | AUCIFP   | AUCINF\_pred          |   2.9764572 |   2.9764572 |          0 |
|       4 | AUCIFP   | AUCINF\_pred          |   2.5245343 |   2.5245343 |          0 |
|       5 | AUCIFP   | AUCINF\_pred          |   2.1173287 |   2.1173287 |          0 |
|       6 | AUCIFP   | AUCINF\_pred          |   3.0004063 |   3.0004063 |          0 |
|       1 | AUCIFPD  | AUCINF\_D\_pred       |   0.0814700 |   0.0814700 |          0 |
|       2 | AUCIFPD  | AUCINF\_D\_pred       |   0.1254580 |   0.1254580 |          0 |
|       3 | AUCIFPD  | AUCINF\_D\_pred       |   0.1190583 |   0.1190583 |          0 |
|       4 | AUCIFPD  | AUCINF\_D\_pred       |   0.1009814 |   0.1009814 |          0 |
|       5 | AUCIFPD  | AUCINF\_D\_pred       |   0.0846931 |   0.0846931 |          0 |
|       6 | AUCIFPD  | AUCINF\_D\_pred       |   0.1200163 |   0.1200163 |          0 |
|       1 | AUCPEP   | AUC\_.Extrap\_pred    |  15.5829013 |  15.5829013 |          0 |
|       2 | AUCPEP   | AUC\_.Extrap\_pred    |   7.8849267 |   7.8849267 |          0 |
|       3 | AUCPEP   | AUC\_.Extrap\_pred    |   3.1831752 |   3.1831752 |          0 |
|       4 | AUCPEP   | AUC\_.Extrap\_pred    |   3.1803265 |   3.1803265 |          0 |
|       5 | AUCPEP   | AUC\_.Extrap\_pred    |   9.2630996 |   9.2630996 |          0 |
|       6 | AUCPEP   | AUC\_.Extrap\_pred    |   5.3023656 |   5.3023656 |          0 |
|       1 | VZFP     | Vz\_F\_pred           |  77.5291765 |  77.5291765 |          0 |
|       2 | VZFP     | Vz\_F\_pred           |  26.3689077 |  26.3689077 |          0 |
|       3 | VZFP     | Vz\_F\_pred           |  19.9084941 |  19.9084941 |          0 |
|       4 | VZFP     | Vz\_F\_pred           |  23.0793917 |  23.0793917 |          0 |
|       5 | VZFP     | Vz\_F\_pred           |  46.7158621 |  46.7158621 |          0 |
|       6 | VZFP     | Vz\_F\_pred           |  23.5692251 |  23.5692252 |          0 |
|       1 | CLFP     | Cl\_F\_pred           |  12.2744566 |  12.2744566 |          0 |
|       2 | CLFP     | Cl\_F\_pred           |   7.9707940 |   7.9707940 |          0 |
|       3 | CLFP     | Cl\_F\_pred           |   8.3992473 |   8.3992473 |          0 |
|       4 | CLFP     | Cl\_F\_pred           |   9.9028165 |   9.9028165 |          0 |
|       5 | CLFP     | Cl\_F\_pred           |  11.8073306 |  11.8073306 |          0 |
|       6 | CLFP     | Cl\_F\_pred           |   8.3322048 |   8.3322048 |          0 |
|       1 | AUMCLST  | AUMClast              |   3.2965543 |   3.2965543 |          0 |
|       2 | AUMCLST  | AUMClast              |   6.4082620 |   6.4082620 |          0 |
|       3 | AUMCLST  | AUMClast              |   5.0353383 |   5.0353382 |          0 |
|       4 | AUMCLST  | AUMClast              |   4.3990453 |   4.3990453 |          0 |
|       5 | AUMCLST  | AUMClast              |   3.7299741 |   3.7299741 |          0 |
|       6 | AUMCLST  | AUMClast              |   5.5775672 |   5.5775672 |          0 |
|       1 | AUMCIFO  | AUMCINF\_obs          |   7.8178588 |   7.8178588 |          0 |
|       2 | AUMCIFO  | AUMCINF\_obs          |   9.4010343 |   9.4010343 |          0 |
|       3 | AUMCIFO  | AUMCINF\_obs          |   7.0017667 |   7.0017667 |          0 |
|       4 | AUMCIFO  | AUMCINF\_obs          |   6.0843899 |   6.0843899 |          0 |
|       5 | AUMCIFO  | AUMCINF\_obs          |   6.5683405 |   6.5683405 |          0 |
|       6 | AUMCIFO  | AUMCINF\_obs          |   8.3343579 |   8.3343579 |          0 |
|       1 | AUMCPEO  | AUMC\_.Extrap\_obs    |  57.8330281 |  57.8330281 |          0 |
|       2 | AUMCPEO  | AUMC\_.Extrap\_obs    |  31.8345005 |  31.8345005 |          0 |
|       3 | AUMCPEO  | AUMC\_.Extrap\_obs    |  28.0847466 |  28.0847466 |          0 |
|       4 | AUMCPEO  | AUMC\_.Extrap\_obs    |  27.6994849 |  27.6994849 |          0 |
|       5 | AUMCPEO  | AUMC\_.Extrap\_obs    |  43.2128382 |  43.2128382 |          0 |
|       6 | AUMCPEO  | AUMC\_.Extrap\_obs    |  33.0774222 |  33.0774222 |          0 |
|       1 | AUMCIFP  | AUMCINF\_pred         |   7.8403303 |   7.8403303 |          0 |
|       2 | AUMCIFP  | AUMCINF\_pred         |   9.2048546 |   9.2048546 |          0 |
|       3 | AUMCIFP  | AUMCINF\_pred         |   6.0178784 |   6.0178784 |          0 |
|       4 | AUMCIFP  | AUMCINF\_pred         |   5.2284721 |   5.2284721 |          0 |
|       5 | AUMCIFP  | AUMCINF\_pred         |   6.0750083 |   6.0750083 |          0 |
|       6 | AUMCIFP  | AUMCINF\_pred         |   7.3003307 |   7.3003307 |          0 |
|       1 | AUMCPEP  | AUMC\_.Extrap\_pred   |  57.9538845 |  57.9538845 |          0 |
|       2 | AUMCPEP  | AUMC\_.Extrap\_pred   |  30.3817147 |  30.3817147 |          0 |
|       3 | AUMCPEP  | AUMC\_.Extrap\_pred   |  16.3270188 |  16.3270188 |          0 |
|       4 | AUMCPEP  | AUMC\_.Extrap\_pred   |  15.8636553 |  15.8636553 |          0 |
|       5 | AUMCPEP  | AUMC\_.Extrap\_pred   |  38.6013327 |  38.6013327 |          0 |
|       6 | AUMCPEP  | AUMC\_.Extrap\_pred   |  23.5984312 |  23.5984312 |          0 |
|       1 | MRTEVLST | MRTlast               |   1.9173089 |   1.9173089 |          0 |
|       2 | MRTEVLST | MRTlast               |   2.2180490 |   2.2180490 |          0 |
|       3 | MRTEVLST | MRTlast               |   1.7473430 |   1.7473430 |          0 |
|       4 | MRTEVLST | MRTlast               |   1.7997556 |   1.7997556 |          0 |
|       5 | MRTEVLST | MRTlast               |   1.9414830 |   1.9414830 |          0 |
|       6 | MRTEVLST | MRTlast               |   1.9630240 |   1.9630240 |          0 |
|       1 | MRTEVIFO | MRTINF\_obs           |   3.8413591 |   3.8413591 |          0 |
|       2 | MRTEVIFO | MRTINF\_obs           |   2.9808604 |   2.9808604 |          0 |
|       3 | MRTEVIFO | MRTINF\_obs           |   2.2797159 |   2.2797159 |          0 |
|       4 | MRTEVIFO | MRTINF\_obs           |   2.3335200 |   2.3335200 |          0 |
|       5 | MRTEVIFO | MRTINF\_obs           |   3.0428858 |   3.0428858 |          0 |
|       6 | MRTEVIFO | MRTINF\_obs           |   2.6920666 |   2.6920666 |          0 |
|       1 | MRTEVIFP | MRTINF\_pred          |   3.8494318 |   3.8494318 |          0 |
|       2 | MRTEVIFP | MRTINF\_pred          |   2.9348000 |   2.9348000 |          0 |
|       3 | MRTEVIFP | MRTINF\_pred          |   2.0218260 |   2.0218260 |          0 |
|       4 | MRTEVIFP | MRTINF\_pred          |   2.0710640 |   2.0710640 |          0 |
|       5 | MRTEVIFP | MRTINF\_pred          |   2.8691853 |   2.8691853 |          0 |
|       6 | MRTEVIFP | MRTINF\_pred          |   2.4331140 |   2.4331140 |          0 |

Indometh (n=6), Log, Extravascular

# Session Information

``` r
devtools::session_info()
```

    ##  setting  value                       
    ##  version  R version 3.5.0 (2018-04-23)
    ##  system   x86_64, mingw32             
    ##  ui       RTerm                       
    ##  language (EN)                        
    ##  collate  Korean_Korea.949            
    ##  tz       Asia/Seoul                  
    ##  date     2018-05-11                  
    ## 
    ##  package     * version date       source        
    ##  assertthat    0.2.0   2017-04-11 CRAN (R 3.4.0)
    ##  backports     1.1.2   2017-12-13 CRAN (R 3.5.0)
    ##  base        * 3.5.0   2018-04-23 local         
    ##  bindr         0.1.1   2018-03-13 CRAN (R 3.4.3)
    ##  bindrcpp    * 0.2.2   2018-03-29 CRAN (R 3.5.0)
    ##  colorspace    1.3-2   2016-12-14 CRAN (R 3.5.0)
    ##  compiler      3.5.0   2018-04-23 local         
    ##  datasets    * 3.5.0   2018-04-23 local         
    ##  devtools      1.13.5  2018-02-18 CRAN (R 3.5.0)
    ##  digest        0.6.15  2018-01-28 CRAN (R 3.5.0)
    ##  dplyr       * 0.7.4   2017-09-28 CRAN (R 3.5.0)
    ##  evaluate      0.10.1  2017-06-24 CRAN (R 3.4.1)
    ##  glue          1.2.0   2017-10-29 CRAN (R 3.5.0)
    ##  graphics    * 3.5.0   2018-04-23 local         
    ##  grDevices   * 3.5.0   2018-04-23 local         
    ##  highr         0.6     2016-05-09 CRAN (R 3.4.0)
    ##  hms           0.4.2   2018-03-10 CRAN (R 3.4.3)
    ##  htmltools     0.3.6   2017-04-28 CRAN (R 3.5.0)
    ##  httr          1.3.1   2017-08-20 CRAN (R 3.4.1)
    ##  kableExtra  * 0.8.0   2018-04-05 CRAN (R 3.4.4)
    ##  knitr       * 1.20    2018-02-20 CRAN (R 3.4.3)
    ##  magrittr      1.5     2014-11-22 CRAN (R 3.4.0)
    ##  memoise       1.1.0   2017-04-21 CRAN (R 3.4.0)
    ##  methods     * 3.5.0   2018-04-23 local         
    ##  munsell       0.4.3   2016-02-13 CRAN (R 3.4.0)
    ##  NonCompart  * 0.4.2   2018-05-11 CRAN (R 3.4.4)
    ##  pillar        1.2.2   2018-04-26 CRAN (R 3.5.0)
    ##  pkgconfig     2.0.1   2017-03-21 CRAN (R 3.4.0)
    ##  plyr          1.8.4   2016-06-08 CRAN (R 3.5.0)
    ##  purrr         0.2.4   2017-10-18 CRAN (R 3.5.0)
    ##  R6            2.2.2   2017-06-17 CRAN (R 3.4.1)
    ##  Rcpp          0.12.16 2018-03-13 CRAN (R 3.5.0)
    ##  readr         1.1.1   2017-05-16 CRAN (R 3.5.0)
    ##  rlang         0.2.0   2018-02-20 CRAN (R 3.5.0)
    ##  rmarkdown     1.9     2018-03-01 CRAN (R 3.4.3)
    ##  rprojroot     1.3-2   2018-01-03 CRAN (R 3.4.3)
    ##  rstudioapi    0.7     2017-09-07 CRAN (R 3.4.1)
    ##  rvest         0.3.2   2016-06-17 CRAN (R 3.4.0)
    ##  scales        0.5.0   2017-08-24 CRAN (R 3.5.0)
    ##  stats       * 3.5.0   2018-04-23 local         
    ##  stringi       1.2.2   2018-05-02 CRAN (R 3.5.0)
    ##  stringr       1.3.1   2018-05-10 CRAN (R 3.5.0)
    ##  tibble        1.4.2   2018-01-22 CRAN (R 3.5.0)
    ##  tidyr       * 0.8.0   2018-01-29 CRAN (R 3.5.0)
    ##  tidyselect    0.2.4   2018-02-26 CRAN (R 3.5.0)
    ##  tools         3.5.0   2018-04-23 local         
    ##  utils       * 3.5.0   2018-04-23 local         
    ##  viridisLite   0.3.0   2018-02-01 CRAN (R 3.4.3)
    ##  withr         2.1.2   2018-03-15 CRAN (R 3.4.4)
    ##  xml2          1.2.0   2018-01-24 CRAN (R 3.5.0)
    ##  yaml          2.1.19  2018-05-01 CRAN (R 3.5.0)

# References

<div id="refs" class="references">

<div id="ref-R-NonCompart">

Bae, Kyun-Seop. 2018. *NonCompart: Noncompartmental Analysis for
Pharmacokinetic Data*. <https://CRAN.R-project.org/package=NonCompart>.

</div>

<div id="ref-wnl">

Certara USA, Inc. 2018. *Phoenix Winnonlin: The Industry Standard for
Nca and Pk/Pd Modeling and Simulation*.
<https://www.certara.com/software/pkpd-modeling-and-simulation/phoenix-winnonlin/>.

</div>

<div id="ref-Kim_2018">

Kim, Hyungsub, Sungpil Han, Yong-Soon Cho, Seok-Kyu Yoon, and Kyun-Seop
Bae. 2018. “Development of R Packages: ‘NonCompart’ and ‘Ncar’ for
Noncompartmental Analysis (Nca).” *Translational and Clinical
Pharmacology* 26 (1). Korean Society for Clinical Pharmacology;
Therapeutics (KAMJE):10–15.
<http://www.tcpharm.org/Data/Journal/2/358.pdf>.

</div>

</div>
