Validation of Noncompartmental Analysis Performed by NonCompart R
package
================
Sungpil Han <shan@acp.kr>
2018-08-07



-----

# Introduction

NonCompart R package (Bae 2018; Kim et al. 2018) can conduct a
noncompartmental analysis as similar as possible to the most widely used
commercial software for pharmacokinetic analysis, i.e.
[Phoenix<sup>®</sup>
WinNonlin<sup>®</sup>](https://www.certara.com/software/pkpd-modeling-and-simulation/phoenix-winnonlin/)
(Certara USA 2018). This document provides validation of
noncompartmental analysis performed by NonCompart R package version
0.4.4 as compared to the results from the commercial software,
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

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Theoph (n=12), Linear,
Extravascular

</caption>

<thead>

<tr>

<th style="border-bottom:hidden" colspan="1">

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Pharmacokinetic
Parameters

</div>

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Values

</div>

</th>

<th style="border-bottom:hidden" colspan="1">

</th>

</tr>

<tr>

<th style="text-align:right;">

Subject

</th>

<th style="text-align:left;">

PPTESTCD

</th>

<th style="text-align:left;">

WNL

</th>

<th style="text-align:right;">

NonCompart

</th>

<th style="text-align:right;">

WinNonlin

</th>

<th style="text-align:right;">

Difference

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9999997

</td>

<td style="text-align:right;">

0.9999997

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9971954

</td>

<td style="text-align:right;">

0.9971954

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9993250

</td>

<td style="text-align:right;">

0.9993250

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9989241

</td>

<td style="text-align:right;">

0.9989241

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9986472

</td>

<td style="text-align:right;">

0.9986472

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9982413

</td>

<td style="text-align:right;">

0.9982413

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9986702

</td>

<td style="text-align:right;">

0.9986702

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9910124

</td>

<td style="text-align:right;">

0.9910124

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9994437

</td>

<td style="text-align:right;">

0.9994437

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9995087

</td>

<td style="text-align:right;">

0.9995087

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9999983

</td>

<td style="text-align:right;">

0.9999983

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9993968

</td>

<td style="text-align:right;">

0.9993968

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9999995

</td>

<td style="text-align:right;">

0.9999995

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9957931

</td>

<td style="text-align:right;">

0.9957931

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9986499

</td>

<td style="text-align:right;">

0.9986499

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9978483

</td>

<td style="text-align:right;">

0.9978483

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9979708

</td>

<td style="text-align:right;">

0.9979708

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9978896

</td>

<td style="text-align:right;">

0.9978896

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9980053

</td>

<td style="text-align:right;">

0.9980053

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9887655

</td>

<td style="text-align:right;">

0.9887655

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9988873

</td>

<td style="text-align:right;">

0.9988873

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9990174

</td>

<td style="text-align:right;">

0.9990174

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9999965

</td>

<td style="text-align:right;">

0.9999965

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9987936

</td>

<td style="text-align:right;">

0.9987936

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9999999

</td>

<td style="text-align:right;">

\-0.9999999

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9985967

</td>

<td style="text-align:right;">

\-0.9985967

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9996624

</td>

<td style="text-align:right;">

\-0.9996624

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9994619

</td>

<td style="text-align:right;">

\-0.9994619

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9993234

</td>

<td style="text-align:right;">

\-0.9993234

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9991203

</td>

<td style="text-align:right;">

\-0.9991203

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9993349

</td>

<td style="text-align:right;">

\-0.9993349

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9954961

</td>

<td style="text-align:right;">

\-0.9954961

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9997218

</td>

<td style="text-align:right;">

\-0.9997218

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9997543

</td>

<td style="text-align:right;">

\-0.9997543

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9999991

</td>

<td style="text-align:right;">

\-0.9999991

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9996984

</td>

<td style="text-align:right;">

\-0.9996984

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

4.0000000

</td>

<td style="text-align:right;">

4.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

4.0000000

</td>

<td style="text-align:right;">

4.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

7.0000000

</td>

<td style="text-align:right;">

7.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

4.0000000

</td>

<td style="text-align:right;">

4.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

6.0000000

</td>

<td style="text-align:right;">

6.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0484570

</td>

<td style="text-align:right;">

0.0484570

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.1040864

</td>

<td style="text-align:right;">

0.1040864

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.1024443

</td>

<td style="text-align:right;">

0.1024443

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0992870

</td>

<td style="text-align:right;">

0.0992870

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0866189

</td>

<td style="text-align:right;">

0.0866189

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0877957

</td>

<td style="text-align:right;">

0.0877957

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0883365

</td>

<td style="text-align:right;">

0.0883365

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0814505

</td>

<td style="text-align:right;">

0.0814505

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0824586

</td>

<td style="text-align:right;">

0.0824586

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0749598

</td>

<td style="text-align:right;">

0.0749598

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0954586

</td>

<td style="text-align:right;">

0.0954586

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.1102595

</td>

<td style="text-align:right;">

0.1102595

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

9.0500000

</td>

<td style="text-align:right;">

9.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

7.0300000

</td>

<td style="text-align:right;">

7.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

9.0200000

</td>

<td style="text-align:right;">

9.0200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

7.0200000

</td>

<td style="text-align:right;">

7.0200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

6.9800000

</td>

<td style="text-align:right;">

6.9800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

3.5300000

</td>

<td style="text-align:right;">

3.5300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

8.8000000

</td>

<td style="text-align:right;">

8.8000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

9.3800000

</td>

<td style="text-align:right;">

9.3800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

9.0300000

</td>

<td style="text-align:right;">

9.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

9.0300000

</td>

<td style="text-align:right;">

9.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.3700000

</td>

<td style="text-align:right;">

24.3700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.3000000

</td>

<td style="text-align:right;">

24.3000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.1700000

</td>

<td style="text-align:right;">

24.1700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.6500000

</td>

<td style="text-align:right;">

24.6500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.3500000

</td>

<td style="text-align:right;">

24.3500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

23.8500000

</td>

<td style="text-align:right;">

23.8500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.2200000

</td>

<td style="text-align:right;">

24.2200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.1200000

</td>

<td style="text-align:right;">

24.1200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.4300000

</td>

<td style="text-align:right;">

24.4300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

23.7000000

</td>

<td style="text-align:right;">

23.7000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.0800000

</td>

<td style="text-align:right;">

24.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.1500000

</td>

<td style="text-align:right;">

24.1500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

14.3043776

</td>

<td style="text-align:right;">

14.3043776

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

6.6593416

</td>

<td style="text-align:right;">

6.6593416

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

6.7660874

</td>

<td style="text-align:right;">

6.7660874

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

6.9812467

</td>

<td style="text-align:right;">

6.9812467

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

8.0022640

</td>

<td style="text-align:right;">

8.0022640

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

7.8949979

</td>

<td style="text-align:right;">

7.8949979

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

7.8466683

</td>

<td style="text-align:right;">

7.8466683

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

8.5100379

</td>

<td style="text-align:right;">

8.5100379

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

8.4059988

</td>

<td style="text-align:right;">

8.4059988

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

9.2469158

</td>

<td style="text-align:right;">

9.2469158

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

7.2612365

</td>

<td style="text-align:right;">

7.2612365

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

6.2865082

</td>

<td style="text-align:right;">

6.2865082

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

1.1200000

</td>

<td style="text-align:right;">

1.1200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

1.9200000

</td>

<td style="text-align:right;">

1.9200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

1.0200000

</td>

<td style="text-align:right;">

1.0200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

1.0700000

</td>

<td style="text-align:right;">

1.0700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

1.1500000

</td>

<td style="text-align:right;">

1.1500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

3.4800000

</td>

<td style="text-align:right;">

3.4800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

2.0200000

</td>

<td style="text-align:right;">

2.0200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.6300000

</td>

<td style="text-align:right;">

0.6300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

3.5500000

</td>

<td style="text-align:right;">

3.5500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.9800000

</td>

<td style="text-align:right;">

0.9800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

3.5200000

</td>

<td style="text-align:right;">

3.5200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

10.5000000

</td>

<td style="text-align:right;">

10.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

8.3300000

</td>

<td style="text-align:right;">

8.3300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

8.2000000

</td>

<td style="text-align:right;">

8.2000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

8.6000000

</td>

<td style="text-align:right;">

8.6000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

11.4000000

</td>

<td style="text-align:right;">

11.4000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

6.4400000

</td>

<td style="text-align:right;">

6.4400000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

7.0900000

</td>

<td style="text-align:right;">

7.0900000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

7.5600000

</td>

<td style="text-align:right;">

7.5600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

9.0300000

</td>

<td style="text-align:right;">

9.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

10.2100000

</td>

<td style="text-align:right;">

10.2100000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

9.7500000

</td>

<td style="text-align:right;">

9.7500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0328125

</td>

<td style="text-align:right;">

0.0328125

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0260312

</td>

<td style="text-align:right;">

0.0260312

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0256250

</td>

<td style="text-align:right;">

0.0256250

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0268750

</td>

<td style="text-align:right;">

0.0268750

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0356250

</td>

<td style="text-align:right;">

0.0356250

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0201250

</td>

<td style="text-align:right;">

0.0201250

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0221562

</td>

<td style="text-align:right;">

0.0221562

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0236250

</td>

<td style="text-align:right;">

0.0236250

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0282188

</td>

<td style="text-align:right;">

0.0282188

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0319063

</td>

<td style="text-align:right;">

0.0319062

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0250000

</td>

<td style="text-align:right;">

0.0250000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0304688

</td>

<td style="text-align:right;">

0.0304688

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.3700000

</td>

<td style="text-align:right;">

24.3700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.3000000

</td>

<td style="text-align:right;">

24.3000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.1700000

</td>

<td style="text-align:right;">

24.1700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.6500000

</td>

<td style="text-align:right;">

24.6500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.3500000

</td>

<td style="text-align:right;">

24.3500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

23.8500000

</td>

<td style="text-align:right;">

23.8500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.2200000

</td>

<td style="text-align:right;">

24.2200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.1200000

</td>

<td style="text-align:right;">

24.1200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.4300000

</td>

<td style="text-align:right;">

24.4300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

23.7000000

</td>

<td style="text-align:right;">

23.7000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.0800000

</td>

<td style="text-align:right;">

24.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.1500000

</td>

<td style="text-align:right;">

24.1500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

3.2800000

</td>

<td style="text-align:right;">

3.2800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.9000000

</td>

<td style="text-align:right;">

0.9000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

1.0500000

</td>

<td style="text-align:right;">

1.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

1.1500000

</td>

<td style="text-align:right;">

1.1500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

1.5700000

</td>

<td style="text-align:right;">

1.5700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.9200000

</td>

<td style="text-align:right;">

0.9200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

1.1500000

</td>

<td style="text-align:right;">

1.1500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

1.2500000

</td>

<td style="text-align:right;">

1.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

1.1200000

</td>

<td style="text-align:right;">

1.1200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

2.4200000

</td>

<td style="text-align:right;">

2.4200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.8600000

</td>

<td style="text-align:right;">

0.8600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

1.1700000

</td>

<td style="text-align:right;">

1.1700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

148.9230500

</td>

<td style="text-align:right;">

148.9230500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

91.5268000

</td>

<td style="text-align:right;">

91.5268000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

99.2865000

</td>

<td style="text-align:right;">

99.2865000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

106.7963000

</td>

<td style="text-align:right;">

106.7963000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

121.2944000

</td>

<td style="text-align:right;">

121.2944000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

73.7755500

</td>

<td style="text-align:right;">

73.7755500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

90.7534000

</td>

<td style="text-align:right;">

90.7534000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

88.5599500

</td>

<td style="text-align:right;">

88.5599500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

86.3261500

</td>

<td style="text-align:right;">

86.3261500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

138.3681000

</td>

<td style="text-align:right;">

138.3681000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

80.0936000

</td>

<td style="text-align:right;">

80.0936000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

119.9775000

</td>

<td style="text-align:right;">

119.9775000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

148.9230500

</td>

<td style="text-align:right;">

148.9230500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

91.5268000

</td>

<td style="text-align:right;">

91.5268000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

99.2865000

</td>

<td style="text-align:right;">

99.2865000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

106.7963000

</td>

<td style="text-align:right;">

106.7963000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

121.2944000

</td>

<td style="text-align:right;">

121.2944000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

73.7755500

</td>

<td style="text-align:right;">

73.7755500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

90.7534000

</td>

<td style="text-align:right;">

90.7534000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

88.5599500

</td>

<td style="text-align:right;">

88.5599500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

86.3261500

</td>

<td style="text-align:right;">

86.3261500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

138.3681000

</td>

<td style="text-align:right;">

138.3681000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

80.0936000

</td>

<td style="text-align:right;">

80.0936000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

119.9775000

</td>

<td style="text-align:right;">

119.9775000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

216.6119330

</td>

<td style="text-align:right;">

216.6119330

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

100.1734591

</td>

<td style="text-align:right;">

100.1734591

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

109.5359707

</td>

<td style="text-align:right;">

109.5359707

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

118.3788814

</td>

<td style="text-align:right;">

118.3788814

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

139.4197778

</td>

<td style="text-align:right;">

139.4197778

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

84.2544183

</td>

<td style="text-align:right;">

84.2544183

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

103.7718018

</td>

<td style="text-align:right;">

103.7718018

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

103.9066868

</td>

<td style="text-align:right;">

103.9066868

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

99.9087179

</td>

<td style="text-align:right;">

99.9087179

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

170.6520606

</td>

<td style="text-align:right;">

170.6520606

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

89.1027449

</td>

<td style="text-align:right;">

89.1027449

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

130.5888316

</td>

<td style="text-align:right;">

130.5888316

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.6769123

</td>

<td style="text-align:right;">

0.6769123

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.3130421

</td>

<td style="text-align:right;">

0.3130421

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.3422999

</td>

<td style="text-align:right;">

0.3422999

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.3699340

</td>

<td style="text-align:right;">

0.3699340

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.4356868

</td>

<td style="text-align:right;">

0.4356868

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.2632951

</td>

<td style="text-align:right;">

0.2632951

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.3242869

</td>

<td style="text-align:right;">

0.3242869

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.3247084

</td>

<td style="text-align:right;">

0.3247084

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.3122147

</td>

<td style="text-align:right;">

0.3122147

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.5332877

</td>

<td style="text-align:right;">

0.5332877

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.2784461

</td>

<td style="text-align:right;">

0.2784461

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.4080901

</td>

<td style="text-align:right;">

0.4080901

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

31.2489169

</td>

<td style="text-align:right;">

31.2489169

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

8.6316867

</td>

<td style="text-align:right;">

8.6316867

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

9.3571734

</td>

<td style="text-align:right;">

9.3571734

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

9.7843309

</td>

<td style="text-align:right;">

9.7843309

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

13.0005786

</td>

<td style="text-align:right;">

13.0005786

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

12.4371737

</td>

<td style="text-align:right;">

12.4371737

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

12.5452209

</td>

<td style="text-align:right;">

12.5452209

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

14.7697297

</td>

<td style="text-align:right;">

14.7697297

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

13.5949777

</td>

<td style="text-align:right;">

13.5949777

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

18.9180022

</td>

<td style="text-align:right;">

18.9180022

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

10.1109623

</td>

<td style="text-align:right;">

10.1109623

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

8.1257573

</td>

<td style="text-align:right;">

8.1257573

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

30.4867482

</td>

<td style="text-align:right;">

30.4867482

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

30.6904416

</td>

<td style="text-align:right;">

30.6904416

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

28.5170999

</td>

<td style="text-align:right;">

28.5170999

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

27.2259641

</td>

<td style="text-align:right;">

27.2259641

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

26.4979947

</td>

<td style="text-align:right;">

26.4979946

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

43.2597345

</td>

<td style="text-align:right;">

43.2597345

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

34.9084408

</td>

<td style="text-align:right;">

34.9084408

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

37.8105081

</td>

<td style="text-align:right;">

37.8105081

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

38.8427934

</td>

<td style="text-align:right;">

38.8427934

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

25.0155401

</td>

<td style="text-align:right;">

25.0155401

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

37.6221852

</td>

<td style="text-align:right;">

37.6221852

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

22.2242936

</td>

<td style="text-align:right;">

22.2242936

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

1.4772963

</td>

<td style="text-align:right;">

1.4772963

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

3.1944589

</td>

<td style="text-align:right;">

3.1944589

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

2.9214147

</td>

<td style="text-align:right;">

2.9214147

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

2.7031849

</td>

<td style="text-align:right;">

2.7031849

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

2.2952267

</td>

<td style="text-align:right;">

2.2952267

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

3.7980204

</td>

<td style="text-align:right;">

3.7980204

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

3.0836893

</td>

<td style="text-align:right;">

3.0836894

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

3.0796863

</td>

<td style="text-align:right;">

3.0796863

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

3.2029237

</td>

<td style="text-align:right;">

3.2029237

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

1.8751605

</td>

<td style="text-align:right;">

1.8751605

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

3.5913596

</td>

<td style="text-align:right;">

3.5913596

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

2.4504393

</td>

<td style="text-align:right;">

2.4504393

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

216.6149558

</td>

<td style="text-align:right;">

216.6149558

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

100.0643176

</td>

<td style="text-align:right;">

100.0643176

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

109.5857218

</td>

<td style="text-align:right;">

109.5857218

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

118.4435586

</td>

<td style="text-align:right;">

118.4435586

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

139.2546304

</td>

<td style="text-align:right;">

139.2546304

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

84.4966986

</td>

<td style="text-align:right;">

84.4966986

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

103.8931470

</td>

<td style="text-align:right;">

103.8931470

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

103.6430515

</td>

<td style="text-align:right;">

103.6430515

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

99.8660677

</td>

<td style="text-align:right;">

99.8660677

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

170.5679125

</td>

<td style="text-align:right;">

170.5679125

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

89.1007190

</td>

<td style="text-align:right;">

89.1007190

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

130.6390680

</td>

<td style="text-align:right;">

130.6390680

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.6769217

</td>

<td style="text-align:right;">

0.6769217

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.3127010

</td>

<td style="text-align:right;">

0.3127010

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.3424554

</td>

<td style="text-align:right;">

0.3424554

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.3701361

</td>

<td style="text-align:right;">

0.3701361

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.4351707

</td>

<td style="text-align:right;">

0.4351707

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.2640522

</td>

<td style="text-align:right;">

0.2640522

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.3246661

</td>

<td style="text-align:right;">

0.3246661

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.3238845

</td>

<td style="text-align:right;">

0.3238845

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.3120815

</td>

<td style="text-align:right;">

0.3120815

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.5330247

</td>

<td style="text-align:right;">

0.5330247

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.2784397

</td>

<td style="text-align:right;">

0.2784397

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.4082471

</td>

<td style="text-align:right;">

0.4082471

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

31.2498763

</td>

<td style="text-align:right;">

31.2498763

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

8.5320300

</td>

<td style="text-align:right;">

8.5320300

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

9.3983245

</td>

<td style="text-align:right;">

9.3983245

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

9.8335939

</td>

<td style="text-align:right;">

9.8335939

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

12.8974027

</td>

<td style="text-align:right;">

12.8974027

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

12.6882455

</td>

<td style="text-align:right;">

12.6882455

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

12.6473665

</td>

<td style="text-align:right;">

12.6473664

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

14.5529307

</td>

<td style="text-align:right;">

14.5529307

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

13.5580763

</td>

<td style="text-align:right;">

13.5580763

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

18.8780012

</td>

<td style="text-align:right;">

18.8780012

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

10.1089184

</td>

<td style="text-align:right;">

10.1089184

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

8.1610870

</td>

<td style="text-align:right;">

8.1610870

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

30.4863228

</td>

<td style="text-align:right;">

30.4863228

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

30.7239161

</td>

<td style="text-align:right;">

30.7239161

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

28.5041534

</td>

<td style="text-align:right;">

28.5041534

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

27.2110972

</td>

<td style="text-align:right;">

27.2110972

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

26.5294196

</td>

<td style="text-align:right;">

26.5294196

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

43.1356944

</td>

<td style="text-align:right;">

43.1356944

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

34.8676684

</td>

<td style="text-align:right;">

34.8676684

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

37.9066862

</td>

<td style="text-align:right;">

37.9066862

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

38.8593822

</td>

<td style="text-align:right;">

38.8593822

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

25.0278813

</td>

<td style="text-align:right;">

25.0278813

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

37.6230406

</td>

<td style="text-align:right;">

37.6230406

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

22.2157473

</td>

<td style="text-align:right;">

22.2157473

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

1.4772757

</td>

<td style="text-align:right;">

1.4772757

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

3.1979432

</td>

<td style="text-align:right;">

3.1979432

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

2.9200884

</td>

<td style="text-align:right;">

2.9200884

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

2.7017088

</td>

<td style="text-align:right;">

2.7017088

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

2.2979487

</td>

<td style="text-align:right;">

2.2979487

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

3.7871302

</td>

<td style="text-align:right;">

3.7871302

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

3.0800877

</td>

<td style="text-align:right;">

3.0800877

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

3.0875201

</td>

<td style="text-align:right;">

3.0875201

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

3.2042916

</td>

<td style="text-align:right;">

3.2042916

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

1.8760856

</td>

<td style="text-align:right;">

1.8760856

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

3.5914413

</td>

<td style="text-align:right;">

3.5914413

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

2.4494970

</td>

<td style="text-align:right;">

2.4494970

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

1459.0711035

</td>

<td style="text-align:right;">

1459.0711040

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

706.5865660

</td>

<td style="text-align:right;">

706.5865660

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

803.1858700

</td>

<td style="text-align:right;">

803.1858700

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

901.0842105

</td>

<td style="text-align:right;">

901.0842105

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

1017.1143165

</td>

<td style="text-align:right;">

1017.1143170

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

609.1523875

</td>

<td style="text-align:right;">

609.1523875

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

782.4198600

</td>

<td style="text-align:right;">

782.4198600

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

739.5345980

</td>

<td style="text-align:right;">

739.5345980

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

705.2296255

</td>

<td style="text-align:right;">

705.2296255

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

1278.1800420

</td>

<td style="text-align:right;">

1278.1800420

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

617.2422125

</td>

<td style="text-align:right;">

617.2422125

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

977.8807235

</td>

<td style="text-align:right;">

977.8807235

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

4505.5348194

</td>

<td style="text-align:right;">

4505.5348190

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

999.7722880

</td>

<td style="text-align:right;">

999.7722880

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

1150.9647687

</td>

<td style="text-align:right;">

1150.9647690

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

1303.2524014

</td>

<td style="text-align:right;">

1303.2524010

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

1667.7216119

</td>

<td style="text-align:right;">

1667.7216120

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

978.4284857

</td>

<td style="text-align:right;">

978.4284857

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

1245.0984083

</td>

<td style="text-align:right;">

1245.0984080

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

1298.1157547

</td>

<td style="text-align:right;">

1298.1157550

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

1201.7715381

</td>

<td style="text-align:right;">

1201.7715380

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

2473.9934274

</td>

<td style="text-align:right;">

2473.9934270

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

928.5599714

</td>

<td style="text-align:right;">

928.5599714

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

1330.3840024

</td>

<td style="text-align:right;">

1330.3840020

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

67.6160287

</td>

<td style="text-align:right;">

67.6160287

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

29.3252499

</td>

<td style="text-align:right;">

29.3252499

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

30.2162940

</td>

<td style="text-align:right;">

30.2162940

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

30.8588107

</td>

<td style="text-align:right;">

30.8588107

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

39.0117446

</td>

<td style="text-align:right;">

39.0117446

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

37.7417567

</td>

<td style="text-align:right;">

37.7417567

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

37.1599984

</td>

<td style="text-align:right;">

37.1599984

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

43.0301500

</td>

<td style="text-align:right;">

43.0301500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

41.3174965

</td>

<td style="text-align:right;">

41.3174965

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

48.3353501

</td>

<td style="text-align:right;">

48.3353501

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

33.5269416

</td>

<td style="text-align:right;">

33.5269415

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

26.4963558

</td>

<td style="text-align:right;">

26.4963558

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

4505.6708646

</td>

<td style="text-align:right;">

4505.6708650

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

996.0715835

</td>

<td style="text-align:right;">

996.0715835

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

1152.6528903

</td>

<td style="text-align:right;">

1152.6528900

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

1305.4981092

</td>

<td style="text-align:right;">

1305.4981090

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

1661.7936744

</td>

<td style="text-align:right;">

1661.7936740

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

986.9664597

</td>

<td style="text-align:right;">

986.9664597

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

1249.4110601

</td>

<td style="text-align:right;">

1249.4110600

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

1288.5201162

</td>

<td style="text-align:right;">

1288.5201160

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

1200.2123597

</td>

<td style="text-align:right;">

1200.2123600

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

2470.8765418

</td>

<td style="text-align:right;">

2470.8765420

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

928.4899636

</td>

<td style="text-align:right;">

928.4899636

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

1332.0528341

</td>

<td style="text-align:right;">

1332.0528340

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

67.6170065

</td>

<td style="text-align:right;">

67.6170065

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

29.0626720

</td>

<td style="text-align:right;">

29.0626720

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

30.3184960

</td>

<td style="text-align:right;">

30.3184960

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

30.9777468

</td>

<td style="text-align:right;">

30.9777468

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

38.7941877

</td>

<td style="text-align:right;">

38.7941877

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

38.2803355

</td>

<td style="text-align:right;">

38.2803355

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

37.3769062

</td>

<td style="text-align:right;">

37.3769062

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

42.6058943

</td>

<td style="text-align:right;">

42.6058943

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

41.2412629

</td>

<td style="text-align:right;">

41.2412629

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

48.2701778

</td>

<td style="text-align:right;">

48.2701778

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

33.5219295

</td>

<td style="text-align:right;">

33.5219295

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

26.5884432

</td>

<td style="text-align:right;">

26.5884432

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

9.7974834

</td>

<td style="text-align:right;">

9.7974834

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

7.7199964

</td>

<td style="text-align:right;">

7.7199964

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.0895778

</td>

<td style="text-align:right;">

8.0895778

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.4374104

</td>

<td style="text-align:right;">

8.4374104

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.3855010

</td>

<td style="text-align:right;">

8.3855010

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.2568329

</td>

<td style="text-align:right;">

8.2568329

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.6213834

</td>

<td style="text-align:right;">

8.6213834

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.3506664

</td>

<td style="text-align:right;">

8.3506664

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.1693626

</td>

<td style="text-align:right;">

8.1693627

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

9.2375341

</td>

<td style="text-align:right;">

9.2375341

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

7.7065110

</td>

<td style="text-align:right;">

7.7065110

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.1505343

</td>

<td style="text-align:right;">

8.1505343

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

20.8000305

</td>

<td style="text-align:right;">

20.8000305

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

9.9804109

</td>

<td style="text-align:right;">

9.9804109

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

10.5076420

</td>

<td style="text-align:right;">

10.5076420

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

11.0091630

</td>

<td style="text-align:right;">

11.0091630

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

11.9618725

</td>

<td style="text-align:right;">

11.9618725

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

11.6127855

</td>

<td style="text-align:right;">

11.6127855

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

11.9984272

</td>

<td style="text-align:right;">

11.9984272

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

12.4930916

</td>

<td style="text-align:right;">

12.4930916

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

12.0286954

</td>

<td style="text-align:right;">

12.0286954

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

14.4972959

</td>

<td style="text-align:right;">

14.4972959

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

10.4212275

</td>

<td style="text-align:right;">

10.4212274

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

10.1875787

</td>

<td style="text-align:right;">

10.1875787

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

20.8003683

</td>

<td style="text-align:right;">

20.8003683

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

9.9543135

</td>

<td style="text-align:right;">

9.9543135

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

10.5182762

</td>

<td style="text-align:right;">

10.5182762

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

11.0221115

</td>

<td style="text-align:right;">

11.0221115

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

11.9334895

</td>

<td style="text-align:right;">

11.9334895

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

11.6805328

</td>

<td style="text-align:right;">

11.6805328

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

12.0259237

</td>

<td style="text-align:right;">

12.0259237

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

12.4322866

</td>

<td style="text-align:right;">

12.4322866

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

12.0182199

</td>

<td style="text-align:right;">

12.0182199

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

14.4861745

</td>

<td style="text-align:right;">

14.4861745

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

10.4206787

</td>

<td style="text-align:right;">

10.4206787

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

10.1964355

</td>

<td style="text-align:right;">

10.1964355

</td>

<td style="text-align:right;">

0

</td>

</tr>

</tbody>

</table>

## Test 2: Theoph (n=12), Log, Extravascular

``` r
table_wres_rres(Wres2, Rres2,
                Caption = 'Theoph (n=12), Log, Extravascular')
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Theoph (n=12), Log,
Extravascular

</caption>

<thead>

<tr>

<th style="border-bottom:hidden" colspan="1">

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Pharmacokinetic
Parameters

</div>

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Values

</div>

</th>

<th style="border-bottom:hidden" colspan="1">

</th>

</tr>

<tr>

<th style="text-align:right;">

Subject

</th>

<th style="text-align:left;">

PPTESTCD

</th>

<th style="text-align:left;">

WNL

</th>

<th style="text-align:right;">

NonCompart

</th>

<th style="text-align:right;">

WinNonlin

</th>

<th style="text-align:right;">

Difference

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9999997

</td>

<td style="text-align:right;">

0.9999997

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9971954

</td>

<td style="text-align:right;">

0.9971954

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9993250

</td>

<td style="text-align:right;">

0.9993250

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9989241

</td>

<td style="text-align:right;">

0.9989241

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9986472

</td>

<td style="text-align:right;">

0.9986472

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9982413

</td>

<td style="text-align:right;">

0.9982413

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9986702

</td>

<td style="text-align:right;">

0.9986702

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9910124

</td>

<td style="text-align:right;">

0.9910124

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9994437

</td>

<td style="text-align:right;">

0.9994437

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9995087

</td>

<td style="text-align:right;">

0.9995087

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9999983

</td>

<td style="text-align:right;">

0.9999983

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9993968

</td>

<td style="text-align:right;">

0.9993968

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9999995

</td>

<td style="text-align:right;">

0.9999995

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9957931

</td>

<td style="text-align:right;">

0.9957931

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9986499

</td>

<td style="text-align:right;">

0.9986499

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9978483

</td>

<td style="text-align:right;">

0.9978483

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9979708

</td>

<td style="text-align:right;">

0.9979708

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9978896

</td>

<td style="text-align:right;">

0.9978896

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9980053

</td>

<td style="text-align:right;">

0.9980053

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9887655

</td>

<td style="text-align:right;">

0.9887655

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9988873

</td>

<td style="text-align:right;">

0.9988873

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9990174

</td>

<td style="text-align:right;">

0.9990174

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9999965

</td>

<td style="text-align:right;">

0.9999965

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9987936

</td>

<td style="text-align:right;">

0.9987936

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9999999

</td>

<td style="text-align:right;">

\-0.9999999

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9985967

</td>

<td style="text-align:right;">

\-0.9985967

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9996624

</td>

<td style="text-align:right;">

\-0.9996624

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9994619

</td>

<td style="text-align:right;">

\-0.9994619

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9993234

</td>

<td style="text-align:right;">

\-0.9993234

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9991203

</td>

<td style="text-align:right;">

\-0.9991203

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9993349

</td>

<td style="text-align:right;">

\-0.9993349

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9954961

</td>

<td style="text-align:right;">

\-0.9954961

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9997218

</td>

<td style="text-align:right;">

\-0.9997218

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9997543

</td>

<td style="text-align:right;">

\-0.9997543

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9999991

</td>

<td style="text-align:right;">

\-0.9999991

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9996984

</td>

<td style="text-align:right;">

\-0.9996984

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

4.0000000

</td>

<td style="text-align:right;">

4.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

4.0000000

</td>

<td style="text-align:right;">

4.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

7.0000000

</td>

<td style="text-align:right;">

7.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

4.0000000

</td>

<td style="text-align:right;">

4.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

6.0000000

</td>

<td style="text-align:right;">

6.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0484570

</td>

<td style="text-align:right;">

0.0484570

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.1040864

</td>

<td style="text-align:right;">

0.1040864

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.1024443

</td>

<td style="text-align:right;">

0.1024443

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0992870

</td>

<td style="text-align:right;">

0.0992870

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0866189

</td>

<td style="text-align:right;">

0.0866189

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0877957

</td>

<td style="text-align:right;">

0.0877957

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0883365

</td>

<td style="text-align:right;">

0.0883365

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0814505

</td>

<td style="text-align:right;">

0.0814505

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0824586

</td>

<td style="text-align:right;">

0.0824586

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0749598

</td>

<td style="text-align:right;">

0.0749598

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.0954586

</td>

<td style="text-align:right;">

0.0954586

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.1102595

</td>

<td style="text-align:right;">

0.1102595

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

9.0500000

</td>

<td style="text-align:right;">

9.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

7.0300000

</td>

<td style="text-align:right;">

7.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

9.0200000

</td>

<td style="text-align:right;">

9.0200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

7.0200000

</td>

<td style="text-align:right;">

7.0200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

6.9800000

</td>

<td style="text-align:right;">

6.9800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

3.5300000

</td>

<td style="text-align:right;">

3.5300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

8.8000000

</td>

<td style="text-align:right;">

8.8000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

9.3800000

</td>

<td style="text-align:right;">

9.3800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

9.0300000

</td>

<td style="text-align:right;">

9.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

9.0300000

</td>

<td style="text-align:right;">

9.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.3700000

</td>

<td style="text-align:right;">

24.3700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.3000000

</td>

<td style="text-align:right;">

24.3000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.1700000

</td>

<td style="text-align:right;">

24.1700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.6500000

</td>

<td style="text-align:right;">

24.6500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.3500000

</td>

<td style="text-align:right;">

24.3500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

23.8500000

</td>

<td style="text-align:right;">

23.8500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.2200000

</td>

<td style="text-align:right;">

24.2200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.1200000

</td>

<td style="text-align:right;">

24.1200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.4300000

</td>

<td style="text-align:right;">

24.4300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

23.7000000

</td>

<td style="text-align:right;">

23.7000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.0800000

</td>

<td style="text-align:right;">

24.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

24.1500000

</td>

<td style="text-align:right;">

24.1500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

14.3043776

</td>

<td style="text-align:right;">

14.3043776

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

6.6593416

</td>

<td style="text-align:right;">

6.6593416

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

6.7660874

</td>

<td style="text-align:right;">

6.7660874

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

6.9812467

</td>

<td style="text-align:right;">

6.9812467

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

8.0022640

</td>

<td style="text-align:right;">

8.0022640

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

7.8949979

</td>

<td style="text-align:right;">

7.8949979

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

7.8466683

</td>

<td style="text-align:right;">

7.8466683

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

8.5100379

</td>

<td style="text-align:right;">

8.5100379

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

8.4059988

</td>

<td style="text-align:right;">

8.4059988

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

9.2469158

</td>

<td style="text-align:right;">

9.2469158

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

7.2612365

</td>

<td style="text-align:right;">

7.2612365

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

6.2865082

</td>

<td style="text-align:right;">

6.2865082

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

1.1200000

</td>

<td style="text-align:right;">

1.1200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

1.9200000

</td>

<td style="text-align:right;">

1.9200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

1.0200000

</td>

<td style="text-align:right;">

1.0200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

1.0700000

</td>

<td style="text-align:right;">

1.0700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

1.1500000

</td>

<td style="text-align:right;">

1.1500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

3.4800000

</td>

<td style="text-align:right;">

3.4800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

2.0200000

</td>

<td style="text-align:right;">

2.0200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.6300000

</td>

<td style="text-align:right;">

0.6300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

3.5500000

</td>

<td style="text-align:right;">

3.5500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.9800000

</td>

<td style="text-align:right;">

0.9800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

3.5200000

</td>

<td style="text-align:right;">

3.5200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

10.5000000

</td>

<td style="text-align:right;">

10.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

8.3300000

</td>

<td style="text-align:right;">

8.3300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

8.2000000

</td>

<td style="text-align:right;">

8.2000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

8.6000000

</td>

<td style="text-align:right;">

8.6000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

11.4000000

</td>

<td style="text-align:right;">

11.4000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

6.4400000

</td>

<td style="text-align:right;">

6.4400000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

7.0900000

</td>

<td style="text-align:right;">

7.0900000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

7.5600000

</td>

<td style="text-align:right;">

7.5600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

9.0300000

</td>

<td style="text-align:right;">

9.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

10.2100000

</td>

<td style="text-align:right;">

10.2100000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

9.7500000

</td>

<td style="text-align:right;">

9.7500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0328125

</td>

<td style="text-align:right;">

0.0328125

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0260312

</td>

<td style="text-align:right;">

0.0260312

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0256250

</td>

<td style="text-align:right;">

0.0256250

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0268750

</td>

<td style="text-align:right;">

0.0268750

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0356250

</td>

<td style="text-align:right;">

0.0356250

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0201250

</td>

<td style="text-align:right;">

0.0201250

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0221562

</td>

<td style="text-align:right;">

0.0221562

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0236250

</td>

<td style="text-align:right;">

0.0236250

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0282188

</td>

<td style="text-align:right;">

0.0282188

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0319063

</td>

<td style="text-align:right;">

0.0319062

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0250000

</td>

<td style="text-align:right;">

0.0250000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0304688

</td>

<td style="text-align:right;">

0.0304688

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.3700000

</td>

<td style="text-align:right;">

24.3700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.3000000

</td>

<td style="text-align:right;">

24.3000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.1700000

</td>

<td style="text-align:right;">

24.1700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.6500000

</td>

<td style="text-align:right;">

24.6500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.3500000

</td>

<td style="text-align:right;">

24.3500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

23.8500000

</td>

<td style="text-align:right;">

23.8500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.2200000

</td>

<td style="text-align:right;">

24.2200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.1200000

</td>

<td style="text-align:right;">

24.1200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.4300000

</td>

<td style="text-align:right;">

24.4300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

23.7000000

</td>

<td style="text-align:right;">

23.7000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.0800000

</td>

<td style="text-align:right;">

24.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

24.1500000

</td>

<td style="text-align:right;">

24.1500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

3.2800000

</td>

<td style="text-align:right;">

3.2800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.9000000

</td>

<td style="text-align:right;">

0.9000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

1.0500000

</td>

<td style="text-align:right;">

1.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

1.1500000

</td>

<td style="text-align:right;">

1.1500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

1.5700000

</td>

<td style="text-align:right;">

1.5700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.9200000

</td>

<td style="text-align:right;">

0.9200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

1.1500000

</td>

<td style="text-align:right;">

1.1500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

1.2500000

</td>

<td style="text-align:right;">

1.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

1.1200000

</td>

<td style="text-align:right;">

1.1200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

2.4200000

</td>

<td style="text-align:right;">

2.4200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.8600000

</td>

<td style="text-align:right;">

0.8600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

1.1700000

</td>

<td style="text-align:right;">

1.1700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

147.2347485

</td>

<td style="text-align:right;">

147.2347485

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

88.7312755

</td>

<td style="text-align:right;">

88.7312755

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

95.8781978

</td>

<td style="text-align:right;">

95.8781978

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

102.6336232

</td>

<td style="text-align:right;">

102.6336232

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

118.1793538

</td>

<td style="text-align:right;">

118.1793538

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

71.6970150

</td>

<td style="text-align:right;">

71.6970150

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

87.9692274

</td>

<td style="text-align:right;">

87.9692274

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

86.8065635

</td>

<td style="text-align:right;">

86.8065635

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

83.9374360

</td>

<td style="text-align:right;">

83.9374360

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

135.5760701

</td>

<td style="text-align:right;">

135.5760701

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

77.8934723

</td>

<td style="text-align:right;">

77.8934723

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

115.2202082

</td>

<td style="text-align:right;">

115.2202082

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

147.2347485

</td>

<td style="text-align:right;">

147.2347485

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

88.7312755

</td>

<td style="text-align:right;">

88.7312755

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

95.8781978

</td>

<td style="text-align:right;">

95.8781978

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

102.6336232

</td>

<td style="text-align:right;">

102.6336232

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

118.1793538

</td>

<td style="text-align:right;">

118.1793538

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

71.6970150

</td>

<td style="text-align:right;">

71.6970150

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

87.9692274

</td>

<td style="text-align:right;">

87.9692274

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

86.8065635

</td>

<td style="text-align:right;">

86.8065635

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

83.9374360

</td>

<td style="text-align:right;">

83.9374360

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

135.5760701

</td>

<td style="text-align:right;">

135.5760701

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

77.8934723

</td>

<td style="text-align:right;">

77.8934723

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

115.2202082

</td>

<td style="text-align:right;">

115.2202082

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

214.9236316

</td>

<td style="text-align:right;">

214.9236316

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

97.3779346

</td>

<td style="text-align:right;">

97.3779346

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

106.1276685

</td>

<td style="text-align:right;">

106.1276685

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

114.2162046

</td>

<td style="text-align:right;">

114.2162046

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

136.3047316

</td>

<td style="text-align:right;">

136.3047316

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

82.1758833

</td>

<td style="text-align:right;">

82.1758833

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

100.9876292

</td>

<td style="text-align:right;">

100.9876292

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

102.1533003

</td>

<td style="text-align:right;">

102.1533003

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

97.5200039

</td>

<td style="text-align:right;">

97.5200039

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

167.8600307

</td>

<td style="text-align:right;">

167.8600307

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

86.9026173

</td>

<td style="text-align:right;">

86.9026173

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

125.8315397

</td>

<td style="text-align:right;">

125.8315397

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.6716363

</td>

<td style="text-align:right;">

0.6716363

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.3043060

</td>

<td style="text-align:right;">

0.3043060

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.3316490

</td>

<td style="text-align:right;">

0.3316490

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.3569256

</td>

<td style="text-align:right;">

0.3569256

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.4259523

</td>

<td style="text-align:right;">

0.4259523

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.2567996

</td>

<td style="text-align:right;">

0.2567996

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.3155863

</td>

<td style="text-align:right;">

0.3155863

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.3192291

</td>

<td style="text-align:right;">

0.3192291

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.3047500

</td>

<td style="text-align:right;">

0.3047500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.5245626

</td>

<td style="text-align:right;">

0.5245626

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.2715707

</td>

<td style="text-align:right;">

0.2715707

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.3932236

</td>

<td style="text-align:right;">

0.3932236

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

31.4943883

</td>

<td style="text-align:right;">

31.4943883

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

8.8794850

</td>

<td style="text-align:right;">

8.8794850

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

9.6576801

</td>

<td style="text-align:right;">

9.6576801

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

10.1409266

</td>

<td style="text-align:right;">

10.1409266

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

13.2976879

</td>

<td style="text-align:right;">

13.2976879

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

12.7517562

</td>

<td style="text-align:right;">

12.7517562

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

12.8910857

</td>

<td style="text-align:right;">

12.8910857

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

15.0232413

</td>

<td style="text-align:right;">

15.0232413

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

13.9279813

</td>

<td style="text-align:right;">

13.9279813

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

19.2326669

</td>

<td style="text-align:right;">

19.2326669

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

10.3669431

</td>

<td style="text-align:right;">

10.3669432

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

8.4329665

</td>

<td style="text-align:right;">

8.4329665

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

30.7262325

</td>

<td style="text-align:right;">

30.7262325

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

31.5715024

</td>

<td style="text-align:right;">

31.5715024

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

29.4329299

</td>

<td style="text-align:right;">

29.4329299

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

28.2182304

</td>

<td style="text-align:right;">

28.2182304

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

27.1035678

</td>

<td style="text-align:right;">

27.1035677

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

44.3539348

</td>

<td style="text-align:right;">

44.3539348

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

35.8708471

</td>

<td style="text-align:right;">

35.8708471

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

38.4594978

</td>

<td style="text-align:right;">

38.4594978

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

39.7942323

</td>

<td style="text-align:right;">

39.7942323

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

25.4316257

</td>

<td style="text-align:right;">

25.4316257

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

38.5746722

</td>

<td style="text-align:right;">

38.5746722

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

23.0645237

</td>

<td style="text-align:right;">

23.0645237

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

1.4889010

</td>

<td style="text-align:right;">

1.4889010

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

3.2861654

</td>

<td style="text-align:right;">

3.2861654

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

3.0152363

</td>

<td style="text-align:right;">

3.0152363

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

2.8017040

</td>

<td style="text-align:right;">

2.8017040

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

2.3476808

</td>

<td style="text-align:right;">

2.3476808

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

3.8940865

</td>

<td style="text-align:right;">

3.8940865

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

3.1687049

</td>

<td style="text-align:right;">

3.1687049

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

3.1325469

</td>

<td style="text-align:right;">

3.1325469

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

3.2813780

</td>

<td style="text-align:right;">

3.2813780

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

1.9063502

</td>

<td style="text-align:right;">

1.9063502

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

3.6822827

</td>

<td style="text-align:right;">

3.6822827

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

2.5430826

</td>

<td style="text-align:right;">

2.5430826

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

214.9266543

</td>

<td style="text-align:right;">

214.9266543

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

97.2687931

</td>

<td style="text-align:right;">

97.2687931

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

106.1774195

</td>

<td style="text-align:right;">

106.1774195

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

114.2808818

</td>

<td style="text-align:right;">

114.2808818

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

136.1395842

</td>

<td style="text-align:right;">

136.1395842

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

82.4181636

</td>

<td style="text-align:right;">

82.4181636

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

101.1089745

</td>

<td style="text-align:right;">

101.1089745

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

101.8896649

</td>

<td style="text-align:right;">

101.8896649

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

97.4773537

</td>

<td style="text-align:right;">

97.4773537

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

167.7758826

</td>

<td style="text-align:right;">

167.7758826

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

86.9005913

</td>

<td style="text-align:right;">

86.9005913

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

125.8817762

</td>

<td style="text-align:right;">

125.8817762

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.6716458

</td>

<td style="text-align:right;">

0.6716458

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.3039650

</td>

<td style="text-align:right;">

0.3039650

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.3318044

</td>

<td style="text-align:right;">

0.3318044

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.3571278

</td>

<td style="text-align:right;">

0.3571278

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.4254362

</td>

<td style="text-align:right;">

0.4254362

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.2575568

</td>

<td style="text-align:right;">

0.2575568

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.3159655

</td>

<td style="text-align:right;">

0.3159655

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.3184052

</td>

<td style="text-align:right;">

0.3184052

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.3046167

</td>

<td style="text-align:right;">

0.3046167

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.5242996

</td>

<td style="text-align:right;">

0.5242996

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.2715643

</td>

<td style="text-align:right;">

0.2715643

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.3933806

</td>

<td style="text-align:right;">

0.3933806

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

31.4953518

</td>

<td style="text-align:right;">

31.4953518

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

8.7772423

</td>

<td style="text-align:right;">

8.7772423

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

9.7000114

</td>

<td style="text-align:right;">

9.7000114

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

10.1917822

</td>

<td style="text-align:right;">

10.1917822

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

13.1925116

</td>

<td style="text-align:right;">

13.1925116

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

13.0082352

</td>

<td style="text-align:right;">

13.0082352

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

12.9956288

</td>

<td style="text-align:right;">

12.9956288

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

14.8033674

</td>

<td style="text-align:right;">

14.8033674

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

13.8903213

</td>

<td style="text-align:right;">

13.8903213

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

19.1921580

</td>

<td style="text-align:right;">

19.1921580

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

10.3648535

</td>

<td style="text-align:right;">

10.3648535

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

8.4695087

</td>

<td style="text-align:right;">

8.4695087

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

30.7258003

</td>

<td style="text-align:right;">

30.7258003

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

31.6069275

</td>

<td style="text-align:right;">

31.6069275

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

29.4191386

</td>

<td style="text-align:right;">

29.4191386

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

28.2022603

</td>

<td style="text-align:right;">

28.2022603

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

27.1364464

</td>

<td style="text-align:right;">

27.1364464

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

44.2235499

</td>

<td style="text-align:right;">

44.2235499

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

35.8277969

</td>

<td style="text-align:right;">

35.8277969

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

38.5590101

</td>

<td style="text-align:right;">

38.5590101

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

39.8116439

</td>

<td style="text-align:right;">

39.8116439

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

25.4443810

</td>

<td style="text-align:right;">

25.4443809

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

38.5755715

</td>

<td style="text-align:right;">

38.5755715

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

23.0553192

</td>

<td style="text-align:right;">

23.0553192

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

1.4888800

</td>

<td style="text-align:right;">

1.4888800

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

3.2898527

</td>

<td style="text-align:right;">

3.2898527

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

3.0138235

</td>

<td style="text-align:right;">

3.0138235

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

2.8001184

</td>

<td style="text-align:right;">

2.8001184

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

2.3505287

</td>

<td style="text-align:right;">

2.3505287

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

3.8826393

</td>

<td style="text-align:right;">

3.8826393

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

3.1649020

</td>

<td style="text-align:right;">

3.1649020

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

3.1406522

</td>

<td style="text-align:right;">

3.1406522

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

3.2828138

</td>

<td style="text-align:right;">

3.2828138

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

1.9073063

</td>

<td style="text-align:right;">

1.9073063

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

3.6823685

</td>

<td style="text-align:right;">

3.6823685

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

2.5420677

</td>

<td style="text-align:right;">

2.5420677

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

1499.1290852

</td>

<td style="text-align:right;">

1499.1290850

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

716.2787279

</td>

<td style="text-align:right;">

716.2787279

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

810.8726830

</td>

<td style="text-align:right;">

810.8726830

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

911.7828093

</td>

<td style="text-align:right;">

911.7828093

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

1038.8799844

</td>

<td style="text-align:right;">

1038.8799840

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

618.6659191

</td>

<td style="text-align:right;">

618.6659191

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

795.6267785

</td>

<td style="text-align:right;">

795.6267785

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

756.3619816

</td>

<td style="text-align:right;">

756.3619816

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

723.3794155

</td>

<td style="text-align:right;">

723.3794155

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

1306.7406149

</td>

<td style="text-align:right;">

1306.7406150

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

626.6357849

</td>

<td style="text-align:right;">

626.6357849

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

982.6343023

</td>

<td style="text-align:right;">

982.6343023

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

4545.5928011

</td>

<td style="text-align:right;">

4545.5928010

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

1009.4644499

</td>

<td style="text-align:right;">

1009.4644500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

1158.6515817

</td>

<td style="text-align:right;">

1158.6515820

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

1313.9510002

</td>

<td style="text-align:right;">

1313.9510000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

1689.4872798

</td>

<td style="text-align:right;">

1689.4872800

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

987.9420173

</td>

<td style="text-align:right;">

987.9420173

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

1258.3053268

</td>

<td style="text-align:right;">

1258.3053270

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

1314.9431383

</td>

<td style="text-align:right;">

1314.9431380

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

1219.9213281

</td>

<td style="text-align:right;">

1219.9213280

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

2502.5540002

</td>

<td style="text-align:right;">

2502.5540000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

937.9535438

</td>

<td style="text-align:right;">

937.9535438

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

1335.1375811

</td>

<td style="text-align:right;">

1335.1375810

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

67.0201632

</td>

<td style="text-align:right;">

67.0201632

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

29.0436897

</td>

<td style="text-align:right;">

29.0436897

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

30.0158308

</td>

<td style="text-align:right;">

30.0158308

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

30.6075486

</td>

<td style="text-align:right;">

30.6075486

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

38.5091562

</td>

<td style="text-align:right;">

38.5091562

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

37.3783169

</td>

<td style="text-align:right;">

37.3783169

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

36.7699745

</td>

<td style="text-align:right;">

36.7699745

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

42.4794913

</td>

<td style="text-align:right;">

42.4794914

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

40.7027815

</td>

<td style="text-align:right;">

40.7027815

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

47.7837196

</td>

<td style="text-align:right;">

47.7837196

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

33.1911704

</td>

<td style="text-align:right;">

33.1911704

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

26.4020191

</td>

<td style="text-align:right;">

26.4020191

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

4545.7288462

</td>

<td style="text-align:right;">

4545.7288460

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

1005.7637454

</td>

<td style="text-align:right;">

1005.7637450

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

1160.3397033

</td>

<td style="text-align:right;">

1160.3397030

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

1316.1967080

</td>

<td style="text-align:right;">

1316.1967080

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

1683.5593423

</td>

<td style="text-align:right;">

1683.5593420

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

996.4799913

</td>

<td style="text-align:right;">

996.4799913

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

1262.6179786

</td>

<td style="text-align:right;">

1262.6179790

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

1305.3474998

</td>

<td style="text-align:right;">

1305.3475000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

1218.3621498

</td>

<td style="text-align:right;">

1218.3621500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

2499.4371146

</td>

<td style="text-align:right;">

2499.4371150

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

937.8835360

</td>

<td style="text-align:right;">

937.8835360

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

1336.8064129

</td>

<td style="text-align:right;">

1336.8064130

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

67.0211503

</td>

<td style="text-align:right;">

67.0211503

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

28.7826061

</td>

<td style="text-align:right;">

28.7826061

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

30.1176474

</td>

<td style="text-align:right;">

30.1176474

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

30.7259467

</td>

<td style="text-align:right;">

30.7259467

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

38.2926424

</td>

<td style="text-align:right;">

38.2926424

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

37.9148679

</td>

<td style="text-align:right;">

37.9148679

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

36.9859457

</td>

<td style="text-align:right;">

36.9859457

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

42.0566568

</td>

<td style="text-align:right;">

42.0566568

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

40.6268969

</td>

<td style="text-align:right;">

40.6268969

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

47.7186040

</td>

<td style="text-align:right;">

47.7186040

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

33.1861835

</td>

<td style="text-align:right;">

33.1861835

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

26.4938967

</td>

<td style="text-align:right;">

26.4938967

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

10.1818973

</td>

<td style="text-align:right;">

10.1818973

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.0724494

</td>

<td style="text-align:right;">

8.0724494

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.4573209

</td>

<td style="text-align:right;">

8.4573209

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.8838607

</td>

<td style="text-align:right;">

8.8838607

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.7907063

</td>

<td style="text-align:right;">

8.7907063

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.6288937

</td>

<td style="text-align:right;">

8.6288937

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

9.0443761

</td>

<td style="text-align:right;">

9.0443761

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.7131889

</td>

<td style="text-align:right;">

8.7131889

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.6180785

</td>

<td style="text-align:right;">

8.6180785

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

9.6384311

</td>

<td style="text-align:right;">

9.6384311

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.0447792

</td>

<td style="text-align:right;">

8.0447792

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

8.5283156

</td>

<td style="text-align:right;">

8.5283156

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

21.1498046

</td>

<td style="text-align:right;">

21.1498045

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

10.3664599

</td>

<td style="text-align:right;">

10.3664599

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

10.9175260

</td>

<td style="text-align:right;">

10.9175260

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

11.5040681

</td>

<td style="text-align:right;">

11.5040681

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

12.3949276

</td>

<td style="text-align:right;">

12.3949276

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

12.0222866

</td>

<td style="text-align:right;">

12.0222866

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

12.4599947

</td>

<td style="text-align:right;">

12.4599947

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

12.8722531

</td>

<td style="text-align:right;">

12.8722531

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

12.5094471

</td>

<td style="text-align:right;">

12.5094471

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

14.9085758

</td>

<td style="text-align:right;">

14.9085758

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

10.7931564

</td>

<td style="text-align:right;">

10.7931564

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

10.6105161

</td>

<td style="text-align:right;">

10.6105161

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

21.1501401

</td>

<td style="text-align:right;">

21.1501401

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

10.3400455

</td>

<td style="text-align:right;">

10.3400455

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

10.9283095

</td>

<td style="text-align:right;">

10.9283095

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

11.5172082

</td>

<td style="text-align:right;">

11.5172082

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

12.3664205

</td>

<td style="text-align:right;">

12.3664205

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

12.0905386

</td>

<td style="text-align:right;">

12.0905386

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

7

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

12.4876944

</td>

<td style="text-align:right;">

12.4876944

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

12.8113828

</td>

<td style="text-align:right;">

12.8113828

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

9

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

12.4989252

</td>

<td style="text-align:right;">

12.4989252

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

10

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

14.8974756

</td>

<td style="text-align:right;">

14.8974756

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

11

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

10.7926025

</td>

<td style="text-align:right;">

10.7926025

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

12

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

10.6195388

</td>

<td style="text-align:right;">

10.6195388

</td>

<td style="text-align:right;">

0

</td>

</tr>

</tbody>

</table>

## Test 3: Indometh (n=6), Linear, IV Bolus

``` r
table_wres_rres(Wres3, Rres3,
                Caption = 'Indometh (n=6), Linear, IV Bolus')
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Indometh (n=6), Linear, IV
Bolus

</caption>

<thead>

<tr>

<th style="border-bottom:hidden" colspan="1">

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Pharmacokinetic
Parameters

</div>

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Values

</div>

</th>

<th style="border-bottom:hidden" colspan="1">

</th>

</tr>

<tr>

<th style="text-align:right;">

Subject

</th>

<th style="text-align:left;">

PPTESTCD

</th>

<th style="text-align:left;">

WNL

</th>

<th style="text-align:right;">

NonCompart

</th>

<th style="text-align:right;">

WinNonlin

</th>

<th style="text-align:right;">

Difference

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9970667

</td>

<td style="text-align:right;">

0.9970667

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9476691

</td>

<td style="text-align:right;">

0.9476691

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8758261

</td>

<td style="text-align:right;">

0.8758261

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8728249

</td>

<td style="text-align:right;">

0.8728249

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8752442

</td>

<td style="text-align:right;">

0.8752442

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9039538

</td>

<td style="text-align:right;">

0.9039538

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9941335

</td>

<td style="text-align:right;">

0.9941335

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9401933

</td>

<td style="text-align:right;">

0.9401933

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8603043

</td>

<td style="text-align:right;">

0.8603043

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8586943

</td>

<td style="text-align:right;">

0.8586943

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8544516

</td>

<td style="text-align:right;">

0.8544516

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8902329

</td>

<td style="text-align:right;">

0.8902329

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9985323

</td>

<td style="text-align:right;">

\-0.9985323

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9734830

</td>

<td style="text-align:right;">

\-0.9734830

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9358558

</td>

<td style="text-align:right;">

\-0.9358558

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9342510

</td>

<td style="text-align:right;">

\-0.9342510

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9355449

</td>

<td style="text-align:right;">

\-0.9355449

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9507649

</td>

<td style="text-align:right;">

\-0.9507649

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

11.0000000

</td>

<td style="text-align:right;">

11.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.1583205

</td>

<td style="text-align:right;">

0.1583205

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.3022800

</td>

<td style="text-align:right;">

0.3022800

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.4218926

</td>

<td style="text-align:right;">

0.4218926

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.4554455

</td>

<td style="text-align:right;">

0.4554455

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.2527478

</td>

<td style="text-align:right;">

0.2527478

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.3535205

</td>

<td style="text-align:right;">

0.3535205

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

5.0000000

</td>

<td style="text-align:right;">

5.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

4.3781270

</td>

<td style="text-align:right;">

4.3781270

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

2.2930632

</td>

<td style="text-align:right;">

2.2930632

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.6429468

</td>

<td style="text-align:right;">

1.6429468

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.5219104

</td>

<td style="text-align:right;">

1.5219104

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

2.7424461

</td>

<td style="text-align:right;">

2.7424461

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.9606986

</td>

<td style="text-align:right;">

1.9606986

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

1.5000000

</td>

<td style="text-align:right;">

1.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.7200000

</td>

<td style="text-align:right;">

2.7200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

1.8500000

</td>

<td style="text-align:right;">

1.8500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.0500000

</td>

<td style="text-align:right;">

2.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.3100000

</td>

<td style="text-align:right;">

2.3100000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0812000

</td>

<td style="text-align:right;">

0.0812000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.1088000

</td>

<td style="text-align:right;">

0.1088000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0740000

</td>

<td style="text-align:right;">

0.0740000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0820000

</td>

<td style="text-align:right;">

0.0820000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0924000

</td>

<td style="text-align:right;">

0.0924000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:right;">

2.3936170

</td>

<td style="text-align:right;">

2.3936170

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:right;">

2.5281595

</td>

<td style="text-align:right;">

2.5281595

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:right;">

4.9653691

</td>

<td style="text-align:right;">

4.9653691

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:right;">

2.4622302

</td>

<td style="text-align:right;">

2.4622302

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:right;">

4.0408654

</td>

<td style="text-align:right;">

4.0408654

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:right;">

3.7056250

</td>

<td style="text-align:right;">

3.7056250

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0500000

</td>

<td style="text-align:right;">

0.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0700000

</td>

<td style="text-align:right;">

0.0700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0900000

</td>

<td style="text-align:right;">

0.0900000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.0404521

</td>

<td style="text-align:right;">

2.0404521

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

3.2485199

</td>

<td style="text-align:right;">

3.2485199

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

3.5544211

</td>

<td style="text-align:right;">

3.5544211

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.7852788

</td>

<td style="text-align:right;">

2.7852788

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.4588582

</td>

<td style="text-align:right;">

2.4588582

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

3.3357031

</td>

<td style="text-align:right;">

3.3357031

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.0404521

</td>

<td style="text-align:right;">

2.0404521

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

3.2485199

</td>

<td style="text-align:right;">

3.2485199

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

3.5544211

</td>

<td style="text-align:right;">

3.5544211

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.7852788

</td>

<td style="text-align:right;">

2.7852788

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.4588582

</td>

<td style="text-align:right;">

2.4588582

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

3.3357031

</td>

<td style="text-align:right;">

3.3357031

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.3562672

</td>

<td style="text-align:right;">

2.3562672

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.5131752

</td>

<td style="text-align:right;">

3.5131752

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.7440428

</td>

<td style="text-align:right;">

3.7440428

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.9389745

</td>

<td style="text-align:right;">

2.9389745

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.6962490

</td>

<td style="text-align:right;">

2.6962490

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.5902852

</td>

<td style="text-align:right;">

3.5902852

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.0942507

</td>

<td style="text-align:right;">

0.0942507

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1405270

</td>

<td style="text-align:right;">

0.1405270

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1497617

</td>

<td style="text-align:right;">

0.1497617

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1175590

</td>

<td style="text-align:right;">

0.1175590

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1078500

</td>

<td style="text-align:right;">

0.1078500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1436114

</td>

<td style="text-align:right;">

0.1436114

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

13.4031956

</td>

<td style="text-align:right;">

13.4031956

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

7.5332215

</td>

<td style="text-align:right;">

7.5332215

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

5.0646241

</td>

<td style="text-align:right;">

5.0646241

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

5.2295685

</td>

<td style="text-align:right;">

5.2295685

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

8.8044838

</td>

<td style="text-align:right;">

8.8044838

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

7.0908603

</td>

<td style="text-align:right;">

7.0908603

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPBEO

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_obs

</td>

<td style="text-align:right;">

20.6556421

</td>

<td style="text-align:right;">

20.6556421

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPBEO

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_obs

</td>

<td style="text-align:right;">

16.2180906

</td>

<td style="text-align:right;">

16.2180906

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPBEO

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_obs

</td>

<td style="text-align:right;">

25.6586578

</td>

<td style="text-align:right;">

25.6586578

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPBEO

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_obs

</td>

<td style="text-align:right;">

18.3407098

</td>

<td style="text-align:right;">

18.3407098

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPBEO

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_obs

</td>

<td style="text-align:right;">

28.2376805

</td>

<td style="text-align:right;">

28.2376805

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPBEO

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_obs

</td>

<td style="text-align:right;">

20.9441054

</td>

<td style="text-align:right;">

20.9441054

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

67.0159780

</td>

<td style="text-align:right;">

67.0159780

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

23.5413171

</td>

<td style="text-align:right;">

23.5413171

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

15.8269504

</td>

<td style="text-align:right;">

15.8269504

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

18.6770303

</td>

<td style="text-align:right;">

18.6770303

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

36.6853493

</td>

<td style="text-align:right;">

36.6853493

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

19.6968341

</td>

<td style="text-align:right;">

19.6968341

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

10.6100020

</td>

<td style="text-align:right;">

10.6100020

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

7.1160698

</td>

<td style="text-align:right;">

7.1160698

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

6.6772740

</td>

<td style="text-align:right;">

6.6772740

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

8.5063686

</td>

<td style="text-align:right;">

8.5063686

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

9.2721407

</td>

<td style="text-align:right;">

9.2721407

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

6.9632351

</td>

<td style="text-align:right;">

6.9632351

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.3578369

</td>

<td style="text-align:right;">

2.3578369

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.4958268

</td>

<td style="text-align:right;">

3.4958268

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.6491670

</td>

<td style="text-align:right;">

3.6491670

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.8554521

</td>

<td style="text-align:right;">

2.8554521

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.6549884

</td>

<td style="text-align:right;">

2.6549884

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.4947956

</td>

<td style="text-align:right;">

3.4947956

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.0943135

</td>

<td style="text-align:right;">

0.0943135

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1398331

</td>

<td style="text-align:right;">

0.1398331

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1459667

</td>

<td style="text-align:right;">

0.1459667

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1142181

</td>

<td style="text-align:right;">

0.1142181

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1061995

</td>

<td style="text-align:right;">

0.1061995

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1397918

</td>

<td style="text-align:right;">

0.1397918

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

13.4608442

</td>

<td style="text-align:right;">

13.4608442

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

7.0743442

</td>

<td style="text-align:right;">

7.0743442

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

2.5963692

</td>

<td style="text-align:right;">

2.5963692

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

2.4575198

</td>

<td style="text-align:right;">

2.4575198

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

7.3872362

</td>

<td style="text-align:right;">

7.3872362

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

4.5522694

</td>

<td style="text-align:right;">

4.5522694

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPBEP

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_pred

</td>

<td style="text-align:right;">

20.6418914

</td>

<td style="text-align:right;">

20.6418914

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPBEP

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_pred

</td>

<td style="text-align:right;">

16.2985748

</td>

<td style="text-align:right;">

16.2985748

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPBEP

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_pred

</td>

<td style="text-align:right;">

26.3257654

</td>

<td style="text-align:right;">

26.3257654

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPBEP

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_pred

</td>

<td style="text-align:right;">

18.8771782

</td>

<td style="text-align:right;">

18.8771782

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPBEP

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_pred

</td>

<td style="text-align:right;">

28.6765156

</td>

<td style="text-align:right;">

28.6765156

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPBEP

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_pred

</td>

<td style="text-align:right;">

21.5163690

</td>

<td style="text-align:right;">

21.5163690

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

66.9713647

</td>

<td style="text-align:right;">

66.9713647

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

23.6581437

</td>

<td style="text-align:right;">

23.6581437

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

16.2384403

</td>

<td style="text-align:right;">

16.2384403

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

19.2233361

</td>

<td style="text-align:right;">

19.2233361

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

37.2554675

</td>

<td style="text-align:right;">

37.2554675

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

20.2350180

</td>

<td style="text-align:right;">

20.2350180

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

10.6029388

</td>

<td style="text-align:right;">

10.6029388

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

7.1513841

</td>

<td style="text-align:right;">

7.1513841

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

6.8508786

</td>

<td style="text-align:right;">

6.8508786

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

8.7551811

</td>

<td style="text-align:right;">

8.7551811

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

9.4162369

</td>

<td style="text-align:right;">

9.4162369

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

7.1534941

</td>

<td style="text-align:right;">

7.1534941

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

3.2712500

</td>

<td style="text-align:right;">

3.2712500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

6.3987500

</td>

<td style="text-align:right;">

6.3987500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

5.0062500

</td>

<td style="text-align:right;">

5.0062500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

4.3818750

</td>

<td style="text-align:right;">

4.3818750

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

3.7075000

</td>

<td style="text-align:right;">

3.7075000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

5.5325000

</td>

<td style="text-align:right;">

5.5325000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

7.7925545

</td>

<td style="text-align:right;">

7.7925545

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

9.3915223

</td>

<td style="text-align:right;">

9.3915223

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

6.9726784

</td>

<td style="text-align:right;">

6.9726784

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

5.9489028

</td>

<td style="text-align:right;">

5.9489028

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

6.5458663

</td>

<td style="text-align:right;">

6.5458663

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

8.2892908

</td>

<td style="text-align:right;">

8.2892908

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

58.0208261

</td>

<td style="text-align:right;">

58.0208261

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

31.8667432

</td>

<td style="text-align:right;">

31.8667432

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

28.2019090

</td>

<td style="text-align:right;">

28.2019090

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

26.3414589

</td>

<td style="text-align:right;">

26.3414589

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

43.3612023

</td>

<td style="text-align:right;">

43.3612023

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

33.2572574

</td>

<td style="text-align:right;">

33.2572574

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

7.8150259

</td>

<td style="text-align:right;">

7.8150259

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

9.1953427

</td>

<td style="text-align:right;">

9.1953427

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

5.9887901

</td>

<td style="text-align:right;">

5.9887901

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

5.0973376

</td>

<td style="text-align:right;">

5.0973376

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

6.0525342

</td>

<td style="text-align:right;">

6.0525342

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

7.2552635

</td>

<td style="text-align:right;">

7.2552635

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

58.1415337

</td>

<td style="text-align:right;">

58.1415337

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

30.4131425

</td>

<td style="text-align:right;">

30.4131425

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

16.4063210

</td>

<td style="text-align:right;">

16.4063210

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

14.0360055

</td>

<td style="text-align:right;">

14.0360055

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

38.7446663

</td>

<td style="text-align:right;">

38.7446663

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

23.7450164

</td>

<td style="text-align:right;">

23.7450164

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.6031986

</td>

<td style="text-align:right;">

1.6031986

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.9697432

</td>

<td style="text-align:right;">

1.9697432

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.4084572

</td>

<td style="text-align:right;">

1.4084572

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.5732267

</td>

<td style="text-align:right;">

1.5732267

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.5078137

</td>

<td style="text-align:right;">

1.5078137

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.6585709

</td>

<td style="text-align:right;">

1.6585709

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

3.3071607

</td>

<td style="text-align:right;">

3.3071607

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.6732291

</td>

<td style="text-align:right;">

2.6732291

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

1.8623394

</td>

<td style="text-align:right;">

1.8623394

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.0241424

</td>

<td style="text-align:right;">

2.0241424

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.4277678

</td>

<td style="text-align:right;">

2.4277678

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.3088112

</td>

<td style="text-align:right;">

2.3088112

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

3.3144897

</td>

<td style="text-align:right;">

3.3144897

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.6303771

</td>

<td style="text-align:right;">

2.6303771

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

1.6411390

</td>

<td style="text-align:right;">

1.6411390

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

1.7851245

</td>

<td style="text-align:right;">

1.7851245

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.2796838

</td>

<td style="text-align:right;">

2.2796838

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.0760194

</td>

<td style="text-align:right;">

2.0760194

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

35.0889819

</td>

<td style="text-align:right;">

35.0889819

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

19.0228851

</td>

<td style="text-align:right;">

19.0228851

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

12.4353504

</td>

<td style="text-align:right;">

12.4353504

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

17.2181012

</td>

<td style="text-align:right;">

17.2181012

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

22.5106044

</td>

<td style="text-align:right;">

22.5106044

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

16.0767951

</td>

<td style="text-align:right;">

16.0767951

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

35.1433309

</td>

<td style="text-align:right;">

35.1433309

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

18.8108371

</td>

<td style="text-align:right;">

18.8108371

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

11.2432438

</td>

<td style="text-align:right;">

11.2432438

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

15.6290886

</td>

<td style="text-align:right;">

15.6290886

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

21.4660427

</td>

<td style="text-align:right;">

21.4660427

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

14.8507925

</td>

<td style="text-align:right;">

14.8507925

</td>

<td style="text-align:right;">

0

</td>

</tr>

</tbody>

</table>

## Test 4: Indometh (n=6), Log, IV Bolus

``` r
table_wres_rres(Wres4, Rres4,
                Caption = 'Indometh (n=6), Log, IV Bolus')
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Indometh (n=6), Log, IV
Bolus

</caption>

<thead>

<tr>

<th style="border-bottom:hidden" colspan="1">

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Pharmacokinetic
Parameters

</div>

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Values

</div>

</th>

<th style="border-bottom:hidden" colspan="1">

</th>

</tr>

<tr>

<th style="text-align:right;">

Subject

</th>

<th style="text-align:left;">

PPTESTCD

</th>

<th style="text-align:left;">

WNL

</th>

<th style="text-align:right;">

NonCompart

</th>

<th style="text-align:right;">

WinNonlin

</th>

<th style="text-align:right;">

Difference

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9970667

</td>

<td style="text-align:right;">

0.9970667

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9476691

</td>

<td style="text-align:right;">

0.9476691

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8758261

</td>

<td style="text-align:right;">

0.8758261

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8728249

</td>

<td style="text-align:right;">

0.8728249

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8752442

</td>

<td style="text-align:right;">

0.8752442

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9039538

</td>

<td style="text-align:right;">

0.9039538

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9941335

</td>

<td style="text-align:right;">

0.9941335

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9401933

</td>

<td style="text-align:right;">

0.9401933

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8603043

</td>

<td style="text-align:right;">

0.8603043

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8586943

</td>

<td style="text-align:right;">

0.8586943

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8544516

</td>

<td style="text-align:right;">

0.8544516

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8902329

</td>

<td style="text-align:right;">

0.8902329

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9985323

</td>

<td style="text-align:right;">

\-0.9985323

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9734830

</td>

<td style="text-align:right;">

\-0.9734830

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9358558

</td>

<td style="text-align:right;">

\-0.9358558

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9342510

</td>

<td style="text-align:right;">

\-0.9342510

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9355449

</td>

<td style="text-align:right;">

\-0.9355449

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9507649

</td>

<td style="text-align:right;">

\-0.9507649

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

11.0000000

</td>

<td style="text-align:right;">

11.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.1583205

</td>

<td style="text-align:right;">

0.1583205

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.3022800

</td>

<td style="text-align:right;">

0.3022800

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.4218926

</td>

<td style="text-align:right;">

0.4218926

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.4554455

</td>

<td style="text-align:right;">

0.4554455

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.2527478

</td>

<td style="text-align:right;">

0.2527478

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.3535205

</td>

<td style="text-align:right;">

0.3535205

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

5.0000000

</td>

<td style="text-align:right;">

5.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

4.3781270

</td>

<td style="text-align:right;">

4.3781270

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

2.2930632

</td>

<td style="text-align:right;">

2.2930632

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.6429468

</td>

<td style="text-align:right;">

1.6429468

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.5219104

</td>

<td style="text-align:right;">

1.5219104

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

2.7424461

</td>

<td style="text-align:right;">

2.7424461

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.9606986

</td>

<td style="text-align:right;">

1.9606986

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

1.5000000

</td>

<td style="text-align:right;">

1.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.7200000

</td>

<td style="text-align:right;">

2.7200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

1.8500000

</td>

<td style="text-align:right;">

1.8500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.0500000

</td>

<td style="text-align:right;">

2.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.3100000

</td>

<td style="text-align:right;">

2.3100000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0812000

</td>

<td style="text-align:right;">

0.0812000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.1088000

</td>

<td style="text-align:right;">

0.1088000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0740000

</td>

<td style="text-align:right;">

0.0740000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0820000

</td>

<td style="text-align:right;">

0.0820000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0924000

</td>

<td style="text-align:right;">

0.0924000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:right;">

2.3936170

</td>

<td style="text-align:right;">

2.3936170

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:right;">

2.5281595

</td>

<td style="text-align:right;">

2.5281595

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:right;">

4.9653691

</td>

<td style="text-align:right;">

4.9653691

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:right;">

2.4622302

</td>

<td style="text-align:right;">

2.4622302

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:right;">

4.0408654

</td>

<td style="text-align:right;">

4.0408654

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:left;">

C0

</td>

<td style="text-align:right;">

3.7056250

</td>

<td style="text-align:right;">

3.7056250

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0500000

</td>

<td style="text-align:right;">

0.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0700000

</td>

<td style="text-align:right;">

0.0700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0900000

</td>

<td style="text-align:right;">

0.0900000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.0098984

</td>

<td style="text-align:right;">

2.0098984

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

3.2028878

</td>

<td style="text-align:right;">

3.2028878

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

3.4743971

</td>

<td style="text-align:right;">

3.4743971

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.7483832

</td>

<td style="text-align:right;">

2.7483832

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.3983736

</td>

<td style="text-align:right;">

2.3983736

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

3.2908266

</td>

<td style="text-align:right;">

3.2908266

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.0098984

</td>

<td style="text-align:right;">

2.0098984

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

3.2028878

</td>

<td style="text-align:right;">

3.2028878

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

3.4743971

</td>

<td style="text-align:right;">

3.4743971

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.7483832

</td>

<td style="text-align:right;">

2.7483832

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.3983736

</td>

<td style="text-align:right;">

2.3983736

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

3.2908266

</td>

<td style="text-align:right;">

3.2908266

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.3257135

</td>

<td style="text-align:right;">

2.3257135

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.4675431

</td>

<td style="text-align:right;">

3.4675431

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.6640188

</td>

<td style="text-align:right;">

3.6640188

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.9020789

</td>

<td style="text-align:right;">

2.9020789

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.6357645

</td>

<td style="text-align:right;">

2.6357645

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.5454087

</td>

<td style="text-align:right;">

3.5454087

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.0930285

</td>

<td style="text-align:right;">

0.0930285

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1387017

</td>

<td style="text-align:right;">

0.1387017

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1465608

</td>

<td style="text-align:right;">

0.1465608

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1160832

</td>

<td style="text-align:right;">

0.1160832

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1054306

</td>

<td style="text-align:right;">

0.1054306

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1418163

</td>

<td style="text-align:right;">

0.1418163

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

13.5792780

</td>

<td style="text-align:right;">

13.5792780

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

7.6323571

</td>

<td style="text-align:right;">

7.6323571

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

5.1752381

</td>

<td style="text-align:right;">

5.1752381

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

5.2960545

</td>

<td style="text-align:right;">

5.2960545

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

9.0065258

</td>

<td style="text-align:right;">

9.0065258

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

7.1806138

</td>

<td style="text-align:right;">

7.1806138

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPBEO

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_obs

</td>

<td style="text-align:right;">

20.5542573

</td>

<td style="text-align:right;">

20.5542573

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPBEO

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_obs

</td>

<td style="text-align:right;">

16.3658871

</td>

<td style="text-align:right;">

16.3658871

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPBEO

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_obs

</td>

<td style="text-align:right;">

25.4552663

</td>

<td style="text-align:right;">

25.4552663

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPBEO

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_obs

</td>

<td style="text-align:right;">

18.4484084

</td>

<td style="text-align:right;">

18.4484084

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPBEO

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_obs

</td>

<td style="text-align:right;">

27.8259014

</td>

<td style="text-align:right;">

27.8259014

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPBEO

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_obs

</td>

<td style="text-align:right;">

20.8230657

</td>

<td style="text-align:right;">

20.8230657

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

67.8963898

</td>

<td style="text-align:right;">

67.8963898

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

23.8511160

</td>

<td style="text-align:right;">

23.8511160

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

16.1726192

</td>

<td style="text-align:right;">

16.1726192

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

18.9144805

</td>

<td style="text-align:right;">

18.9144805

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

37.5271908

</td>

<td style="text-align:right;">

37.5271908

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

19.9461495

</td>

<td style="text-align:right;">

19.9461495

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

10.7493892

</td>

<td style="text-align:right;">

10.7493892

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

7.2097158

</td>

<td style="text-align:right;">

7.2097158

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

6.8231092

</td>

<td style="text-align:right;">

6.8231092

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

8.6145142

</td>

<td style="text-align:right;">

8.6145142

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

9.4849143

</td>

<td style="text-align:right;">

9.4849143

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

7.0513732

</td>

<td style="text-align:right;">

7.0513732

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.3272832

</td>

<td style="text-align:right;">

2.3272832

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.4501946

</td>

<td style="text-align:right;">

3.4501946

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.5691429

</td>

<td style="text-align:right;">

3.5691429

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.8185565

</td>

<td style="text-align:right;">

2.8185565

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.5945039

</td>

<td style="text-align:right;">

2.5945039

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.4499191

</td>

<td style="text-align:right;">

3.4499191

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.0930913

</td>

<td style="text-align:right;">

0.0930913

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1380078

</td>

<td style="text-align:right;">

0.1380078

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1427657

</td>

<td style="text-align:right;">

0.1427657

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1127423

</td>

<td style="text-align:right;">

0.1127423

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1037802

</td>

<td style="text-align:right;">

0.1037802

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1379968

</td>

<td style="text-align:right;">

0.1379968

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

13.6375646

</td>

<td style="text-align:right;">

13.6375646

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

7.1679092

</td>

<td style="text-align:right;">

7.1679092

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

2.6545826

</td>

<td style="text-align:right;">

2.6545826

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

2.4896893

</td>

<td style="text-align:right;">

2.4896893

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

7.5594516

</td>

<td style="text-align:right;">

7.5594516

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

4.6114853

</td>

<td style="text-align:right;">

4.6114853

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPBEP

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_pred

</td>

<td style="text-align:right;">

20.5403945

</td>

<td style="text-align:right;">

20.5403945

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPBEP

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_pred

</td>

<td style="text-align:right;">

16.4481790

</td>

<td style="text-align:right;">

16.4481790

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPBEP

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_pred

</td>

<td style="text-align:right;">

26.1319245

</td>

<td style="text-align:right;">

26.1319245

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPBEP

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_pred

</td>

<td style="text-align:right;">

18.9950907

</td>

<td style="text-align:right;">

18.9950907

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPBEP

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_pred

</td>

<td style="text-align:right;">

28.2684182

</td>

<td style="text-align:right;">

28.2684182

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPBEP

</td>

<td style="text-align:left;">

AUC\_.Back\_Ext\_pred

</td>

<td style="text-align:right;">

21.3994230

</td>

<td style="text-align:right;">

21.3994230

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

67.8505969

</td>

<td style="text-align:right;">

67.8505969

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

23.9710455

</td>

<td style="text-align:right;">

23.9710455

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

16.6025238

</td>

<td style="text-align:right;">

16.6025238

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

19.4749739

</td>

<td style="text-align:right;">

19.4749739

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

38.1239878

</td>

<td style="text-align:right;">

38.1239878

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

20.4982349

</td>

<td style="text-align:right;">

20.4982349

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

10.7421392

</td>

<td style="text-align:right;">

10.7421392

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

7.2459681

</td>

<td style="text-align:right;">

7.2459681

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

7.0044827

</td>

<td style="text-align:right;">

7.0044827

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

8.8697884

</td>

<td style="text-align:right;">

8.8697884

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

9.6357534

</td>

<td style="text-align:right;">

9.6357534

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

7.2465467

</td>

<td style="text-align:right;">

7.2465467

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

3.3047961

</td>

<td style="text-align:right;">

3.3047961

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

6.4131687

</td>

<td style="text-align:right;">

6.4131687

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

5.0552993

</td>

<td style="text-align:right;">

5.0552993

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

4.4049718

</td>

<td style="text-align:right;">

4.4049718

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

3.7472994

</td>

<td style="text-align:right;">

3.7472994

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

5.5904206

</td>

<td style="text-align:right;">

5.5904206

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

7.8261005

</td>

<td style="text-align:right;">

7.8261005

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

9.4059410

</td>

<td style="text-align:right;">

9.4059410

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

7.0217278

</td>

<td style="text-align:right;">

7.0217278

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

5.9719996

</td>

<td style="text-align:right;">

5.9719996

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

6.5856658

</td>

<td style="text-align:right;">

6.5856658

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

8.3472113

</td>

<td style="text-align:right;">

8.3472113

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

57.7721236

</td>

<td style="text-align:right;">

57.7721236

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

31.8178935

</td>

<td style="text-align:right;">

31.8178935

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

28.0049084

</td>

<td style="text-align:right;">

28.0049084

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

26.2395827

</td>

<td style="text-align:right;">

26.2395827

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

43.0991557

</td>

<td style="text-align:right;">

43.0991557

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

33.0264883

</td>

<td style="text-align:right;">

33.0264883

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

7.8485720

</td>

<td style="text-align:right;">

7.8485720

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

9.2097614

</td>

<td style="text-align:right;">

9.2097614

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

6.0378395

</td>

<td style="text-align:right;">

6.0378395

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

5.1204344

</td>

<td style="text-align:right;">

5.1204344

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

6.0923336

</td>

<td style="text-align:right;">

6.0923336

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

7.3131841

</td>

<td style="text-align:right;">

7.3131841

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

57.8930274

</td>

<td style="text-align:right;">

57.8930274

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

30.3655279

</td>

<td style="text-align:right;">

30.3655279

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

16.2730417

</td>

<td style="text-align:right;">

16.2730417

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

13.9726930

</td>

<td style="text-align:right;">

13.9726930

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

38.4915588

</td>

<td style="text-align:right;">

38.4915588

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

23.5569554

</td>

<td style="text-align:right;">

23.5569554

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.6442602

</td>

<td style="text-align:right;">

1.6442602

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

2.0023083

</td>

<td style="text-align:right;">

2.0023083

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.4550148

</td>

<td style="text-align:right;">

1.4550148

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.6027502

</td>

<td style="text-align:right;">

1.6027502

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.5624335

</td>

<td style="text-align:right;">

1.5624335

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.6987892

</td>

<td style="text-align:right;">

1.6987892

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

3.3650320

</td>

<td style="text-align:right;">

3.3650320

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.7125665

</td>

<td style="text-align:right;">

2.7125665

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

1.9164006

</td>

<td style="text-align:right;">

1.9164006

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.0578350

</td>

<td style="text-align:right;">

2.0578350

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.4985790

</td>

<td style="text-align:right;">

2.4985790

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.3543721

</td>

<td style="text-align:right;">

2.3543721

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

3.3724181

</td>

<td style="text-align:right;">

3.3724181

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.6693455

</td>

<td style="text-align:right;">

2.6693455

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

1.6916777

</td>

<td style="text-align:right;">

1.6916777

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

1.8166868

</td>

<td style="text-align:right;">

1.8166868

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.3481690

</td>

<td style="text-align:right;">

2.3481690

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.1198132

</td>

<td style="text-align:right;">

2.1198132

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

36.1720388

</td>

<td style="text-align:right;">

36.1720388

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

19.5568334

</td>

<td style="text-align:right;">

19.5568334

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

13.0758105

</td>

<td style="text-align:right;">

13.0758105

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

17.7272490

</td>

<td style="text-align:right;">

17.7272490

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

23.6988080

</td>

<td style="text-align:right;">

23.6988080

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

16.6015562

</td>

<td style="text-align:right;">

16.6015562

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

36.2269851

</td>

<td style="text-align:right;">

36.2269851

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

19.3419923

</td>

<td style="text-align:right;">

19.3419923

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

11.8493272

</td>

<td style="text-align:right;">

11.8493272

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

16.1136274

</td>

<td style="text-align:right;">

16.1136274

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

22.6263772

</td>

<td style="text-align:right;">

22.6263772

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

15.3613252

</td>

<td style="text-align:right;">

15.3613252

</td>

<td style="text-align:right;">

0

</td>

</tr>

</tbody>

</table>

## Test 5: Indometh (n=6), Linear, IV Infusion (0.25hr)

``` r
table_wres_rres(Wres5, Rres5,
                Caption = 'Indometh (n=6), Linear, IV Infusion (0.25hr)')
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Indometh (n=6), Linear, IV Infusion
(0.25hr)

</caption>

<thead>

<tr>

<th style="border-bottom:hidden" colspan="1">

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Pharmacokinetic
Parameters

</div>

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Values

</div>

</th>

<th style="border-bottom:hidden" colspan="1">

</th>

</tr>

<tr>

<th style="text-align:right;">

Subject

</th>

<th style="text-align:left;">

PPTESTCD

</th>

<th style="text-align:left;">

WNL

</th>

<th style="text-align:right;">

NonCompart

</th>

<th style="text-align:right;">

WinNonlin

</th>

<th style="text-align:right;">

Difference

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9970667

</td>

<td style="text-align:right;">

0.9970667

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9476691

</td>

<td style="text-align:right;">

0.9476691

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8758261

</td>

<td style="text-align:right;">

0.8758261

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8671179

</td>

<td style="text-align:right;">

0.8671179

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8752442

</td>

<td style="text-align:right;">

0.8752442

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9039538

</td>

<td style="text-align:right;">

0.9039538

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9941335

</td>

<td style="text-align:right;">

0.9941335

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9401933

</td>

<td style="text-align:right;">

0.9401933

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8603043

</td>

<td style="text-align:right;">

0.8603043

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8505077

</td>

<td style="text-align:right;">

0.8505077

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8544516

</td>

<td style="text-align:right;">

0.8544516

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8902329

</td>

<td style="text-align:right;">

0.8902329

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9985323

</td>

<td style="text-align:right;">

\-0.9985323

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9734830

</td>

<td style="text-align:right;">

\-0.9734830

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9358558

</td>

<td style="text-align:right;">

\-0.9358558

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9311917

</td>

<td style="text-align:right;">

\-0.9311917

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9355449

</td>

<td style="text-align:right;">

\-0.9355449

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9507649

</td>

<td style="text-align:right;">

\-0.9507649

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.1583205

</td>

<td style="text-align:right;">

0.1583205

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.3022800

</td>

<td style="text-align:right;">

0.3022800

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.4218926

</td>

<td style="text-align:right;">

0.4218926

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.4290762

</td>

<td style="text-align:right;">

0.4290762

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.2527478

</td>

<td style="text-align:right;">

0.2527478

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.3535205

</td>

<td style="text-align:right;">

0.3535205

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

5.0000000

</td>

<td style="text-align:right;">

5.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

4.3781270

</td>

<td style="text-align:right;">

4.3781270

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

2.2930632

</td>

<td style="text-align:right;">

2.2930632

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.6429468

</td>

<td style="text-align:right;">

1.6429468

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.6154409

</td>

<td style="text-align:right;">

1.6154409

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

2.7424461

</td>

<td style="text-align:right;">

2.7424461

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.9606986

</td>

<td style="text-align:right;">

1.9606986

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

1.5000000

</td>

<td style="text-align:right;">

1.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.7200000

</td>

<td style="text-align:right;">

2.7200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

1.8500000

</td>

<td style="text-align:right;">

1.8500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.0500000

</td>

<td style="text-align:right;">

2.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.3100000

</td>

<td style="text-align:right;">

2.3100000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0812000

</td>

<td style="text-align:right;">

0.0812000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.1088000

</td>

<td style="text-align:right;">

0.1088000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0740000

</td>

<td style="text-align:right;">

0.0740000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0820000

</td>

<td style="text-align:right;">

0.0820000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0924000

</td>

<td style="text-align:right;">

0.0924000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0500000

</td>

<td style="text-align:right;">

0.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0700000

</td>

<td style="text-align:right;">

0.0700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0900000

</td>

<td style="text-align:right;">

0.0900000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

1.7412500

</td>

<td style="text-align:right;">

1.7412500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.9325000

</td>

<td style="text-align:right;">

2.9325000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.9337500

</td>

<td style="text-align:right;">

2.9337500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.4775000

</td>

<td style="text-align:right;">

2.4775000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

1.9537500

</td>

<td style="text-align:right;">

1.9537500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.8725000

</td>

<td style="text-align:right;">

2.8725000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

1.7412500

</td>

<td style="text-align:right;">

1.7412500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.9325000

</td>

<td style="text-align:right;">

2.9325000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.9337500

</td>

<td style="text-align:right;">

2.9337500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.4775000

</td>

<td style="text-align:right;">

2.4775000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

1.9537500

</td>

<td style="text-align:right;">

1.9537500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.8725000

</td>

<td style="text-align:right;">

2.8725000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.0570651

</td>

<td style="text-align:right;">

2.0570651

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.1971553

</td>

<td style="text-align:right;">

3.1971553

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.1233717

</td>

<td style="text-align:right;">

3.1233717

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.6406412

</td>

<td style="text-align:right;">

2.6406412

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.1911408

</td>

<td style="text-align:right;">

2.1911408

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.1270821

</td>

<td style="text-align:right;">

3.1270821

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.0822826

</td>

<td style="text-align:right;">

0.0822826

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1278862

</td>

<td style="text-align:right;">

0.1278862

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1249349

</td>

<td style="text-align:right;">

0.1249349

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1056256

</td>

<td style="text-align:right;">

0.1056256

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.0876456

</td>

<td style="text-align:right;">

0.0876456

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1250833

</td>

<td style="text-align:right;">

0.1250833

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

15.3527035

</td>

<td style="text-align:right;">

15.3527035

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

8.2778360

</td>

<td style="text-align:right;">

8.2778360

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

6.0710577

</td>

<td style="text-align:right;">

6.0710577

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

6.1780905

</td>

<td style="text-align:right;">

6.1780905

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

10.8341191

</td>

<td style="text-align:right;">

10.8341191

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

8.1412032

</td>

<td style="text-align:right;">

8.1412032

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

76.7635175

</td>

<td style="text-align:right;">

76.7635175

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

25.8682374

</td>

<td style="text-align:right;">

25.8682374

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

18.9720552

</td>

<td style="text-align:right;">

18.9720552

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

22.0646091

</td>

<td style="text-align:right;">

22.0646091

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

45.1421631

</td>

<td style="text-align:right;">

45.1421631

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

22.6144534

</td>

<td style="text-align:right;">

22.6144534

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

12.1532371

</td>

<td style="text-align:right;">

12.1532371

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

7.8194513

</td>

<td style="text-align:right;">

7.8194513

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

8.0041706

</td>

<td style="text-align:right;">

8.0041706

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

9.4673975

</td>

<td style="text-align:right;">

9.4673975

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

11.4095817

</td>

<td style="text-align:right;">

11.4095817

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

7.9946733

</td>

<td style="text-align:right;">

7.9946733

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.0586347

</td>

<td style="text-align:right;">

2.0586347

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.1798068

</td>

<td style="text-align:right;">

3.1798068

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.0284958

</td>

<td style="text-align:right;">

3.0284958

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.5577884

</td>

<td style="text-align:right;">

2.5577884

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.1498803

</td>

<td style="text-align:right;">

2.1498803

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.0315925

</td>

<td style="text-align:right;">

3.0315925

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.0823454

</td>

<td style="text-align:right;">

0.0823454

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1271923

</td>

<td style="text-align:right;">

0.1271923

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1211398

</td>

<td style="text-align:right;">

0.1211398

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1023115

</td>

<td style="text-align:right;">

0.1023115

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.0859952

</td>

<td style="text-align:right;">

0.0859952

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1212637

</td>

<td style="text-align:right;">

0.1212637

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

15.4172443

</td>

<td style="text-align:right;">

15.4172443

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

7.7774164

</td>

<td style="text-align:right;">

7.7774164

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

3.1284787

</td>

<td style="text-align:right;">

3.1284787

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

3.1389787

</td>

<td style="text-align:right;">

3.1389787

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

9.1228460

</td>

<td style="text-align:right;">

9.1228460

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

5.2478198

</td>

<td style="text-align:right;">

5.2478198

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

76.7049878

</td>

<td style="text-align:right;">

76.7049878

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

26.0093699

</td>

<td style="text-align:right;">

26.0093699

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

19.5664063

</td>

<td style="text-align:right;">

19.5664063

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

22.7793336

</td>

<td style="text-align:right;">

22.7793336

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

46.0085322

</td>

<td style="text-align:right;">

46.0085322

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

23.3267671

</td>

<td style="text-align:right;">

23.3267671

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

12.1439707

</td>

<td style="text-align:right;">

12.1439707

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

7.8621128

</td>

<td style="text-align:right;">

7.8621128

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

8.2549230

</td>

<td style="text-align:right;">

8.2549230

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

9.7740688

</td>

<td style="text-align:right;">

9.7740688

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

11.6285546

</td>

<td style="text-align:right;">

11.6285546

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

8.2464909

</td>

<td style="text-align:right;">

8.2464909

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

3.2712500

</td>

<td style="text-align:right;">

3.2712500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

6.3987500

</td>

<td style="text-align:right;">

6.3987500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

5.0062500

</td>

<td style="text-align:right;">

5.0062500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

4.3818750

</td>

<td style="text-align:right;">

4.3818750

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

3.7075000

</td>

<td style="text-align:right;">

3.7075000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

5.5325000

</td>

<td style="text-align:right;">

5.5325000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

7.7925545

</td>

<td style="text-align:right;">

7.7925545

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

9.3915223

</td>

<td style="text-align:right;">

9.3915223

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

6.9726784

</td>

<td style="text-align:right;">

6.9726784

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

6.0672197

</td>

<td style="text-align:right;">

6.0672197

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

6.5458663

</td>

<td style="text-align:right;">

6.5458663

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

8.2892908

</td>

<td style="text-align:right;">

8.2892908

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

58.0208261

</td>

<td style="text-align:right;">

58.0208261

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

31.8667432

</td>

<td style="text-align:right;">

31.8667432

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

28.2019090

</td>

<td style="text-align:right;">

28.2019090

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

27.7778746

</td>

<td style="text-align:right;">

27.7778746

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

43.3612023

</td>

<td style="text-align:right;">

43.3612023

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

33.2572574

</td>

<td style="text-align:right;">

33.2572574

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

7.8150259

</td>

<td style="text-align:right;">

7.8150259

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

9.1953427

</td>

<td style="text-align:right;">

9.1953427

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

5.9887901

</td>

<td style="text-align:right;">

5.9887901

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

5.2113018

</td>

<td style="text-align:right;">

5.2113018

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

6.0525342

</td>

<td style="text-align:right;">

6.0525342

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

7.2552635

</td>

<td style="text-align:right;">

7.2552635

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

58.1415337

</td>

<td style="text-align:right;">

58.1415337

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

30.4131425

</td>

<td style="text-align:right;">

30.4131425

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

16.4063210

</td>

<td style="text-align:right;">

16.4063210

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

15.9159231

</td>

<td style="text-align:right;">

15.9159231

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

38.7446663

</td>

<td style="text-align:right;">

38.7446663

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

23.7450164

</td>

<td style="text-align:right;">

23.7450164

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.7536791

</td>

<td style="text-align:right;">

1.7536791

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

2.0570119

</td>

<td style="text-align:right;">

2.0570119

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.5814337

</td>

<td style="text-align:right;">

1.5814337

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.6436680

</td>

<td style="text-align:right;">

1.6436680

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.7726328

</td>

<td style="text-align:right;">

1.7726328

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.8010226

</td>

<td style="text-align:right;">

1.8010226

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

3.6631905

</td>

<td style="text-align:right;">

3.6631905

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.8124621

</td>

<td style="text-align:right;">

2.8124621

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.1074203

</td>

<td style="text-align:right;">

2.1074203

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.1726312

</td>

<td style="text-align:right;">

2.1726312

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.8624239

</td>

<td style="text-align:right;">

2.8624239

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.5258069

</td>

<td style="text-align:right;">

2.5258069

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

3.6712178

</td>

<td style="text-align:right;">

3.6712178

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.7667929

</td>

<td style="text-align:right;">

2.7667929

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

1.8524801

</td>

<td style="text-align:right;">

1.8524801

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

1.9124249

</td>

<td style="text-align:right;">

1.9124249

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.6902890

</td>

<td style="text-align:right;">

2.6902890

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.2682186

</td>

<td style="text-align:right;">

2.2682186

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

44.5196227

</td>

<td style="text-align:right;">

44.5196227

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

21.9919102

</td>

<td style="text-align:right;">

21.9919102

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

16.8681518

</td>

<td style="text-align:right;">

16.8681518

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

20.5691634

</td>

<td style="text-align:right;">

20.5691634

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

32.6590590

</td>

<td style="text-align:right;">

32.6590590

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

20.1930009

</td>

<td style="text-align:right;">

20.1930009

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

44.5831617

</td>

<td style="text-align:right;">

44.5831617

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

21.7528377

</td>

<td style="text-align:right;">

21.7528377

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

15.2920802

</td>

<td style="text-align:right;">

15.2920802

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

18.6921722

</td>

<td style="text-align:right;">

18.6921722

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

31.2841719

</td>

<td style="text-align:right;">

31.2841719

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

18.7048438

</td>

<td style="text-align:right;">

18.7048438

</td>

<td style="text-align:right;">

0

</td>

</tr>

</tbody>

</table>

## Test 6: Indometh (n=6), Log, IV Infusion (0.25hr)

``` r
table_wres_rres(Wres6, Rres6,
                Caption = 'Indometh (n=6), Log, IV Infusion (0.25hr)')
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Indometh (n=6), Log, IV Infusion
(0.25hr)

</caption>

<thead>

<tr>

<th style="border-bottom:hidden" colspan="1">

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Pharmacokinetic
Parameters

</div>

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Values

</div>

</th>

<th style="border-bottom:hidden" colspan="1">

</th>

</tr>

<tr>

<th style="text-align:right;">

Subject

</th>

<th style="text-align:left;">

PPTESTCD

</th>

<th style="text-align:left;">

WNL

</th>

<th style="text-align:right;">

NonCompart

</th>

<th style="text-align:right;">

WinNonlin

</th>

<th style="text-align:right;">

Difference

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9970667

</td>

<td style="text-align:right;">

0.9970667

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9476691

</td>

<td style="text-align:right;">

0.9476691

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8758261

</td>

<td style="text-align:right;">

0.8758261

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8671179

</td>

<td style="text-align:right;">

0.8671179

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8752442

</td>

<td style="text-align:right;">

0.8752442

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9039538

</td>

<td style="text-align:right;">

0.9039538

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9941335

</td>

<td style="text-align:right;">

0.9941335

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9401933

</td>

<td style="text-align:right;">

0.9401933

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8603043

</td>

<td style="text-align:right;">

0.8603043

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8505077

</td>

<td style="text-align:right;">

0.8505077

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8544516

</td>

<td style="text-align:right;">

0.8544516

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8902329

</td>

<td style="text-align:right;">

0.8902329

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9985323

</td>

<td style="text-align:right;">

\-0.9985323

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9734830

</td>

<td style="text-align:right;">

\-0.9734830

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9358558

</td>

<td style="text-align:right;">

\-0.9358558

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9311917

</td>

<td style="text-align:right;">

\-0.9311917

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9355449

</td>

<td style="text-align:right;">

\-0.9355449

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9507649

</td>

<td style="text-align:right;">

\-0.9507649

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.1583205

</td>

<td style="text-align:right;">

0.1583205

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.3022800

</td>

<td style="text-align:right;">

0.3022800

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.4218926

</td>

<td style="text-align:right;">

0.4218926

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.4290762

</td>

<td style="text-align:right;">

0.4290762

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.2527478

</td>

<td style="text-align:right;">

0.2527478

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.3535205

</td>

<td style="text-align:right;">

0.3535205

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

5.0000000

</td>

<td style="text-align:right;">

5.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

4.3781270

</td>

<td style="text-align:right;">

4.3781270

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

2.2930632

</td>

<td style="text-align:right;">

2.2930632

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.6429468

</td>

<td style="text-align:right;">

1.6429468

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.6154409

</td>

<td style="text-align:right;">

1.6154409

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

2.7424461

</td>

<td style="text-align:right;">

2.7424461

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.9606986

</td>

<td style="text-align:right;">

1.9606986

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

1.5000000

</td>

<td style="text-align:right;">

1.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.7200000

</td>

<td style="text-align:right;">

2.7200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

1.8500000

</td>

<td style="text-align:right;">

1.8500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.0500000

</td>

<td style="text-align:right;">

2.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.3100000

</td>

<td style="text-align:right;">

2.3100000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0812000

</td>

<td style="text-align:right;">

0.0812000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.1088000

</td>

<td style="text-align:right;">

0.1088000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0740000

</td>

<td style="text-align:right;">

0.0740000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0820000

</td>

<td style="text-align:right;">

0.0820000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0924000

</td>

<td style="text-align:right;">

0.0924000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0500000

</td>

<td style="text-align:right;">

0.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0700000

</td>

<td style="text-align:right;">

0.0700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0900000

</td>

<td style="text-align:right;">

0.0900000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

1.7193653

</td>

<td style="text-align:right;">

1.7193653

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.8891436

</td>

<td style="text-align:right;">

2.8891436

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.8817113

</td>

<td style="text-align:right;">

2.8817113

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.4442459

</td>

<td style="text-align:right;">

2.4442459

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

1.9211984

</td>

<td style="text-align:right;">

1.9211984

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.8413138

</td>

<td style="text-align:right;">

2.8413138

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

1.7193653

</td>

<td style="text-align:right;">

1.7193653

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.8891436

</td>

<td style="text-align:right;">

2.8891436

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.8817113

</td>

<td style="text-align:right;">

2.8817113

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.4442459

</td>

<td style="text-align:right;">

2.4442459

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

1.9211984

</td>

<td style="text-align:right;">

1.9211984

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.8413138

</td>

<td style="text-align:right;">

2.8413138

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.0351804

</td>

<td style="text-align:right;">

2.0351804

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.1537989

</td>

<td style="text-align:right;">

3.1537989

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.0713330

</td>

<td style="text-align:right;">

3.0713330

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.6073871

</td>

<td style="text-align:right;">

2.6073871

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.1585892

</td>

<td style="text-align:right;">

2.1585892

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.0958959

</td>

<td style="text-align:right;">

3.0958959

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.0814072

</td>

<td style="text-align:right;">

0.0814072

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1261520

</td>

<td style="text-align:right;">

0.1261520

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1228533

</td>

<td style="text-align:right;">

0.1228533

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1042955

</td>

<td style="text-align:right;">

0.1042955

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.0863436

</td>

<td style="text-align:right;">

0.0863436

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1238358

</td>

<td style="text-align:right;">

0.1238358

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

15.5177942

</td>

<td style="text-align:right;">

15.5177942

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

8.3916343

</td>

<td style="text-align:right;">

8.3916343

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

6.1739217

</td>

<td style="text-align:right;">

6.1739217

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

6.2568848

</td>

<td style="text-align:right;">

6.2568848

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

10.9974979

</td>

<td style="text-align:right;">

10.9974979

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

8.2232127

</td>

<td style="text-align:right;">

8.2232127

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

77.5889712

</td>

<td style="text-align:right;">

77.5889712

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

26.2238573

</td>

<td style="text-align:right;">

26.2238573

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

19.2935053

</td>

<td style="text-align:right;">

19.2935053

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

22.3460171

</td>

<td style="text-align:right;">

22.3460171

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

45.8229078

</td>

<td style="text-align:right;">

45.8229078

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZO

</td>

<td style="text-align:left;">

Vz\_obs

</td>

<td style="text-align:right;">

22.8422576

</td>

<td style="text-align:right;">

22.8422576

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

12.2839234

</td>

<td style="text-align:right;">

12.2839234

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

7.9269481

</td>

<td style="text-align:right;">

7.9269481

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

8.1397881

</td>

<td style="text-align:right;">

8.1397881

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

9.5881430

</td>

<td style="text-align:right;">

9.5881430

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

11.5816384

</td>

<td style="text-align:right;">

11.5816384

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLO

</td>

<td style="text-align:left;">

Cl\_obs

</td>

<td style="text-align:right;">

8.0752068

</td>

<td style="text-align:right;">

8.0752068

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.0367500

</td>

<td style="text-align:right;">

2.0367500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.1364504

</td>

<td style="text-align:right;">

3.1364504

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.9764572

</td>

<td style="text-align:right;">

2.9764572

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.5245343

</td>

<td style="text-align:right;">

2.5245343

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.1173287

</td>

<td style="text-align:right;">

2.1173287

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.0004063

</td>

<td style="text-align:right;">

3.0004063

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.0814700

</td>

<td style="text-align:right;">

0.0814700

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1254580

</td>

<td style="text-align:right;">

0.1254580

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1190583

</td>

<td style="text-align:right;">

0.1190583

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1009814

</td>

<td style="text-align:right;">

0.1009814

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.0846931

</td>

<td style="text-align:right;">

0.0846931

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1200163

</td>

<td style="text-align:right;">

0.1200163

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

15.5829013

</td>

<td style="text-align:right;">

15.5829013

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

7.8849267

</td>

<td style="text-align:right;">

7.8849267

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

3.1831752

</td>

<td style="text-align:right;">

3.1831752

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

3.1803265

</td>

<td style="text-align:right;">

3.1803265

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

9.2630996

</td>

<td style="text-align:right;">

9.2630996

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

5.3023656

</td>

<td style="text-align:right;">

5.3023656

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

77.5291765

</td>

<td style="text-align:right;">

77.5291765

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

26.3689077

</td>

<td style="text-align:right;">

26.3689077

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

19.9084941

</td>

<td style="text-align:right;">

19.9084941

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

23.0793917

</td>

<td style="text-align:right;">

23.0793917

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

46.7158621

</td>

<td style="text-align:right;">

46.7158621

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZP

</td>

<td style="text-align:left;">

Vz\_pred

</td>

<td style="text-align:right;">

23.5692251

</td>

<td style="text-align:right;">

23.5692251

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

12.2744566

</td>

<td style="text-align:right;">

12.2744566

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

7.9707940

</td>

<td style="text-align:right;">

7.9707940

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

8.3992473

</td>

<td style="text-align:right;">

8.3992473

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

9.9028165

</td>

<td style="text-align:right;">

9.9028165

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

11.8073306

</td>

<td style="text-align:right;">

11.8073306

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLP

</td>

<td style="text-align:left;">

Cl\_pred

</td>

<td style="text-align:right;">

8.3322048

</td>

<td style="text-align:right;">

8.3322048

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

3.2965543

</td>

<td style="text-align:right;">

3.2965543

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

6.4082620

</td>

<td style="text-align:right;">

6.4082620

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

5.0353383

</td>

<td style="text-align:right;">

5.0353383

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

4.3990453

</td>

<td style="text-align:right;">

4.3990453

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

3.7299741

</td>

<td style="text-align:right;">

3.7299741

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

5.5775672

</td>

<td style="text-align:right;">

5.5775672

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

7.8178588

</td>

<td style="text-align:right;">

7.8178588

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

9.4010343

</td>

<td style="text-align:right;">

9.4010343

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

7.0017667

</td>

<td style="text-align:right;">

7.0017667

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

6.0843899

</td>

<td style="text-align:right;">

6.0843899

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

6.5683405

</td>

<td style="text-align:right;">

6.5683405

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

8.3343579

</td>

<td style="text-align:right;">

8.3343579

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

57.8330281

</td>

<td style="text-align:right;">

57.8330281

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

31.8345005

</td>

<td style="text-align:right;">

31.8345005

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

28.0847466

</td>

<td style="text-align:right;">

28.0847466

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

27.6994849

</td>

<td style="text-align:right;">

27.6994849

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

43.2128382

</td>

<td style="text-align:right;">

43.2128382

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

33.0774222

</td>

<td style="text-align:right;">

33.0774222

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

7.8403303

</td>

<td style="text-align:right;">

7.8403303

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

9.2048546

</td>

<td style="text-align:right;">

9.2048546

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

6.0178784

</td>

<td style="text-align:right;">

6.0178784

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

5.2284721

</td>

<td style="text-align:right;">

5.2284721

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

6.0750083

</td>

<td style="text-align:right;">

6.0750083

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

7.3003307

</td>

<td style="text-align:right;">

7.3003307

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

57.9538845

</td>

<td style="text-align:right;">

57.9538845

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

30.3817147

</td>

<td style="text-align:right;">

30.3817147

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

16.3270188

</td>

<td style="text-align:right;">

16.3270188

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

15.8636553

</td>

<td style="text-align:right;">

15.8636553

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

38.6013327

</td>

<td style="text-align:right;">

38.6013327

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

23.5984312

</td>

<td style="text-align:right;">

23.5984312

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.7923089

</td>

<td style="text-align:right;">

1.7923089

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

2.0930490

</td>

<td style="text-align:right;">

2.0930490

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.6223430

</td>

<td style="text-align:right;">

1.6223430

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.6747556

</td>

<td style="text-align:right;">

1.6747556

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.8164830

</td>

<td style="text-align:right;">

1.8164830

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTIVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.8380240

</td>

<td style="text-align:right;">

1.8380240

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

3.7163591

</td>

<td style="text-align:right;">

3.7163591

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.8558604

</td>

<td style="text-align:right;">

2.8558604

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.1547159

</td>

<td style="text-align:right;">

2.1547159

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.2085200

</td>

<td style="text-align:right;">

2.2085200

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.9178858

</td>

<td style="text-align:right;">

2.9178858

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTIVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.5670666

</td>

<td style="text-align:right;">

2.5670666

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

3.7244318

</td>

<td style="text-align:right;">

3.7244318

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.8098000

</td>

<td style="text-align:right;">

2.8098000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

1.8968260

</td>

<td style="text-align:right;">

1.8968260

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

1.9460640

</td>

<td style="text-align:right;">

1.9460640

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.7441853

</td>

<td style="text-align:right;">

2.7441853

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTIVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.3081140

</td>

<td style="text-align:right;">

2.3081140

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

45.6514707

</td>

<td style="text-align:right;">

45.6514707

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

22.6382575

</td>

<td style="text-align:right;">

22.6382575

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

17.5389306

</td>

<td style="text-align:right;">

17.5389306

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

21.1756058

</td>

<td style="text-align:right;">

21.1756058

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

33.7938980

</td>

<td style="text-align:right;">

33.7938980

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VSSO

</td>

<td style="text-align:left;">

Vss\_obs

</td>

<td style="text-align:right;">

20.7295934

</td>

<td style="text-align:right;">

20.7295934

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

45.7153760

</td>

<td style="text-align:right;">

45.7153760

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

22.3963368

</td>

<td style="text-align:right;">

22.3963368

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

15.9319103

</td>

<td style="text-align:right;">

15.9319103

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

19.2715146

</td>

<td style="text-align:right;">

19.2715146

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

32.4015028

</td>

<td style="text-align:right;">

32.4015028

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VSSP

</td>

<td style="text-align:left;">

Vss\_pred

</td>

<td style="text-align:right;">

19.2316785

</td>

<td style="text-align:right;">

19.2316785

</td>

<td style="text-align:right;">

0

</td>

</tr>

</tbody>

</table>

## Test 7: Indometh (n=6), Linear, Extravascular

``` r
table_wres_rres(Wres7, Rres7,
                Caption = 'Indometh (n=6), Linear, Extravascular')
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Indometh (n=6), Linear,
Extravascular

</caption>

<thead>

<tr>

<th style="border-bottom:hidden" colspan="1">

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Pharmacokinetic
Parameters

</div>

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Values

</div>

</th>

<th style="border-bottom:hidden" colspan="1">

</th>

</tr>

<tr>

<th style="text-align:right;">

Subject

</th>

<th style="text-align:left;">

PPTESTCD

</th>

<th style="text-align:left;">

WNL

</th>

<th style="text-align:right;">

NonCompart

</th>

<th style="text-align:right;">

WinNonlin

</th>

<th style="text-align:right;">

Difference

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9970667

</td>

<td style="text-align:right;">

0.9970667

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9476691

</td>

<td style="text-align:right;">

0.9476691

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8758261

</td>

<td style="text-align:right;">

0.8758261

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8671179

</td>

<td style="text-align:right;">

0.8671179

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8752442

</td>

<td style="text-align:right;">

0.8752442

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9039538

</td>

<td style="text-align:right;">

0.9039538

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9941335

</td>

<td style="text-align:right;">

0.9941335

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9401933

</td>

<td style="text-align:right;">

0.9401933

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8603043

</td>

<td style="text-align:right;">

0.8603043

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8505077

</td>

<td style="text-align:right;">

0.8505077

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8544516

</td>

<td style="text-align:right;">

0.8544516

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8902329

</td>

<td style="text-align:right;">

0.8902329

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9985323

</td>

<td style="text-align:right;">

\-0.9985323

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9734830

</td>

<td style="text-align:right;">

\-0.9734830

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9358558

</td>

<td style="text-align:right;">

\-0.9358558

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9311917

</td>

<td style="text-align:right;">

\-0.9311917

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9355449

</td>

<td style="text-align:right;">

\-0.9355449

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9507649

</td>

<td style="text-align:right;">

\-0.9507649

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.1583205

</td>

<td style="text-align:right;">

0.1583205

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.3022800

</td>

<td style="text-align:right;">

0.3022800

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.4218926

</td>

<td style="text-align:right;">

0.4218926

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.4290762

</td>

<td style="text-align:right;">

0.4290762

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.2527478

</td>

<td style="text-align:right;">

0.2527478

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.3535205

</td>

<td style="text-align:right;">

0.3535205

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

5.0000000

</td>

<td style="text-align:right;">

5.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

4.3781270

</td>

<td style="text-align:right;">

4.3781270

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

2.2930632

</td>

<td style="text-align:right;">

2.2930632

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.6429468

</td>

<td style="text-align:right;">

1.6429468

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.6154409

</td>

<td style="text-align:right;">

1.6154409

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

2.7424461

</td>

<td style="text-align:right;">

2.7424461

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.9606986

</td>

<td style="text-align:right;">

1.9606986

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

1.5000000

</td>

<td style="text-align:right;">

1.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.7200000

</td>

<td style="text-align:right;">

2.7200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

1.8500000

</td>

<td style="text-align:right;">

1.8500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.0500000

</td>

<td style="text-align:right;">

2.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.3100000

</td>

<td style="text-align:right;">

2.3100000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0812000

</td>

<td style="text-align:right;">

0.0812000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.1088000

</td>

<td style="text-align:right;">

0.1088000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0740000

</td>

<td style="text-align:right;">

0.0740000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0820000

</td>

<td style="text-align:right;">

0.0820000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0924000

</td>

<td style="text-align:right;">

0.0924000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0500000

</td>

<td style="text-align:right;">

0.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0700000

</td>

<td style="text-align:right;">

0.0700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0900000

</td>

<td style="text-align:right;">

0.0900000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

1.7412500

</td>

<td style="text-align:right;">

1.7412500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.9325000

</td>

<td style="text-align:right;">

2.9325000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.9337500

</td>

<td style="text-align:right;">

2.9337500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.4775000

</td>

<td style="text-align:right;">

2.4775000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

1.9537500

</td>

<td style="text-align:right;">

1.9537500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.8725000

</td>

<td style="text-align:right;">

2.8725000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

1.7412500

</td>

<td style="text-align:right;">

1.7412500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.9325000

</td>

<td style="text-align:right;">

2.9325000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.9337500

</td>

<td style="text-align:right;">

2.9337500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.4775000

</td>

<td style="text-align:right;">

2.4775000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

1.9537500

</td>

<td style="text-align:right;">

1.9537500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.8725000

</td>

<td style="text-align:right;">

2.8725000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.0570651

</td>

<td style="text-align:right;">

2.0570651

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.1971553

</td>

<td style="text-align:right;">

3.1971553

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.1233717

</td>

<td style="text-align:right;">

3.1233717

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.6406412

</td>

<td style="text-align:right;">

2.6406412

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.1911408

</td>

<td style="text-align:right;">

2.1911408

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.1270821

</td>

<td style="text-align:right;">

3.1270821

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.0822826

</td>

<td style="text-align:right;">

0.0822826

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1278862

</td>

<td style="text-align:right;">

0.1278862

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1249349

</td>

<td style="text-align:right;">

0.1249349

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1056256

</td>

<td style="text-align:right;">

0.1056256

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.0876456

</td>

<td style="text-align:right;">

0.0876456

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1250833

</td>

<td style="text-align:right;">

0.1250833

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

15.3527035

</td>

<td style="text-align:right;">

15.3527035

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

8.2778360

</td>

<td style="text-align:right;">

8.2778360

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

6.0710577

</td>

<td style="text-align:right;">

6.0710577

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

6.1780905

</td>

<td style="text-align:right;">

6.1780905

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

10.8341191

</td>

<td style="text-align:right;">

10.8341191

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

8.1412032

</td>

<td style="text-align:right;">

8.1412032

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

76.7635175

</td>

<td style="text-align:right;">

76.7635175

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

25.8682374

</td>

<td style="text-align:right;">

25.8682374

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

18.9720552

</td>

<td style="text-align:right;">

18.9720552

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

22.0646091

</td>

<td style="text-align:right;">

22.0646091

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

45.1421631

</td>

<td style="text-align:right;">

45.1421631

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

22.6144534

</td>

<td style="text-align:right;">

22.6144534

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

12.1532371

</td>

<td style="text-align:right;">

12.1532371

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

7.8194513

</td>

<td style="text-align:right;">

7.8194513

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

8.0041706

</td>

<td style="text-align:right;">

8.0041706

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

9.4673975

</td>

<td style="text-align:right;">

9.4673975

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

11.4095817

</td>

<td style="text-align:right;">

11.4095817

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

7.9946733

</td>

<td style="text-align:right;">

7.9946733

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.0586347

</td>

<td style="text-align:right;">

2.0586347

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.1798068

</td>

<td style="text-align:right;">

3.1798068

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.0284958

</td>

<td style="text-align:right;">

3.0284958

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.5577884

</td>

<td style="text-align:right;">

2.5577884

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.1498803

</td>

<td style="text-align:right;">

2.1498803

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.0315925

</td>

<td style="text-align:right;">

3.0315925

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.0823454

</td>

<td style="text-align:right;">

0.0823454

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1271923

</td>

<td style="text-align:right;">

0.1271923

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1211398

</td>

<td style="text-align:right;">

0.1211398

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1023115

</td>

<td style="text-align:right;">

0.1023115

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.0859952

</td>

<td style="text-align:right;">

0.0859952

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1212637

</td>

<td style="text-align:right;">

0.1212637

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

15.4172443

</td>

<td style="text-align:right;">

15.4172443

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

7.7774164

</td>

<td style="text-align:right;">

7.7774164

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

3.1284787

</td>

<td style="text-align:right;">

3.1284787

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

3.1389787

</td>

<td style="text-align:right;">

3.1389787

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

9.1228460

</td>

<td style="text-align:right;">

9.1228460

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

5.2478198

</td>

<td style="text-align:right;">

5.2478198

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

76.7049878

</td>

<td style="text-align:right;">

76.7049878

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

26.0093699

</td>

<td style="text-align:right;">

26.0093699

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

19.5664063

</td>

<td style="text-align:right;">

19.5664063

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

22.7793336

</td>

<td style="text-align:right;">

22.7793336

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

46.0085322

</td>

<td style="text-align:right;">

46.0085322

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

23.3267671

</td>

<td style="text-align:right;">

23.3267671

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

12.1439707

</td>

<td style="text-align:right;">

12.1439707

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

7.8621128

</td>

<td style="text-align:right;">

7.8621128

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

8.2549230

</td>

<td style="text-align:right;">

8.2549230

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

9.7740688

</td>

<td style="text-align:right;">

9.7740687

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

11.6285546

</td>

<td style="text-align:right;">

11.6285546

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

8.2464909

</td>

<td style="text-align:right;">

8.2464909

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

3.2712500

</td>

<td style="text-align:right;">

3.2712500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

6.3987500

</td>

<td style="text-align:right;">

6.3987500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

5.0062500

</td>

<td style="text-align:right;">

5.0062500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

4.3818750

</td>

<td style="text-align:right;">

4.3818750

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

3.7075000

</td>

<td style="text-align:right;">

3.7075000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

5.5325000

</td>

<td style="text-align:right;">

5.5325000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

7.7925545

</td>

<td style="text-align:right;">

7.7925545

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

9.3915223

</td>

<td style="text-align:right;">

9.3915223

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

6.9726784

</td>

<td style="text-align:right;">

6.9726784

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

6.0672197

</td>

<td style="text-align:right;">

6.0672197

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

6.5458663

</td>

<td style="text-align:right;">

6.5458663

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

8.2892908

</td>

<td style="text-align:right;">

8.2892908

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

58.0208261

</td>

<td style="text-align:right;">

58.0208261

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

31.8667432

</td>

<td style="text-align:right;">

31.8667432

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

28.2019090

</td>

<td style="text-align:right;">

28.2019090

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

27.7778746

</td>

<td style="text-align:right;">

27.7778746

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

43.3612023

</td>

<td style="text-align:right;">

43.3612023

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

33.2572574

</td>

<td style="text-align:right;">

33.2572574

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

7.8150259

</td>

<td style="text-align:right;">

7.8150259

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

9.1953427

</td>

<td style="text-align:right;">

9.1953427

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

5.9887901

</td>

<td style="text-align:right;">

5.9887901

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

5.2113018

</td>

<td style="text-align:right;">

5.2113018

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

6.0525342

</td>

<td style="text-align:right;">

6.0525342

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

7.2552635

</td>

<td style="text-align:right;">

7.2552635

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

58.1415337

</td>

<td style="text-align:right;">

58.1415337

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

30.4131425

</td>

<td style="text-align:right;">

30.4131425

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

16.4063210

</td>

<td style="text-align:right;">

16.4063210

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

15.9159231

</td>

<td style="text-align:right;">

15.9159230

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

38.7446663

</td>

<td style="text-align:right;">

38.7446663

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

23.7450164

</td>

<td style="text-align:right;">

23.7450164

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.8786791

</td>

<td style="text-align:right;">

1.8786791

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

2.1820119

</td>

<td style="text-align:right;">

2.1820119

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.7064337

</td>

<td style="text-align:right;">

1.7064337

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.7686680

</td>

<td style="text-align:right;">

1.7686680

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.8976328

</td>

<td style="text-align:right;">

1.8976328

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.9260226

</td>

<td style="text-align:right;">

1.9260226

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

3.7881905

</td>

<td style="text-align:right;">

3.7881905

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.9374621

</td>

<td style="text-align:right;">

2.9374621

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.2324203

</td>

<td style="text-align:right;">

2.2324203

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.2976312

</td>

<td style="text-align:right;">

2.2976312

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.9874239

</td>

<td style="text-align:right;">

2.9874239

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.6508069

</td>

<td style="text-align:right;">

2.6508069

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

3.7962178

</td>

<td style="text-align:right;">

3.7962178

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.8917929

</td>

<td style="text-align:right;">

2.8917929

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

1.9774801

</td>

<td style="text-align:right;">

1.9774801

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.0374249

</td>

<td style="text-align:right;">

2.0374249

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.8152890

</td>

<td style="text-align:right;">

2.8152890

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.3932186

</td>

<td style="text-align:right;">

2.3932186

</td>

<td style="text-align:right;">

0

</td>

</tr>

</tbody>

</table>

## Test 8: Indometh (n=6), Log, Extravascular

``` r
table_wres_rres(Wres8, Rres8,
                Caption = 'Indometh (n=6), Log, Extravascular')
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Indometh (n=6), Log,
Extravascular

</caption>

<thead>

<tr>

<th style="border-bottom:hidden" colspan="1">

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Pharmacokinetic
Parameters

</div>

</th>

<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px;">

Values

</div>

</th>

<th style="border-bottom:hidden" colspan="1">

</th>

</tr>

<tr>

<th style="text-align:right;">

Subject

</th>

<th style="text-align:left;">

PPTESTCD

</th>

<th style="text-align:left;">

WNL

</th>

<th style="text-align:right;">

NonCompart

</th>

<th style="text-align:right;">

WinNonlin

</th>

<th style="text-align:right;">

Difference

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9970667

</td>

<td style="text-align:right;">

0.9970667

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9476691

</td>

<td style="text-align:right;">

0.9476691

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8758261

</td>

<td style="text-align:right;">

0.8758261

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8671179

</td>

<td style="text-align:right;">

0.8671179

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.8752442

</td>

<td style="text-align:right;">

0.8752442

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2

</td>

<td style="text-align:left;">

Rsq

</td>

<td style="text-align:right;">

0.9039538

</td>

<td style="text-align:right;">

0.9039538

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9941335

</td>

<td style="text-align:right;">

0.9941335

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.9401933

</td>

<td style="text-align:right;">

0.9401933

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8603043

</td>

<td style="text-align:right;">

0.8603043

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8505077

</td>

<td style="text-align:right;">

0.8505077

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8544516

</td>

<td style="text-align:right;">

0.8544516

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

R2ADJ

</td>

<td style="text-align:left;">

Rsq\_adjusted

</td>

<td style="text-align:right;">

0.8902329

</td>

<td style="text-align:right;">

0.8902329

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9985323

</td>

<td style="text-align:right;">

\-0.9985323

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9734830

</td>

<td style="text-align:right;">

\-0.9734830

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9358558

</td>

<td style="text-align:right;">

\-0.9358558

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9311917

</td>

<td style="text-align:right;">

\-0.9311917

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9355449

</td>

<td style="text-align:right;">

\-0.9355449

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CORRXY

</td>

<td style="text-align:left;">

Corr\_XY

</td>

<td style="text-align:right;">

\-0.9507649

</td>

<td style="text-align:right;">

\-0.9507649

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

3.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

10.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZNPT

</td>

<td style="text-align:left;">

No\_points\_lambda\_z

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

9.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.1583205

</td>

<td style="text-align:right;">

0.1583205

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.3022800

</td>

<td style="text-align:right;">

0.3022800

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.4218926

</td>

<td style="text-align:right;">

0.4218926

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.4290762

</td>

<td style="text-align:right;">

0.4290762

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.2527478

</td>

<td style="text-align:right;">

0.2527478

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZ

</td>

<td style="text-align:left;">

Lambda\_z

</td>

<td style="text-align:right;">

0.3535205

</td>

<td style="text-align:right;">

0.3535205

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

5.0000000

</td>

<td style="text-align:right;">

5.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

1.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZLL

</td>

<td style="text-align:left;">

Lambda\_z\_lower

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0.7500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZUL

</td>

<td style="text-align:left;">

Lambda\_z\_upper

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

4.3781270

</td>

<td style="text-align:right;">

4.3781270

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

2.2930632

</td>

<td style="text-align:right;">

2.2930632

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.6429468

</td>

<td style="text-align:right;">

1.6429468

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.6154409

</td>

<td style="text-align:right;">

1.6154409

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

2.7424461

</td>

<td style="text-align:right;">

2.7424461

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

LAMZHL

</td>

<td style="text-align:left;">

HL\_Lambda\_z

</td>

<td style="text-align:right;">

1.9606986

</td>

<td style="text-align:right;">

1.9606986

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TLAG

</td>

<td style="text-align:left;">

Tlag

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TMAX

</td>

<td style="text-align:left;">

Tmax

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0.2500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

1.5000000

</td>

<td style="text-align:right;">

1.5000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

2.0300000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.7200000

</td>

<td style="text-align:right;">

2.7200000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

1.8500000

</td>

<td style="text-align:right;">

1.8500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.0500000

</td>

<td style="text-align:right;">

2.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAX

</td>

<td style="text-align:left;">

Cmax

</td>

<td style="text-align:right;">

2.3100000

</td>

<td style="text-align:right;">

2.3100000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0812000

</td>

<td style="text-align:right;">

0.0812000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.1088000

</td>

<td style="text-align:right;">

0.1088000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0740000

</td>

<td style="text-align:right;">

0.0740000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0820000

</td>

<td style="text-align:right;">

0.0820000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CMAXD

</td>

<td style="text-align:left;">

Cmax\_D

</td>

<td style="text-align:right;">

0.0924000

</td>

<td style="text-align:right;">

0.0924000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

TLST

</td>

<td style="text-align:left;">

Tlast

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

8.0000000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0500000

</td>

<td style="text-align:right;">

0.0500000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0.0800000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0700000

</td>

<td style="text-align:right;">

0.0700000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0.0600000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLST

</td>

<td style="text-align:left;">

Clast

</td>

<td style="text-align:right;">

0.0900000

</td>

<td style="text-align:right;">

0.0900000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

1.7193653

</td>

<td style="text-align:right;">

1.7193653

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.8891436

</td>

<td style="text-align:right;">

2.8891436

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.8817113

</td>

<td style="text-align:right;">

2.8817113

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.4442459

</td>

<td style="text-align:right;">

2.4442459

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

1.9211984

</td>

<td style="text-align:right;">

1.9211984

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCLST

</td>

<td style="text-align:left;">

AUClast

</td>

<td style="text-align:right;">

2.8413138

</td>

<td style="text-align:right;">

2.8413138

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

1.7193653

</td>

<td style="text-align:right;">

1.7193653

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.8891436

</td>

<td style="text-align:right;">

2.8891436

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.8817113

</td>

<td style="text-align:right;">

2.8817113

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.4442459

</td>

<td style="text-align:right;">

2.4442459

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

1.9211984

</td>

<td style="text-align:right;">

1.9211984

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCALL

</td>

<td style="text-align:left;">

AUCall

</td>

<td style="text-align:right;">

2.8413138

</td>

<td style="text-align:right;">

2.8413138

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.0351804

</td>

<td style="text-align:right;">

2.0351804

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.1537989

</td>

<td style="text-align:right;">

3.1537989

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.0713330

</td>

<td style="text-align:right;">

3.0713330

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.6073871

</td>

<td style="text-align:right;">

2.6073871

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

2.1585892

</td>

<td style="text-align:right;">

2.1585892

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFO

</td>

<td style="text-align:left;">

AUCINF\_obs

</td>

<td style="text-align:right;">

3.0958959

</td>

<td style="text-align:right;">

3.0958959

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.0814072

</td>

<td style="text-align:right;">

0.0814072

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1261520

</td>

<td style="text-align:right;">

0.1261520

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1228533

</td>

<td style="text-align:right;">

0.1228533

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1042955

</td>

<td style="text-align:right;">

0.1042955

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.0863436

</td>

<td style="text-align:right;">

0.0863436

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFOD

</td>

<td style="text-align:left;">

AUCINF\_D\_obs

</td>

<td style="text-align:right;">

0.1238358

</td>

<td style="text-align:right;">

0.1238358

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

15.5177942

</td>

<td style="text-align:right;">

15.5177942

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

8.3916343

</td>

<td style="text-align:right;">

8.3916343

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

6.1739217

</td>

<td style="text-align:right;">

6.1739217

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

6.2568848

</td>

<td style="text-align:right;">

6.2568848

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

10.9974979

</td>

<td style="text-align:right;">

10.9974979

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEO

</td>

<td style="text-align:left;">

AUC\_.Extrap\_obs

</td>

<td style="text-align:right;">

8.2232127

</td>

<td style="text-align:right;">

8.2232127

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

77.5889712

</td>

<td style="text-align:right;">

77.5889712

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

26.2238573

</td>

<td style="text-align:right;">

26.2238573

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

19.2935053

</td>

<td style="text-align:right;">

19.2935053

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

22.3460171

</td>

<td style="text-align:right;">

22.3460171

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

45.8229078

</td>

<td style="text-align:right;">

45.8229078

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZFO

</td>

<td style="text-align:left;">

Vz\_F\_obs

</td>

<td style="text-align:right;">

22.8422576

</td>

<td style="text-align:right;">

22.8422576

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

12.2839234

</td>

<td style="text-align:right;">

12.2839234

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

7.9269481

</td>

<td style="text-align:right;">

7.9269481

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

8.1397881

</td>

<td style="text-align:right;">

8.1397881

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

9.5881430

</td>

<td style="text-align:right;">

9.5881430

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

11.5816384

</td>

<td style="text-align:right;">

11.5816384

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLFO

</td>

<td style="text-align:left;">

Cl\_F\_obs

</td>

<td style="text-align:right;">

8.0752068

</td>

<td style="text-align:right;">

8.0752068

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.0367500

</td>

<td style="text-align:right;">

2.0367500

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.1364504

</td>

<td style="text-align:right;">

3.1364504

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.9764572

</td>

<td style="text-align:right;">

2.9764572

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.5245343

</td>

<td style="text-align:right;">

2.5245343

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

2.1173287

</td>

<td style="text-align:right;">

2.1173287

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFP

</td>

<td style="text-align:left;">

AUCINF\_pred

</td>

<td style="text-align:right;">

3.0004063

</td>

<td style="text-align:right;">

3.0004063

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.0814700

</td>

<td style="text-align:right;">

0.0814700

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1254580

</td>

<td style="text-align:right;">

0.1254580

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1190583

</td>

<td style="text-align:right;">

0.1190583

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1009814

</td>

<td style="text-align:right;">

0.1009814

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.0846931

</td>

<td style="text-align:right;">

0.0846931

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCIFPD

</td>

<td style="text-align:left;">

AUCINF\_D\_pred

</td>

<td style="text-align:right;">

0.1200163

</td>

<td style="text-align:right;">

0.1200163

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

15.5829013

</td>

<td style="text-align:right;">

15.5829013

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

7.8849267

</td>

<td style="text-align:right;">

7.8849267

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

3.1831752

</td>

<td style="text-align:right;">

3.1831752

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

3.1803265

</td>

<td style="text-align:right;">

3.1803265

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

9.2630996

</td>

<td style="text-align:right;">

9.2630996

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUCPEP

</td>

<td style="text-align:left;">

AUC\_.Extrap\_pred

</td>

<td style="text-align:right;">

5.3023656

</td>

<td style="text-align:right;">

5.3023656

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

77.5291765

</td>

<td style="text-align:right;">

77.5291765

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

26.3689077

</td>

<td style="text-align:right;">

26.3689077

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

19.9084941

</td>

<td style="text-align:right;">

19.9084941

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

23.0793917

</td>

<td style="text-align:right;">

23.0793917

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

46.7158621

</td>

<td style="text-align:right;">

46.7158621

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

VZFP

</td>

<td style="text-align:left;">

Vz\_F\_pred

</td>

<td style="text-align:right;">

23.5692251

</td>

<td style="text-align:right;">

23.5692252

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

12.2744566

</td>

<td style="text-align:right;">

12.2744566

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

7.9707940

</td>

<td style="text-align:right;">

7.9707940

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

8.3992473

</td>

<td style="text-align:right;">

8.3992473

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

9.9028165

</td>

<td style="text-align:right;">

9.9028165

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

11.8073306

</td>

<td style="text-align:right;">

11.8073306

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

CLFP

</td>

<td style="text-align:left;">

Cl\_F\_pred

</td>

<td style="text-align:right;">

8.3322048

</td>

<td style="text-align:right;">

8.3322048

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

3.2965543

</td>

<td style="text-align:right;">

3.2965543

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

6.4082620

</td>

<td style="text-align:right;">

6.4082620

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

5.0353383

</td>

<td style="text-align:right;">

5.0353382

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

4.3990453

</td>

<td style="text-align:right;">

4.3990453

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

3.7299741

</td>

<td style="text-align:right;">

3.7299741

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCLST

</td>

<td style="text-align:left;">

AUMClast

</td>

<td style="text-align:right;">

5.5775672

</td>

<td style="text-align:right;">

5.5775672

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

7.8178588

</td>

<td style="text-align:right;">

7.8178588

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

9.4010343

</td>

<td style="text-align:right;">

9.4010343

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

7.0017667

</td>

<td style="text-align:right;">

7.0017667

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

6.0843899

</td>

<td style="text-align:right;">

6.0843899

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

6.5683405

</td>

<td style="text-align:right;">

6.5683405

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFO

</td>

<td style="text-align:left;">

AUMCINF\_obs

</td>

<td style="text-align:right;">

8.3343579

</td>

<td style="text-align:right;">

8.3343579

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

57.8330281

</td>

<td style="text-align:right;">

57.8330281

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

31.8345005

</td>

<td style="text-align:right;">

31.8345005

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

28.0847466

</td>

<td style="text-align:right;">

28.0847466

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

27.6994849

</td>

<td style="text-align:right;">

27.6994849

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

43.2128382

</td>

<td style="text-align:right;">

43.2128382

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEO

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_obs

</td>

<td style="text-align:right;">

33.0774222

</td>

<td style="text-align:right;">

33.0774222

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

7.8403303

</td>

<td style="text-align:right;">

7.8403303

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

9.2048546

</td>

<td style="text-align:right;">

9.2048546

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

6.0178784

</td>

<td style="text-align:right;">

6.0178784

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

5.2284721

</td>

<td style="text-align:right;">

5.2284721

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

6.0750083

</td>

<td style="text-align:right;">

6.0750083

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCIFP

</td>

<td style="text-align:left;">

AUMCINF\_pred

</td>

<td style="text-align:right;">

7.3003307

</td>

<td style="text-align:right;">

7.3003307

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

57.9538845

</td>

<td style="text-align:right;">

57.9538845

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

30.3817147

</td>

<td style="text-align:right;">

30.3817147

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

16.3270188

</td>

<td style="text-align:right;">

16.3270188

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

15.8636553

</td>

<td style="text-align:right;">

15.8636553

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

38.6013327

</td>

<td style="text-align:right;">

38.6013327

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

AUMCPEP

</td>

<td style="text-align:left;">

AUMC\_.Extrap\_pred

</td>

<td style="text-align:right;">

23.5984312

</td>

<td style="text-align:right;">

23.5984312

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.9173089

</td>

<td style="text-align:right;">

1.9173089

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

2.2180490

</td>

<td style="text-align:right;">

2.2180490

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.7473430

</td>

<td style="text-align:right;">

1.7473430

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.7997556

</td>

<td style="text-align:right;">

1.7997556

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.9414830

</td>

<td style="text-align:right;">

1.9414830

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTEVLST

</td>

<td style="text-align:left;">

MRTlast

</td>

<td style="text-align:right;">

1.9630240

</td>

<td style="text-align:right;">

1.9630240

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

3.8413591

</td>

<td style="text-align:right;">

3.8413591

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.9808604

</td>

<td style="text-align:right;">

2.9808604

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.2797159

</td>

<td style="text-align:right;">

2.2797159

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.3335200

</td>

<td style="text-align:right;">

2.3335200

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

3.0428858

</td>

<td style="text-align:right;">

3.0428858

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTEVIFO

</td>

<td style="text-align:left;">

MRTINF\_obs

</td>

<td style="text-align:right;">

2.6920666

</td>

<td style="text-align:right;">

2.6920666

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

3.8494318

</td>

<td style="text-align:right;">

3.8494318

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.9348000

</td>

<td style="text-align:right;">

2.9348000

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.0218260

</td>

<td style="text-align:right;">

2.0218260

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.0710640

</td>

<td style="text-align:right;">

2.0710640

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.8691853

</td>

<td style="text-align:right;">

2.8691853

</td>

<td style="text-align:right;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:left;">

MRTEVIFP

</td>

<td style="text-align:left;">

MRTINF\_pred

</td>

<td style="text-align:right;">

2.4331140

</td>

<td style="text-align:right;">

2.4331140

</td>

<td style="text-align:right;">

0

</td>

</tr>

</tbody>

</table>

# Session Information

``` r
devtools::session_info()
```

    ##  setting  value                       
    ##  version  R version 3.5.1 (2018-07-02)
    ##  system   x86_64, mingw32             
    ##  ui       RTerm                       
    ##  language (EN)                        
    ##  collate  Korean_Korea.949            
    ##  tz       Asia/Seoul                  
    ##  date     2018-08-07                  
    ## 
    ##  package     * version date       source                             
    ##  assertthat    0.2.0   2017-04-11 CRAN (R 3.5.0)                     
    ##  backports     1.1.2   2017-12-13 CRAN (R 3.5.0)                     
    ##  base        * 3.5.1   2018-07-02 local                              
    ##  bindr         0.1.1   2018-03-13 CRAN (R 3.5.0)                     
    ##  bindrcpp    * 0.2.2   2018-03-29 CRAN (R 3.5.0)                     
    ##  colorspace    1.3-2   2016-12-14 CRAN (R 3.5.0)                     
    ##  compiler      3.5.1   2018-07-02 local                              
    ##  crayon        1.3.4   2018-06-08 Github (gaborcsardi/crayon@3e751fb)
    ##  datasets    * 3.5.1   2018-07-02 local                              
    ##  devtools      1.13.6  2018-06-27 CRAN (R 3.5.0)                     
    ##  digest        0.6.15  2018-01-28 CRAN (R 3.5.0)                     
    ##  dplyr       * 0.7.6   2018-06-29 CRAN (R 3.5.0)                     
    ##  evaluate      0.11    2018-07-17 CRAN (R 3.5.1)                     
    ##  glue          1.3.0   2018-07-17 CRAN (R 3.5.1)                     
    ##  graphics    * 3.5.1   2018-07-02 local                              
    ##  grDevices   * 3.5.1   2018-07-02 local                              
    ##  highr         0.7     2018-06-09 CRAN (R 3.5.0)                     
    ##  hms           0.4.2   2018-03-10 CRAN (R 3.5.0)                     
    ##  htmltools     0.3.6   2017-04-28 CRAN (R 3.5.0)                     
    ##  httr          1.3.1   2017-08-20 CRAN (R 3.5.0)                     
    ##  kableExtra  * 0.9.0   2018-05-21 CRAN (R 3.5.0)                     
    ##  knitr       * 1.20    2018-02-20 CRAN (R 3.5.0)                     
    ##  magrittr      1.5     2014-11-22 CRAN (R 3.5.0)                     
    ##  memoise       1.1.0   2017-04-21 CRAN (R 3.5.0)                     
    ##  methods     * 3.5.1   2018-07-02 local                              
    ##  munsell       0.5.0   2018-06-12 CRAN (R 3.5.0)                     
    ##  NonCompart  * 0.4.4   2018-08-06 CRAN (R 3.5.1)                     
    ##  pillar        1.3.0   2018-07-14 CRAN (R 3.5.1)                     
    ##  pkgconfig     2.0.1   2017-03-21 CRAN (R 3.5.0)                     
    ##  plyr          1.8.4   2016-06-08 CRAN (R 3.5.0)                     
    ##  purrr         0.2.5   2018-05-29 CRAN (R 3.5.0)                     
    ##  R6            2.2.2   2017-06-17 CRAN (R 3.5.0)                     
    ##  Rcpp          0.12.18 2018-07-23 CRAN (R 3.5.1)                     
    ##  readr         1.1.1   2017-05-16 CRAN (R 3.5.0)                     
    ##  rlang         0.2.1   2018-05-30 CRAN (R 3.5.0)                     
    ##  rmarkdown     1.10    2018-06-11 CRAN (R 3.5.0)                     
    ##  rprojroot     1.3-2   2018-01-03 CRAN (R 3.5.0)                     
    ##  rstudioapi    0.7     2017-09-07 CRAN (R 3.5.0)                     
    ##  rvest         0.3.2   2016-06-17 CRAN (R 3.5.0)                     
    ##  scales        0.5.0   2017-08-24 CRAN (R 3.5.0)                     
    ##  stats       * 3.5.1   2018-07-02 local                              
    ##  stringi       1.2.4   2018-07-20 CRAN (R 3.5.1)                     
    ##  stringr       1.3.1   2018-05-10 CRAN (R 3.5.0)                     
    ##  tibble        1.4.2   2018-01-22 CRAN (R 3.5.0)                     
    ##  tidyr       * 0.8.1   2018-05-18 CRAN (R 3.5.0)                     
    ##  tidyselect    0.2.4   2018-02-26 CRAN (R 3.5.0)                     
    ##  tools         3.5.1   2018-07-02 local                              
    ##  utils       * 3.5.1   2018-07-02 local                              
    ##  viridisLite   0.3.0   2018-02-01 CRAN (R 3.5.0)                     
    ##  withr         2.1.2   2018-03-15 CRAN (R 3.5.0)                     
    ##  xml2          1.2.0   2018-01-24 CRAN (R 3.5.0)                     
    ##  yaml          2.2.0   2018-07-25 CRAN (R 3.5.1)

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
