# SPRE.R

SPRE.R is a Semiparanetric Data Estimation model to predict outcomes for single-subject and small group designs for use in 'R'.   

SPRE predicts residuals, normal Q-Q, predictions from a change point, effective clinical stability, SPRE Weibull probability, computations for the highest F statistic, all predictions, and the probability curve for effective clinical stability.

## SPRE.R running in R

The data must be named: Dataset <- read.csv("~/Documents/NameOfYourFile.csv"). There are 2 Excel data columns: FData  Session. This is a quick and simpler way to run SPRE.R.

```
Weissman-Miller, D., Shotwell, M.P. & Miller, R.J. (2012).  New single-subject and small-n design in occupational therapy:  Application to weight loss in obesity.  American Journal of Occupational Therapy, 66, 455-462. DOI: http://dx.doi.org/10.5014/ajot.2012.004788

Weissman-Miller, D. (2013). Novel point estimation from a Semiparametric Ratio Estimator (SPRE): Long-term health outcomes from short-term linear data, with application to weight loss in obesity. International Journal of Biostatistics, 9(2): 175-184 DOI: http://dx.doi.org/10.1515/ijb-2012-0049

```

## Getting Started

These instructions will get you a copy of the 'R' project up and running on your local machine for use of SPRE in Occupational Therapy (OT) and Personalized Medicine research and therapy.

### Prerequisites

* Install the R computer language from www.r-project.org on your machine.

### Installing

A step by step series of examples that tell you how to get a development env running

* GOTO www.github.com/Debbiewm5/SPRE.R and download the code files onto your computer.
* Open the download SPRE.R-master to view files.
* Open the \data folder, then copy and paste the 2 sample files into your Documents folder. On linux in a Chromebook, paste into Downloads folder. Note that FData, the y axis, is given as an order statistic to analyse data in SPRE. Some of the data referenced below is naturally ordered.  If not, the FData column is ordered together with the asocciated session number in Excel and saved as '.csv'.  Read about this in the wiki.

* In Windows:
	- Open the 'R' program
        - Set the R session to the Documents folder or wherever your data is located:
	- Windows: setwd("~/")  
	- Mac: setwd("~/Documents")  
	- Linux: setwd("~/Downloads")

		- To Enter Code in the Windows console type:				
		- setwd("~/") and Enter  
		- Dataset <- read.csv("~/NameOfYourFile.csv") and Enter   
		- source("test191.R") and Enter  
	
### Run the EASY WAY in Windows:

```
		- File -> Change dir... -> Documents
		- Dataset <- read.csv("~/NameOfYourFile.csv") and Enter
		- File -> Source R Code -> test192

	- Use the Dataset given in 'Dataset' for practice, substitute in NameOfYourFile.csv
	- - Use the Dataset given in the file downloaded in 'data' for practice, substitute in NameOfYourFile.csv as... 	
		- Dataset <- read.csv("~/HEAT_LILLY_Heat.Stat.14_2016.csv")
	- The results from the dataset appear in the R console 
	- The resulting 'change point', health predictions from the change point, the effective clinical stability and the comparative values for the 	  	  	  Weibull predictive curve are printed on the console.
	
	- There are 2 pages of plots: 1a) The SPRE residuals from the initial OLS linear regression from the client or patients' data; 1b)The SPRE 	  	  	  Q-Q plot often useful to compare the pattern of these responses to another set of responses in SPRE; 1c)Finally the predictive curve for 	  	          that participant.
	- The 2nd. plot sheet indicates the session number beyond the change point where the participant becomes effectively clinically stable for the 	  	  	  treatment.
	- The probability for the Weibull CDF used for prediction is given on the console.

```

* In Mac:
	- To Enter Code in the console type:
		- setwd("~/Documents") 
		- Dataset <- read.csv("~/Documents/NameOfYourFile.csv") and Enter
		- source("test191.R")

### Run the EASY WAY in Mac:

```
		- Misc -> Change Working Directory... -> Documents
		- Dataset <- read.csv("~/NameOfYourFile.csv") and Enter
		- File -> Source File... -> test192.R  (or Documents - iCloud)

	- The rest of the instructions above are the same.

```

## Downloaded from 'http://www.github.com/Debbiewm5/SPRE.R'

```
Examples: In Windows downloaded from GitHub

> In the console type: Dataset <- read.csv("~/HEAT_LILLY_Heat.Stat.14_2016.csv".csv") and Enter
    FData Session
1  4.8462       4
2  4.6429       3
3  4.6250      14
4  4.4375      11
5  4.3792      13
6  4.3750       6
7  4.1765      10
8  4.1429       8
9  4.1333      12
10 4.0667       2
11 4.0000       1
12 4.0000       5
13 3.9375       7 
14 3.8750       9


```

### Results to compute the 'change point'where the participant adapts to the treatment.

```
Backwards stepwise OLS linear regression 
 
        Rsquare Fstatistic FstatNumer. FstatDenom.   Pvalue Session
 [1,] 0.0000000         NA           4           0       NA       4
 [2,] 1.0000000        NaN           1           0      NaN       3
 [3,] 0.2419160  0.3191150           1           1 0.672642      14
 [4,] 0.3279450  0.9759470           1           2 0.427335      11
 [5,] 0.4544550  2.4990800           1           3 0.212049      13
 [6,] 0.2109580  1.0694400           1           4 0.359494       6
 [7,] 0.1874930  1.1537900           1           5 0.331837      10
 [8,] 0.1131720  0.7656880           1           6 0.415204       8
 [9,] 0.1705230  1.4390500           1           7 0.269317      12
[10,] 0.0107366  0.0868249           1           8 0.775755       2
[11,] 0.0136682  0.1247190           1           9 0.732103       1
[12,] 0.0283044  0.2912880           1          10 0.601207       5
[13,] 0.0276055  0.3122820           1          11 0.587479       7
[14,] 0.0131730  0.1601870           1          12 0.696018       9
 
  
 slope parameter: betaHat1 = [1] 0.008139
 
 shape parameter: k_value -Session 
 0.1118 
 
 scale parameter is the value (tau) at the change point - given as the highest F statistic

Session number for - Session 
     13 
 
 Session value for - Pvalue 
 0.212 
 
 If R^2 is less than 0.45 at the change point, do not predict for inference 

```

### An example if getting some predictive data, from the 'change point', out of the system.
The first fifteen rows of data from the 'change point' for this participant.

```
     Session    Ratio tau k.factor Predictions
 [1,]      NA       NA  NA       NA          NA
 [2,]      15 0.995800  13 0.111829     4.36081
 [3,]      16 0.996054  13 0.111829     4.34360
 [4,]      17 0.996280  13 0.111829     4.32744
 [5,]      18 0.996481  13 0.111829     4.31222
 [6,]      19 0.996662  13 0.111829     4.29782
 [7,]      20 0.996825  13 0.111829     4.28417
 [8,]      21 0.996972  13 0.111829     4.27120
 [9,]      22 0.997107  13 0.111829     4.25885
[10,]      23 0.997230  13 0.111829     4.24705
[11,]      24 0.997343  13 0.111829     4.23577
[12,]      25 0.997447  13 0.111829     4.22495
[13,]      26 0.997544  13 0.111829     4.21458
[14,]      27 0.997633  13 0.111829     4.20460
[15,]      28 0.997716  13 0.111829     4.19500

```

## Add final data plots in R for Windows and Mac 
### Treatment prognosis and evaluation can be determined from these plots. 
For the first run of the data on your computer, click 'Windows' and last plot number, usually 3. Then under the 'Windows' tab, click the next R Graphics device 2 (Active). This enables you to see the previous plot. Plots can be switched in Windows. 
Start with plots 1a), 1b), and 1c), then view plot 2 to analyse the data.

## Saving and Comparing Participant Plot Data
### Saving R Plots
Display the plot and then goto File. Then click 'save as' for options. 
### Comparing Participant Plot Data
	* For only one participant, save the 2 plots and then exit R. This will clear plot data.
	* To compare plots, do not exit R but start a new analysis in the console. Then, for example, 4 plots can be directly compared.
		- In the Windows console type:
		- setwd("~/") and Enter
		- Dataset <- read.csv("~/NameOfFirstFile.csv") and Enter
		- source("test191.R") and Enter

		- In the Windows console type:
		- setwd("~/") and Enter
		- Dataset <- read.csv("~/NameOfAnotherFile.csv") and Enter
		- source("test191.R") and Enter
	* All 4 plots will be available for comparison in 'Windows'.

## Deployment and transforming Excel data
Additional notes on transforming the data in Excel to analyse in SPRE are given in a Github Wiki.

## Built With
* R programming language

## Versioning

We use [Git] for versioning. 

## Authors

* **Deborah Weissman-Miller** - *Initial work* - (https://github.com/Debbiewm5/SPRE.R)

## License

This project is licensed under GPL-2

## Acknowledgments

* Inspiration - the HEAT measure in OT by Dr. Mary Shotwell, Brenau University, Gainesville, GA

## Further References

```
Weissman-Miller, D., Miller, R., & Shotwell, M. (2017). Translational Research for Occupational Therapy: Using SPRE in Hippotherapy for Children with Developmental Disabilities. Occupational Therapy International.  https://www.hindawi.com/journals/oti/2017/2305402/

Weissman-Miller, D. (2016). On Predicting Survival in Prostate Cancer: Using an extended maximum spacing method at the change point of the semiparametric ratio estimator (SPRE).  International Journal of Statistics and Probability, 5, No. 2, 2016. http://www.ccsenet.org/journal/index.php/ijsp/article/view/57236)
DOI: http://dx.doi.org/10.5539/ijsp.v5n2p19

Weissman-Miller, D. Holmes, W. (2015). Novel Low-Error Interpolation Method for a Fall Prevention Program Using the Single- Subject Design Statistical Model SPRE. International Journal of Current Research in Life Sciences. Vol. 4, No. 10, pp.429-435 October 21015. http://www.ijcrls.com

Weissman-Miller, D. and Graham, K. 2015. Novel scale development for fear of falling and falls: Analysed using a Semiparametric Ratio Estimator (SPRE). International Journal of Statistics and Probability, 4, No. 3, 2015. DOI: http://dx.doi.org/10.5539/ijsp.v4n3p161

Weissman-Miller, D., Graham, K. C. (2014).  Novel Scale Development for Fear of Falling and Falls: Analysed using a Semiparametric Ratio Estimator (SPRE). In the Joint Statistical Meeting (JSM) Proceedings, Non-Parametric Section.  Boston, MA: American Statistical Association (2014). 

```