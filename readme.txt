DATE: May 13, 2020
TITLE: Code and sample data to accompany "The Shape of Educational Inequality" by Quarles, Budak & Resnick, published in Science Advances
AUTHOR: Christopher L. Quarles, chrisquarles@gmail.com

This repository contains four files:
-- readme.txt: The file you're reading now.
-- mlecens.R: This file contains the R code for estimating student capital in a population of students. It contains one function, mlecens, which performs right-censored maximum likelihood estimation to fit a distribution to a data set. 
-- code from QBR paper.R: This contains the R code used to make (most of) the images and tables in the paper. Because our data is unavailable, all of the code will run on the sample_data.csv. If you want to make an image from the paper with your own data, you can just format your data like in sample_data.csv and then run the code in this file.
-- sample_data.csv: For privacy reasons, the data used in the paper is not available to the public. This dataset mimics the type of data used for the analysis. The dataset has 4 variables:
      - credits_earned = # of credits earned by a given student, rounded to the nearest positive integer
      - droppedout = FALSE if the student graduated or transferred, TRUE otherwise 
      - transferred = TRUE iff the student transferred to a 4-year college 
      - transnograd = TRUE iff the student transferred but didn't graduate


IF YOU JUST WANT TO CALCULATE THE STUDENT CAPITAL OF A GROUP OF STUDENTS:

 - Make sure that your cohort is large enough. In simulations based on real data, the standard error of the estimated average student capital was roughly: SE = 150/sqrt(sample size).
 - Also, make sure that a middling number of your students dropped out. Otherwise, you won't observe enough students' capital to make an accurate inference. I don't have a good rule of thumb here, but 20% or fewer dropouts probably won't work. Nor will >90% dropouts.
 - Save your data in the same format as sample_data.csv, or you can copy and paste over the sample data. You only need two variables: credits_earned and droppedout. droppedout can be either TRUE/FALSE or 1/0. 
 - Make sure all the files are in the same directory. 
 - Make sure that you have the VGAM package installed. You can run install.packages("VGAM") or install it through Tools menu in RStudio.
 - Run the following lines in R. You'll have to change the file name to match your file. (The sample file should give q=.9917 and mu_s=120.2.) 

source("mlecens.R")
coldat <- read.csv("sample_data.csv")
q = mlecens(x=coldat$credits_earned, yc=coldat$droppedout)  # This returns the "per-credit retention rate"
mu_s = 1/(1-q)  # This returns the "average student capital", measured in credits.
