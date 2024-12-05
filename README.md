The R file that is attached to this repository uses Nimble to build an MCMC using a data that details Duchene Muscular Dystrophy, a sex-linked genetic disease mostly affecting boys.
Girls can be carriers without experiencing symptoms and thus it aims to develop a test to detect female carriers.
The study involved an independent sample of 38 known carriers and 82 non-carriers. The results showed that 20 carriers and 1 non carrier tested positive. 


## Question 1

Let Se denote the sensitivity of the test, i.e., the probability of a woman testing positive
given she is a carrier, and Sp the specificity of the test, i.e., the probability of a woman
testing negative given that she is not a carrier. Assuming Se ∼ Beta(a_Se; b_Se) and
Sp ∼ Beta(a_Sp; b_Sp), derive the posterior distributions for Se and Sp.

## Question 2

Using hyperparameter values for aSe, bSe, aSp, bSp that reflect constant prior probability
on the parameter space of Se and Sp, report suitable summary statistics for the posterior
distributions of Se and Sp. Briefly comment on the test performance after observing the
data of this experiment.

## Question 3

The probability that a woman is a carrier given that she tests positive is called the
Positive Predictive Value, PPV. Use Bayes theorem to write down the formula for the
PPV in terms of Se, Sp and π where π denotes the unconditional probability of being a
carrier. Note: the unconditional probability of being a carrier is called the prevalence of
the disease.

## Question 4

Assume that the prior distribution for the prevalence of the disease, expressed on what
is called a logit scale, log(π/(1 − π)), is N(0; 102). Use your results from question 2 on
the posterior distributions for Se and Sp and simple Monte Carlo methods to generate
a sample from the posterior distribution of PPV. Summarise the results, examining the
mean value and a histogram of the PPV values. Hint: note that since this is a case
control study (i.e., the number of women that are carriers and not carriers was fixed
beforehand), the data provide no information on π, and so its posterior distribution is
identical to the prior distribution.

## Question 5

Now assume that we have access to data on the number of carriers from 15 hospitals
where genetic studies on the prevalence of the DMD-causing allele were conducted. In each study, genetic tests were administered to young female patients enrolled in an HPV screening (ni) and the resulting
number of carriers recorded (ncarriers). In addition, the following information is also
available from each hospital:

• X hp: proportion of females enrolled in the HPV screening with diagnosed heart
disorders.
• X hi: proportion of Hispanic females enrolled in the HPV screening.
Model the number of carriers from each hospital i using a Binomial distribution with
parameter πi whose prior depends on the following relation:
logit(πi) = β0 + βW hZhpi + βHiZhii,

where Z are the normalised covariates. All beta parameters are given Normal priors with
mean 0 and standard deviation 5.


## Question 6

Do the inferences from the model fitted in Question 5 (with both
covariates) make sense? If the model fits well, the observed data should look plausible
under the posterior predictive distribution

