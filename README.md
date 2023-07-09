# Interval Censored Survival Models in Stan

This repository contains the Stan codes that I wrote for various Interval Censored Survival Models.  

## Motivation

Interval censored models are not well-covered in the Stan manuel, despite its popularity in real-life applications.

## Prerequisite 

I assume users are familiar enough with [Stan](https://mc-stan.org/), a Bayesian Probabilistic programming language.  

## Data
The data are taked from a 1980 paper in the Organizational Behavior and Human Performance journal, written by Morrison and Schmittlein.  A copy is provided [here](https://github.com/kaloklee/Interval_surv/blob/main/morrison_schmittlein_obhp_80.pdf).

## Where to begin

`Weibull_interval.stan` and `WeibullGamma_interval.stan` are the two main Stan files that I wrote for the Weibull and the Weibull-Gamma distributions.  

The user may run `Interval_censored setup.R` to follow what I am doing.  I am sure it is straight-forward to follow.

## Changes

Initially, two Weibull models are considered for this exercise . A standard Weibull model can leverage the built-in Stan functions. A Weibull-Gamma model requires customized functions.  

In a follow-up exercise I consider finite mixture models.  For the Morrison and Schmittlein data at hand, these models are not necessary.  I nevertheless share what I did as a way to summarize the finite mixture coding, which may be useful for future use.

A noteworth point is that in the Stan code's generating quantities block, I consider both parameter uncertainty and process uncertainy.  

I use cmdstanR as interface with Stan. Users are welcome to use the codes as they see fit.  

Please cite this repository if you do use the codes.  Thank you.

BONUS: I coded a mixture model of Weibull and Exponential in MLE.  Putting them in this repository for future use.
