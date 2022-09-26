# Interval Censored Survival Models in Stan

This repository contains the Stan codes that I wrote for various Interval Censored Survival Models.

The data used are from a 1980 paper in the Organizational Behavior and Human Performance journal, written by Morrison and Schmittlein.  A copy is attached.

Interval censored models are not well-covered in the Stan manuel, despite its popularity in real-life applications.

Initially, two Weibull models are considered for this exercise . A standard Weibull model can leverage the built-in Stan functions. A Weibull-Gamma model requires customized functions.  

In a follow-up exercise I consider finite mixture models.  For the Morrison and Schmittlein data at hand, these models are not necessary.  I nevertheless share what I did as a way to summarize the finite mixture coding for future use.

A noteworth point is that in the Stan code's generating quantities block, I consider both parameter uncertainty and process uncertainy.  

I use cmdstanR as interface with Stan. Users are welcome to use the codes as they see fit.  

Please cite this repository if you do use the codes.  Thank you.
