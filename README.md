# Interval Censored Survival Models in Stan

This repository contains the Stan codes that I wrote for various Interval Censored Survival Models.

The data used are from a 1980 paper in the Organizational Behavior and Human Performance journal, written by Morrison and Schmittlein.  A copy is attached.

Interval censored models are not well-covered in the Stan manuel, despite its popularity in real-life applications.

Two Weibull models are considered for this exercise as of September 25th, 2022. A standard Weibull model can leverage the built-in Stan functions. A Weibull-Gamma model requires customized functions. 
I plan to add more models in the future.

A noteworth point is that in the Stan code's generating quantities block, I consider both parameter uncertainty and process uncertainy.  

I use cmdstanR as interface with Stan. Users are welcome to use the codes as they see fit.  

Please cite this repository if you do use the codes.  Thank you.
