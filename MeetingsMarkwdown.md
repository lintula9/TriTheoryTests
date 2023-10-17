---
title: "Kuuma Peruna (KP) simulation results"
author: "Sakari"
date: "2023-10-17"
output: ioslides_presentation
fig_caption: no
---



## General idea

Three models of how psychopathology is conceputalized to explain positive manifold.  

1. Latent variable -model (LV)  

2. Network model  

3. Secondary Impairments (SI)  


Latent process is a 1-factor causal model, with measurement errors.  

Network model is a item-level causal model.  

SI is a additive confounding causal model that explains the positive manifold by stating that psychopathology is highly correlated because of the non-specific impairments that stress individuals.  
Widiger, Oltmans and others suggested this [I believe at least here](https://doi.org/10.1521/pedi_2021_35_530) as a response to p factor theory.  

## General idea
![A: simple (causal) network, B: Simple LV model, C: SI model represented with bonds that are synonymous to impairments.](figure/applicationFigure.png)

## General idea
They often make the same predictions at least in cross-sectional data:  

- [Bartholomew, Deary & Law (2009)](doi.org/10.1037/a0016262) explain how a factor model can be represented as caused by a random sample of, essentially confounding, variables.  

- [Golino & Epskamp (2017)](doi.org/10.1371/journal.pone.0174035) estimated latent variables using networks.  

- [Marsman et al., 2015](https://www.nature.com/articles/srep09050) estimated networks using latent variables.  

- [Christensen & Golino (2021)](https://doi.org/10.3758/s13428-020-01500-6) showed near singularity between factor and 'network loadings'.  

- ['Swiss cheese' Fried also commented on this](https://doi.org/10.1080/1047840X.2020.1853461).  

-And many others, for example [Molenaar](https://doi.org/10.1017/S0140525X10000798) was early on this.  


## Scenarios where they might differ are unclear.
- Longitudinal repeated cross-sectional (i.e., measurement invariance).  

- Intervention data distinguishes in some scenarios (but 'fat hand' etc problems).  

- Temporal analysis (e.g., CT-VAR models, lv-VAR models).  

- Independent vs. common pathways (e.g., genetic studies).  

- Within person data?  


## Longitudinal data, temporal analysis.

For example, sampling from VAR(1) process with symmetrical drift matrix does not create measurement invariance for 1 latent factor.  

Vars: 3.  

Drift matrix: symmetric (.3 in diagonal and off-diagonal), undirected.  

Covariance matrix stays the same over time (as does the mean by stationarity; ergodic process).  

![Within person VAR(1) model.](figure/plot of within level VAR model-1.png)


## Longitudinal data, temporal analysis.
![Example VAR(1) process.](figure/unnamed-chunk-1-1.png)

The plot already hints that you could as well view this as being due to one underlying latent process.  

## Longitudinal data, temporal analysis.

Between person covariances are equal, when taking equidistant samples at random timepoints of individuals.  
I.e., we get measurement invariance.  

(Between person) covariance at random time 1, of 10000 samples.  

```
##      [,1] [,2] [,3]
## [1,] 0.25 0.15 0.14
## [2,] 0.15 0.25 0.14
## [3,] 0.14 0.14 0.24
```
Covariance at random time 2 =  time 1 + 100, of 10000 samples.  

```
##      [,1] [,2] [,3]
## [1,] 0.24 0.14 0.14
## [2,] 0.14 0.25 0.14
## [3,] 0.14 0.14 0.24
```
## Longitudinal data, temporal analysis.

AR processes are of course limited in real world scenarios without time-varying parameters.  

![Image of a white noise, true AR(1), and AR(1) model fit to true stepwise linear model.](figure/ARProcess.png)


## Intervention designs.
Definitely different implications for interventions.  

- SI model predicts that effects should be very specific when intervening on symptoms/disorders/illness, but general effects should be expected when intervening on the secondary impairments - i.e., joblessness, basic needs etc. I.e., adding beneficiary resources to the overall pool of problems/resources produces transdiagnostic healing.  

- Network predicts that intervening on a symptom causes changes in adjacent symptoms. In directed networks intervening on a parent node produces healing of nodes that parent has edges to.  

- LP predicts, that most effective intervention is directly aimed at the underlying process. E.g., attacking emotion regulation might be one such process, that then creates (transdiagnostic) healing.  

In practice, the difficulty might be 'fat handedness' [, e.g. Eronen (2020)](https://doi.org/10.1016/j.newideapsych.2020.100785), of interventions.    

## Intervention designs.

![Some differences in how interventions might affect different models. A: Intervention hitting a 'central node' (although this idea has already been critized. B: Most effective intervention for an LV model. C: Conceputalization of a most (transdiagnostically) effective intervetion to SI model. D: Sketch of absolutely fat handed itervention.)](figure/applicationFigure4.png)

- This could've been used in a Bayesian averaging setting, with each of the A-D comprising an individual process through which intervention (e.g., psychotherapy) affects. Then you could've seen if any model is superior in terms of weight in the Bayesian model.  

## Ergodicity

- I suspect that most ergodic data-generating processes won't be able to distinguish between Network, LV or SI models.  

This is because we already know that the cross-sectional data does not distinguish between the models (up to second order moments at least). Then it follows that these processes won't also distinguish between models, when we're inspecting them over time instead. By definition of ergodicity, they are the same processess when observed cross-sectionally over multiple subjects, or observed over time.  


## 
