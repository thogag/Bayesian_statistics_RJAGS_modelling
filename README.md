---
title: "Bayesian mini-project: Mixture model for alternative splicing"
author: "Thomas Gagnieu"
date: "2025-11-05"
output: html_document
---

# Bayesian Mini-Project: Mixture Model for Alternative Splicing

This project implements a **Bayesian mixture model** to study the hypothesis of *transcriptional noise* in alternative splicing (Saudemont et al., *Genome Biology*, 2017).  
The goal is to model the relationship between **gene expression** and **alternative splicing rate**, and to distinguish between **functional** and **non-functional (noisy)** splicing events.

---

## Biological Context

According to Saudemont *et al.* (2017), a large fraction of observed alternative splicing (AS) might result from **non-adaptive transcriptional noise** rather than functional regulation.  
Under this hypothesis:

- The fraction of *minor isoforms* (AS rate) should be **higher in weakly expressed genes**.
- Highly expressed genes should have **lower splicing noise** due to stronger quality control mechanisms (e.g., NMD).

The aim of this Bayesian project is to **model this dependency** and to **quantify the relative contribution** of noise vs. functional splicing.

---

## Model Overview

A **Bayesian mixture model** is fitted to the binomial counts *(n₂ successes among n₁+n₂ trials)* per intron.

| Component | Meaning | Parameters |
|------------|----------|------------|
| **Noise (θ_b)** | Random splicing events (transcriptional noise) | logit(θ_b) = β₀ + β₁·expr + β₂·introns |
| **Functional (θ_f)** | Regulated, biologically relevant AS | θ_f ~ Beta(a, b) |
| **p** | Proportion of functional AS among all introns | p ~ Uniform(0,1) |

---

## Code Summary

### 1. Data Preparation
```r
expr_df <- read.csv("data/gene_data.csv", sep = ";")
splice_df <- read.csv("data/intron_data_v2.csv", sep = ";")

df <- splice_df %>%
  left_join(expr_df %>% select(gene_id, weighted_fpkm, type), by = "gene_id") %>%
  mutate(
    y = n2, n = n1 + n2,
    expr = log1p(weighted_fpkm),
    nmd = ifelse(type == "pseudogene", 1, 0)
  ) %>%
  filter(n > 0, !is.na(expr), y <= n)
  ``` 
A subsample of 1000 introns is used for demonstration.

### 2. Bayesian Mixture Model (JAGS)
```{r}
model_string_mix <- "
model {
  for (i in 1:N) {
    y[i] ~ dbin(theta[i], n[i])
    Z[i] ~ dbern(p)
    theta_f[i] ~ dbeta(a, b)
    logit(theta_b[i]) <- beta0 + beta1*expr[i] + beta2*introns[i]
    theta[i] <- Z[i]*theta_f[i] + (1 - Z[i])*theta_b[i]
  }

  beta0 ~ dnorm(0, 0.001)
  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0, 0.001)
  p ~ dunif(0, 1)
  a ~ dunif(0, 10)
  b ~ dunif(0, 10)
}
```
Sampling is performed using rjags and coda over 3 chains with convergence diagnostics (Gelman-Rubin, autocorrelation).

### 3. Posterior Estimates
Key results:

| Parameter        | Median           | 95% CI                                         | Interpretation |
| ---------------- | ---------------- | ---------------------------------------------- | -------------- |
| β₀ = -7.14       | [-7.49, -6.84]   | Baseline log-odds of splicing noise (very low) |                |
| β₁ = -0.66       | [-0.87, -0.47]   | Negative effect of expression on noise         |                |
| β₂ = +0.29       | [+0.10, +0.46]   | Positive effect of intron count                |                |
| p = 0.78         | [0.72, 0.83]     | 78% of introns classified as functional        |                |
| a,b = 0.11, 3.11 | → E[θ_f] ≈ 0.034 | Mean splicing rate for functional AS           |                |

### 4. Posterior Classification
Each intron is assigned a posterior probability of being functional, P(Zi=1∣data),
then categorized as:

| Class                  | Criterion     | Meaning         |
| ---------------------- | ------------- | --------------- |
| Very likely functional | P ≥ 0.9       | Regulated AS    |
| Likely functional      | 0.5 ≤ P < 0.9 | Intermediate    |
| Uncertain              | 0.1 ≤ P < 0.5 | Ambiguous       |
| Likely noise           | P < 0.1       | Random splicing |

### 5. Posterior Predictive Checks
+ Global PPC: Observed splicing rate (2.5%) lies within model’s simulated 95% interval [1.6%, 4.8%].
+ Class-level PPC:
Functional: Model reproduces 4% observed rate.
Noise: Model overestimates slightly (predicted ~2% vs observed 0.04%), suggesting an unmodelled zero-excess component.

### 6. Marginal Effects
| Component            | Behavior                  | Interpretation                                   |
| -------------------- | ------------------------- | ------------------------------------------------ |
| **Functional (θ_f)** | Flat (~3–4%)              | Regulated AS independent of expression           |
| **Noise (θ_b)**      | Decreases with expression | Strong quality control in highly expressed genes |

Graphs show that as expression increases:
+ Noise probability drops from ~0.006 to <0.001.
+ Functional AS remains constant.


## Empirical Comparison
| Component  | total_y | total_n   | Proportion |
| ---------- | ------- | --------- | ---------- |
| Functional | 109,156 | 2,694,746 | 0.0405     |
| Noise      | 709     | 1,721,394 | 0.00041    |

→ Functional AS ≈ 100× more frequent than noise.


## Interpretation Summary
| Aspect                    | Result   | Biological Meaning                      |
| ------------------------- | -------- | --------------------------------------- |
| Functional proportion (p) | 0.78     | Majority of AS is functional            |
| Mean θ_f                  | 0.034    | Stable, regulated splicing              |
| Mean θ_b                  | <0.001   | Rare, random AS                         |
| Expression effect         | Negative | Fewer noisy events at high expression   |
| Introns effect            | Positive | More introns → more noise opportunities |
| Correlation (obs vs pred) | ~0.99    | Excellent model fit                     |
| R² (Bayesian)             | ~0.46    | Moderate uncertainty across introns     |


## Biological Insight
The Bayesian mixture model confirms that alternative splicing arises from two distinct regimes:
+ A functional, expression-independent component (~3–4% rate)
+ A stochastic noise component, strongly decreasing with expression

This supports the hypothesis that most observed alternative splicing events are biologically functional, and that splicing noise is rare and tightly controlled in highly expressed genes.

If you want to see all the results in html format, click [here](https://thogag.github.io/Bayesian_statistics_RJAGS_modelling/Results.html)


## Dependencies
R packages: `dplyr`, `ggplot2`, `rjags`, `coda`, `bayesplot`
Input data:
+ `data/gene_data.csv`
+ `data/intron_data_v2.csv`


### References
Saudemont B. et al.,
Nonadaptive Alternative Splicing Is the Main Determinant of the Transcriptome Evolution in Fungi.
Genome Biology (2017) 18:208.
DOI: [10.1186/s13059-017-1344-6](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1344-6)


Author: Thomas Gagnieu
Date: 05 November 2025
