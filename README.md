---
title: Plasticity of the gastrocnemius elastic system in response to decreased work
  and power demand during growth
author: "Zanne Cox"
date: '2023-01-29'
output:
  pdf_document: default
  html_document: default
---
## Plasticity of the gastrocnemius elastic system in response to decreased work and power demand during growth



####  Abstract

Elastic energy storage and release can enhance performance
that would otherwise be limited by the force–velocity constraints of
muscle. Although functional influence of a biological spring depends
on tuning between components of an elastic system (the muscle,
spring-driven mass and lever system), we do not know whether elastic
systems systematically adapt to functional demand. To test whether
altering work and power generation during maturation alters the
morphology of an elastic system, we prevented growing guinea
fowl (Numida meleagris) from jumping. We compared the jump
performance of our treatment group at maturity with that of controls
and measured the morphology of the gastrocnemius elastic system.
We found that restricted birds jumped with lower jump power and
work, yet there were no significant between-group differences in the
components of the elastic system. Further, subject-specific models
revealed no difference in energy storage capacity between groups,
though energy storage was most sensitive to variations in muscle
properties (most significantly operating length and least dependent
on tendon stiffness). We conclude that the gastrocnemius elastic
system in the guinea fowl displays little to no plastic response to
decreased demand during growth and hypothesize that neural
plasticity may explain performance variation


###  Data Analysis pipeline

1. **Collect and Analyze Force data from jumping tests of 16 Guinea Fowl**

    + As described previously (Cox et al., 2020), at skeletal maturity
(between 29 and 31 weeks old) jump performance was measured by
placing each bird in turn on 6×6 inch (15.24×15.24 cm) force plates
(AMTI HE6x6; Watertown, MA, USA) enclosed in a tapered box
and encouraging the birds to jump
    +  Calculated:
      1. Jump power from instantaneous net vertical ground reaction and the vertical center of mass velocity
      2. Instantaneous Velocity: integrating center of mass acceleration
      3. Jump work: integrating power with respec to time

2. **Build subject Specific Musculoskeletal models of individuals**
    + Collect Morphology Measurements
    
        1. Muscle Analysis of lateral and medial gastroc
        2. Moment Arm of Achillies about the ankle
        3. Tendon force/Length 
        4. Tendon Slack Length
        5. Leg bone segment lengths
        
    + Fit model parameters to match experimentally measured values via particle swarm optimization

3. **Simulate muscle activation and tendon strain across a range of postures**
    + Extract Elastic Energy Storage
  
4. **Statistical Tests**
1. *Do components of the gastrocnemius elastic system change systematically in response to changes in power and work demand during growth?*

    + To determine whether components of the gastrocnemius elastic
system change systematically in response to changes in demand,
we evaluated the influence of treatment group (restricted versus
control) on each element of morphology measured. This was
accomplished using t-tests if the homogeneity of variance
assumption test was passed, and using a Kruskal–Wallis test by
ranks when this criterion was not met. The threshold for statistical
significance was set at 0.005 after a Bonferroni correction for
multiple comparisons.

2.  *Is the energy storage capacity reduced in individuals that
did not jump during growth?*

    + The relationship between treatment group and elastic energy storage capacity was evaluated with a t-test after data passed tests for normality and homogeneity of variance, as described above for evaluation of
differences between groups of individual elastic system components


3.  *Which type of morphological variation has the greatest
influence on energy storage capacity?*

    + We used stepwise comparison of Akaike information criterion
(AIC) values (stepAIC R Mass package; Venables and Ripley, 2002) to determine the parameters and coefficients of the full model that best predicted elastic energy storage potential across natural variation of joint postures in preparation for jumps. The full statistical model evaluated included stored strain energy (PE) as a dependent factor and, as potential independent variables, tendon stiffness (tendonK), the summed maximum isometric force capacity of LG and MG along the tendon (sumFMax), the average LG and MG optimal fascicle length (avOFL) and starting muscle length.  We included possible interaction terms between muscle force capacity, tendon stiffness and muscle start length (sumMaxF×tendonK×avLenA0c) and between optimal fascicle length and tendon stiffness (avOFL×tendonK) following recommendations by Zajac (1989) of functional equivalent muscle tendon joint properties at zero activation of the LG and MG in the pre-jump posture (avLenA0c). 

4.  *Does elastic energy storage capacity predict peak jump
powers and work?*
    + The relationship between Achilles tendon elastic energy storage
capacity and experimentally measured muscle-mass-normalized peak power output and jump work were both tested with a linear model with elastic energy storage as the dependent variable and peak power or jump work as the independent variable.






