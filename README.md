# Active-Inference-in-R

 R conversion of matlab code from the supplemantary material
 of the paper "A Tutorial on Active Inference Modelling and its Application to Empirical Data" 
 by Ryan Smith, Karl Friston and Christopher J. Whyte.

 Link to paper: https://www.sciencedirect.com/science/article/pii/S0022249621000973
 
 Original Matlab code also available here (Ryan Smith's github):  https://github.com/rssmith33/Active-Inference-Tutorial-Scripts

 Conversion to R by Steffen Schwerdtfeger

 Feedback or contributions always welcome. Still fairly new to R on that level of complexity.
 
 
# CONVERTED SO FAR:

- Pencil and Paper Example
- Prediction Error Example
- Message Passing Example
- EFE Precision Updating Example 
- Simplified Simuloation Script (incomplete, scruffy attempt from last year; changed to vector lists - a lot easier)
- SimpSimScript.R (incomplete, second attempt 2022)
- Replication2022.m corresponding Matlab script (goes as far as the SimpSimScript for comparison). 

# Next:
SimpSimScript is my second attempt to replicate the Simpliefied_Simulation_Script.m 
- Status: Got to the loop in 2 days, so looking good so far. ones() and zeros() from the pracma package was used, but could replaced.
- got to VFE, now working on replications of G_epistemic_value() and spm_cross
- note that the object chosen_action is modified to all ones, in order to get the loop running at this stage
 
Very basic Active inference Tutorial is completed, just going through some Math from Shannon's paper "A mathematical theory of communication", 
Jensens's inequality, Bayes' rule, Thermodynamics... Will be released in the near future along with some other stat tuts via BEM... 

# Current struggle: 
 - [[SOLVED]]get plot running for Figure 3 and 1 of Example 2 in the script for message passing example. 
      - the function spm_cat could be translated into R in the future - or some equivalent code written for it... 

# Issues:
- the message passing example and EFE conflicts with each other. Clear environment seems not to be enough, so restarting the 
  session might be necessary for correct results


 
