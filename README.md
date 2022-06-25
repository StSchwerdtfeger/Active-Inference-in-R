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
- SimpSimScript.R (incomplete, second attempt 2022)
- Replication2022.m corresponding Matlab script (goes as far as the SimpSimScript for comparison). 

Outdated:
- OLD_ATTEMPT_2021Simplified Simuloation Script (incomplete, scruffy attempt from last year; changed to vector lists - a lot easier)
- SimpSimScript is as far as my first attempt by now (19.06.2022)

# Next:
SimpSimScript is my second attempt to replicate the Simpliefied_Simulation_Script.m 
- Note that the object chosen_action is modified to all ones, in order to get the loop running at this stage
- Script now runs VFE correctly. Issue was that list objects are not always treated as double()... as.matrix() 
  within the message passing part helped. 
 
Very basic Active inference Tutorial is completed, just going through some Math from Shannon's paper "A mathematical theory of communication", 
Jensens's inequality, Bayes' rule, Thermodynamics... Will be released in the near future along with some other stat tuts via BEM... 



 
