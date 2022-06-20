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
- I will probably switch to using arrays for the most of the script, but leave it for now (change has a lot of impact on the script, so Ill do it later)
- going well so far, but plotting routines will probably be more difficult to reproduce, will see. Some of the functions can be improved as well. 
 
Very basic Active inference Tutorial is completed, just going through some Math from Shannon's paper "A mathematical theory of communication", 
Jensens's inequality, Bayes' rule, Thermodynamics... Will be released in the near future along with some other stat tuts via BEM... 

# Current struggle: 
 - [[SOLVED]]get plot running for Figure 3 and 1 of Example 2 in the script for message passing example. 
      - the function spm_cat could be translated into R in the future - or some equivalent code written for it... 

# Issues:
- the message passing example and EFE conflicts with each other. Clear environment seems not to be enough, so restarting the 
  session might be necessary for correct results


 
