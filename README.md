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
- SimpSimScript.R (incomplete, second attempt 2022; now further than first attempt...)
- Replication2022.m corresponding Matlab script (goes as far as the SimpSimScript for comparison). 

Outdated:
- OLD_ATTEMPT_2021Simplified Simuloation Script (incomplete, scruffy attempt from last year; changed to vector lists /arrays 
  which is a lot easier, but failed for b and B_norm, where the old attempt ran correct results for VFE (up to the modified .m script)
- SimpSimScript is as further than my first attempt by now (27.06.2022)

# Next:
SimpSimScript is my second attempt to replicate the Simpliefied_Simulation_Script.m 
- Note that the object chosen_action is modified to all ones, in order to get the loop running at this stage
- Got as far as policy_posterior. cell_md_dot for the Comb input needs external sum() so far (not the most elegant
  solution, but not as bad as with B_norm() - both have to be checked for flexibility on different models based on the script
  in the future. For the ATUT model it works so far.
 
Very basic Active inference Tutorial is completed, just going through some Math from Shannon's paper "A mathematical theory of communication", 
Jensens's inequality, Bayes' rule, Thermodynamics and some Semiotics/Linguistics... Will be released in the near future along with some other stat tuts via BEM... 



 
