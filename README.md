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
- SimpSimScript now further than my first attempt (27.06.2022)

# Next:
SimpSimScript is my second attempt to replicate the Simpliefied_Simulation_Script.m 
- Finsihed the big loop! Appears to be fine. Only issue appears to be in G_error, related to G_epistemic_value (again...)
- cell_md_dot for the Comb input needs external sum(), but works for Gen_model=1 and 2. Still not a very professional solution for cell_md_dot.
 
Very basic Active inference Tutorial is completed, just going through some Math from Shannon's paper "A mathematical theory of communication", 
Jensens's inequality, Bayes' rule, Thermodynamics and some Semiotics/Linguistics... Will be released in the near future along with some other stat tuts via BEM... 



 
