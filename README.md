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
- SimpSimScript.R (completed, except of spm plotting routines)

# Next:
SimpSimScript is my second attempt to replicate the Simpliefied_Simulation_Script.m 
- Script is done!! Just plotting is missing (spm plotting routines). I have replicated some of them in other scripts - hope it wont be to tough. 
- cell_md_dot for the Comb input needs external sum(), but works for Gen_model=1 and 2. Still not a very professional solution for cell_md_dot.
- The functions B_norm, spm_KL_dir also need some improvement, to be flexible for other POMDP inputs. 
 
I will try to dig into scripts such as spm_VB_X and look into the way it was carried out in the pymdp package (and will go through the great ActInf Lab Sessions 
on the ATUT material). The learning curve might be a bit steep again, but I will definitly give it a shot right after I finished the plotting of the SimpSimScript, 
as it is also a great opportunity to learn about structuring and writing packages with generalized and flexible functions and intuitive handling (interface
qualities). The pymdp approach also appears to avoid too much code complexity all at once and may be easier to achieve (also refers to recommendations I got when 
translating the SimpSimScript). I will see.  

# Tutorial for some basic math/physics (information theory) of active inference finished: https://journal.medicine.berlinexchange.de/pub/bi7ccdsa 



 
