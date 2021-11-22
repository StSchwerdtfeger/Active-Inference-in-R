# Active-Inference-in-R

 R conversion of matlab code from the supplemantary material
 of the paper "A Tutorial on Active Inference Modelling and its Application to Empirical Data" 
 by Ryan Smith, Karl Friston and Christopher J. Whyte.

 Link to paper: https://psyarxiv.com/b4jm6/
 
 Original Matlab code also available here (Ryan Smith's github):  https://github.com/rssmith33/Active-Inference-Tutorial-Scripts

 Conversion to R by Steffen Schwerdtfeger

 I am still new to all of this. Feedback or contributions always welcome.
 
 
# CONVERTED SO FAR:

- Pencil and Paper Example
- Prediction Error Example
- Message Passing Example
- EFE Precision Updating Example 

# Current struggle: 
 - [[SOLVED]]get plot running for Figure 3 and 1 of Example 2 in the script for message passing example. 
      - the function spm_cat could be translated into R in the future - or some equivalent code written for it... 

# Issues:
- the message passing example and EFE conflicts with each other. Clear environment seems not to be enough, so restarting the 
  session might be necessary for correct results

# Next:
- Simpl. Simulation script. Around 85% converted so far. Got stuck with the functions smp_cross + G_epistemic_value but will 
  get back to it soon again, as I nearly solved it... 
