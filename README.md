# Active-Inference-in-R

 R conversion of matlab code from the supplemantary material
 of the paper "A Tutorial on Active Inference Modelling and its Application to Empirical Data" 
 by Ryan Smith, Karl Friston and Christopher J. Whyte.

 Link to paper: https://www.sciencedirect.com/science/article/pii/S0022249621000973
 
 Original Matlab code also available here (Ryan Smith's github):  https://github.com/rssmith33/Active-Inference-Tutorial-Scripts

 Conversion to R by Steffen Schwerdtfeger

 I am still new to all of this. Feedback or contributions always welcome.
 
 
# CONVERTED SO FAR:

- Pencil and Paper Example
- Prediction Error Example
- Message Passing Example
- EFE Precision Updating Example 
- Simplified Simuloation Script (incomplete, see below for details)

# Current struggle: 
 - [[SOLVED]]get plot running for Figure 3 and 1 of Example 2 in the script for message passing example. 
      - the function spm_cat could be translated into R in the future - or some equivalent code written for it... 

# Issues:
- the message passing example and EFE conflicts with each other. Clear environment seems not to be enough, so restarting the 
  session might be necessary for correct results

# Next:
- Simpl. Simulation script. Around 87,5% converted so far. I uploaded the incomplete script. I think I really got far, and
  learned a lot from doing this. Currently working on the md_dot function, but I am busy until end of January and would like 
  to finish this together with other people and focus more on playing around with the actual matlab scripts. 
  The R script is in general very slow and alot of things could be optimized, probably by using functionals from the apply family
  instead of loops... Other things are concerning functions that could be fit together (spm_cross1 and 2). I was mostly interested 
  to get the code runnoing, i.e. creating equivalent results to the matlab code. 
  
  UPDATE: Note that you can run the script up to the current point. The object chosen_action has been altered to make the loop do so. 
  This of course just delivers abstract results, but it can be compared with a in parallel prepared matlab script. This helped me
  to fully understand what the results of some of the parts are. Some problems that I have encountered were hard to track, as it sometimes
  just involved the class() of an object that was changed by some function (from packages such as pracma). Packages in R that deliver 
  matlab functions are oft just equipped with the basic input possibilities, so other functions had to be found. E.g., I did not find a 
  reasonable solution to cover bsxfun... I will write a detailed status report on remaining problems in the near future.
  
# Roughly missing in simp. sim. script so far: 
  Line 351-516 of the big nested loop, the respective functions that are needed for the missing part of the loop, as well as all the spm 
  plotting functions that are used in the script - most of rest is done, just needs some cosmetics, to make it run faster.
  
  
- Information theory for active inference tutorial. I want to put together all the code and tutorials on e.g., mutual information,
  Jensen's inequality etc. that I found for R so far. 
