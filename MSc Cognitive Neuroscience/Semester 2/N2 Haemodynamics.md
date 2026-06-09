### Lecture 5
Echo Planar Imaging 

Beckmann
Using bet to remove non brains areas
Multichannel segementation 
All should be pre registered 

mixture of Gaussians into the histogram of intensities 

Bias field correction
**histograms**
Bias field correction helps with the histogram detections
Segmentations 

Local information - Markof Random Field
Higher b value, more of trusting the neighbour 

PVE - partial volume estimates 

Segmentation of sub structures 

Deformable model, mesh wrapping 
 Jacobian map. Is multiplied with the normal image, this helps to make sure it is standardised after the non linear transformation 

Principle of bold signal 
	
Static dephasing vs. dynamic dephasing 
Depends on the size of the blood vessel 

Distortion field moves as well as



# Notes to keep in mind

EPI spiral makes more sense 
Spin echo, straight makes more sense
Distortions - susceptibility gradient, 
Frequency dependent ofset, sheering effect. 
All appear in the phase encoding direction, dropout. 
Distortion and bandwidth. increase in acquisition bandwidth in read direction, reduces the time it takes. Shorter time between successive datapoints. 
Bandwidth, frequency, gradient, noise
Epi- accelerating epi, limiting factor of this acc, 2d epi - all acc in the slice direction. Multislice. If worried about regions with short t2 star, multiecho, accelerate in plane. 
3d epi for higher resolution, 

DWI not based on relaxation time for uses epi. 

Motion leads to a signal increase, 
Spin history effects increase or decrease signal. Depends on the movement. 
Spin echo, gradient echo across all 4 contrasts
1 question per lecture/topic. 


# Statistic in fMRI
Null hypthesis
	Opposite of what I hope to see
	What happens if there is no effect of your intervention or no difference between the groups
Test statistic 




Family wise error

Using resel-which is the number of voxels that….

Gaussian random field theory 

Cluster based, always reject the null hypothesis of the largest cluster. 


FDR - false discovery rate. 
FDR vs FWE


# Advance linear modelling 
Reasons for finding the temporal derivative. 
Motion parameter confounds - 

Outlier timepoint detection 
	Removes all time points that are outliers, however it does not go away completely, it just doesn’t take into account any influences at that time point. 

ICA denoising 

Location of certain effects - cardiacs in more central areas, respiratory everywhere. 
	Cant be controlled by recording heart rate and respiration. 

Demeaning EVs
	demeaning your regressors 

Parametric designs 
Either different levels of pain are modelled in the same mode, this assumes that lower pain and higher pain is proportional 
Or you can seperate both stimulus and look at them as 2 seperate designs and not in a gradient. 
Linear trend(model it using multiple EVs) vs 

F-contrasts - 

Modelling positive and negative interaction. 

Correlation of EVs

## Beckmann summary class
first few slide minus artefacts
No BET in first exam, bet in the resit
Registration important(Anatomy)
![[Pasted image 20260605085037.png|286]]
Cost function not as important
interpolation important
	pros and cons of types of interpolation
masking, transfering binary masks
Anatomical structre not as important but read,
Probability model iportant 
	Bias field correction
	spatial neighbourhood information!

Optional use of prior important for resit(modelling shape and intensity)

Voxel based analysis of gm volume, what are the disadvantaged
no atrohpy
![[Pasted image 20260605090408.png|243]]

What order should you do motion correction and slice timing(ideally should be done in 4d)
use the temporal derevitive as an additional regressor in the GLM
4d would be good but its numerically challenging but doesnt do much. Instead we skip slice timing and add temporal derivative to GLM (Taylor expansion) to save ourselves the interpolation. The derivative can also account for different Haemodynamic delays within a slice.


intensity normalization - why its done

![[Pasted image 20260605091213.png|240]]
Multiple comparison problems
cluster thresholding
False discovery rate(important)
Group analysis very important
Paired t test
Demeaning
Parametric variation - linear trends
How can you test when you are close to 
Structured noise and GLM

Structural connectivity - diffusion


DTI important