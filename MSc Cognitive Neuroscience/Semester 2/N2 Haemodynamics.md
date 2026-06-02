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



# 

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

