# Electrophysiology

# MRI 

# fMRI Analysis

## Preprocessing

What is fMRI data, grey scale pixel intensity data
DICOMs
4dimensional data. 3 Dimensions of space in which voxel/part of the brain followed by a time dimension.
##### Steps of preprocessing
* Reconstruction  
	Converts k-space data into images that we see. 
	Corrects in case of any errors during the scanning process
* Registration  
	fd
* Motion correction  
* Slice timing correction  
* Spatial filtering  
* Temporal filtering  
* Denoising


## fMRI Analysis
### Single session analysis

![[Pasted image 20251124103852.png]]

Model is a mixture of the HRF and an instantaneous response 
Compare the modelled/predicted response with the actual response to check for correlations

Different beta for every single voxel




### Design types 
![[Pasted image 20251124113847.png]]
##### Block model 
is good for detecting the activity but doesnt state anything about the underlying estimation. 
Not ideal for single event analysis

#### Event related design
Inter stimulus interval is required
Longer than 8s ideally


#### Jittered design
having mini spaces in between helps with the processing better
Need to ensure no multicolinearity
Including null trials
Allowing the HRF to return to baseline

## Contrast


Popular frame used here is the Tripartite network
[[Triple network theory]]

t test on the correlation on each voxel
## Brain Stimulation
