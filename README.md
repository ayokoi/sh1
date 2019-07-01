# Neural Organization of Hierarchical Motor Sequence Representations in the Human Neocortex: MATLAB codes and related data

This repository contains MATLAB codes and related data to reproduce the Figures (except for schematics) and related analyses in Yokoi and Diedrichsen (2019) paper "Neural organization of hierarchical motor sequence representations in the human neocortex" published in Neuron (doi here). Each sub-folders contains codes (`/code`), data (`/data`), and some resultant figures (`/figure`), respectively. To make the codes work, all the contents (i.e., `/code`, `/data`, and `/figure`) of this repository need to be downloaded.

## Dependency
The following softwares/toolboxes need to be correctly installed to run the present codes.
- MATLAB R2015b (consistency is not guaranteed for other versions)
- SPM12 (for `spm_BMS.m`, https://www.fil.ion.ucl.ac.uk/spm)
- Dataframe toolbox (https://github.com/jdiedrichsen/dataframe)
- PCM toolbox (https://githum.com/jdiedrichsen/pcm_toolbox)
- Caret (http://brainvis.wustl.edu/wiki/index.php/Caret:About)

Note also that I haven't tested if the code correctly works on PCs.

## Main functions to reproduce figures
Here is a short summary of two main codes for reproducing behabioral and imaging figures. For more detail, see help messages for each function. 
 
|Function |Description |Example |
|----|--------|----|
|sh1_behavior.m |Reproduces figures for behavioral analyses. | `sh1_behavior('Figure_all')` |
|sh1_imaging.m |Reproduces figures for imaging analyses. | `sh1_imaging('Figure4c')` |

- The repository also contains some MATLAB functions that utilize the Caret functionalities (`/code/helper`) in the background for visualization of imaging anlayses on flattened cortical surface (fsaverage_sym template).
