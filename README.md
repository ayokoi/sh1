# Hierarchical representations of motor sequences in the human neocortex: MATLAB codes and related data

This repository contains MATLAB codes and related data to reproduce the Figures (except for schematics) and related analyses in Yokoi and Diedrichsen (2019) paper. Each sub-folders contains codes (`/code`), data (`/data`), and some resultant figures (`/figure`), respectively.

## Dependency
- MATLAB R2015b (consistency is not guaranteed for other versions)
- Dataframe toolbox (https://github.com/jdiedrichsen/dataframe)
- PCM toolbox (https://githum.com/jdiedrichsen/pcm_toolbox)
- Caret (http://brainvis.wustl.edu/wiki/index.php/Caret:About)
 
## Main functions to reproduce figures
Here is a short summary of two main codes for reproducing behabioral and imaging figures. For more detail, see help messages for each function. 
 
|Function |Description |Example |
|----|--------|----|
|sh1_behavior.m |Reproduces figures for behavioral analyses. | `sh1_behavior('Figure_all')` |
|sh1_imaging.m |Reproduces figures for imaging analyses. | `sh1_imaging('Figure4c')` |

- The repository also contains some MATLAB functions which utilizes the Caret for visualization.