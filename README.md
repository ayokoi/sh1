# Hierarchical representations of motor sequences in the human neocortex: MATLAB codes and related data

This repository contains MATLAB codes and related data to reproduce the Figures (except for schematics) and related analyses in Yokoi and Diedrichsen (2019) paper. Each sub-folders contains data (`/data`), codes (`/code`), and some resultant figures (`/figure`), respectively.

## Dependency
- MATLAB R2015b (consistency is not guaranteed for other versions)
- Dataframe toolbox (https://github.com/jdiedrichsen/dataframe)
- PCM toolbox (https://githum.com/jdiedrichsen/pcm_toolbox)
- Caret (http://brainvis.wustl.edu/wiki/index.php/Caret:About)
 
## Main functions to reproduce figures
Here is a short summary of two main codes reproducing behabioral and imaging figures. For more detail, see help messages for each function. 
 
|Function |Description |Example |
|sh1_behabior.m |Reproduces figures for behaioral analyses. | `sh1_behabior('Figure_all')` |
|sh1_imaging.m |Reproduces figures for imaging analyses. | `sh1_imaging('Figure4c')` |

- The repository also contains some MATLAB functions which utilizes Caret for visualization.