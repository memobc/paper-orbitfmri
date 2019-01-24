# The Paper

This repository includes scripts and data for the following paper:

[**Cooper, R. & Ritchey, M. (preprint). Cortico-hippocampal network connections support the multidimensional quality of episodic memory.**](https://www.biorxiv.org/content/10.1101/526657v1?rss=1)

# Abstract

Episodic memories reflect a bound representation of multimodal features that can be reinstated with varying levels of precision. Yet little is known about how brain networks involved in memory, including the hippocampus and posterior-medial (PM) and anterior-temporal (AT) cortical systems, functionally interact to support the quality and the content of recollection. Participants learned color, spatial, and emotion associations of objects, later reconstructing the visual features using a continuous color spectrum and 360° panorama scenes. Behaviorally, dependencies in memory were observed for the gist but not precision of these event associations. Supporting this integration, hippocampus, AT, and PM regions showed increased inter-network connectivity and reduced modularity during retrieval compared to encoding. These network connections, particularly to hippocampus, tracked a multidimensional, continuous measure of objective memory quality. Moreover, distinct patterns of connectivity tracked item color precision and spatial memory precision. These findings demonstrate not only how hippocampal-cortical connections reconfigure during episodic retrieval, but how such dynamic interactions might flexibly support the multidimensional quality of remembered events.

# Resources

Psychtoolbox task scripts are included in the `task` folder. The stimuli used in the experiment can be obtained from the following links: [objects](https://bradylab.ucsd.edu/stimuli.html), [panorama scenes](http://people.csail.mit.edu/jxiao/SUN360/index.html), and [sounds](https://csea.phhp.ufl.edu/media/iadsmessage.html). 

I have also shared a few key analysis scripts in the `analysis` folder along with some corresponding `data` files and `reports`.

The general flow of the included analysis scripts is as follows:
- Analyze **behavioral data**: `Orbit-fMRI-Behavior_Paper.Rmd`
    - The formatted report with code can be found [here](http://www.thememolab.org/paper-orbitfmri/reports/Orbit-fMRI-Behavior_Paper.nb.html). This contains all analysis output from behavioral data in `Behavioral_data.csv`.
- Analyze **univariate data**: `Orbit-fMRI-Univariate_Paper.Rmd`. This script analyzes first level betas, reflecting the change in mean ROI activity with increasing memory quality.
    - The analysis output and code can be found [here](http://www.thememolab.org/paper-orbitfmri/reports/Orbit-fMRI-Univariate_Paper.nb.html).
- Analyze **background connectivity data**: 
    - Run the first level analysis using the [CONN toobox](https://sites.google.com/view/conn/): `conn_batch_firstlevel_background.m`. This script requires that all task regressors have already been generated.
    - Analyze the first level connectivity data: `Orbit-fMRI-BackgroundConnectivity_Paper.Rmd`. This script calls functions in `background_functions_paper.R` to format, analyze, and visualize the connectivity matrices at the group level. 
    - The analysis output and code can be found [here](http://www.thememolab.org/paper-orbitfmri/reports/Orbit-fMRI-BackgroundConnectivity_Paper.nb.html).
    - First level ROIxROI connectivity matrices for encoding and retrieval tasks can be found in `Background_connectivity_data.RData`. 
- Analyze **memory-modulated connectivity (gPPI) data**: 
    - Run the first level analysis using CONN: `conn_batch_firstlevel_memorygPPI.m`. This script requires that all task regressors have already been generated.
    - Analyze the first level connectivity data: `Orbit-fMRI-MemorygPPI_Paper.Rmd`. This script calls functions in `gPPI_functions_paper.R` to format, analyze, and visualize the network and ROI connectivity matrices at the group level.
    - The analysis output and code can be found [here](http://www.thememolab.org/paper-orbitfmri/reports/Orbit-fMRI-MemorygPPI_Paper.nb.html), which tests changes in connectivity with i) *overall memory quality*, ii) *color memory precision*, and iii) *scene memory precision*. 
    - First level ROIxROI connectivity matrices (beta estimates for the PPI regressor) for each memory modulator can be found in `Memory_gPPI_data.RData`. 
    - The group-level results of the hippocampus seed to voxel analysis (change in whole-brain hippocampal connectivity with increasing memory quality) are also provided as an spmT.nii file: `HippSeed_wholebrain_MemoryQuality_spmT.nii`.

# Comments?

Please direct any comments to Rose Cooper, rose.cooper at bc.edu. Please feel free to use any of these scripts. Unfortunately I cannot provide support for you to adapt them to your own data. Notice a bug? Please tell me.
