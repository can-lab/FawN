# FawN
FSL analysis with NiPype

by [Florian Krause](https://www.floriankrause.org/)


## Introduction
FawN is a collection of NiPype workflows for building FSL-style fMRI analyisis pipelines.

Currently available workflows:

* **"resampling"**  
  Resample images into size of reference image
* **"timecourse_extraction"**  
  Extract timecourses from binary masks
* **"first_level"**  
  Analysis of single functional runs
* **"session_level"**  
  Analysis across all functional runs (average) of a single session (convenience workflow)
* **"higher_level"**  
  Analysis across runs/sessions/subject
* **"thresholding"**  
  Thresholding of results on voxel-level (FWE corrected) and cluster-level
  
The workflows expect preprocessed images (see also https://github.com/can-lab/finish-the-job).

## Prerequisites
1. Install [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/)
2. Install nipype with
   ```
   pip3 install nipype fslpy
   ```
3. Download [FawN](https://github.com/can-lab/FawN/master.zip)
4. Install with
   ```
   pip3 install FawN-X.X.X.zip
   ```
   (replace X.X.X with latest release version)

### Donders cluster
If you are working on the compute cluster of the Donders Institute, please follow the following steps:
1. Load Anaconda3 module by running command: `module load anaconda3`
2. Create new environment in home directory by running command: `cd && python3 -m venv fawn_env`
4. Activate new environment by running command: `source fawn_env/bin/activate`
5. Install Nipype into environment by running command: `pip3 install nipype fslpy`
6. Download [FawN](https://github.com/can-lab/fawn/master.zip)
7. Install with
   ```
   pip3 install FawN-X.X.X.zip
   ```
   (replace X.X.X with latest release version)
   
## Usage

See [examples](https://github.com/can-lab/FawN/tree/master/examples).

### Donders cluster
If you are working on the compute cluster of the Donders Institute, please follow the following steps:
1. Start a new interactive job by running command: `qsub -I -l 'procs=8, mem=64gb, walltime=24:00:00'`
2. Load Anaconda3 module by running command: `module load anaconda3`
3. Activate environment by running command: `source fawn_env/bin/activate`
4. Write script `mystudy_fawn.py` (see [examples](https://github.com/can-lab/FawN/tree/master/examples))
5. Run script by running command: `python3 mystudy_fawn.py`
