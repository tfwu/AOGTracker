# AOGTracker 

## Implementation of the papers
Tianfu Wu , Yang Lu and Song-Chun Zhu, Online Object Tracking, Learning and Parsing with And-Or Graphs, under review for TPAMI (arXiv 1509.08067, a short version appeared in CVPR14).

Written by Matt Tianfu Wu (tfwu@stat.ucla.edu)

## 0. System requirements:

### OS
We tested our tracker on Ubuntu 14.04 LTS. Other OS will be supported later on.

### Third-party libraries
sudo apt-get install build-essential libboost1.55-all-dev libopencv-dev libeigen3-dev libfftw3-dev mpich2

###  The TRAX library for VOT
It is needed for integrating AOGTracker4VOT into vot-toolkit. Please follow https://github.com/votchallenge/trax.git.


## 1. Datasets
TB100/50 is available at http://cvlab.hanyang.ac.kr/tracker_benchmark/. Please download all the data to PATH_TO_YOUR_DATA (e.g. /home/matt/Data/TB100/). Please also download "TB100-occ" (provided by the authors of TB-100) for running TRE. 

VOT datasets are vailable at http://www.votchallenge.net/. vot-toolkit will download the data automatically.

## 2. Run AOGTracker in the termial
change setting in the downloaded configuratin file, tracker_config_release.xml (e.g., specify your data directory, etc)

cd to PATH_TO_YOUR_DOWNLOAD_AOGTracker

./entry Tracking ./tracker_config_release.xml


## 3. Run AOGTrackerMPI in the terminal
If you have mulitple workstations available, please follow https://help.ubuntu.com/community/MpichCluster to set up the cluster.

change setting in the downloaded configuratin file, tracker_config_release.xml (e.g., specify your data directory, etc)

cd to PATH_TO_YOUR_DOWNLOAD_AOGTracker

/usr/bin/mpiexec.mpich -f machine ./entryMPI Tracking ./tracker_config_release.xml

Note: "machine" is a txt file specifying the cluster machines, and the executable and data directory should be shared among all cluster machines.

## 4. Evaluate the performance in TB-100
Run ./evaluation/AOGTracker_TB100_GenPerfMat.m first and then ./evaluation/my_TB100_performance_comp.m to generate the plots. 


## 5. Run AOGTracker4VOT
a) Follow the vot-toolkit tutorial to set up the testing environment using matlab.

b) Modify the configuration.m file by adding to the end: set_global_variable('trax_timeout', 10*60); 

c) Modify the tracker_AOGTracker.m, e.g., 

tracker_label = ['AOGTracker'];

> for VOT2013, VOT2014 and VOT2015
tracker_command = 'your_path_to/AOGTracker4VOT your_path_to/vot_config_release.xml';  % change your_path_to 

> for VOT-TIR2015
tracker_command = 'your_path_to/AOGTracker4VOT your_path_to/vottir_config_release.xml';  % change your_path_to 

