# Background Studies using **FairShip**

This repository contains scripts that were used to produce the **DIS background studies** presented in *[link to presentation]*.  

If you do not already have FairShip installed, please refer to the [FairShip README](https://github.com/ShipSoft/FairShip?tab=readme-ov-file#build-instructions-using-cvmfs) (Steps 1–4) for detailed installation and setup instructions on **lxplus**.

---

## 0. Preliminary Setup

Before launching large-scale analyses, it is a good idea to verify that all scripts run correctly and that there are no missing dependencies.  
Please install the following Python packages (required for the `sbtveto` module used in the analysis):

```bash
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
pip install scikit-learn dm-tree torch_geometric
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv \
     -f https://data.pyg.org/whl/torch-2.4.0+cpu.html
```

## 1. Setting up the FairShip Environment for HTCondor

Source the setup script (replace \$SHIP_RELEASE with the release you want to use):
```bash
source /cvmfs/ship.cern.ch/$SHIP_RELEASE/setUp.sh
```
To load the FairShip environment for use within HTCondor, execute:
```bash
alienv load FairShip/latest-<version>-release > /afs/cern.ch/user/<path>/<to_your_folder>/HTCondor_scripts/config_<version>.sh
```
Alternatively, you can also use `eval` as described in [Fairship README](https://github.com/ShipSoft/FairShip?tab=readme-ov-file#build-instructions-using-cvmfs) (Step 5).


## 2. Modifying the Submission Scripts

All submission files can be found inside the `condor_scripts` directory.

### 2.1  Submission Files: (`condor_muonDIS_complete.sub`,`condor_neuDIS_complete.sub`)

The example below shows the changes to be done case for muonDIS, but these can be generalised for the other scripts (condor_<>.sub) as well.

#### Email Configuration

Change the `notify_user` field to your respective email address if you would like to be notified of the job status. 

#### Modify Folder Paths:

Replace the folder paths in **line 3** with your own folder paths:

```bash
arguments = $(ClusterId) $(ProcId) $(Runnumber) /eos/user/<path>/<to_your_output_folder>/ /afs/cern.ch/user/<path>/<to_the_cloned_repository>/ muonDIS
```
Also update the input list path in **line 21**:

```bash
queue Runnumber from joblists_muDIS_ECN3_2024.csv
```

if you wish to print the logs, you may uncomment the following lines, but ensure that the necessary directories exist for HTCondor logs, outputs, and errors.

```bash
#error   = /afs/cern.ch/work/<path>/<to_your_Condor_folder>/error/muDIS_$(ClusterId).$(ProcId).err
#log     = /afs/cern.ch/work/<path>/<to_your_Condor_folder>/log/muDIS_$(ClusterId).$(ProcId).log
#output  = /afs/cern.ch/work/<path>/<to_your_Condor_folder>/output/neuDIS_$(ClusterId).$(ProcId).out
```
It is advised to use the afs work space (`/afs/cern.ch/work/`) for these files since the [quota is much larger](https://resources.web.cern.ch/resources/Manage/ListServices.aspx).


### 2.2 Executable (`condor_complete.sh`)

Change the \$SHiP_RELEASE version on **line 13** with the release you want to use):
```bash
source /cvmfs/ship.cern.ch/$SHIP_RELEASE/setUp.sh
```
Change the path to the `config_<version>.sh file on **line 18** with the one produced in [Step 1]():
```bash
source /afs/cern.ch/user/<path>/<to_your_folder>/HTCondor_scripts/config_<version>.sh
```
Replace all mentions of  `root://eospublic.cern.ch` with `root://eosuser.cern.ch` to save the files to your personal EOS space.
---

## 3. Submitting to HTCondor
To submit jobs, simply run:

```bash
condor_submit condor_muonDIS_complete.sub
```
To monitor the job queue, use `condor_q`. 

For more useful commands and adavnced usage, please refer to the [official HTCondor documentation](https://htcondor.readthedocs.io/en/latest/)


## 4. Merging Results

After the Condor jobs finish successfully, the results can be combined using the provided merge submission file.
This step collects and merges the outputs from all individual jobs.

Similar to Step 2, you may configure email notifications and error logs in the same way as for the main submission scripts.

```bash
arguments = /eos/user/<path>/<to_your_output_folder>/ /afs/cern.ch/user/<path>/<to_the_cloned_repository>/ /afs/cern.ch/<path>/<to_directory_for_merged_outputs>/
```

To combine the results from the successful Condor jobs,

launch
```bash
condor_submit merge_all_studies.sub
```



