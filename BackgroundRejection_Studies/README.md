# DIS Background Studies using **FairShip**

This repository contains scripts that were used to produce the DIS background studies presented at the [34th SHiP Collaboration meeting](https://indico.cern.ch/event/1578670/contributions/6689107/attachments/3137319/5568518/Status_MuonDIS_nuDIS.pdf). 

If you do not already have FairShip installed, please refer to the [FairShip README](https://github.com/ShipSoft/FairShip?tab=readme-ov-file#build-instructions-using-cvmfs) (Steps 1–4) for detailed installation and setup instructions on **lxplus**. 
Commit: `d2b67f8318` can be used for the existing DIS simulations in `/eos/experiment/ship/simulation/bkg/`.

---

## 0. Preliminary Setup

Please install the following Python packages (required for the `sbtveto` module used in the analysis):

```bash
python3 -m pip install torch==2.5.1+cpu torchvision torchaudio \
    --index-url https://download.pytorch.org/whl/cpu
python3 -m pip install pyg_lib torch_scatter torch_sparse torch_cluster \
    torch_spline_conv -f https://data.pyg.org/whl/torch-2.5.1+cpu.html
python3 -m pip install torch_geometric scikit-learn dm-tree 
```

## 1. Setting up the FairShip Environment for HTCondor

Source the setup script:
```bash
source /cvmfs/ship.cern.ch/24.10/setUp.sh
```

To load the FairShip environment for use within HTCondor, execute:
```bash
alienv load FairShip/latest-d2b67f8318-release > /afs/cern.ch/user/<path>/<to_your_folder>/config_<version>.sh
```
Alternatively, you can also use `eval` as described in [Fairship README](https://github.com/ShipSoft/FairShip?tab=readme-ov-file#build-instructions-using-cvmfs) (Step 5).


## 2. Modifying the Submission Scripts

After successfully cloning the repository:

```bash
git clone https://github.com/anupama-reghunath/FairShip_Analysis.git
cd BackgroundRejection_Studies
```

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
#output  = /afs/cern.ch/work/<path>/<to_your_Condor_folder>/output/muDIS_$(ClusterId).$(ProcId).out
```
It is advised to use the afs work space (`/afs/cern.ch/work/`) for these files since the [quota is much larger](https://resources.web.cern.ch/resources/Manage/ListServices.aspx).

### 2.2 Executable (`condor_complete.sh`)

Change the path to the `config_<version>.sh file on **line 18** with the one produced in [Step 1]():
```bash
source /afs/cern.ch/user/<path>/<to_your_folder>/config_<version>.sh
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

Update the folder paths in `merge_all_studies.sub`:

```bash
arguments = /eos/user/<path>/<to_your_output_folder>/ /afs/cern.ch/user/<path>/<to_the_cloned_repository>/ /afs/cern.ch/<path>/<to_directory_for_merged_outputs>/
```

and simply run:

```bash
condor_submit merge_all_studies.sub
```
You can also make use of the `consolidate_tables.py` to get the tables in a matching format as to the presentation.



