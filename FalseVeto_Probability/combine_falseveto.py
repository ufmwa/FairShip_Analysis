import os
import csv
from tabulate import tabulate
import sys

sys.path.insert(1, '/afs/cern.ch/user/a/anupamar/Analysis/Tools')# caution: path[0] is reserved for script path (or '' in REPL)


#threshold_list=[0,10,20,30,45,60,90]
Summation={}
path='/afs/cern.ch/work/a/anupamar/FalseVetoProbability/'#100batch'

threshold_probabilities = {}  # { threshold: [prob1, prob2, ...] }

for job_folder in os.listdir(path):
    
    print(job_folder)
    folder_index=job_folder.split('_')[-1]
    filename=f"{path}/{job_folder}/falseveto_results_{folder_index}.csv"
    
    try:
        with open(filename, 'r') as file:
            reader=csv.DictReader(file, delimiter=',')
            for line in reader:
                
                threshold=int(line[ 'Threshold'])
                    
                if threshold not in Summation:
                    Summation[threshold]={'nEvents checked':0,'nEvents vetoed':0}
                    threshold_probabilities[threshold] = []


                checked = float(line['nEvents checked'])
                vetoed = float(line['nEvents vetoed'])

                Summation[threshold]['nEvents checked']+=checked#float(line['nEvents checked'])
                Summation[threshold]['nEvents vetoed']+=vetoed#float(line['nEvents vetoed'])

                prob = (vetoed / checked * 100) if checked != 0 else 0
                threshold_probabilities[threshold].append(prob)
    
    except Exception as e:
        print(e)            


print('-----------------------------------------SUMMARY OF ANALYSIS-----------------------------------------\n')
tabulate_header = ['Threshold','nEvents checked','nEvents vetoed','Probability(%)']
table_data=[]

for threshold in Summation:        
        
        Probability=Summation[threshold]['nEvents vetoed']/Summation[threshold]['nEvents checked']
        table_data.append([threshold,Summation[threshold]['nEvents checked'],Summation[threshold]['nEvents vetoed'],Probability*100])

print(tabulate(table_data,headers=tabulate_header,tablefmt="pretty"))
#print(tabulate(table_data, headers=tabulate_header, tablefmt="pretty"))
##########################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
from matplotlib.gridspec import GridSpec


# Use the path you're inserting into sys.path
tools_dir = Path('/afs/cern.ch/user/a/anupamar/Analysis/Tools')
style_path = tools_dir / 'thesis.mplstyle'

plt.style.use(str(style_path))
"""

thresholds = sorted(threshold_probabilities.keys())
fig, axs = plt.subplots(2, 4, figsize=(18, 10), constrained_layout=True)
axs = axs.flatten()

for i, threshold in enumerate(thresholds):
    probs = threshold_probabilities[threshold]
    cumulative_avg = np.cumsum(probs) / np.arange(1, len(probs) + 1)

    axs[i].plot(range(len(probs)), probs, marker='.', label='Iteration Probability')
    axs[i].plot(range(len(cumulative_avg)), cumulative_avg, 'r--', label='Cumulative Average')

    axs[i].set_title(f'Threshold@{threshold} MeV', fontsize=12)
    axs[i].set_xlabel('Iteration Index', labelpad=0)  # closer
    axs[i].set_ylabel('False Veto Probability (%)', labelpad=0)
    axs[i].grid(True)


# Use the last subplot as legend container
legend_ax = axs[len(thresholds)]
legend_ax.axis('off')  # Hide the plot frame

# Create custom legend
legend_elements = [
    Line2D([0], [0], color='blue', marker='.', label='Iteration Probability'),
    Line2D([0], [0], color='red', linestyle='--', label='Cumulative Average')
]

legend_ax.legend(handles=legend_elements, loc='center', fontsize=12)

# Title on top, layout adjusted to fit both
fig.suptitle('False Veto Probability Across Iterations for Different Thresholds', fontsize=16)

plt.show()  # Uncomment to view inline
"""


# Set global style
#plt.rcParams.update({'font.size': 18})
#plt.rc('font', family='serif')
# Prepare x and y data
x = sorted(threshold_probabilities.keys())  # Thresholds
y = [np.cumsum(threshold_probabilities[t])[-1] / len(threshold_probabilities[t]) for t in x]  # Final cumulative average

# Create plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(1, 1, 1)

# Plot line and points
plt.plot(x, y, color='0.50', ls='dashed')
plt.scatter(x, y, s=200, color='red', marker='s')

# Annotate each point with its value
for i in range(len(x)):
    plt.text(x[i] + 1, y[i] + 1, f"{y[i]:.1f}%", fontsize=15)

# Axes labels and limits
ax.set_xlabel("Energy threshold for SBT cell (MeV)")
ax.set_ylabel("False veto probability (%)")
ax.set_xlim(-5, 105)
ax.set_ylim(30, 105)
ax.set_xticks(x)

# Save to file
plt.savefig("falsevetoprobability.svg")
#####################

thresholds = sorted(threshold_probabilities.keys())

fig = plt.figure(figsize=(8, 14))
gs = GridSpec(len(thresholds), 1, figure=fig, hspace=0)  # hspace=0 = no vertical gaps

axs = [fig.add_subplot(gs[i, 0]) for i in range(len(thresholds))]

for i, (ax, threshold) in enumerate(zip(axs, thresholds)):
    probs = threshold_probabilities[threshold]
    cumulative_avg = np.cumsum(probs) / np.arange(1, len(probs) + 1)

    ax.plot(probs, marker='.', linestyle='-', color='gray', alpha=0.6)
    ax.plot(cumulative_avg, 'r--', linewidth=2)

    # Internal title
    ax.text(
        0.05, 0.95, f"@{threshold} MeV",
        transform=ax.transAxes,
        fontsize=11, ha='left', va='top',
        bbox=dict(facecolor='none', edgecolor='none', alpha=0.8)
    )

    
    #ax.set_yticks(np.arange(0, 111, 10))  # up to 110 for padding
    ax.set_ylim(-20, 140)
    ax.grid(True, axis='both', linestyle=':', linewidth=0.5)
    # Show ticks, but no y-labels
    ax.set_ylabel('')
    ax.tick_params(labelleft=True, labelsize=9)

    # Remove x-ticks from all except bottom
    if i != len(thresholds) - 1:
        ax.set_xticklabels([])

# Global x and y labels
fig.text(0.04, 0.5, 'False Veto Probability (%)', va='center', rotation='vertical', fontsize=13)
fig.text(0.5, 0.04, 'Iteration Index', ha='center', fontsize=13)

# Global legend
legend_elements = [
    Line2D([0], [0], color='gray', marker='.', label='Iteration Probability', alpha=0.6),
    Line2D([0], [0], color='red', linestyle='--', label='Cumulative Average')
]
fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=2, fontsize=11, frameon=False)

plt.show()    