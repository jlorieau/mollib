"""Analysis of partially aligned data for ubiquitin
"""
import glob

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from mollib import Molecule
from mollib.pa import read_pa_file, calc_pa_SVD, Process
from mollib.utils.interactions import interaction_type


# Plot parameters
figsize = (7, 1.85)  # size of each subplot in inches
title_fontsize = 9  # size of the font for subplottitles
title_font = 'Arial'
axes_fontsize = 9  # font size of the axes
axes_font = 'Arial'
label_fontsize = 8  # font size of the tick labels
annotation_fontsize = 7

matplotlib.rc('text', usetex=False)
matplotlib.rc('font', family='sans-serif')
matplotlib.rc('axes', linewidth=0.5)

for i in ('xtick', 'ytick'):
    matplotlib.rcParams[i + '.major.size'] = 2.0
    matplotlib.rcParams[i + '.major.width'] = 0.75
    matplotlib.rcParams[i + '.minor.size'] = 1.0
    matplotlib.rcParams[i + '.minor.width'] = 0.5


# Load the molecule
mol = Molecule('2MJB')


# Load the partial alignment RDC and RACS data
data = {}
for filename in glob.glob('./*.pa'):
    data_subset = read_pa_file(filename)
    data.update(data_subset)


# Process the magnetic interactions
process = Process(mol)
magnetic_interactions = process.process(labels=data.keys())


# Conduct the SVD fit
data_pred, Saupe_components, stats = calc_pa_SVD(magnetic_interactions, data)

# Print the overall statistics
print('Q (%)', stats['Overall']['Q (%)'])

# Prepare a subplot for each interaction type
int_types = ('N-H', 'C', 'N', 'H')  # The different interaction types to plot


f, axarr = plt.subplots(nrows=1, ncols=len(int_types),
                        figsize=figsize)


for col, int_type in enumerate(int_types):
    # Prepare the data for each interaction type
    x_obs = []
    y_pred = []
    for k, v in data.items():
        # Only proceed if there's a predicted data value and it matches the
        # interaction type
        if k not in data_pred or not interaction_type(k) == int_type:
            continue

        x_obs.append(v.value)
        y_pred.append(data_pred[k].value)

    x_obs = np.array(x_obs)
    y_pred = np.array(y_pred)

    # Get the fit stats
    Q = stats[int_type]['Q (%)']
    RMS = stats[int_type]['RMS']
    N = len(x_obs)

    # Set the ticks outside the grid, and only put ticks on the bottom and left
    axarr[col].get_xaxis().set_ticks_position('bottom')
    axarr[col].get_xaxis().set_tick_params(direction='out')
    axarr[col].get_yaxis().set_tick_params(direction='out')
    axarr[col].get_yaxis().set_ticks_position('left')

    # Set the title
    if '-' in int_type:
        axarr[col].set_title(int_type + ' RDC (Hz)',
                             size=title_fontsize,
                             fontname=title_font)
    else:
        axarr[col].set_title(int_type + ' RACS (ppb)',
                             size=title_fontsize,
                             fontname=title_font)

    # Set the x- and y-axis label
    if col == 0:
        axarr[0].set_xlabel('Observed',
                            fontsize=axes_fontsize, fontname=axes_font)
        axarr[0].set_ylabel('Predicted',
                            fontsize=axes_fontsize, fontname=axes_font)

    # Set the axis tick label size
    axarr[col].tick_params(labelsize=label_fontsize)

    # Set the number of ticks for the x- and y-axes. The x- and y-axes should
    # have the same ticks and ranges
    locator = MaxNLocator(steps=[1, 2, 8, 10])
    axarr[col].set_aspect('equal', 'datalim')
    axarr[col].get_xaxis().set_major_locator(locator)
    axarr[col].get_yaxis().set_major_locator(locator)

    # Set the axis tick label font
    labels = axarr[col].get_xticklabels() + axarr[col].get_yticklabels()
    for label in labels:
        label.set_fontname(axes_font)

    # Add the annotations
    if '-' in int_type:
        text = "Q={}%\nRMS={}Hz\nN={}".format(Q, RMS, N)  # RDC
    else:
        text = "Q={}%\nRMS={}ppb\nN={}".format(Q, RMS, N)  # RACS
    axarr[col].text(0.05, 0.95, text,
                    verticalalignment='top',
                    horizontalalignment='left',
                    transform=axarr[col].transAxes,
                    fontsize=annotation_fontsize)

    # Plot the subplot
    axarr[col].scatter(x_obs, y_pred)


# Save the figure
plt.tight_layout(pad=0.0)
plt.savefig('ubq_2mjb_racs_lowres.png', format='PNG',
            dpi=220, bbox_inches='tight', pad_inches=0.02)
plt.savefig('ubq_2mjb_racs.png', format='PNG',
            dpi=1200, bbox_inches='tight', pad_inches=0.02)
plt.savefig('ubq_2mjb_racs.svg', format='SVG',
            bbox_inches='tight', pad_inches=0.02)


