"""Analysis of the Hydrogen bond datasets.
"""
import tarfile
import json
from itertools import izip_longest as zip_longest

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# Plot parameters
figsize = (7, 1.5)  # size of each subplot in inches
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

# Histogram parameters
d_hist_range = (1.5, 3.0)
# d_hist_range = (-180, 180)
theta_hist_range = (60, 180.)
bins = 16

# import the dataset
filename = '../../mollib/data/hbondstatistics/measurements.tar'
tfile = tarfile.open(name=filename, mode='r')
measurement_dict = {}

with tfile:
    # Extract the files
    for member in tfile.getmembers():
        f = tfile.extractfile(member)
        try:
            identifier = member.name.strip('.json')
            string = f.read().decode()
            return_dict = json.loads(string)
            measurement_dict[identifier] = return_dict
        except KeyError:
            continue
        finally:
            f.close()

d1a1_theta = {}
for identifier, return_dict in measurement_dict.items():
    for classification, hbond_dict_list in return_dict.items():
        major, minor, modifier = classification.split('__')
        if major != 'bb-bb amide':
            continue
        values = d1a1_theta.setdefault(minor, list())
        for d in hbond_dict_list:
            # values += [(d['angles']['theta'], d['angles']['phi'])]
            values += [(d['distances']['d1a1'], d['angles']['theta'])]


# Define the plots. These are organized by keys in results and the
# corresponding color map.
labels = (('alpha-helix', plt.cm.Greens_r),
          ('310-helix', plt.cm.Greens_r),
          ('pi-helix', plt.cm.Greens_r),
          ('sheet', plt.cm.Blues_r),
          ('type I turn', plt.cm.Reds_r),
          ('type II turn', plt.cm.Reds_r),
          ("type I' turn", plt.cm.Reds_r),
          ("type II' turn", plt.cm.Reds_r),
          ('isolated', plt.cm.Greys_r),
          )


def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


subfigure_groups = list(grouper(4, labels))

for row, subfigure_labels in enumerate(subfigure_groups, 1):

    # Create the subplots in this row
    f, axarr = plt.subplots(1, len(subfigure_labels), sharey=True,
                            figsize=figsize)

    # Set the shared y-axis label
    axarr[0].set_ylabel(r'$\theta$ (deg)', fontsize=axes_fontsize,
                        fontname=axes_font)

    for count, values in enumerate(subfigure_labels):
        if values is None:
            axarr[count].axis('off')
            continue
        title, cmap = values

        # Prepare the data
        xy = d1a1_theta[title]

        x = np.array([i[0] for i in xy])
        y = np.array([i[1] for i in xy])

        hist2d, y, x = np.histogram2d(y, x, bins=36, range=np.array(
            [theta_hist_range, d_hist_range]))
        hist2d = -1. * np.log(hist2d + 0.1)

        minimum = np.min(hist2d)
        hist2d -= minimum

        # Set the plot title
        axarr[count].set_title(title.replace('__', ' ')
                               .replace('alpha', "$\\alpha$")
                               .replace('pi', "$\\pi$")
                               .replace('No', 'no'),
                               size=title_fontsize,
                               fontname=title_font)

        # Set the x-axis label
        axarr[count].set_xlabel(r'$d_{H...O}$ (A)', fontsize=axes_fontsize,
                                fontname=axes_font)

        # Set the x-axis ticks
        # axarr[count].xaxis.set_ticks(np.arange(1.5, 3.0, 0.5))
        # axarr[count].set_xlim(1.5, 3.0)

        # Set the axis tick label size
        axarr[count].tick_params(labelsize=label_fontsize)

        # Set the axis tick label font
        labels = axarr[count].get_xticklabels() + axarr[count].get_yticklabels()
        for label in labels:
            label.set_fontname(axes_font)

        # Place ticks on the bottom and left side only
        axarr[count].get_xaxis().tick_bottom()
        axarr[count].get_yaxis().tick_left()

        # Set the contour levels
        levels = np.arange(0.0, 5.1, 1.0)

        axarr[count].contourf(x[:-1], y[:-1], hist2d, levels, cmap=cmap)

    plt.savefig('hbond_countour_{}.png'.format(row), format='PNG',
                dpi=1200, bbox_inches='tight', pad_inches=0.05)
    plt.savefig('hbond_countour_{}.svg'.format(row), format='SVG',
                bbox_inches='tight', pad_inches=0.05)
