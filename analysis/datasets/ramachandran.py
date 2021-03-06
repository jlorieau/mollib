"""Analysis of the Ramachandran datasets.
"""

import tarfile
import json
from itertools import izip_longest as zip_longest

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import ndimage

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

# import the dataset
filename = '../../mollib/data/ramachandranstatistics/measurements.tar'
tfile = tarfile.open(name=filename, mode='r')
measurement_dict = {}

with tfile:
    # Extract the files
    for member in tfile.getmembers():
        f= tfile.extractfile(member)
        try:
            identifier = member.name.strip('.json')
            string = f.read().decode()
            return_dict = json.loads(string)
            measurement_dict[identifier] = return_dict
        except KeyError:
            continue
        finally:
            f.close()


# Import the phi/psi angles
results = {}
for identifier, return_dict in measurement_dict.items():
    for classification, phi_psi_list in return_dict.items():
        l = results.setdefault(classification, list())
        phi_psi_list = [(i,j) for i,j in phi_psi_list if isinstance(i, float)
                        and isinstance(j, float)]
        l.extend(phi_psi_list)


# Create and Overall dataset
phi_list = []
psi_list = []
for classification, phi_psi in results.items():
    phi, psi = zip(*phi_psi)
    phi_list += phi
    psi_list += psi

results['Overall'] = zip(phi_list, psi_list)


# Prepare the plots

# Define the plots. These are organized by keys in results and the
# corresponding color map.
# The first item is the dataset name, the second item is the color map, and the
# third item is the position of the text label for the number of items used to
# calculate the density map
labels = (('Overall', plt.cm.plasma, (0.95, 0.01)),
          ('alpha-helix', plt.cm.Greens_r, (0.95, 0.01)),
          ('alpha-helix__N-term', plt.cm.Greens_r, (0.95, 0.01)),
          ('alpha-helix__C-term', plt.cm.Greens_r, (0.95, 0.01)),
          ('310-helix', plt.cm.Greens_r, (0.95, 0.01)),
          ('pi-helix', plt.cm.Greens_r, (0.95, 0.01)),
          ('sheet', plt.cm.Blues_r, (0.95, 0.01)),
          ('sheet__N-term', plt.cm.Blues_r, (0.95, 0.01)),
          ('sheet__C-term', plt.cm.Blues_r, (0.95, 0.01)),
          ('type I turn', plt.cm.Reds_r, (0.95, 0.01)),
          ('type II turn', plt.cm.Reds_r, (0.95, 0.01)),
          ("type I' turn", plt.cm.Reds_r, (0.95, 0.01)),
          ("type II' turn", plt.cm.Reds_r, (0.95, 0.89)),
          ('Gly', plt.cm.Greys_r, (0.95, 0.70)),
          ('No classification', plt.cm.Greys_r, (0.95, 0.01)),
          )


def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


subfigure_groups = list(grouper(4, labels))


# Make rows of plots
for row, subfigure_labels in enumerate(subfigure_groups, 1):
    f, axarr = plt.subplots(nrows=1, ncols=len(subfigure_labels), sharey=True,
                            figsize=figsize)

    # Set the shared y-axis label
    axarr[0].set_ylabel(r'$\psi$ (deg)', fontsize=axes_fontsize,
                        fontname=axes_font)

    for count, values in enumerate(subfigure_labels):

        # Skip empty plots
        if values is None:
            axarr[count].axis('off')
            # axarr[count].set_visible(False)
            continue
        label, cmap, annot_xy = values

        # Prepare the data and create 2d histograms
        phi_psi = results[label]
        x, y = zip(*phi_psi)

        x = np.array(x)
        y = np.array(y)
        N = len(x)
        hist2d, x, y = np.histogram2d(y, x, bins=36, range=np.array(
            [(-180., 190.), (-180., 190.)]))
        hist2d = -1. * np.log(hist2d + 0.1)

        minimum = np.min(hist2d)
        hist2d -= minimum

        levels = np.arange(0.0, 5.1, 1.0)

        # Set the title and x-axis label
        title = (label.replace('__', ' ')
                      .replace('alpha', "$\\alpha$")
                      .replace('pi', "$\\pi$")
                      .replace("\\'", "'")
                      .replace('No', 'no'))

        axarr[count].set_title(title,
                               size=title_fontsize,
                               fontname=title_font)

        # Set the x-axis label
        axarr[count].set_xlabel(r'$\phi$ (deg)', fontsize=axes_fontsize,
                                fontname=axes_font)

        # Set the axis tick spacing
        axarr[count].xaxis.set_ticks(np.arange(-180, 181, 90))
        axarr[count].set_xlim(-180, 180)
        axarr[count].yaxis.set_ticks(np.arange(-180, 181, 90))
        axarr[count].set_ylim(-180, 180)

        # Set the axis tick label size
        axarr[count].tick_params(labelsize=label_fontsize)

        # Set the axis tick label font
        labels = axarr[count].get_xticklabels() + axarr[count].get_yticklabels()
        for label in labels:
            label.set_fontname(axes_font)

        # Annotate the number of measurements on the plot
        if annot_xy is not None:
            axarr[count].text(annot_xy[0], annot_xy[1],
                              "N={:,.0f}".format(N),
                              verticalalignment='bottom',
                              horizontalalignment='right',
                              transform=axarr[count].transAxes,
                              fontsize=annotation_fontsize)

        # Create the 2d contour plot
        axarr[count].contourf(x[:-1], y[:-1], hist2d, levels, cmap=cmap)

    # Save the figures
    plt.savefig('ramachandran/ramachandran_countour_{}.png'.format(row),
                format='PNG', dpi=1200, bbox_inches='tight', pad_inches=0.05)
    plt.savefig('ramachandran/ramachandran_countour_{}_lowres.png'.format(row),
                format='PNG', dpi=220, bbox_inches='tight', pad_inches=0.02)
    plt.savefig('ramachandran/ramachandran_countour_{}.svg'.format(row),
                format='SVG', bbox_inches='tight', pad_inches=0.05)

# Prepare an overall contour plot with lines
phi_psi = results['Overall']
x, y = zip(*phi_psi)

x = np.array(x)
y = np.array(y)
N = len(x)

# Convert to a histogram
hist2d, x, y = np.histogram2d(y, x, bins=36, range=np.array(
        [(-180., 185.), (-180., 185.)]))
hist2d = -1. * np.log(hist2d + 0.1)

minimum = np.min(hist2d)
hist2d -= minimum

# Optionally smooth the data
# hist2d = ndimage.gaussian_filter(hist2d, sigma=0.25, order=0)

# Contour levels at 98% (4.0) and 99.8%
levels = np.array([4.1, 6.1])

f, axarr = plt.subplots(nrows=1, ncols=1, sharey=True,
                            figsize=(4,4))
# Set the x-axis and y-axis labels
axarr.set_xlabel(r'$\phi$ (deg)', fontsize=axes_fontsize,
                 fontname=axes_font)
axarr.set_ylabel(r'$\psi$ (deg)', fontsize=axes_fontsize,
                 fontname=axes_font)

# Set the axis tick spacing
axarr.xaxis.set_ticks(np.arange(-180, 181, 90))
axarr.set_xlim(-180, 180)
axarr.yaxis.set_ticks(np.arange(-180, 181, 90))
axarr.set_ylim(-180, 180)

# Set the axis tick label size
axarr.tick_params(labelsize=label_fontsize)

# Outer ticks
axarr.get_yaxis().set_tick_params(direction='out')
axarr.get_xaxis().set_tick_params(direction='out')

# Create the 2d contour plot
axarr.contour(x[:-1], y[:-1], hist2d, levels, cmap=plt.cm.Blues_r)

plt.savefig('ramachandran/ramachandran_line_countour_overall_lowres.png',
            format='PNG', dpi=220, bbox_inches='tight', pad_inches=0.02)
plt.savefig('ramachandran/ramachandran_line_countour_overall.svg',
            format='SVG', bbox_inches='tight', pad_inches=0.05)
