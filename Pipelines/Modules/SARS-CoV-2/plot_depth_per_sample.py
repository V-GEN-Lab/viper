# !/usr/bin/env python3
# -*- coding: utf-8 -*-
 
import sys

from matplotlib import (pyplot as plt,
                        lines)
 
 
def parse_depth(depth_input, genome_size):
    """Parse depth file.
 
    Args:
        depth_input (str): Path to depth file.
        genome_size (int): Genome size.
 
    Returns:
        list: List with depth.
 
    """
    depth = [0] * genome_size
    references = set()
 
    with open(depth_input) as depth_object:
        for row in depth_object:
            genome_id, position, depth_count = row.split()
            position_index = int(position) - 1

            references.add(genome_id)
 
            if len(references) > 1:
                raise Exception(' This script only handles one genome - contig.')

            if 0 <= position_index < genome_size:
                depth[position_index] = int(depth_count)
 
    return depth
 
 
def plot_depth(depth_report, output_name, plot_title, genome_size, normalize=False, depth_cut_off=20):
    """Plot genome Depth across genome.
 
    Args:
        depth_report (str): Path to samtool's depth file.
        output_name (str): Path to output PNG image.
        plot_title (str): Plot title.
        genome_size (int): Genome size.
        normalize (bool): If `True`, normalizes the depth by the largest depth (default = `False`).
        depth_cut_off (int): Plot a line to represent a targeted depth (default = 20).
 
    """
    data = parse_depth(depth_report, genome_size)
 
    y_label = "Normalized Depth" if normalize else "Depth"
    max_depth = max(data) if data else 0
    if normalize and max_depth > 0:
        data = [xx / max_depth for xx in data]
 
    figure, ax = plt.subplots()
    ax.set_title(plot_title)
 
    ax.plot(range(len(data)), data, linewidth=1.5)
    ax.set_xlabel('Genome Position (bp)')
    ax.set_ylabel(y_label)
 
    if not normalize:
        ax.add_line(lines.Line2D([0, genome_size + 1], [depth_cut_off], color="r"))
 
    figure.set_size_inches(14, 10)
    plt.savefig(output_name, bbox_inches='tight', dpi=400)
    plt.close()
 
    print("Done :)")


# Getting files from sys.argv
depth_file = sys.argv[1]
depth_plot = sys.argv[2]
sample_name = sys.argv[3]

plot_depth(depth_file, depth_plot, sample_name, 29903, False, 100)
