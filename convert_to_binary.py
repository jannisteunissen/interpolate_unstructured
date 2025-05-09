#!/usr/bin/env python3

import os
import argparse
from collections import defaultdict
import numpy as np
import meshio


def write_array_to_binary(array, fh):
    """Write array to binary file"""
    # Use fixed size string of length 64
    dtype_ascii = f'{str(array.dtype):<64}'.encode('ascii')
    shape_int32 = np.array(array.shape, dtype=np.int32)
    ndim_int32 = np.asarray(array.ndim, dtype=np.int32)
    fh.write(dtype_ascii)
    fh.write(ndim_int32)
    fh.write(shape_int32)
    fh.write(array)


def get_cell_neighbors_2d(cells):
    """Determine neighbors for each cell, with cell_neighbor[i_cell, i_vertex]
    being the index of the neighbor connected to the edge between i_vertex and
    i_vertex+1, and -1 if there is no neighbor
    """
    edge_to_cells = defaultdict(list)

    # Create array of neighbors per cell, -1 indicates no neighbor
    cell_neighb = np.full_like(cells, -1, dtype=np.int32)

    # Store cells per edge
    for cell_id, cell in enumerate(cells):
        # Create edges for the current cell
        n_vertices = len(cell)
        for i in range(n_vertices):
            edge = tuple(sorted((cell[i], cell[(i + 1) % n_vertices])))
            edge_to_cells[edge].append(cell_id)

    for cell_id, cell in enumerate(cells):
        n_vertices = len(cell)
        for i in range(n_vertices):
            edge = tuple(sorted((cell[i], cell[(i + 1) % n_vertices])))
            all_cells = edge_to_cells[edge]

            if len(all_cells) == 2:
                # There can only be neighbor if the edge connect two cells
                neighbor = all_cells[1] if all_cells[0] == cell_id \
                    else all_cells[0]
                cell_neighb[cell_id, i] = neighbor

    return cell_neighb


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Convert unstructured grid to binary files')
parser.add_argument('infile', type=str, help='Input file')
parser.add_argument('-output_basename', type=str, help='Basename for output')
args = parser.parse_args()

mesh = meshio.read(args.infile)

if args.output_basename is None:
    args.output_basename = os.path.splitext(args.infile)[0]

with open(args.output_basename + '_points.bin', 'wb') as f:
    write_array_to_binary(mesh.points, f)

with open(args.output_basename + '_cell_type.txt', 'w') as f:
    # Assume single cellblock
    f.write(mesh.cells[0].type + '\n')

with open(args.output_basename + '_cells.bin', 'wb') as f:
    # Assume single cellblock
    write_array_to_binary(mesh.cells[0].data, f)

cell_neighbors = get_cell_neighbors_2d(mesh.cells[0].data)
with open(args.output_basename + '_neighbors.bin', 'wb') as f:
    # Assume single cellblock
    write_array_to_binary(cell_neighbors, f)

for var in mesh.point_data:
    with open(args.output_basename + f'_point_data_{var}.bin', 'wb') as f:
        write_array_to_binary(mesh.point_data[var], f)
