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


def get_cell_neighbors(cells, points_per_face):
    """Determine neighbors for each cell, with cell_neighbor[i_cell, i_vertex]
    being the index of the neighbor connected to the face with index i_vertex,
    and -1 if there is no neighbor
    """
    face_to_cells = defaultdict(list)

    # Create array of neighbors per cell, -1 indicates no neighbor
    cell_neighb = np.full_like(cells, -1, dtype=np.int32)

    # Store cells per face
    for cell_id, cell in enumerate(cells):
        # Create faces for the current cell
        n_vertices = len(cell)
        for i in range(n_vertices):
            face_points = [cell[k % n_vertices] for
                           k in range(points_per_face)]
            face = tuple(sorted(face_points))
            face_to_cells[face].append(cell_id)

    for cell_id, cell in enumerate(cells):
        n_vertices = len(cell)
        for i in range(n_vertices):
            face_points = [cell[k % n_vertices] for
                           k in range(points_per_face)]
            face = tuple(sorted(face_points))
            all_cells = face_to_cells[face]

            if len(all_cells) == 2:
                # There can only be a neighbor if the face connect two cells
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
    if len(mesh.cells) > 1:
        raise ValueError('Mixed cell types not yet implemented')

    f.write(mesh.cells[0].type + '\n')

with open(args.output_basename + '_cells.bin', 'wb') as f:
    write_array_to_binary(mesh.cells[0].data, f)

if mesh.cells[0].type in ['triangle', 'quad']:
    n_points_per_face = 2
elif mesh.cells[0].type in ['tetra']:
    n_points_per_face = 3
else:
    raise ValueError(f'Cell type {mesh.cells[0].type} not implemented')

cell_neighbors = get_cell_neighbors(mesh.cells[0].data, n_points_per_face)

with open(args.output_basename + '_neighbors.bin', 'wb') as f:
    write_array_to_binary(cell_neighbors, f)

for var in mesh.point_data:
    with open(args.output_basename + f'_point_data_{var}.bin', 'wb') as f:
        write_array_to_binary(mesh.point_data[var], f)
