#!/usr/bin/env python3

import struct
import os
import argparse
from collections import defaultdict
import numpy as np
import meshio


class BinaryDataStore:
    """
    A class to store multiple binary data entries along with their metadata.

    Each entry consists of a name, data type, number of dimensions,
    shape, and an offset to the binary data within the file.
    """

    def __init__(self):
        """Initialize the BinaryDataStore."""
        self.entries = []
        self.binary_data_storage = bytearray()

    def add_entry(self, name, data):
        """
        Add a data entry to the store.

        Args:
            name (str): An identifier for the data (max length 128).
            data (np.ndarray): The data to store (must be a numpy array).

        Raises:
            ValueError: If the name exceeds 128 characters or
                         if the data is not a numpy array
                         or if the array dimensionality exceeds 8.
        """
        if len(name) > 128:
            raise ValueError("Name must be at most 128 characters.")

        if not isinstance(data, np.ndarray):
            raise ValueError("Data must be a numpy array.")

        if data.ndim > 8:
            raise ValueError("Number of dimensions cannot exceed 8.")

        # Convert numpy array to binary data
        binary_data = data.tobytes()
        self.binary_data_storage.extend(binary_data)

        # Compute offset for the current data
        offset = len(self.binary_data_storage) - len(binary_data)

        # Append entry info (name, data type, ndim, shape, offset)
        self.entries.append((
            name.ljust(128).encode('ascii'),
            str(data.dtype).ljust(128).encode('ascii'),
            data.ndim,
            data.shape,
            offset
        ))

    def write_to_file(self, filename):
        """
        Write the stored entries and binary data to a file.

        Args:
            filename (str): Path to the file where the data will be written.
        """
        header_size = (
            struct.calcsize('128s') +   # Size for name
            struct.calcsize('128s') +   # Size for data type
            struct.calcsize('q') +      # Size for ndim
            struct.calcsize('8q') +     # Size for shape, 8 dimensions max
            struct.calcsize('q')        # Size for offset
        )

        n_entries = len(self.entries)
        total_header_size = 2 * struct.calcsize('q') + n_entries * header_size

        with open(filename, 'wb') as f:
            # Write the number of entries and total header size
            f.write(struct.pack('q', n_entries))
            f.write(struct.pack('q', total_header_size))

            # Write each entry's metadata
            for entry in self.entries:
                name, data_type, ndim, shape, offset = entry
                offset += total_header_size  # Adjust for header size
                f.write(struct.pack('128s', name))
                f.write(struct.pack('128s', data_type))
                f.write(struct.pack('q', ndim))
                # Pad `shape` with zeros if its length is smaller than 8
                f.write(struct.pack('8q', *(shape + (0,) * (8 - len(shape)))))
                f.write(struct.pack('q', offset))

            # Write all binary data at once at the end of the file
            f.write(self.binary_data_storage)


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


def get_cell_neighbors(cells, n_points_face):
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
            face_points = [cell[(i + k) % n_vertices]
                           for k in range(n_points_face)]
            face = tuple(sorted(face_points))
            face_to_cells[face].append(cell_id)

    for cell_id, cell in enumerate(cells):
        n_vertices = len(cell)
        for i in range(n_vertices):
            face_points = [cell[(i + k) % n_vertices]
                           for k in range(n_points_face)]
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
binstore = BinaryDataStore()

if args.output_basename is None:
    args.output_basename = os.path.splitext(args.infile)[0]

with open(args.output_basename + '_points.bin', 'wb') as f:
    write_array_to_binary(mesh.points, f)
    binstore.add_entry('points', mesh.points)

with open(args.output_basename + '_cell_type.txt', 'w') as f:
    if len(mesh.cells) > 1:
        raise ValueError('Mixed cell types not yet implemented')

    f.write(mesh.cells[0].type + '\n')

with open(args.output_basename + '_cells.bin', 'wb') as f:
    write_array_to_binary(mesh.cells[0].data, f)
    binstore.add_entry('cells_' + mesh.cells[0].type, mesh.cells[0].data)

if mesh.cells[0].type in ['triangle', 'quad']:
    n_points_per_face = 2
elif mesh.cells[0].type in ['tetra']:
    n_points_per_face = 3
else:
    raise ValueError(f'Cell type {mesh.cells[0].type} not implemented')

cell_neighbors = get_cell_neighbors(mesh.cells[0].data, n_points_per_face)

with open(args.output_basename + '_neighbors.bin', 'wb') as f:
    write_array_to_binary(cell_neighbors, f)
    binstore.add_entry('cell_neighbors', cell_neighbors)

for var in mesh.point_data:
    with open(args.output_basename + f'_point_data_{var}.bin', 'wb') as f:
        write_array_to_binary(mesh.point_data[var], f)
        binstore.add_entry(f'point_data_{var}', mesh.point_data[var])

binstore.write_to_file(args.output_basename + '.bin')
