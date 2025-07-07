#!/usr/bin/env python3

import struct
import os
import argparse
from collections import defaultdict
import numpy as np
import meshio


class BindaWriter:
    """
    A class to store multiple binary data entries along with their metadata.

    Each entry consists of a name, data type, number of dimensions,
    shape, and an offset to the binary data within the file.
    """

    def __init__(self):
        """Initialize the BindaWriter."""
        self.entries = []
        self.binary_data_storage = bytearray()

    def add_entry(self, name, data, metadata=""):
        """
        Add a data entry to the store.

        Args:
            name (str): An identifier for the data (max length 128).
            data (np.ndarray): The data to store (must be a numpy array).
            metadata (str): optional metadata (max length 128).

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

        # Append entry info (name, data type, metadata, ndim, shape, offset)
        self.entries.append((
            name.ljust(128).encode('ascii'),
            str(data.dtype).ljust(128).encode('ascii'),
            metadata.ljust(128).encode('ascii'),
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

        # Size of header per entry
        header_size = (
            struct.calcsize('128s') +   # Size for name
            struct.calcsize('128s') +   # Size for data type
            struct.calcsize('128s') +   # Size for metadata
            struct.calcsize('q') +      # Size for ndim
            struct.calcsize('8q') +     # Size for shape, 8 dimensions max
            struct.calcsize('q')        # Size for offset
        )

        n_entries = len(self.entries)
        total_header_size = (
            struct.calcsize('8s') +  # size for 'BINDA' identifier
            struct.calcsize('q') +  # size for n_entries
            struct.calcsize('q') +  # size for total_header_size
            n_entries * header_size  # size of rest of header
        )

        with open(filename, 'wb') as f:
            identifier = 'BINDA'.ljust(8).encode('ascii')
            f.write(struct.pack('8s', identifier))

            # Write the number of entries and total header size
            f.write(struct.pack('q', n_entries))
            f.write(struct.pack('q', total_header_size))

            # Write each entry's metadata
            for entry in self.entries:
                name, data_type, metadata, ndim, shape, offset = entry
                offset += total_header_size  # Adjust for header size
                f.write(struct.pack('128s', name))
                f.write(struct.pack('128s', data_type))
                f.write(struct.pack('128s', metadata))
                f.write(struct.pack('q', ndim))
                # Pad `shape` with zeros if its length is smaller than 8
                f.write(struct.pack('8q', *(shape + (0,) * (8 - len(shape)))))
                f.write(struct.pack('q', offset))

            # Write all binary data at once at the end of the file
            f.write(self.binary_data_storage)


def get_cell_neighbors(cells, points, n_points_face):
    """Determine neighbors for each cell, with cell_neighbor[i_cell, i_vertex]
    being the index of the neighbor connected to the face with index i_vertex,
    and -1 if there is no neighbor. This method first 'removes' duplicate
    points to more robustly find neighbors.

    """
    face_to_cells = defaultdict(list)

    # Create array of neighbors per cell, -1 indicates no neighbor
    cell_neighb = np.full_like(cells, -1, dtype=np.int32)

    points_uniq, idx = np.unique(points, axis=0, return_inverse=True)

    if len(points_uniq) < len(points):
        print(f'Found {len(points) - len(points_uniq)} duplicate points')

    # Cells but now with index to unique points
    cells_uniq = idx[cells.reshape(-1)].reshape(cells.shape)

    # Store cells per face
    for cell_id, cell in enumerate(cells_uniq):
        # Create faces for the current cell
        n_vertices = len(cell)
        for i in range(n_vertices):
            face_points = [cell[(i + k) % n_vertices]
                           for k in range(n_points_face)]
            face = tuple(sorted(face_points))
            face_to_cells[face].append(cell_id)

    for cell_id, cell in enumerate(cells_uniq):
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

if args.output_basename is None:
    args.output_basename = os.path.splitext(args.infile)[0]

mesh = meshio.read(args.infile)

if len(mesh.cells) > 1:
    raise ValueError('Mixed cell types not yet implemented')

if mesh.cells[0].type in ['triangle', 'quad']:
    n_points_per_face = 2
elif mesh.cells[0].type in ['tetra']:
    n_points_per_face = 3
else:
    raise ValueError(f'Cell type {mesh.cells[0].type} not implemented')

cell_neighbors = get_cell_neighbors(mesh.cells[0].data, mesh.points,
                                    n_points_per_face)

binstore = BindaWriter()

binstore.add_entry('points', mesh.points)
binstore.add_entry('cells', mesh.cells[0].data, mesh.cells[0].type)
binstore.add_entry('cell_neighbors', cell_neighbors)

for var in mesh.point_data:
    clean_name = var.replace(',', '')
    binstore.add_entry('point_data', mesh.point_data[var], clean_name)
    print('Storing point data:', clean_name)

for var in mesh.cell_data:
    clean_name = var.replace(',', '')
    binstore.add_entry('cell_data', mesh.cell_data[var], clean_name)
    print('Storing cell data: ', clean_name)

fname = args.output_basename + '.binda'
binstore.write_to_file(fname)
print(f'Stored {fname}')
