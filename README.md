# Efficient interpolation on unstructured grids

This library can be used to efficiently interpolate on unstructured grids. Two approaches are implemented to interpolate at a location **r**:

1. If a previous (nearby) position and cell are known, traverse the grid to find the cell that contains **r**
2. Otherwise, use a kd-tree to locate a nearby cell, followed by step 1.

# Supported cell types

* Triangles (2D)
* Quads (2D)
* Tetrahedra (3D)

# Installation

    git clone --recurse-submodules <url>
    cd interpolate_unstructured
    make

The `convert_to_binary.py` script depends on the [meshio](https://github.com/nschloe/meshio/) library.

# Usage

See the included examples. The basic step are to:

1. Convert an unstructured grid to binary format with the `convert_to_binary.py` script. Point and cell data will be included.
2. Call `iu_read_grid` to read the binary file.
3. Call `iu_interpolate_at` or `iu_interpolate_scalar_at` to interpolate point data at a position.

The following methods are available:

* `iu_get_point_data_index`: Find index of point data with a given name
* `iu_get_cell_data_index`: Find index of cell data with a given name
* `iu_get_cell_center`: Compute center of a cell
* `iu_get_cell`: Determine the cell that contains a point **r**
* `iu_get_cell_scalar_at`: Get value of cell data at a point **r**
* `iu_interpolate_at`: Interpolate multiple point-data variables at **r** simultaneously
* `iu_interpolate_scalar_at`: Interpolate a single point-data variable at **r**

