# Efficient interpolation on unstructured grids

This library can be used to efficiently interpolate on unstructured grids. Two approaches are implemented to interpolate at a location **r**:

1. If a previous (nearby) position and cell are known, traverse the grid to find the cell that contains **r**
2. Otherwise, use a kd-tree to locate a nearby cell, followed by step 1.

# Supported cell types

* Triangles (2d)

# Compilation

TODO
