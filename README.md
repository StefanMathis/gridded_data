gridded_data
============

A lightweight library for interpolating on a regular / rectilinear multidimensional grid.

[`GriddedData`]: https://docs.rs/gridded_data/0.1.1/gridded_data/struct.GriddedData.html

This library revolves around the struct [`GriddedData`], which stores data values on a regular /
rectilinear grid of arbitrary dimensions. This data can then be used for [multivariate interpolation](https://en.wikipedia.org/wiki/Multivariate_interpolation).
Currently, the following algorithms are available:
- Nearest-neighbor interpolation
- n-linear interpolation

If more algorithms are needed, please do not hesistate to open an issue on the repository
website: [https://github.com/StefanMathis/gridded_data](https://github.com/StefanMathis/gridded_data)

The full API documentation is available at [https://docs.rs/gridded_data/0.1.1/gridded_data/](https://docs.rs/gridded_data/0.1.1/gridded_data/).

# Concept

The [`GriddedData`] struct defines a n-dimensional grid via a corresponding number of axes. Each axis is defined as a vector of "nodes". For example, this axis
has the three nodes 1, 2 and 5: `ax = [1, 2, 5]`. An axis vector must be strictly monotonically increasing (each node must be smaller than its successor) and it must have at least two nodes.

One or more axis define a grid. As an example, the three axes `x = [1, 2, 5]`; `y = [0, 1]`; `z = [-1, 3]` result in a 3-dimensional cube with 12 [vertices](https://en.wikipedia.org/wiki/Vertex_(geometry)). The grid vertices are created via the Cartesian product of the axes:
`[1, 0, -1]`; `[1, 0, 3]`; `[1, 1, -1]`; `[1, 1, 3]` and so on.

The individual grid vertices form n-dimensional hypercuboids which are called "cells". These can be identified via the n index pairs of the axes.
The index pairs `[1, 2], [0, 1], [0, 1]` identify the hypercuboid formed by the vertices `[2, 0, -1]`; `[2, 0, 3]`; `[2, 1, -1]`; `[2, 1, 3]`; `[5, 0, -1]`; `[5, 0, 3]`; `[5, 1, -1]`; `[5, 1, 3]`. The method [`GriddedData::cell_bounds`](https://docs.rs/gridded_data/0.1.1/gridded_data/struct.GriddedData.html#method.cell_bounds) can be used to find the index pairs of cells.

For each vertex, a corresponding value needs to be given during the construction of [`GriddedData`]. This is done via a `data` vector which provides the data in [row-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order). For the three-dimensional grid defined above, `data` therefore needs 12 entries, e.g.: `data = [1, 2, ..., 11, 12]`. Examples can be found in the docstring of [`GriddedData::new`](https://docs.rs/gridded_data/0.1.1/gridded_data/struct.GriddedData.html#method.new)

This results in the following vertex–value pairs:
| Vertex      | Data value |
| ----------- | ---------- |
| `[1, 0, -1]`| `0`        |
| `[1, 0, 3]` | `1`        |
| `[1, 1, -1]`| `2`        |
| ...         | ...        |
| `[5, 1, 3]` | `12`       |

For interpolation purposes, the data values are treated as the output of an underlying (unknown) function of the corresponding vertices.

# Examples

This section contains some examples on how to use [`GriddedData`] for inter- and extrapolation:

```rust
use gridded_data::GriddedData;

/*
1-dimensional example:
The underlying function is known to have the following vertex–value pairs:
vertex : 0 1 2 3 4
         ---------
data:    0 2 4 2 0
*/
let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
let data = vec![0.0, 2.0, 4.0, 2.0, 0.0];
let grid1 = GriddedData::new([x], data).unwrap();

// Access the data at the grid points
assert_eq!(*grid1.get_cart(&[3]).unwrap(), 2.0);

// Get values between grid points via different interpolation methods
assert_eq!(grid1.nearest_neighbor_interp(&[0.2]), 0.0);
assert_eq!(grid1.nearest_neighbor_interp(&[1.5]), 2.0);
assert_eq!(grid1.linear_interp(&[0.2]), 0.4);
assert_eq!(grid1.linear_interp(&[1.5]), 3.0);

/*
2-dimensional example:
The underlying function is known to have the following vertex–value pairs:
  x 0 1 2
y -------
0 | 0 1 2
1 | 3 4 5
*/
let x = vec![0.0, 1.0];
let y = vec![0.0, 1.0, 2.0];
let data = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
let grid2 = GriddedData::new([x, y], data).unwrap();

// Access the data at the grid points
assert_eq!(*grid2.get_cart(&[0, 2]).unwrap(), 2.0);

// Get values between grid points via different interpolation methods
assert_eq!(grid2.nearest_neighbor_interp(&[0.75, 0.4]), 3.0);
assert_eq!(grid2.linear_interp(&[0.5, 0.25]), 1.75);
```

# Feature flags

All features are disabled by default.

## Serialization and deserialization

The [`GriddedData`] struct can be serialized and deserialized via the [serde](https://crates.io/crates/serde) crate.
Unfortunately, it is currently not possible to implement this feature for arbitrary dimensions
due to limitations within Rust itself. Therefore, manual implementations for the dimensions
1 to 16 exist.

This functionality is gated behind the **serde** feature flag.

## Providing data via matrix libraries

### nalgebra integration

The function [`from_nalgebra_matrix`]() allows to provide the data values
for a 2-dimensional `GriddedData<2>` via a [nalgebra](https://crates.io/crates/nalgebra) matrix.
See the function docstring for an example.

This function is gated behind the **nalgebra** feature flag.