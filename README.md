A lightweight library for interpolating on a regular / rectilinear multidimensional grid.

This library revolves around the struct `GriddedData`, which stores data values on a regular /
rectilinear grid of arbitrary dimensions. This data can then be used for [multivariate interpolation](https://en.wikipedia.org/wiki/Multivariate_interpolation).
Currently, the following algorithms are available:
- Nearest-neighbor interpolation
- n-linear interpolation

If more algorithms are needed, please do not hesistate to open an issue on the repository
website: [https://github.com/StefanMathis/gridded_data]

```rust
use gridded_data::GriddedData;

/*
1-dimensional example:
The underlying function is known to have the following data / node pairs:
node : 0 1 2 3 4
       ---------
data:  0 2 4 2 0
*/
let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
let data = vec![0.0, 2.0, 4.0, 2.0, 0.0];
let grid1 = GriddedData::<1>::new([x], data).unwrap();

// Access the data at the grid points
assert_eq!(*grid1.get_cart(&[3]).unwrap(), 2.0);

// Get values between grid points via different interpolation methods
assert_eq!(grid1.nearest_neighbor_interp(&[0.2]), 0.0);
assert_eq!(grid1.nearest_neighbor_interp(&[1.5]), 2.0);
assert_eq!(grid1.linear_interp(&[0.2]), 0.4);
assert_eq!(grid1.linear_interp(&[1.5]), 3.0);

/*
2-dimensional example:
The underlying function is known to have the following data / node pairs (with a node defined by both its x- and y-value):
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

# Serialization and deserialization

`GriddedData` can be serialized and deserialized via the [serde](https://crates.io/crates/serde) crate.
Unfortunately, it is currently not possible to implement this feature for arbitrary dimensions
due to limitations within Rust itself. Therefore, manual implementations for the dimensions
1 to 16 exist.

This functionality is gated behind the `serde` feature flag.

# Providing data via matrix libraries

## nalgebra integration

The function [`from_nalgebra_matrix`] allows to provide the data values
for a 2-dimensional `GriddedData<2>` via a [nalgebra](https://crates.io/crates/nalgebra) matrix.
See the function docstring for an example.

This function is gated behind the `nalgebra` feature flag.