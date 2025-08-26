A lightweight library for interpolating on a regular / rectilinear multidimensional grid.

This library revolves around the struct `GriddedData`, which stores data values on a regular /
rectilinear grid of arbitrary dimensions. This data can then be used for multivariate interpolation.
Currently, the following algorithms are available:
- Nearest-neighbor interpolation
- n-linear interpolation

If more algorithms are needed, please do not hesistate to open an issue on the repository
website: [https://github.com/StefanMathis/gridded_data]

```rust
use gridded_data::GriddedData;

// 1-dimensional example
let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
let data = vec![0.0, 2.0, 4.0, 2.0, 0.0];
let grid1 = GriddedData::<1>::new([x], data).unwrap();

assert_eq!(grid1.nearest_neighbor(&[0.2]), 0.0);
assert_eq!(grid1.nearest_neighbor(&[1.5]), 2.0);
assert_eq!(grid1.linear_interpolation(&[0.2]), 0.4);
assert_eq!(grid1.linear_interpolation(&[1.5]), 3.0);

// 2-dimensional example (3 x 3 matrix)
// [0 3 6]
// [1 4 7]
// [2 5 8]
let x = vec![0.0, 1.0];
let y = vec![0.0, 1.0];
let data = vec![0.0, 1.0, 2.0, 3.0];
let grid2 = GriddedData::new([x, y], data).unwrap();

assert_eq!(grid2.nearest_neighbor(&[0.75, 0.4]), 2.0);
assert_eq!(grid2.linear_interpolation(&[0.75, 0.4]), 1.9);

// 3-dimensional example (1 x 2 x 4 matrix)
// [0 4] [1 5]
// [2 6] [3 7]
let x = vec![0.0, 1.0];
let y = vec![2.0, 3.0];
let z = vec![4.0, 5.0];
let data = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
let grid3 = GriddedData::new([x, y, z], data).unwrap();

assert_eq!(grid3.nearest_neighbor(&[0.6, 2.3, 4.5]), 4.0);
assert_eq!(grid3.linear_interpolation(&[0.6, 2.3, 4.5]),3.5);

// Points outside the grid are extrapolated via the nearest neighbor method:

assert_eq!(grid3.nearest_neighbor(&[2.0, 3.6, 4.8]), 7.0);
assert_eq!(grid3.linear_interpolation(&[2.0, 3.6, 4.8]), 7.0);
```

# Serialization and deserialization

`GriddedData` can be serialized and deserialized via the (serde)[https://crates.io/crates/serde] crate.
Unfortunately, it is currently not possible to implement this feature for arbitrary dimensions
due to limitations within Rust itself. Therefore, manual implementations for the dimensions
1 to 16 exist.
This functionality is gated behind the `serde` feature flag.

# Providing data via matrix libraries

## nalgebra integration

The function `from_nalgebra_matrix` allows to provide the data values
for a 2-dimensional `GriddedData<2>` via a (nalgebra)[https://crates.io/crates/nalgebra] matrix.
See the function docstring for an example.s
This function is gated behind the `nalgebra` feature flag.