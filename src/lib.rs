/*!
A lightweight library for interpolating on a regular / rectilinear multidimensional grid.

This library revolves around the struct [`GriddedData`], which stores data values on a regular /
rectilinear grid of arbitrary dimensions. This data can then be used for [multivariate interpolation](https://en.wikipedia.org/wiki/Multivariate_interpolation).
Currently, the following algorithms are available:
- Nearest-neighbor interpolation
- n-linear interpolation

If more algorithms are needed, please do not hesistate to open an issue on the repository
website: (https://github.com/StefanMathis/gridded_data)

```
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
*/

use cart_lin::{CartesianIndices, cart_to_lin, cart_to_lin_unchecked};

#[cfg(feature = "serde")]
mod serde;

/**
This struct contains the N-dimensional grid data and provides interpolation methods which work on this data.
 */
#[derive(Debug, Clone)]
pub struct GriddedData<const N: usize> {
    pub(crate) axes: [Vec<f64>; N],
    pub(crate) data: Vec<f64>,
}

impl<const N: usize> GriddedData<N> {
    /**
    Creates a new n-dimensional [`GriddedData`] with the number of dimensions defined by the input.

    If the arguments are valid, this function creates a new instance of `GriddedData` with the number
    of dimension equal to the length of `axes`. The given arguments  must fulfil the following criteria:
    - Each vector in the `axes` array must be strictly monotonic increasing
    - Each vector in the `axes` must have at least two values (if it doesn't, the corresponding
    dimension collapses anyway and the vector can be omitted entirely).
    - The number of datapoints (`data.len()`) must be equal to the product of the lengths of the
    axis vectors within `axes`.

    The `data` must be given in row-major fashion (index of last dimension changes fastest). Some examples:

    - The 2D-matrix
        ```text
        1.0 3.0 5.0
        2.0 4.0 6.0
        ```
        is represented by `data = [1.0 3.0 5.0 2.0 4.0 6.0]` (first dimension is the row number, second one the column number).

    - The 3D-matrix (| denotes a page / 3rd dimension separator)
        ```text
        1.0 3.0 5.0 | 7.0  9.0 11.0
        2.0 4.0 6.0 | 8.0 10.0 12.0
        ```
        is represented by `data = [1.0 7.0 3.0 9.0 5.0 11.0 2.0 8.0 4.0 10.0 6.0 12.0]`
        (first dimension is the row number, second one the column number, third one the page number).

    # Examples

    ```
    use gridded_data::GriddedData;

    /*
    2 x 3 x 3 matrix
      z = 1.0     z = 2.0
      y 3 4 6     3  4  6
    x -------     -------
    0 | 0 3 5     7  9 11
    1 | 2 4 6     8 10 12
     */
    let data = vec![1.0, 7.0, 3.0, 9.0, 5.0, 11.0, 2.0, 8.0, 4.0, 10.0, 6.0, 12.0];

    let x = vec![0.0, 1.0];
    let y = vec![3.0, 4.0, 6.0]; // Spacing does not need to be uniform!
    let z = vec![1.0, 2.0];
    assert!(GriddedData::new([x, y, z], data).is_ok());
    ```

    If the data is available in matrix form, consider using the `CartesianIndices` iterator from the [`cart_lin`](https://crates.io/crates/cart_lin)
    library (which is a dependency of this library anyway):

    ```
    use cart_lin::CartesianIndices;
    use gridded_data::GriddedData;

    // 2 x 2 matrix
    // 1.0 3.0
    // 2.0 4.0
    let matrix = [[1.0, 3.0], [2.0, 4.0]];
    let mut data = Vec::with_capacity(4);
    for idx in CartesianIndices::new([2, 2]) {
        data.push(matrix[idx[0]][idx[1]]);
    }
    let x = vec![0.0, 1.0];
    let y = vec![1.0, 2.0];
    assert_eq!(data, vec![1.0, 3.0, 2.0, 4.0]);
    assert!(GriddedData::new([x, y], data).is_ok());
    ```

    As an alternative for two-dimensional data, consider using [`GriddedData::from_slice_of_slices`].
     */
    pub fn new(axes: [Vec<f64>; N], data: Vec<f64>) -> Result<Self, Error> {
        for (idx, axis) in axes.iter().enumerate() {
            // Assert that the input data for the axes is strictly monotonic increasing
            if !(axis.windows(2).all(|w| w[0] < w[1])) {
                return Err(Error(format!(
                    "data for axis {idx} is not strictly monotonic increasing"
                )));
            }

            // Assert that the axis length is at least 2
            if axis.len() < 2 {
                return Err(Error(format!("axis must contain at least two values")));
            }
        }

        // The length of the data must be equal to the product of all axes lengths
        if data.len() != axes.iter().fold(1, |acc, axis| acc * axis.len()) {
            return Err(Error(
                "length of the data must be equal to the product of all axes lengths".into(),
            ));
        }

        return Ok(Self { axes, data });
    }

    /**
    Access the raw data of the underlying axes.
     */
    pub fn axes(&self) -> &[Vec<f64>; N] {
        return &self.axes;
    }

    /**
    Access the raw underlying data.
     */
    pub fn data(&self) -> &[f64] {
        return self.data.as_slice();
    }

    /**
    Access the raw underlying data mutably.
     */
    pub fn data_mut(&mut self) -> &mut [f64] {
        return self.data.as_mut_slice();
    }

    /**
    Access the raw underlying data via cartesian indices.

    This function uses `cart_to_lin` from the [cart_lin](https://crates.io/crates/cart_lin) crate to convert the given cartesian
    indices into a linear index and use that to access [`GriddedData::data`]. If the given cartesian index would
    result in an out-of bounds access, `None` is returned instead.
     */
    pub fn get_cart(&self, indices: &[usize; N]) -> Option<&f64> {
        let index = cart_to_lin(indices, &self.axes_len())?;
        return self.data.get(index);
    }

    /**
    Access the raw underlying data mutably via cartesian indices.

    This function uses `cart_to_lin` from the [cart_lin](https://crates.io/crates/cart_lin) crate to convert the given cartesian
    indices into a linear index and use that to access [`GriddedData::data_mut`]. If the given cartesian index would
    result in an out-of bounds access, `None` is returned instead.
     */
    pub fn get_cart_mut(&mut self, indices: &[usize; N]) -> Option<&mut f64> {
        let index = cart_to_lin(indices, &self.axes_len())?;
        return self.data.get_mut(index);
    }

    /**
    Like [`GriddedData::get_cart`], but does not perform any bound checking and is therefore unsafe.
     */
    pub unsafe fn get_cart_unchecked(&self, indices: &[usize; N]) -> &f64 {
        let index = cart_to_lin_unchecked(indices, &self.axes_len());
        return unsafe { self.data.get_unchecked(index) };
    }

    /**
    Like [`GriddedData::get_cart_mut`], but does not perform any bound checking and is therefore unsafe.
     */
    pub unsafe fn get_cart_unchecked_mut(&mut self, indices: &[usize; N]) -> &mut f64 {
        let index = cart_to_lin_unchecked(indices, &self.axes_len());
        return unsafe { self.data.get_unchecked_mut(index) };
    }

    /**
    Returns the lengths of all axes.
     */
    pub fn axes_len(&self) -> [usize; N] {
        let mut bounds = [0; N];
        for (bound, axis) in bounds.iter_mut().zip(self.axes()) {
            *bound = axis.len();
        }
        return bounds;
    }

    /**
    Access the raw underlying data via a linear index.

    This is a convenience wrapper around `self.data().get()`.
     */
    pub fn get_lin(&self, index: usize) -> Option<&f64> {
        return self.data.get(index);
    }

    /**
    Access the raw underlying data mutably via a linear index.

    This is a convenience wrapper around `self.get_mut().get()`.
     */
    pub fn get_lin_mut(&mut self, index: usize) -> Option<&mut f64> {
        return self.data.get_mut(index);
    }

    /**
    Like [`GriddedData::get_lin`], but does not perform bounds checking and is therefore unsafe.
     */
    pub unsafe fn get_lin_unchecked(&self, index: usize) -> &f64 {
        return unsafe { self.data.get_unchecked(index) };
    }

    /**
    Like [`GriddedData::get_lin_mut`], but does not perform bounds checking and is therefore unsafe.
     */
    pub unsafe fn get_lin_unchecked_mut(&mut self, index: usize) -> &mut f64 {
        return unsafe { self.data.get_unchecked_mut(index) };
    }

    /**
    Returns the axes values which correspond to the given cartesian index.

    ```
    use gridded_data::GriddedData;

    let grid = GriddedData::new([vec![0.0, 1.0], vec![1.0, 2.0]], vec![1.0, 3.0, 2.0, 4.0]).expect("valid inputs");
    assert_eq!(grid.axes_values_cart(&[1, 0]).expect("valid index"), [1.0, 1.0]);
    ```
     */
    pub fn axes_values_cart(&self, indices: &[usize; N]) -> Option<[f64; N]> {
        let mut point = [0.0; N];
        for ((pt, index), axis) in point.iter_mut().zip(indices.iter()).zip(self.axes.iter()) {
            *pt = *axis.get(*index)?;
        }
        return Some(point);
    }

    /**
    Like [`GriddedData::axes_values_cart`], but does not perform bounds checking and is therefore unsafe.
     */
    pub unsafe fn axes_values_cart_unchecked(&self, indices: &[usize; N]) -> [f64; N] {
        let mut point = [0.0; N];
        for ((pt, index), axis) in point.iter_mut().zip(indices.iter()).zip(self.axes.iter()) {
            *pt = *unsafe { axis.get_unchecked(*index) };
        }
        return point;
    }

    /**
    Returns the grid cell containing the given point.

    If the given point is within the grid, this function finds the indices of the node interval of each axis which contain
    the corresponding point value.

    # Examples

    ```
    use gridded_data::GriddedData;

    let data = vec![1.0, 7.0, 3.0, 9.0, 5.0, 11.0, 2.0, 8.0, 4.0, 10.0, 6.0, 12.0];
    let x = vec![0.0, 1.0];
    let y = vec![3.0, 4.0, 6.0]; // Spacing does not need to be uniform!
    let z = vec![1.0, 2.0];
    let grid = GriddedData::new([x, y, z], data).expect("valid input");

    /*
    The x-coordinate of the given point is between the 0th and the 1st node
    of the x-axis, the y-coordinate is between the 1st and the 2nd node and the
    z-coordinate is between the 0th and the 1st node.
     */
    assert_eq!(grid.cell_bounds(&[0.5, 4.1, 1.5]), Some([[0, 1], [1, 2], [0, 1]]));

    /*
    x-coordinate directly on 0th node
    y-coordinate directly on 1st node -> Two different intervals could be returned:
    [0, 1] or [1, 2]. This function returns the "right-side" interval
    z-coordinate directly on 1st node
     */
    assert_eq!(grid.cell_bounds(&[0.0, 4.0, 2.0]), Some([[0, 1], [1, 2], [0, 1]]));

    // y-coordinate outside grid -> Return None.
    assert_eq!(grid.cell_bounds(&[0.5, 6.1, 1.5]), None);
    ```
     */
    pub fn cell_bounds(&self, point: &[f64; N]) -> Option<[[usize; 2]; N]> {
        if !self.contains(point) {
            return None;
        }
        return Some(self.cell_bounds_raw(point));
    }

    fn cell_bounds_raw(&self, point: &[f64; N]) -> [[usize; 2]; N] {
        let mut bounds = [[0, 0]; N];
        for (index, (p, axis)) in bounds.iter_mut().zip(point.iter().zip(self.axes().iter())) {
            let partition_point = axis.as_slice().partition_point(|&val| val <= *p);
            match partition_point.checked_sub(1) {
                Some(pp_minus_one) => {
                    // SAFETY: axis contains at least one element (checked in constructor)
                    if partition_point == axis.len() {
                        if unsafe { *axis.last().unwrap_unchecked() } == *p {
                            *index = [pp_minus_one - 1, pp_minus_one]
                        } else {
                            *index = [pp_minus_one, pp_minus_one]
                        }
                    } else {
                        *index = [pp_minus_one, partition_point]
                    }
                }
                None => *index = [partition_point, partition_point],
            }
        }
        return bounds;
    }

    /**
    Returns whether the given point is inside the grid or not.

    # Examples
    ```
    use gridded_data::GriddedData;

    let grid = GriddedData::new([vec![0.0, 1.0], vec![1.0, 2.0]], vec![1.0, 3.0, 2.0, 4.0]).expect("valid inputs");

    assert!(grid.contains(&[0.5, 1.2]));
    assert!(!grid.contains(&[1.5, 1.2]));
    ```
     */
    pub fn contains(&self, point: &[f64; N]) -> bool {
        return self.axes().iter().zip(point.iter()).all(|(axis, val)| {
            // SAFETY: axis has at least one element -> This is checked in the grid constructor
            unsafe {
                axis.first().unwrap_unchecked() <= val && axis.last().unwrap_unchecked() >= val
            }
        });
    }

    /**
    Performs a nearest-neighbor interpolation for the given point.

    This function first finds the "cell" (see [`GriddedData::cell_bounds`]) which contains
    the given point. Afterwards, it calculates the [Euclidian distance](https://en.wikipedia.org/wiki/Euclidean_distance)
    of the point to each cell node and identifies the cartesian indices of the nearest node.
    Finally, it returns the corresponding `data` value.

    In case the point is outside the grid, the only difference is that there is just one node per dimension for each
    point value outside of the respective axis.

    # Examples
    ```
    use gridded_data::GriddedData;

    let data = vec![1.0, 7.0, 3.0, 9.0, 5.0, 11.0, 2.0, 8.0, 4.0, 10.0, 6.0, 12.0];
    let x = vec![0.0, 1.0];
    let y = vec![3.0, 4.0, 6.0]; // Spacing does not need to be uniform!
    let z = vec![1.0, 2.0];
    let grid = GriddedData::new([x, y, z], data).expect("valid input");

    /*
    This point is inside the cell with the bounds x = [0.0, 1.0], y = [4.0, 6.0], z = [1.0, 2.0]
    Within this cell, it is closest to the node [1.0, 4.0, 2.0], which corresponds to the data value 10.0
     */
    assert_eq!(grid.nearest_neighbor_interp(&[0.9, 4.2, 1.7]), 10.0);

    /*
    The x-coordinate of this point is outside the grid -> Extrapolation by only considering the x-node with value 0.0
     */
    assert_eq!(grid.nearest_neighbor_interp(&[-1.0, 4.2, 1.7]), 9.0);
    ```
     */
    pub fn nearest_neighbor_interp(&self, point: &[f64; N]) -> f64 {
        let mut cell_bounds = self.cell_bounds_raw(point);

        // CartesianIndices::from_bounds_unchecked expects a half-open interval, hence we need to add 1 to the upper limits
        add_one_to_upper(&mut cell_bounds);

        // Temporary buffer for the current corner point
        let mut corner_point_indices = [0; N];

        // We are looking for the closest corner point
        let mut closest_corner_point = [0; N];

        // Start values of the iteration
        let mut minimum_distance = std::f64::INFINITY;

        // Look up the cell corner points and calculate the distance to point
        for cart_prod in CartesianIndices::from_bounds_unchecked(cell_bounds).into_iter() {
            // Populate the corner point indices for all axes
            for (corner_point_index, cart_prod_index) in
                corner_point_indices.iter_mut().zip(cart_prod.iter())
            {
                *corner_point_index = *cart_prod_index;
            }

            // Now calculate the distance between this point and the input
            let distance: f64 = self
                .axes()
                .iter()
                .zip(corner_point_indices.iter())
                .zip(point.iter())
                .map(|((axis, corner_point_index), point_value)| {
                    // SAFETY: This value exists for sure since it is build from cell_bounds
                    let corner_point_value = unsafe { *axis.get_unchecked(*corner_point_index) };
                    (corner_point_value - *point_value).powi(2)
                })
                .sum();
            if distance < minimum_distance {
                closest_corner_point = corner_point_indices;
                minimum_distance = distance;
            }
        }

        // SAFETY: The closest corner point is one of the cell corner points and therefore exists for sure.
        return unsafe { *self.get_cart_unchecked(&closest_corner_point) };
    }

    /**
    Performs a multivariate linear interpolation for the given point.

    This function performs a generalized bilinear interpolation and calculates the
    return value as
    ```math
    output = num / denom
    ```
    with
    ```math
    num = [ P11..11 * (a2-a) * (b2-b) * ... (n2-n) + P21..11 * (a-a1) * (b2-b) * ... (n2-n) + ... + P22..22 * (a-a1) * (b-b1) * ... (n-n1) ]
    denom = [ (a2-a1) * (b2-b1) * ... (n2-n1) ]
    ```
    - 1 is the index of the smaller value, 2 is the index of the larger value
    - P11..11 to P22..22 are the node values of the multidimensional grid
    - a1, a2 are the two data values along the first axis, b1, b2 along the second axis and so on.
    - a, b, ..., n are the values of the input point for the corresponding dimensions.

    If the point is outside the grid, a [nearest-neighbor extrapolation](GriddedData::nearest_neighbor_interp) is performed instead.

    # Examples
    ```
    use gridded_data::GriddedData;

    let data = vec![1.0, 7.0, 3.0, 9.0, 5.0, 11.0, 2.0, 8.0, 4.0, 10.0, 6.0, 12.0];
    let x = vec![0.0, 1.0];
    let y = vec![3.0, 4.0, 6.0]; // Spacing does not need to be uniform!
    let z = vec![1.0, 2.0];
    let grid = GriddedData::new([x, y, z], data).expect("valid input");

    // Values within the grid (tolerance accounts for limited floating-point accuracy)
    assert!((grid.linear_interp(&[0.9, 4.2, 1.7]) - 8.3).abs() < 1e-10);
    assert!((grid.linear_interp(&[0.8, 3.7, 1.2]) - 4.4).abs() < 1e-10);

    // Value outside grid -> Nearest-neighbor extrapolation
    assert_eq!(grid.linear_interp(&[-1.0, 4.2, 1.7]), 9.0);
    ```
     */
    pub fn linear_interp(&self, point: &[f64; N]) -> f64 {
        let cell_bounds = match self.cell_bounds(&point) {
            Some(limits) => limits,
            None => {
                return {
                    // Extrapolation with the nearest-neighbor method
                    self.nearest_neighbor_interp(&point)
                };
            }
        };

        // Calculate the denominator [ (a2-a1) * (b2-b1) * ... (n2-n1) ]
        let mut denominator = 1.0;
        for (limit, axis) in cell_bounds.iter().zip(self.axes()) {
            // SAFETY: the function contains checked that all dimensions of point are inside the grid.
            // This also makes sure that the bracket value is always positive, since limit[0] < limit[1]
            // and the axis is strictly monotonic increasing.
            let lower = unsafe { *axis.get_unchecked(limit[0]) };
            let upper = unsafe { *axis.get_unchecked(limit[1]) };
            denominator = denominator * (upper - lower);
        }

        // CartesianIndices::from_bounds_unchecked expects a half-open interval, hence we need to add 1 to the upper limits
        let mut adj_cell_bounds = cell_bounds.clone();
        add_one_to_upper(&mut adj_cell_bounds);

        // Calculate the numerator
        let mut numerator = 0.0;
        for grid_corner in CartesianIndices::from_bounds_unchecked(adj_cell_bounds).into_iter() {
            // Get the indices of the diagonally opposing point
            // Those values are needed to calculate the brackets.
            let mut opposing_grid_corner = [0; N];
            opposing_grid_corner
                .iter_mut()
                .zip(cell_bounds.iter())
                .zip(grid_corner.iter())
                .for_each(|((op_corner, limit), corner)| {
                    *op_corner = limit[0] * usize::from(*corner == limit[1])
                        + limit[1] * usize::from(*corner == limit[0]);
                });

            // SAFETY: grid_corner is generated from cell_bounds which is fully inside the grid.
            let mut product = unsafe { *self.get_cart_unchecked(&grid_corner) };

            // Calculate the bracket values
            for (((pt_val, corner_idx), limit), axis) in point
                .iter()
                .zip(grid_corner.iter())
                .zip(cell_bounds.iter())
                .zip(self.axes().iter())
            {
                // If this is true, we need to calculate (x2-x), otherwise (x-x1)
                // SAFETY: All indices of cell_bounds are guaranteed to be in bounds for their respective axis.
                if *corner_idx == limit[0] {
                    product = product * (unsafe { axis.get_unchecked(limit[1]) } - *pt_val);
                } else {
                    product = product * (*pt_val - unsafe { axis.get_unchecked(limit[0]) });
                }
            }

            numerator += product;
        }

        return numerator / denominator;
    }
}

/**
Add an offset of + 1 to the upper bound of each dimension
*/
fn add_one_to_upper<const N: usize>(bounds: &mut [[usize; 2]; N]) {
    for [_, upper] in bounds.iter_mut() {
        *upper += 1;
    }
}

impl<const N: usize> std::ops::Index<&[usize; N]> for GriddedData<N> {
    type Output = f64;

    fn index(&self, index: &[usize; N]) -> &Self::Output {
        return self.get_cart(index).expect("out-of-bounds access");
    }
}

impl<const N: usize> std::ops::IndexMut<&[usize; N]> for GriddedData<N> {
    fn index_mut(&mut self, index: &[usize; N]) -> &mut Self::Output {
        return self.get_cart_mut(index).expect("out-of-bounds access");
    }
}

impl<const N: usize> std::ops::Index<usize> for GriddedData<N> {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        return &self.data[index];
    }
}

impl<const N: usize> std::ops::IndexMut<usize> for GriddedData<N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        return &mut self.data[index];
    }
}

impl GriddedData<2> {
    /**
    Create a 2-dimensional [`GriddedData`] using a slices-of-slices matrix as data.

    # Examples
    ```
    use gridded_data::GriddedData;

    let row1 = [1.0, 3.0];
    let row2 = [2.0, 4.0];
    let data = [row1.as_slice(), row2.as_slice()];
    let grid = GriddedData::from_slice_of_slices([vec![0.0, 1.0], vec![1.0, 2.0]], data.as_slice()).expect("valid inputs");
    assert_eq!(grid.data(), &[1.0, 2.0, 3.0, 4.0]);
    ```
     */
    pub fn from_slice_of_slices(axes: [Vec<f64>; 2], data: &[&[f64]]) -> Result<Self, Error> {
        let length_second_axis = axes[1].len();
        if data.len() != length_second_axis {
            return Err(Error(format!(
                "length {} of outer slice is not equal to the length {} of the second axis",
                data.len(),
                length_second_axis
            )));
        }

        let length_first_axis = axes[0].len();
        for (idx, inner_slice) in data.iter().enumerate() {
            if inner_slice.len() != axes[0].len() {
                return Err(Error(format!(
                    "length {} of inner slice {} is not equal to the length {} of the first axis",
                    data.len(),
                    idx,
                    length_first_axis
                )));
            }
        }

        // Concatenate the data into a vector
        let mut vec = Vec::with_capacity(length_first_axis * length_second_axis);
        for first in 0..length_first_axis {
            for second in 0..length_second_axis {
                // SAFETY: We checked before that the length of the outer slice is equal to length_second_axis
                // and that the length of each inner slice is equal to length_first_axis
                vec.push(unsafe { *data.get_unchecked(second).get_unchecked(first) });
            }
        }

        return Self::new(axes, vec);
    }
}

/**
Error for gridded data.
 */
#[derive(Debug, Clone)]
pub struct Error(pub String);

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        return self.0.fmt(f);
    }
}

impl std::error::Error for Error {}

#[cfg(feature = "nalgebra")]
mod nalgebra_impl {
    use super::*;
    use nalgebra::{Dim, Matrix, RawStorage};

    impl GriddedData<2> {
        /**
        Create a 2-dimensional [`GriddedData`] using a [nalgebra](https://crates.io/crates/nalgebra) matrix.

        # Examples
        ```
        use gridded_data::GriddedData;
        use nalgebra::Matrix2;

        let row1 = [1.0, 3.0];
        let row2 = [2.0, 4.0];
        let data = Matrix2::from_vec(vec![1.0, 3.0, 2.0, 4.0]); // nalgebra is column-major!
        let grid = GriddedData::from_nalgebra_matrix([vec![0.0, 1.0], vec![1.0, 2.0]], &data).expect("valid inputs");
        assert_eq!(grid.data(), &[1.0, 2.0, 3.0, 4.0]); // The data within grid is stored row-major
        ```
        */
        pub fn from_nalgebra_matrix<R: Dim, C: Dim, S: RawStorage<f64, R, C>>(
            axes: [Vec<f64>; 2],
            matrix: &Matrix<f64, R, C, S>,
        ) -> Result<Self, Error> {
            if axes[0].len() != matrix.nrows() {
                return Err(Error(
                    "number of rows must be equal to the length of the first axis".to_string(),
                ));
            }
            if axes[1].len() != matrix.ncols() {
                return Err(Error(
                    "number of columns must be equal to the length of the second axis".to_string(),
                ));
            }

            let mut vector = Vec::with_capacity(matrix.nrows() * matrix.ncols());
            for row in matrix.row_iter() {
                for value in row.iter() {
                    vector.push(*value);
                }
            }

            return Self::new(axes, vector);
        }
    }
}
