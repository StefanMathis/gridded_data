use gridded_data::GriddedData;

#[test]
fn test_from_slice_of_slices() {
    {
        // From vec of vecs
        let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let y = vec![0.0, 1.0];
        let first = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let second = vec![-1.0, -2.0, -3.0, -4.0, -5.0];
        let data = [first.as_slice(), second.as_slice()];

        let grid = GriddedData::from_slice_of_slices([x, y], &data).unwrap();
        assert_eq!(*grid.get_cart(&[0, 0]).unwrap(), 1.0);
        assert_eq!(*grid.get_cart(&[0, 1]).unwrap(), -1.0);
        assert_eq!(*grid.get_cart(&[1, 0]).unwrap(), 2.0);
        assert_eq!(*grid.get_cart(&[1, 1]).unwrap(), -2.0);
    }
}

#[test]
fn test_cell_bounds_1d() {
    let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let data = vec![0.0; x.len()];
    let grid = GriddedData::<1>::new([x], data).unwrap();
    assert!(grid.cell_bounds(&[-1.0]).is_none());
    assert_eq!(grid.cell_bounds(&[0.5]).unwrap()[0], [0, 1]);
    assert_eq!(grid.cell_bounds(&[3.5]).unwrap()[0], [3, 4]);
    assert!(grid.cell_bounds(&[4.5]).is_none());
}

#[test]
fn test_cell_bounds_2d() {
    {
        // 3 x 3 matrix
        let x = vec![0.0, 1.0, 2.0];
        let y = vec![0.0, 1.0, 2.0];
        let data = vec![0.0; x.len() * y.len()];
        let grid = GriddedData::<2>::new([x, y], data).unwrap();
        assert!(grid.cell_bounds(&[-1.0, 1.5]).is_none());
        assert_eq!(grid.cell_bounds(&[0.0, 2.0]).unwrap(), [[0, 1], [1, 2]]);
        assert!(grid.cell_bounds(&[1.1, 5.0]).is_none());
        assert_eq!(grid.cell_bounds(&[0.5, 0.5]).unwrap(), [[0, 1], [0, 1]]);
    }
    {
        // 2 x 5 matrix with five columns and two rows
        let x = vec![0.0, 1.0];
        let y = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let data = vec![0.0; x.len() * y.len()];
        let grid = GriddedData::<2>::new([x, y], data).unwrap();

        assert!(grid.cell_bounds(&[1.1, 3.0]).is_none());
        assert!(grid.cell_bounds(&[3.0, 1.1]).is_none());
        assert!(grid.cell_bounds(&[3.1, -1.1]).is_none());
        assert!(grid.cell_bounds(&[2.9, 0.5]).is_none());
        assert!(grid.cell_bounds(&[-0.5, 5.5]).is_none());
    }
}

#[test]
fn test_nearest_neighbor_interp_1d() {
    let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let data = vec![0.0, 2.0, 4.0, 2.0, 0.0];
    let grid = GriddedData::<1>::new([x], data).unwrap();
    assert_eq!(grid.nearest_neighbor_interp(&[-1.0]), 0.0);
    assert_eq!(grid.nearest_neighbor_interp(&[0.2]), 0.0);
    assert_eq!(grid.nearest_neighbor_interp(&[0.6]), 2.0);
    assert_eq!(grid.nearest_neighbor_interp(&[1.6]), 4.0);
    assert_eq!(grid.nearest_neighbor_interp(&[3.1]), 2.0);
    assert_eq!(grid.nearest_neighbor_interp(&[3.5]), 2.0);
    assert_eq!(grid.nearest_neighbor_interp(&[3.8]), 0.0);
    assert_eq!(grid.nearest_neighbor_interp(&[4.5]), 0.0);
}

#[test]
fn test_nearest_neighbor_interp_2d() {
    // 3 x 3 matrix
    // [0 3 6]
    // [1 4 7]
    // [2 5 8]
    let x = vec![0.0, 1.0, 2.0];
    let y = vec![0.0, 1.0, 2.0];
    let data = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let grid = GriddedData::new([x, y], data).unwrap();

    // Values at the grid corners
    assert_eq!(grid.nearest_neighbor_interp(&[0.0, 0.0]), 0.0);
    assert_eq!(grid.nearest_neighbor_interp(&[1.0, 0.0]), 3.0);
    assert_eq!(grid.nearest_neighbor_interp(&[2.0, 0.0]), 6.0);
    assert_eq!(grid.nearest_neighbor_interp(&[0.0, 1.0]), 1.0);
    assert_eq!(grid.nearest_neighbor_interp(&[1.0, 1.0]), 4.0);
    assert_eq!(grid.nearest_neighbor_interp(&[2.0, 1.0]), 7.0);
    assert_eq!(grid.nearest_neighbor_interp(&[0.0, 2.0]), 2.0);
    assert_eq!(grid.nearest_neighbor_interp(&[1.0, 2.0]), 5.0);
    assert_eq!(grid.nearest_neighbor_interp(&[2.0, 2.0]), 8.0);

    // Inter- and extrapolation
    assert_eq!(grid.nearest_neighbor_interp(&[-1.0, 0.0]), 0.0);
    assert_eq!(grid.nearest_neighbor_interp(&[0.8, 0.8]), 4.0);
    assert_eq!(grid.nearest_neighbor_interp(&[1.2, 1.2]), 4.0);
    assert_eq!(grid.nearest_neighbor_interp(&[1.2, 0.4]), 3.0);
    assert_eq!(grid.nearest_neighbor_interp(&[1.2, 1.7]), 5.0);
    assert_eq!(grid.nearest_neighbor_interp(&[2.5, 1.7]), 8.0);
    assert_eq!(grid.nearest_neighbor_interp(&[1.5, 1.7]), 5.0);
    assert_eq!(grid.nearest_neighbor_interp(&[0.1, 1.7]), 2.0);
    assert_eq!(grid.nearest_neighbor_interp(&[1.7, 0.1]), 6.0);
}

#[test]
fn test_nearest_neighbor_interp_3d() {
    // Cube
    // [0 4] [1 5]
    // [2 6] [3 7]
    let x = vec![0.0, 1.0];
    let y = vec![2.0, 3.0];
    let z = vec![4.0, 5.0];
    let data = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
    let grid = GriddedData::new([x, y, z], data).unwrap();

    // Values at the cube corners
    assert_eq!(grid.nearest_neighbor_interp(&[0.0, 2.0, 4.0]), 0.0);
    assert_eq!(grid.nearest_neighbor_interp(&[1.0, 2.0, 4.0]), 4.0);
    assert_eq!(grid.nearest_neighbor_interp(&[0.0, 3.0, 4.0]), 2.0);
    assert_eq!(grid.nearest_neighbor_interp(&[1.0, 3.0, 4.0]), 6.0);
    assert_eq!(grid.nearest_neighbor_interp(&[1.0, 3.0, 5.0]), 7.0);
}

#[test]
fn test_linear_interp_1d() {
    let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let data = vec![0.0, 2.0, 4.0, 2.0, 0.0];
    let grid = GriddedData::<1>::new([x], data).unwrap();
    approx::assert_abs_diff_eq!(grid.linear_interp(&[-1.0]), 0.0, epsilon = 1e-15);
    approx::assert_abs_diff_eq!(grid.linear_interp(&[0.2]), 0.4, epsilon = 1e-15);
    approx::assert_abs_diff_eq!(grid.linear_interp(&[0.6]), 1.2, epsilon = 1e-15);
    approx::assert_abs_diff_eq!(grid.linear_interp(&[1.6]), 3.2, epsilon = 1e-15);
    approx::assert_abs_diff_eq!(grid.linear_interp(&[3.1]), 1.8, epsilon = 1e-15);
    approx::assert_abs_diff_eq!(grid.linear_interp(&[3.5]), 1.0, epsilon = 1e-15);
    approx::assert_abs_diff_eq!(grid.linear_interp(&[3.8]), 0.4, epsilon = 1e-15);
    approx::assert_abs_diff_eq!(grid.linear_interp(&[4.5]), 0.0, epsilon = 1e-15);
}

#[test]
fn test_linear_interp_2d() {
    // 3 x 3 matrix
    // [0 3 0]
    // [1 4 1]
    // [2 4 2]
    let x = vec![0.0, 1.0, 2.0];
    let y = vec![0.0, 1.0, 2.0];
    let data = vec![0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 0.0, 1.0, 2.0];
    let grid = GriddedData::new([x, y], data).unwrap();

    // Values at the grid corners
    assert_eq!(grid.linear_interp(&[0.0, 0.0]), 0.0);
    assert_eq!(grid.linear_interp(&[1.0, 0.0]), 3.0);
    assert_eq!(grid.linear_interp(&[2.0, 0.0]), 0.0);
    assert_eq!(grid.linear_interp(&[0.0, 1.0]), 1.0);
    assert_eq!(grid.linear_interp(&[1.0, 1.0]), 4.0);
    assert_eq!(grid.linear_interp(&[2.0, 1.0]), 1.0);
    assert_eq!(grid.linear_interp(&[0.0, 2.0]), 2.0);
    assert_eq!(grid.linear_interp(&[1.0, 2.0]), 4.0);
    assert_eq!(grid.linear_interp(&[2.0, 2.0]), 2.0);

    // Inter- and extrapolation
    approx::assert_abs_diff_eq!(grid.linear_interp(&[-1.0, 0.0]), 0.0, epsilon = 1e-15);

    approx::assert_abs_diff_eq!(grid.linear_interp(&[0.8, 0.8]), 3.2, epsilon = 1e-15);
    approx::assert_abs_diff_eq!(grid.linear_interp(&[1.2, 0.8]), 3.2, epsilon = 1e-15);

    approx::assert_abs_diff_eq!(grid.linear_interp(&[0.0, 0.5]), 0.5, epsilon = 1e-15);
    approx::assert_abs_diff_eq!(grid.linear_interp(&[0.0, 1.5]), 1.5, epsilon = 1e-15);

    approx::assert_abs_diff_eq!(grid.linear_interp(&[2.0, 0.5]), 0.5, epsilon = 1e-15);
    approx::assert_abs_diff_eq!(grid.linear_interp(&[2.0, 1.5]), 1.5, epsilon = 1e-15);

    approx::assert_abs_diff_eq!(grid.linear_interp(&[1.0, 0.5]), 3.5, epsilon = 1e-15);
    approx::assert_abs_diff_eq!(grid.linear_interp(&[1.0, 1.5]), 4.0, epsilon = 1e-15);

    approx::assert_abs_diff_eq!(grid.linear_interp(&[0.5, 0.0]), 1.5, epsilon = 1e-15);
    approx::assert_abs_diff_eq!(grid.linear_interp(&[1.5, 0.0]), 1.5, epsilon = 1e-15);

    approx::assert_abs_diff_eq!(grid.linear_interp(&[1.5, 1.5]), 2.75, epsilon = 1e-15);
    approx::assert_abs_diff_eq!(grid.linear_interp(&[0.5, 1.5]), 2.75, epsilon = 1e-15);

    approx::assert_abs_diff_eq!(grid.linear_interp(&[0.25, 0.75]), 1.5, epsilon = 1e-15);
}

#[test]
fn test_linear_interp_3d() {
    // Cube
    // [0 1] [1 2]
    // [0 1] [1 2]
    let x = vec![0.0, 1.0];
    let y = vec![0.0, 1.0];
    let z = vec![0.0, 1.0];
    let data = vec![0.0, 1.0, 0.0, 1.0, 1.0, 2.0, 1.0, 2.0];
    let grid = GriddedData::new([x, y, z], data).unwrap();

    assert_eq!(grid.linear_interp(&[0.5, 0.5, 0.5]), 1.0);
    assert_eq!(grid.linear_interp(&[0.0, 0.5, 0.5]), 0.5);
    assert_eq!(grid.linear_interp(&[0.5, 0.0, 0.5]), 1.0);
    assert_eq!(grid.linear_interp(&[0.0, 0.0, 0.5]), 0.5);
    approx::assert_abs_diff_eq!(grid.linear_interp(&[0.9, 0.9, 0.9]), 1.8, epsilon = 1e-15);
    assert_eq!(grid.linear_interp(&[0.0, 0.0, -1.0]), 0.0);
    assert_eq!(grid.linear_interp(&[0.0, 2.0, -1.0]), 0.0);
    assert_eq!(grid.linear_interp(&[2.0, 2.0, 2.0]), 2.0);
    assert_eq!(grid.linear_interp(&[2.0, 2.0, -1.0]), 1.0);
}
