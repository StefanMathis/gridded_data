use gridded_data::GriddedData;
use nalgebra::DMatrix;

#[test]
fn test_nearest_neighbor_interp_2d() {
    {
        let x = vec![0.0, 1.0];
        let y = vec![0.0, 1.0, 2.0, 3.0, 4.0];

        // [1.5  1.5 ]
        // [2.5  2.5 ]
        // [3.5  3.5 ]
        // [1.25 1.25]
        // [0.0  0.0 ]
        let mut matrix = DMatrix::from_element(x.len(), y.len(), 0.0);
        matrix[(0, 0)] = 1.5;
        matrix[(1, 0)] = 1.5;
        matrix[(0, 1)] = 2.5;
        matrix[(1, 1)] = 2.5;
        matrix[(0, 2)] = 3.5;
        matrix[(1, 2)] = 3.5;
        matrix[(0, 3)] = 1.25;
        matrix[(1, 3)] = 1.25;

        let grid = GriddedData::from_nalgebra_matrix([x, y], &matrix).unwrap();

        assert_eq!(grid.nearest_neighbor_interp(&[0.0, 0.0]), 1.5);
        assert_eq!(grid.nearest_neighbor_interp(&[1.0, 0.0]), 1.5);

        assert_eq!(grid.nearest_neighbor_interp(&[0.0, 1.0]), 2.5);
        assert_eq!(grid.nearest_neighbor_interp(&[1.0, 1.0]), 2.5);

        assert_eq!(grid.nearest_neighbor_interp(&[0.0, 2.0]), 3.5);
        assert_eq!(grid.nearest_neighbor_interp(&[1.0, 2.0]), 3.5);

        assert_eq!(grid.nearest_neighbor_interp(&[0.0, 3.0]), 1.25);
        assert_eq!(grid.nearest_neighbor_interp(&[1.0, 3.0]), 1.25);

        assert_eq!(grid.nearest_neighbor_interp(&[0.0, 4.0]), 0.0);
        assert_eq!(grid.nearest_neighbor_interp(&[1.0, 4.0]), 0.0);
    }
}
