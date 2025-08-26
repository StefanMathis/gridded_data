use gridded_data::GriddedData;
use indoc::indoc;

#[test]
fn serialize_and_deserialize_1d() {
    let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let data = vec![0.0, 2.0, 4.0, 2.0, 0.0];
    let grid = GriddedData::<1>::new([x], data).unwrap();

    let txt = serde_yaml::to_string(&grid).unwrap();
    let cmp = indoc! {"
        ---
        axes:
          - - 0.0
            - 1.0
            - 2.0
            - 3.0
            - 4.0
        data:
          - 0.0
          - 2.0
          - 4.0
          - 2.0
          - 0.0
        "};
    assert_eq!(txt, cmp);

    let de_grid: GriddedData<1> = serde_yaml::from_str(&txt).unwrap();
    assert_eq!(
        grid.nearest_neighbor_interp(&[1.5]),
        de_grid.nearest_neighbor_interp(&[1.5])
    );
}

#[test]
fn serialize_and_deserialize_2d() {
    let x = vec![0.0, 1.0, 2.0];
    let y = vec![0.0, 1.0, 2.0];
    let data = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let grid = GriddedData::new([x, y], data).unwrap();

    let txt = serde_yaml::to_string(&grid).unwrap();
    let cmp = indoc! {"
        ---
        axes:
          - - 0.0
            - 1.0
            - 2.0
          - - 0.0
            - 1.0
            - 2.0
        data:
          - 0.0
          - 1.0
          - 2.0
          - 3.0
          - 4.0
          - 5.0
          - 6.0
          - 7.0
          - 8.0
        "};
    assert_eq!(txt, cmp);

    let de_grid: GriddedData<2> = serde_yaml::from_str(&txt).unwrap();
    assert_eq!(
        grid.nearest_neighbor_interp(&[1.5, 1.5]),
        de_grid.nearest_neighbor_interp(&[1.5, 1.5])
    );
}
