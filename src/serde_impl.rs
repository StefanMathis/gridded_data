//! Serialization and deserialization of [`GriddedData`] using the [serde](https://crates.io/crates/serde) library.
//! Since serde currently does not support serialization / deserialization of generic arrays `[T; N]`, this module
//! contains manual implementations for 1 up to 16 dimensions.

use serde::{Deserialize, Serialize};

use crate::GriddedData;

macro_rules! impl_gridded_data {
    ($n:literal) => {
        impl Serialize for GriddedData<$n> {
            fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
            where
                S: serde::Serializer,
            {
                #[derive(Serialize)]
                #[serde(rename = "GriddedData")]
                struct GriddedDataSer<'a> {
                    axes: &'a [Vec<f64>; $n],
                    data: &'a Vec<f64>,
                }

                let grid = GriddedDataSer {
                    axes: &self.axes,
                    data: &self.data,
                };

                grid.serialize(serializer)
            }
        }

        impl<'de> Deserialize<'de> for GriddedData<$n> {
            fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
            where
                D: serde::Deserializer<'de>,
            {
                #[derive(Deserialize)]
                #[serde(rename = "GriddedData")]
                struct GriddedDataDe {
                    axes: [Vec<f64>; $n],
                    data: Vec<f64>,
                }

                let data = GriddedDataDe::deserialize(deserializer)?;
                GriddedData::new(data.axes, data.data).map_err(serde::de::Error::custom)
            }
        }
    };
}

impl_gridded_data!(1);
impl_gridded_data!(2);
impl_gridded_data!(3);
impl_gridded_data!(4);
impl_gridded_data!(5);
impl_gridded_data!(6);
impl_gridded_data!(7);
impl_gridded_data!(8);
impl_gridded_data!(9);
impl_gridded_data!(10);
impl_gridded_data!(11);
impl_gridded_data!(12);
impl_gridded_data!(13);
impl_gridded_data!(14);
impl_gridded_data!(15);
impl_gridded_data!(16);
