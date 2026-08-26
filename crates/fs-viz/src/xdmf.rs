//! Deterministic XDMF 3.0 / 2.0 descriptor emission for heavy dataset visualization.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.8`
//!
//! Generates clean, standard XDMF XML files with raw binary heavy-data containers
//! (strictly pure Rust, zero C / HDF5 dependency) for large-scale solution fields.

use super::vtu::{CellType, DataAssociation, DataValues, UnstructuredGrid, VtuError};
use std::fmt::Write as _;

/// Writer for deterministic XDMF representations.
pub struct XdmfWriter;

impl XdmfWriter {
    /// Generate the XDMF XML text and the raw binary companion buffer.
    pub fn write_xdmf_with_binary(
        grid: &UnstructuredGrid,
        binary_file_name: &str,
    ) -> Result<(String, Vec<u8>), VtuError> {
        grid.validate()?;

        let mut bin = Vec::with_capacity(16 * 1024);
        let mut xmf = String::with_capacity(4096);

        xmf.push_str("<?xml version=\"1.0\" ?>\n");
        xmf.push_str("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
        xmf.push_str("<Xdmf Version=\"3.0\">\n");
        xmf.push_str("  <Domain>\n");
        xmf.push_str("    <Grid Name=\"FrankenSimGrid\" GridType=\"Uniform\">\n");

        // 1. Topology
        let num_cells = grid.num_cells();
        let topology_type = if grid
            .cells_types
            .iter()
            .all(|&t| t == CellType::Tetra.type_id())
        {
            "Tetrahedron"
        } else if grid
            .cells_types
            .iter()
            .all(|&t| t == CellType::Triangle.type_id())
        {
            "Triangle"
        } else {
            "Mixed"
        };

        let conn_offset = bin.len();
        for &c in &grid.cells_connectivity {
            bin.extend_from_slice(&(c as u64).to_le_bytes());
        }
        let conn_len = grid.cells_connectivity.len();

        writeln!(
            xmf,
            "      <Topology TopologyType=\"{topology_type}\" NumberOfElements=\"{num_cells}\">"
        )
        .ok();
        writeln!(
            xmf,
            "        <DataItem DataType=\"UInt\" Precision=\"8\" Dimensions=\"{conn_len}\" Format=\"Binary\" Seek=\"{conn_offset}\">"
        )
        .ok();
        writeln!(xmf, "          {binary_file_name}").ok();
        xmf.push_str("        </DataItem>\n");
        xmf.push_str("      </Topology>\n");

        // 2. Geometry (Points)
        let num_points = grid.num_points();
        let geom_offset = bin.len();
        for p in &grid.points {
            bin.extend_from_slice(&p[0].to_le_bytes());
            bin.extend_from_slice(&p[1].to_le_bytes());
            bin.extend_from_slice(&p[2].to_le_bytes());
        }

        xmf.push_str("      <Geometry GeometryType=\"XYZ\">\n");
        writeln!(
            xmf,
            "        <DataItem DataType=\"Float\" Precision=\"8\" Dimensions=\"{num_points} 3\" Format=\"Binary\" Seek=\"{geom_offset}\">"
        )
        .ok();
        writeln!(xmf, "          {binary_file_name}").ok();
        xmf.push_str("        </DataItem>\n");
        xmf.push_str("      </Geometry>\n");

        // 3. Attributes (Arrays)
        for arr in &grid.arrays {
            let center = match arr.association {
                DataAssociation::PointData => "Node",
                DataAssociation::CellData => "Cell",
            };
            let attr_type = if arr.components == 1 {
                "Scalar"
            } else if arr.components == 3 {
                "Vector"
            } else {
                "Matrix"
            };

            let arr_offset = bin.len();
            let dim_str = match arr.association {
                DataAssociation::PointData => {
                    if arr.components == 1 {
                        format!("{num_points}")
                    } else {
                        format!("{num_points} {}", arr.components)
                    }
                }
                DataAssociation::CellData => {
                    if arr.components == 1 {
                        format!("{num_cells}")
                    } else {
                        format!("{num_cells} {}", arr.components)
                    }
                }
            };

            match &arr.values {
                DataValues::Float64(vals) => {
                    for &v in vals {
                        bin.extend_from_slice(&v.to_le_bytes());
                    }
                }
                DataValues::Float32(vals) => {
                    for &v in vals {
                        bin.extend_from_slice(&v.to_le_bytes());
                    }
                }
                DataValues::Int32(vals) => {
                    for &v in vals {
                        bin.extend_from_slice(&v.to_le_bytes());
                    }
                }
                DataValues::Int64(vals) => {
                    for &v in vals {
                        bin.extend_from_slice(&v.to_le_bytes());
                    }
                }
                DataValues::UInt8(vals) => {
                    bin.extend_from_slice(vals);
                }
            }

            writeln!(
                xmf,
                "      <Attribute Name=\"{}\" AttributeType=\"{}\" Center=\"{center}\">",
                arr.name, attr_type
            )
            .ok();
            writeln!(
                xmf,
                "        <DataItem DataType=\"Float\" Precision=\"8\" Dimensions=\"{dim_str}\" Format=\"Binary\" Seek=\"{arr_offset}\">"
            )
            .ok();
            writeln!(xmf, "          {binary_file_name}").ok();
            xmf.push_str("        </DataItem>\n");
            xmf.push_str("      </Attribute>\n");
        }

        xmf.push_str("    </Grid>\n");
        xmf.push_str("  </Domain>\n");
        xmf.push_str("</Xdmf>\n");

        Ok((xmf, bin))
    }
}
