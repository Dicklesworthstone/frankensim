//! Deterministic XDMF 3.0 / 2.0 descriptor emission for heavy dataset visualization.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.8`
//!
//! Generates clean, standard XDMF XML files with raw binary heavy-data containers
//! (strictly pure Rust, zero C / HDF5 dependency) for large-scale solution fields.

use super::vtu::{
    CellType, DataAssociation, DataValues, UnstructuredGrid, VtuError, push_xml_attribute_escaped,
    push_xml_text_escaped, validate_xml_chars,
};
use std::fmt::Write as _;

struct TopologyData {
    topology_type: &'static str,
    nodes_per_element: Option<usize>,
    dimensions: String,
    values: Vec<u64>,
}

fn append_cell_indices(values: &mut Vec<u64>, cell_type: CellType, indices: &[i64]) {
    let order: &[usize] = match cell_type {
        // VTK_PIXEL and VTK_VOXEL use axis-fastest ordering, unlike XDMF's
        // perimeter/facet ordering for quadrilaterals and hexahedra.
        CellType::Pixel => &[0, 1, 3, 2],
        CellType::Voxel => &[0, 1, 3, 2, 4, 5, 7, 6],
        _ => {
            values.extend(indices.iter().map(|&index| index as u64));
            return;
        }
    };
    values.extend(order.iter().map(|&index| indices[index] as u64));
}

fn homogeneous_topology(cell_type: CellType) -> Option<(&'static str, bool)> {
    match cell_type {
        CellType::Vertex | CellType::PolyVertex => Some(("Polyvertex", true)),
        CellType::Line | CellType::PolyLine => Some(("Polyline", true)),
        CellType::Triangle => Some(("Triangle", false)),
        CellType::Polygon => Some(("Polygon", true)),
        CellType::Pixel | CellType::Quad => Some(("Quadrilateral", false)),
        CellType::Tetra => Some(("Tetrahedron", false)),
        CellType::Voxel | CellType::Hexahedron => Some(("Hexahedron", false)),
        CellType::Wedge => Some(("Wedge", false)),
        CellType::Pyramid => Some(("Pyramid", false)),
        CellType::TriangleStrip => None,
    }
}

fn mixed_cell_code(cell_type: CellType) -> Option<(u64, bool)> {
    match cell_type {
        CellType::Vertex | CellType::PolyVertex => Some((1, true)),
        CellType::Line | CellType::PolyLine => Some((2, true)),
        CellType::Polygon => Some((3, true)),
        CellType::Triangle => Some((4, false)),
        CellType::Pixel | CellType::Quad => Some((5, false)),
        CellType::Tetra => Some((6, false)),
        CellType::Pyramid => Some((7, false)),
        CellType::Wedge => Some((8, false)),
        CellType::Voxel | CellType::Hexahedron => Some((9, false)),
        CellType::TriangleStrip => None,
    }
}

fn topology_data(grid: &UnstructuredGrid) -> Result<TopologyData, VtuError> {
    let first_type_id = *grid
        .cells_types
        .first()
        .ok_or(VtuError::EmptyXdmfTopology)?;
    let first_type = CellType::from_type_id(first_type_id).ok_or(VtuError::InvalidCellType {
        cell: 0,
        type_id: first_type_id,
    })?;

    let mut cells = Vec::with_capacity(grid.num_cells());
    let mut start = 0usize;
    for (cell, (&end, &type_id)) in grid.cells_offsets.iter().zip(&grid.cells_types).enumerate() {
        let end = end as usize;
        let cell_type =
            CellType::from_type_id(type_id).ok_or(VtuError::InvalidCellType { cell, type_id })?;
        cells.push((cell_type, &grid.cells_connectivity[start..end]));
        start = end;
    }

    let homogeneous = cells.iter().all(|(cell_type, _)| *cell_type == first_type);
    if homogeneous
        && let Some((topology_type, needs_nodes_per_element)) = homogeneous_topology(first_type)
    {
        let nodes_per_element = cells[0].1.len();
        if cells
            .iter()
            .all(|(_, indices)| indices.len() == nodes_per_element)
        {
            let mut values = Vec::with_capacity(grid.cells_connectivity.len());
            for (_, indices) in &cells {
                append_cell_indices(&mut values, first_type, indices);
            }
            return Ok(TopologyData {
                topology_type,
                nodes_per_element: needs_nodes_per_element.then_some(nodes_per_element),
                dimensions: format!("{} {nodes_per_element}", grid.num_cells()),
                values,
            });
        }
    }

    let mut values = Vec::with_capacity(grid.cells_connectivity.len() + grid.num_cells() * 2);
    for (cell, (cell_type, indices)) in cells.iter().enumerate() {
        let (code, needs_node_count) =
            mixed_cell_code(*cell_type).ok_or(VtuError::UnsupportedXdmfCellType {
                cell,
                type_id: cell_type.type_id(),
            })?;
        values.push(code);
        if needs_node_count {
            values.push(indices.len() as u64);
        }
        append_cell_indices(&mut values, *cell_type, indices);
    }
    Ok(TopologyData {
        topology_type: "Mixed",
        nodes_per_element: None,
        dimensions: values.len().to_string(),
        values,
    })
}

/// Writer for deterministic XDMF representations.
pub struct XdmfWriter;

impl XdmfWriter {
    /// Generate the XDMF XML text and the raw binary companion buffer.
    pub fn write_xdmf_with_binary(
        grid: &UnstructuredGrid,
        binary_file_name: &str,
    ) -> Result<(String, Vec<u8>), VtuError> {
        grid.validate()?;
        validate_xml_chars("XDMF binary file name", binary_file_name)?;

        let mut bin = Vec::with_capacity(16 * 1024);
        let mut xmf = String::with_capacity(4096);

        xmf.push_str("<?xml version=\"1.0\" ?>\n");
        xmf.push_str("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
        xmf.push_str("<Xdmf Version=\"3.0\">\n");
        xmf.push_str("  <Domain>\n");
        xmf.push_str("    <Grid Name=\"FrankenSimGrid\" GridType=\"Uniform\">\n");

        // 1. Topology
        let num_cells = grid.num_cells();
        let topology = topology_data(grid)?;

        let conn_offset = bin.len();
        for &value in &topology.values {
            bin.extend_from_slice(&value.to_le_bytes());
        }

        write!(
            xmf,
            "      <Topology TopologyType=\"{}\" NumberOfElements=\"{num_cells}\"",
            topology.topology_type
        )
        .ok();
        if let Some(nodes_per_element) = topology.nodes_per_element {
            write!(xmf, " NodesPerElement=\"{nodes_per_element}\"").ok();
        }
        xmf.push_str(">\n");
        writeln!(
            xmf,
            "        <DataItem NumberType=\"UInt\" Precision=\"8\" Dimensions=\"{}\" Format=\"Binary\" Endian=\"Little\" Seek=\"{conn_offset}\">",
            topology.dimensions
        )
        .ok();
        xmf.push_str("          ");
        push_xml_text_escaped(&mut xmf, binary_file_name);
        xmf.push('\n');
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
            "        <DataItem NumberType=\"Float\" Precision=\"8\" Dimensions=\"{num_points} 3\" Format=\"Binary\" Endian=\"Little\" Seek=\"{geom_offset}\">"
        )
        .ok();
        xmf.push_str("          ");
        push_xml_text_escaped(&mut xmf, binary_file_name);
        xmf.push('\n');
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

            let (number_type, precision) = match &arr.values {
                DataValues::Float64(vals) => {
                    for &v in vals {
                        bin.extend_from_slice(&v.to_le_bytes());
                    }
                    ("Float", 8)
                }
                DataValues::Float32(vals) => {
                    for &v in vals {
                        bin.extend_from_slice(&v.to_le_bytes());
                    }
                    ("Float", 4)
                }
                DataValues::Int32(vals) => {
                    for &v in vals {
                        bin.extend_from_slice(&v.to_le_bytes());
                    }
                    ("Int", 4)
                }
                DataValues::Int64(vals) => {
                    for &v in vals {
                        bin.extend_from_slice(&v.to_le_bytes());
                    }
                    ("Int", 8)
                }
                DataValues::UInt8(vals) => {
                    bin.extend_from_slice(vals);
                    ("UInt", 1)
                }
            };

            xmf.push_str("      <Attribute Name=\"");
            push_xml_attribute_escaped(&mut xmf, &arr.name);
            writeln!(xmf, "\" AttributeType=\"{attr_type}\" Center=\"{center}\">").ok();
            writeln!(
                xmf,
                "        <DataItem NumberType=\"{number_type}\" Precision=\"{precision}\" Dimensions=\"{dim_str}\" Format=\"Binary\" Endian=\"Little\" Seek=\"{arr_offset}\">"
            )
            .ok();
            xmf.push_str("          ");
            push_xml_text_escaped(&mut xmf, binary_file_name);
            xmf.push('\n');
            xmf.push_str("        </DataItem>\n");
            xmf.push_str("      </Attribute>\n");
        }

        xmf.push_str("    </Grid>\n");
        xmf.push_str("  </Domain>\n");
        xmf.push_str("</Xdmf>\n");

        Ok((xmf, bin))
    }
}
