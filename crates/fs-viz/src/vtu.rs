//! Deterministic VTU (VTK Unstructured Grid XML) export and independent verification.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.8`
//!
//! Provides canonical, bit-reproducible VTU XML emission for 3D unstructured meshes
//! and solution fields (conduction temperatures, heat flux vectors, material IDs,
//! pressure, velocity) compatible with ParaView, VisIt, and standard CAE post-processors.
//! Includes an independent reader/checker that validates mesh topology, array bounds,
//! finite values, and extrema.

use fs_blake3::{ContentHash, hash_domain};
use std::fmt::Write as _;

/// Standard VTK Cell Types (VTK 4.2+ / XML).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum CellType {
    Vertex = 1,
    PolyVertex = 2,
    Line = 3,
    PolyLine = 4,
    Triangle = 5,
    TriangleStrip = 6,
    Polygon = 7,
    Pixel = 8,
    Quad = 9,
    Tetra = 10,
    Voxel = 11,
    Hexahedron = 12,
    Wedge = 13,
    Pyramid = 14,
}

impl CellType {
    /// Return the integer type ID.
    #[must_use]
    pub const fn type_id(self) -> u8 {
        self as u8
    }

    /// Parse a cell type from its VTK integer ID.
    #[must_use]
    pub const fn from_type_id(id: u8) -> Option<Self> {
        match id {
            1 => Some(Self::Vertex),
            2 => Some(Self::PolyVertex),
            3 => Some(Self::Line),
            4 => Some(Self::PolyLine),
            5 => Some(Self::Triangle),
            6 => Some(Self::TriangleStrip),
            7 => Some(Self::Polygon),
            8 => Some(Self::Pixel),
            9 => Some(Self::Quad),
            10 => Some(Self::Tetra),
            11 => Some(Self::Voxel),
            12 => Some(Self::Hexahedron),
            13 => Some(Self::Wedge),
            14 => Some(Self::Pyramid),
            _ => None,
        }
    }

    /// Expected number of points for fixed-size cell types.
    #[must_use]
    pub const fn fixed_point_count(self) -> Option<usize> {
        match self {
            Self::Vertex => Some(1),
            Self::Line => Some(2),
            Self::Triangle => Some(3),
            Self::Pixel | Self::Quad | Self::Tetra => Some(4),
            Self::Pyramid => Some(5),
            Self::Wedge => Some(6),
            Self::Voxel | Self::Hexahedron => Some(8),
            Self::PolyVertex | Self::PolyLine | Self::Polygon | Self::TriangleStrip => None,
        }
    }
}

/// Association of a data array with mesh elements.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataAssociation {
    PointData,
    CellData,
}

/// Data types supported in VTK DataArray elements.
#[derive(Debug, Clone, PartialEq)]
pub enum DataValues {
    Float64(Vec<f64>),
    Float32(Vec<f32>),
    Int32(Vec<i32>),
    Int64(Vec<i64>),
    UInt8(Vec<u8>),
}

impl DataValues {
    /// Number of scalar elements in the array.
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            Self::Float64(v) => v.len(),
            Self::Float32(v) => v.len(),
            Self::Int32(v) => v.len(),
            Self::Int64(v) => v.len(),
            Self::UInt8(v) => v.len(),
        }
    }

    /// Whether the array is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// VTK data type name.
    #[must_use]
    pub const fn type_name(&self) -> &'static str {
        match self {
            Self::Float64(_) => "Float64",
            Self::Float32(_) => "Float32",
            Self::Int32(_) => "Int32",
            Self::Int64(_) => "Int64",
            Self::UInt8(_) => "UInt8",
        }
    }
}

/// A named data array associated with points or cells.
#[derive(Debug, Clone, PartialEq)]
pub struct DataArray {
    pub name: String,
    pub association: DataAssociation,
    pub components: usize,
    pub unit: Option<String>,
    pub values: DataValues,
}

impl DataArray {
    /// Create a new scalar float64 array.
    #[must_use]
    pub fn new_point_scalar(name: impl Into<String>, values: Vec<f64>) -> Self {
        Self {
            name: name.into(),
            association: DataAssociation::PointData,
            components: 1,
            unit: None,
            values: DataValues::Float64(values),
        }
    }

    /// Create a new 3-component vector float64 array on points.
    #[must_use]
    pub fn new_point_vector(name: impl Into<String>, vectors: &[[f64; 3]]) -> Self {
        let mut flat = Vec::with_capacity(vectors.len() * 3);
        for v in vectors {
            flat.extend_from_slice(v);
        }
        Self {
            name: name.into(),
            association: DataAssociation::PointData,
            components: 3,
            unit: None,
            values: DataValues::Float64(flat),
        }
    }

    /// Create a new cell scalar float64 array.
    #[must_use]
    pub fn new_cell_scalar(name: impl Into<String>, values: Vec<f64>) -> Self {
        Self {
            name: name.into(),
            association: DataAssociation::CellData,
            components: 1,
            unit: None,
            values: DataValues::Float64(values),
        }
    }

    /// Create a new cell int32 array (e.g. RegionId / MaterialIndex).
    #[must_use]
    pub fn new_cell_int32(name: impl Into<String>, values: Vec<i32>) -> Self {
        Self {
            name: name.into(),
            association: DataAssociation::CellData,
            components: 1,
            unit: None,
            values: DataValues::Int32(values),
        }
    }

    /// Attach a physical unit label.
    #[must_use]
    pub fn with_unit(mut self, unit: impl Into<String>) -> Self {
        self.unit = Some(unit.into());
        self
    }
}

/// An unstructured grid representing 3D geometry and solution fields.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct UnstructuredGrid {
    /// 3D Coordinates of points [x, y, z].
    pub points: Vec<[f64; 3]>,
    /// Connectivity list of vertex indices for all cells.
    pub cells_connectivity: Vec<i64>,
    /// Cumulative offsets into connectivity (one per cell).
    pub cells_offsets: Vec<i64>,
    /// Cell type IDs (one per cell).
    pub cells_types: Vec<u8>,
    /// Data arrays (point data and cell data).
    pub arrays: Vec<DataArray>,
}

/// Errors during VTU validation or serialization.
/// Errors arising from VTU grid validation or parsing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VtuError {
    /// The grid contains zero points.
    EmptyGrid,
    /// A coordinate component is NaN or infinite.
    NonFinitePointCoordinate {
        /// Point index.
        index: usize,
        /// Coordinate axis (0=x, 1=y, 2=z).
        component: usize,
    },
    /// A float field value is NaN or infinite.
    NonFiniteFieldValue {
        /// Name of the array.
        array: String,
        /// Element index.
        index: usize,
    },
    /// A cell references a point index out of bounds.
    InvalidCellIndex {
        /// Cell index.
        cell: usize,
        /// Point index referenced.
        point_index: i64,
        /// Maximum points in the mesh.
        max_points: usize,
    },
    /// Unrecognized VTK cell type ID.
    InvalidCellType {
        /// Cell index.
        cell: usize,
        /// VTK type ID byte.
        type_id: u8,
    },
    /// Number of points in cell connectivity does not match the cell type requirement.
    CellPointCountMismatch {
        /// Cell index.
        cell: usize,
        /// Expected points for the type.
        expected: usize,
        /// Points found in connectivity.
        found: usize,
    },
    /// Length of cell offsets array does not match cell count.
    OffsetsMismatch {
        /// Expected cell count.
        expected_cells: usize,
        /// Offsets array length.
        found_offsets: usize,
    },
    /// Length of cell types array does not match cell count.
    TypesMismatch {
        /// Expected cell count.
        expected_cells: usize,
        /// Types array length.
        found_types: usize,
    },
    /// Offsets are not strictly monotonically increasing.
    OffsetsNotMonotonic {
        /// Cell index.
        cell: usize,
        /// Previous offset value.
        prev: i64,
        /// Current offset value.
        current: i64,
    },
    /// Array length does not match expected point or cell count.
    ArrayLengthMismatch {
        /// Array name.
        array: String,
        /// Expected scalar items.
        expected: usize,
        /// Actual scalar items found.
        found: usize,
    },
    /// Duplicate array name within PointData or CellData.
    DuplicateArrayName {
        /// Name of the duplicate array.
        name: String,
    },
    /// Array declared with zero components per item.
    ZeroComponents {
        /// Array name.
        array: String,
    },
    /// XML syntax or formatting parse error.
    ParseError {
        /// Error details.
        detail: String,
    },
}

impl std::fmt::Display for VtuError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyGrid => write!(f, "unstructured grid has zero points"),
            Self::NonFinitePointCoordinate { index, component } => {
                write!(
                    f,
                    "point {index} component {component} is non-finite (NaN or Inf)"
                )
            }
            Self::NonFiniteFieldValue { array, index } => {
                write!(
                    f,
                    "array `{array}` at index {index} contains non-finite float"
                )
            }
            Self::InvalidCellIndex {
                cell,
                point_index,
                max_points,
            } => {
                write!(
                    f,
                    "cell {cell} references point index {point_index} >= {max_points}"
                )
            }
            Self::InvalidCellType { cell, type_id } => {
                write!(f, "cell {cell} has unknown VTK cell type ID {type_id}")
            }
            Self::CellPointCountMismatch {
                cell,
                expected,
                found,
            } => {
                write!(f, "cell {cell} expects {expected} points, found {found}")
            }
            Self::OffsetsMismatch {
                expected_cells,
                found_offsets,
            } => {
                write!(
                    f,
                    "expected {expected_cells} offsets, found {found_offsets}"
                )
            }
            Self::TypesMismatch {
                expected_cells,
                found_types,
            } => {
                write!(
                    f,
                    "expected {expected_cells} cell types, found {found_types}"
                )
            }
            Self::OffsetsNotMonotonic {
                cell,
                prev,
                current,
            } => {
                write!(f, "offset for cell {cell} ({current}) <= previous ({prev})")
            }
            Self::ArrayLengthMismatch {
                array,
                expected,
                found,
            } => {
                write!(
                    f,
                    "array `{array}` length mismatch: expected {expected}, found {found}"
                )
            }
            Self::DuplicateArrayName { name } => {
                write!(f, "duplicate data array name `{name}`")
            }
            Self::ZeroComponents { array } => {
                write!(f, "array `{array}` has zero components")
            }
            Self::ParseError { detail } => {
                write!(f, "VTU parse error: {detail}")
            }
        }
    }
}

impl std::error::Error for VtuError {}

impl UnstructuredGrid {
    /// Create a new empty grid.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a 3D point. Returns its index.
    pub fn add_point(&mut self, x: f64, y: f64, z: f64) -> usize {
        let index = self.points.len();
        self.points.push([x, y, z]);
        index
    }

    /// Add a cell of given type with its vertex indices.
    pub fn add_cell(&mut self, cell_type: CellType, indices: &[usize]) {
        for &idx in indices {
            self.cells_connectivity.push(idx as i64);
        }
        self.cells_offsets
            .push(self.cells_connectivity.len() as i64);
        self.cells_types.push(cell_type.type_id());
    }

    /// Add a 4-node tetrahedral cell.
    pub fn add_tetra(&mut self, i0: usize, i1: usize, i2: usize, i3: usize) {
        self.add_cell(CellType::Tetra, &[i0, i1, i2, i3]);
    }

    /// Add a 3-node triangle cell.
    pub fn add_triangle(&mut self, i0: usize, i1: usize, i2: usize) {
        self.add_cell(CellType::Triangle, &[i0, i1, i2]);
    }

    /// Add a data array.
    pub fn add_array(&mut self, array: DataArray) {
        self.arrays.push(array);
    }

    /// Number of cells.
    #[must_use]
    pub fn num_cells(&self) -> usize {
        self.cells_types.len()
    }

    /// Number of points.
    #[must_use]
    pub fn num_points(&self) -> usize {
        self.points.len()
    }

    /// Validate internal structural consistency and bounds.
    pub fn validate(&self) -> Result<(), VtuError> {
        if self.points.is_empty() {
            return Err(VtuError::EmptyGrid);
        }

        // Validate point coordinates
        for (i, p) in self.points.iter().enumerate() {
            for (c, &val) in p.iter().enumerate() {
                if !val.is_finite() {
                    return Err(VtuError::NonFinitePointCoordinate {
                        index: i,
                        component: c,
                    });
                }
            }
        }

        let num_cells = self.num_cells();
        if self.cells_offsets.len() != num_cells {
            return Err(VtuError::OffsetsMismatch {
                expected_cells: num_cells,
                found_offsets: self.cells_offsets.len(),
            });
        }

        // Validate offsets monotonicity and cell connectivity
        let mut prev_offset = 0i64;
        for (c_idx, (&offset, &type_id)) in
            self.cells_offsets.iter().zip(&self.cells_types).enumerate()
        {
            if offset <= prev_offset {
                return Err(VtuError::OffsetsNotMonotonic {
                    cell: c_idx,
                    prev: prev_offset,
                    current: offset,
                });
            }
            let cell_type = CellType::from_type_id(type_id).ok_or(VtuError::InvalidCellType {
                cell: c_idx,
                type_id,
            })?;
            let pts_count = (offset - prev_offset) as usize;
            if let Some(expected) = cell_type.fixed_point_count()
                && pts_count != expected
            {
                return Err(VtuError::CellPointCountMismatch {
                    cell: c_idx,
                    expected,
                    found: pts_count,
                });
            }
            for p_idx in (prev_offset as usize)..(offset as usize) {
                let p = self.cells_connectivity[p_idx];
                if p < 0 || (p as usize) >= self.points.len() {
                    return Err(VtuError::InvalidCellIndex {
                        cell: c_idx,
                        point_index: p,
                        max_points: self.points.len(),
                    });
                }
            }
            prev_offset = offset;
        }

        // Validate data arrays
        let mut seen_names = std::collections::HashSet::new();
        for arr in &self.arrays {
            if !seen_names.insert(&arr.name) {
                return Err(VtuError::DuplicateArrayName {
                    name: arr.name.clone(),
                });
            }
            if arr.components == 0 {
                return Err(VtuError::ZeroComponents {
                    array: arr.name.clone(),
                });
            }
            let expected_items = match arr.association {
                DataAssociation::PointData => self.points.len() * arr.components,
                DataAssociation::CellData => num_cells * arr.components,
            };
            if arr.values.len() != expected_items {
                return Err(VtuError::ArrayLengthMismatch {
                    array: arr.name.clone(),
                    expected: expected_items,
                    found: arr.values.len(),
                });
            }
            // Check finite
            match &arr.values {
                DataValues::Float64(vals) => {
                    for (i, &v) in vals.iter().enumerate() {
                        if !v.is_finite() {
                            return Err(VtuError::NonFiniteFieldValue {
                                array: arr.name.clone(),
                                index: i,
                            });
                        }
                    }
                }
                DataValues::Float32(vals) => {
                    for (i, &v) in vals.iter().enumerate() {
                        if !v.is_finite() {
                            return Err(VtuError::NonFiniteFieldValue {
                                array: arr.name.clone(),
                                index: i,
                            });
                        }
                    }
                }
                _ => {}
            }
        }

        Ok(())
    }
}

/// Writer for deterministic VTU XML representation.
pub struct VtuWriter;

impl VtuWriter {
    /// Export an unstructured grid to a bit-reproducible VTU XML string.
    pub fn write_ascii(grid: &UnstructuredGrid) -> Result<String, VtuError> {
        grid.validate()?;

        let mut out = String::with_capacity(16 * 1024);
        out.push_str("<?xml version=\"1.0\"?>\n");
        out.push_str("<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
        out.push_str("  <UnstructuredGrid>\n");
        let _ = writeln!(
            out,
            "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
            grid.points.len(),
            grid.num_cells()
        );

        // 1. Points
        out.push_str("      <Points>\n");
        out.push_str("        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n          ");
        for (i, p) in grid.points.iter().enumerate() {
            if i > 0 {
                out.push(' ');
            }
            let _ = write!(out, "{:.17e} {:.17e} {:.17e}", p[0], p[1], p[2]);
        }
        out.push_str("\n        </DataArray>\n");
        out.push_str("      </Points>\n");

        // 2. Cells
        out.push_str("      <Cells>\n");
        // Connectivity
        out.push_str(
            "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n          ",
        );
        for (i, &c) in grid.cells_connectivity.iter().enumerate() {
            if i > 0 {
                out.push(' ');
            }
            let _ = write!(out, "{c}");
        }
        out.push_str("\n        </DataArray>\n");

        // Offsets
        out.push_str(
            "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n          ",
        );
        for (i, &off) in grid.cells_offsets.iter().enumerate() {
            if i > 0 {
                out.push(' ');
            }
            let _ = write!(out, "{off}");
        }
        out.push_str("\n        </DataArray>\n");

        // Types
        out.push_str(
            "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ",
        );
        for (i, &t) in grid.cells_types.iter().enumerate() {
            if i > 0 {
                out.push(' ');
            }
            let _ = write!(out, "{t}");
        }
        out.push_str("\n        </DataArray>\n");
        out.push_str("      </Cells>\n");

        // 3. PointData
        let point_arrays: Vec<&DataArray> = grid
            .arrays
            .iter()
            .filter(|a| a.association == DataAssociation::PointData)
            .collect();
        if !point_arrays.is_empty() {
            out.push_str("      <PointData>\n");
            for arr in point_arrays {
                Self::write_data_array(&mut out, arr);
            }
            out.push_str("      </PointData>\n");
        }

        // 4. CellData
        let cell_arrays: Vec<&DataArray> = grid
            .arrays
            .iter()
            .filter(|a| a.association == DataAssociation::CellData)
            .collect();
        if !cell_arrays.is_empty() {
            out.push_str("      <CellData>\n");
            for arr in cell_arrays {
                Self::write_data_array(&mut out, arr);
            }
            out.push_str("      </CellData>\n");
        }

        out.push_str("    </Piece>\n");
        out.push_str("  </UnstructuredGrid>\n");
        out.push_str("</VTKFile>\n");

        Ok(out)
    }

    fn write_data_array(out: &mut String, arr: &DataArray) {
        let type_name = arr.values.type_name();
        let unit_attr = arr
            .unit
            .as_ref()
            .map_or(String::new(), |u| format!(" Unit=\"{u}\""));
        let _ = write!(
            out,
            "        <DataArray type=\"{}\" Name=\"{}\" NumberOfComponents=\"{}\" format=\"ascii\"{}>\n          ",
            type_name, arr.name, arr.components, unit_attr
        );
        match &arr.values {
            DataValues::Float64(vals) => {
                for (i, &v) in vals.iter().enumerate() {
                    if i > 0 {
                        out.push(' ');
                    }
                    let _ = write!(out, "{v:.17e}");
                }
            }
            DataValues::Float32(vals) => {
                for (i, &v) in vals.iter().enumerate() {
                    if i > 0 {
                        out.push(' ');
                    }
                    let _ = write!(out, "{v:.9e}");
                }
            }
            DataValues::Int32(vals) => {
                for (i, &v) in vals.iter().enumerate() {
                    if i > 0 {
                        out.push(' ');
                    }
                    let _ = write!(out, "{v}");
                }
            }
            DataValues::Int64(vals) => {
                for (i, &v) in vals.iter().enumerate() {
                    if i > 0 {
                        out.push(' ');
                    }
                    let _ = write!(out, "{v}");
                }
            }
            DataValues::UInt8(vals) => {
                for (i, &v) in vals.iter().enumerate() {
                    if i > 0 {
                        out.push(' ');
                    }
                    let _ = write!(out, "{v}");
                }
            }
        }
        out.push_str("\n        </DataArray>\n");
    }

    /// Compute the BLAKE3 digest of the canonical VTU serialization.
    pub fn content_hash(grid: &UnstructuredGrid) -> Result<ContentHash, VtuError> {
        let xml = Self::write_ascii(grid)?;
        Ok(hash_domain("org.frankensim.vtu.ascii.v1", xml.as_bytes()))
    }
}

/// An independent offline VTU reader and checker for round-trip verification.
pub struct VtuChecker;

/// Summary report produced by the VTU checker.
#[derive(Debug, Clone, PartialEq)]
pub struct VtuReport {
    /// Total number of mesh vertices/points.
    pub num_points: usize,
    /// Total number of mesh elements/cells.
    pub num_cells: usize,
    /// Number of PointData and CellData arrays attached.
    pub array_count: usize,
    /// Axis-aligned bounding box [[x_min, x_max], [y_min, y_max], [z_min, z_max]].
    pub point_bounds: [[f64; 2]; 3],
    /// Min/max scalar extrema for each attached field array.
    pub array_extrema: Vec<(String, [f64; 2])>,
    /// Cryptographic content digest of the verified mesh.
    pub content_hash: ContentHash,
}

impl VtuChecker {
    /// Parse and independently verify a VTU XML string.
    pub fn check(xml: &str) -> Result<VtuReport, VtuError> {
        let grid = Self::parse_ascii(xml)?;
        grid.validate()?;

        let mut x_min = f64::INFINITY;
        let mut x_max = f64::NEG_INFINITY;
        let mut y_min = f64::INFINITY;
        let mut y_max = f64::NEG_INFINITY;
        let mut z_min = f64::INFINITY;
        let mut z_max = f64::NEG_INFINITY;

        for p in &grid.points {
            x_min = x_min.min(p[0]);
            x_max = x_max.max(p[0]);
            y_min = y_min.min(p[1]);
            y_max = y_max.max(p[1]);
            z_min = z_min.min(p[2]);
            z_max = z_max.max(p[2]);
        }

        let mut array_extrema = Vec::new();
        for arr in &grid.arrays {
            let mut v_min = f64::INFINITY;
            let mut v_max = f64::NEG_INFINITY;
            match &arr.values {
                DataValues::Float64(vals) => {
                    for &v in vals {
                        v_min = v_min.min(v);
                        v_max = v_max.max(v);
                    }
                }
                DataValues::Float32(vals) => {
                    for &v in vals {
                        v_min = v_min.min(f64::from(v));
                        v_max = v_max.max(f64::from(v));
                    }
                }
                DataValues::Int32(vals) => {
                    for &v in vals {
                        v_min = v_min.min(f64::from(v));
                        v_max = v_max.max(f64::from(v));
                    }
                }
                DataValues::Int64(vals) => {
                    for &v in vals {
                        v_min = v_min.min(v as f64);
                        v_max = v_max.max(v as f64);
                    }
                }
                DataValues::UInt8(vals) => {
                    for &v in vals {
                        v_min = v_min.min(f64::from(v));
                        v_max = v_max.max(f64::from(v));
                    }
                }
            }
            array_extrema.push((arr.name.clone(), [v_min, v_max]));
        }

        let content_hash = hash_domain("org.frankensim.vtu.ascii.v1", xml.as_bytes());

        Ok(VtuReport {
            num_points: grid.num_points(),
            num_cells: grid.num_cells(),
            array_count: grid.arrays.len(),
            point_bounds: [[x_min, x_max], [y_min, y_max], [z_min, z_max]],
            array_extrema,
            content_hash,
        })
    }

    /// Independent parser for ASCII VTU format.
    pub fn parse_ascii(xml: &str) -> Result<UnstructuredGrid, VtuError> {
        let mut grid = UnstructuredGrid::new();

        // 1. Extract Points
        if let Some(points_content) = extract_tag_content(xml, "<Points>", "</Points>")
            && let Some(data) = extract_data_array_content(&points_content)
        {
            let coords: Vec<f64> = data
                .split_whitespace()
                .map(|s| {
                    s.parse::<f64>().map_err(|e| VtuError::ParseError {
                        detail: format!("invalid point float `{s}`: {e}"),
                    })
                })
                .collect::<Result<_, _>>()?;
            if !coords.len().is_multiple_of(3) {
                return Err(VtuError::ParseError {
                    detail: format!(
                        "points coordinate count {} is not multiple of 3",
                        coords.len()
                    ),
                });
            }
            for chunk in coords.as_chunks::<3>().0 {
                grid.points.push([chunk[0], chunk[1], chunk[2]]);
            }
        }

        // 2. Extract Cells
        if let Some(cells_content) = extract_tag_content(xml, "<Cells>", "</Cells>") {
            // Connectivity
            if let Some(conn_tag) = find_data_array_by_name(&cells_content, "connectivity") {
                let indices: Vec<i64> = conn_tag
                    .split_whitespace()
                    .map(|s| {
                        s.parse::<i64>().map_err(|e| VtuError::ParseError {
                            detail: format!("invalid connectivity int `{s}`: {e}"),
                        })
                    })
                    .collect::<Result<_, _>>()?;
                grid.cells_connectivity = indices;
            }

            // Offsets
            if let Some(offsets_tag) = find_data_array_by_name(&cells_content, "offsets") {
                let offsets: Vec<i64> = offsets_tag
                    .split_whitespace()
                    .map(|s| {
                        s.parse::<i64>().map_err(|e| VtuError::ParseError {
                            detail: format!("invalid offset int `{s}`: {e}"),
                        })
                    })
                    .collect::<Result<_, _>>()?;
                grid.cells_offsets = offsets;
            }

            // Types
            if let Some(types_tag) = find_data_array_by_name(&cells_content, "types") {
                let types: Vec<u8> = types_tag
                    .split_whitespace()
                    .map(|s| {
                        s.parse::<u8>().map_err(|e| VtuError::ParseError {
                            detail: format!("invalid type u8 `{s}`: {e}"),
                        })
                    })
                    .collect::<Result<_, _>>()?;
                grid.cells_types = types;
            }
        }

        // 3. Extract PointData
        if let Some(pd_content) = extract_tag_content(xml, "<PointData>", "</PointData>") {
            parse_arrays_in_section(&pd_content, DataAssociation::PointData, &mut grid)?;
        }

        // 4. Extract CellData
        if let Some(cd_content) = extract_tag_content(xml, "<CellData>", "</CellData>") {
            parse_arrays_in_section(&cd_content, DataAssociation::CellData, &mut grid)?;
        }

        Ok(grid)
    }
}

fn extract_tag_content(source: &str, open_tag: &str, close_tag: &str) -> Option<String> {
    let start = source.find(open_tag)? + open_tag.len();
    let end = source[start..].find(close_tag)?;
    Some(source[start..start + end].to_string())
}

fn extract_data_array_content(source: &str) -> Option<String> {
    let open_end = source.find("<DataArray")?;
    let content_start = source[open_end..].find('>')? + open_end + 1;
    let content_end = source[content_start..].find("</DataArray>")? + content_start;
    Some(source[content_start..content_end].trim().to_string())
}

fn find_data_array_by_name(section: &str, name: &str) -> Option<String> {
    let target = format!("Name=\"{name}\"");
    let target_pos = section.find(&target)?;
    let array_start = section[..target_pos].rfind("<DataArray")?;
    let content_start = section[array_start..].find('>')? + array_start + 1;
    let content_end = section[content_start..].find("</DataArray>")? + content_start;
    Some(section[content_start..content_end].trim().to_string())
}

fn parse_arrays_in_section(
    section: &str,
    assoc: DataAssociation,
    grid: &mut UnstructuredGrid,
) -> Result<(), VtuError> {
    let mut cursor = section;
    while let Some(start_pos) = cursor.find("<DataArray") {
        let tag_header_end = cursor[start_pos..]
            .find('>')
            .ok_or_else(|| VtuError::ParseError {
                detail: "unterminated DataArray tag".to_string(),
            })?
            + start_pos;
        let tag_header = &cursor[start_pos..tag_header_end];

        let content_end = cursor[tag_header_end..]
            .find("</DataArray>")
            .ok_or_else(|| VtuError::ParseError {
                detail: "missing </DataArray> closing tag".to_string(),
            })?
            + tag_header_end;
        let content = cursor[tag_header_end + 1..content_end].trim();

        // Extract attributes
        let name = extract_attr(tag_header, "Name").unwrap_or_else(|| "unnamed".to_string());
        let type_name = extract_attr(tag_header, "type").unwrap_or_else(|| "Float64".to_string());
        let components = extract_attr(tag_header, "NumberOfComponents")
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(1);
        let unit = extract_attr(tag_header, "Unit");

        let values = match type_name.as_str() {
            "Float64" => {
                let vals: Vec<f64> = content
                    .split_whitespace()
                    .map(|s| {
                        s.parse::<f64>().map_err(|e| VtuError::ParseError {
                            detail: format!("invalid float in array `{name}`: {e}"),
                        })
                    })
                    .collect::<Result<_, _>>()?;
                DataValues::Float64(vals)
            }
            "Int32" => {
                let vals: Vec<i32> = content
                    .split_whitespace()
                    .map(|s| {
                        s.parse::<i32>().map_err(|e| VtuError::ParseError {
                            detail: format!("invalid int in array `{name}`: {e}"),
                        })
                    })
                    .collect::<Result<Vec<i32>, VtuError>>()?;
                DataValues::Int32(vals)
            }
            _ => {
                let vals: Vec<f64> = content
                    .split_whitespace()
                    .map(|s| {
                        s.parse::<f64>().map_err(|e| VtuError::ParseError {
                            detail: format!("invalid float in array `{name}`: {e}"),
                        })
                    })
                    .collect::<Result<_, _>>()?;
                DataValues::Float64(vals)
            }
        };

        grid.add_array(DataArray {
            name,
            association: assoc,
            components,
            unit,
            values,
        });

        cursor = &cursor[content_end + "</DataArray>".len()..];
    }
    Ok(())
}

fn extract_attr(header: &str, attr: &str) -> Option<String> {
    let pattern = format!("{attr}=\"");
    let start = header.find(&pattern)? + pattern.len();
    let end = header[start..].find('"')? + start;
    Some(header[start..end].to_string())
}
