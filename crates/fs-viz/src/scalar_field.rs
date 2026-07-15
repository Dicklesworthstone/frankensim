//! Versioned, bounded scalar-field artifacts shared by L5 visualization paths.

use crate::{Grid3, Grid3Error};

/// Ledger artifact kind for [`ScalarField3`] schema v1 bytes.
pub const SCALAR_FIELD3_ARTIFACT_KIND: &str = "frankensim.scalar-field3";
/// Current scalar-field byte-schema version.
pub const SCALAR_FIELD3_SCHEMA_VERSION: u32 = 1;

const MAGIC: [u8; 8] = *b"FSVIZF3\0";
const F64_LE_ENCODING: u32 = 1;
const FIXED_HEADER_BYTES: usize = 116;
const MAX_SEMANTIC_BYTES: usize = 64;

/// Meaning of grid coordinates in a scalar-field artifact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScalarLayout3 {
    /// `origin` is the first sample node and `spacing` separates nodes.
    NodeCentered,
    /// `origin` is the minimum cell corner and samples lie at cell centers.
    CellCentered,
}

impl ScalarLayout3 {
    const fn code(self) -> u32 {
        match self {
            Self::NodeCentered => 1,
            Self::CellCentered => 2,
        }
    }
}

/// Unit and quantity semantics carried inside scalar-field bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ScalarFieldSemantics {
    /// Physical or mathematical quantity, such as `density` or `temperature`.
    pub quantity: String,
    /// Unit for `origin` and `spacing`, such as `m`.
    pub coordinate_unit: String,
    /// Unit for scalar samples, such as `kg/m^3`.
    pub value_unit: String,
}

/// A validated owned 3-D scalar field ready for artifact encoding.
#[derive(Debug, Clone, PartialEq)]
pub struct ScalarField3 {
    layout: ScalarLayout3,
    dimensions: [usize; 3],
    origin: [f64; 3],
    spacing: [f64; 3],
    semantics: ScalarFieldSemantics,
    values: Vec<f64>,
}

/// Fail-closed scalar-field construction, codec, and conversion errors.
#[derive(Debug, Clone, PartialEq)]
pub enum ScalarField3Error {
    /// Every dimension must contain at least one sample.
    InvalidDimensions {
        /// Rejected dimensions.
        dimensions: [usize; 3],
    },
    /// The dimension product does not fit in `usize`.
    SampleCountOverflow {
        /// Rejected dimensions.
        dimensions: [usize; 3],
    },
    /// The field exceeds its explicit sample budget.
    SampleBudgetExceeded {
        /// Required samples.
        required: usize,
        /// Caller-provided limit.
        limit: usize,
    },
    /// Origin, spacing, or the implied world bound is invalid.
    InvalidGeometry {
        /// Cartesian axis, `0..3`.
        axis: usize,
    },
    /// Supplied storage does not match the dimension product.
    ValueCountMismatch {
        /// Required values.
        expected: usize,
        /// Supplied values.
        actual: usize,
    },
    /// A scalar sample is not finite.
    NonFiniteValue {
        /// Linear x-fastest sample index.
        index: usize,
    },
    /// A quantity or unit label is empty, too long, or contains controls.
    InvalidSemantic {
        /// Rejected semantic field.
        field: &'static str,
    },
    /// Encoded-size arithmetic overflowed.
    EncodedSizeOverflow,
    /// Encoded bytes exceed an explicit read/write budget.
    ByteBudgetExceeded {
        /// Required bytes.
        required: usize,
        /// Caller-provided limit.
        limit: usize,
    },
    /// A required allocation could not be reserved in full.
    AllocationFailed {
        /// Allocation purpose.
        context: &'static str,
        /// Requested elements or bytes, as named by `context`.
        requested: usize,
    },
    /// Artifact bytes are truncated, inconsistent, or have trailing data.
    Malformed {
        /// Stable parser diagnosis.
        what: &'static str,
    },
    /// The byte schema is newer or otherwise unsupported.
    UnsupportedSchema {
        /// Encoded schema version.
        found: u32,
    },
    /// The scalar encoding is unsupported.
    UnsupportedEncoding {
        /// Encoded representation tag.
        found: u32,
    },
    /// The sample-layout tag is unsupported.
    UnsupportedLayout {
        /// Encoded layout tag.
        found: u32,
    },
    /// An encoded dimension cannot be represented by this target's `usize`.
    DimensionOutOfRange {
        /// Cartesian axis, `0..3`.
        axis: usize,
    },
    /// Only node-centered artifacts can become an isosurface [`Grid3`].
    NotNodeCentered,
    /// Node-centered conversion failed the stricter isosurface-grid contract.
    Grid(Grid3Error),
}

impl core::fmt::Display for ScalarField3Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidDimensions { dimensions } => {
                write!(
                    f,
                    "scalar-field dimensions must be positive (got {dimensions:?})"
                )
            }
            Self::SampleCountOverflow { dimensions } => {
                write!(f, "scalar-field sample count overflows for {dimensions:?}")
            }
            Self::SampleBudgetExceeded { required, limit } => write!(
                f,
                "scalar field requires {required} samples, exceeding limit {limit}"
            ),
            Self::InvalidGeometry { axis } => {
                write!(f, "scalar-field geometry is invalid on axis {axis}")
            }
            Self::ValueCountMismatch { expected, actual } => write!(
                f,
                "scalar field requires {expected} values but received {actual}"
            ),
            Self::NonFiniteValue { index } => {
                write!(f, "scalar-field value {index} is non-finite")
            }
            Self::InvalidSemantic { field } => {
                write!(f, "scalar-field semantic {field} is invalid")
            }
            Self::EncodedSizeOverflow => write!(f, "scalar-field encoded size overflows"),
            Self::ByteBudgetExceeded { required, limit } => write!(
                f,
                "scalar field requires {required} encoded bytes, exceeding limit {limit}"
            ),
            Self::AllocationFailed { context, requested } => {
                write!(f, "could not reserve {requested} {context}")
            }
            Self::Malformed { what } => write!(f, "malformed scalar-field artifact: {what}"),
            Self::UnsupportedSchema { found } => {
                write!(f, "unsupported scalar-field schema {found}")
            }
            Self::UnsupportedEncoding { found } => {
                write!(f, "unsupported scalar-field encoding {found}")
            }
            Self::UnsupportedLayout { found } => {
                write!(f, "unsupported scalar-field layout {found}")
            }
            Self::DimensionOutOfRange { axis } => {
                write!(f, "scalar-field dimension {axis} does not fit this target")
            }
            Self::NotNodeCentered => {
                write!(
                    f,
                    "cell-centered scalar fields cannot become isosurface node grids"
                )
            }
            Self::Grid(error) => write!(f, "scalar-field node-grid conversion failed: {error}"),
        }
    }
}

impl std::error::Error for ScalarField3Error {}

impl ScalarField3 {
    /// Validate and adopt x-fastest scalar samples.
    ///
    /// # Errors
    /// [`ScalarField3Error`] on dimensions, budget, geometry, semantics,
    /// length, or sample finiteness.
    pub fn new(
        layout: ScalarLayout3,
        dimensions: [usize; 3],
        origin: [f64; 3],
        spacing: [f64; 3],
        semantics: ScalarFieldSemantics,
        sample_limit: usize,
        values: Vec<f64>,
    ) -> Result<Self, ScalarField3Error> {
        let sample_count = validate_layout(layout, dimensions, origin, spacing, sample_limit)?;
        validate_semantics(&semantics)?;
        if values.len() != sample_count {
            return Err(ScalarField3Error::ValueCountMismatch {
                expected: sample_count,
                actual: values.len(),
            });
        }
        if let Some(index) = values.iter().position(|value| !value.is_finite()) {
            return Err(ScalarField3Error::NonFiniteValue { index });
        }
        Ok(Self {
            layout,
            dimensions,
            origin,
            spacing,
            semantics,
            values,
        })
    }

    /// Sample layout.
    #[must_use]
    pub const fn layout(&self) -> ScalarLayout3 {
        self.layout
    }

    /// `[nx, ny, nz]` sample dimensions.
    #[must_use]
    pub const fn dimensions(&self) -> [usize; 3] {
        self.dimensions
    }

    /// First node for node-centered fields, or lower cell corner otherwise.
    #[must_use]
    pub const fn origin(&self) -> [f64; 3] {
        self.origin
    }

    /// Per-axis node or cell spacing.
    #[must_use]
    pub const fn spacing(&self) -> [f64; 3] {
        self.spacing
    }

    /// Quantity and unit semantics stored in the artifact.
    #[must_use]
    pub const fn semantics(&self) -> &ScalarFieldSemantics {
        &self.semantics
    }

    /// Read-only x-fastest scalar values.
    #[must_use]
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// World-space coverage bounds.
    ///
    /// Node-centered bounds end at the final node; cell-centered bounds end
    /// at the upper corner of the final cell.
    #[must_use]
    pub fn world_bounds(&self) -> [[f64; 3]; 2] {
        [
            self.origin,
            upper_bound(self.layout, self.dimensions, self.origin, self.spacing),
        ]
    }

    /// Encode schema-v1 bytes after admitting the complete byte budget.
    pub fn encode(&self, byte_limit: usize) -> Result<Vec<u8>, ScalarField3Error> {
        let required = self.encoded_len()?;
        if required > byte_limit {
            return Err(ScalarField3Error::ByteBudgetExceeded {
                required,
                limit: byte_limit,
            });
        }
        let mut out = Vec::new();
        out.try_reserve_exact(required)
            .map_err(|_| ScalarField3Error::AllocationFailed {
                context: "encoded bytes",
                requested: required,
            })?;
        out.extend_from_slice(&MAGIC);
        out.extend_from_slice(&SCALAR_FIELD3_SCHEMA_VERSION.to_le_bytes());
        out.extend_from_slice(&F64_LE_ENCODING.to_le_bytes());
        out.extend_from_slice(&self.layout.code().to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes());
        for dimension in self.dimensions {
            let dimension =
                u64::try_from(dimension).map_err(|_| ScalarField3Error::EncodedSizeOverflow)?;
            out.extend_from_slice(&dimension.to_le_bytes());
        }
        for value in self.origin.into_iter().chain(self.spacing) {
            out.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        let sample_count =
            u64::try_from(self.values.len()).map_err(|_| ScalarField3Error::EncodedSizeOverflow)?;
        out.extend_from_slice(&sample_count.to_le_bytes());
        for value in [
            self.semantics.quantity.len(),
            self.semantics.coordinate_unit.len(),
            self.semantics.value_unit.len(),
        ] {
            let value = u32::try_from(value).map_err(|_| ScalarField3Error::EncodedSizeOverflow)?;
            out.extend_from_slice(&value.to_le_bytes());
        }
        out.extend_from_slice(self.semantics.quantity.as_bytes());
        out.extend_from_slice(self.semantics.coordinate_unit.as_bytes());
        out.extend_from_slice(self.semantics.value_unit.as_bytes());
        for value in &self.values {
            out.extend_from_slice(&value.to_bits().to_le_bytes());
        }
        debug_assert_eq!(out.len(), required);
        Ok(out)
    }

    /// Decode schema-v1 bytes after admitting byte and sample budgets.
    pub fn decode(
        bytes: &[u8],
        sample_limit: usize,
        byte_limit: usize,
    ) -> Result<Self, ScalarField3Error> {
        if bytes.len() > byte_limit {
            return Err(ScalarField3Error::ByteBudgetExceeded {
                required: bytes.len(),
                limit: byte_limit,
            });
        }
        if bytes.len() < FIXED_HEADER_BYTES {
            return Err(ScalarField3Error::Malformed {
                what: "truncated fixed header",
            });
        }
        let mut cursor = Cursor::new(bytes);
        if cursor.read_array::<8>()? != MAGIC {
            return Err(ScalarField3Error::Malformed { what: "bad magic" });
        }
        let version = cursor.read_u32()?;
        if version != SCALAR_FIELD3_SCHEMA_VERSION {
            return Err(ScalarField3Error::UnsupportedSchema { found: version });
        }
        let encoding = cursor.read_u32()?;
        if encoding != F64_LE_ENCODING {
            return Err(ScalarField3Error::UnsupportedEncoding { found: encoding });
        }
        let layout_code = cursor.read_u32()?;
        let layout = match layout_code {
            1 => ScalarLayout3::NodeCentered,
            2 => ScalarLayout3::CellCentered,
            found => return Err(ScalarField3Error::UnsupportedLayout { found }),
        };
        if cursor.read_u32()? != 0 {
            return Err(ScalarField3Error::Malformed {
                what: "reserved header bits are nonzero",
            });
        }
        let mut dimensions = [0usize; 3];
        for (axis, dimension) in dimensions.iter_mut().enumerate() {
            *dimension = usize::try_from(cursor.read_u64()?)
                .map_err(|_| ScalarField3Error::DimensionOutOfRange { axis })?;
        }
        let mut origin = [0.0; 3];
        let mut spacing = [0.0; 3];
        for value in origin.iter_mut().chain(&mut spacing) {
            *value = f64::from_bits(cursor.read_u64()?);
        }
        let encoded_sample_count = usize::try_from(cursor.read_u64()?)
            .map_err(|_| ScalarField3Error::SampleCountOverflow { dimensions })?;
        let semantic_lengths = [
            cursor.read_u32()? as usize,
            cursor.read_u32()? as usize,
            cursor.read_u32()? as usize,
        ];

        let sample_count = validate_layout(layout, dimensions, origin, spacing, sample_limit)?;
        if encoded_sample_count != sample_count {
            return Err(ScalarField3Error::Malformed {
                what: "sample count disagrees with dimensions",
            });
        }
        for (field, length) in ["quantity", "coordinate_unit", "value_unit"]
            .into_iter()
            .zip(semantic_lengths)
        {
            if length == 0 || length > MAX_SEMANTIC_BYTES {
                return Err(ScalarField3Error::InvalidSemantic { field });
            }
        }
        let expected = encoded_len(sample_count, semantic_lengths)?;
        if bytes.len() != expected {
            return Err(ScalarField3Error::Malformed {
                what: "encoded length mismatch",
            });
        }
        let semantics = ScalarFieldSemantics {
            quantity: cursor.read_string(semantic_lengths[0], "quantity")?,
            coordinate_unit: cursor.read_string(semantic_lengths[1], "coordinate_unit")?,
            value_unit: cursor.read_string(semantic_lengths[2], "value_unit")?,
        };
        validate_semantics(&semantics)?;
        let mut values = Vec::new();
        values.try_reserve_exact(sample_count).map_err(|_| {
            ScalarField3Error::AllocationFailed {
                context: "scalar samples",
                requested: sample_count,
            }
        })?;
        for _ in 0..sample_count {
            values.push(f64::from_bits(cursor.read_u64()?));
        }
        debug_assert_eq!(cursor.position, bytes.len());
        Self::new(
            layout,
            dimensions,
            origin,
            spacing,
            semantics,
            sample_limit,
            values,
        )
    }

    /// Convert a node-centered artifact into the isosurface grid type.
    pub fn into_node_grid(self, node_limit: usize) -> Result<Grid3, ScalarField3Error> {
        if self.layout != ScalarLayout3::NodeCentered {
            return Err(ScalarField3Error::NotNodeCentered);
        }
        let upper = self.world_bounds()[1];
        Grid3::from_values(self.dimensions, self.origin, upper, node_limit, self.values)
            .map_err(ScalarField3Error::Grid)
    }

    fn encoded_len(&self) -> Result<usize, ScalarField3Error> {
        encoded_len(
            self.values.len(),
            [
                self.semantics.quantity.len(),
                self.semantics.coordinate_unit.len(),
                self.semantics.value_unit.len(),
            ],
        )
    }
}

fn validate_layout(
    layout: ScalarLayout3,
    dimensions: [usize; 3],
    origin: [f64; 3],
    spacing: [f64; 3],
    sample_limit: usize,
) -> Result<usize, ScalarField3Error> {
    if dimensions.contains(&0) {
        return Err(ScalarField3Error::InvalidDimensions { dimensions });
    }
    let sample_count = dimensions
        .into_iter()
        .try_fold(1usize, |count, dimension| count.checked_mul(dimension))
        .ok_or(ScalarField3Error::SampleCountOverflow { dimensions })?;
    if sample_count > sample_limit {
        return Err(ScalarField3Error::SampleBudgetExceeded {
            required: sample_count,
            limit: sample_limit,
        });
    }
    let upper = upper_bound(layout, dimensions, origin, spacing);
    for axis in 0..3 {
        let multiplier = match layout {
            ScalarLayout3::NodeCentered => dimensions[axis] - 1,
            ScalarLayout3::CellCentered => dimensions[axis],
        };
        if !origin[axis].is_finite()
            || !spacing[axis].is_finite()
            || spacing[axis] <= 0.0
            || !upper[axis].is_finite()
            || (multiplier > 0 && upper[axis] <= origin[axis])
        {
            return Err(ScalarField3Error::InvalidGeometry { axis });
        }
    }
    Ok(sample_count)
}

fn upper_bound(
    layout: ScalarLayout3,
    dimensions: [usize; 3],
    origin: [f64; 3],
    spacing: [f64; 3],
) -> [f64; 3] {
    let multiplier = |axis: usize| match layout {
        ScalarLayout3::NodeCentered => dimensions[axis] - 1,
        ScalarLayout3::CellCentered => dimensions[axis],
    };
    [
        spacing[0].mul_add(multiplier(0) as f64, origin[0]),
        spacing[1].mul_add(multiplier(1) as f64, origin[1]),
        spacing[2].mul_add(multiplier(2) as f64, origin[2]),
    ]
}

fn validate_semantics(semantics: &ScalarFieldSemantics) -> Result<(), ScalarField3Error> {
    for (field, value) in [
        ("quantity", semantics.quantity.as_str()),
        ("coordinate_unit", semantics.coordinate_unit.as_str()),
        ("value_unit", semantics.value_unit.as_str()),
    ] {
        if value.is_empty()
            || value.len() > MAX_SEMANTIC_BYTES
            || value.chars().any(char::is_control)
        {
            return Err(ScalarField3Error::InvalidSemantic { field });
        }
    }
    Ok(())
}

fn encoded_len(
    sample_count: usize,
    semantic_lengths: [usize; 3],
) -> Result<usize, ScalarField3Error> {
    let semantic_bytes = semantic_lengths
        .into_iter()
        .try_fold(0usize, |count, length| count.checked_add(length))
        .ok_or(ScalarField3Error::EncodedSizeOverflow)?;
    let sample_bytes = sample_count
        .checked_mul(core::mem::size_of::<f64>())
        .ok_or(ScalarField3Error::EncodedSizeOverflow)?;
    FIXED_HEADER_BYTES
        .checked_add(semantic_bytes)
        .and_then(|size| size.checked_add(sample_bytes))
        .ok_or(ScalarField3Error::EncodedSizeOverflow)
}

struct Cursor<'a> {
    bytes: &'a [u8],
    position: usize,
}

impl<'a> Cursor<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, position: 0 }
    }

    fn read_array<const N: usize>(&mut self) -> Result<[u8; N], ScalarField3Error> {
        let end = self
            .position
            .checked_add(N)
            .ok_or(ScalarField3Error::Malformed {
                what: "cursor overflow",
            })?;
        let value = self
            .bytes
            .get(self.position..end)
            .ok_or(ScalarField3Error::Malformed {
                what: "truncated payload",
            })?;
        self.position = end;
        Ok(value.try_into().expect("slice length is N"))
    }

    fn read_u32(&mut self) -> Result<u32, ScalarField3Error> {
        Ok(u32::from_le_bytes(self.read_array()?))
    }

    fn read_u64(&mut self) -> Result<u64, ScalarField3Error> {
        Ok(u64::from_le_bytes(self.read_array()?))
    }

    fn read_string(
        &mut self,
        length: usize,
        field: &'static str,
    ) -> Result<String, ScalarField3Error> {
        let end = self
            .position
            .checked_add(length)
            .ok_or(ScalarField3Error::Malformed {
                what: "semantic range overflow",
            })?;
        let value = self
            .bytes
            .get(self.position..end)
            .ok_or(ScalarField3Error::Malformed {
                what: "truncated semantics",
            })?;
        self.position = end;
        std::str::from_utf8(value)
            .map(str::to_owned)
            .map_err(|_| ScalarField3Error::InvalidSemantic { field })
    }
}
