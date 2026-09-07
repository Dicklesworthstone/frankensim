//! Declared elastic, strain and stress tensor coordinates. Numerical rotation belongs to
//! the constitutive/operator layer; this module only preserves meaning.

use fs_blake3::ContentHash;

use crate::MatDbError;

/// Stress/strain vector convention for an elastic matrix [Pa].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElasticTensorNotation {
    /// Physical tensor shear for both vectors; shear columns include factor two.
    Tensor,
    /// Tensor stress shear and engineering strain shear `2 epsilon_ij`.
    Engineering,
    /// Both vectors store `sqrt(2)` times physical tensor shear.
    Mandel,
}

/// Shared row/column order of an elastic matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElasticTensorOrder {
    /// `[xx, yy, zz, xy, yz, zx]`.
    XxYyZzXyYzZx,
    /// `[xx, yy, zz, yz, zx, xy]`.
    XxYyZzYzZxXy,
}

/// Source basis; no frame or shear convention is inferred from a material name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ElasticTensorBasis {
    /// Stress/strain shear convention.
    pub notation: ElasticTensorNotation,
    /// Component order.
    pub order: ElasticTensorOrder,
    /// Nonzero source coordinate-frame identity.
    pub frame: ContentHash,
}

/// Shear coordinates of a symmetric infinitesimal STRAIN tensor.
/// These are not stress-vector or fourth-order stiffness conventions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrainTensorNotation {
    /// Store the physical off-diagonal strain `epsilon_ij` once.
    Tensor,
    /// Store engineering shear strain `2 epsilon_ij`.
    Engineering,
    /// Store `sqrt(2) epsilon_ij`; Euclidean norm equals Frobenius norm.
    Mandel,
}

/// Coordinates of a symmetric second-order strain tensor. Its six entries
/// declare symmetry by construction; no displacement gradient, antisymmetric
/// part or coefficient of thermal expansion is inferred from this descriptor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StrainTensorBasis {
    /// Strain-specific shear scaling.
    pub notation: StrainTensorNotation,
    /// Same six symmetric coordinate positions used by elastic rows/columns.
    pub order: ElasticTensorOrder,
    /// Nonzero source-frame identity, validated by the numerical consumer.
    pub frame: ContentHash,
}

/// Coordinates of a symmetric small-strain Cauchy stress tensor [Pa].
/// Engineering stress stores physical shear; only engineering STRAIN doubles it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StressTensorNotation {
    /// Physical off-diagonal stress `sigma_ij`, stored once.
    Tensor,
    /// Physical shear stress conjugate to engineering strain `2 epsilon_ij`.
    Engineering,
    /// `sqrt(2) sigma_ij`, conjugate to Mandel strain.
    Mandel,
}

/// Explicit coordinates for six symmetric stress components. This does not
/// describe first/second Piola stress or provide a finite-deformation mapping.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StressTensorBasis {
    /// Stress-specific shear scaling.
    pub notation: StressTensorNotation,
    /// Order of the six symmetric stress coordinates.
    pub order: ElasticTensorOrder,
    /// Nonzero coordinate-frame identity, checked by the numerical consumer.
    pub frame: ContentHash,
}

/// One coordinate of a source-declared symmetric Cauchy stress [Pa]. Grouping
/// six components does not establish their correlation or material calibration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StressTensorComponent {
    basis: StressTensorBasis,
    source_tensor: ContentHash,
    index: u8,
}

impl StressTensorComponent {
    /// Three convention/index bytes and two complete identities.
    pub const ENCODED_LEN: usize = 67;

    /// Declare a bounded coordinate with nonzero frame and source identities.
    pub fn new(
        basis: StressTensorBasis,
        source_tensor: ContentHash,
        index: u8,
    ) -> Result<Self, MatDbError> {
        if basis.frame == ContentHash([0; 32]) || source_tensor == ContentHash([0; 32]) {
            return Err(MatDbError::InvalidTensorContext {
                reason: "stress frame and source tensor identities must be nonzero",
            });
        }
        if index >= 6 {
            return Err(MatDbError::InvalidTensorContext {
                reason: "stress component index must lie in 0..6",
            });
        }
        Ok(Self {
            basis,
            source_tensor,
            index,
        })
    }

    /// Stress convention, coordinate order and source frame.
    #[must_use]
    pub const fn basis(self) -> StressTensorBasis {
        self.basis
    }

    /// Source-declared tensor grouping, distinct from a component claim ID.
    #[must_use]
    pub const fn source_tensor(self) -> ContentHash {
        self.source_tensor
    }

    /// Position in the explicitly declared order.
    #[must_use]
    pub const fn index(self) -> usize {
        self.index as usize
    }

    /// Canonical coordinates; the enclosing claim/pack owns versioning and
    /// distinguishes this descriptor from strain or stiffness descriptors.
    #[must_use]
    pub fn canonical_bytes(self) -> [u8; Self::ENCODED_LEN] {
        let mut bytes = [0; Self::ENCODED_LEN];
        bytes[0] = match self.basis.notation {
            StressTensorNotation::Tensor => 0,
            StressTensorNotation::Engineering => 1,
            StressTensorNotation::Mandel => 2,
        };
        bytes[1] = match self.basis.order {
            ElasticTensorOrder::XxYyZzXyYzZx => 0,
            ElasticTensorOrder::XxYyZzYzZxXy => 1,
        };
        bytes[2] = self.index;
        bytes[3..35].copy_from_slice(self.basis.frame.as_bytes());
        bytes[35..].copy_from_slice(self.source_tensor.as_bytes());
        bytes
    }

    /// Decode exact-length bytes, refusing unknown conventions and identities.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, MatDbError> {
        let invalid = || MatDbError::InvalidTensorContext {
            reason: "invalid stress component descriptor length or convention tag",
        };
        if bytes.len() != Self::ENCODED_LEN {
            return Err(invalid());
        }
        let notation = match bytes[0] {
            0 => StressTensorNotation::Tensor,
            1 => StressTensorNotation::Engineering,
            2 => StressTensorNotation::Mandel,
            _ => return Err(invalid()),
        };
        let order = match bytes[1] {
            0 => ElasticTensorOrder::XxYyZzXyYzZx,
            1 => ElasticTensorOrder::XxYyZzYzZxXy,
            _ => return Err(invalid()),
        };
        Self::new(
            StressTensorBasis {
                notation,
                order,
                frame: ContentHash(bytes[3..35].try_into().map_err(|_| invalid())?),
            },
            ContentHash(bytes[35..].try_into().map_err(|_| invalid())?),
            bytes[2],
        )
    }
}

/// One of six explicitly supplied symmetric strain coordinates. The source
/// tensor identity groups components; it does not establish their correlation
/// or identify a thermal expansion coefficient or temperature integration law.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StrainTensorComponent {
    basis: StrainTensorBasis,
    source_tensor: ContentHash,
    index: u8,
}

impl StrainTensorComponent {
    /// Three convention/index bytes and two complete identities.
    pub const ENCODED_LEN: usize = 67;

    pub fn new(
        basis: StrainTensorBasis,
        source_tensor: ContentHash,
        index: u8,
    ) -> Result<Self, MatDbError> {
        if basis.frame == ContentHash([0; 32]) || source_tensor == ContentHash([0; 32]) {
            return Err(MatDbError::InvalidTensorContext {
                reason: "strain frame and source tensor identities must be nonzero",
            });
        }
        if index >= 6 {
            return Err(MatDbError::InvalidTensorContext {
                reason: "strain component index must lie in 0..6",
            });
        }
        Ok(Self {
            basis,
            source_tensor,
            index,
        })
    }

    #[must_use]
    pub const fn basis(self) -> StrainTensorBasis {
        self.basis
    }

    #[must_use]
    pub const fn source_tensor(self) -> ContentHash {
        self.source_tensor
    }

    /// Coordinate position in the explicitly declared order.
    #[must_use]
    pub const fn index(self) -> usize {
        self.index as usize
    }

    /// Symmetry is implicit in the six-coordinate type. The enclosing format
    /// owns versioning; the frozen descriptor retains exact conventions.
    #[must_use]
    pub fn canonical_bytes(self) -> [u8; Self::ENCODED_LEN] {
        let mut bytes = [0; Self::ENCODED_LEN];
        bytes[0] = match self.basis.notation {
            StrainTensorNotation::Tensor => 0,
            StrainTensorNotation::Engineering => 1,
            StrainTensorNotation::Mandel => 2,
        };
        bytes[1] = match self.basis.order {
            ElasticTensorOrder::XxYyZzXyYzZx => 0,
            ElasticTensorOrder::XxYyZzYzZxXy => 1,
        };
        bytes[2] = self.index;
        bytes[3..35].copy_from_slice(self.basis.frame.as_bytes());
        bytes[35..].copy_from_slice(self.source_tensor.as_bytes());
        bytes
    }

    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, MatDbError> {
        let invalid = || MatDbError::InvalidTensorContext {
            reason: "invalid strain component descriptor length or convention tag",
        };
        if bytes.len() != Self::ENCODED_LEN {
            return Err(invalid());
        }
        let notation = match bytes[0] {
            0 => StrainTensorNotation::Tensor,
            1 => StrainTensorNotation::Engineering,
            2 => StrainTensorNotation::Mandel,
            _ => return Err(invalid()),
        };
        let order = match bytes[1] {
            0 => ElasticTensorOrder::XxYyZzXyYzZx,
            1 => ElasticTensorOrder::XxYyZzYzZxXy,
            _ => return Err(invalid()),
        };
        Self::new(
            StrainTensorBasis {
                notation,
                order,
                frame: ContentHash(bytes[3..35].try_into().map_err(|_| invalid())?),
            },
            ContentHash(bytes[35..].try_into().map_err(|_| invalid())?),
            bytes[2],
        )
    }
}

/// Declared symmetry of the complete law, checked numerically by its consumer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElasticTensorSymmetry {
    /// Small-strain elastic law with both minor symmetries and major reciprocity.
    MajorMinor,
}

/// One coefficient of a named complete source tensor. All 36 coefficients are
/// explicit, including zeros. A common source identity groups their declarations;
/// it is not proof of measured correlation, symmetry or positive definiteness.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ElasticTensorComponent {
    basis: ElasticTensorBasis,
    symmetry: ElasticTensorSymmetry,
    source_tensor: ContentHash,
    row: u8,
    column: u8,
}

impl ElasticTensorComponent {
    /// Fixed v1 descriptor size: five tags/indices and two full identities.
    pub const ENCODED_LEN: usize = 69;

    /// Declare one bounded component; neither identity may be zero.
    pub fn new(
        basis: ElasticTensorBasis,
        symmetry: ElasticTensorSymmetry,
        source_tensor: ContentHash,
        row: u8,
        column: u8,
    ) -> Result<Self, MatDbError> {
        if basis.frame == ContentHash([0; 32]) || source_tensor == ContentHash([0; 32]) {
            return Err(MatDbError::InvalidTensorContext {
                reason: "elastic frame and source tensor identities must be nonzero",
            });
        }
        if row >= 6 || column >= 6 {
            return Err(MatDbError::InvalidTensorContext {
                reason: "elastic component indices must lie in 0..6",
            });
        }
        Ok(Self {
            basis,
            symmetry,
            source_tensor,
            row,
            column,
        })
    }

    /// Source basis shared by every coefficient of the complete tensor.
    #[must_use]
    pub const fn basis(self) -> ElasticTensorBasis {
        self.basis
    }

    /// Declared symmetry class, not a numerical certificate.
    #[must_use]
    pub const fn symmetry(self) -> ElasticTensorSymmetry {
        self.symmetry
    }

    /// Source-declared tensor group; independent of each coefficient's claim id.
    #[must_use]
    pub const fn source_tensor(self) -> ContentHash {
        self.source_tensor
    }

    /// Zero-based matrix row and column in the declared order.
    #[must_use]
    pub const fn indices(self) -> (usize, usize) {
        (self.row as usize, self.column as usize)
    }

    /// Canonical descriptor bytes. The enclosing claim/pack owns versioning.
    #[must_use]
    pub fn canonical_bytes(self) -> [u8; Self::ENCODED_LEN] {
        let mut bytes = [0; Self::ENCODED_LEN];
        bytes[0] = match self.basis.notation {
            ElasticTensorNotation::Tensor => 0,
            ElasticTensorNotation::Engineering => 1,
            ElasticTensorNotation::Mandel => 2,
        };
        bytes[1] = match self.basis.order {
            ElasticTensorOrder::XxYyZzXyYzZx => 0,
            ElasticTensorOrder::XxYyZzYzZxXy => 1,
        };
        bytes[2] = match self.symmetry {
            ElasticTensorSymmetry::MajorMinor => 0,
        };
        bytes[3] = self.row;
        bytes[4] = self.column;
        bytes[5..37].copy_from_slice(self.basis.frame.as_bytes());
        bytes[37..].copy_from_slice(self.source_tensor.as_bytes());
        bytes
    }

    /// Decode exact-length canonical bytes, refusing unknown conventions.
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, MatDbError> {
        let invalid = || MatDbError::InvalidTensorContext {
            reason: "invalid elastic component descriptor length or convention tag",
        };
        if bytes.len() != Self::ENCODED_LEN {
            return Err(invalid());
        }
        let notation = match bytes[0] {
            0 => ElasticTensorNotation::Tensor,
            1 => ElasticTensorNotation::Engineering,
            2 => ElasticTensorNotation::Mandel,
            _ => return Err(invalid()),
        };
        let order = match bytes[1] {
            0 => ElasticTensorOrder::XxYyZzXyYzZx,
            1 => ElasticTensorOrder::XxYyZzYzZxXy,
            _ => return Err(invalid()),
        };
        let symmetry = match bytes[2] {
            0 => ElasticTensorSymmetry::MajorMinor,
            _ => return Err(invalid()),
        };
        Self::new(
            ElasticTensorBasis {
                notation,
                order,
                frame: ContentHash(bytes[5..37].try_into().map_err(|_| invalid())?),
            },
            symmetry,
            ContentHash(bytes[37..].try_into().map_err(|_| invalid())?),
            bytes[3],
            bytes[4],
        )
    }
}
