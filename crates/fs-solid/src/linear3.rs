//! Three-dimensional small-strain elasticity on body-fitted linear
//! tetrahedra.
//!
//! This module owns the generic volume operator needed by solid vibration,
//! transient deformation, and coupled acoustics. It assembles the standard
//! Galerkin stiffness and consistent mass matrices
//!
//! `K = sum_e integral B^T C B dV`, `M = sum_e integral rho N^T N dV`
//!
//! for arbitrary symmetric positive-definite linear-elastic material states.
//! Isotropic and oriented-orthotropic constructors lower their constitutive
//! laws into the same Mandel-basis tensor, so no object class or material name
//! is encoded: a squat metal disc, a wooden body, a crystal, or a glass
//! support are geometry plus resolved material fields. The output is the
//! ordinary `(K,M)` pencil consumed by `fs-modal`.
//!
//! The v1 applicability boundary is intentionally narrow and explicit:
//! infinitesimal strain, symmetric positive-definite linear elasticity,
//! undeformed geometry, and conforming non-degenerate P1 tetrahedra.
//! A separately resolved total thermal eigenstrain field can be assembled into
//! its work-conjugate equivalent nodal load. The caller must still solve the
//! constrained displacement field and explicitly update or remesh geometry.
//! Plasticity, finite strain, phase change, and evolving topology must select a
//! different constitutive/kinematic rung rather than silently passing through
//! this operator.

use std::collections::BTreeSet;

use crate::linear::Jacobi;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_material::OrthotropicElastic;
use fs_material::state_point::{
    CorrelationUnknownReason, ElasticTensorStatePoint, IntegratedIsotropicThermalExpansion,
    IsotropicElasticStatePoint, IsotropicSolidStatePoint, JointCorrelation,
    JointStrainTensorStatePoint, OrthotropicElasticStatePoint, StrainTensorStatePoint,
};
pub use fs_material::state_point::{ElasticTensorBasis, ElasticTensorNotation, ElasticTensorOrder};
pub use fs_material::tensor::{
    StrainTensorBasis, StrainTensorNotation, StressTensorBasis, StressTensorNotation,
};
use fs_solver::krylov::CgState;
use fs_solver::op::{CsrOp, LinearOp};
use fs_sparse::{Coo, Csr};

/// Bounded assembly request. Limits are checked before matrix allocation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TetAssemblyBudget {
    /// Maximum admitted mesh vertices.
    pub maximum_nodes: usize,
    /// Maximum admitted tetrahedra.
    pub maximum_tetrahedra: usize,
    /// Maximum admitted free displacement degrees of freedom.
    pub maximum_free_dofs: usize,
    /// Minimum `abs(det(J)) / h_max^3` admitted for any tetrahedron.
    pub minimum_scaled_jacobian: f64,
}

impl TetAssemblyBudget {
    /// A conservative interactive/default envelope.
    #[must_use]
    pub const fn standard() -> Self {
        Self {
            maximum_nodes: 1_000_000,
            maximum_tetrahedra: 6_000_000,
            maximum_free_dofs: 3_000_000,
            minimum_scaled_jacobian: 1.0e-12,
        }
    }
}

impl Default for TetAssemblyBudget {
    fn default() -> Self {
        Self::standard()
    }
}

/// Mandel ordering used by the three-dimensional operator:
/// `[xx, yy, zz, sqrt(2)xy, sqrt(2)yz, sqrt(2)zx]`.
///
/// Mandel scaling makes the ordinary six-vector dot product equal to tensor
/// double contraction. Consequently an elastic matrix is symmetric positive
/// definite in the usual Euclidean sense, including after material-frame
/// rotations.
pub type MandelStiffness6 = [[f64; 6]; 6];

/// Retained numerical transform and its executable, admitted material.
///
/// This records a major/minor-symmetric, positive-definite elasticity law,
/// not an arbitrary coupled operator. It does not certify observations or
/// frame calibration, and currently has no portable byte codec.
#[derive(Debug, Clone, PartialEq)]
pub struct ElasticTensorTransformReceipt {
    source_stiffness_pa: [[f64; 6]; 6],
    source_basis: ElasticTensorBasis,
    target_frame: ContentHash,
    source_to_target: [[f64; 3]; 3],
    source_material_identity: ContentHash,
    material: TetElasticMaterial,
}

impl ElasticTensorTransformReceipt {
    /// Complete source matrix, preserving its declared order and shear scaling.
    #[must_use]
    pub const fn source_stiffness_pa(&self) -> &[[f64; 6]; 6] {
        &self.source_stiffness_pa
    }

    /// Declared source notation, component ordering and frame identity.
    #[must_use]
    pub const fn source_basis(&self) -> ElasticTensorBasis {
        self.source_basis
    }

    /// Target coordinate-frame identity.
    #[must_use]
    pub const fn target_frame(&self) -> ContentHash {
        self.target_frame
    }

    /// Proper rotation: columns are source axes expressed in the target frame.
    #[must_use]
    pub const fn source_to_target(&self) -> &[[f64; 3]; 3] {
        &self.source_to_target
    }

    /// Upstream identity supplied for the admitted material data.
    #[must_use]
    pub const fn source_material_identity(&self) -> ContentHash {
        self.source_material_identity
    }

    /// Target-frame material consumed directly by the existing tet operator.
    /// Its identity binds all source/transform fields, density and output bits.
    #[must_use]
    pub const fn material(&self) -> &TetElasticMaterial {
        &self.material
    }
}

/// One uniform linear-elastic state carried into the element operator.
#[derive(Debug, Clone, PartialEq)]
pub struct TetElasticMaterial {
    /// Density [kg/m^3].
    pub density_kg_m3: f64,
    /// Symmetric positive-definite constitutive tensor [Pa] in the global
    /// frame and the module's declared Mandel ordering.
    stiffness_mandel_pa: MandelStiffness6,
    /// Raw bytes of the evidence-bearing `fs-material` state identity.
    pub material_state_identity: [u8; 32],
}

impl TetElasticMaterial {
    /// Construct an isotropic state from scalar properties and an upstream
    /// material-state identity.
    ///
    /// The identity is mandatory because otherwise changing a temperature-
    /// dependent property could leave no trace in downstream modal or
    /// acoustic artifacts.
    pub fn try_new(
        density_kg_m3: f64,
        young_modulus_pa: f64,
        poisson_ratio: f64,
        material_state_identity: [u8; 32],
    ) -> Result<Self, TetElasticError> {
        if !(young_modulus_pa.is_finite() && young_modulus_pa > 0.0) {
            return Err(TetElasticError::InvalidMaterial {
                what: "young_modulus_pa must be finite and positive",
            });
        }
        if !(poisson_ratio.is_finite() && poisson_ratio > -1.0 && poisson_ratio < 0.5) {
            return Err(TetElasticError::InvalidMaterial {
                what: "poisson_ratio must lie in (-1, 0.5)",
            });
        }
        let lambda = young_modulus_pa * poisson_ratio
            / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
        let mu = young_modulus_pa / (2.0 * (1.0 + poisson_ratio));
        let mut stiffness = [[0.0; 6]; 6];
        for (row, values) in stiffness.iter_mut().enumerate().take(3) {
            for value in values.iter_mut().take(3) {
                *value = lambda;
            }
            values[row] += 2.0 * mu;
        }
        for (index, row) in stiffness.iter_mut().enumerate().skip(3) {
            row[index] = 2.0 * mu;
        }
        Self::try_new_mandel(density_kg_m3, stiffness, material_state_identity)
    }

    /// Construct from a complete global-frame Mandel stiffness tensor.
    ///
    /// This is the material-name-independent admission seam for anisotropic
    /// solids. The matrix must be finite, exactly symmetric, and positive
    /// definite. Callers that obtain paired tensor entries from noisy data
    /// must perform and identify their chosen symmetry projection upstream;
    /// this operator never repairs constitutive evidence silently.
    pub fn try_new_mandel(
        density_kg_m3: f64,
        stiffness_mandel_pa: MandelStiffness6,
        material_state_identity: [u8; 32],
    ) -> Result<Self, TetElasticError> {
        let material = Self {
            density_kg_m3,
            stiffness_mandel_pa,
            material_state_identity,
        };
        material.validate()?;
        Ok(material)
    }

    /// Admit a full anisotropic elastic matrix with explicit basis and frames.
    ///
    /// Converts notation/order before applying the fourth-order rotation
    /// `C'_{ijkl} = Q_{ip} Q_{jq} Q_{kr} Q_{ls} C_{pqrs}` through its Mandel
    /// representation. Major symmetry is checked in work-conjugate coordinates,
    /// not by assuming a tensor-shear matrix is ordinarily symmetric. Both the
    /// complete source law and transformed law must be positive definite.
    /// No noisy-data projection or missing-coefficient completion is performed.
    pub fn try_new_oriented_tensor(
        density_kg_m3: f64,
        stiffness_pa: [[f64; 6]; 6],
        source_basis: ElasticTensorBasis,
        target_frame: ContentHash,
        source_to_target: [[f64; 3]; 3],
        source_material_identity: ContentHash,
    ) -> Result<ElasticTensorTransformReceipt, TetElasticError> {
        validate_rotation(source_to_target)?;
        if source_basis.frame == ContentHash([0; 32]) || target_frame == ContentHash([0; 32]) {
            return Err(TetElasticError::InvalidMaterial {
                what: "source and target tensor frame identities must be nonzero",
            });
        }
        if source_basis.frame == target_frame
            && source_to_target != [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        {
            return Err(TetElasticError::InvalidMaterial {
                what: "a frame mapped to itself requires the identity rotation",
            });
        }
        let principal = normalize_elastic_tensor(stiffness_pa, source_basis)?;
        Self::try_new_mandel(density_kg_m3, principal, source_material_identity.0)?;
        let world = rotate_mandel_stiffness(principal, source_to_target);
        let mut identity =
            DomainHasher::new("org.frankensim.fs-solid.elastic-tensor-frame-transform.v1");
        identity.update(source_material_identity.as_bytes());
        identity.update(source_basis.frame.as_bytes());
        identity.update(target_frame.as_bytes());
        identity.update(&[match source_basis.notation {
            ElasticTensorNotation::Tensor => 0,
            ElasticTensorNotation::Engineering => 1,
            ElasticTensorNotation::Mandel => 2,
        }]);
        identity.update(&[match source_basis.order {
            ElasticTensorOrder::XxYyZzXyYzZx => 0,
            ElasticTensorOrder::XxYyZzYzZxXy => 1,
        }]);
        identity.update(&density_kg_m3.to_bits().to_le_bytes());
        for value in stiffness_pa
            .iter()
            .flatten()
            .chain(source_to_target.iter().flatten())
            .chain(world.iter().flatten())
        {
            identity.update(&value.to_bits().to_le_bytes());
        }
        let material = Self::try_new_mandel(density_kg_m3, world, identity.finalize().0)?;
        Ok(ElasticTensorTransformReceipt {
            source_stiffness_pa: stiffness_pa,
            source_basis,
            target_frame,
            source_to_target,
            source_material_identity,
            material,
        })
    }

    /// Admit a complete material-card tensor and rotate it into the solver frame.
    /// The upstream bundle identity binds all selected coefficient receipts,
    /// density and the exact physical query point. Numerical admission checks
    /// the whole law; component uncertainty is retained upstream, not promoted
    /// to a certified rotated uncertainty bound.
    pub fn from_resolved_elastic_tensor(
        state: &ElasticTensorStatePoint,
        target_frame: ContentHash,
        source_to_target: [[f64; 3]; 3],
    ) -> Result<ElasticTensorTransformReceipt, TetElasticError> {
        Self::try_new_oriented_tensor(
            state.density_kg_m3(),
            *state.stiffness_pa(),
            state.basis(),
            target_frame,
            source_to_target,
            state.resolved().identity(),
        )
    }

    /// Construct an oriented orthotropic state.
    ///
    /// `principal_to_world` is a proper orthonormal rotation whose columns are
    /// the material principal axes expressed in world coordinates. The
    /// constitutive law is rotated exactly once at admission; element hot
    /// loops consume the resulting global-frame tensor without dispatch.
    pub fn try_new_oriented_orthotropic(
        density_kg_m3: f64,
        law: &OrthotropicElastic,
        principal_to_world: [[f64; 3]; 3],
        material_state_identity: [u8; 32],
        orientation_identity: [u8; 32],
    ) -> Result<Self, TetElasticError> {
        validate_rotation(principal_to_world)?;
        if orientation_identity == [0; 32] {
            return Err(TetElasticError::InvalidMaterial {
                what: "orientation_identity must not be the zero identity",
            });
        }
        let principal = law.stiffness();
        let world = rotate_mandel_stiffness(principal, principal_to_world);
        let mut identity =
            DomainHasher::new("org.frankensim.fs-solid.oriented-orthotropic-constitutive-state.v1");
        identity.update(&material_state_identity);
        identity.update(&orientation_identity);
        for row in principal_to_world {
            for value in row {
                identity.update(&value.to_bits().to_le_bytes());
            }
        }
        Self::try_new_mandel(density_kg_m3, world, *identity.finalize().as_bytes())
    }

    /// Bind directly to the canonical temperature/environment-resolved
    /// isotropic material state.
    #[must_use]
    pub fn from_resolved_state(state: &IsotropicSolidStatePoint) -> Self {
        Self::try_new(
            state.density_kg_m3(),
            state.young_modulus_pa(),
            state.poisson_ratio(),
            *state.resolved().identity().as_bytes(),
        )
        .expect("an admitted isotropic material state satisfies tet invariants")
    }

    /// Bind the minimal evidence-bearing isotropic-elastic state used by
    /// vibration analysis, without requiring an unrelated yield datum.
    #[must_use]
    pub fn from_resolved_elastic_state(state: &IsotropicElasticStatePoint) -> Self {
        Self::try_new(
            state.density_kg_m3(),
            state.young_modulus_pa(),
            state.poisson_ratio(),
            *state.resolved().identity().as_bytes(),
        )
        .expect("an admitted isotropic-elastic state satisfies tet invariants")
    }

    /// Bind an evidence-bearing orthotropic state and an independently
    /// identified material-axis orientation to the global operator frame.
    pub fn try_from_resolved_orthotropic_state(
        state: &OrthotropicElasticStatePoint,
        principal_to_world: [[f64; 3]; 3],
        orientation_identity: ContentHash,
    ) -> Result<Self, TetElasticError> {
        Self::try_new_oriented_orthotropic(
            state.density_kg_m3(),
            state.law(),
            principal_to_world,
            *state.resolved().identity().as_bytes(),
            *orientation_identity.as_bytes(),
        )
    }

    /// Complete global-frame Mandel stiffness tensor [Pa].
    #[must_use]
    pub const fn stiffness_mandel_pa(&self) -> &MandelStiffness6 {
        &self.stiffness_mandel_pa
    }

    fn validate(&self) -> Result<(), TetElasticError> {
        if !(self.density_kg_m3.is_finite() && self.density_kg_m3 > 0.0) {
            return Err(TetElasticError::InvalidMaterial {
                what: "density_kg_m3 must be finite and positive",
            });
        }
        if self
            .stiffness_mandel_pa
            .iter()
            .flatten()
            .any(|value| !value.is_finite())
        {
            return Err(TetElasticError::InvalidMaterial {
                what: "stiffness_mandel_pa must be finite",
            });
        }
        for row in 0..6 {
            for column in 0..row {
                if self.stiffness_mandel_pa[row][column].to_bits()
                    != self.stiffness_mandel_pa[column][row].to_bits()
                {
                    return Err(TetElasticError::InvalidMaterial {
                        what: "stiffness_mandel_pa must be exactly symmetric",
                    });
                }
            }
        }
        if !is_positive_definite(self.stiffness_mandel_pa) {
            return Err(TetElasticError::InvalidMaterial {
                what: "stiffness_mandel_pa must be positive definite",
            });
        }
        if self.material_state_identity == [0; 32] {
            return Err(TetElasticError::InvalidMaterial {
                what: "material_state_identity must not be the zero identity",
            });
        }
        Ok(())
    }
}

/// Material field over tetrahedra. A one-card solid uses `Uniform`; a
/// multi-material or spatially varying state uses one entry per element.
#[derive(Debug, Clone, Copy)]
pub enum TetMaterialField<'a> {
    /// One state for every tetrahedron.
    Uniform(&'a TetElasticMaterial),
    /// One state per tetrahedron, in connectivity order.
    PerElement(&'a [TetElasticMaterial]),
}

impl TetMaterialField<'_> {
    fn validate(&self, element_count: usize) -> Result<(), TetElasticError> {
        match self {
            Self::Uniform(material) => material.validate(),
            Self::PerElement(materials) => {
                if materials.len() != element_count {
                    return Err(TetElasticError::MaterialCountMismatch {
                        expected: element_count,
                        got: materials.len(),
                    });
                }
                for material in *materials {
                    material.validate()?;
                }
                Ok(())
            }
        }
    }

    fn at(&self, element: usize) -> &TetElasticMaterial {
        match self {
            Self::Uniform(material) => material,
            Self::PerElement(materials) => &materials[element],
        }
    }
}

/// Total stress-free thermal strain at one resolved material state.
///
/// The six components use this module's Mandel ordering.  The strain is the
/// already-integrated change from `reference_temperature_k` to
/// `temperature_k`; this operator does not assume a constant coefficient of
/// thermal expansion.  An upstream material law may therefore integrate a
/// nonlinear or anisotropic expansion tensor before constructing this value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TetThermalStrainState {
    temperature_k: f64,
    reference_temperature_k: f64,
    free_strain_mandel: [f64; 6],
    elastic_material_state_identity: [u8; 32],
    thermal_law_identity: ContentHash,
    maximum_absolute_mandel_strain: f64,
    identity: ContentHash,
}

/// A second-order thermal strain transform and its directly executable state.
/// The input is already integrated over the temperature path. This receipt
/// retains nominal values and frame declarations, not calibrated orientation
/// or propagated uncertainty, and currently has no portable byte codec.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalStrainTransformReceipt {
    source_free_strain: [f64; 6],
    source_basis: StrainTensorBasis,
    target_frame: ContentHash,
    source_to_target: [[f64; 3]; 3],
    state: TetThermalStrainState,
}

impl ThermalStrainTransformReceipt {
    /// Source values in the declared coordinate order and shear convention.
    #[must_use]
    pub const fn source_free_strain(&self) -> &[f64; 6] {
        &self.source_free_strain
    }

    /// Symmetric strain coordinates and source-frame identity.
    #[must_use]
    pub const fn source_basis(&self) -> StrainTensorBasis {
        self.source_basis
    }

    /// The admitted elastic material's target frame.
    #[must_use]
    pub const fn target_frame(&self) -> ContentHash {
        self.target_frame
    }

    /// Proper rotation: columns are source axes expressed in the target frame.
    #[must_use]
    pub const fn source_to_target(&self) -> &[[f64; 3]; 3] {
        &self.source_to_target
    }

    /// Target-frame strain accepted directly by the tet thermal operator.
    /// Its identity binds this complete transform and both physical states.
    #[must_use]
    pub const fn state(&self) -> &TetThermalStrainState {
        &self.state
    }
}

/// Conditional second moments in the target frame, with source correlation
/// retained. Stiffness, frame rotation and temperatures are treated as fixed;
/// covariance alone does not bound realizations inside the material domain.
#[derive(Debug, Clone, PartialEq)]
pub enum ThermalStrainUncertainty {
    Covariance {
        /// Dimensionless squared; target Mandel strain coordinate order.
        strain_mandel: [[f64; 6]; 6],
        /// Pa squared; covariance of C * epsilon_free in target Mandel order.
        free_stress_mandel_pa2: [[f64; 6]; 6],
    },
    Unknown {
        reason: CorrelationUnknownReason,
    },
}

/// An executable nominal thermal state and its conditional strain/stress
/// covariance. This is a linear second-moment transform, not a confidence band
/// or a claim of independently measured, calibrated material behavior.
#[derive(Debug, Clone, PartialEq)]
pub struct JointThermalStrainTransformReceipt {
    nominal: ThermalStrainTransformReceipt,
    uncertainty: ThermalStrainUncertainty,
    source_joint_identity: ContentHash,
    identity: ContentHash,
}

impl JointThermalStrainTransformReceipt {
    #[must_use]
    pub const fn nominal(&self) -> &ThermalStrainTransformReceipt {
        &self.nominal
    }
    #[must_use]
    pub const fn uncertainty(&self) -> &ThermalStrainUncertainty {
        &self.uncertainty
    }
    #[must_use]
    pub const fn source_joint_identity(&self) -> ContentHash {
        self.source_joint_identity
    }
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

// Fixed 6x6 tile, no unbounded allocation or iteration. Source covariance has
// already been admitted by fs-matdb. Pairwise publication retains exact output
// symmetry; no eigenvalue clipping or independence assumption repairs input.
fn transform_covariance6(
    map: &[[f64; 6]; 6],
    source: &[[f64; 6]; 6],
) -> Result<[[f64; 6]; 6], TetElasticError> {
    let mut left = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            for k in 0..6 {
                left[i][j] = map[i][k].mul_add(source[k][j], left[i][j]);
            }
        }
    }
    let mut output = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in i..6 {
            let mut value = 0.0;
            for k in 0..6 {
                value = left[i][k].mul_add(map[j][k], value);
            }
            if !value.is_finite() || (i == j && value < 0.0) {
                return Err(TetElasticError::InvalidThermalStrain {
                    what: "covariance transform is nonfinite or has a negative variance",
                });
            }
            output[i][j] = value;
            output[j][i] = value;
        }
    }
    if source.iter().flatten().any(|v| *v != 0.0) && output.iter().flatten().all(|v| *v == 0.0) {
        return Err(TetElasticError::InvalidThermalStrain {
            what: "covariance transform loses the complete nonzero covariance",
        });
    }
    Ok(output)
}

impl TetThermalStrainState {
    /// Propagate source joint strain statistics through a fixed proper rotation
    /// and fixed admitted stiffness. This computes A Sigma A^T, preserving all
    /// correlations. Unknown source correlation remains unknown. Upstream must
    /// still declare/admit the integrated stress-free temperature path.
    pub fn try_from_joint_strain(
        elastic: &ElasticTensorTransformReceipt,
        strain: &JointStrainTensorStatePoint,
        temperature_k: f64,
        reference_temperature_k: f64,
        source_to_target: [[f64; 3]; 3],
        maximum_absolute_mandel_strain: f64,
    ) -> Result<JointThermalStrainTransformReceipt, TetElasticError> {
        let nominal = Self::try_from_resolved_strain(
            elastic,
            strain.nominal(),
            temperature_k,
            reference_temperature_k,
            source_to_target,
            maximum_absolute_mandel_strain,
        )?;
        let uncertainty = match &strain.joint_receipt().correlation {
            JointCorrelation::Unknown { reason } => {
                ThermalStrainUncertainty::Unknown { reason: *reason }
            }
            JointCorrelation::Covariance { covariance, .. } => {
                let source = core::array::from_fn(|i| {
                    core::array::from_fn(|j| {
                        let hi = i.max(j);
                        let lo = i.min(j);
                        covariance[hi * (hi + 1) / 2 + lo]
                    })
                });
                let basis = strain.nominal().basis();
                let order = match basis.order {
                    ElasticTensorOrder::XxYyZzXyYzZx => [0, 1, 2, 3, 4, 5],
                    ElasticTensorOrder::XxYyZzYzZxXy => [0, 1, 2, 5, 3, 4],
                };
                let shear = match basis.notation {
                    StrainTensorNotation::Tensor => core::f64::consts::SQRT_2,
                    StrainTensorNotation::Engineering => core::f64::consts::FRAC_1_SQRT_2,
                    StrainTensorNotation::Mandel => 1.0,
                };
                let rotation =
                    if source_to_target == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
                        core::array::from_fn(|i| {
                            core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 })
                        })
                    } else {
                        mandel_rotation(source_to_target)
                    };
                let mut map = [[0.0; 6]; 6];
                for row in 0..6 {
                    for column in 0..6 {
                        map[row][order[column]] =
                            rotation[row][column] * if column < 3 { 1.0 } else { shear };
                    }
                }
                let strain_mandel = transform_covariance6(&map, &source)?;
                let free_stress_mandel_pa2 = transform_covariance6(
                    elastic.material().stiffness_mandel_pa(),
                    &strain_mandel,
                )?;
                ThermalStrainUncertainty::Covariance {
                    strain_mandel,
                    free_stress_mandel_pa2,
                }
            }
        };
        let source_joint_identity = strain.identity();
        let mut identity =
            DomainHasher::new("org.frankensim.fs-solid.joint-thermal-strain-transform.v1");
        identity.update(nominal.state().identity().as_bytes());
        identity.update(source_joint_identity.as_bytes());
        match &uncertainty {
            ThermalStrainUncertainty::Unknown { reason } => {
                identity.update(&[0]);
                identity.update(reason.tag().as_bytes());
            }
            ThermalStrainUncertainty::Covariance {
                strain_mandel,
                free_stress_mandel_pa2,
            } => {
                identity.update(&[1]);
                for value in strain_mandel
                    .iter()
                    .flatten()
                    .chain(free_stress_mandel_pa2.iter().flatten())
                {
                    identity.update(&value.to_bits().to_le_bytes());
                }
            }
        }
        Ok(JointThermalStrainTransformReceipt {
            nominal,
            uncertainty,
            source_joint_identity,
            identity: identity.finalize(),
        })
    }

    /// Use all six source-backed strain components and bind their selected
    /// claims/query point to the thermal state. The caller declares these to be
    /// stress-free strain already integrated over the supplied temperature path;
    /// point support alone does not establish that integration or path support.
    pub fn try_from_resolved_strain(
        elastic: &ElasticTensorTransformReceipt,
        strain: &StrainTensorStatePoint,
        temperature_k: f64,
        reference_temperature_k: f64,
        source_to_target: [[f64; 3]; 3],
        maximum_absolute_mandel_strain: f64,
    ) -> Result<ThermalStrainTransformReceipt, TetElasticError> {
        Self::try_from_oriented_strain(
            elastic,
            temperature_k,
            reference_temperature_k,
            *strain.strain(),
            strain.basis(),
            source_to_target,
            strain.resolved().identity(),
            maximum_absolute_mandel_strain,
        )
    }

    /// Rotate an upstream-integrated symmetric strain into the elastic frame.
    ///
    /// This applies the SECOND-order law `epsilon' = Q epsilon Q^T`.
    /// Engineering inputs double only the strain shear coordinates. The
    /// admitted elastic receipt supplies the target frame and material identity;
    /// no caller-provided target label can detach the strain from that material.
    /// `maximum_absolute_mandel_strain` bounds the TARGET Mandel components.
    /// Temperatures and `thermal_law_identity` describe the supplied integrated
    /// path; the upstream law must integrate expansion coefficients and admit
    /// the path against its support. This adapter does not check that support.
    #[allow(clippy::too_many_arguments)]
    pub fn try_from_oriented_strain(
        elastic: &ElasticTensorTransformReceipt,
        temperature_k: f64,
        reference_temperature_k: f64,
        source_free_strain: [f64; 6],
        source_basis: StrainTensorBasis,
        source_to_target: [[f64; 3]; 3],
        thermal_law_identity: ContentHash,
        maximum_absolute_mandel_strain: f64,
    ) -> Result<ThermalStrainTransformReceipt, TetElasticError> {
        validate_rotation(source_to_target)?;
        let identity_rotation = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        if source_basis.frame == ContentHash([0; 32]) {
            return Err(TetElasticError::InvalidThermalStrain {
                what: "source strain frame identity must be nonzero",
            });
        }
        if source_basis.frame == elastic.target_frame() && source_to_target != identity_rotation {
            return Err(TetElasticError::InvalidThermalStrain {
                what: "a strain frame mapped to itself requires the identity rotation",
            });
        }
        if source_free_strain.iter().any(|value| !value.is_finite()) {
            return Err(TetElasticError::InvalidThermalStrain {
                what: "source_free_strain must be finite",
            });
        }
        let order = match source_basis.order {
            ElasticTensorOrder::XxYyZzXyYzZx => [0, 1, 2, 3, 4, 5],
            ElasticTensorOrder::XxYyZzYzZxXy => [0, 1, 2, 5, 3, 4],
        };
        let shear_scale = match source_basis.notation {
            StrainTensorNotation::Tensor => core::f64::consts::SQRT_2,
            StrainTensorNotation::Engineering => core::f64::consts::FRAC_1_SQRT_2,
            StrainTensorNotation::Mandel => 1.0,
        };
        let source_mandel: [f64; 6] = core::array::from_fn(|index| {
            source_free_strain[order[index]] * if index < 3 { 1.0 } else { shear_scale }
        });
        if source_mandel.iter().any(|value| !value.is_finite()) {
            return Err(TetElasticError::InvalidThermalStrain {
                what: "strain shear conversion overflows",
            });
        }
        let target_mandel = if source_to_target == identity_rotation {
            source_mandel
        } else {
            let transform = mandel_rotation(source_to_target);
            core::array::from_fn(|row| {
                (0..6).fold(0.0, |sum, column| {
                    transform[row][column].mul_add(source_mandel[column], sum)
                })
            })
        };
        if source_mandel.iter().any(|value| *value != 0.0)
            && target_mandel.iter().all(|value| *value == 0.0)
        {
            return Err(TetElasticError::InvalidThermalStrain {
                what: "strain rotation underflows the complete nonzero tensor",
            });
        }
        let mut state = Self::try_new(
            temperature_k,
            reference_temperature_k,
            target_mandel,
            elastic.material().material_state_identity,
            thermal_law_identity,
            maximum_absolute_mandel_strain,
        )?;
        let mut identity =
            DomainHasher::new("org.frankensim.fs-solid.thermal-strain-frame-transform.v1");
        // The ordinary state identity already binds temperatures, upstream law,
        // exact elastic state, all target values and the target validity bound.
        identity.update(state.identity().as_bytes());
        identity.update(source_basis.frame.as_bytes());
        identity.update(elastic.target_frame().as_bytes());
        identity.update(&[match source_basis.notation {
            StrainTensorNotation::Tensor => 0,
            StrainTensorNotation::Engineering => 1,
            StrainTensorNotation::Mandel => 2,
        }]);
        identity.update(&[match source_basis.order {
            ElasticTensorOrder::XxYyZzXyYzZx => 0,
            ElasticTensorOrder::XxYyZzYzZxXy => 1,
        }]);
        for value in source_free_strain
            .iter()
            .chain(source_to_target.iter().flatten())
        {
            identity.update(&value.to_bits().to_le_bytes());
        }
        state.identity = identity.finalize();
        Ok(ThermalStrainTransformReceipt {
            source_free_strain,
            source_basis,
            target_frame: elastic.target_frame(),
            source_to_target,
            state,
        })
    }

    /// Admit an upstream-integrated thermal strain tensor and its validity bound.
    pub fn try_new(
        temperature_k: f64,
        reference_temperature_k: f64,
        free_strain_mandel: [f64; 6],
        elastic_material_state_identity: [u8; 32],
        thermal_law_identity: ContentHash,
        maximum_absolute_mandel_strain: f64,
    ) -> Result<Self, TetElasticError> {
        if !(temperature_k.is_finite() && temperature_k > 0.0) {
            return Err(TetElasticError::InvalidThermalStrain {
                what: "temperature_k must be finite and positive",
            });
        }
        if !(reference_temperature_k.is_finite() && reference_temperature_k > 0.0) {
            return Err(TetElasticError::InvalidThermalStrain {
                what: "reference_temperature_k must be finite and positive",
            });
        }
        if !(maximum_absolute_mandel_strain.is_finite() && maximum_absolute_mandel_strain > 0.0) {
            return Err(TetElasticError::InvalidThermalStrain {
                what: "maximum_absolute_mandel_strain must be finite and positive",
            });
        }
        if free_strain_mandel.iter().any(|value| !value.is_finite()) {
            return Err(TetElasticError::InvalidThermalStrain {
                what: "free_strain_mandel must be finite",
            });
        }
        if free_strain_mandel
            .iter()
            .any(|value| value.abs() > maximum_absolute_mandel_strain)
        {
            return Err(TetElasticError::InvalidThermalStrain {
                what: "free_strain_mandel exceeds its declared small-strain validity bound",
            });
        }
        if elastic_material_state_identity == [0; 32] {
            return Err(TetElasticError::InvalidThermalStrain {
                what: "elastic_material_state_identity must not be zero",
            });
        }
        if thermal_law_identity
            .as_bytes()
            .iter()
            .all(|byte| *byte == 0)
        {
            return Err(TetElasticError::InvalidThermalStrain {
                what: "thermal_law_identity must not be zero",
            });
        }
        let mut identity = DomainHasher::new("org.frankensim.fs-solid.tet-thermal-strain.v1");
        identity.update(&temperature_k.to_bits().to_le_bytes());
        identity.update(&reference_temperature_k.to_bits().to_le_bytes());
        for value in free_strain_mandel {
            identity.update(&value.to_bits().to_le_bytes());
        }
        identity.update(&elastic_material_state_identity);
        identity.update(thermal_law_identity.as_bytes());
        identity.update(&maximum_absolute_mandel_strain.to_bits().to_le_bytes());
        Ok(Self {
            temperature_k,
            reference_temperature_k,
            free_strain_mandel,
            elastic_material_state_identity,
            thermal_law_identity,
            maximum_absolute_mandel_strain,
            identity: identity.finalize(),
        })
    }

    /// Bind an evidence-bearing isotropic expansion path to its elastic state.
    ///
    /// The current endpoint must be the exact same card and query point used
    /// for tangent elasticity. This prevents a strain integrated for one
    /// temperature, pressure, phase, or processing state from loading another.
    pub fn try_from_isotropic_expansion(
        elastic: &IsotropicElasticStatePoint,
        expansion: &IntegratedIsotropicThermalExpansion,
        maximum_absolute_mandel_strain: f64,
    ) -> Result<Self, TetElasticError> {
        let current = expansion.current().resolved();
        if current.card_identity() != elastic.resolved().card_identity()
            || current.query_point() != elastic.resolved().query_point()
        {
            return Err(TetElasticError::ThermalElasticStatePointMismatch);
        }
        let temperature_k = current
            .query_point()
            .binary_search_by(|(axis, _)| axis.as_str().cmp("T"))
            .ok()
            .and_then(|index| current.query_point().get(index))
            .map(|(_, temperature)| *temperature)
            .ok_or(TetElasticError::InvalidThermalStrain {
                what: "integrated expansion path has no temperature coordinate",
            })?;
        let reference_temperature_k = expansion
            .reference()
            .resolved()
            .query_point()
            .binary_search_by(|(axis, _)| axis.as_str().cmp("T"))
            .ok()
            .and_then(|index| expansion.reference().resolved().query_point().get(index))
            .map(|(_, temperature)| *temperature)
            .ok_or(TetElasticError::InvalidThermalStrain {
                what: "integrated expansion reference has no temperature coordinate",
            })?;
        let strain = expansion.free_linear_strain();
        Self::try_new(
            temperature_k,
            reference_temperature_k,
            [strain, strain, strain, 0.0, 0.0, 0.0],
            *elastic.resolved().identity().as_bytes(),
            expansion.identity(),
            maximum_absolute_mandel_strain,
        )
    }

    /// Current absolute temperature [K].
    #[must_use]
    pub const fn temperature_k(self) -> f64 {
        self.temperature_k
    }

    /// Stress-free reference temperature [K].
    #[must_use]
    pub const fn reference_temperature_k(self) -> f64 {
        self.reference_temperature_k
    }

    /// Total stress-free thermal strain in Mandel ordering.
    #[must_use]
    pub const fn free_strain_mandel(self) -> [f64; 6] {
        self.free_strain_mandel
    }

    /// Elastic material state to which this thermal strain was resolved.
    #[must_use]
    pub const fn elastic_material_state_identity(self) -> [u8; 32] {
        self.elastic_material_state_identity
    }

    /// Upstream expansion-law or integrated-data identity.
    #[must_use]
    pub const fn thermal_law_identity(self) -> ContentHash {
        self.thermal_law_identity
    }

    /// Declared small-strain applicability ceiling.
    #[must_use]
    pub const fn maximum_absolute_mandel_strain(self) -> f64 {
        self.maximum_absolute_mandel_strain
    }

    /// Complete state identity.
    #[must_use]
    pub const fn identity(self) -> ContentHash {
        self.identity
    }
}

/// Uniform or per-element thermal strain resolved on the elastic mesh.
#[derive(Debug, Clone, Copy)]
pub enum TetThermalStrainField<'a> {
    /// One thermal state for every tetrahedron.
    Uniform(&'a TetThermalStrainState),
    /// One thermal state per tetrahedron, in connectivity order.
    PerElement(&'a [TetThermalStrainState]),
}

impl TetThermalStrainField<'_> {
    fn validate(&self, element_count: usize) -> Result<(), TetElasticError> {
        match self {
            Self::Uniform(_) => Ok(()),
            Self::PerElement(states) if states.len() == element_count => Ok(()),
            Self::PerElement(states) => Err(TetElasticError::ThermalStrainCountMismatch {
                expected: element_count,
                got: states.len(),
            }),
        }
    }

    fn at(&self, element: usize) -> &TetThermalStrainState {
        match self {
            Self::Uniform(state) => state,
            Self::PerElement(states) => &states[element],
        }
    }
}

/// A body-fitted linear-tetrahedron elasticity/mass assembly request.
pub struct TetLinearElasticProblem<'a> {
    /// Vertex coordinates [m].
    pub nodes_m: &'a [[f64; 3]],
    /// Four vertex indices per conforming tetrahedron.
    pub tetrahedra: &'a [[usize; 4]],
    /// Uniform or per-element resolved material field.
    pub materials: TetMaterialField<'a>,
    /// Strongly fixed global displacement DOFs (`3*node + component`).
    pub fixed_dofs: &'a [usize],
    /// Explicit work/quality envelope.
    pub budget: TetAssemblyBudget,
}

/// Canonical reduced matrices and the map back to full displacement DOFs.
#[derive(Debug, Clone)]
pub struct TetElasticAssembly {
    /// Reduced symmetric stiffness matrix [N/m].
    pub stiffness: Csr,
    /// Reduced symmetric consistent mass matrix [kg].
    pub mass: Csr,
    /// `free_dofs[i]` is the full displacement DOF represented by reduced row i.
    pub free_dofs: Vec<usize>,
    /// Total physical mass `sum rho_e V_e` [kg].
    pub total_mass_kg: f64,
    /// Smallest admitted element quality `abs(det(J))/h_max^3`.
    pub minimum_scaled_jacobian: f64,
    /// Positive element volumes [m^3], in connectivity order.
    pub element_volumes_m3: Vec<f64>,
}

/// Equivalent nodal load induced by a stress-free thermal strain field.
///
/// The sign convention is `K u = f_external + equivalent_force`: an
/// unconstrained affine displacement matching a uniform free strain satisfies
/// `K u = equivalent_force`. The full vector retains load contributions at
/// constrained DOFs; actual support reactions follow from the solved stress.
/// `reduced_force_n` follows the elastic assembly's free-DOF order.
#[derive(Debug, Clone, PartialEq)]
pub struct TetThermalLoad {
    /// Full nodal force vector in node-major xyz order [N].
    pub full_equivalent_force_n: Vec<f64>,
    /// Equivalent force on unconstrained DOFs [N].
    pub reduced_force_n: Vec<f64>,
    /// Full DOF represented by each reduced-force entry.
    pub free_dofs: Vec<usize>,
    /// Identity binding mesh, elastic states, thermal states, and constraints.
    pub identity: ContentHash,
}

/// Bounded convergence controls for a static thermal-displacement solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TetStaticSolveConfig {
    /// Required true Euclidean relative residual.
    pub relative_residual_tolerance: f64,
    /// Maximum preconditioned-CG iterations.
    pub maximum_iterations: usize,
}

impl Default for TetStaticSolveConfig {
    fn default() -> Self {
        Self {
            relative_residual_tolerance: 1.0e-10,
            maximum_iterations: 40_000,
        }
    }
}

/// Converged small-strain displacement generated by a thermal strain field.
#[derive(Debug, Clone, PartialEq)]
pub struct TetThermalDisplacementSolution {
    /// Full nodal displacement vectors [m]; strongly fixed DOFs are zero.
    displacement_m: Vec<[f64; 3]>,
    /// Iterations used by the reduced SPD solve.
    iterations: usize,
    /// Independently recomputed true Euclidean relative residual.
    true_relative_residual: f64,
    /// Identity binding thermal load, solve controls, and result.
    identity: ContentHash,
    reference_geometry_identity: ContentHash,
    response_input_identity: ContentHash,
}

impl TetThermalDisplacementSolution {
    /// Full nodal displacement vectors [m].
    #[must_use]
    pub fn displacement_m(&self) -> &[[f64; 3]] {
        &self.displacement_m
    }

    /// Iterations used by the reduced solve.
    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.iterations
    }

    /// Independently recomputed true Euclidean relative residual.
    #[must_use]
    pub const fn true_relative_residual(&self) -> f64 {
        self.true_relative_residual
    }

    /// Complete solution identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Constant P1 element response in the reference mesh's Mandel coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct TetElementThermalStress {
    elastic_strain_mandel: [f64; 6],
    stress_mandel_pa: [f64; 6],
    elastic_energy_j: f64,
}

impl TetElementThermalStress {
    /// Compatible total strain minus the integrated stress-free strain.
    #[must_use]
    pub const fn elastic_strain_mandel(&self) -> &[f64; 6] {
        &self.elastic_strain_mandel
    }
    /// Actual `C * (B u - epsilon_free)`, not the equivalent thermal-load stress.
    #[must_use]
    pub const fn stress_mandel_pa(&self) -> &[f64; 6] {
        &self.stress_mandel_pa
    }
    /// Stored energy `V/2 epsilon_elastic:sigma` for this element [J].
    #[must_use]
    pub const fn elastic_energy_j(&self) -> f64 {
        self.elastic_energy_j
    }
}

/// Physical stress and nodal equilibrium of an accepted thermal solution.
/// The caller declares the mesh frame; this is not a calibrated frame claim.
#[derive(Debug, Clone, PartialEq)]
pub struct TetThermalStressReport {
    elements: Vec<TetElementThermalStress>,
    nodal_internal_force_n: Vec<[f64; 3]>,
    elastic_energy_j: f64,
    mesh_frame: ContentHash,
    solution_identity: ContentHash,
    identity: ContentHash,
}

impl TetThermalStressReport {
    /// Per-element response in connectivity order.
    #[must_use]
    pub fn elements(&self) -> &[TetElementThermalStress] {
        &self.elements
    }
    /// `integral B^T sigma dV`: support forces on fixed DOFs, equilibrium
    /// residual on free DOFs. No residual entries are forcibly zeroed.
    #[must_use]
    pub fn nodal_internal_force_n(&self) -> &[[f64; 3]] {
        &self.nodal_internal_force_n
    }
    /// Sum of element stored elastic energies [J].
    #[must_use]
    pub const fn elastic_energy_j(&self) -> f64 {
        self.elastic_energy_j
    }
    /// Caller-declared coordinate frame of the reference mesh.
    #[must_use]
    pub const fn mesh_frame(&self) -> ContentHash {
        self.mesh_frame
    }
    /// Accepted displacement solution supplying this response.
    #[must_use]
    pub const fn solution_identity(&self) -> ContentHash {
        self.solution_identity
    }
    /// Identity of the complete input-bound recovered response.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Express an element's actual stress as `Q sigma Q^T`. Tensor and
    /// engineering output use identical physical shear values, distinct from
    /// engineering strain. The fixed 6x6 calculation has bounded work.
    pub fn stress_in_frame(
        &self,
        element: usize,
        target_basis: StressTensorBasis,
        mesh_to_target: [[f64; 3]; 3],
    ) -> Result<StressTensorTransformReceipt, TetElasticError> {
        let source = self
            .elements
            .get(element)
            .ok_or(TetElasticError::InvalidStressObservation {
                what: "stress element index is outside the recovered field",
            })?
            .stress_mandel_pa;
        if target_basis.frame == ContentHash([0; 32]) {
            return Err(TetElasticError::InvalidStressObservation {
                what: "target stress frame must not be zero",
            });
        }
        validate_rotation(mesh_to_target)?;
        let identity_rotation =
            mesh_to_target == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        if target_basis.frame == self.mesh_frame && !identity_rotation {
            return Err(TetElasticError::InvalidStressObservation {
                what: "a stress frame mapped to itself requires identity rotation",
            });
        }
        let rotated = if identity_rotation {
            source
        } else {
            let map = mandel_rotation(mesh_to_target);
            core::array::from_fn(|i| (0..6).fold(0.0, |sum, j| map[i][j].mul_add(source[j], sum)))
        };
        let order = match target_basis.order {
            ElasticTensorOrder::XxYyZzXyYzZx => [0, 1, 2, 3, 4, 5],
            ElasticTensorOrder::XxYyZzYzZxXy => [0, 1, 2, 4, 5, 3],
        };
        let shear = match target_basis.notation {
            StressTensorNotation::Tensor | StressTensorNotation::Engineering => {
                core::f64::consts::FRAC_1_SQRT_2
            }
            StressTensorNotation::Mandel => 1.0,
        };
        let stress_pa: [f64; 6] =
            core::array::from_fn(|i| rotated[order[i]] * if i < 3 { 1.0 } else { shear });
        if stress_pa.iter().any(|v| !v.is_finite())
            || (source.iter().any(|v| *v != 0.0) && stress_pa.iter().all(|v| *v == 0.0))
        {
            return Err(TetElasticError::InvalidStressObservation {
                what: "stress rotation is nonfinite or loses the complete nonzero stress",
            });
        }
        let mut identity = DomainHasher::new("org.frankensim.fs-solid.stress-frame-transform.v1");
        identity.update(self.identity.as_bytes());
        hash_usize(&mut identity, element)?;
        identity.update(target_basis.frame.as_bytes());
        identity.update(&[
            match target_basis.notation {
                StressTensorNotation::Tensor => 0,
                StressTensorNotation::Engineering => 1,
                StressTensorNotation::Mandel => 2,
            },
            match target_basis.order {
                ElasticTensorOrder::XxYyZzXyYzZx => 0,
                ElasticTensorOrder::XxYyZzYzZxXy => 1,
            },
        ]);
        for value in mesh_to_target.iter().flatten().chain(stress_pa.iter()) {
            identity.update(&value.to_bits().to_le_bytes());
        }
        Ok(StressTensorTransformReceipt {
            source_report_identity: self.identity,
            element,
            source_frame: self.mesh_frame,
            source_stress_mandel_pa: source,
            target_basis,
            source_to_target: mesh_to_target,
            stress_pa,
            identity: identity.finalize(),
        })
    }
}

/// Exact source stress, output coordinates and the second-order transform.
#[derive(Debug, Clone, PartialEq)]
pub struct StressTensorTransformReceipt {
    source_report_identity: ContentHash,
    element: usize,
    source_frame: ContentHash,
    source_stress_mandel_pa: [f64; 6],
    target_basis: StressTensorBasis,
    source_to_target: [[f64; 3]; 3],
    stress_pa: [f64; 6],
    identity: ContentHash,
}

impl StressTensorTransformReceipt {
    /// Identity of the recovered stress field supplying this element.
    #[must_use]
    pub const fn source_report_identity(&self) -> ContentHash {
        self.source_report_identity
    }
    /// Element index in the source report's connectivity order.
    #[must_use]
    pub const fn element(&self) -> usize {
        self.element
    }
    /// Declared reference-mesh coordinate frame.
    #[must_use]
    pub const fn source_frame(&self) -> ContentHash {
        self.source_frame
    }
    /// Source stress before rotation, in canonical Mandel coordinates.
    #[must_use]
    pub const fn source_stress_mandel_pa(&self) -> &[f64; 6] {
        &self.source_stress_mandel_pa
    }
    /// Output stress convention, component order and coordinate frame.
    #[must_use]
    pub const fn target_basis(&self) -> StressTensorBasis {
        self.target_basis
    }
    /// Proper rotation; columns are source axes in the target frame.
    #[must_use]
    pub const fn source_to_target(&self) -> &[[f64; 3]; 3] {
        &self.source_to_target
    }
    /// Transformed stress [Pa] in the declared target coordinates.
    #[must_use]
    pub const fn stress_pa(&self) -> &[f64; 6] {
        &self.stress_pa
    }
    /// Identity binding the source element, transform and output stress.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Updated fixed-topology tetrahedral geometry after an admitted displacement.
#[derive(Debug, Clone, PartialEq)]
pub struct TetDeformedMesh {
    nodes_m: Vec<[f64; 3]>,
    tetrahedra: Vec<[usize; 4]>,
    minimum_scaled_jacobian: f64,
    identity: ContentHash,
}

impl TetDeformedMesh {
    /// Updated vertex coordinates [m].
    #[must_use]
    pub fn nodes_m(&self) -> &[[f64; 3]] {
        &self.nodes_m
    }

    /// Unchanged fixed-topology connectivity.
    #[must_use]
    pub fn tetrahedra(&self) -> &[[usize; 4]] {
        &self.tetrahedra
    }

    /// Smallest scaled Jacobian after the update.
    #[must_use]
    pub const fn minimum_scaled_jacobian(&self) -> f64 {
        self.minimum_scaled_jacobian
    }

    /// Identity binding reference geometry, displacement, and updated chart.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Typed, stable refusal surface for 3-D elastic assembly.
#[derive(Debug, Clone, PartialEq)]
pub enum TetElasticError {
    /// Empty node or element set.
    EmptyMesh,
    /// A coordinate is NaN or infinite.
    NonFiniteNode {
        /// Zero-based node index.
        node: usize,
    },
    /// A connectivity entry is outside the node array.
    NodeIndexOutOfRange {
        /// Zero-based tetrahedron index.
        element: usize,
        /// Invalid node index.
        node: usize,
    },
    /// A tetrahedron repeats one of its four vertices.
    RepeatedTetNode {
        /// Zero-based tetrahedron index.
        element: usize,
    },
    /// A mesh node belongs to no tetrahedron.
    UnreferencedNode {
        /// Zero-based node index.
        node: usize,
    },
    /// A tetrahedron is singular or below the declared quality floor.
    DegenerateTet {
        /// Zero-based tetrahedron index.
        element: usize,
        /// Measured dimensionless quality.
        scaled_jacobian: f64,
        /// Declared quality floor.
        minimum: f64,
    },
    /// A material scalar or identity is inadmissible.
    InvalidMaterial {
        /// Failed material invariant.
        what: &'static str,
    },
    /// A per-element material array has the wrong length.
    MaterialCountMismatch {
        /// Required number of cards.
        expected: usize,
        /// Supplied number of cards.
        got: usize,
    },
    /// A per-element thermal-state array has the wrong length.
    ThermalStrainCountMismatch {
        /// Required number of states.
        expected: usize,
        /// Supplied number of states.
        got: usize,
    },
    /// A thermal strain state or its authority was inadmissible.
    InvalidThermalStrain {
        /// Failed invariant.
        what: &'static str,
    },
    /// Thermal strain and elastic stiffness were resolved for different states.
    ThermalElasticStateMismatch {
        /// Zero-based tetrahedron index.
        element: usize,
    },
    /// Expansion path and elastic tangent were resolved at different state points.
    ThermalElasticStatePointMismatch,
    /// Stress recovery controls, transformed values or applicability failed.
    InvalidStressObservation {
        /// Failed recovery, frame or numerical invariant.
        what: &'static str,
    },
    /// Recovery was requested with a different mesh/material/thermal/load state.
    StressRecoveryInputMismatch,
    /// Finite inputs produced a non-finite stress or nodal force.
    NonFiniteThermalLoad {
        /// Zero-based tetrahedron index.
        element: usize,
        /// Local displacement DOF when force evaluation failed, if applicable.
        local_dof: Option<usize>,
    },
    /// Static solve controls were malformed.
    InvalidStaticSolveConfig {
        /// Failed invariant.
        what: &'static str,
    },
    /// The static displacement solve missed its true-residual gate.
    StaticSolveFailed {
        /// Iterations performed.
        iterations: usize,
        /// Independently recomputed true relative residual.
        true_relative_residual: f64,
    },
    /// Dirichlet constraints do not remove all six rigid modes of a component.
    UnconstrainedRigidModes {
        /// Canonical smallest node of the connected component.
        component: usize,
        /// Rank of the six-column rigid-mode constraint matrix.
        constrained_rank: usize,
    },
    /// Displacement and reference mesh identities do not match.
    GeometryUpdateIdentityMismatch,
    /// Geometry-update controls or displacements were inadmissible.
    InvalidGeometryUpdate {
        /// Failed invariant.
        what: &'static str,
    },
    /// An accepted displacement would invert a tetrahedron.
    InvertedTet {
        /// Zero-based tetrahedron index.
        element: usize,
    },
    /// A fixed displacement DOF is outside the full vector.
    FixedDofOutOfRange {
        /// Invalid full displacement DOF.
        dof: usize,
        /// Total full displacement DOFs.
        full_dofs: usize,
    },
    /// A fixed displacement DOF was declared more than once.
    DuplicateFixedDof {
        /// Repeated full displacement DOF.
        dof: usize,
    },
    /// Matrix/modal assembly has no free displacement DOFs.
    NoFreeDofs,
    /// A checked integer size operation overflowed.
    SizeOverflow,
    /// A declared work or quality envelope was exceeded.
    BudgetExceeded {
        /// Bounded resource.
        what: &'static str,
        /// Requested count.
        requested: usize,
        /// Admitted count.
        maximum: usize,
    },
    /// Cancellation was observed at an element boundary.
    Cancelled,
}

impl core::fmt::Display for TetElasticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyMesh => write!(f, "FS-SOLID-TET-EMPTY-MESH"),
            Self::NonFiniteNode { node } => {
                write!(f, "FS-SOLID-TET-NONFINITE-NODE: node {node}")
            }
            Self::NodeIndexOutOfRange { element, node } => write!(
                f,
                "FS-SOLID-TET-NODE-RANGE: element {element} names node {node}"
            ),
            Self::RepeatedTetNode { element } => {
                write!(f, "FS-SOLID-TET-REPEATED-NODE: element {element}")
            }
            Self::UnreferencedNode { node } => {
                write!(f, "FS-SOLID-TET-UNREFERENCED-NODE: node {node}")
            }
            Self::DegenerateTet {
                element,
                scaled_jacobian,
                minimum,
            } => write!(
                f,
                "FS-SOLID-TET-DEGENERATE: element {element} scaled Jacobian \
                 {scaled_jacobian:.3e} below {minimum:.3e}"
            ),
            Self::InvalidMaterial { what } => {
                write!(f, "FS-SOLID-TET-MATERIAL: {what}")
            }
            Self::MaterialCountMismatch { expected, got } => write!(
                f,
                "FS-SOLID-TET-MATERIAL-COUNT: expected {expected}, got {got}"
            ),
            Self::ThermalStrainCountMismatch { expected, got } => write!(
                f,
                "FS-SOLID-TET-THERMAL-COUNT: expected {expected}, got {got}"
            ),
            Self::InvalidThermalStrain { what } => {
                write!(f, "FS-SOLID-TET-THERMAL-STRAIN: {what}")
            }
            Self::ThermalElasticStateMismatch { element } => write!(
                f,
                "FS-SOLID-TET-THERMAL-ELASTIC-STATE: element {element} resolves different states"
            ),
            Self::ThermalElasticStatePointMismatch => write!(
                f,
                "FS-SOLID-TET-THERMAL-ELASTIC-POINT: expansion and elasticity resolve different state points"
            ),
            Self::InvalidStressObservation { what } => write!(f, "FS-SOLID-TET-STRESS: {what}"),
            Self::StressRecoveryInputMismatch => write!(
                f,
                "FS-SOLID-TET-STRESS-INPUT: solution belongs to different physical inputs"
            ),
            Self::NonFiniteThermalLoad { element, local_dof } => write!(
                f,
                "FS-SOLID-TET-THERMAL-NONFINITE: element {element}, local dof {local_dof:?}"
            ),
            Self::InvalidStaticSolveConfig { what } => {
                write!(f, "FS-SOLID-TET-STATIC-CONFIG: {what}")
            }
            Self::StaticSolveFailed {
                iterations,
                true_relative_residual,
            } => write!(
                f,
                "FS-SOLID-TET-STATIC-SOLVE: residual {true_relative_residual:.6e} after {iterations} iterations"
            ),
            Self::UnconstrainedRigidModes {
                component,
                constrained_rank,
            } => write!(
                f,
                "FS-SOLID-TET-RIGID-MODES: component {component} constrains rank {constrained_rank} of 6 rigid modes"
            ),
            Self::GeometryUpdateIdentityMismatch => write!(
                f,
                "FS-SOLID-TET-GEOMETRY-IDENTITY: displacement belongs to another reference mesh"
            ),
            Self::InvalidGeometryUpdate { what } => {
                write!(f, "FS-SOLID-TET-GEOMETRY-UPDATE: {what}")
            }
            Self::InvertedTet { element } => {
                write!(f, "FS-SOLID-TET-INVERTED: element {element}")
            }
            Self::FixedDofOutOfRange { dof, full_dofs } => write!(
                f,
                "FS-SOLID-TET-FIXED-DOF-RANGE: {dof} outside 0..{full_dofs}"
            ),
            Self::DuplicateFixedDof { dof } => {
                write!(f, "FS-SOLID-TET-DUPLICATE-FIXED-DOF: {dof}")
            }
            Self::NoFreeDofs => write!(f, "FS-SOLID-TET-NO-FREE-DOFS"),
            Self::SizeOverflow => write!(f, "FS-SOLID-TET-SIZE-OVERFLOW"),
            Self::BudgetExceeded {
                what,
                requested,
                maximum,
            } => write!(
                f,
                "FS-SOLID-TET-BUDGET: {what} requested {requested}, maximum {maximum}"
            ),
            Self::Cancelled => write!(f, "FS-SOLID-TET-CANCELLED"),
        }
    }
}

impl std::error::Error for TetElasticError {}

impl TetLinearElasticProblem<'_> {
    /// Assemble reduced stiffness and consistent mass matrices.
    ///
    /// Cancellation is polled before allocation and at every element
    /// boundary. No partial matrix is published on refusal.
    pub fn assemble(&self, cx: &Cx<'_>) -> Result<TetElasticAssembly, TetElasticError> {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        self.validate_budget()?;
        self.materials.validate(self.tetrahedra.len())?;

        let (free_dofs, reduced_of) = free_dof_map(
            self.nodes_m.len(),
            self.fixed_dofs,
            self.budget.maximum_free_dofs,
        )?;
        if free_dofs.is_empty() {
            return Err(TetElasticError::NoFreeDofs);
        }

        let mut incident = vec![false; self.nodes_m.len()];
        let n = free_dofs.len();
        let mut stiffness = Coo::new(n, n);
        let mut mass = Coo::new(n, n);
        let mut total_mass_kg = 0.0;
        let mut minimum_scaled_jacobian = f64::INFINITY;
        let mut element_volumes_m3 = Vec::with_capacity(self.tetrahedra.len());

        for (element, tet) in self.tetrahedra.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            validate_connectivity(element, tet, self.nodes_m.len())?;
            for &node in tet {
                incident[node] = true;
            }
            let points = tet.map(|node| self.nodes_m[node]);
            let geometry =
                TetGeometry::try_new(element, points, self.budget.minimum_scaled_jacobian)?;
            minimum_scaled_jacobian = minimum_scaled_jacobian.min(geometry.scaled_jacobian);
            element_volumes_m3.push(geometry.volume);

            let material = self.materials.at(element);
            let (ke, me) = element_matrices(&geometry, material);
            total_mass_kg += material.density_kg_m3 * geometry.volume;
            for local_row in 0..12 {
                let full_row = 3 * tet[local_row / 3] + local_row % 3;
                let Some(row) = reduced_of[full_row] else {
                    continue;
                };
                for local_col in 0..12 {
                    let full_col = 3 * tet[local_col / 3] + local_col % 3;
                    let Some(col) = reduced_of[full_col] else {
                        continue;
                    };
                    let stiffness_value = ke[local_row][local_col];
                    if stiffness_value != 0.0 {
                        stiffness.push(row, col, stiffness_value);
                    }
                    let mass_value = me[local_row][local_col];
                    if mass_value != 0.0 {
                        mass.push(row, col, mass_value);
                    }
                }
            }
        }
        if let Some(node) = incident.iter().position(|seen| !seen) {
            return Err(TetElasticError::UnreferencedNode { node });
        }
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok(TetElasticAssembly {
            stiffness: stiffness.assemble(),
            mass: mass.assemble(),
            free_dofs,
            total_mass_kg,
            minimum_scaled_jacobian,
            element_volumes_m3,
        })
    }

    /// Assemble the equivalent force generated by an admitted thermal strain field.
    ///
    /// This is the generic small-strain thermomechanical bridge. It accepts an
    /// arbitrary integrated symmetric thermal strain per element and uses the
    /// same element geometry, elastic tensor, constraints, cancellation, and
    /// quality budget as [`Self::assemble`]. It does not solve the displacement
    /// field or update mesh coordinates; those remain explicit downstream acts.
    pub fn assemble_thermal_load(
        &self,
        thermal_strains: TetThermalStrainField<'_>,
        cx: &Cx<'_>,
    ) -> Result<TetThermalLoad, TetElasticError> {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        self.validate_budget()?;
        self.materials.validate(self.tetrahedra.len())?;
        thermal_strains.validate(self.tetrahedra.len())?;
        let (free_dofs, _) = free_dof_map(
            self.nodes_m.len(),
            self.fixed_dofs,
            self.budget.maximum_free_dofs,
        )?;
        let full_dofs = self
            .nodes_m
            .len()
            .checked_mul(3)
            .ok_or(TetElasticError::SizeOverflow)?;
        let mut full_equivalent_force_n = vec![0.0; full_dofs];
        let mut incident = vec![false; self.nodes_m.len()];
        let mut identity = DomainHasher::new("org.frankensim.fs-solid.tet-thermal-load.v1");
        hash_usize(&mut identity, self.nodes_m.len())?;
        for point in self.nodes_m {
            for coordinate in point {
                identity.update(&coordinate.to_bits().to_le_bytes());
            }
        }
        hash_usize(&mut identity, self.tetrahedra.len())?;
        hash_usize(&mut identity, free_dofs.len())?;
        for &dof in &free_dofs {
            hash_usize(&mut identity, dof)?;
        }

        for (element, tet) in self.tetrahedra.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            validate_connectivity(element, tet, self.nodes_m.len())?;
            for &node in tet {
                incident[node] = true;
                hash_usize(&mut identity, node)?;
            }
            let geometry = TetGeometry::try_new(
                element,
                tet.map(|node| self.nodes_m[node]),
                self.budget.minimum_scaled_jacobian,
            )?;
            let material = self.materials.at(element);
            let thermal = thermal_strains.at(element);
            if thermal.elastic_material_state_identity != material.material_state_identity {
                return Err(TetElasticError::ThermalElasticStateMismatch { element });
            }
            identity.update(&material.material_state_identity);
            identity.update(thermal.identity.as_bytes());

            let b = strain_displacement_matrix(&geometry);
            let mut free_stress_mandel_pa = [0.0; 6];
            for (row, stress) in free_stress_mandel_pa.iter_mut().enumerate() {
                for column in 0..6 {
                    *stress = material.stiffness_mandel_pa[row][column]
                        .mul_add(thermal.free_strain_mandel[column], *stress);
                }
            }
            if free_stress_mandel_pa
                .iter()
                .any(|stress| !stress.is_finite())
            {
                return Err(TetElasticError::NonFiniteThermalLoad {
                    element,
                    local_dof: None,
                });
            }
            for local_dof in 0..12 {
                let mut force_n = 0.0;
                for component in 0..6 {
                    force_n =
                        b[component][local_dof].mul_add(free_stress_mandel_pa[component], force_n);
                }
                force_n *= geometry.volume;
                if !force_n.is_finite() {
                    return Err(TetElasticError::NonFiniteThermalLoad {
                        element,
                        local_dof: Some(local_dof),
                    });
                }
                let full_dof = 3 * tet[local_dof / 3] + local_dof % 3;
                full_equivalent_force_n[full_dof] += force_n;
                if !full_equivalent_force_n[full_dof].is_finite() {
                    return Err(TetElasticError::NonFiniteThermalLoad {
                        element,
                        local_dof: Some(local_dof),
                    });
                }
            }
        }
        if let Some(node) = incident.iter().position(|seen| !seen) {
            return Err(TetElasticError::UnreferencedNode { node });
        }
        let reduced_force_n = free_dofs
            .iter()
            .map(|&full_dof| full_equivalent_force_n[full_dof])
            .collect();
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok(TetThermalLoad {
            full_equivalent_force_n,
            reduced_force_n,
            free_dofs,
            identity: identity.finalize(),
        })
    }

    /// Assemble and solve the fixed-topology small-strain thermal response.
    ///
    /// PCG is advanced in bounded chunks so cancellation is observed during a
    /// long solve. Acceptance uses a freshly recomputed `||f-Ku||/||f||`, not
    /// the recursively updated Krylov estimate. Geometry update or remeshing is
    /// deliberately a separate admitted operation.
    /// A fully fixed mesh has the exact prescribed zero displacement and an
    /// empty reduced system: zero iterations/residual, with full thermal forces
    /// retained for subsequent stress and support-reaction recovery.
    pub fn solve_thermal_displacement(
        &self,
        thermal_strains: TetThermalStrainField<'_>,
        config: TetStaticSolveConfig,
        cx: &Cx<'_>,
    ) -> Result<TetThermalDisplacementSolution, TetElasticError> {
        if !(config.relative_residual_tolerance.is_finite()
            && config.relative_residual_tolerance > 0.0
            && config.relative_residual_tolerance < 1.0)
        {
            return Err(TetElasticError::InvalidStaticSolveConfig {
                what: "relative_residual_tolerance must lie in (0,1)",
            });
        }
        if config.maximum_iterations == 0 {
            return Err(TetElasticError::InvalidStaticSolveConfig {
                what: "maximum_iterations must be positive",
            });
        }
        let load = self.assemble_thermal_load(thermal_strains, cx)?;
        let (reduced_displacement, iterations, true_relative_residual) =
            if load.free_dofs.is_empty() {
                // Every displacement is prescribed zero. The empty reduced system
                // is solved exactly; nonzero full forces remain support reactions.
                (Vec::new(), 0, 0.0)
            } else {
                let assembly = self.assemble(cx)?;
                validate_rigid_mode_constraints(self.nodes_m, self.tetrahedra, self.fixed_dofs)?;
                let preconditioner = Jacobi::new(&assembly.stiffness);
                let operator = CsrOp::symmetric(assembly.stiffness);
                let mut state = CgState::new(&operator, &preconditioner, &load.reduced_force_n);
                while state.iters < config.maximum_iterations
                    && state.rel_residual() >= config.relative_residual_tolerance
                {
                    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
                    let remaining = config.maximum_iterations - state.iters;
                    let _ = state.run(
                        &operator,
                        &preconditioner,
                        config.relative_residual_tolerance,
                        remaining.min(64),
                    );
                }
                let mut applied = vec![0.0; state.x.len()];
                operator.apply(&state.x, &mut applied);
                let residual_norm = euclidean_norm(
                    load.reduced_force_n
                        .iter()
                        .zip(&applied)
                        .map(|(force, ku)| force - ku),
                );
                let force_norm = euclidean_norm(load.reduced_force_n.iter().copied());
                let true_relative_residual = residual_norm / force_norm.max(f64::MIN_POSITIVE);
                if !true_relative_residual.is_finite()
                    || true_relative_residual > config.relative_residual_tolerance
                {
                    return Err(TetElasticError::StaticSolveFailed {
                        iterations: state.iters,
                        true_relative_residual,
                    });
                }
                (state.x, state.iters, true_relative_residual)
            };
        let mut flat = vec![0.0; self.nodes_m.len() * 3];
        for (&full_dof, &value) in load.free_dofs.iter().zip(&reduced_displacement) {
            flat[full_dof] = value;
        }
        let displacement_m = flat
            .chunks_exact(3)
            .map(|value| [value[0], value[1], value[2]])
            .collect::<Vec<_>>();
        let mut identity = DomainHasher::new("org.frankensim.fs-solid.tet-thermal-solve.v1");
        identity.update(load.identity.as_bytes());
        identity.update(&config.relative_residual_tolerance.to_bits().to_le_bytes());
        hash_usize(&mut identity, config.maximum_iterations)?;
        hash_usize(&mut identity, iterations)?;
        identity.update(&true_relative_residual.to_bits().to_le_bytes());
        for displacement in &displacement_m {
            for value in displacement {
                identity.update(&value.to_bits().to_le_bytes());
            }
        }
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok(TetThermalDisplacementSolution {
            displacement_m,
            iterations,
            true_relative_residual,
            identity: identity.finalize(),
            reference_geometry_identity: tet_geometry_identity(self.nodes_m, self.tetrahedra)?,
            response_input_identity: self.thermal_response_input_identity(load.identity, cx)?,
        })
    }

    // A caller-supplied material identity alone cannot detect changed numerical
    // stiffness. Retain the actual law consumed by the solve, without changing
    // the existing v1 thermal-load or displacement identity formats.
    fn thermal_response_input_identity(
        &self,
        load_identity: ContentHash,
        cx: &Cx<'_>,
    ) -> Result<ContentHash, TetElasticError> {
        let mut identity = DomainHasher::new("org.frankensim.fs-solid.thermal-response-input.v1");
        identity.update(load_identity.as_bytes());
        for element in 0..self.tetrahedra.len() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            let material = self.materials.at(element);
            identity.update(&material.density_kg_m3.to_bits().to_le_bytes());
            for value in material.stiffness_mandel_pa.iter().flatten() {
                identity.update(&value.to_bits().to_le_bytes());
            }
        }
        Ok(identity.finalize())
    }

    /// Recover actual small-strain Cauchy stress, stored elastic energy and
    /// nodal internal force from an accepted thermal solve. The caller names the
    /// mesh coordinate frame and bounds both total and elastic Mandel strain.
    /// No extrapolation to yield, finite strain or a confidence bound is implied.
    /// Cancellation is polled before allocation and at every element; refusal
    /// never publishes a partial field. The mesh/material/thermal/constraint
    /// input must match the solved state, including the actual stiffness bits.
    pub fn recover_thermal_stress(
        &self,
        solution: &TetThermalDisplacementSolution,
        thermal_strains: TetThermalStrainField<'_>,
        mesh_frame: ContentHash,
        maximum_absolute_mandel_strain: f64,
        cx: &Cx<'_>,
    ) -> Result<TetThermalStressReport, TetElasticError> {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        if mesh_frame == ContentHash([0; 32])
            || !(maximum_absolute_mandel_strain.is_finite() && maximum_absolute_mandel_strain > 0.0)
        {
            return Err(TetElasticError::InvalidStressObservation {
                what: "mesh frame must be nonzero and strain ceiling must be finite and positive",
            });
        }
        let load = self.assemble_thermal_load(thermal_strains, cx)?;
        let inputs = self.thermal_response_input_identity(load.identity, cx)?;
        if inputs != solution.response_input_identity
            || solution.displacement_m.len() != self.nodes_m.len()
        {
            return Err(TetElasticError::StressRecoveryInputMismatch);
        }
        let mut elements = Vec::with_capacity(self.tetrahedra.len());
        let mut nodal_internal_force_n = vec![[0.0; 3]; self.nodes_m.len()];
        let mut elastic_energy_j = 0.0;
        let mut identity = DomainHasher::new("org.frankensim.fs-solid.thermal-stress-recovery.v1");
        identity.update(inputs.as_bytes());
        identity.update(solution.identity.as_bytes());
        identity.update(mesh_frame.as_bytes());
        identity.update(&maximum_absolute_mandel_strain.to_bits().to_le_bytes());
        for (element, tet) in self.tetrahedra.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            let geometry = TetGeometry::try_new(
                element,
                tet.map(|node| self.nodes_m[node]),
                self.budget.minimum_scaled_jacobian,
            )?;
            let b = strain_displacement_matrix(&geometry);
            let total: [f64; 6] = core::array::from_fn(|i| {
                (0..12).fold(0.0, |sum, d| {
                    b[i][d].mul_add(solution.displacement_m[tet[d / 3]][d % 3], sum)
                })
            });
            let elastic_strain_mandel: [f64; 6] = core::array::from_fn(|i| {
                total[i] - thermal_strains.at(element).free_strain_mandel[i]
            });
            if total
                .iter()
                .chain(&elastic_strain_mandel)
                .any(|v| !v.is_finite() || v.abs() > maximum_absolute_mandel_strain)
            {
                return Err(TetElasticError::InvalidStressObservation {
                    what: "recovered total or elastic strain exceeds the admitted finite small-strain ceiling",
                });
            }
            let material = self.materials.at(element);
            let stress_mandel_pa: [f64; 6] = core::array::from_fn(|i| {
                (0..6).fold(0.0, |sum, j| {
                    material.stiffness_mandel_pa[i][j].mul_add(elastic_strain_mandel[j], sum)
                })
            });
            let energy = (0..6).fold(0.0, |sum, i| {
                elastic_strain_mandel[i].mul_add(stress_mandel_pa[i], sum)
            }) * (0.5 * geometry.volume);
            elastic_energy_j += energy;
            if stress_mandel_pa.iter().any(|v| !v.is_finite())
                || !elastic_energy_j.is_finite()
                || !energy.is_finite()
                || energy < 0.0
            {
                return Err(TetElasticError::InvalidStressObservation {
                    what: "recovered stress or stored elastic energy is nonfinite or energy is negative",
                });
            }
            for d in 0..12 {
                let force = (0..6).fold(0.0, |sum, i| b[i][d].mul_add(stress_mandel_pa[i], sum))
                    * geometry.volume;
                let accumulated = &mut nodal_internal_force_n[tet[d / 3]][d % 3];
                *accumulated += force;
                if !accumulated.is_finite() {
                    return Err(TetElasticError::InvalidStressObservation {
                        what: "recovered nodal internal force is nonfinite",
                    });
                }
            }
            for value in elastic_strain_mandel
                .iter()
                .chain(&stress_mandel_pa)
                .chain(core::iter::once(&energy))
            {
                identity.update(&value.to_bits().to_le_bytes());
            }
            elements.push(TetElementThermalStress {
                elastic_strain_mandel,
                stress_mandel_pa,
                elastic_energy_j: energy,
            });
        }
        for (node, force) in nodal_internal_force_n.iter().enumerate() {
            if node % 4096 == 0 {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            }
            for value in force {
                identity.update(&value.to_bits().to_le_bytes());
            }
        }
        identity.update(&elastic_energy_j.to_bits().to_le_bytes());
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok(TetThermalStressReport {
            elements,
            nodal_internal_force_n,
            elastic_energy_j,
            mesh_frame,
            solution_identity: solution.identity,
            identity: identity.finalize(),
        })
    }

    /// Apply an admitted displacement to the shared fixed-topology mesh.
    ///
    /// This is a transactional geometry update: it refuses foreign solutions,
    /// non-finite or excessive nodal motion, degenerate elements, and every
    /// orientation flip before publishing any updated chart.
    pub fn update_geometry_from_displacement(
        &self,
        solution: &TetThermalDisplacementSolution,
        maximum_nodal_displacement_m: f64,
        cx: &Cx<'_>,
    ) -> Result<TetDeformedMesh, TetElasticError> {
        if !(maximum_nodal_displacement_m.is_finite() && maximum_nodal_displacement_m > 0.0) {
            return Err(TetElasticError::InvalidGeometryUpdate {
                what: "maximum_nodal_displacement_m must be finite and positive",
            });
        }
        if solution.reference_geometry_identity
            != tet_geometry_identity(self.nodes_m, self.tetrahedra)?
        {
            return Err(TetElasticError::GeometryUpdateIdentityMismatch);
        }
        if solution.displacement_m.len() != self.nodes_m.len() {
            return Err(TetElasticError::InvalidGeometryUpdate {
                what: "displacement node count differs from reference mesh",
            });
        }
        let mut nodes_m = Vec::with_capacity(self.nodes_m.len());
        for (node, (point, displacement)) in self
            .nodes_m
            .iter()
            .zip(&solution.displacement_m)
            .enumerate()
        {
            if node % 4096 == 0 {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            }
            let magnitude = euclidean_norm(displacement.iter().copied());
            if !magnitude.is_finite() || magnitude > maximum_nodal_displacement_m {
                return Err(TetElasticError::InvalidGeometryUpdate {
                    what: "nodal displacement exceeds its admitted bound",
                });
            }
            let updated = [
                point[0] + displacement[0],
                point[1] + displacement[1],
                point[2] + displacement[2],
            ];
            if updated.iter().any(|coordinate| !coordinate.is_finite()) {
                return Err(TetElasticError::InvalidGeometryUpdate {
                    what: "updated coordinate is non-finite",
                });
            }
            nodes_m.push(updated);
        }
        let mut minimum_scaled_jacobian = f64::INFINITY;
        for (element, tet) in self.tetrahedra.iter().enumerate() {
            if element % 4096 == 0 {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            }
            let reference = tet.map(|node| self.nodes_m[node]);
            let updated = tet.map(|node| nodes_m[node]);
            if signed_tet_determinant(reference).is_sign_positive()
                != signed_tet_determinant(updated).is_sign_positive()
            {
                return Err(TetElasticError::InvertedTet { element });
            }
            let geometry =
                TetGeometry::try_new(element, updated, self.budget.minimum_scaled_jacobian)?;
            minimum_scaled_jacobian = minimum_scaled_jacobian.min(geometry.scaled_jacobian);
        }
        let mut identity = DomainHasher::new("org.frankensim.fs-solid.tet-deformed-mesh.v1");
        identity.update(solution.reference_geometry_identity.as_bytes());
        identity.update(solution.identity.as_bytes());
        identity.update(&maximum_nodal_displacement_m.to_bits().to_le_bytes());
        for point in &nodes_m {
            for coordinate in point {
                identity.update(&coordinate.to_bits().to_le_bytes());
            }
        }
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok(TetDeformedMesh {
            nodes_m,
            tetrahedra: self.tetrahedra.to_vec(),
            minimum_scaled_jacobian,
            identity: identity.finalize(),
        })
    }

    fn validate_budget(&self) -> Result<(), TetElasticError> {
        if self.nodes_m.is_empty() || self.tetrahedra.is_empty() {
            return Err(TetElasticError::EmptyMesh);
        }
        if self.nodes_m.len() > self.budget.maximum_nodes {
            return Err(TetElasticError::BudgetExceeded {
                what: "nodes",
                requested: self.nodes_m.len(),
                maximum: self.budget.maximum_nodes,
            });
        }
        if self.tetrahedra.len() > self.budget.maximum_tetrahedra {
            return Err(TetElasticError::BudgetExceeded {
                what: "tetrahedra",
                requested: self.tetrahedra.len(),
                maximum: self.budget.maximum_tetrahedra,
            });
        }
        if !(self.budget.minimum_scaled_jacobian.is_finite()
            && self.budget.minimum_scaled_jacobian > 0.0)
        {
            return Err(TetElasticError::DegenerateTet {
                element: 0,
                scaled_jacobian: f64::NAN,
                minimum: self.budget.minimum_scaled_jacobian,
            });
        }
        for (node, point) in self.nodes_m.iter().enumerate() {
            if point.iter().any(|value| !value.is_finite()) {
                return Err(TetElasticError::NonFiniteNode { node });
            }
        }
        Ok(())
    }
}

fn free_dof_map(
    node_count: usize,
    fixed_dofs: &[usize],
    maximum_free_dofs: usize,
) -> Result<(Vec<usize>, Vec<Option<usize>>), TetElasticError> {
    let full_dofs = node_count
        .checked_mul(3)
        .ok_or(TetElasticError::SizeOverflow)?;
    let mut fixed = BTreeSet::new();
    for &dof in fixed_dofs {
        if dof >= full_dofs {
            return Err(TetElasticError::FixedDofOutOfRange { dof, full_dofs });
        }
        if !fixed.insert(dof) {
            return Err(TetElasticError::DuplicateFixedDof { dof });
        }
    }
    let free_dofs: Vec<usize> = (0..full_dofs).filter(|dof| !fixed.contains(dof)).collect();
    if free_dofs.len() > maximum_free_dofs {
        return Err(TetElasticError::BudgetExceeded {
            what: "free_dofs",
            requested: free_dofs.len(),
            maximum: maximum_free_dofs,
        });
    }
    let mut reduced_of = vec![None; full_dofs];
    for (reduced, &full) in free_dofs.iter().enumerate() {
        reduced_of[full] = Some(reduced);
    }
    Ok((free_dofs, reduced_of))
}

fn hash_usize(hasher: &mut DomainHasher, value: usize) -> Result<(), TetElasticError> {
    let value = u64::try_from(value).map_err(|_| TetElasticError::SizeOverflow)?;
    hasher.update(&value.to_le_bytes());
    Ok(())
}

fn euclidean_norm(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .fold(0.0, |norm, value| norm.hypot(value))
}

fn signed_tet_determinant(points: [[f64; 3]; 4]) -> f64 {
    let columns = [
        sub(points[1], points[0]),
        sub(points[2], points[0]),
        sub(points[3], points[0]),
    ];
    determinant3([
        [columns[0][0], columns[1][0], columns[2][0]],
        [columns[0][1], columns[1][1], columns[2][1]],
        [columns[0][2], columns[1][2], columns[2][2]],
    ])
}

fn tet_geometry_identity(
    nodes_m: &[[f64; 3]],
    tetrahedra: &[[usize; 4]],
) -> Result<ContentHash, TetElasticError> {
    let mut identity = DomainHasher::new("org.frankensim.fs-solid.tet-geometry.v1");
    hash_usize(&mut identity, nodes_m.len())?;
    for point in nodes_m {
        for coordinate in point {
            identity.update(&coordinate.to_bits().to_le_bytes());
        }
    }
    hash_usize(&mut identity, tetrahedra.len())?;
    for tet in tetrahedra {
        for &node in tet {
            hash_usize(&mut identity, node)?;
        }
    }
    Ok(identity.finalize())
}

fn validate_rigid_mode_constraints(
    nodes_m: &[[f64; 3]],
    tetrahedra: &[[usize; 4]],
    fixed_dofs: &[usize],
) -> Result<(), TetElasticError> {
    let mut parent = (0..nodes_m.len()).collect::<Vec<_>>();
    fn root(parent: &mut [usize], mut node: usize) -> usize {
        while parent[node] != node {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        node
    }
    for tet in tetrahedra {
        let first = root(&mut parent, tet[0]);
        for &node in &tet[1..] {
            let other = root(&mut parent, node);
            if first != other {
                parent[other] = first;
            }
        }
    }
    let roots = (0..nodes_m.len())
        .map(|node| root(&mut parent, node))
        .collect::<Vec<_>>();
    let mut components = BTreeSet::new();
    components.extend(roots.iter().copied());
    for component in components {
        let members = roots
            .iter()
            .enumerate()
            .filter_map(|(node, &candidate)| (candidate == component).then_some(node))
            .collect::<Vec<_>>();
        let member_count =
            u32::try_from(members.len()).map_err(|_| TetElasticError::SizeOverflow)?;
        let inverse_count = 1.0 / f64::from(member_count);
        let mut centroid = [0.0; 3];
        for &node in &members {
            for (coordinate, sum) in nodes_m[node].iter().zip(&mut centroid) {
                *sum += coordinate * inverse_count;
            }
        }
        let member_set = members.iter().copied().collect::<BTreeSet<_>>();
        let mut rows = Vec::new();
        for &dof in fixed_dofs {
            let node = dof / 3;
            if !member_set.contains(&node) {
                continue;
            }
            let component_axis = dof % 3;
            let r = [
                nodes_m[node][0] - centroid[0],
                nodes_m[node][1] - centroid[1],
                nodes_m[node][2] - centroid[2],
            ];
            let translations = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
            let rotations = [[0.0, -r[2], r[1]], [r[2], 0.0, -r[0]], [-r[1], r[0], 0.0]];
            rows.push([
                translations[0][component_axis],
                translations[1][component_axis],
                translations[2][component_axis],
                rotations[0][component_axis],
                rotations[1][component_axis],
                rotations[2][component_axis],
            ]);
        }
        let constrained_rank = normalized_rank_six(&mut rows);
        if constrained_rank < 6 {
            return Err(TetElasticError::UnconstrainedRigidModes {
                component: *members.first().unwrap_or(&component),
                constrained_rank,
            });
        }
    }
    Ok(())
}

fn normalized_rank_six(rows: &mut [[f64; 6]]) -> usize {
    for column in 0..6 {
        let norm = euclidean_norm(rows.iter().map(|row| row[column]));
        if norm > 0.0 {
            for row in rows.iter_mut() {
                row[column] /= norm;
            }
        }
    }
    let mut rank = 0;
    for column in 0..6 {
        let Some((offset, _)) = rows[rank..]
            .iter()
            .enumerate()
            .max_by(|(_, left), (_, right)| left[column].abs().total_cmp(&right[column].abs()))
        else {
            break;
        };
        let pivot = rank + offset;
        if rows[pivot][column].abs() <= 1.0e-12 {
            continue;
        }
        rows.swap(rank, pivot);
        let divisor = rows[rank][column];
        for value in &mut rows[rank][column..] {
            *value /= divisor;
        }
        let pivot_row = rows[rank];
        for row in rows.iter_mut().skip(rank + 1) {
            let factor = row[column];
            for index in column..6 {
                row[index] = (-factor).mul_add(pivot_row[index], row[index]);
            }
        }
        rank += 1;
        if rank == 6 {
            break;
        }
    }
    rank
}

fn validate_connectivity(
    element: usize,
    tet: &[usize; 4],
    node_count: usize,
) -> Result<(), TetElasticError> {
    for &node in tet {
        if node >= node_count {
            return Err(TetElasticError::NodeIndexOutOfRange { element, node });
        }
    }
    let unique: BTreeSet<usize> = tet.iter().copied().collect();
    if unique.len() != 4 {
        return Err(TetElasticError::RepeatedTetNode { element });
    }
    Ok(())
}

struct TetGeometry {
    gradients: [[f64; 3]; 4],
    volume: f64,
    scaled_jacobian: f64,
}

impl TetGeometry {
    fn try_new(element: usize, p: [[f64; 3]; 4], minimum: f64) -> Result<Self, TetElasticError> {
        let j = [sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0])];
        // `j` stores columns; spell the determinant and inverse in coordinate
        // rows to keep the gradient convention auditable.
        let a = [
            [j[0][0], j[1][0], j[2][0]],
            [j[0][1], j[1][1], j[2][1]],
            [j[0][2], j[1][2], j[2][2]],
        ];
        let det = determinant3(a);
        let h_max = maximum_edge_length(p);
        let scaled_jacobian = det.abs() / (h_max * h_max * h_max);
        if !scaled_jacobian.is_finite() || scaled_jacobian < minimum {
            return Err(TetElasticError::DegenerateTet {
                element,
                scaled_jacobian,
                minimum,
            });
        }
        let inv = inverse3(a, det);
        // grad_x N_i = J^{-T} grad_ref N_i. Rows of J^{-1} are the
        // physical gradients of reference coordinates xi_i.
        let gradients = [
            [
                -inv[0][0] - inv[1][0] - inv[2][0],
                -inv[0][1] - inv[1][1] - inv[2][1],
                -inv[0][2] - inv[1][2] - inv[2][2],
            ],
            inv[0],
            inv[1],
            inv[2],
        ];
        Ok(Self {
            gradients,
            volume: det.abs() / 6.0,
            scaled_jacobian,
        })
    }
}

fn element_matrices(
    geometry: &TetGeometry,
    material: &TetElasticMaterial,
) -> ([[f64; 12]; 12], [[f64; 12]; 12]) {
    let d = material.stiffness_mandel_pa;
    let b = strain_displacement_matrix(geometry);
    let mut db = [[0.0; 12]; 6];
    for i in 0..6 {
        for j in 0..12 {
            for k in 0..6 {
                db[i][j] = d[i][k].mul_add(b[k][j], db[i][j]);
            }
        }
    }
    let mut ke = [[0.0; 12]; 12];
    for i in 0..12 {
        for j in 0..12 {
            for k in 0..6 {
                ke[i][j] = b[k][i].mul_add(db[k][j], ke[i][j]);
            }
            ke[i][j] *= geometry.volume;
        }
    }

    let mut me = [[0.0; 12]; 12];
    let mass_scale = material.density_kg_m3 * geometry.volume / 20.0;
    for a in 0..4 {
        for b_node in 0..4 {
            let value = mass_scale * if a == b_node { 2.0 } else { 1.0 };
            for component in 0..3 {
                me[3 * a + component][3 * b_node + component] = value;
            }
        }
    }
    (ke, me)
}

fn strain_displacement_matrix(geometry: &TetGeometry) -> [[f64; 12]; 6] {
    let mut b = [[0.0; 12]; 6];
    let inverse_sqrt_two = core::f64::consts::FRAC_1_SQRT_2;
    for (node, g) in geometry.gradients.iter().enumerate() {
        let column = 3 * node;
        b[0][column] = g[0];
        b[1][column + 1] = g[1];
        b[2][column + 2] = g[2];
        b[3][column] = g[1] * inverse_sqrt_two;
        b[3][column + 1] = g[0] * inverse_sqrt_two;
        b[4][column + 1] = g[2] * inverse_sqrt_two;
        b[4][column + 2] = g[1] * inverse_sqrt_two;
        b[5][column] = g[2] * inverse_sqrt_two;
        b[5][column + 2] = g[0] * inverse_sqrt_two;
    }
    b
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn is_positive_definite(matrix: MandelStiffness6) -> bool {
    let mut lower = [[0.0_f64; 6]; 6];
    for row in 0..6 {
        for column in 0..=row {
            let mut remainder = matrix[row][column];
            for inner in 0..column {
                remainder = lower[row][inner].mul_add(-lower[column][inner], remainder);
            }
            if row == column {
                if !(remainder.is_finite() && remainder > 0.0) {
                    return false;
                }
                lower[row][column] = remainder.sqrt();
            } else {
                lower[row][column] = remainder / lower[column][column];
                if !lower[row][column].is_finite() {
                    return false;
                }
            }
        }
    }
    true
}

fn validate_rotation(rotation: [[f64; 3]; 3]) -> Result<(), TetElasticError> {
    if rotation.iter().flatten().any(|value| !value.is_finite()) {
        return Err(TetElasticError::InvalidMaterial {
            what: "principal_to_world must be finite",
        });
    }
    const ORTHONORMAL_TOLERANCE: f64 = 1.0e-12;
    for column_a in 0..3 {
        for column_b in 0..3 {
            let dot = (0..3)
                .map(|row| rotation[row][column_a] * rotation[row][column_b])
                .sum::<f64>();
            let target = if column_a == column_b { 1.0 } else { 0.0 };
            if (dot - target).abs() > ORTHONORMAL_TOLERANCE {
                return Err(TetElasticError::InvalidMaterial {
                    what: "principal_to_world must be orthonormal",
                });
            }
        }
    }
    if (determinant3(rotation) - 1.0).abs() > ORTHONORMAL_TOLERANCE {
        return Err(TetElasticError::InvalidMaterial {
            what: "principal_to_world must be a proper rotation",
        });
    }
    Ok(())
}

fn normalize_elastic_tensor(
    source: [[f64; 6]; 6],
    basis: ElasticTensorBasis,
) -> Result<MandelStiffness6, TetElasticError> {
    if source.iter().flatten().any(|value| !value.is_finite()) {
        return Err(TetElasticError::InvalidMaterial {
            what: "source elastic tensor must be finite",
        });
    }
    let order = match basis.order {
        ElasticTensorOrder::XxYyZzXyYzZx => [0, 1, 2, 3, 4, 5],
        ElasticTensorOrder::XxYyZzYzZxXy => [0, 1, 2, 5, 3, 4],
    };
    let mut conjugate = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let value = source[order[i]][order[j]];
            // Tensor stress = C_tensor * epsilon_tensor. Engineering strain
            // is W epsilon_tensor, so C_engineering = C_tensor W^-1.
            let normalized = if basis.notation == ElasticTensorNotation::Tensor && j >= 3 {
                value * 0.5
            } else {
                value
            };
            if value != 0.0 && normalized == 0.0 {
                return Err(TetElasticError::InvalidMaterial {
                    what: "elastic tensor shear conversion underflows",
                });
            }
            conjugate[i][j] = normalized;
        }
    }
    let mut mandel = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in i..6 {
            if conjugate[i][j].to_bits() != conjugate[j][i].to_bits() {
                return Err(TetElasticError::InvalidMaterial {
                    what: "elastic tensor must have exact major symmetry in its declared notation",
                });
            }
            // C_mandel = S C_engineering S, S=diag(1,1,1,sqrt(2),...).
            // One product per pair avoids manufacturing a final-bit asymmetry.
            let scale = if basis.notation == ElasticTensorNotation::Mandel {
                1.0
            } else {
                match (i >= 3, j >= 3) {
                    (false, false) => 1.0,
                    (true, true) => 2.0,
                    _ => core::f64::consts::SQRT_2,
                }
            };
            let value = conjugate[i][j] * scale;
            mandel[i][j] = value;
            mandel[j][i] = value;
        }
    }
    Ok(mandel)
}

fn rotate_mandel_stiffness(
    principal: MandelStiffness6,
    rotation: [[f64; 3]; 3],
) -> MandelStiffness6 {
    let transform = mandel_rotation(rotation);
    let mut world = [[0.0; 6]; 6];
    for a in 0..6 {
        for b in 0..6 {
            for p in 0..6 {
                for q in 0..6 {
                    world[a][b] =
                        transform[a][p].mul_add(principal[p][q] * transform[b][q], world[a][b]);
                }
            }
        }
    }
    // Rotation is mathematically symmetric; retain the original orthotropic
    // operation tree and its deterministic triangle convention. Input symmetry
    // is checked before this step; this is not a projection of noisy source data.
    for row in 0..6 {
        for column in 0..row {
            world[row][column] = world[column][row];
        }
    }
    world
}

fn mandel_rotation(rotation: [[f64; 3]; 3]) -> [[f64; 6]; 6] {
    let basis = mandel_basis();
    let mut transform = [[0.0; 6]; 6];
    for principal in 0..6 {
        let mut rotated = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for p in 0..3 {
                    for q in 0..3 {
                        rotated[i][j] = rotation[i][p]
                            .mul_add(basis[principal][p][q] * rotation[j][q], rotated[i][j]);
                    }
                }
            }
        }
        for world in 0..6 {
            for i in 0..3 {
                for j in 0..3 {
                    transform[world][principal] =
                        basis[world][i][j].mul_add(rotated[i][j], transform[world][principal]);
                }
            }
        }
    }
    transform
}

fn mandel_basis() -> [[[f64; 3]; 3]; 6] {
    let mut basis = [[[0.0; 3]; 3]; 6];
    basis[0][0][0] = 1.0;
    basis[1][1][1] = 1.0;
    basis[2][2][2] = 1.0;
    let shear = core::f64::consts::FRAC_1_SQRT_2;
    basis[3][0][1] = shear;
    basis[3][1][0] = shear;
    basis[4][1][2] = shear;
    basis[4][2][1] = shear;
    basis[5][2][0] = shear;
    basis[5][0][2] = shear;
    basis
}

fn determinant3(a: [[f64; 3]; 3]) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

fn inverse3(a: [[f64; 3]; 3], det: f64) -> [[f64; 3]; 3] {
    let inv_det = det.recip();
    [
        [
            (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * inv_det,
            (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * inv_det,
            (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * inv_det,
        ],
        [
            (a[1][2] * a[2][0] - a[1][0] * a[2][2]) * inv_det,
            (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * inv_det,
            (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * inv_det,
        ],
        [
            (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * inv_det,
            (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * inv_det,
            (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * inv_det,
        ],
    ]
}

fn maximum_edge_length(p: [[f64; 3]; 4]) -> f64 {
    let mut max_squared: f64 = 0.0;
    for i in 0..4 {
        for j in i + 1..4 {
            let d = sub(p[i], p[j]);
            max_squared = max_squared.max(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        }
    }
    max_squared.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_mesh::{RoundedCylinderMeshSpec, rounded_cylinder_tet_mesh};
    use fs_modal::eigh_gen_dense;

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 3,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn material(e: f64, rho: f64) -> TetElasticMaterial {
        TetElasticMaterial::try_new(rho, e, 0.25, [7; 32]).unwrap()
    }

    fn reference_problem<'a>(
        nodes: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        mat: &'a TetElasticMaterial,
        fixed: &'a [usize],
    ) -> TetLinearElasticProblem<'a> {
        TetLinearElasticProblem {
            nodes_m: nodes,
            tetrahedra,
            materials: TetMaterialField::Uniform(mat),
            fixed_dofs: fixed,
            budget: TetAssemblyBudget::standard(),
        }
    }

    #[test]
    fn g0_single_tet_has_rigid_translation_nullspace_and_exact_mass() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let mat = material(120.0, 6.0);
        let assembly =
            with_cx(|cx| reference_problem(&nodes, &tets, &mat, &[]).assemble(cx)).unwrap();
        assert_eq!(assembly.total_mass_kg, 1.0);
        let k = assembly.stiffness.to_dense();
        let m = assembly.mass.to_dense();
        for component in 0..3 {
            let u: Vec<f64> = (0..12).map(|dof| f64::from(dof % 3 == component)).collect();
            for row in 0..12 {
                let ku: f64 = (0..12).map(|col| k[row * 12 + col] * u[col]).sum();
                assert!(ku.abs() < 1.0e-12, "row {row}: {ku}");
            }
            let directional_mass: f64 = (0..12)
                .map(|row| {
                    (0..12)
                        .map(|col| u[row] * m[row * 12 + col] * u[col])
                        .sum::<f64>()
                })
                .sum();
            assert!((directional_mass - 1.0).abs() < 1.0e-14);
        }
    }

    #[test]
    fn g0_affine_strain_energy_and_rigid_rotation_are_exact() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let mat = material(120.0, 6.0);
        let assembly =
            with_cx(|cx| reference_problem(&nodes, &tets, &mat, &[]).assemble(cx)).unwrap();
        let k = assembly.stiffness.to_dense();

        // Infinitesimal rotation about z: u=(-y,x,0) has zero symmetric strain.
        let rotation: Vec<f64> = nodes.iter().flat_map(|p| [-p[1], p[0], 0.0]).collect();
        let rotation_energy = quadratic(&k, &rotation);
        assert!(rotation_energy.abs() < 1.0e-12, "{rotation_energy}");

        // u=(alpha*x,0,0): eps_xx=alpha. Analytic energy is
        // 1/2 * (lambda+2mu) * alpha^2 * V.
        let alpha = 0.125;
        let extension: Vec<f64> = nodes
            .iter()
            .flat_map(|p| [alpha * p[0], 0.0, 0.0])
            .collect();
        let expected = 0.5 * mat.stiffness_mandel_pa()[0][0] * alpha * alpha / 6.0;
        let actual = 0.5 * quadratic(&k, &extension);
        assert!(
            (actual - expected).abs() < 1.0e-13,
            "{actual} != {expected}"
        );
    }

    fn portable_strain_pack(
        basis: StrainTensorBasis,
        strain: [f64; 6],
        uncertainty: fs_matdb::UncertaintyModel,
        covariance: Option<([[f64; 6]; 6], usize)>,
    ) -> (
        fs_matdb::NormalizedMaterialCardPack,
        [fs_matdb::PropertyKey; 6],
    ) {
        use fs_matdb::{
            ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack,
            NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue,
            Provenance, StrainTensorComponent,
        };
        let mut claims = ClaimSet::new();
        let provenance = Provenance {
            source: "synthetic integrated strain".into(),
            license: "synthetic test data".into(),
            artifact: Some(ContentHash([5; 32])),
        };
        let observation = claims
            .register_observation(ObservationDataset {
                specimen: "synthetic".into(),
                method: "analytic free strain fixture".into(),
                artifact: ContentHash([5; 32]),
                caveats: "no empirical or path-integration claim".into(),
                provenance: provenance.clone(),
            })
            .unwrap();
        let mut claim_indices = Vec::new();
        let keys = core::array::from_fn(|index| {
            let key = PropertyKey::new("strain", fs_qty::Dims::NONE)
                .with_strain_component(
                    StrainTensorComponent::new(basis, ContentHash([4; 32]), index as u8).unwrap(),
                )
                .unwrap();
            let id = claims
                .insert_claim(PropertyClaim {
                    key: key.clone(),
                    value: PropertyValue::Scalar {
                        value: strain[index],
                        dims: fs_qty::Dims::NONE,
                    },
                    validity: fs_evidence::ValidityDomain::unconstrained()
                        .with("T", 400.0, 400.0)
                        .with("Tref", 293.15, 293.15),
                    uncertainty: uncertainty.clone(),
                    interpolation: InterpolationPolicy::TabulatedOnly,
                    observations: vec![observation],
                    provenance: provenance.clone(),
                })
                .unwrap();
            claim_indices.push((id, index));
            key
        });
        let mut blocks = Vec::new();
        if let Some((covariance, count)) = covariance {
            claim_indices.truncate(count);
            claim_indices.sort_by_key(|(id, _)| *id);
            let mut packed = Vec::new();
            for (row, (_, source_i)) in claim_indices.iter().enumerate() {
                for (_, source_j) in &claim_indices[..=row] {
                    packed.push(covariance[*source_i][*source_j]);
                }
            }
            blocks.push(fs_matdb::JointStatistics::new(
                observation,
                "synthetic correlated strains",
                claim_indices
                    .iter()
                    .map(|(id, _)| fs_matdb::StatisticMember::scalar(*id))
                    .collect(),
                packed,
                None,
            ));
        }
        let pack = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: "synthetic".into(),
                phase: "solid".into(),
                process: "integrated strain fixture".into(),
                revision: 0,
            },
            NormalizedPack::new(
                "synthetic-strain",
                "fixture-v6",
                ContentHash([5; 32]),
                "synthetic redistribution permitted",
                claims,
                blocks,
                Vec::new(),
            )
            .unwrap(),
        )
        .unwrap();
        let portable =
            NormalizedMaterialCardPack::from_bytes_verified(pack.content_hash(), &pack.to_bytes())
                .unwrap();
        (portable, keys)
    }

    fn portable_strain_state(basis: StrainTensorBasis, strain: [f64; 6]) -> StrainTensorStatePoint {
        use fs_matdb::{QueryPoint, UncertaintyModel};
        use fs_material::state_point::{
            MaterialPropertySelection, resolve_strain_tensor_state_point,
        };
        let (portable, keys) =
            portable_strain_pack(basis, strain, UncertaintyModel::Unstated, None);
        let point = QueryPoint::new()
            .with("T", 400.0)
            .unwrap()
            .with("Tref", 293.15)
            .unwrap();
        let state = resolve_strain_tensor_state_point(
            portable.card(),
            &point,
            &keys,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert_eq!(state.strain(), &strain);
        assert_eq!(state.basis(), basis);
        assert_eq!(state.resolved().properties().len(), 6);
        state
    }

    #[test]
    fn g3_joint_strain_covariance_matches_correlated_ensemble_and_thermal_forces() {
        use fs_matdb::{QueryPoint, SelectionPolicy, UncertaintyModel};
        use fs_material::state_point::resolve_joint_strain_tensor_state_point;
        let q = fs_material::tensor::rotation([2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0], 0.73);
        let elastic = TetElasticMaterial::try_new_oriented_tensor(
            6.0,
            coupled_engineering_stiffness(),
            ElasticTensorBasis {
                notation: ElasticTensorNotation::Engineering,
                order: ElasticTensorOrder::XxYyZzXyYzZx,
                frame: ContentHash([1; 32]),
            },
            ContentHash([2; 32]),
            q,
            ContentHash([3; 32]),
        )
        .unwrap();
        let mean = [0.01, -0.004, 0.002, 0.003, -0.001, 0.002];
        let factors = [
            [1.0, 0.5, -0.25, 0.125, 0.0, 0.25],
            [-0.25, 0.25, 1.0, -0.25, 0.5, 0.0],
            [0.0, 0.0, 0.5, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.5, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.5],
        ];
        let scale = 1.0 / 4096.0;
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let problem = reference_problem(&nodes, &tets, elastic.material(), &[0, 1, 2, 4, 5, 8]);
        let point = QueryPoint::new()
            .with("T", 400.0)
            .unwrap()
            .with("Tref", 293.15)
            .unwrap();
        // Twelve equally likely atoms: +/-sqrt(6) times each physical tensor
        // factor. This has full-rank covariance L L^T; no sampling or covariance
        // transform under test is used to construct the reference ensemble.
        let mut strain_atoms = Vec::new();
        let mut force_atoms = Vec::new();
        for factor in factors {
            for sign in [-1.0, 1.0] {
                let physical =
                    core::array::from_fn(|i| mean[i] + sign * 6.0_f64.sqrt() * scale * factor[i]);
                let rotated = fs_material::tensor::rotate(&physical, &q);
                let mandel = core::array::from_fn(|i| {
                    rotated[i]
                        * if i < 3 {
                            1.0
                        } else {
                            core::f64::consts::SQRT_2
                        }
                });
                let state = TetThermalStrainState::try_new(
                    400.0,
                    293.15,
                    mandel,
                    elastic.material().material_state_identity,
                    ContentHash([7; 32]),
                    0.05,
                )
                .unwrap();
                strain_atoms.push(mandel);
                force_atoms.push(
                    with_cx(|cx| {
                        problem.assemble_thermal_load(TetThermalStrainField::Uniform(&state), cx)
                    })
                    .unwrap()
                    .full_equivalent_force_n,
                );
            }
        }
        assert_eq!(strain_atoms.len(), 12);
        assert_eq!(force_atoms.len(), 12);
        let atom_covariance = |samples: &[Vec<f64>], i: usize, j: usize| {
            let count = samples.len() as f64;
            let mi = samples.iter().map(|v| v[i]).sum::<f64>() / count;
            let mj = samples.iter().map(|v| v[j]).sum::<f64>() / count;
            samples
                .iter()
                .map(|v| (v[i] - mi) * (v[j] - mj))
                .sum::<f64>()
                / count
        };
        let strain_atoms: Vec<_> = strain_atoms.iter().map(|v| v.to_vec()).collect();
        let gradients = [[-1.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let mut force_map = [[0.0; 6]; 12];
        for (node, [x, y, z]) in gradients.into_iter().enumerate() {
            let h = core::f64::consts::FRAC_1_SQRT_2;
            force_map[3 * node] = [x / 6.0, 0.0, 0.0, y * h / 6.0, 0.0, z * h / 6.0];
            force_map[3 * node + 1] = [0.0, y / 6.0, 0.0, x * h / 6.0, z * h / 6.0, 0.0];
            force_map[3 * node + 2] = [0.0, 0.0, z / 6.0, 0.0, y * h / 6.0, x * h / 6.0];
        }
        let mut identities = BTreeSet::new();
        for notation in [
            StrainTensorNotation::Tensor,
            StrainTensorNotation::Engineering,
            StrainTensorNotation::Mandel,
        ] {
            for order in [
                ElasticTensorOrder::XxYyZzXyYzZx,
                ElasticTensorOrder::XxYyZzYzZxXy,
            ] {
                let indices = if order == ElasticTensorOrder::XxYyZzXyYzZx {
                    [0, 1, 2, 3, 4, 5]
                } else {
                    [0, 1, 2, 4, 5, 3]
                };
                let shear = match notation {
                    StrainTensorNotation::Tensor => 1.0,
                    StrainTensorNotation::Engineering => 2.0,
                    StrainTensorNotation::Mandel => core::f64::consts::SQRT_2,
                };
                let multipliers: [f64; 6] =
                    core::array::from_fn(|i| if i < 3 { 1.0 } else { shear });
                let input = core::array::from_fn(|i| mean[indices[i]] * multipliers[i]);
                let covariance = core::array::from_fn(|i| {
                    core::array::from_fn(|j| {
                        factors
                            .iter()
                            .map(|l| l[indices[i]] * l[indices[j]])
                            .sum::<f64>()
                            * scale
                            * scale
                            * multipliers[i]
                            * multipliers[j]
                    })
                });
                let basis = StrainTensorBasis {
                    notation,
                    order,
                    frame: ContentHash([1; 32]),
                };
                let (pack, keys) = portable_strain_pack(
                    basis,
                    input,
                    UncertaintyModel::HalfWidth {
                        half_width: 0.001,
                        confidence: 0.95,
                    },
                    Some((covariance, 6)),
                );
                for policy in [
                    SelectionPolicy::SingleClaimOnly,
                    SelectionPolicy::PreferObservationBacked,
                ] {
                    let state =
                        resolve_joint_strain_tensor_state_point(&pack, &point, &keys, policy)
                            .unwrap();
                    pack.claims_pack()
                        .verify_joint_receipt(state.joint_receipt())
                        .unwrap();
                    let result = TetThermalStrainState::try_from_joint_strain(
                        &elastic, &state, 400.0, 293.15, q, 0.05,
                    )
                    .unwrap();
                    assert_eq!(
                        result,
                        TetThermalStrainState::try_from_joint_strain(
                            &elastic, &state, 400.0, 293.15, q, 0.05
                        )
                        .unwrap()
                    );
                    assert_eq!(result.source_joint_identity(), state.identity());
                    assert!(identities.insert(result.identity()));
                    let ThermalStrainUncertainty::Covariance {
                        strain_mandel,
                        free_stress_mandel_pa2,
                    } = result.uncertainty()
                    else {
                        panic!("admitted full covariance must survive")
                    };
                    for i in 0..6 {
                        for j in 0..6 {
                            let expected = atom_covariance(&strain_atoms, i, j);
                            assert!(
                                (strain_mandel[i][j] - expected).abs() < 1.0e-18,
                                "{notation:?}/{order:?} ({i},{j})"
                            );
                            assert_eq!(
                                strain_mandel[i][j].to_bits(),
                                strain_mandel[j][i].to_bits()
                            );
                        }
                    }
                    for i in 0..12 {
                        for j in 0..12 {
                            let mut actual = 0.0;
                            for a in 0..6 {
                                for b in 0..6 {
                                    actual += force_map[i][a]
                                        * free_stress_mandel_pa2[a][b]
                                        * force_map[j][b];
                                }
                            }
                            let expected = atom_covariance(&force_atoms, i, j);
                            assert!(
                                (actual - expected).abs() < 1.0e-14,
                                "nodal covariance ({i},{j}): {actual} != {expected}"
                            );
                        }
                    }
                    let mut tampered = state.joint_receipt().clone();
                    if let JointCorrelation::Covariance { covariance, .. } =
                        &mut tampered.correlation
                    {
                        covariance[0] *= 2.0;
                    }
                    assert!(pack.claims_pack().verify_joint_receipt(&tampered).is_err());
                }
            }
        }
        assert_eq!(identities.len(), 12);
    }

    #[test]
    fn g0_joint_strain_preserves_unknown_correlation_and_refuses_numeric_loss() {
        use fs_matdb::{QueryPoint, SelectionPolicy, UncertaintyModel};
        use fs_material::state_point::resolve_joint_strain_tensor_state_point;
        let basis = StrainTensorBasis {
            notation: StrainTensorNotation::Mandel,
            order: ElasticTensorOrder::XxYyZzXyYzZx,
            frame: ContentHash([1; 32]),
        };
        let q = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let elastic = TetElasticMaterial::try_new_oriented_tensor(
            6.0,
            coupled_engineering_stiffness(),
            ElasticTensorBasis {
                notation: ElasticTensorNotation::Engineering,
                order: basis.order,
                frame: basis.frame,
            },
            ContentHash([2; 32]),
            q,
            ContentHash([3; 32]),
        )
        .unwrap();
        let point = QueryPoint::new()
            .with("T", 400.0)
            .unwrap()
            .with("Tref", 293.15)
            .unwrap();
        let covariance =
            core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0e-6 } else { 0.0 }));
        let stated = UncertaintyModel::HalfWidth {
            half_width: 0.001,
            confidence: 0.95,
        };
        for (uncertainty, block, reason) in [
            (stated.clone(), None, CorrelationUnknownReason::NoBlock),
            (
                stated,
                Some((covariance, 5)),
                CorrelationUnknownReason::PartialMembership,
            ),
            (
                UncertaintyModel::Unstated,
                Some((covariance, 6)),
                CorrelationUnknownReason::UnstatedMarginal,
            ),
        ] {
            let (pack, keys) = portable_strain_pack(basis, [0.01; 6], uncertainty, block);
            let state = resolve_joint_strain_tensor_state_point(
                &pack,
                &point,
                &keys,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
            let result = TetThermalStrainState::try_from_joint_strain(
                &elastic, &state, 400.0, 293.15, q, 0.05,
            )
            .unwrap();
            assert_eq!(
                result.uncertainty(),
                &ThermalStrainUncertainty::Unknown { reason }
            );
            assert_eq!(
                result.nominal(),
                &TetThermalStrainState::try_from_resolved_strain(
                    &elastic,
                    state.nominal(),
                    400.0,
                    293.15,
                    q,
                    0.05
                )
                .unwrap()
            );
        }
        let identity =
            core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
        assert_eq!(
            transform_covariance6(&identity, &[[0.0; 6]; 6]).unwrap(),
            [[0.0; 6]; 6]
        );
        for scale in [f64::MAX, 1.0e-200] {
            let map = identity.map(|row| row.map(|v| v * scale));
            assert!(transform_covariance6(&map, &identity).is_err());
        }
    }

    #[test]
    fn g1_g3_oriented_strain_conventions_drive_anisotropic_thermal_expansion() {
        let q = fs_material::tensor::rotation([2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0], 0.73);
        let elastic = TetElasticMaterial::try_new_oriented_tensor(
            6.0,
            coupled_engineering_stiffness(),
            ElasticTensorBasis {
                notation: ElasticTensorNotation::Engineering,
                order: ElasticTensorOrder::XxYyZzXyYzZx,
                frame: ContentHash([1; 32]),
            },
            ContentHash([2; 32]),
            q,
            ContentHash([3; 32]),
        )
        .unwrap();
        let source = [
            [0.010, 0.003, 0.002],
            [0.003, -0.004, -0.001],
            [0.002, -0.001, 0.002],
        ];
        // Independent second-order oracle, with no six-vector transform.
        let mut target = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for p in 0..3 {
                    for r in 0..3 {
                        target[i][j] += q[i][p] * source[p][r] * q[j][r];
                    }
                }
            }
        }
        let expected_mandel = [
            target[0][0],
            target[1][1],
            target[2][2],
            core::f64::consts::SQRT_2 * target[0][1],
            core::f64::consts::SQRT_2 * target[1][2],
            core::f64::consts::SQRT_2 * target[2][0],
        ];
        let source_norm2: f64 = source.iter().flatten().map(|v| v * v).sum();
        let source_trace = source[0][0] + source[1][1] + source[2][2];
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let problem = reference_problem(&nodes, &tets, elastic.material(), &[0, 1, 2, 4, 5, 8]);
        // A skew displacement gradient removes exactly the constrained rigid
        // rotation while preserving sym(grad u) = the independent target strain.
        let expected_displacement = [
            [0.0; 3],
            [target[0][0], 0.0, 0.0],
            [2.0 * target[0][1], target[1][1], 0.0],
            [2.0 * target[0][2], 2.0 * target[1][2], target[2][2]],
        ];
        let c = four_index_rotated_stiffness(q);
        let mut stress = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        stress[i][j] += c[i][j][k][l] * target[k][l];
                    }
                }
            }
        }
        let gradients = [[-1.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let mut identities = BTreeSet::new();
        for notation in [
            StrainTensorNotation::Tensor,
            StrainTensorNotation::Engineering,
            StrainTensorNotation::Mandel,
        ] {
            let shear_scale = match notation {
                StrainTensorNotation::Tensor => 1.0,
                StrainTensorNotation::Engineering => 2.0,
                StrainTensorNotation::Mandel => core::f64::consts::SQRT_2,
            };
            let canonical = [
                source[0][0],
                source[1][1],
                source[2][2],
                source[0][1] * shear_scale,
                source[1][2] * shear_scale,
                source[2][0] * shear_scale,
            ];
            for order in [
                ElasticTensorOrder::XxYyZzXyYzZx,
                ElasticTensorOrder::XxYyZzYzZxXy,
            ] {
                let input = match order {
                    ElasticTensorOrder::XxYyZzXyYzZx => canonical,
                    ElasticTensorOrder::XxYyZzYzZxXy => [
                        canonical[0],
                        canonical[1],
                        canonical[2],
                        canonical[4],
                        canonical[5],
                        canonical[3],
                    ],
                };
                let basis = StrainTensorBasis {
                    notation,
                    order,
                    frame: ContentHash([1; 32]),
                };
                let receipt = TetThermalStrainState::try_from_oriented_strain(
                    &elastic,
                    400.0,
                    293.15,
                    input,
                    basis,
                    q,
                    ContentHash([0x41; 32]),
                    0.05,
                )
                .unwrap();
                assert_eq!(receipt.source_free_strain(), &input);
                assert_eq!(receipt.source_basis(), basis);
                assert_eq!(receipt.target_frame(), elastic.target_frame());
                assert_eq!(receipt.source_to_target(), &q);
                assert_eq!(
                    receipt.state().thermal_law_identity(),
                    ContentHash([0x41; 32])
                );
                let actual = receipt.state().free_strain_mandel();
                for (actual, expected) in actual.iter().zip(expected_mandel) {
                    assert!(
                        (actual - expected).abs() < 1.0e-16,
                        "{notation:?}/{order:?}: {actual} != {expected}"
                    );
                }
                assert!((actual[..3].iter().sum::<f64>() - source_trace).abs() < 1.0e-16);
                assert!((actual.iter().map(|v| v * v).sum::<f64>() - source_norm2).abs() < 1.0e-18);
                assert_eq!(
                    receipt,
                    TetThermalStrainState::try_from_oriented_strain(
                        &elastic,
                        400.0,
                        293.15,
                        input,
                        basis,
                        q,
                        ContentHash([0x41; 32]),
                        0.05,
                    )
                    .unwrap()
                );
                assert!(
                    identities.insert(receipt.state().identity()),
                    "source convention stays identity-bearing"
                );
                let resolved = portable_strain_state(basis, input);
                let sourced = TetThermalStrainState::try_from_resolved_strain(
                    &elastic, &resolved, 400.0, 293.15, q, 0.05,
                )
                .unwrap();
                assert_eq!(sourced.source_free_strain(), receipt.source_free_strain());
                assert_eq!(
                    sourced.state().free_strain_mandel(),
                    receipt.state().free_strain_mandel()
                );
                assert_eq!(
                    sourced.state().thermal_law_identity(),
                    resolved.resolved().identity()
                );
                assert_ne!(sourced.state().identity(), receipt.state().identity());
                let field = TetThermalStrainField::Uniform(sourced.state());
                let load = with_cx(|cx| problem.assemble_thermal_load(field, cx)).unwrap();
                for (node, gradient) in gradients.iter().enumerate() {
                    for component in 0..3 {
                        let expected = (0..3)
                            .map(|j| stress[component][j] * gradient[j])
                            .sum::<f64>()
                            / 6.0;
                        assert!(
                            (load.full_equivalent_force_n[3 * node + component] - expected).abs()
                                < 2.0e-14
                        );
                    }
                }
                let solution = with_cx(|cx| {
                    problem.solve_thermal_displacement(field, TetStaticSolveConfig::default(), cx)
                })
                .unwrap();
                assert!(solution.true_relative_residual() < 1.0e-12);
                let updated =
                    with_cx(|cx| problem.update_geometry_from_displacement(&solution, 0.1, cx))
                        .unwrap();
                for node in 0..4 {
                    for component in 0..3 {
                        assert!(
                            (solution.displacement_m()[node][component]
                                - expected_displacement[node][component])
                                .abs()
                                < 1.0e-12
                        );
                        assert!(
                            (updated.nodes_m()[node][component]
                                - nodes[node][component]
                                - expected_displacement[node][component])
                                .abs()
                                < 1.0e-12
                        );
                    }
                }
                assert_eq!(updated.tetrahedra(), &tets);
            }
        }
        assert_eq!(identities.len(), 6);
    }

    #[test]
    fn g0_oriented_thermal_strain_refuses_bad_frames_and_lost_values() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let elastic = TetElasticMaterial::try_new_oriented_tensor(
            6.0,
            coupled_engineering_stiffness(),
            ElasticTensorBasis {
                notation: ElasticTensorNotation::Engineering,
                order: ElasticTensorOrder::XxYyZzXyYzZx,
                frame: ContentHash([1; 32]),
            },
            ContentHash([2; 32]),
            identity,
            ContentHash([3; 32]),
        )
        .unwrap();
        let basis = StrainTensorBasis {
            notation: StrainTensorNotation::Mandel,
            order: ElasticTensorOrder::XxYyZzXyYzZx,
            frame: ContentHash([1; 32]),
        };
        let make = |input, basis, q, limit| {
            TetThermalStrainState::try_from_oriented_strain(
                &elastic,
                400.0,
                293.15,
                input,
                basis,
                q,
                ContentHash([0x41; 32]),
                limit,
            )
        };
        let input = [0.01, -0.003, 0.002, -0.001, 0.005, 0.007];
        assert_eq!(
            make(input, basis, identity, 0.01)
                .unwrap()
                .state()
                .free_strain_mandel()
                .map(f64::to_bits),
            input.map(f64::to_bits)
        );
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut values = input;
            values[4] = bad;
            assert!(matches!(
                make(values, basis, identity, 0.05),
                Err(TetElasticError::InvalidThermalStrain { .. })
            ));
        }
        assert!(matches!(
            make(
                input,
                StrainTensorBasis {
                    frame: ContentHash([0; 32]),
                    ..basis
                },
                identity,
                0.05
            ),
            Err(TetElasticError::InvalidThermalStrain { .. })
        ));
        for q in [
            [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[1.0, 0.1, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[f64::NAN, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        ] {
            assert!(matches!(
                make(input, basis, q, 0.05),
                Err(TetElasticError::InvalidMaterial { .. })
            ));
        }
        let q = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        assert!(matches!(
            make(
                input,
                StrainTensorBasis {
                    frame: elastic.target_frame(),
                    ..basis
                },
                q,
                0.05
            ),
            Err(TetElasticError::InvalidThermalStrain { .. })
        ));
        for limit in [0.005, 0.0, f64::INFINITY] {
            assert!(matches!(
                make(input, basis, identity, limit),
                Err(TetElasticError::InvalidThermalStrain { .. })
            ));
        }
        let mut overflowing = input;
        overflowing[3] = f64::MAX;
        assert!(matches!(
            make(
                overflowing,
                StrainTensorBasis {
                    notation: StrainTensorNotation::Tensor,
                    ..basis
                },
                identity,
                f64::MAX
            ),
            Err(TetElasticError::InvalidThermalStrain { .. })
        ));
        // Rank-one strain along an equally inclined axis has all six Mandel
        // coefficients below half a minimum subnormal: refuse an invented zero.
        let a = 1.0 / 3.0_f64.sqrt();
        let b = 1.0 / 2.0_f64.sqrt();
        let c = 1.0 / 6.0_f64.sqrt();
        let inclined = [[a, -b, -c], [a, b, -c], [a, 0.0, 2.0 * c]];
        let tiny = [f64::from_bits(1), 0.0, 0.0, 0.0, 0.0, 0.0];
        assert!(matches!(
            make(tiny, basis, inclined, 0.05),
            Err(TetElasticError::InvalidThermalStrain {
                what: "strain rotation underflows the complete nonzero tensor"
            })
        ));
        assert_eq!(
            make(tiny, basis, identity, 0.05)
                .unwrap()
                .state()
                .free_strain_mandel()
                .map(f64::to_bits),
            tiny.map(f64::to_bits)
        );
        // Even zero expansion cannot discard the declared frame transform.
        let zero = make([0.0; 6], basis, identity, 0.05).unwrap();
        let rotated_zero = make([0.0; 6], basis, q, 0.05).unwrap();
        assert_eq!(
            zero.state().free_strain_mandel(),
            rotated_zero.state().free_strain_mandel()
        );
        assert_ne!(zero.state().identity(), rotated_zero.state().identity());
        let relabelled = make(
            [0.0; 6],
            StrainTensorBasis {
                frame: ContentHash([9; 32]),
                ..basis
            },
            identity,
            0.05,
        )
        .unwrap();
        assert_ne!(zero.state().identity(), relabelled.state().identity());
    }

    #[test]
    fn g0_uniform_thermal_eigenstrain_reproduces_the_affine_free_expansion_load() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let mat = material(120.0, 6.0);
        let sqrt_two = core::f64::consts::SQRT_2;
        let tensor = [
            0.010,
            -0.004,
            0.002,
            sqrt_two * 0.003,
            sqrt_two * -0.001,
            sqrt_two * 0.002,
        ];
        let thermal = TetThermalStrainState::try_new(
            400.0,
            293.15,
            tensor,
            mat.material_state_identity,
            ContentHash([0x31; 32]),
            0.05,
        )
        .unwrap();
        let problem = reference_problem(&nodes, &tets, &mat, &[]);
        let assembly = with_cx(|cx| problem.assemble(cx)).unwrap();
        let load = with_cx(|cx| {
            problem.assemble_thermal_load(TetThermalStrainField::Uniform(&thermal), cx)
        })
        .unwrap();

        let exy = tensor[3] / sqrt_two;
        let eyz = tensor[4] / sqrt_two;
        let ezx = tensor[5] / sqrt_two;
        let displacement: Vec<f64> = nodes
            .iter()
            .flat_map(|point| {
                [
                    tensor[0] * point[0] + exy * point[1] + ezx * point[2],
                    exy * point[0] + tensor[1] * point[1] + eyz * point[2],
                    ezx * point[0] + eyz * point[1] + tensor[2] * point[2],
                ]
            })
            .collect();
        let stiffness = assembly.stiffness.to_dense();
        for row in 0..12 {
            let ku = (0..12)
                .map(|column| stiffness[row * 12 + column] * displacement[column])
                .sum::<f64>();
            assert!(
                (ku - load.full_equivalent_force_n[row]).abs() < 1.0e-12,
                "row {row}: {ku} != {}",
                load.full_equivalent_force_n[row]
            );
        }
        for component in 0..3 {
            let resultant: f64 = load
                .full_equivalent_force_n
                .iter()
                .skip(component)
                .step_by(3)
                .sum();
            assert!(
                resultant.abs() < 1.0e-12,
                "component {component}: {resultant}"
            );
        }
        assert_eq!(load.reduced_force_n, load.full_equivalent_force_n);
        assert_eq!(load.free_dofs, (0..12).collect::<Vec<_>>());
    }

    #[test]
    fn g1_static_thermal_solve_recovers_the_admitted_free_expansion() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let mat = material(120.0, 6.0);
        let alpha = 0.0125;
        let thermal = TetThermalStrainState::try_new(
            400.0,
            293.15,
            [alpha, alpha, alpha, 0.0, 0.0, 0.0],
            mat.material_state_identity,
            ContentHash([0x41; 32]),
            0.05,
        )
        .unwrap();
        // Remove exactly the six rigid modes without suppressing diagonal
        // free expansion: origin xyz, node-x yz, and node-y z.
        let problem = reference_problem(&nodes, &tets, &mat, &[0, 1, 2, 4, 5, 8]);
        let solution = with_cx(|cx| {
            problem.solve_thermal_displacement(
                TetThermalStrainField::Uniform(&thermal),
                TetStaticSolveConfig::default(),
                cx,
            )
        })
        .expect("constrained SPD thermal solve");
        let expected = [
            [0.0, 0.0, 0.0],
            [alpha, 0.0, 0.0],
            [0.0, alpha, 0.0],
            [0.0, 0.0, alpha],
        ];
        for (node, (actual, expected)) in solution
            .displacement_m()
            .iter()
            .zip(expected.iter())
            .enumerate()
        {
            for component in 0..3 {
                assert!(
                    (actual[component] - expected[component]).abs() < 1.0e-12,
                    "node {node} component {component}: {} != {}",
                    actual[component],
                    expected[component]
                );
            }
        }
        assert!(solution.true_relative_residual() < 1.0e-12);
        let stress = with_cx(|cx| {
            problem.recover_thermal_stress(
                &solution,
                TetThermalStrainField::Uniform(&thermal),
                ContentHash([0x51; 32]),
                0.05,
                cx,
            )
        })
        .unwrap();
        assert!(
            stress.elements()[0]
                .stress_mandel_pa()
                .iter()
                .all(|v| v.abs() < 1.0e-11)
        );
        assert!(
            stress
                .nodal_internal_force_n()
                .iter()
                .flatten()
                .all(|v| v.abs() < 1.0e-11)
        );
        assert!(stress.elastic_energy_j() < 1.0e-24);
        let deformed = with_cx(|cx| problem.update_geometry_from_displacement(&solution, 0.1, cx))
            .expect("non-inverting fixed-topology update");
        for (node, point) in deformed.nodes_m().iter().enumerate() {
            for component in 0..3 {
                assert!(
                    (point[component] - (nodes[node][component] + expected[node][component])).abs()
                        < 1.0e-12
                );
            }
        }
        assert_eq!(deformed.tetrahedra(), &tets);

        let mut inverted = solution.clone();
        inverted.displacement_m[3][2] = -2.0;
        assert!(matches!(
            with_cx(|cx| problem.update_geometry_from_displacement(&inverted, 3.0, cx)),
            Err(TetElasticError::InvertedTet { element: 0 })
        ));

        let free_body = reference_problem(&nodes, &tets, &mat, &[]);
        assert!(matches!(
            with_cx(|cx| free_body.solve_thermal_displacement(
                TetThermalStrainField::Uniform(&thermal),
                TetStaticSolveConfig::default(),
                cx,
            )),
            Err(TetElasticError::UnconstrainedRigidModes {
                component: 0,
                constrained_rank: 0,
            })
        ));
    }

    #[test]
    fn g1_g3_thermal_stress_reports_constraints_materials_reactions_and_frames() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [2.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let materials = [material(1200.0, 6.0), material(2400.0, 3.0)];
        let alphas = [0.01, 0.02];
        let thermal: Vec<_> = materials
            .iter()
            .zip(alphas)
            .map(|(mat, a)| {
                TetThermalStrainState::try_new(
                    400.0,
                    293.15,
                    [a, a, a, 0.0, 0.0, 0.0],
                    mat.material_state_identity,
                    ContentHash([0x41; 32]),
                    0.05,
                )
                .unwrap()
            })
            .collect();
        // Two supported bodies; only x displacement of each x-axis vertex is
        // free. Analytically sigma_xx=0, sigma_yy=sigma_zz=-E alpha/(1-nu),
        // and epsilon_xx=(1+nu)/(1-nu) alpha. Each volume is 1/6 m^3.
        let fixed: Vec<_> = (0..24).filter(|i| ![3, 15].contains(i)).collect();
        let problem = TetLinearElasticProblem {
            nodes_m: &nodes,
            tetrahedra: &tets,
            materials: TetMaterialField::PerElement(&materials),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
        };
        let solution = with_cx(|cx| {
            problem.solve_thermal_displacement(
                TetThermalStrainField::PerElement(&thermal),
                TetStaticSolveConfig::default(),
                cx,
            )
        })
        .unwrap();
        let frame = ContentHash([0x51; 32]);
        let report = with_cx(|cx| {
            problem.recover_thermal_stress(
                &solution,
                TetThermalStrainField::PerElement(&thermal),
                frame,
                0.05,
                cx,
            )
        })
        .unwrap();
        assert_eq!(
            report,
            with_cx(|cx| problem.recover_thermal_stress(
                &solution,
                TetThermalStrainField::PerElement(&thermal),
                frame,
                0.05,
                cx,
            ))
            .unwrap()
        );
        assert_eq!(report.solution_identity(), solution.identity());
        assert_eq!(report.mesh_frame(), frame);
        assert_eq!(report.elements().len(), 2);
        let mut expected_energy = 0.0;
        let q = fs_material::tensor::rotation([2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0], 0.73);
        let mut identities = BTreeSet::new();
        for element in 0..2 {
            let e = [1200.0, 2400.0][element];
            let a = alphas[element];
            let s = -e * a / 0.75;
            let expected_stress = [0.0, s, s, 0.0, 0.0, 0.0];
            let expected_elastic = [2.0 * a / 3.0, -a, -a, 0.0, 0.0, 0.0];
            assert!(
                (solution.displacement_m()[4 * element + 1][0] - 5.0 * a / 3.0).abs() < 1.0e-13
            );
            for i in 0..6 {
                assert!(
                    (report.elements()[element].stress_mandel_pa()[i] - expected_stress[i]).abs()
                        < 1.0e-10
                );
                assert!(
                    (report.elements()[element].elastic_strain_mandel()[i] - expected_elastic[i])
                        .abs()
                        < 1.0e-13
                );
            }
            let energy = -a * s / 6.0;
            expected_energy += energy;
            assert!((report.elements()[element].elastic_energy_j() - energy).abs() < 1.0e-13);
            let expected_forces = [
                [0.0, -s / 6.0, -s / 6.0],
                [0.0; 3],
                [0.0, s / 6.0, 0.0],
                [0.0, 0.0, s / 6.0],
            ];
            for node in 0..4 {
                for d in 0..3 {
                    assert!(
                        (report.nodal_internal_force_n()[4 * element + node][d]
                            - expected_forces[node][d])
                            .abs()
                            < 1.0e-11
                    );
                }
            }
            // Independent physical 3x3 tensor rotation, without Mandel-map code.
            let rotated = fs_material::tensor::rotate(&expected_stress, &q);
            let rotated_strain = fs_material::tensor::rotate(&expected_elastic, &q);
            for notation in [
                StressTensorNotation::Tensor,
                StressTensorNotation::Engineering,
                StressTensorNotation::Mandel,
            ] {
                for order in [
                    ElasticTensorOrder::XxYyZzXyYzZx,
                    ElasticTensorOrder::XxYyZzYzZxXy,
                ] {
                    let target = StressTensorBasis {
                        notation,
                        order,
                        frame: ContentHash([0x52; 32]),
                    };
                    let receipt = report.stress_in_frame(element, target, q).unwrap();
                    assert_eq!(receipt, report.stress_in_frame(element, target, q).unwrap());
                    assert_eq!(receipt.source_report_identity(), report.identity());
                    assert_eq!(receipt.source_frame(), frame);
                    assert_eq!(receipt.element(), element);
                    assert_eq!(receipt.target_basis(), target);
                    assert_eq!(receipt.source_to_target(), &q);
                    assert_eq!(
                        receipt.source_stress_mandel_pa(),
                        report.elements()[element].stress_mandel_pa()
                    );
                    assert!(identities.insert(receipt.identity()));
                    let indices = if order == ElasticTensorOrder::XxYyZzXyYzZx {
                        [0, 1, 2, 3, 4, 5]
                    } else {
                        [0, 1, 2, 4, 5, 3]
                    };
                    let stress_shear = if notation == StressTensorNotation::Mandel {
                        core::f64::consts::SQRT_2
                    } else {
                        1.0
                    };
                    let strain_shear = if notation == StressTensorNotation::Mandel {
                        core::f64::consts::SQRT_2
                    } else {
                        2.0
                    };
                    let mut work_density = 0.0;
                    for i in 0..6 {
                        let expected = rotated[indices[i]] * if i < 3 { 1.0 } else { stress_shear };
                        assert!(
                            (receipt.stress_pa()[i] - expected).abs() < 1.0e-10,
                            "{notation:?}/{order:?} {i}"
                        );
                        work_density += receipt.stress_pa()[i]
                            * rotated_strain[indices[i]]
                            * if i < 3 { 1.0 } else { strain_shear };
                    }
                    assert!((0.5 * work_density / 6.0 - energy).abs() < 1.0e-12);
                }
            }
        }
        assert_eq!(identities.len(), 12);
        assert!((report.elastic_energy_j() - expected_energy).abs() < 1.0e-13);
    }

    #[test]
    fn g1_fully_clamped_heating_has_zero_motion_and_physical_stress_reactions() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let mat = material(1200.0, 6.0);
        let fixed: Vec<_> = (0..12).collect();
        let problem = reference_problem(&nodes, &tets, &mat, &fixed);
        let thermal = TetThermalStrainState::try_new(
            400.0,
            293.15,
            [0.01, 0.01, 0.01, 0.0, 0.0, 0.0],
            mat.material_state_identity,
            ContentHash([0x41; 32]),
            0.05,
        )
        .unwrap();
        let field = TetThermalStrainField::Uniform(&thermal);
        // Empty modal/matrix spaces remain inadmissible. Thermal equilibrium
        // instead has prescribed u=0 and the nonzero clamp force -f_thermal.
        assert!(matches!(
            with_cx(|cx| problem.assemble(cx)),
            Err(TetElasticError::NoFreeDofs)
        ));
        let load = with_cx(|cx| problem.assemble_thermal_load(field, cx)).unwrap();
        assert!(load.free_dofs.is_empty());
        assert!(load.reduced_force_n.is_empty());
        let solution = with_cx(|cx| {
            problem.solve_thermal_displacement(field, TetStaticSolveConfig::default(), cx)
        })
        .unwrap();
        assert_eq!(solution.displacement_m(), &[[0.0; 3]; 4]);
        assert_eq!(solution.iterations(), 0);
        assert_eq!(solution.true_relative_residual(), 0.0);
        let report = with_cx(|cx| {
            problem.recover_thermal_stress(&solution, field, ContentHash([0x51; 32]), 0.05, cx)
        })
        .unwrap();
        // sigma=-E alpha/(1-2 nu) I=-24 I Pa, and V=1/6 m^3.
        for i in 0..6 {
            let expected = if i < 3 { -24.0 } else { 0.0 };
            assert!((report.elements()[0].stress_mandel_pa()[i] - expected).abs() < 1.0e-12);
        }
        let forces = [
            [4.0, 4.0, 4.0],
            [-4.0, 0.0, 0.0],
            [0.0, -4.0, 0.0],
            [0.0, 0.0, -4.0],
        ];
        for node in 0..4 {
            for d in 0..3 {
                let actual = report.nodal_internal_force_n()[node][d];
                assert!((actual - forces[node][d]).abs() < 1.0e-12);
                assert!((actual + load.full_equivalent_force_n[3 * node + d]).abs() < 1.0e-12);
            }
        }
        assert!((report.elastic_energy_j() - 0.06).abs() < 1.0e-14);
        assert_eq!(
            with_cx(|cx| problem.update_geometry_from_displacement(&solution, 0.1, cx))
                .unwrap()
                .nodes_m(),
            &nodes
        );
        // The empty solve still admits geometry, materials and thermal states.
        assert!(matches!(
            with_cx(|cx| problem.solve_thermal_displacement(
                TetThermalStrainField::PerElement(&[]),
                TetStaticSolveConfig::default(),
                cx
            )),
            Err(TetElasticError::ThermalStrainCountMismatch { .. })
        ));
        let collapsed = [[0.0; 3]; 4];
        let bad = reference_problem(&collapsed, &tets, &mat, &fixed);
        assert!(
            with_cx(|cx| bad.solve_thermal_displacement(
                field,
                TetStaticSolveConfig::default(),
                cx
            ))
            .is_err()
        );
    }

    #[test]
    fn g0_g4_thermal_stress_refuses_stale_inputs_bounds_frames_and_cancellation() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let mat = material(1200.0, 6.0);
        let fixed = [0, 1, 2, 4, 5, 8];
        let problem = reference_problem(&nodes, &tets, &mat, &fixed);
        let thermal = TetThermalStrainState::try_new(
            400.0,
            293.15,
            [0.01, 0.01, 0.01, 0.0, 0.0, 0.0],
            mat.material_state_identity,
            ContentHash([0x41; 32]),
            0.05,
        )
        .unwrap();
        let field = TetThermalStrainField::Uniform(&thermal);
        let solution = with_cx(|cx| {
            problem.solve_thermal_displacement(field, TetStaticSolveConfig::default(), cx)
        })
        .unwrap();
        let frame = ContentHash([0x51; 32]);
        let changed_mat = material(2400.0, 6.0);
        assert_eq!(
            mat.material_state_identity,
            changed_mat.material_state_identity
        );
        let changed = reference_problem(&nodes, &tets, &changed_mat, &fixed);
        // The old load identity alone cannot distinguish this reused material ID.
        assert_eq!(
            with_cx(|cx| problem.assemble_thermal_load(field, cx))
                .unwrap()
                .identity,
            with_cx(|cx| changed.assemble_thermal_load(field, cx))
                .unwrap()
                .identity
        );
        assert!(matches!(
            with_cx(|cx| changed.recover_thermal_stress(&solution, field, frame, 0.05, cx)),
            Err(TetElasticError::StressRecoveryInputMismatch)
        ));
        let hotter = TetThermalStrainState::try_new(
            410.0,
            293.15,
            [0.011, 0.011, 0.011, 0.0, 0.0, 0.0],
            mat.material_state_identity,
            ContentHash([0x41; 32]),
            0.05,
        )
        .unwrap();
        assert!(matches!(
            with_cx(|cx| problem.recover_thermal_stress(
                &solution,
                TetThermalStrainField::Uniform(&hotter),
                frame,
                0.05,
                cx
            )),
            Err(TetElasticError::StressRecoveryInputMismatch)
        ));
        let changed_constraints = reference_problem(&nodes, &tets, &mat, &[0, 1, 2, 4, 5, 6, 8]);
        assert!(matches!(
            with_cx(
                |cx| changed_constraints.recover_thermal_stress(&solution, field, frame, 0.05, cx)
            ),
            Err(TetElasticError::StressRecoveryInputMismatch)
        ));
        for bound in [0.0, f64::NAN, f64::INFINITY, 0.001] {
            assert!(matches!(
                with_cx(|cx| problem.recover_thermal_stress(&solution, field, frame, bound, cx)),
                Err(TetElasticError::InvalidStressObservation { .. })
            ));
        }
        assert!(matches!(
            with_cx(|cx| problem.recover_thermal_stress(
                &solution,
                field,
                ContentHash([0; 32]),
                0.05,
                cx
            )),
            Err(TetElasticError::InvalidStressObservation { .. })
        ));
        let report =
            with_cx(|cx| problem.recover_thermal_stress(&solution, field, frame, 0.05, cx))
                .unwrap();
        let identity_q = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let target = StressTensorBasis {
            notation: StressTensorNotation::Mandel,
            order: ElasticTensorOrder::XxYyZzXyYzZx,
            frame,
        };
        assert_eq!(
            report
                .stress_in_frame(0, target, identity_q)
                .unwrap()
                .stress_pa(),
            report.elements()[0].stress_mandel_pa()
        );
        assert!(matches!(
            report.stress_in_frame(1, target, identity_q),
            Err(TetElasticError::InvalidStressObservation { .. })
        ));
        assert!(
            report
                .stress_in_frame(
                    0,
                    StressTensorBasis {
                        frame: ContentHash([0; 32]),
                        ..target
                    },
                    identity_q
                )
                .is_err()
        );
        let q = fs_material::tensor::rotation([0.0, 0.0, 1.0], 0.25);
        assert!(matches!(
            report.stress_in_frame(0, target, q),
            Err(TetElasticError::InvalidStressObservation { .. })
        ));
        assert!(report.stress_in_frame(0, target, [[0.0; 3]; 3]).is_err());
        assert!(
            report
                .stress_in_frame(
                    0,
                    target,
                    [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
                )
                .is_err()
        );
        let gate = CancelGate::new();
        gate.request();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 3,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(matches!(
                problem.recover_thermal_stress(&solution, field, frame, 0.05, &cx),
                Err(TetElasticError::Cancelled)
            ));
        });
    }

    #[test]
    fn g0_thermal_strain_refuses_invalid_bounds_counts_and_stale_elastic_state() {
        let invalid = TetThermalStrainState::try_new(
            400.0,
            293.15,
            [0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
            [7; 32],
            ContentHash([0x32; 32]),
            0.05,
        );
        assert!(matches!(
            invalid,
            Err(TetElasticError::InvalidThermalStrain { .. })
        ));

        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let mat = material(120.0, 6.0);
        let problem = reference_problem(&nodes, &tets, &mat, &[]);
        assert!(matches!(
            with_cx(|cx| problem.assemble_thermal_load(TetThermalStrainField::PerElement(&[]), cx)),
            Err(TetElasticError::ThermalStrainCountMismatch {
                expected: 1,
                got: 0
            })
        ));
        let stale = TetThermalStrainState::try_new(
            400.0,
            293.15,
            [0.01, 0.01, 0.01, 0.0, 0.0, 0.0],
            [8; 32],
            ContentHash([0x33; 32]),
            0.05,
        )
        .unwrap();
        assert!(matches!(
            with_cx(|cx| problem.assemble_thermal_load(TetThermalStrainField::Uniform(&stale), cx)),
            Err(TetElasticError::ThermalElasticStateMismatch { element: 0 })
        ));

        let mut extreme_stiffness = [[0.0; 6]; 6];
        for (index, row) in extreme_stiffness.iter_mut().enumerate() {
            row[index] = 1.0e308;
        }
        let extreme = TetElasticMaterial::try_new_mandel(1.0, extreme_stiffness, [9; 32]).unwrap();
        let extreme_problem = reference_problem(&nodes, &tets, &extreme, &[]);
        let extreme_strain = TetThermalStrainState::try_new(
            400.0,
            293.15,
            [2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            extreme.material_state_identity,
            ContentHash([0x34; 32]),
            2.0,
        )
        .unwrap();
        assert!(matches!(
            with_cx(|cx| extreme_problem
                .assemble_thermal_load(TetThermalStrainField::Uniform(&extreme_strain), cx)),
            Err(TetElasticError::NonFiniteThermalLoad { element: 0, .. })
        ));
    }

    #[test]
    fn g0_oriented_orthotropy_rotates_energy_and_preserves_mass() {
        let law = OrthotropicElastic::new(
            [120.0, 60.0, 30.0],
            [0.2, 0.1, 0.15],
            [20.0, 15.0, 10.0],
            0.01,
        )
        .unwrap();
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let quarter_turn_about_z = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let principal = TetElasticMaterial::try_new_oriented_orthotropic(
            6.0, &law, identity, [8; 32], [0x18; 32],
        )
        .unwrap();
        let rotated = TetElasticMaterial::try_new_oriented_orthotropic(
            6.0,
            &law,
            quarter_turn_about_z,
            [9; 32],
            [0x19; 32],
        )
        .unwrap();
        assert!(
            (principal.stiffness_mandel_pa()[0][0] - rotated.stiffness_mandel_pa()[1][1]).abs()
                < 1.0e-12
        );
        assert!(
            (principal.stiffness_mandel_pa()[1][1] - rotated.stiffness_mandel_pa()[0][0]).abs()
                < 1.0e-12
        );

        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let principal_assembly =
            with_cx(|cx| reference_problem(&nodes, &tets, &principal, &[]).assemble(cx)).unwrap();
        let rotated_assembly =
            with_cx(|cx| reference_problem(&nodes, &tets, &rotated, &[]).assemble(cx)).unwrap();
        assert_eq!(principal_assembly.total_mass_kg, 1.0);
        assert_eq!(rotated_assembly.total_mass_kg, 1.0);

        let alpha = 0.125;
        let x_extension: Vec<f64> = nodes
            .iter()
            .flat_map(|point| [alpha * point[0], 0.0, 0.0])
            .collect();
        let principal_energy =
            0.5 * quadratic(&principal_assembly.stiffness.to_dense(), &x_extension);
        let rotated_energy = 0.5 * quadratic(&rotated_assembly.stiffness.to_dense(), &x_extension);
        assert!(principal_energy > rotated_energy);
        assert!(
            (principal_energy / rotated_energy
                - principal.stiffness_mandel_pa()[0][0] / rotated.stiffness_mandel_pa()[0][0])
                .abs()
                < 1.0e-12
        );
    }

    fn coupled_engineering_stiffness() -> [[f64; 6]; 6] {
        // Synthetic, strictly diagonally dominant SPD law with normal/shear
        // coupling and signed off-diagonal coefficients; no material preset.
        [
            [120.0, 20.0, 10.0, 4.0, -3.0, 2.0],
            [20.0, 110.0, 15.0, -2.0, 5.0, -1.0],
            [10.0, 15.0, 100.0, 3.0, -4.0, 6.0],
            [4.0, -2.0, 3.0, 50.0, 2.0, -1.0],
            [-3.0, 5.0, -4.0, 2.0, 40.0, 3.0],
            [2.0, -1.0, 6.0, -1.0, 3.0, 30.0],
        ]
    }

    // Independent oracle: expand C_ijkl directly from engineering coefficients
    // and contract all four rotation indices. No Mandel rotation helper or
    // production notation converter is used here.
    fn four_index_rotated_stiffness(q: [[f64; 3]; 3]) -> [[[[f64; 3]; 3]; 3]; 3] {
        let engineering = coupled_engineering_stiffness();
        let component = |i, j| match (i, j) {
            (0, 0) => 0,
            (1, 1) => 1,
            (2, 2) => 2,
            (0, 1) | (1, 0) => 3,
            (1, 2) | (2, 1) => 4,
            (0, 2) | (2, 0) => 5,
            _ => unreachable!(),
        };
        let mut out = [[[[0.0; 3]; 3]; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        for p in 0..3 {
                            for r in 0..3 {
                                for s in 0..3 {
                                    for t in 0..3 {
                                        out[i][j][k][l] += q[i][p]
                                            * q[j][r]
                                            * q[k][s]
                                            * q[l][t]
                                            * engineering[component(p, r)][component(s, t)];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        out
    }

    #[test]
    fn g3_full_tensor_notations_orders_and_four_index_rotation_agree() {
        let q = fs_material::tensor::rotation([2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0], 0.73);
        let expected = four_index_rotated_stiffness(q);
        let pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (2, 0)];
        let engineering = coupled_engineering_stiffness();
        let mut identities = BTreeSet::new();
        for notation in [
            ElasticTensorNotation::Engineering,
            ElasticTensorNotation::Tensor,
            ElasticTensorNotation::Mandel,
        ] {
            for order in [
                ElasticTensorOrder::XxYyZzXyYzZx,
                ElasticTensorOrder::XxYyZzYzZxXy,
            ] {
                let permutation = if order == ElasticTensorOrder::XxYyZzXyYzZx {
                    [0, 1, 2, 3, 4, 5]
                } else {
                    [0, 1, 2, 4, 5, 3]
                };
                let source = core::array::from_fn(|i| {
                    core::array::from_fn(|j| {
                        let (a, b) = (permutation[i], permutation[j]);
                        let scale = match notation {
                            ElasticTensorNotation::Engineering => 1.0,
                            ElasticTensorNotation::Tensor => {
                                if b >= 3 {
                                    2.0
                                } else {
                                    1.0
                                }
                            }
                            ElasticTensorNotation::Mandel => match (a >= 3, b >= 3) {
                                (false, false) => 1.0,
                                (true, true) => 2.0,
                                _ => core::f64::consts::SQRT_2,
                            },
                        };
                        engineering[a][b] * scale
                    })
                });
                let basis = ElasticTensorBasis {
                    notation,
                    order,
                    frame: ContentHash([1; 32]),
                };
                let receipt = TetElasticMaterial::try_new_oriented_tensor(
                    6.0,
                    source,
                    basis,
                    ContentHash([2; 32]),
                    q,
                    ContentHash([3; 32]),
                )
                .unwrap();
                assert_eq!(receipt.source_stiffness_pa(), &source);
                assert_eq!(receipt.source_basis(), basis);
                assert_eq!(receipt.target_frame(), ContentHash([2; 32]));
                assert_eq!(receipt.source_to_target(), &q);
                assert_eq!(receipt.source_material_identity(), ContentHash([3; 32]));
                assert!(identities.insert(receipt.material().material_state_identity));
                for (a, &(i, j)) in pairs.iter().enumerate() {
                    for (b, &(k, l)) in pairs.iter().enumerate() {
                        let scale = match (a >= 3, b >= 3) {
                            (false, false) => 1.0,
                            (true, true) => 2.0,
                            _ => core::f64::consts::SQRT_2,
                        };
                        let oracle = expected[i][j][k][l] * scale;
                        assert!(
                            (receipt.material().stiffness_mandel_pa()[a][b] - oracle).abs()
                                < 2.0e-12,
                            "{notation:?} {order:?} coefficient {a},{b}"
                        );
                    }
                }
                assert_eq!(
                    receipt,
                    TetElasticMaterial::try_new_oriented_tensor(
                        6.0,
                        source,
                        basis,
                        ContentHash([2; 32]),
                        q,
                        ContentHash([3; 32]),
                    )
                    .unwrap()
                );
            }
        }
        assert_eq!(
            identities.len(),
            6,
            "equivalent responses retain distinct source conventions"
        );
    }

    #[test]
    fn g1_full_tensor_drives_tet_forces_and_energy_in_the_target_frame() {
        use fs_matdb::{
            ClaimSet, ElasticTensorComponent, ElasticTensorSymmetry, InterpolationPolicy,
            MaterialStateId, NormalizedMaterialCardPack, NormalizedPack, PropertyClaim,
            PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
        };
        use fs_material::state_point::{
            MaterialPropertySelection, resolve_elastic_tensor_state_point,
        };
        let q = fs_material::tensor::rotation([2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0], 0.73);
        let receipt = TetElasticMaterial::try_new_oriented_tensor(
            6.0,
            coupled_engineering_stiffness(),
            ElasticTensorBasis {
                notation: ElasticTensorNotation::Engineering,
                order: ElasticTensorOrder::XxYyZzXyYzZx,
                frame: ContentHash([1; 32]),
            },
            ContentHash([2; 32]),
            q,
            ContentHash([3; 32]),
        )
        .unwrap();
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let assembly =
            with_cx(|cx| reference_problem(&nodes, &tets, receipt.material(), &[]).assemble(cx))
                .unwrap();
        assert_eq!(assembly.total_mass_kg, 1.0);
        // Independently exercise the portable-card -> 37 selected receipts ->
        // complete tensor admission -> same numerical element path.
        let basis = receipt.source_basis();
        let keys: [[PropertyKey; 6]; 6] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                PropertyKey::new("stiffness", fs_qty::Pressure::DIMS)
                    .with_elastic_component(
                        ElasticTensorComponent::new(
                            basis,
                            ElasticTensorSymmetry::MajorMinor,
                            ContentHash([4; 32]),
                            i as u8,
                            j as u8,
                        )
                        .unwrap(),
                    )
                    .unwrap()
            })
        });
        let mut claims = ClaimSet::new();
        let provenance = Provenance {
            source: "synthetic full-tensor software fixture".into(),
            license: "synthetic test data".into(),
            artifact: Some(ContentHash([5; 32])),
        };
        let observation = claims
            .register_observation(fs_matdb::ObservationDataset {
                specimen: "synthetic tensor specimen".into(),
                method: "synthetic software fixture".into(),
                artifact: ContentHash([5; 32]),
                caveats: "no empirical or calibration claim".into(),
                provenance: provenance.clone(),
            })
            .unwrap();
        let coefficients = coupled_engineering_stiffness();
        let mut values: Vec<_> = keys
            .iter()
            .flatten()
            .enumerate()
            .map(|(i, key)| (key.clone(), coefficients[i / 6][i % 6]))
            .collect();
        values.push((PropertyKey::new("density", fs_qty::Density::DIMS), 6.0));
        for (key, value) in values {
            claims
                .insert_claim(PropertyClaim {
                    value: PropertyValue::Scalar {
                        value,
                        dims: key.dims(),
                    },
                    key,
                    validity: fs_evidence::ValidityDomain::unconstrained().with("T", 300.0, 300.0),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::TabulatedOnly,
                    observations: vec![observation],
                    provenance: provenance.clone(),
                })
                .unwrap();
        }
        let pack = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: "synthetic".into(),
                phase: "solid".into(),
                process: "numerical fixture".into(),
                revision: 0,
            },
            NormalizedPack::new(
                "synthetic",
                "fixture-v5",
                ContentHash([5; 32]),
                "synthetic redistribution permitted",
                claims,
                Vec::new(),
                Vec::new(),
            )
            .unwrap(),
        )
        .unwrap();
        let portable =
            NormalizedMaterialCardPack::from_bytes_verified(pack.content_hash(), &pack.to_bytes())
                .unwrap();
        let state = resolve_elastic_tensor_state_point(
            portable.card(),
            &QueryPoint::new().with("T", 300.0).unwrap(),
            &keys,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert_eq!(state.resolved().properties().len(), 37);
        assert!(state.resolved().property("stiffness").is_none());
        for (index, key) in keys.iter().flatten().enumerate() {
            assert_eq!(
                state.resolved().property_by_key(key).unwrap().value_si(),
                coefficients[index / 6][index % 6]
            );
        }
        let mut pins: Vec<_> = state
            .resolved()
            .properties()
            .iter()
            .map(|property| {
                (
                    property.requirement().name().to_owned(),
                    property.answer().receipt.selected,
                )
            })
            .collect();
        pins.reverse();
        let pinned = resolve_elastic_tensor_state_point(
            portable.card(),
            &QueryPoint::new().with("T", 300.0).unwrap(),
            &keys,
            MaterialPropertySelection::PinnedByProperty(pins.clone()),
        )
        .unwrap();
        assert_eq!(pinned.stiffness_pa(), &coefficients);
        let stiffness_pin = pins
            .iter()
            .position(|(name, _)| name == "stiffness")
            .unwrap();
        let another = pins
            .iter()
            .enumerate()
            .find(|(index, (name, _))| *index != stiffness_pin && name == "stiffness")
            .unwrap()
            .0;
        pins[another].1 = pins[stiffness_pin].1;
        assert!(matches!(
            resolve_elastic_tensor_state_point(
                portable.card(),
                &QueryPoint::new().with("T", 300.0).unwrap(),
                &keys,
                MaterialPropertySelection::PinnedByProperty(pins),
            ),
            Err(fs_material::state_point::MaterialStatePointError::InvalidSelectionPlan)
        ));
        let from_card =
            TetElasticMaterial::from_resolved_elastic_tensor(&state, ContentHash([2; 32]), q)
                .unwrap();
        assert_eq!(
            from_card.source_material_identity(),
            state.resolved().identity()
        );
        let card_assembly =
            with_cx(|cx| reference_problem(&nodes, &tets, from_card.material(), &[]).assemble(cx))
                .unwrap();
        assert_eq!(
            card_assembly.stiffness.to_dense(),
            assembly.stiffness.to_dense()
        );
        assert_eq!(card_assembly.total_mass_kg, assembly.total_mass_kg);
        let epsilon = [
            [0.001, 0.0002, -0.0003],
            [0.0002, -0.0008, 0.0004],
            [-0.0003, 0.0004, 0.0005],
        ];
        let c = four_index_rotated_stiffness(q);
        let mut stress = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        stress[i][j] += c[i][j][k][l] * epsilon[k][l];
                    }
                }
            }
        }
        let displacement: Vec<f64> = nodes
            .iter()
            .flat_map(|x| epsilon.map(|row| row.iter().zip(x).map(|(a, b)| a * b).sum::<f64>()))
            .collect();
        let k = assembly.stiffness.to_dense();
        let gradients = [
            [-1.0, -1.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        for (node, gradient) in gradients.iter().enumerate() {
            for (a, row) in stress.iter().enumerate() {
                let actual: f64 = (0..12)
                    .map(|dof| k[(3 * node + a) * 12 + dof] * displacement[dof])
                    .sum();
                let expected: f64 = row.iter().zip(gradient).map(|(s, g)| s * g / 6.0).sum();
                assert!((actual - expected).abs() < 2.0e-14);
            }
        }
        let expected_energy: f64 = epsilon
            .iter()
            .flatten()
            .zip(stress.iter().flatten())
            .map(|(e, s)| e * s / 12.0)
            .sum();
        assert!((0.5 * quadratic(&k, &displacement) - expected_energy).abs() < 2.0e-16);
    }

    #[test]
    fn g0_full_tensor_admission_checks_whole_law_frames_and_transform_identity() {
        let source = coupled_engineering_stiffness();
        let basis = ElasticTensorBasis {
            notation: ElasticTensorNotation::Engineering,
            order: ElasticTensorOrder::XxYyZzXyYzZx,
            frame: ContentHash([1; 32]),
        };
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let admit = |matrix, basis, frame, rotation, state| {
            TetElasticMaterial::try_new_oriented_tensor(6.0, matrix, basis, frame, rotation, state)
        };
        let receipt = admit(
            source,
            basis,
            ContentHash([2; 32]),
            identity,
            ContentHash([3; 32]),
        )
        .unwrap();
        // Misdeclaring engineering data as tensor shear breaks reciprocity;
        // do not accept it merely because its raw array is symmetric.
        assert!(
            admit(
                source,
                ElasticTensorBasis {
                    notation: ElasticTensorNotation::Tensor,
                    ..basis
                },
                ContentHash([2; 32]),
                identity,
                ContentHash([3; 32])
            )
            .is_err()
        );
        let mut asymmetric = source;
        asymmetric[0][3] += 1.0;
        let mut unstable = source;
        unstable[0][3] = 1000.0;
        unstable[3][0] = 1000.0; // positive diagonal blocks, indefinite whole law
        let mut nonfinite = source;
        nonfinite[0][0] = f64::NAN;
        let mut overflow = source;
        overflow[3][3] = f64::MAX;
        for matrix in [asymmetric, unstable, nonfinite, overflow] {
            assert!(
                admit(
                    matrix,
                    basis,
                    ContentHash([2; 32]),
                    identity,
                    ContentHash([3; 32])
                )
                .is_err()
            );
        }
        for (source_frame, target_frame, state) in [(0, 2, 3), (1, 0, 3), (1, 2, 0)] {
            assert!(
                admit(
                    source,
                    ElasticTensorBasis {
                        frame: ContentHash([source_frame; 32]),
                        ..basis
                    },
                    ContentHash([target_frame; 32]),
                    identity,
                    ContentHash([state; 32])
                )
                .is_err()
            );
        }
        let quarter_turn = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        assert!(
            admit(
                source,
                basis,
                basis.frame,
                quarter_turn,
                ContentHash([3; 32])
            )
            .is_err()
        );
        for q in [
            [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        ] {
            assert!(admit(source, basis, ContentHash([2; 32]), q, ContentHash([3; 32])).is_err());
        }
        let mut changed = source;
        changed[0][0] += 1.0;
        for different in [
            admit(
                changed,
                basis,
                ContentHash([2; 32]),
                identity,
                ContentHash([3; 32]),
            )
            .unwrap(),
            admit(
                source,
                basis,
                ContentHash([4; 32]),
                identity,
                ContentHash([3; 32]),
            )
            .unwrap(),
            admit(
                source,
                ElasticTensorBasis {
                    frame: ContentHash([4; 32]),
                    ..basis
                },
                ContentHash([2; 32]),
                identity,
                ContentHash([3; 32]),
            )
            .unwrap(),
            admit(
                source,
                basis,
                ContentHash([2; 32]),
                identity,
                ContentHash([4; 32]),
            )
            .unwrap(),
            admit(
                source,
                basis,
                ContentHash([2; 32]),
                quarter_turn,
                ContentHash([3; 32]),
            )
            .unwrap(),
        ] {
            assert_ne!(
                receipt.material().material_state_identity,
                different.material().material_state_identity
            );
        }
    }

    #[test]
    fn g0_constitutive_admission_refuses_asymmetry_indefiniteness_and_bad_orientation() {
        let mut asymmetric = [[0.0; 6]; 6];
        for (index, row) in asymmetric.iter_mut().enumerate() {
            row[index] = 1.0;
        }
        asymmetric[0][1] = 0.25;
        assert_eq!(
            TetElasticMaterial::try_new_mandel(1.0, asymmetric, [1; 32]).unwrap_err(),
            TetElasticError::InvalidMaterial {
                what: "stiffness_mandel_pa must be exactly symmetric"
            }
        );

        let mut indefinite = [[0.0; 6]; 6];
        for (index, row) in indefinite.iter_mut().enumerate() {
            row[index] = 1.0;
        }
        indefinite[5][5] = -1.0;
        assert_eq!(
            TetElasticMaterial::try_new_mandel(1.0, indefinite, [1; 32]).unwrap_err(),
            TetElasticError::InvalidMaterial {
                what: "stiffness_mandel_pa must be positive definite"
            }
        );

        let law = OrthotropicElastic::new([3.0, 2.0, 1.0], [0.1, 0.1, 0.1], [0.8, 0.7, 0.6], 0.01)
            .unwrap();
        let reflection = [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert_eq!(
            TetElasticMaterial::try_new_oriented_orthotropic(
                1.0, &law, reflection, [1; 32], [2; 32],
            )
            .unwrap_err(),
            TetElasticError::InvalidMaterial {
                what: "principal_to_world must be a proper rotation"
            }
        );
    }

    #[test]
    fn g1_frequency_scaling_follows_geometry_modulus_and_density() {
        let base = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        // Minimal exact rigid-body suppression: node 0 xyz, node 1 yz,
        // node 2 z. Six elastic coordinates remain.
        let fixed = [0, 1, 2, 4, 5, 8];
        let lowest = |nodes: &[[f64; 3]], e: f64, rho: f64| {
            let mat = material(e, rho);
            let a =
                with_cx(|cx| reference_problem(nodes, &tets, &mat, &fixed).assemble(cx)).unwrap();
            let modes = eigh_gen_dense(
                &a.stiffness.to_dense(),
                &a.mass.to_dense(),
                a.free_dofs.len(),
            )
            .unwrap();
            modes[0].lambda
        };
        let lambda = lowest(&base, 120.0, 6.0);
        let scaled = base.map(|p| p.map(|x| 2.0 * x));
        let lambda_size = lowest(&scaled, 120.0, 6.0);
        let lambda_e = lowest(&base, 480.0, 6.0);
        let lambda_rho = lowest(&base, 120.0, 24.0);
        assert!((lambda_size / lambda - 0.25).abs() < 1.0e-10);
        assert!((lambda_e / lambda - 4.0).abs() < 1.0e-10);
        assert!((lambda_rho / lambda - 0.25).abs() < 1.0e-10);
    }

    #[test]
    fn malformed_topology_materials_and_quality_refuse() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0e-16, 1.0e-16, 0.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let mat = material(1.0, 1.0);
        let error =
            with_cx(|cx| reference_problem(&nodes, &tets, &mat, &[]).assemble(cx)).unwrap_err();
        assert!(matches!(error, TetElasticError::DegenerateTet { .. }));

        let bad = TetElasticMaterial::try_new(1.0, 1.0, 0.25, [0; 32]).unwrap_err();
        assert_eq!(
            bad,
            TetElasticError::InvalidMaterial {
                what: "material_state_identity must not be the zero identity"
            }
        );
    }

    #[test]
    fn g1_rounded_cylinder_geometry_composes_with_elastic_mass() {
        let mesh = with_cx(|cx| {
            rounded_cylinder_tet_mesh(
                RoundedCylinderMeshSpec {
                    outer_radius_m: 0.038,
                    thickness_m: 0.006,
                    fillet_radius_m: 0.001,
                    core_radial_segments: 3,
                    fillet_radial_segments: 2,
                    azimuthal_segments: 12,
                    axial_segments: 2,
                    maximum_vertices: 10_000,
                    maximum_tetrahedra: 50_000,
                },
                cx,
            )
        })
        .unwrap();
        let density = 7_800.0;
        let mat = material(200.0e9, density);
        let assembly = with_cx(|cx| {
            TetLinearElasticProblem {
                nodes_m: &mesh.nodes_m,
                tetrahedra: &mesh.tetrahedra,
                materials: TetMaterialField::Uniform(&mat),
                fixed_dofs: &[],
                budget: TetAssemblyBudget::standard(),
            }
            .assemble(cx)
        })
        .unwrap();
        let geometric_volume: f64 = mesh
            .tetrahedra
            .iter()
            .map(|tet| {
                let [a, b, c, d] = tet.map(|vertex| mesh.nodes_m[vertex]);
                determinant_columns([sub(b, a), sub(c, a), sub(d, a)]).abs() / 6.0
            })
            .sum();
        assert!((assembly.total_mass_kg - density * geometric_volume).abs() < 1.0e-12);
        assert_eq!(assembly.free_dofs.len(), 3 * mesh.nodes_m.len());
        assert!(!mesh.boundary.triangles.is_empty());
    }

    fn quadratic(matrix: &[f64], vector: &[f64]) -> f64 {
        let n = vector.len();
        (0..n)
            .map(|row| {
                (0..n)
                    .map(|col| vector[row] * matrix[row * n + col] * vector[col])
                    .sum::<f64>()
            })
            .sum()
    }

    fn determinant_columns(columns: [[f64; 3]; 3]) -> f64 {
        columns[0][0] * (columns[1][1] * columns[2][2] - columns[1][2] * columns[2][1])
            - columns[1][0] * (columns[0][1] * columns[2][2] - columns[0][2] * columns[2][1])
            + columns[2][0] * (columns[0][1] * columns[1][2] - columns[0][2] * columns[1][1])
    }
}
