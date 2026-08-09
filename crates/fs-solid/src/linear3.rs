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
    IntegratedIsotropicThermalExpansion, IsotropicElasticStatePoint, IsotropicSolidStatePoint,
    OrthotropicElasticStatePoint,
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
        let transform = mandel_rotation(principal_to_world);
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
        // The operation tree above is symmetric mathematically, but separate
        // floating-point accumulation paths can differ by a final bit. Use one
        // deterministic triangle as the exact stored constitutive tensor.
        for row in 0..6 {
            for column in 0..row {
                world[row][column] = world[column][row];
            }
        }
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

impl TetThermalStrainState {
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
/// `K u = equivalent_force`.  The full vector retains reactions at constrained
/// DOFs; `reduced_force_n` follows the elastic assembly's free-DOF order.
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
    /// Every displacement DOF was constrained.
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
        let assembly = self.assemble(cx)?;
        validate_rigid_mode_constraints(self.nodes_m, self.tetrahedra, self.fixed_dofs)?;
        let load = self.assemble_thermal_load(thermal_strains, cx)?;
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
        let mut flat = vec![0.0; self.nodes_m.len() * 3];
        for (&full_dof, &value) in load.free_dofs.iter().zip(&state.x) {
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
        hash_usize(&mut identity, state.iters)?;
        identity.update(&true_relative_residual.to_bits().to_le_bytes());
        for displacement in &displacement_m {
            for value in displacement {
                identity.update(&value.to_bits().to_le_bytes());
            }
        }
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        Ok(TetThermalDisplacementSolution {
            displacement_m,
            iterations: state.iters,
            true_relative_residual,
            identity: identity.finalize(),
            reference_geometry_identity: tet_geometry_identity(self.nodes_m, self.tetrahedra)?,
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
    if free_dofs.is_empty() {
        return Err(TetElasticError::NoFreeDofs);
    }
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
