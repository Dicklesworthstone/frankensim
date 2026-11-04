//! Small-strain tensor-Prony solids (IRP M1.3-W2a).
//!
//! The material owner advances every quadrature point by its exact linear
//! strain-ramp update. The element owner applies the same work-conjugate B/Bᵀ
//! pair to displacement and stress, and sums the material energy densities
//! with positive volume weights. Thus nonnegative constitutive heat stays
//! nonnegative after assembly. A pointwise relative quadratic-form bound on
//! C′ and C″ also bounds the corresponding assembled quadratic forms: integrate
//! the pointwise inequalities with these same positive weights. This statement
//! does not certify modal truncation or solve the dispersive eigenproblem.
//!
//! Histories borrow their immutable prepared system, so neither geometry nor
//! laws can change underneath internal variables. Updates are transactional,
//! including cancellation and nonfinite-result refusals. Prescribed nodal
//! histories and load-controlled quasistatic creep share the same material
//! histories. Dynamics couples their exact work-mean force to consistent-mass
//! midpoint inertia on the original spatial operators.
//!
//! HEX8 also admits the existing volumetric B-bar strain projection. Every
//! equilibrium, memory, Kelvin, mean-force and harmonic operator uses those
//! same columns, while compatible strain independently respects the source
//! domain. Inertia keeps the elastic owner's consistent reference quadrature.
//! This is small-strain, elementwise-material B-bar, not exact incompressibility
//! or an evolving-temperature constitutive model.
//! Prescribed thermal expansion at the original quadrature is a separate work input. It acts
//! on these same histories with fixed constitutive coefficients and reference
//! mass; neither a reciprocal caloric law nor coefficient transport is inferred.

use super::hex::{
    HexConnectivity, HexElasticProblem, HexFamily, HexPoint, HexQuadrature, hex20_point,
    hex27_point,
};
use super::hex8::Hex8ElasticProblem;
use super::tet10::{QuadraticTetProblem, TetQuadrature, tet10_point};
use super::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticError, TetElasticMaterial, TetGeometry,
    TetLinearElasticProblem, TetMaterialField,
};
use crate::linear::Formulation;
use fs_exec::Cx;
use fs_feec::highorder::simplex::duffy_quadrature;
use fs_material::visco::ViscoError;
use fs_material::visco::tensor::{Mandel6, TensorProny, TensorPronyState};
use fs_sparse::{Coo, Csr};

mod hex8;
mod project;
mod thermal;
pub use project::{ProjectedViscoSolid, SolidProjectionBudget};
mod quasistatic;
pub use quasistatic::{ViscoQuasistaticOperator, ViscoQuasistaticStep};
mod dynamics;
pub use dynamics::{ViscoDynamicConfig, ViscoDynamicOperator, ViscoDynamicState, ViscoDynamicStep};
mod harmonic;
pub use harmonic::{ViscoHarmonicError, ViscoHarmonicOperator, ViscoHarmonicResponse};
mod poles;
pub use poles::{
    PreparedViscoPoles, ViscoJointHarmonicBudget, ViscoJointHarmonicError,
    ViscoJointHarmonicFactor, ViscoJointHarmonicLoad, ViscoJointHarmonicResponse,
    ViscoJointHarmonicStateSpaceLoad, ViscoMemorySite, ViscoPole, ViscoPoleBudget,
    ViscoPoleDescriptorLoad, ViscoPoleDynamicBudget, ViscoPoleDynamicLoad,
    ViscoPoleDynamicResponse, ViscoPoleError, ViscoPolePortBudget, ViscoPolePortLoad,
    ViscoPolePortResponse, ViscoPoleSpectrum, ViscoPoleStateSpaceLoad, ViscoPoleStateSpaceResponse,
};

/// Geometry and explicit quadrature; material frames must already be resolved.
#[derive(Clone, Copy)]
pub enum SolidMeshRef<'a> {
    /// Constant-strain affine tetrahedra.
    Tet4 {
        /// Complete node positions [m].
        nodes_m: &'a [[f64; 3]],
        /// Four corner IDs per element.
        elements: &'a [[usize; 4]],
    },
    /// Isoparametric quadratic tetrahedra, in canonical TET10 order.
    Tet10 {
        /// Complete corner and shared midside positions [m].
        nodes_m: &'a [[f64; 3]],
        /// Ten node IDs per element, in canonical TET10 order.
        elements: &'a [[usize; 10]],
        /// Required positive-weight stiffness and consistent-mass rules.
        quadrature: TetQuadrature,
    },
    /// Trilinear eight-node bricks using the native elastic B-bar owner.
    Hex8 {
        /// Complete corner positions [m].
        nodes_m: &'a [[f64; 3]],
        /// Canonical corners: four around z-, then corresponding corners at z+.
        elements: &'a [[usize; 8]],
        /// Compatible strain or reference-volume-mean volumetric B-bar.
        formulation: Formulation,
        /// Positive tensor Gauss rule, in 2..=12. Inertia uses at least three
        /// points per axis; every stiffness point retains its own history.
        stiffness_points_per_dir: usize,
    },
    /// Fully integrated isoparametric HEX20 or HEX27, in canonical node order.
    /// The original elastic owner admits topology and sampled Jacobians.
    Hex {
        /// Complete corner, edge and (for HEX27) face/interior positions [m].
        nodes_m: &'a [[f64; 3]],
        /// Existing hexahedral family and connectivity, without remeshing.
        elements: HexConnectivity<'a>,
        /// Positive tensor-product rule, with at least three points per axis.
        quadrature: HexQuadrature,
    },
}

/// One law already rotated into the mesh frame; SI density and caller identity.
#[derive(Clone)]
pub struct TetViscoMaterial {
    /// Constitutive owner, including its applicability band.
    pub law: TensorProny,
    /// Density [kg/m³].
    pub density_kg_m3: f64,
    /// Resolved material or explicitly named test-fixture identity.
    pub material_state_identity: [u8; 32],
}

impl TetViscoMaterial {
    /// Original relaxed tensor, reference density and coefficient-source identity.
    ///
    /// Pass this material to the existing `ThermalExpansionModel::resolve` to
    /// bind imposed expansion to these frozen Prony coefficients. It does not
    /// requery stiffness or density at the expansion's target temperature.
    pub fn equilibrium_material(&self) -> Result<TetElasticMaterial, TetElasticError> {
        TetElasticMaterial::try_new_mandel(
            self.density_kg_m3,
            *self.law.c_eq(),
            self.material_state_identity,
        )
    }

    /// Lower one same-card constitutive/inertial source into the original
    /// tetrahedral owner. Every tensor uses the explicit material-to-mesh
    /// rotation; density, branch rates and applicability are unchanged. Retain
    /// `source` alongside the specimen when its full receipts are needed.
    pub fn from_source(
        source: &fs_material::visco_tensor::source::ResolvedTensorMaterial,
        material_to_mesh: [[f64; 3]; 3],
    ) -> Result<Self, fs_material::visco_tensor::source::TensorSourceError> {
        Ok(Self {
            law: source.rotated_law(material_to_mesh)?,
            density_kg_m3: source.density_kg_m3(),
            material_state_identity: source.state().identity().0,
        })
    }
}

/// Uniform or explicitly per-element constitutive data.
#[derive(Clone, Copy)]
pub enum TetViscoField<'a> {
    /// One immutable law for the mesh.
    Uniform(&'a TetViscoMaterial),
    /// One law per element, in connectivity order.
    PerElement(&'a [TetViscoMaterial]),
}

impl TetViscoField<'_> {
    fn at(&self, element: usize) -> &TetViscoMaterial {
        match self {
            Self::Uniform(m) => m,
            Self::PerElement(m) => &m[element],
        }
    }
}

/// Mandatory history allocation cap, counted as six scalars per point/branch.
/// Optional thermal endpoints retain one state per admitted element or stiffness point;
/// those source snapshots are separate from this branch-strain scalar count.
#[derive(Debug, Clone, Copy)]
pub struct ViscoStateBudget {
    /// Maximum total internal strain scalars, not a byte count.
    pub max_internal_variables: usize,
}

/// Named admission or compute failure; no partial state/operator escapes.
#[derive(Debug, Clone, PartialEq)]
pub enum ViscoElementError {
    /// Existing mesh, material-density or assembly admission failed.
    Elastic(TetElasticError),
    /// Constitutive owner refused (including out-of-band queries).
    Material(ViscoError),
    /// The original work-conjugate tensor-material projection refused.
    Projection(fs_material::visco_tensor::TensorPronyError),
    /// Required internal strain scalars exceeded the explicit cap.
    InternalStateBudget {
        /// Internal scalar strains required by the mesh and law field.
        needed: usize,
        /// Caller-admitted scalar count.
        cap: usize,
    },
    /// Wrong displacement size, nonfinite values, or zero/negative duration.
    InvalidHistory {
        /// Violated history or observation invariant.
        what: &'static str,
    },
    /// Nonfinite stress, force or energy; candidate state was discarded.
    NonFiniteResult,
    /// Independently evaluated inertial/material work balance missed its gate.
    DynamicEnergy {
        /// Signed work minus kinetic/storage change and physical heat [J].
        residual_j: f64,
        /// Explicit absolute plus relative numerical allowance [J].
        allowance_j: f64,
    },
    /// Explicit mode-count, matrix-size, iteration or frequency controls failed.
    InvalidModalBudget,
    /// The existing modal eigensolver refused the pencil.
    Modal(fs_modal::ModalError),
    /// The caller's search interval does not bracket the selected ordered mode.
    UnbracketedMode {
        /// Zero-based ascending mode index.
        mode: usize,
    },
    /// Frequency iteration exhausted its budget or hit algebraic resolution.
    ModalConvergence {
        /// Zero-based ascending mode index.
        mode: usize,
    },
    /// The algebraic eigensolve cannot resolve the requested root accuracy.
    EigenResolution {
        /// Zero-based ascending mode index.
        mode: usize,
        /// Frequency at which the eigen residual exceeded the budget [rad/s].
        omega_rad_s: f64,
        /// Certified eigenvalue radius divided by omega².
        relative_radius: f64,
        /// Maximum admitted relative eigenvalue radius.
        maximum: f64,
    },
}

impl core::fmt::Display for ViscoElementError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-SOLID-VISCO: {self:?}")
    }
}
impl std::error::Error for ViscoElementError {}
impl From<TetElasticError> for ViscoElementError {
    fn from(e: TetElasticError) -> Self {
        Self::Elastic(e)
    }
}
impl From<ViscoError> for ViscoElementError {
    fn from(e: ViscoError) -> Self {
        Self::Material(e)
    }
}

struct IntegrationPoint {
    b: Vec<[f64; 6]>,
    // B-bar changes only the normal rows. Retain 3N raw gradient scalars,
    // instead of a second 6-by-3N matrix, for independent strain admission.
    compatible_gradients: Option<Vec<[f64; 3]>>,
    volume_m3: f64,
}

impl IntegrationPoint {
    fn strain(&self, dofs: &[usize], displacement: &[f64]) -> [f64; 6] {
        core::array::from_fn(|i| {
            self.b
                .iter()
                .zip(dofs)
                .map(|(b, &d)| b[i] * displacement[d])
                .sum()
        })
    }

    fn validate_compatible_strain(
        &self,
        law: &TensorProny,
        dofs: &[usize],
        displacement: &[f64],
        strain: [f64; 6],
    ) -> Result<(), ViscoElementError> {
        if self.compatible_gradients.is_some() {
            law.validate_mandel_strain(&self.compatible_strain(dofs, displacement, strain))?;
        }
        Ok(())
    }

    fn compatible_strain(
        &self,
        dofs: &[usize],
        displacement: &[f64],
        mut strain: [f64; 6],
    ) -> [f64; 6] {
        if let Some(gradients) = &self.compatible_gradients {
            for i in 0..3 {
                strain[i] = gradients
                    .iter()
                    .zip(dofs.chunks_exact(3))
                    .map(|(g, d)| g[i] * displacement[d[i]])
                    .sum();
            }
        }
        strain
    }
}
struct Element {
    dofs: Vec<usize>,
    law: TensorProny,
    density_kg_m3: f64,
    material_state_identity: [u8; 32],
    points: Vec<IntegrationPoint>,
}

/// Prepared geometry/laws and their canonical equilibrium stiffness/mass.
pub struct ViscoElements {
    elements: Vec<Element>,
    equilibrium: TetElasticAssembly,
    full_dofs: usize,
    internal_variables: usize,
    // Free-body modal/history use remains admitted; load control requires the
    // existing geometric rigid-mode constraint admission on every component.
    static_constraint_error: Option<TetElasticError>,
}

/// Storage/loss matrices in the equilibrium assembly's free-DOF order.
#[derive(Debug, Clone)]
pub struct ComplexStiffness {
    /// K′(ω) [N/m].
    pub real: Csr,
    /// K″(ω) [N/m].
    pub imag: Csr,
    /// Angular frequency [rad/s], inside every element law's band.
    pub omega_rad_s: f64,
}

/// Required controls for the bounded dense reference dispersive eigensolve.
#[derive(Debug, Clone, Copy)]
pub struct ViscoModeBudget {
    /// Number of lowest ordered modes to solve.
    pub modes: usize,
    /// Dense matrix dimension cap, at most 256.
    pub max_free_dofs: usize,
    /// Maximum bisection steps per mode.
    pub max_iterations: usize,
    /// Required search interval [rad/s], inside every law's band.
    pub frequency_bracket_rad_s: (f64, f64),
    /// Relative half-width of each returned angular-frequency bracket.
    pub relative_frequency_tolerance: f64,
}

/// A storage-dispersive mode; its loss quotient is not a complex-mode certificate.
#[derive(Debug, Clone)]
pub struct ViscoMode {
    /// Angular frequency [rad/s].
    pub omega_rad_s: f64,
    /// Numerical bisection interval [rad/s], not a roundoff enclosure.
    pub frequency_bracket_rad_s: (f64, f64),
    /// Mass-normalized free-DOF shape at the returned frequency.
    pub phi: Vec<f64>,
    /// φᵀK″(ω)φ / φᵀK′(ω)φ; no diagonal-damping adequacy claim.
    pub loss_factor: f64,
    /// Storage-pencil Ritz value at the returned frequency [s⁻²].
    pub storage_eigenvalue_per_s2: f64,
    /// M⁻¹-norm eigen residual from fs-modal [s⁻²], separate from root error.
    pub eigen_residual_per_s2: f64,
    /// Certified eigenvalue radius [s⁻²]: quadratic (Kato–Temple) where the
    /// complete spectrum certifies the gap, else the residual.
    pub eigen_radius_per_s2: f64,
    /// Relative residual of λ(K′(ω),M) = ω² at the returned frequency.
    pub relative_frequency_residual: f64,
}

impl ViscoElements {
    /// Solve the storage-dispersive eigenproblem at each mode's own frequency.
    ///
    /// For passive Prony laws, K′(ω)/ω² = K_eq/ω² +
    /// Σ_j τ_j² K_j/(1+ω²τ_j²) is decreasing in Loewner order for ω>0.
    /// Hence each positive ordered eigenvalue crosses one at most once. This
    /// brackets the nonlinear root without assuming a contractive fixed-point
    /// iteration. KV contributes loss only. The dense reference eigensolver is
    /// the existing fs-modal owner; its residual must be smaller than the
    /// requested frequency accuracy. Returned brackets are numerical estimates,
    /// not outward-rounded eigenvalue/root certificates. Complex-mode damping
    /// adequacy is a separate fs-modal certification obligation.
    pub fn modes(
        &self,
        budget: ViscoModeBudget,
        cx: &Cx<'_>,
    ) -> Result<Vec<ViscoMode>, ViscoElementError> {
        checkpoint(cx)?;
        let n = self.equilibrium.free_dofs.len();
        let (low, high) = budget.frequency_bracket_rad_s;
        if budget.modes == 0
            || budget.modes > n
            || n > budget.max_free_dofs
            || budget.max_free_dofs > 256
            || budget.max_iterations == 0
            || !(budget.relative_frequency_tolerance.is_finite()
                && budget.relative_frequency_tolerance > 0.0
                && budget.relative_frequency_tolerance < 1.0)
            || !(low.is_finite()
                && low > 0.0
                && high.is_finite()
                && high > low
                && (high * high).is_finite())
        {
            return Err(ViscoElementError::InvalidModalBudget);
        }
        let dense = |a: &Csr| {
            let mut result = vec![0.0; n * n];
            for i in 0..n {
                let (cols, values) = a.row(i);
                for (&j, &v) in cols.iter().zip(values) {
                    result[i * n + j] = v;
                }
            }
            result
        };
        let mass = dense(&self.equilibrium.mass);
        // A dense f64 eigenvector of a stiff mesh keeps an M⁻¹ residual near
        // eps*lambda_max however it is solved, so the residual is no measure
        // of the eigenvalue. Recentre each value on its own Rayleigh quotient
        // and certify it with the quadratic radius of the complete spectrum.
        type Spectrum = (ComplexStiffness, Vec<fs_modal::ModePair>, Vec<f64>);
        let evaluate = |omega| -> Result<Spectrum, ViscoElementError> {
            checkpoint(cx)?;
            let k = self.complex_stiffness(omega, cx)?;
            let pairs = fs_modal::eigh_gen_dense(&dense(&k.real), &mass, n)
                .map_err(ViscoElementError::Modal)?;
            checkpoint(cx)?;
            let (pairs, radii) =
                fs_modal::recentred_dense_radii(&pairs, &k.real, &self.equilibrium.mass)
                    .map_err(ViscoElementError::Modal)?;
            checkpoint(cx)?;
            Ok((k, pairs, radii))
        };
        let (_, at_low, _) = evaluate(low)?;
        let (_, at_high, _) = evaluate(high)?;
        let mut result = Vec::with_capacity(budget.modes);
        for mode in 0..budget.modes {
            if !(at_low[mode].interval.0 > low * low && at_high[mode].interval.1 < high * high) {
                return Err(ViscoElementError::UnbracketedMode { mode });
            }
            let (mut lo, mut hi) = (low, high);
            let mut solved = None;
            for _ in 0..budget.max_iterations {
                let omega = f64::midpoint(lo, hi);
                if omega == lo || omega == hi {
                    break;
                }
                let (k, pairs, radii) = evaluate(omega)?;
                let pair = &pairs[mode];
                let radius = radii[mode];
                let relative = (pair.lambda / (omega * omega) - 1.0).abs();
                if !pair.lambda.is_finite() || !radius.is_finite() || pair.lambda <= 0.0 {
                    return Err(ViscoElementError::ModalConvergence { mode });
                }
                let relative_radius = radius / (omega * omega);
                let maximum = budget.relative_frequency_tolerance / 16.0;
                let frequency_converged = (hi - lo) / (2.0 * omega)
                    <= budget.relative_frequency_tolerance
                    && relative <= 4.0 * budget.relative_frequency_tolerance;
                // A trial far from the root needs a reliable sign, not the
                // final eigen accuracy. Charge that accuracy when returning
                // a mode or when the eigen interval cannot decide the sign.
                let sign_unresolved = (pair.lambda - omega * omega).abs() <= radius;
                if relative_radius > maximum && (frequency_converged || sign_unresolved) {
                    return Err(ViscoElementError::EigenResolution {
                        mode,
                        omega_rad_s: omega,
                        relative_radius,
                        maximum,
                    });
                }
                if frequency_converged {
                    solved = Some(ViscoMode {
                        omega_rad_s: omega,
                        frequency_bracket_rad_s: (lo, hi),
                        phi: pair.phi.clone(),
                        loss_factor: modal_loss_energy_quotient(&k, &pair.phi)?,
                        storage_eigenvalue_per_s2: pair.lambda,
                        eigen_residual_per_s2: pair.residual,
                        eigen_radius_per_s2: radius,
                        relative_frequency_residual: relative,
                    });
                    break;
                }
                if pair.lambda > omega * omega {
                    lo = omega;
                } else {
                    hi = omega;
                }
            }
            result.push(solved.ok_or(ViscoElementError::ModalConvergence { mode })?);
        }
        Ok(result)
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ViscoElementError> {
    cx.checkpoint()
        .map_err(|_| TetElasticError::Cancelled.into())
}

// The two existing hex geometry kernels differ only in their node count.
// Keep the Mandel pullback, quadrature order and source copying identical.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn prepare_hex<const N: usize, const D: usize>(
    prepared: &mut Vec<Element>,
    nodes: &[[f64; 3]],
    cells: &[[usize; N]],
    field: TetViscoField<'_>,
    rule: &[([f64; 3], f64)],
    minimum_scaled_jacobian: f64,
    point: fn([[f64; 3]; N], [f64; 3], f64, usize, usize) -> Result<HexPoint<N>, TetElasticError>,
    cx: &Cx<'_>,
) -> Result<(), ViscoElementError> {
    for (element, conn) in cells.iter().enumerate() {
        checkpoint(cx)?;
        let coordinates = conn.map(|node| nodes[node]);
        let mut points = Vec::with_capacity(rule.len());
        for (index, &(location, weight)) in rule.iter().enumerate() {
            checkpoint(cx)?;
            let geometry = point(
                coordinates,
                location,
                minimum_scaled_jacobian,
                element,
                index,
            )?;
            let b = super::isoparametric::mandel_b::<N, D>(geometry.gradients);
            points.push(IntegrationPoint {
                b: (0..D)
                    .map(|column| core::array::from_fn(|row| b[row][column]))
                    .collect(),
                compatible_gradients: None,
                volume_m3: weight * geometry.determinant_m3,
            });
        }
        prepared.push(Element {
            dofs: conn
                .iter()
                .flat_map(|node| [3 * node, 3 * node + 1, 3 * node + 2])
                .collect(),
            law: field.at(element).law.clone(),
            density_kg_m3: field.at(element).density_kg_m3,
            material_state_identity: field.at(element).material_state_identity,
            points,
        });
    }
    Ok(())
}

impl ViscoElements {
    /// Prepare P1/TET10/HEX8/HEX20/HEX27 using existing elastic admission/assembly.
    /// Caps are checked before allocating per-element laws or histories.
    pub fn new(
        mesh: SolidMeshRef<'_>,
        field: TetViscoField<'_>,
        fixed: &[usize],
        budget: TetAssemblyBudget,
        state_budget: ViscoStateBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, ViscoElementError> {
        checkpoint(cx)?;
        let (nodes, count, points) = match mesh {
            SolidMeshRef::Tet4 { nodes_m, elements } => (nodes_m, elements.len(), 1),
            SolidMeshRef::Tet10 {
                nodes_m,
                elements,
                quadrature,
            } => {
                quadrature.validate()?;
                (
                    nodes_m,
                    elements.len(),
                    quadrature.stiffness_points_per_dir.pow(3),
                )
            }
            SolidMeshRef::Hex8 {
                nodes_m,
                elements,
                stiffness_points_per_dir,
                ..
            } => {
                super::hex8::validate_stiffness_order(stiffness_points_per_dir)?;
                (nodes_m, elements.len(), stiffness_points_per_dir.pow(3))
            }
            SolidMeshRef::Hex {
                nodes_m,
                elements,
                quadrature,
            } => {
                quadrature.validate()?;
                let count = match elements {
                    HexConnectivity::Hex20(cells) => cells.len(),
                    HexConnectivity::Hex27(cells) => cells.len(),
                };
                (nodes_m, count, quadrature.points_per_dir.pow(3))
            }
        };
        for (what, requested, maximum) in [
            ("nodes", nodes.len(), budget.maximum_nodes),
            ("tetrahedra", count, budget.maximum_tetrahedra),
        ] {
            if requested > maximum {
                return Err(TetElasticError::BudgetExceeded {
                    what,
                    requested,
                    maximum,
                }
                .into());
            }
        }
        if let TetViscoField::PerElement(m) = field {
            if m.len() != count {
                return Err(TetElasticError::MaterialCountMismatch {
                    expected: count,
                    got: m.len(),
                }
                .into());
            }
        }
        let mut needed = 0usize;
        for e in 0..count {
            checkpoint(cx)?;
            needed = needed
                .checked_add(
                    field
                        .at(e)
                        .law
                        .branches()
                        .len()
                        .checked_mul(points)
                        .and_then(|v| v.checked_mul(6))
                        .ok_or(TetElasticError::SizeOverflow)?,
                )
                .ok_or(TetElasticError::SizeOverflow)?;
        }
        if needed > state_budget.max_internal_variables {
            return Err(ViscoElementError::InternalStateBudget {
                needed,
                cap: state_budget.max_internal_variables,
            });
        }
        let mut materials = Vec::with_capacity(count);
        for e in 0..count {
            checkpoint(cx)?;
            let m = field.at(e);
            materials.push(m.equilibrium_material()?);
        }
        let eq = match mesh {
            SolidMeshRef::Tet4 { nodes_m, elements } => TetLinearElasticProblem {
                nodes_m,
                tetrahedra: elements,
                materials: TetMaterialField::PerElement(&materials),
                fixed_dofs: fixed,
                budget,
            }
            .assemble(cx)?,
            SolidMeshRef::Tet10 {
                nodes_m,
                elements,
                quadrature,
            } => QuadraticTetProblem {
                nodes_m,
                elements,
                materials: TetMaterialField::PerElement(&materials),
                fixed_dofs: fixed,
                budget,
                quadrature,
            }
            .assemble(cx)?,
            SolidMeshRef::Hex8 {
                nodes_m,
                elements,
                formulation,
                stiffness_points_per_dir,
            } => Hex8ElasticProblem {
                nodes_m,
                hexahedra: elements,
                materials: TetMaterialField::PerElement(&materials),
                formulation,
                fixed_dofs: fixed,
                budget,
                stiffness_points_per_dir,
            }
            .assemble(cx)?,
            SolidMeshRef::Hex {
                nodes_m,
                elements,
                quadrature,
            } => HexElasticProblem {
                nodes_m,
                family: match elements {
                    HexConnectivity::Hex20(_) => HexFamily::Serendipity20,
                    HexConnectivity::Hex27(_) => HexFamily::Lagrange27,
                },
                elements,
                materials: TetMaterialField::PerElement(&materials),
                fixed_dofs: fixed,
                budget,
                quadrature,
            }
            .assemble(cx)?,
        };
        let mut prepared = Vec::with_capacity(count);
        match mesh {
            SolidMeshRef::Tet4 { elements, .. } => {
                for (e, conn) in elements.iter().enumerate() {
                    checkpoint(cx)?;
                    let g = TetGeometry::try_new(
                        e,
                        conn.map(|n| nodes[n]),
                        budget.minimum_scaled_jacobian,
                    )?;
                    let b = super::isoparametric::mandel_b::<4, 12>(g.gradients);
                    prepared.push(Element {
                        dofs: conn
                            .iter()
                            .flat_map(|n| [3 * n, 3 * n + 1, 3 * n + 2])
                            .collect(),
                        law: field.at(e).law.clone(),
                        density_kg_m3: field.at(e).density_kg_m3,
                        material_state_identity: field.at(e).material_state_identity,
                        points: vec![IntegrationPoint {
                            b: (0..12).map(|j| core::array::from_fn(|i| b[i][j])).collect(),
                            compatible_gradients: None,
                            volume_m3: g.volume,
                        }],
                    });
                }
            }
            SolidMeshRef::Tet10 {
                elements,
                quadrature,
                ..
            } => {
                let rule = duffy_quadrature(quadrature.stiffness_points_per_dir);
                for (e, conn) in elements.iter().enumerate() {
                    checkpoint(cx)?;
                    let mut p = Vec::with_capacity(rule.len());
                    for (q, &(lambda, weight)) in rule.iter().enumerate() {
                        checkpoint(cx)?;
                        let g = tet10_point(
                            conn.map(|n| nodes[n]),
                            lambda,
                            budget.minimum_scaled_jacobian,
                            e,
                            q,
                        )?;
                        let b = super::isoparametric::mandel_b::<10, 30>(g.gradients);
                        p.push(IntegrationPoint {
                            b: (0..30).map(|j| core::array::from_fn(|i| b[i][j])).collect(),
                            compatible_gradients: None,
                            volume_m3: weight * g.determinant_m3,
                        });
                    }
                    prepared.push(Element {
                        dofs: conn
                            .iter()
                            .flat_map(|n| [3 * n, 3 * n + 1, 3 * n + 2])
                            .collect(),
                        law: field.at(e).law.clone(),
                        density_kg_m3: field.at(e).density_kg_m3,
                        material_state_identity: field.at(e).material_state_identity,
                        points: p,
                    });
                }
            }
            SolidMeshRef::Hex8 {
                nodes_m,
                elements,
                formulation,
                stiffness_points_per_dir,
            } => hex8::prepare(
                &mut prepared,
                &Hex8ElasticProblem {
                    nodes_m,
                    hexahedra: elements,
                    materials: TetMaterialField::PerElement(&materials),
                    formulation,
                    fixed_dofs: fixed,
                    budget,
                    stiffness_points_per_dir,
                },
                field,
                cx,
            )?,
            SolidMeshRef::Hex {
                elements,
                quadrature,
                ..
            } => {
                let rule = quadrature.points()?;
                match elements {
                    HexConnectivity::Hex20(cells) => prepare_hex::<20, 60>(
                        &mut prepared,
                        nodes,
                        cells,
                        field,
                        &rule,
                        budget.minimum_scaled_jacobian,
                        hex20_point,
                        cx,
                    )?,
                    HexConnectivity::Hex27(cells) => prepare_hex::<27, 81>(
                        &mut prepared,
                        nodes,
                        cells,
                        field,
                        &rule,
                        budget.minimum_scaled_jacobian,
                        hex27_point,
                        cx,
                    )?,
                }
            }
        }
        let static_constraint_error = match mesh {
            SolidMeshRef::Tet4 { elements, .. } => {
                super::validate_rigid_mode_constraints(nodes, elements, fixed)
            }
            SolidMeshRef::Tet10 { elements, .. } => {
                super::validate_rigid_mode_constraints(nodes, elements, fixed)
            }
            SolidMeshRef::Hex8 { elements, .. } => {
                super::validate_rigid_mode_constraints(nodes, elements, fixed)
            }
            SolidMeshRef::Hex {
                elements: HexConnectivity::Hex20(elements),
                ..
            } => super::validate_rigid_mode_constraints(nodes, elements, fixed),
            SolidMeshRef::Hex {
                elements: HexConnectivity::Hex27(elements),
                ..
            } => super::validate_rigid_mode_constraints(nodes, elements, fixed),
        }
        .err();
        checkpoint(cx)?;
        Ok(Self {
            elements: prepared,
            equilibrium: eq,
            full_dofs: nodes.len() * 3,
            internal_variables: needed,
            static_constraint_error,
        })
    }

    /// Canonical equilibrium elastic pencil, including consistent mass.
    pub fn equilibrium(&self) -> &TetElasticAssembly {
        &self.equilibrium
    }

    /// Original material, reference density and source identity for an element.
    /// This is the material that assembled the retained reference mass matrix.
    pub fn reference_material(&self, element: usize) -> Option<TetViscoMaterial> {
        self.elements.get(element).map(|e| TetViscoMaterial {
            law: e.law.clone(),
            density_kg_m3: e.density_kg_m3,
            material_state_identity: e.material_state_identity,
        })
    }

    /// Count charged to the caller's state budget.
    pub const fn internal_variables(&self) -> usize {
        self.internal_variables
    }

    /// Assemble ∫Bᵀ C*(ω) B dV without diagonalizing or dropping cross terms.
    pub fn complex_stiffness(
        &self,
        omega_rad_s: f64,
        cx: &Cx<'_>,
    ) -> Result<ComplexStiffness, ViscoElementError> {
        checkpoint(cx)?;
        // Check all bands before starting the assembly; zero-frequency queries
        // are admitted only if zero belongs to every constitutive band.
        let tensors = self
            .elements
            .iter()
            .map(|e| {
                checkpoint(cx)?;
                Ok(e.law.complex_modulus(omega_rad_s)?)
            })
            .collect::<Result<Vec<_>, ViscoElementError>>()?;
        let n = self.equilibrium.free_dofs.len();
        let mut reduced = vec![None; self.full_dofs];
        for (i, &dof) in self.equilibrium.free_dofs.iter().enumerate() {
            reduced[dof] = Some(i);
        }
        let mut real = Coo::new(n, n);
        let mut imag = Coo::new(n, n);
        let mut equilibrium_storage = true;
        for (element, (storage, loss)) in self.elements.iter().zip(tensors) {
            checkpoint(cx)?;
            equilibrium_storage &= storage == *element.law.c_eq();
            let d = element.dofs.len();
            let mut kr = vec![0.0; d * d];
            let mut ki = vec![0.0; d * d];
            for p in &element.points {
                checkpoint(cx)?;
                for i in 0..d {
                    for j in i..d {
                        kr[i * d + j] += p.volume_m3 * bilinear(&p.b[i], &storage, &p.b[j]);
                        ki[i * d + j] += p.volume_m3 * bilinear(&p.b[i], &loss, &p.b[j]);
                    }
                }
            }
            for i in 0..d {
                for j in i..d {
                    let (r, l) = (kr[i * d + j], ki[i * d + j]);
                    if !r.is_finite() || !l.is_finite() {
                        return Err(ViscoElementError::NonFiniteResult);
                    }
                    if let (Some(row), Some(col)) =
                        (reduced[element.dofs[i]], reduced[element.dofs[j]])
                    {
                        real.push(row, col, r);
                        imag.push(row, col, l);
                        if row != col {
                            real.push(col, row, r);
                            imag.push(col, row, l);
                        }
                    }
                }
            }
        }
        let real = if equilibrium_storage {
            self.equilibrium.stiffness.clone()
        } else {
            real.assemble()
        };
        let imag = imag.assemble();
        checkpoint(cx)?;
        for a in [&real, &imag] {
            if (0..n).any(|i| a.row(i).1.iter().any(|v| !v.is_finite())) {
                return Err(ViscoElementError::NonFiniteResult);
            }
        }
        Ok(ComplexStiffness {
            real,
            imag,
            omega_rad_s,
        })
    }

    /// Check the complete displacement against every original material point.
    /// This reuses the element's exact Mandel strain map and source domain;
    /// it evolves no history and introduces no alternate stiffness or stress.
    /// High-order elements check declared quadrature points, not unsampled extrema.
    /// B-bar checks both assumed and raw compatible strain at those points.
    pub fn validate_displacement(
        &self,
        displacement_m: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ViscoElementError> {
        if cancelled() {
            return Err(TetElasticError::Cancelled.into());
        }
        if displacement_m.len() != self.full_dofs
            || displacement_m.iter().any(|value| !value.is_finite())
        {
            return Err(ViscoElementError::InvalidHistory {
                what: "complete finite node-major displacement required",
            });
        }
        for element in &self.elements {
            for point in &element.points {
                if cancelled() {
                    return Err(TetElasticError::Cancelled.into());
                }
                let strain = point.strain(&element.dofs, displacement_m);
                element.law.validate_mandel_strain(&strain)?;
                point.validate_compatible_strain(
                    &element.law,
                    &element.dofs,
                    displacement_m,
                    strain,
                )?;
            }
        }
        if cancelled() {
            return Err(TetElasticError::Cancelled.into());
        }
        Ok(())
    }

    /// Start a relaxed zero-displacement history, bound to this exact system.
    pub fn state(&self, cx: &Cx<'_>) -> Result<ViscoElementState<'_>, ViscoElementError> {
        let mut histories = Vec::with_capacity(self.elements.len());
        for e in &self.elements {
            let mut points = Vec::with_capacity(e.points.len());
            for _ in &e.points {
                checkpoint(cx)?;
                points.push(e.law.state());
            }
            histories.push(points);
        }
        Ok(ViscoElementState {
            system: self,
            histories,
            energy: ViscoEnergy::default(),
            displacement_m: vec![0.0; self.full_dofs],
            thermal_strains: None,
            eigenstrain_work_j: 0.0,
        })
    }
}

fn bilinear(x: &[f64; 6], c: &Mandel6, y: &[f64; 6]) -> f64 {
    (0..6)
        .map(|i| x[i] * (0..6).map(|j| c[i][j] * y[j]).sum::<f64>())
        .sum()
}

/// Integrated constitutive ledger [J], excluding inertia/support/radiation.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ViscoEnergy {
    /// Current Helmholtz storage.
    pub stored_j: f64,
    /// Cumulative exact mechanical-strain work: nodal work plus imposed
    /// eigenstrain work. It equals storage change plus constitutive heat.
    pub work_j: f64,
    /// Cumulative dashpot heat, integrated as sums of squares, without clamping.
    pub dissipated_j: f64,
}

/// Bound quadrature histories. Only successful steps publish state changes.
#[derive(Clone)]
pub struct ViscoElementState<'a> {
    system: &'a ViscoElements,
    histories: Vec<Vec<TensorPronyState>>,
    energy: ViscoEnergy,
    displacement_m: Vec<f64>,
    thermal_strains: Option<std::sync::Arc<thermal::ThermalSnapshot>>,
    eigenstrain_work_j: f64,
}

/// Work-conjugate force observation for one admitted prescribed-history step.
pub struct ViscoForceStep {
    /// Endpoint internal force at every DOF [N]. Inertia and applied-load
    /// subtraction belong to the consuming mechanics owner.
    pub full_force_n: Vec<f64>,
    /// Endpoint internal force in canonical free-DOF order [N].
    pub reduced_force_n: Vec<f64>,
    /// Exact interval-mean internal force, work-conjugate to the displacement
    /// increment; this is generally not the average of endpoint forces [N].
    pub full_mean_force_n: Vec<f64>,
    /// Exact interval-mean internal force in canonical free-DOF order [N].
    pub reduced_mean_force_n: Vec<f64>,
    /// Total integrated energy after this step [J].
    pub energy: ViscoEnergy,
    /// Constitutive work during this step, not endpoint force times displacement.
    pub step_work_j: f64,
    /// Actual interval-mean nodal internal force times displacement increment,
    /// including prescribed supports [J].
    pub step_mechanical_work_j: f64,
    /// Signed imposed-expansion work `-integral sigma_mean : delta epsilon_T`
    /// [J]. Together with mechanical work it supplies `step_work_j`. This is
    /// an external eigenstrain work port, not caloric heat or added damping.
    pub step_eigenstrain_work_j: f64,
    /// Nonnegative integrated heat during this step, without clamping.
    pub step_dissipated_j: f64,
}

impl ViscoElementState<'_> {
    /// Accepted full node-major displacement, including prescribed supports [m].
    pub fn displacement_m(&self) -> &[f64] {
        &self.displacement_m
    }

    /// Read the current constitutive energy ledger.
    pub const fn energy(&self) -> ViscoEnergy {
        self.energy
    }

    /// Prescribe complete node-major xyz displacements [m] after `dt_s > 0`.
    /// Fixed DOFs are prescribed inputs too; their force rows retain reactions.
    /// A finite history is not a certification that its spectrum lies inside
    /// the law's band. Frequency-domain queries enforce that band explicitly.
    /// Any previously accepted thermal strain is held, including its original
    /// reference temperature and source identity.
    pub fn step_internal_force(
        &mut self,
        u_new: &[f64],
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<ViscoForceStep, ViscoElementError> {
        self.step_internal_force_impl(u_new, dt_s, None, cx)
    }

    fn step_internal_force_impl(
        &mut self,
        u_new: &[f64],
        dt_s: f64,
        thermal: Option<&std::sync::Arc<thermal::ThermalSnapshot>>,
        cx: &Cx<'_>,
    ) -> Result<ViscoForceStep, ViscoElementError> {
        checkpoint(cx)?;
        if u_new.len() != self.system.full_dofs || u_new.iter().any(|v| !v.is_finite()) {
            return Err(ViscoElementError::InvalidHistory {
                what: "complete finite node-major displacement required",
            });
        }
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(ViscoElementError::InvalidHistory {
                what: "finite positive dt_s required",
            });
        }
        let mut histories = Vec::with_capacity(self.histories.len());
        for states in &self.histories {
            let mut candidate = Vec::with_capacity(states.len());
            for state in states {
                checkpoint(cx)?;
                candidate.push(state.clone());
            }
            histories.push(candidate);
        }
        let mut force = vec![0.0; self.system.full_dofs];
        let mut mean_force = vec![0.0; self.system.full_dofs];
        let mut energy = ViscoEnergy::default();
        let mut step_work_j = 0.0;
        let mut step_dissipated_j = 0.0;
        let mut step_eigenstrain_work_j = 0.0;
        let target = thermal
            .map(|s| s.as_ref())
            .or(self.thermal_strains.as_deref());
        for (element, (e, states)) in self.system.elements.iter().zip(&mut histories).enumerate() {
            for (point, (p, state)) in e.points.iter().zip(states).enumerate() {
                checkpoint(cx)?;
                let delta_free =
                    self.thermal_increment(element, point, thermal.map(|s| s.as_ref()));
                let strain = p.strain(&e.dofs, u_new);
                let strain = if let Some(thermal) = target {
                    p.thermal_mechanical_strain(
                        e,
                        u_new,
                        strain,
                        *thermal.at(element, point),
                        thermal.free_strain(element, point),
                    )?
                } else {
                    p.validate_compatible_strain(&e.law, &e.dofs, u_new, strain)?;
                    strain
                };
                let response = e.law.step_with_mean(state, strain, dt_s)?;
                let stress = response.stress;
                if stress
                    .iter()
                    .chain(state.q.iter().flatten())
                    .any(|v| !v.is_finite())
                    || !state.stored.is_finite()
                    || !state.work.is_finite()
                    || !state.dissipated.is_finite()
                {
                    return Err(ViscoElementError::NonFiniteResult);
                }
                for (b, &d) in p.b.iter().zip(&e.dofs) {
                    force[d] += p.volume_m3 * b.iter().zip(stress).map(|(a, s)| a * s).sum::<f64>();
                    mean_force[d] += p.volume_m3
                        * b.iter()
                            .zip(response.mean_stress)
                            .map(|(a, s)| a * s)
                            .sum::<f64>();
                }
                energy.stored_j += p.volume_m3 * state.stored;
                energy.work_j += p.volume_m3 * state.work;
                energy.dissipated_j += p.volume_m3 * state.dissipated;
                step_work_j += p.volume_m3 * response.step_work;
                step_dissipated_j += p.volume_m3 * response.step_dissipated;
                step_eigenstrain_work_j -= p.volume_m3
                    * response
                        .mean_stress
                        .iter()
                        .zip(delta_free)
                        .map(|(s, e)| s * e)
                        .sum::<f64>();
            }
        }
        let step_mechanical_work_j = mean_force
            .iter()
            .zip(u_new.iter().zip(&self.displacement_m))
            .map(|(f, (new, old))| f * (new - old))
            .sum::<f64>();
        let eigenstrain_work_j = self.eigenstrain_work_j + step_eigenstrain_work_j;
        if force.iter().chain(&mean_force).any(|v| !v.is_finite())
            || [
                energy.stored_j,
                energy.work_j,
                energy.dissipated_j,
                step_work_j,
                step_dissipated_j,
                step_mechanical_work_j,
                step_eigenstrain_work_j,
                eigenstrain_work_j,
            ]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(ViscoElementError::NonFiniteResult);
        }
        checkpoint(cx)?;
        let reduced_force_n = self
            .system
            .equilibrium
            .free_dofs
            .iter()
            .map(|&d| force[d])
            .collect();
        let reduced_mean_force_n = self
            .system
            .equilibrium
            .free_dofs
            .iter()
            .map(|&d| mean_force[d])
            .collect();
        self.histories = histories;
        self.energy = energy;
        self.displacement_m.copy_from_slice(u_new);
        if let Some(strains) = thermal {
            self.thermal_strains = Some(std::sync::Arc::clone(strains));
        }
        self.eigenstrain_work_j = eigenstrain_work_j;
        Ok(ViscoForceStep {
            full_force_n: force,
            reduced_force_n,
            full_mean_force_n: mean_force,
            reduced_mean_force_n,
            energy,
            step_work_j,
            step_mechanical_work_j,
            step_eigenstrain_work_j,
            step_dissipated_j,
        })
    }
}

/// φᵀK″φ / φᵀK′φ, retaining all modal cross terms. Invalid/zero energy refuses.
pub fn modal_loss_energy_quotient(
    k: &ComplexStiffness,
    phi: &[f64],
) -> Result<f64, ViscoElementError> {
    let n = k.real.nrows();
    if k.real.ncols() != n
        || k.imag.nrows() != n
        || k.imag.ncols() != n
        || phi.len() != n
        || phi.iter().any(|v| !v.is_finite())
    {
        return Err(ViscoElementError::InvalidHistory {
            what: "finite free-DOF modal vector required",
        });
    }
    let mut real = vec![0.0; n];
    let mut imag = vec![0.0; n];
    k.real.spmv(phi, &mut real);
    k.imag.spmv(phi, &mut imag);
    let storage: f64 = phi.iter().zip(real).map(|(a, b)| a * b).sum();
    let loss: f64 = phi.iter().zip(imag).map(|(a, b)| a * b).sum();
    if !storage.is_finite()
        || storage <= 0.0
        || !loss.is_finite()
        || loss < 0.0
        || !(loss / storage).is_finite()
    {
        return Err(ViscoElementError::InvalidHistory {
            what: "positive finite storage and nonnegative loss required",
        });
    }
    Ok(loss / storage)
}
