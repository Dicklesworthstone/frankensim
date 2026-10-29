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
//! histories are inputs; this module is not a displacement time integrator.

use super::tet10::{QuadraticTetProblem, TetQuadrature, tet10_point};
use super::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticError, TetElasticMaterial, TetGeometry,
    TetLinearElasticProblem, TetMaterialField,
};
use fs_exec::Cx;
use fs_feec::highorder::simplex::duffy_quadrature;
use fs_material::visco::ViscoError;
use fs_material::visco::tensor::{Mandel6, TensorProny, TensorPronyState};
use fs_sparse::{Coo, Csr};

mod project;
pub use project::{ProjectedViscoSolid, SolidProjectionBudget};

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
    volume_m3: f64,
}
struct Element {
    dofs: Vec<usize>,
    law: TensorProny,
    points: Vec<IntegrationPoint>,
}

/// Prepared geometry/laws and their canonical equilibrium stiffness/mass.
pub struct ViscoElements {
    elements: Vec<Element>,
    equilibrium: TetElasticAssembly,
    full_dofs: usize,
    internal_variables: usize,
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

impl ViscoElements {
    /// Prepare P1/TET10 geometry using existing elastic admission/assembly.
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
            materials.push(TetElasticMaterial::try_new_mandel(
                m.density_kg_m3,
                *m.law.c_eq(),
                m.material_state_identity,
            )?);
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
                        points: vec![IntegrationPoint {
                            b: (0..12).map(|j| core::array::from_fn(|i| b[i][j])).collect(),
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
                            volume_m3: weight * g.determinant_m3,
                        });
                    }
                    prepared.push(Element {
                        dofs: conn
                            .iter()
                            .flat_map(|n| [3 * n, 3 * n + 1, 3 * n + 2])
                            .collect(),
                        law: field.at(e).law.clone(),
                        points: p,
                    });
                }
            }
        }
        Ok(Self {
            elements: prepared,
            equilibrium: eq,
            full_dofs: nodes.len() * 3,
            internal_variables: needed,
        })
    }

    /// Canonical equilibrium elastic pencil, including consistent mass.
    pub fn equilibrium(&self) -> &TetElasticAssembly {
        &self.equilibrium
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
    /// TET10 checks its declared quadrature points, not unsampled extrema.
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
                let strain = core::array::from_fn(|i| {
                    point.b.iter().zip(&element.dofs)
                        .map(|(b, &dof)| b[i] * displacement_m[dof]).sum()
                });
                element.law.validate_mandel_strain(&strain)?;
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
    /// Cumulative exact strain-ramp work.
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
}

/// Work-conjugate force observation for one admitted prescribed-history step.
pub struct ViscoForceStep {
    /// Endpoint internal force at every DOF, including support reactions [N].
    pub full_force_n: Vec<f64>,
    /// Endpoint internal force in canonical free-DOF order [N].
    pub reduced_force_n: Vec<f64>,
    /// Total integrated energy after this step [J].
    pub energy: ViscoEnergy,
    /// Constitutive work during this step, not endpoint force times displacement.
    pub step_work_j: f64,
    /// Nonnegative integrated heat during this step, without clamping.
    pub step_dissipated_j: f64,
}

impl ViscoElementState<'_> {
    /// Read the current constitutive energy ledger.
    pub const fn energy(&self) -> ViscoEnergy {
        self.energy
    }

    /// Prescribe complete node-major xyz displacements [m] after `dt_s > 0`.
    /// Fixed DOFs are prescribed inputs too; their force rows retain reactions.
    /// A finite history is not a certification that its spectrum lies inside
    /// the law's band. Frequency-domain queries enforce that band explicitly.
    pub fn step_internal_force(
        &mut self,
        u_new: &[f64],
        dt_s: f64,
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
        let mut energy = ViscoEnergy::default();
        let mut step_work_j = 0.0;
        let mut step_dissipated_j = 0.0;
        for (e, states) in self.system.elements.iter().zip(&mut histories) {
            for (p, state) in e.points.iter().zip(states) {
                checkpoint(cx)?;
                let strain = core::array::from_fn(|i| {
                    p.b.iter().zip(&e.dofs).map(|(b, &d)| b[i] * u_new[d]).sum()
                });
                let old_work = state.work;
                let old_heat = state.dissipated;
                let stress = e.law.step(state, strain, dt_s)?;
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
                }
                energy.stored_j += p.volume_m3 * state.stored;
                energy.work_j += p.volume_m3 * state.work;
                energy.dissipated_j += p.volume_m3 * state.dissipated;
                step_work_j += p.volume_m3 * (state.work - old_work);
                step_dissipated_j += p.volume_m3 * (state.dissipated - old_heat);
            }
        }
        if force.iter().any(|v| !v.is_finite())
            || [
                energy.stored_j,
                energy.work_j,
                energy.dissipated_j,
                step_work_j,
                step_dissipated_j,
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
        self.histories = histories;
        self.energy = energy;
        Ok(ViscoForceStep {
            full_force_n: force,
            reduced_force_n,
            energy,
            step_work_j,
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
