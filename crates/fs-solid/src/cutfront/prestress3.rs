//! Signed initial-stress modes on the original Cartesian or adaptive CutFEM grid.
//!
//! The native owners apply (K + G) phi = lambda M phi. Only the consistent
//! physical mass is assembled and factored sparsely; positive LDLᵀ whitening
//! gives the symmetric operator Sᵀ(K + G)S, with Sᵀ M S = I. Neither K nor G
//! is probed on coordinate vectors or assembled into a global dense matrix.
//! The existing spectral service targets the lowest algebraic Ritz values.
//! Negative and zero values are retained, never clamped or square-rooted.
//!
//! The result contains residual-admitted numerical estimates, not a certificate
//! that every lower mode was found. No inertia/count theorem, exhaustive
//! buckling factor, nonlinear stability, finite-deformation tangent or
//! follower-load correction is inferred. The prepared stress must represent
//! the caller's intended reference equilibrium; its owner admits the material
//! state but does not prove equilibrium. Geometry, supports, mass fractions,
//! elastic scales and geometric prestress stay fixed throughout this query.

use std::cell::{Cell, RefCell};
use std::ops::ControlFlow;

use fs_cutfem::elastic3::mass::MaterialMass3;
use fs_cutfem::elastic3::{ElasticityError3, PreparedPrestress3};
use fs_exec::CancelGate;
use fs_sparse::direct::DirectLimits;
use fs_sparse::direct::positive::{PositiveFactorError, PositiveWhitening};
use fs_sparse::{Coo, DirectOrdering, LdltError, LdltFactor, LdltOptions, SymbolicLdlt};
use fs_spectral::service::{
    CertifiedEigenvalue, EigenBackend, EigenQuery, EigenService, ServiceError, SymmetricOp,
};

use super::dynamics3::CutElasticityRef3;

/// Explicit native, sparse-factor and spectral-work allowances.
#[derive(Debug, Clone, Copy)]
pub struct CutPrestressLimits3 {
    /// Bounds independent DOFs and the adaptive physical reconstruction size.
    pub max_native_dofs: usize,
    /// Bounds the original, indivisible native stiffness action.
    pub max_cells: usize,
    /// Conservative sixfold bound on native scalar mass entries, before
    /// duplicate summation, xyz replication and off-diagonal mirroring.
    pub max_mass_entries: usize,
    /// Maximum retained scalar geometric-stiffness entries.
    pub max_geometric_entries: usize,
    /// Native sparse symbolic, frontal and scalar-work admission.
    pub factor: DirectLimits,
    /// Numeric sparse factor plus its front/scratch payload; excludes sources.
    pub max_factor_bytes: usize,
    /// Conservative scalar/index workspace for native actions, LOBPCG blocks,
    /// rollback snapshots and returned modes. Sparse mass assembly and factors
    /// have their separate limits above; allocator overhead is not process RSS.
    pub max_workspace_components: usize,
    /// Total native K+G actions per query, including original-pencil checks.
    pub max_operator_applications: usize,
}

impl Default for CutPrestressLimits3 {
    fn default() -> Self {
        Self {
            max_native_dofs: 100_000,
            max_cells: 32_768,
            max_mass_entries: 4_000_000,
            max_geometric_entries: 2_000_000,
            factor: DirectLimits {
                max_symbolic_entries: 4_000_000,
                max_analysis_work: 400_000_000,
                max_numeric_work: 4_000_000_000,
                max_front_dim: 512,
            },
            max_factor_bytes: 256 * 1024 * 1024,
            max_workspace_components: 8_000_000,
            max_operator_applications: 100_000,
        }
    }
}

/// Fixed query and numerical admission policy for signed Ritz estimates.
#[derive(Debug, Clone, Copy)]
pub struct CutPrestressModeOptions3 {
    /// Block size; the native LOBPCG owner requires 3*modes <= free DOFs.
    pub modes: usize,
    /// Maximum one-iteration service ticks in this query.
    pub max_ticks: usize,
    /// Absolute residual tolerance of the mass-whitened spectral operator.
    pub spectral_tolerance: f64,
    /// Original-pencil whitened residual must satisfy absolute + relative*scale.
    pub residual_relative_tolerance: f64,
    /// Absolute allowance is necessary near zero eigenvalues. This has the
    /// same units as the whitened eigenvalues (s^-2 with SI native sources).
    pub residual_absolute_tolerance: f64,
    /// Maximum entrywise defect of the returned modes' physical M Gram matrix.
    pub mass_orthogonality_tolerance: f64,
}

impl Default for CutPrestressModeOptions3 {
    fn default() -> Self {
        Self {
            modes: 1,
            max_ticks: 256,
            spectral_tolerance: 1e-10,
            residual_relative_tolerance: 1e-8,
            residual_absolute_tolerance: 1e-10,
            mass_orthogonality_tolerance: 1e-8,
        }
    }
}

/// Typed preparation, numerical convergence and cancellation refusals.
#[derive(Debug)]
pub enum CutPrestressError3 {
    /// Invalid query policy, inertia or unrepresentable numerical arithmetic.
    Invalid(&'static str),
    /// An explicit native storage, workspace or application cap was exceeded.
    Budget(&'static str),
    /// Original material, geometry or prepared-state owner refused.
    Physics(ElasticityError3),
    /// Sparse mass analysis or factorization refused without regularization.
    Factor(LdltError),
    /// The complete mass factor or one triangular action was not admissible.
    Whitening(PositiveFactorError),
    /// The shared spectral owner refused its query or one bounded tick.
    Spectral(ServiceError),
    /// No partial mode packet is published on any failed numerical gate.
    Unresolved {
        /// Service iterations consumed.
        ticks: usize,
        /// Largest available whitened residual, or infinity before a packet exists.
        whitened_residual: f64,
        /// Largest checked physical mass-Gram defect, or infinity if unchecked.
        mass_orthogonality_error: f64,
    },
    /// Cooperative cancellation observed before result publication.
    Cancelled,
}

impl std::fmt::Display for CutPrestressError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "CutFEM prestress modes refused: {self:?}")
    }
}
impl std::error::Error for CutPrestressError3 {}
impl From<ElasticityError3> for CutPrestressError3 {
    fn from(error: ElasticityError3) -> Self {
        match error {
            ElasticityError3::Cancelled => Self::Cancelled,
            other => Self::Physics(other),
        }
    }
}
impl From<LdltError> for CutPrestressError3 {
    fn from(error: LdltError) -> Self {
        match error {
            LdltError::Cancelled => Self::Cancelled,
            other => Self::Factor(other),
        }
    }
}
impl From<PositiveFactorError> for CutPrestressError3 {
    fn from(error: PositiveFactorError) -> Self {
        match error {
            PositiveFactorError::Cancelled => Self::Cancelled,
            other => Self::Whitening(other),
        }
    }
}
impl From<ServiceError> for CutPrestressError3 {
    fn from(error: ServiceError) -> Self {
        match error {
            ServiceError::Cancelled => Self::Cancelled,
            other => Self::Spectral(other),
        }
    }
}

/// One signed eigenvalue of the frozen initial-stress approximation.
#[derive(Debug, Clone)]
pub struct CutPrestressMode3 {
    /// lambda, in s^-2 with SI sources. A nonpositive value is not a real
    /// vibration frequency. Values are not rounded to zero or absolutized.
    pub eigenvalue: f64,
    /// Native independent coordinates, with exact zero strong-support entries.
    /// Adaptive hanging reconstruction remains with the original native owner.
    /// Physical normalization is phi^T M phi = 1 to the reported tolerance.
    pub vector: Vec<f64>,
    /// ||Sᵀ[(K+G)phi - lambda M phi]||, recomputed through native K/G/M.
    pub whitened_residual_norm: f64,
    /// ||Sᵀ(K+G)phi|| + |lambda| ||Sᵀ M phi|| for the relative gate.
    pub whitened_residual_scale: f64,
    /// Euclidean norm of the same original-coordinate residual, diagnostic only.
    pub native_residual_norm: f64,
}

/// Complete numerically admitted packet at one frozen native configuration.
#[derive(Debug, Clone)]
pub struct CutPrestressModeReport3 {
    /// Ascending signed estimates; this does not certify spectral completeness.
    pub modes: Vec<CutPrestressMode3>,
    /// The requested search policy and gates used for this packet.
    pub options: CutPrestressModeOptions3,
    /// Spectral service iterations consumed.
    pub ticks: usize,
    /// Native K+G actions, including the independent residual checks.
    pub operator_applications: usize,
    /// Maximum entrywise error in phi_i^T M phi_j = delta_ij.
    pub maximum_mass_orthogonality_error: f64,
}

struct Scratch {
    physical: Vec<f64>,
    force: Vec<f64>,
    triangular: Vec<f64>,
    native: Vec<f64>,
    elastic: Vec<f64>,
    geometric: Vec<f64>,
    mass: Vec<f64>,
}

/// Bound native pencil and sparse consistent-mass factor, reusable for queries.
/// The immutable source borrow prevents scale changes under cached prestress.
/// No modal amplitude is interpreted as a physical strain-domain request.
pub struct CutPrestressModes3<'a> {
    source: CutElasticityRef3<'a>,
    prestress: &'a PreparedPrestress3,
    mass: &'a MaterialMass3,
    occupancy: &'a [f64],
    gate: &'a CancelGate,
    free: Vec<usize>,
    physical_dofs: usize,
    factor: LdltFactor,
    limits: CutPrestressLimits3,
    scratch: RefCell<Scratch>,
}

fn poll(gate: &CancelGate) -> Result<(), CutPrestressError3> {
    if gate.is_requested() {
        Err(CutPrestressError3::Cancelled)
    } else {
        Ok(())
    }
}
fn checkpoint(gate: &CancelGate) -> ControlFlow<()> {
    if gate.is_requested() {
        ControlFlow::Break(())
    } else {
        ControlFlow::Continue(())
    }
}
fn norm(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |n, value| n.hypot(*value))
}
fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

impl<'a> CutPrestressModes3<'a> {
    /// Prepare only the physical mass factor. Scalar native mass contributions
    /// are replicated/mirrored and summed by the existing sparse assembler;
    /// there is no stiffness probing, mass lumping or regularization.
    ///
    /// The complete sparse factor must have positive scalar pivots. Positive
    /// diagonals alone do not establish an admissible inertia metric. Sparse
    /// fill/front/work and numeric payload are bounded before publication.
    /// Mass sorting and the native stiffness action remain bounded indivisible
    /// kernels; cancellation polls surround them and run inside native G/M,
    /// sparse factor/triangular actions and each spectral tick.
    #[allow(clippy::too_many_lines)]
    pub fn new(
        source: CutElasticityRef3<'a>,
        prestress: &'a PreparedPrestress3,
        mass: &'a MaterialMass3,
        occupancy: &'a [f64],
        limits: CutPrestressLimits3,
        gate: &'a CancelGate,
    ) -> Result<Self, CutPrestressError3> {
        poll(gate)?;
        source.validate_mass(mass)?;
        let (cells, physical_dofs) = match source {
            CutElasticityRef3::Cartesian(op) => {
                op.validate_prestress(prestress)?;
                (op.cells(), source.operator().n())
            }
            CutElasticityRef3::Adaptive(op) => {
                op.validate_prestress(prestress)?;
                (
                    op.cells(),
                    op.physical_nodes()
                        .len()
                        .checked_mul(3)
                        .ok_or(CutPrestressError3::Budget("physical DOF count overflow"))?,
                )
            }
        };
        let native_n = source.operator().n();
        if native_n == 0
            || native_n > limits.max_native_dofs
            || physical_dofs > limits.max_native_dofs
            || cells > limits.max_cells
            || prestress.scalar_entries() > limits.max_geometric_entries
        {
            return Err(CutPrestressError3::Budget(
                "native prestress pencil exceeds its bounds",
            ));
        }
        if mass
            .scalar_terms()
            .checked_mul(6)
            .is_none_or(|n| n > limits.max_mass_entries)
        {
            return Err(CutPrestressError3::Budget(
                "consistent mass triplet allowance exhausted",
            ));
        }
        if limits.max_factor_bytes == 0 || limits.max_operator_applications == 0 {
            return Err(CutPrestressError3::Invalid(
                "positive factor and application allowances required",
            ));
        }
        // Bound diagonal, maps and action scratch before their allocation.
        // Native dimension is conservative when many coordinates are clamped.
        let base_work = physical_dofs
            .checked_mul(8)
            .and_then(|v| native_n.checked_mul(8).and_then(|w| v.checked_add(w)));
        if base_work.is_none_or(|v| v > limits.max_workspace_components) {
            return Err(CutPrestressError3::Budget(
                "native modal scratch allowance exhausted",
            ));
        }
        let diagonal = mass.diagonal(occupancy, || checkpoint(gate))?;
        let mut free = Vec::new();
        let mut index = vec![None; native_n];
        for (dof, &value) in diagonal.iter().enumerate() {
            if dof % 192 == 0 {
                poll(gate)?;
            }
            if !source.fixed()[dof / 3] {
                if value <= 0.0 || !value.is_finite() {
                    return Err(CutPrestressError3::Invalid(
                        "every free mode coordinate needs positive physical inertia",
                    ));
                }
                index[dof] = Some(free.len());
                free.push(dof);
            }
        }
        let n = free.len();
        if n == 0 {
            return Err(CutPrestressError3::Invalid(
                "prestress modes require a free coordinate",
            ));
        }
        let mut staged = Coo::new(n, n);
        for (entry, (cell, [a, b], reference)) in mass.reference_scalar_entries().enumerate() {
            if entry % 64 == 0 {
                poll(gate)?;
            }
            let value = occupancy[cell] * reference;
            if value == 0.0 {
                continue;
            }
            for component in 0..3 {
                if let (Some(i), Some(j)) = (index[3 * a + component], index[3 * b + component]) {
                    staged.push(i, j, value);
                    if i != j {
                        staged.push(j, i, value);
                    }
                }
            }
        }
        poll(gate)?;
        let matrix = staged.assemble();
        drop(staged);
        drop(index);
        drop(diagonal);
        poll(gate)?;
        let cancelled = || gate.is_requested();
        let mut factor_limits = limits.factor;
        factor_limits.max_front_dim = factor_limits.max_front_dim.min(n.saturating_sub(1));
        let symbolic = SymbolicLdlt::analyze_bounded(
            &matrix,
            DirectOrdering::Amd,
            &factor_limits,
            &cancelled,
        )?;
        let factor = symbolic.factor_bounded(
            &matrix,
            &LdltOptions {
                allow_2x2: false,
                max_factor_bytes: Some(limits.max_factor_bytes),
                ..Default::default()
            },
            &factor_limits,
            &cancelled,
        )?;
        // Admit the entire factor now. Each query prepares its O(n) inverse
        // square roots once, without self-referential ownership or refactoring.
        PositiveWhitening::new(&factor, &cancelled)?;
        poll(gate)?;
        Ok(Self {
            source,
            prestress,
            mass,
            occupancy,
            gate,
            free,
            physical_dofs,
            factor,
            limits,
            scratch: RefCell::new(Scratch {
                physical: vec![0.0; n],
                force: vec![0.0; n],
                triangular: vec![0.0; n],
                native: vec![0.0; native_n],
                elastic: vec![0.0; native_n],
                geometric: vec![0.0; native_n],
                mass: vec![0.0; native_n],
            }),
        })
    }

    /// Sorted native independent coordinates retained in the physical pencil.
    #[must_use]
    pub fn free_dofs(&self) -> &[usize] {
        &self.free
    }

    /// Query the smallest algebraic Ritz values of the fixed pencil. The
    /// existing LOBPCG owner handles indefinite operators and bounds its own
    /// scalar work between polls. No inverse of K+G or positivity is assumed.
    /// Returned vectors pass both the native mass Gram gate and a recomputed
    /// mass-whitened residual of the ORIGINAL K/G/M actions. Failure publishes
    /// no partial spectrum and does not mutate the native owners or factor.
    #[allow(clippy::too_many_lines)]
    pub fn lowest_modes(
        &self,
        options: CutPrestressModeOptions3,
    ) -> Result<CutPrestressModeReport3, CutPrestressError3> {
        poll(self.gate)?;
        let n = self.free.len();
        let tolerances = [
            options.spectral_tolerance,
            options.residual_relative_tolerance,
            options.residual_absolute_tolerance,
            options.mass_orthogonality_tolerance,
        ];
        if options.modes == 0
            || options.max_ticks == 0
            || tolerances.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || options.modes.checked_mul(3).is_none_or(|v| v > n)
        {
            return Err(CutPrestressError3::Invalid(
                "positive tolerances/ticks and 3*modes <= free DOFs required",
            ));
        }
        // Above the native 24*n*k + 27*k^2 block estimate, also account for
        // native reconstruction, retained scratch, whitening and result modes.
        let workspace = n
            .checked_mul(options.modes)
            .and_then(|v| v.checked_mul(32))
            .and_then(|v| {
                options
                    .modes
                    .checked_mul(options.modes)
                    .and_then(|w| w.checked_mul(36))
                    .and_then(|w| v.checked_add(w))
            })
            .and_then(|v| {
                options
                    .modes
                    .checked_add(8)
                    .and_then(|m| self.physical_dofs.checked_mul(m))
                    .and_then(|w| v.checked_add(w))
            })
            .and_then(|v| n.checked_mul(12).and_then(|w| v.checked_add(w)));
        if workspace.is_none_or(|v| v > self.limits.max_workspace_components) {
            return Err(CutPrestressError3::Budget(
                "LOBPCG and native modal workspace allowance exhausted",
            ));
        }
        let cancelled = || self.gate.is_requested();
        let whitening = PositiveWhitening::new(&self.factor, &cancelled)?;
        let mut service = EigenService::new(
            EigenBackend::Lobpcg,
            n,
            EigenQuery {
                k: options.modes,
                largest: false,
                tol: options.spectral_tolerance,
                steps_per_tick: 1,
            },
        )?;
        let operator = WhitenedOp {
            problem: self,
            whitening: &whitening,
            applications: Cell::new(0),
            failure: RefCell::new(None),
        };
        let mut last = CutPrestressError3::Unresolved {
            ticks: 0,
            whitened_residual: f64::INFINITY,
            mass_orthogonality_error: f64::INFINITY,
        };
        for tick in 0..options.max_ticks {
            poll(self.gate)?;
            let progress = service.tick_cancellable(&operator, &cancelled);
            if let Some(error) = operator.failure.borrow_mut().take() {
                return Err(error);
            }
            let progress = progress?;
            if progress.converged || tick + 1 == options.max_ticks {
                match operator.physical_result(&progress.pairs, options, progress.ticks) {
                    Ok(report) => return Ok(report),
                    Err(error @ CutPrestressError3::Unresolved { .. }) => last = error,
                    Err(error) => return Err(error),
                }
            }
        }
        Err(last)
    }
}

struct WhitenedOp<'a, 'source, 'factor> {
    problem: &'a CutPrestressModes3<'source>,
    whitening: &'a PositiveWhitening<'factor>,
    applications: Cell<usize>,
    failure: RefCell<Option<CutPrestressError3>>,
}

impl WhitenedOp<'_, '_, '_> {
    fn native_stiffness(
        &self,
        input: &[f64],
        elastic: &mut [f64],
        geometric: &mut [f64],
    ) -> Result<(), CutPrestressError3> {
        let problem = self.problem;
        poll(problem.gate)?;
        if self.applications.get() >= problem.limits.max_operator_applications {
            return Err(CutPrestressError3::Budget(
                "native stiffness application allowance exhausted",
            ));
        }
        self.applications.set(self.applications.get() + 1);
        problem.source.operator().apply(input, elastic);
        poll(problem.gate)?;
        problem
            .prestress
            .apply_into(input, geometric, || checkpoint(problem.gate))?;
        if elastic.iter().any(|v| !v.is_finite()) {
            return Err(CutPrestressError3::Invalid(
                "nonfinite native elastic action",
            ));
        }
        Ok(())
    }

    #[allow(clippy::too_many_lines)]
    fn physical_result(
        &self,
        pairs: &[CertifiedEigenvalue],
        options: CutPrestressModeOptions3,
        ticks: usize,
    ) -> Result<CutPrestressModeReport3, CutPrestressError3> {
        let problem = self.problem;
        let cancelled = || problem.gate.is_requested();
        let mut modes: Vec<CutPrestressMode3> = Vec::with_capacity(options.modes);
        let mut gram_error = 0.0_f64;
        let mut worst_residual = 0.0_f64;
        if pairs.len() != options.modes {
            return Err(CutPrestressError3::Unresolved {
                ticks,
                whitened_residual: f64::INFINITY,
                mass_orthogonality_error: f64::INFINITY,
            });
        }
        let spectral_residual = pairs.iter().fold(0.0_f64, |r, pair| r.max(pair.residual()));
        if spectral_residual > options.spectral_tolerance {
            return Err(CutPrestressError3::Unresolved {
                ticks,
                whitened_residual: spectral_residual,
                mass_orthogonality_error: f64::INFINITY,
            });
        }
        for pair in pairs {
            poll(problem.gate)?;
            let mut scratch = problem.scratch.borrow_mut();
            let Scratch {
                physical,
                force,
                triangular,
                native,
                elastic,
                geometric,
                mass,
            } = &mut *scratch;
            self.whitening
                .apply(pair.vector(), physical, triangular, &cancelled)?;
            let amplitude = physical.iter().fold(0.0_f64, |scale, v| scale.max(v.abs()));
            if !amplitude.is_finite() || amplitude <= 0.0 {
                return Err(CutPrestressError3::Invalid("unrepresentable physical mode"));
            }
            native.fill(0.0);
            for (&dof, value) in problem.free.iter().zip(physical.iter()) {
                native[dof] = value / amplitude;
            }
            problem
                .mass
                .apply_into(problem.occupancy, native, mass, || checkpoint(problem.gate))?;
            let work = dot(native, mass);
            if !work.is_finite() || work <= 0.0 {
                return Err(CutPrestressError3::Invalid(
                    "mode has no representable positive mass norm",
                ));
            }
            let normalization = work.sqrt();
            for value in native.iter_mut() {
                *value /= normalization;
            }
            // Reapply after normalization: gates describe the actual rounded
            // vector being returned, not a rescaled stale operator image.
            problem
                .mass
                .apply_into(problem.occupancy, native, mass, || checkpoint(problem.gate))?;
            let mass_norm = dot(native, mass);
            if !mass_norm.is_finite() || mass_norm <= 0.0 {
                return Err(CutPrestressError3::Invalid(
                    "normalized physical mode overflow",
                ));
            }
            gram_error = gram_error.max((mass_norm - 1.0).abs());
            for previous in &modes {
                let cross = dot(&previous.vector, mass);
                if !cross.is_finite() {
                    return Err(CutPrestressError3::Invalid("mode mass Gram overflow"));
                }
                gram_error = gram_error.max(cross.abs());
            }
            self.native_stiffness(native, elastic, geometric)?;
            let eigenvalue = (dot(native, elastic) + dot(native, geometric)) / mass_norm;
            if !eigenvalue.is_finite() {
                return Err(CutPrestressError3::Invalid(
                    "physical modal Rayleigh quotient overflow",
                ));
            }
            for (value, &dof) in force.iter_mut().zip(&problem.free) {
                *value = elastic[dof] + geometric[dof];
            }
            self.whitening
                .apply_transpose(force, physical, triangular, &cancelled)?;
            let stiffness_norm = norm(physical);
            for (value, &dof) in force.iter_mut().zip(&problem.free) {
                *value = mass[dof];
            }
            self.whitening
                .apply_transpose(force, physical, triangular, &cancelled)?;
            let residual_scale = stiffness_norm + eigenvalue.abs() * norm(physical);
            for (value, &dof) in force.iter_mut().zip(&problem.free) {
                *value = (-eigenvalue).mul_add(mass[dof], elastic[dof] + geometric[dof]);
            }
            let native_residual_norm = norm(force);
            self.whitening
                .apply_transpose(force, physical, triangular, &cancelled)?;
            let residual = norm(physical);
            let allowed = options.residual_absolute_tolerance
                + options.residual_relative_tolerance * residual_scale;
            if !residual.is_finite() || !residual_scale.is_finite() || !allowed.is_finite() {
                return Err(CutPrestressError3::Invalid("modal residual scale overflow"));
            }
            worst_residual = worst_residual.max(residual);
            if residual > allowed || gram_error > options.mass_orthogonality_tolerance {
                return Err(CutPrestressError3::Unresolved {
                    ticks,
                    whitened_residual: worst_residual,
                    mass_orthogonality_error: gram_error,
                });
            }
            modes.push(CutPrestressMode3 {
                eigenvalue,
                vector: native.clone(),
                whitened_residual_norm: residual,
                whitened_residual_scale: residual_scale,
                native_residual_norm,
            });
        }
        modes.sort_by(|a, b| a.eigenvalue.total_cmp(&b.eigenvalue));
        poll(problem.gate)?;
        Ok(CutPrestressModeReport3 {
            modes,
            options,
            ticks,
            operator_applications: self.applications.get(),
            maximum_mass_orthogonality_error: gram_error,
        })
    }
}

impl SymmetricOp for WhitenedOp<'_, '_, '_> {
    fn dim(&self) -> usize {
        self.problem.free.len()
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) {
        if self.failure.borrow().is_some() {
            output.fill(f64::NAN);
            return;
        }
        let result = (|| -> Result<(), CutPrestressError3> {
            let problem = self.problem;
            let cancelled = || problem.gate.is_requested();
            let mut scratch = problem.scratch.borrow_mut();
            let Scratch {
                physical,
                force,
                triangular,
                native,
                elastic,
                geometric,
                ..
            } = &mut *scratch;
            self.whitening
                .apply(input, physical, triangular, &cancelled)?;
            native.fill(0.0);
            for (&dof, &value) in problem.free.iter().zip(physical.iter()) {
                native[dof] = value;
            }
            self.native_stiffness(native, elastic, geometric)?;
            for (value, &dof) in force.iter_mut().zip(&problem.free) {
                *value = elastic[dof] + geometric[dof];
            }
            self.whitening
                .apply_transpose(force, output, triangular, &cancelled)?;
            Ok(())
        })();
        if let Err(error) = result {
            *self.failure.borrow_mut() = Some(error);
            output.fill(f64::NAN);
        }
    }
}
