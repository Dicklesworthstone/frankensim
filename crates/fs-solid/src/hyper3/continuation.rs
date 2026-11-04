//! Pseudo-arclength paths on the original three-dimensional hyperelastic solid.
//!
//! The native prepared tetrahedra own reference geometry, constitutive domains,
//! internal force and dP/dF. The existing bordered continuation solver owns
//! prediction, correction and step control. Only free coordinates enter its
//! system: `q = u / displacement_m`, `r = (f_int - lambda f_ref) / force_n`.
//! The reduced Jacobian is exactly `K L / F`.
//!
//! All supports are homogeneous and remain fixed; the load factor scales the
//! original dead nodal and reference-volume body forces. A state records
//! dimensionless FREE coordinates in [`HyperTetContinuation::free_dofs`] order.
//! Use [`HyperTetContinuation::response`] for physical displacement, reactions
//! and energy. A resumed state must satisfy this same problem's force equation.
//! Each advance returns a complete candidate and cannot mutate its input,
//! including when a later requested step or cancellation refuses.
//!
//! This adds no follower loading, moving supports, thermal/history evolution,
//! mixed incompressibility or spatial stability certificate. P1 locking and
//! discretization limitations remain those of the original solid. Path events
//! are the generic solver's diagnostics, not certified bifurcation locations.

use std::cell::RefCell;

use fs_exec::Cx;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_sparse::Csr;

use super::prepared::PreparedHyperTetProblem;
use super::{HyperTetError, HyperTetProblem, checkpoint};
use crate::SolidError;
use crate::continuation::{self, ArcSettings, PathResidual, PathState};

/// Positive physical scales defining the dimensionless continuation equation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HyperTetPathScales {
    /// Displacement coordinate scale, m; free coordinates are u / this value.
    pub displacement_m: f64,
    /// Force-residual scale, N; `ArcSettings::tol` multiplies this scale.
    pub force_n: f64,
}

/// Preserve the original native-solid refusal separately from path failures.
#[derive(Clone, Debug, PartialEq)]
pub enum HyperTetContinuationError {
    /// Original mesh, constitutive point, assembly or cancellation diagnosis.
    Solid(HyperTetError),
    /// Continuation input, corrector or bordered linear-solver diagnosis.
    Path(SolidError),
}

impl From<HyperTetError> for HyperTetContinuationError {
    fn from(error: HyperTetError) -> Self {
        Self::Solid(error)
    }
}
impl From<SolidError> for HyperTetContinuationError {
    fn from(error: SolidError) -> Self {
        Self::Path(error)
    }
}
impl core::fmt::Display for HyperTetContinuationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Solid(error) => error.fmt(f),
            Self::Path(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for HyperTetContinuationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(match self {
            Self::Solid(error) => error,
            Self::Path(error) => error,
        })
    }
}

fn invalid(what: &str) -> HyperTetContinuationError {
    SolidError::InvalidInput { what: what.into() }.into()
}

/// Physical observation of one path state on the original complete mesh.
#[derive(Clone, Debug, PartialEq)]
pub struct HyperTetPathResponse {
    /// Multiplier of the original dead reference loads.
    pub load_factor: f64,
    /// Node-major displacement in metres, including exact zero supports.
    pub displacement_m: Vec<[f64; 3]>,
    /// Internal minus applied force on original prescribed DOFs, N.
    pub reactions_n: Vec<(usize, f64)>,
    /// Physical force residual in `free_dofs()` order, N.
    pub free_force_residual_n: Vec<f64>,
    /// Original native stored energy, J.
    pub strain_energy_j: f64,
    /// Smallest actual det F among the original tetrahedra.
    pub min_det_f: f64,
}

/// Immutable original material/geometry binding and exact support elimination.
#[derive(Debug)]
pub struct HyperTetContinuation<'a, L: HyperelasticMaterial = Hyperelastic> {
    prepared: PreparedHyperTetProblem<'a, L>,
    prescribed: &'a [(usize, f64)],
    free_dofs: Vec<usize>,
    full_to_free: Vec<usize>,
    reference_force_n: Vec<f64>,
    reference_force_scaled: Vec<f64>,
    scales: HyperTetPathScales,
    tangent_scale: f64,
}

impl<'a, L: HyperelasticMaterial> HyperTetProblem<'a, L> {
    /// Prepare a conservative dead-load path with fixed homogeneous supports.
    /// Native admission and budgets apply. Scales and their L/F ratio must be
    /// positive finite representable values; the scaled free load must be
    /// finite and nonzero. Inelastic incremental potentials do not implement
    /// this entry's `HyperelasticMaterial` bound.
    pub fn prepare_continuation(
        &self,
        scales: HyperTetPathScales,
        cx: &Cx<'_>,
    ) -> Result<HyperTetContinuation<'a, L>, HyperTetContinuationError> {
        checkpoint(cx)?;
        let tangent_scale = scales.displacement_m / scales.force_n;
        if ![scales.displacement_m, scales.force_n, tangent_scale]
            .iter()
            .all(|value| value.is_finite() && *value > 0.0)
        {
            return Err(invalid(
                "positive finite displacement, force and L/F scales required",
            ));
        }
        if self.prescribed_m.iter().any(|&(_, value)| value != 0.0) {
            return Err(invalid(
                "three-dimensional continuation requires homogeneous fixed supports",
            ));
        }
        let prepared = self.prepare(cx)?;
        let reference_force_n =
            prepared.external_force_n(self.nodal_forces_n, self.body_force_n_m3, cx)?;
        let mut full_to_free = vec![0; reference_force_n.len()];
        for &(dof, _) in self.prescribed_m {
            full_to_free[dof] = usize::MAX;
        }
        let mut free_dofs = Vec::new();
        for (dof, index) in full_to_free.iter_mut().enumerate() {
            if *index != usize::MAX {
                *index = free_dofs.len();
                free_dofs.push(dof);
            }
        }
        let reference_force_scaled: Vec<_> = free_dofs
            .iter()
            .map(|&dof| reference_force_n[dof] / scales.force_n)
            .collect();
        if reference_force_scaled
            .iter()
            .any(|value| !value.is_finite())
            || !reference_force_scaled.iter().any(|&value| value != 0.0)
        {
            return Err(invalid(
                "finite nonzero normalized dead force on free coordinates required",
            ));
        }
        checkpoint(cx)?;
        Ok(HyperTetContinuation {
            prepared,
            prescribed: self.prescribed_m,
            free_dofs,
            full_to_free,
            reference_force_n,
            reference_force_scaled,
            scales,
            tangent_scale,
        })
    }
}

impl<L: HyperelasticMaterial> HyperTetContinuation<'_, L> {
    /// Original node-major xyz DOFs represented by the path coordinates.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }

    /// Physical scales retained throughout this path.
    pub const fn scales(&self) -> HyperTetPathScales {
        self.scales
    }

    /// Original complete dead load, including loads applied at fixed supports.
    pub fn reference_force_n(&self) -> &[f64] {
        &self.reference_force_n
    }

    fn displacement(&self, q: &[f64], lambda: f64) -> Result<Vec<f64>, HyperTetContinuationError> {
        if q.len() != self.free_dofs.len()
            || !lambda.is_finite()
            || q.iter().any(|value| !value.is_finite())
            || self
                .reference_force_n
                .iter()
                .any(|force| !(lambda * force).is_finite())
        {
            return Err(invalid(
                "finite dimension-matched path coordinates and representable applied loads required",
            ));
        }
        let mut u = vec![0.0; self.full_to_free.len()];
        for (&dof, &value) in self.free_dofs.iter().zip(q) {
            u[dof] = value * self.scales.displacement_m;
            if !u[dof].is_finite() {
                return Err(invalid("physical path displacement is not representable"));
            }
        }
        Ok(u)
    }

    /// Original free-force equation divided by the retained force scale.
    pub fn residual(
        &self,
        q: &[f64],
        lambda: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, HyperTetContinuationError> {
        checkpoint(cx)?;
        let u = self.displacement(q, lambda)?;
        let internal = self.prepared.internal_force(&u, cx)?;
        let result: Vec<_> = self
            .free_dofs
            .iter()
            .zip(&self.reference_force_scaled)
            .map(|(&dof, force)| internal[dof] / self.scales.force_n - lambda * force)
            .collect();
        if result.iter().any(|value| !value.is_finite()) {
            return Err(invalid("normalized free-force residual is not finite"));
        }
        checkpoint(cx)?;
        Ok(result)
    }

    /// Exact derivative of [`Self::residual`] with respect to q.
    /// Fixed rows/columns are eliminated; no support identities or reaction
    /// forces participate in path normalization or the arclength metric.
    pub fn tangent(
        &self,
        q: &[f64],
        lambda: f64,
        cx: &Cx<'_>,
    ) -> Result<Csr, HyperTetContinuationError> {
        checkpoint(cx)?;
        let u = self.displacement(q, lambda)?;
        let full = self.prepared.tangent_matrix(&u, cx)?;
        let n = self.free_dofs.len();
        let mut rows = Vec::with_capacity(n + 1);
        let mut columns = Vec::new();
        let mut values = Vec::new();
        rows.push(0);
        for (row, &dof) in self.free_dofs.iter().enumerate() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            let (full_columns, full_values) = full.row(dof);
            for (&column, &value) in full_columns.iter().zip(full_values) {
                let free_column = self.full_to_free[column];
                if free_column != usize::MAX {
                    let scaled = value * self.tangent_scale;
                    if !scaled.is_finite() {
                        return Err(invalid("normalized native tangent is not finite"));
                    }
                    columns.push(free_column);
                    values.push(scaled);
                }
            }
            rows.push(values.len());
        }
        // Original rows are canonical and the free map is increasing. Restrict
        // directly in O(nnz), without another COO buffer, duplicate reduction
        // or sort at every predictor/corrector evaluation.
        let tangent =
            Csr::try_from_parts_with_checkpoint(n, n, rows, columns, values, || checkpoint(cx))?
                .ok_or_else(|| SolidError::InternalInvariant {
                    what: "canonical native tangent lost ordering during free-DOF restriction"
                        .into(),
                })?;
        checkpoint(cx)?;
        Ok(tangent)
    }

    /// Start at the unloaded reference configuration, admitted as equilibrium.
    pub fn start(
        &self,
        settings: &ArcSettings,
        cx: &Cx<'_>,
    ) -> Result<PathState, HyperTetContinuationError> {
        self.start_from(&vec![0.0; self.full_to_free.len()], 0.0, settings, cx)
    }

    /// Start from a supplied physical equilibrium without projecting it.
    /// `displacement_m` is the complete node-major xyz field. Supports must be
    /// exactly zero and the original normalized free-force norm must meet tol.
    /// The initial predictor follows increasing lambda; a retained checkpoint
    /// instead preserves its own tangent through [`Self::advance`].
    pub fn start_from(
        &self,
        displacement_m: &[f64],
        lambda: f64,
        settings: &ArcSettings,
        cx: &Cx<'_>,
    ) -> Result<PathState, HyperTetContinuationError> {
        checkpoint(cx)?;
        if displacement_m.len() != self.full_to_free.len()
            || displacement_m.iter().any(|value| !value.is_finite())
            || self
                .prescribed
                .iter()
                .any(|&(dof, _)| displacement_m[dof] != 0.0)
        {
            return Err(invalid(
                "finite complete initial displacement with exact homogeneous supports required",
            ));
        }
        let mut state = PathState::start(self.free_dofs.len(), settings);
        state.u = self
            .free_dofs
            .iter()
            .map(|&dof| displacement_m[dof] / self.scales.displacement_m)
            .collect();
        state.lambda = lambda;
        let deflection = state
            .u
            .iter()
            .fold(0.0_f64, |max, value| max.max(value.abs()));
        state.trace = vec![(lambda, deflection)];
        self.advance(&state, settings, 0, cx)
    }

    /// Return a complete candidate after the requested number of arc steps.
    /// A failure at any step leaves `previous` untouched. Before resuming, the
    /// supplied checkpoint is checked against the original force equation and
    /// material domain. Tightening tol can therefore require a better starting
    /// equilibrium, rather than silently projecting an old state.
    pub fn advance(
        &self,
        previous: &PathState,
        settings: &ArcSettings,
        steps: usize,
        cx: &Cx<'_>,
    ) -> Result<PathState, HyperTetContinuationError> {
        checkpoint(cx)?;
        let adapter = PathAdapter {
            owner: self,
            cx,
            native_error: RefCell::new(None),
        };
        let mut next = previous.clone();
        continuation::advance(&adapter, &mut next, settings, 0)
            .map_err(|error| adapter.recover(error))?;
        let residual = self.residual(&previous.u, previous.lambda, cx)?;
        let norm = residual
            .iter()
            .fold(0.0_f64, |norm, &value| norm.hypot(value));
        if norm > settings.tol {
            return Err(invalid(
                "initial path state does not satisfy the original scaled force equation",
            ));
        }
        continuation::advance(&adapter, &mut next, settings, steps)
            .map_err(|error| adapter.recover(error))?;
        checkpoint(cx)?;
        Ok(next)
    }

    /// Recover physical fields with the original native material evaluation.
    /// This observation does not assert that a caller-mutated state converged;
    /// the returned free-force residual remains explicit.
    pub fn response(
        &self,
        state: &PathState,
        cx: &Cx<'_>,
    ) -> Result<HyperTetPathResponse, HyperTetContinuationError> {
        checkpoint(cx)?;
        let u = self.displacement(&state.u, state.lambda)?;
        let evaluated = self.prepared.evaluate(&u, cx)?;
        let residual: Vec<_> = evaluated
            .internal
            .iter()
            .zip(&self.reference_force_n)
            .map(|(internal, force)| internal - state.lambda * force)
            .collect();
        if residual.iter().any(|value| !value.is_finite()) {
            return Err(invalid("physical path force residual is not finite"));
        }
        let response = HyperTetPathResponse {
            load_factor: state.lambda,
            displacement_m: u.chunks_exact(3).map(|v| [v[0], v[1], v[2]]).collect(),
            reactions_n: self
                .prescribed
                .iter()
                .map(|&(dof, _)| (dof, residual[dof]))
                .collect(),
            free_force_residual_n: self.free_dofs.iter().map(|&dof| residual[dof]).collect(),
            strain_energy_j: evaluated.energy,
            min_det_f: evaluated.min_det_f,
        };
        checkpoint(cx)?;
        Ok(response)
    }
}

/// Scope-local error bridge; a final material refusal retains its element ID.
struct PathAdapter<'s, 'a, 'c, 'cx, L: HyperelasticMaterial> {
    owner: &'s HyperTetContinuation<'a, L>,
    cx: &'c Cx<'cx>,
    native_error: RefCell<Option<HyperTetError>>,
}

impl<L: HyperelasticMaterial> PathAdapter<'_, '_, '_, '_, L> {
    fn bridge<T>(&self, result: Result<T, HyperTetContinuationError>) -> Result<T, SolidError> {
        *self.native_error.borrow_mut() = None;
        result.map_err(|error| match error {
            HyperTetContinuationError::Path(error) => error,
            HyperTetContinuationError::Solid(error) => {
                *self.native_error.borrow_mut() = Some(error.clone());
                match error {
                    HyperTetError::Cancelled => SolidError::Cancelled,
                    HyperTetError::MaterialRefused { element, what } => {
                        SolidError::MaterialRefused {
                            what: format!("element {element}: {what}"),
                        }
                    }
                    other => SolidError::InvalidInput {
                        what: other.to_string(),
                    },
                }
            }
        })
    }

    fn recover(&self, error: SolidError) -> HyperTetContinuationError {
        self.native_error.borrow_mut().take().map_or(
            HyperTetContinuationError::Path(error),
            HyperTetContinuationError::Solid,
        )
    }
}

impl<L: HyperelasticMaterial> PathResidual for PathAdapter<'_, '_, '_, '_, L> {
    fn checkpoint(&self) -> Result<(), SolidError> {
        self.bridge(checkpoint(self.cx).map_err(Into::into))
    }
    fn ndof(&self) -> usize {
        self.owner.free_dofs.len()
    }
    fn residual(&self, u: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        self.bridge(self.owner.residual(u, lambda, self.cx))
    }
    fn tangent(&self, u: &[f64], lambda: f64) -> Result<Csr, SolidError> {
        self.bridge(self.owner.tangent(u, lambda, self.cx))
    }
    fn load_vector(&self) -> Vec<f64> {
        self.owner.reference_force_scaled.clone()
    }
}
