//! Pseudo-arclength paths on the original finite-strain Taylor–Hood solid.
//!
//! Free displacements and independent pressures share the original mixed
//! potential, reference quadrature and consistent saddle tangent. With diagonal
//! physical variable scale D and work scale W = stress_scale * reference_volume,
//! q = D^-1 x, R = D r / W and R_q = D K D / W. The load parameter ramps the
//! original dead forces, moving supports, pressure gauge datum and any native
//! held pressure/gravity. Its negative derivative includes -D K_fc g / W;
//! pressure-constraint rows participate in this lifting as well.
//!
//! The native bordered solver crosses simple folds without inverting K or
//! requiring positive elastic inertia. Exact incompressibility remains the
//! original weak P1 constraint. Every accepted step checks ALL pressure rows,
//! including the omitted gauge equation, and validates the gauge null mode.
//! The mixed saddle potential is not relabeled as physical stored energy.
//!
//! No material/reference/history evolution, branch switching, mesh convergence,
//! injectivity or stability certificate is supplied. The saddle tangent is
//! normally indefinite even on a stable branch; its generic definiteness
//! diagnostic is not a mechanical bifurcation test. Only the load-turn events
//! have their original continuation meaning.

use std::cell::RefCell;

use fs_exec::Cx;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_material::hyper::mixed::MixedFiniteStrainPotential;
use fs_sparse::{Coo, Csr};

use super::{
    MixedUpError, MixedUpEvaluation, MixedUpProblem, PreparedMixedUpProblem, checkpoint, invalid,
    point_deformation, pressure,
};
use crate::SolidError;
use crate::continuation::{self, ArcSettings, PathState};

/// Preserve native material, geometry, gauge and cancellation diagnoses.
#[derive(Clone, Debug, PartialEq)]
pub enum MixedContinuationError {
    Mixed(MixedUpError),
    Path(SolidError),
}

impl From<MixedUpError> for MixedContinuationError {
    fn from(error: MixedUpError) -> Self {
        Self::Mixed(error)
    }
}
impl From<SolidError> for MixedContinuationError {
    fn from(error: SolidError) -> Self {
        Self::Path(error)
    }
}
impl core::fmt::Display for MixedContinuationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Mixed(e) => e.fmt(f),
            Self::Path(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for MixedContinuationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(match self {
            Self::Mixed(e) => e,
            Self::Path(e) => e,
        })
    }
}

/// Physical observation of a path point, without claiming its equilibrium.
#[derive(Clone, Debug, PartialEq)]
pub struct MixedPathResponse {
    pub load_factor: f64,
    pub displacement_m: Vec<[f64; 3]>,
    /// Original continuous corner-pressure order from pressure_nodes().
    pub pressure_pa: Vec<f64>,
    pub reactions_n: Vec<(usize, f64)>,
    pub free_force_residual_n: Vec<f64>,
    /// ALL original pressure equations, including a dropped gauge equation.
    pub pressure_residual_m3: Vec<f64>,
    /// Original mixed material saddle potential; external loads are excluded.
    pub mixed_potential_j: f64,
    pub min_det_f: f64,
}

struct CachedPoint {
    q: Vec<f64>,
    lambda: f64,
    displacement_m: Vec<[f64; 3]>,
    pressure_pa: Vec<f64>,
    physical: MixedUpEvaluation,
    residual: Vec<f64>,
    tangent: Csr,
    negative_parameter_derivative: Vec<f64>,
}

/// Immutable source-bound mixed problem with a shared exact-point cache.
/// Only hyperelastic owners are admitted; frozen incremental inelastic
/// potentials cannot be used as a history-independent equilibrium path.
pub struct MixedContinuation<
    'a,
    L: HyperelasticMaterial + MixedFiniteStrainPotential = Hyperelastic,
> {
    prepared: PreparedMixedUpProblem<'a, L>,
    prescribed_rate: Vec<f64>,
    cache: RefCell<Option<CachedPoint>>,
}

impl<'a, L: HyperelasticMaterial + MixedFiniteStrainPotential> MixedUpProblem<'a, L> {
    /// Prepare with the original length/stress/work scales and dead loads.
    pub fn prepare_continuation(
        &self,
        cx: &Cx<'_>,
    ) -> Result<MixedContinuation<'a, L>, MixedContinuationError> {
        self.prepare(cx)?.into_continuation(cx)
    }
}

impl<'a, L: HyperelasticMaterial + MixedFiniteStrainPotential> PreparedMixedUpProblem<'a, L> {
    /// Retain this original prepared geometry, constitutive field and held loads.
    /// Both support displacement and a supplied gauge value mean their values
    /// at lambda=1, as in the original mixed endpoint solver.
    pub fn into_continuation(
        self,
        cx: &Cx<'_>,
    ) -> Result<MixedContinuation<'a, L>, MixedContinuationError> {
        checkpoint(cx)?;
        let nu = 3 * self.mesh.nodes_count;
        let np = self.mesh.pressure_nodes.len();
        let mut path = MixedContinuation {
            prepared: self,
            prescribed_rate: Vec::new(),
            cache: RefCell::new(None),
        };
        path.admit_storage(0, 0, 0)?;
        path.prescribed_rate = vec![0.0; nu + np];
        for &(d, value) in path.prepared.problem.prescribed_m {
            path.prescribed_rate[d] = value;
        }
        if let (Some(pin), Some((_, value))) = (
            path.prepared.gauge_pressure,
            path.prepared.problem.pressure_gauge,
        ) {
            path.prescribed_rate[nu + pin] = value;
        }
        checkpoint(cx)?;
        Ok(path)
    }
}

impl<L: HyperelasticMaterial + MixedFiniteStrainPotential> MixedContinuation<'_, L> {
    /// Global mixed DOFs: original node-major xyz, then corner pressures.
    pub fn free_dofs(&self) -> &[usize] {
        &self.prepared.free
    }
    pub fn pressure_nodes(&self) -> &[usize] {
        &self.prepared.mesh.pressure_nodes
    }
    /// Original mesh/material/dead-load definition; held loads remain separate.
    pub fn problem(&self) -> &MixedUpProblem<'_, L> {
        &self.prepared.problem
    }
    pub fn boundary_pressure(&self) -> Option<super::pressure::UniformBoundaryPressure> {
        self.prepared.boundary_pressure
    }
    pub fn gravity(&self) -> Option<super::gravity::HeldMixedGravity<'_>> {
        self.prepared.gravity
    }
    pub fn work_scale_j(&self) -> f64 {
        self.prepared.work_scale_j
    }

    // Bounds retain the original payload caps. Account for a full assembly,
    // retained/restricted CSR and callback clones. The original path solver
    // uses <=1024 dense bordered coordinates or restart<=64, <=320-step FGMRES.
    fn admit_storage(&self, steps: usize, trace: usize, events: usize) -> Result<(), MixedUpError> {
        for scale in [
            self.prepared.problem.settings.length_scale_m,
            self.prepared.problem.settings.stress_scale_pa,
        ] {
            let row_scale = scale / self.prepared.work_scale_j;
            if !row_scale.is_finite() || row_scale <= 0.0 {
                return Err(invalid(
                    "mixed path work-conjugate row scales are unrepresentable",
                ));
            }
        }
        let n = self.prepared.free.len() as u128;
        let full = 3 * self.prepared.mesh.nodes_count as u128
            + self.prepared.mesh.pressure_nodes.len() as u128;
        let entries = self.prepared.mesh.elements.len() as u128 * 34 * 34;
        let assembly = entries * 192 + full * 512;
        if assembly > self.prepared.problem.settings.maximum_assembly_bytes as u128 {
            return Err(invalid(
                "mixed continuation assembly/cache exceeds maximum_assembly_bytes",
            ));
        }
        let dense = if n < 1024 { 24 * (n + 1) * (n + 1) } else { 0 };
        let krylov = 2048 * (n + 1) + 65_536;
        let records =
            (trace as u128 + steps as u128) * 64 + (events as u128 + 2 * steps as u128) * 128;
        let scratch = dense.max(krylov) + 512 * (n + 1) + records;
        if scratch > self.prepared.problem.settings.maximum_factor_bytes as u128 {
            return Err(invalid(
                "mixed continuation border/path exceeds maximum_factor_bytes",
            ));
        }
        Ok(())
    }

    fn physical_state(
        &self,
        q: &[f64],
        lambda: f64,
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 3]>, Vec<f64>), MixedUpError> {
        if q.len() != self.prepared.free.len()
            || !lambda.is_finite()
            || q.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "finite dimension-matched mixed path coordinates required",
            ));
        }
        let nu = self.prepared.mesh.nodes_count * 3;
        let mut u = vec![[0.0; 3]; self.prepared.mesh.nodes_count];
        let mut p = vec![0.0; self.prepared.mesh.pressure_nodes.len()];
        for (d, &rate) in self.prescribed_rate.iter().enumerate() {
            if d % 256 == 0 {
                checkpoint(cx)?;
            }
            let value = lambda * rate;
            if d < nu {
                u[d / 3][d % 3] = value;
            } else {
                p[d - nu] = value;
            }
        }
        for (i, (&d, &value)) in self.prepared.free.iter().zip(q).enumerate() {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            let physical = value * self.prepared.variable_scale(d);
            if d < nu {
                u[d / 3][d % 3] = physical;
            } else {
                p[d - nu] = physical;
            }
        }
        self.prepared.validate_state(&u, &p)?;
        if self
            .prepared
            .external_n
            .iter()
            .any(|f| !(lambda * f).is_finite())
        {
            return Err(invalid("mixed path applied force is unrepresentable"));
        }
        Ok((u, p))
    }

    // Exact negative explicit parameter derivative at fixed complete u,p.
    // Native pressure/gravity owners supply original current geometry and
    // quadrature. This visits loads only, never repeats constitutive assembly.
    fn reference_force(&self, u: &[[f64; 3]], cx: &Cx<'_>) -> Result<Vec<f64>, MixedUpError> {
        let prepared = &self.prepared;
        let mut force = prepared.external_n.clone();
        if prepared.boundary_pressure.is_some() || prepared.gravity.is_some() {
            for (e, element) in prepared.mesh.elements.iter().enumerate() {
                for (q, point) in element.points.iter().enumerate() {
                    checkpoint(cx)?;
                    if let Some(gravity) = prepared.gravity {
                        let (mass, sample) =
                            gravity.sample(e, q, element, point, prepared.problem.nodes_m, u)?;
                        let mut local = [0.0; 34];
                        super::gravity::add_point(&point.shapes, mass, sample, &mut local, None);
                        for a in 0..10 {
                            for i in 0..3 {
                                force[3 * element.connectivity[a] + i] -= local[3 * a + i];
                            }
                        }
                    }
                    if let Some(boundary) = prepared.boundary_pressure {
                        let f = point_deformation(element, point, u);
                        let mut piola = [0.0; 9];
                        pressure::add_volume_gradient(&f, boundary.pressure_pa(), &mut piola);
                        for a in 0..10 {
                            for i in 0..3 {
                                force[3 * element.connectivity[a] + i] -= point.weight_m3
                                    * (0..3)
                                        .map(|j| piola[3 * i + j] * point.gradients[a][j])
                                        .sum::<f64>();
                            }
                        }
                    }
                }
            }
        }
        if force.iter().any(|v| !v.is_finite()) {
            return Err(invalid("mixed path force derivative is unrepresentable"));
        }
        Ok(force)
    }

    fn ensure(&self, q: &[f64], lambda: f64, cx: &Cx<'_>) -> Result<(), MixedUpError> {
        checkpoint(cx)?;
        if self
            .cache
            .borrow()
            .as_ref()
            .is_some_and(|c| c.q == q && c.lambda == lambda)
        {
            return Ok(());
        }
        *self.cache.borrow_mut() = None;
        let (u, p) = self.physical_state(q, lambda, cx)?;
        let mut physical = self
            .prepared
            .evaluate_at_load(&u, &p, true, None, lambda, cx)?;
        let full = physical.tangent.take().expect("requested mixed tangent");
        self.prepared.validate_pressure_gauge(&full, cx)?;
        let force = self.reference_force(&u, cx)?;
        let nu = u.len() * 3;
        let n = self.prepared.free.len();
        let mut residual = Vec::with_capacity(n);
        let mut parameter = Vec::with_capacity(n);
        let mut coo = Coo::new(n, n);
        for (i, &d) in self.prepared.free.iter().enumerate() {
            if i % 128 == 0 {
                checkpoint(cx)?;
            }
            let row_scale = self.prepared.variable_scale(d) / self.prepared.work_scale_j;
            let raw = if d < nu {
                physical.internal_forces_n[d] - lambda * self.prepared.external_n[d]
            } else {
                physical.pressure_residual_m3[d - nu]
            };
            residual.push(row_scale * raw);
            let mut load = if d < nu { force[d] } else { 0.0 };
            let (cols, values) = full.row(d);
            for (&col, &value) in cols.iter().zip(values) {
                if let Some(j) = self.prepared.free_index[col] {
                    let entry = row_scale * value * self.prepared.variable_scale(col);
                    if !entry.is_finite() {
                        return Err(invalid("mixed path tangent is unrepresentable"));
                    }
                    coo.push(i, j, entry);
                } else {
                    load -= value * self.prescribed_rate[col];
                }
            }
            parameter.push(row_scale * load);
        }
        if residual.iter().chain(&parameter).any(|v| !v.is_finite()) {
            return Err(invalid(
                "mixed path residual/load column is unrepresentable",
            ));
        }
        let tangent = coo.assemble();
        checkpoint(cx)?;
        *self.cache.borrow_mut() = Some(CachedPoint {
            q: q.to_vec(),
            lambda,
            displacement_m: u,
            pressure_pa: p,
            physical,
            residual,
            tangent,
            negative_parameter_derivative: parameter,
        });
        Ok(())
    }

    pub fn residual(
        &self,
        q: &[f64],
        lambda: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, MixedContinuationError> {
        self.ensure(q, lambda, cx)?;
        Ok(self
            .cache
            .borrow()
            .as_ref()
            .expect("prepared point")
            .residual
            .clone())
    }
    pub fn tangent(
        &self,
        q: &[f64],
        lambda: f64,
        cx: &Cx<'_>,
    ) -> Result<Csr, MixedContinuationError> {
        self.ensure(q, lambda, cx)?;
        Ok(self
            .cache
            .borrow()
            .as_ref()
            .expect("prepared point")
            .tangent
            .clone())
    }
    /// Exact -dR/dlambda, including moving support and pressure-datum lifting.
    pub fn load_vector_at(
        &self,
        q: &[f64],
        lambda: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, MixedContinuationError> {
        self.ensure(q, lambda, cx)?;
        Ok(self
            .cache
            .borrow()
            .as_ref()
            .expect("prepared point")
            .negative_parameter_derivative
            .clone())
    }

    fn admit_equilibrium(
        &self,
        state: &PathState,
        tol: f64,
        cx: &Cx<'_>,
    ) -> Result<(), MixedContinuationError> {
        self.ensure(&state.u, state.lambda, cx)?;
        let cache = self.cache.borrow();
        let point = cache.as_ref().expect("prepared point");
        let norm = point.residual.iter().fold(0.0_f64, |n, &v| n.hypot(v));
        let pressure_scale =
            self.prepared.problem.settings.stress_scale_pa / self.prepared.work_scale_j;
        let pressure_norm = point
            .physical
            .pressure_residual_m3
            .iter()
            .fold(0.0_f64, |n, &v| n.hypot(v * pressure_scale));
        if pressure_norm > tol {
            if let Some((node, _)) = self.prepared.problem.pressure_gauge {
                return Err(MixedUpError::IncompatiblePressureGauge { node }.into());
            }
            return Err(
                invalid("mixed continuation pressure equilibrium exceeds tolerance").into(),
            );
        }
        if norm > tol {
            return Err(invalid("mixed continuation equilibrium exceeds tolerance").into());
        }
        Ok(())
    }

    /// Start at the original unloaded configuration and zero pressure.
    /// A thermally prestretched source may require start_from instead.
    pub fn start(
        &self,
        settings: &ArcSettings,
        cx: &Cx<'_>,
    ) -> Result<PathState, MixedContinuationError> {
        self.start_from(
            &vec![[0.0; 3]; self.prepared.mesh.nodes_count],
            &vec![0.0; self.prepared.mesh.pressure_nodes.len()],
            0.0,
            settings,
            cx,
        )
    }

    /// Start from a supplied physical equilibrium without projecting it.
    /// Authored supports and any gauge datum must equal lambda times their
    /// full-load values exactly; a resumed checkpoint retains its direction.
    pub fn start_from(
        &self,
        displacement_m: &[[f64; 3]],
        pressure_pa: &[f64],
        lambda: f64,
        settings: &ArcSettings,
        cx: &Cx<'_>,
    ) -> Result<PathState, MixedContinuationError> {
        checkpoint(cx)?;
        self.prepared.validate_state(displacement_m, pressure_pa)?;
        if !lambda.is_finite() {
            return Err(invalid("finite mixed path load factor required").into());
        }
        let nu = displacement_m.len() * 3;
        for (d, index) in self.prepared.free_index.iter().enumerate() {
            if d % 256 == 0 {
                checkpoint(cx)?;
            }
            if index.is_none() {
                let value = if d < nu {
                    displacement_m[d / 3][d % 3]
                } else {
                    pressure_pa[d - nu]
                };
                if value != lambda * self.prescribed_rate[d] {
                    return Err(invalid(
                        "initial mixed path support/gauge differs from the authored datum",
                    )
                    .into());
                }
            }
        }
        let mut state = PathState::start(self.prepared.free.len(), settings);
        for (i, (q, &d)) in state.u.iter_mut().zip(&self.prepared.free).enumerate() {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            let value = if d < nu {
                displacement_m[d / 3][d % 3]
            } else {
                pressure_pa[d - nu]
            };
            *q = value / self.prepared.variable_scale(d);
        }
        state.lambda = lambda;
        state.last_definite = false; // A pressure saddle is not an SPD elastic pencil.
        state.trace = vec![(lambda, state.u.iter().fold(0.0_f64, |m, v| m.max(v.abs())))];
        self.advance(&state, settings, 0, cx)
    }

    /// Return a complete candidate; failure/cancellation leaves previous intact.
    /// All pressure equations are checked after EACH accepted step, so a dropped
    /// gauge equation cannot hide an incompatible prescribed volume history.
    pub fn advance(
        &self,
        previous: &PathState,
        settings: &ArcSettings,
        steps: usize,
        cx: &Cx<'_>,
    ) -> Result<PathState, MixedContinuationError> {
        checkpoint(cx)?;
        if previous.u.len() != self.prepared.free.len()
            || previous.tangent.0.len() != self.prepared.free.len()
            || previous.step.checked_add(steps).is_none()
        {
            return Err(invalid(
                "dimension-matched mixed checkpoint and representable step count required",
            )
            .into());
        }
        self.admit_storage(steps, previous.trace.len(), previous.events.len())?;
        if previous.pending_switch {
            return Err(invalid(
                "mixed continuation does not admit an unverified branch-switch direction",
            )
            .into());
        }
        let adapter = Adapter {
            owner: self,
            cx,
            native_error: RefCell::new(None),
        };
        let mut next = previous.clone();
        continuation::advance(&adapter, &mut next, settings, 0).map_err(|e| adapter.recover(e))?;
        self.admit_equilibrium(&next, settings.tol, cx)?;
        for _ in 0..steps {
            continuation::advance(&adapter, &mut next, settings, 1)
                .map_err(|e| adapter.recover(e))?;
            self.admit_equilibrium(&next, settings.tol, cx)?;
        }
        checkpoint(cx)?;
        Ok(next)
    }

    pub fn response(
        &self,
        state: &PathState,
        cx: &Cx<'_>,
    ) -> Result<MixedPathResponse, MixedContinuationError> {
        self.ensure(&state.u, state.lambda, cx)?;
        let cache = self.cache.borrow();
        let point = cache.as_ref().expect("prepared point");
        let nu = point.displacement_m.len() * 3;
        let net = |d: usize| {
            point.physical.internal_forces_n[d] - state.lambda * self.prepared.external_n[d]
        };
        Ok(MixedPathResponse {
            load_factor: state.lambda,
            displacement_m: point.displacement_m.clone(),
            pressure_pa: point.pressure_pa.clone(),
            reactions_n: self
                .prepared
                .problem
                .prescribed_m
                .iter()
                .map(|&(d, _)| (d, net(d)))
                .collect(),
            free_force_residual_n: self
                .prepared
                .free
                .iter()
                .copied()
                .take_while(|&d| d < nu)
                .map(net)
                .collect(),
            pressure_residual_m3: point.physical.pressure_residual_m3.clone(),
            mixed_potential_j: point.physical.potential_j,
            min_det_f: point.physical.min_det_f,
        })
    }
}

struct Adapter<'s, 'a, 'c, 'cx, L: HyperelasticMaterial + MixedFiniteStrainPotential> {
    owner: &'s MixedContinuation<'a, L>,
    cx: &'c Cx<'cx>,
    native_error: RefCell<Option<MixedUpError>>,
}
impl<L: HyperelasticMaterial + MixedFiniteStrainPotential> Adapter<'_, '_, '_, '_, L> {
    fn bridge<T>(&self, result: Result<T, MixedContinuationError>) -> Result<T, SolidError> {
        *self.native_error.borrow_mut() = None;
        result.map_err(|e| match e {
            MixedContinuationError::Path(e) => e,
            MixedContinuationError::Mixed(e) => {
                *self.native_error.borrow_mut() = Some(e.clone());
                match e {
                    MixedUpError::Cancelled => SolidError::Cancelled,
                    MixedUpError::Material {
                        element,
                        point,
                        error,
                    } => SolidError::MaterialRefused {
                        what: format!("element {element}, point {point}: {error}"),
                    },
                    other => SolidError::InvalidInput {
                        what: other.to_string(),
                    },
                }
            }
        })
    }
    fn recover(&self, error: SolidError) -> MixedContinuationError {
        self.native_error.borrow_mut().take().map_or(
            MixedContinuationError::Path(error),
            MixedContinuationError::Mixed,
        )
    }
}
impl<L: HyperelasticMaterial + MixedFiniteStrainPotential> continuation::PathResidual
    for Adapter<'_, '_, '_, '_, L>
{
    fn checkpoint(&self) -> Result<(), SolidError> {
        self.bridge(checkpoint(self.cx).map_err(Into::into))
    }
    fn ndof(&self) -> usize {
        self.owner.free_dofs().len()
    }
    fn residual(&self, q: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        self.bridge(self.owner.residual(q, lambda, self.cx))
    }
    fn tangent(&self, q: &[f64], lambda: f64) -> Result<Csr, SolidError> {
        self.bridge(self.owner.tangent(q, lambda, self.cx))
    }
    fn load_vector(&self) -> Vec<f64> {
        vec![0.0; self.ndof()]
    }
    fn load_vector_at(&self, q: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        self.bridge(self.owner.load_vector_at(q, lambda, self.cx))
    }
}
