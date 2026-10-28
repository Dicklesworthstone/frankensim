//! Capillary material-surface motion on the existing structural midpoint clock.
//!
//! Fixed reference nodal masses and triangular material connectivity supply
//! inertia; the original Gram area discrete gradient supplies surface stress.
//! `fs-time` owns Newmark kinematics and its matrix-free Newton--Krylov solve.
//! With rho_inf=1 its stage is the midpoint: delta x=dt*v_mean and
//! M*delta v=dt*(F_surface+F_external). Their contraction closes kinetic plus
//! surface free energy, with prescribed coefficient supply accounted separately.
//!
//! This is an inertial tension-only material surface, not bulk incompressible
//! liquid or an elastic-solid constitutive replacement. Mass never follows
//! current area. The caller owns which physical interfaces each face represents
//! (a two-sided film needs both tensions), supports, and the validity of the
//! fixed topology. No bending, viscosity, contact, remeshing, surface entropy or
//! inferred thermal heat is supplied. Time/mesh refinement is still required.

use super::*;
use fs_time::galpha::{
    ImplicitSolveConfig, ImplicitStepTelemetry, OperatorGeneralizedAlpha,
    SecondOrderOperatorWeights, SecondOrderProblem, SecondOrderState, TimeSolveError,
};
use std::cell::{Cell, RefCell};

/// Nonlinear controls and independently checked physical acceptance bounds.
#[derive(Clone, Copy, Debug)]
pub struct SurfaceDynamicsControls {
    /// Original structural solver policy; its residual has force units [N].
    pub solve: ImplicitSolveConfig,
    /// Allocation cap on the original FGMRES restart dimension.
    pub maximum_krylov_restart: usize,
    /// Euclidean norm of all nodal impulse residuals [N s].
    pub impulse_tolerance_n_s: f64,
    /// Euclidean norm of midpoint displacement defects [m].
    pub kinematic_tolerance_m: f64,
    /// Absolute kinetic+surface-minus-supply residual [J].
    pub energy_tolerance_j: f64,
}

/// A complete accepted step, including the original structural solver report.
#[derive(Clone, Debug)]
pub struct SurfaceDynamicsReport {
    /// The native accepted time/step/Newton telemetry.
    pub native: ImplicitStepTelemetry,
    /// Original geometric energy, work and prescribed-sigma supply.
    pub surface: SurfaceEnergyChange,
    /// Initial kinetic storage [J].
    pub kinetic_before_j: f64,
    /// Accepted kinetic storage [J].
    pub kinetic_after_j: f64,
    /// Kinetic change evaluated from momentum increments [J].
    pub kinetic_change_j: f64,
    /// Held external nodal forces times actual displacement [J].
    pub external_work_j: f64,
    /// Fresh momentum balance norm, independent of the Newton report [N s].
    pub impulse_residual_n_s: f64,
    /// Fresh midpoint kinematic balance norm [m].
    pub kinematic_residual_m: f64,
    /// Delta kinetic + delta surface - external work - sigma supply [J].
    pub energy_residual_j: f64,
}

/// Refusal publishes no position, velocity, acceleration, clock or history.
#[derive(Debug)]
pub enum SurfaceDynamicsError {
    /// Invalid dimensions, finite data, topology budget or solve controls.
    Invalid(&'static str),
    /// Original geometric surface-force refusal.
    Surface(CapillarySurfaceError),
    /// Original structural/Newton--Krylov refusal.
    Time(TimeSolveError),
    /// Cancellation at a bounded callback or face boundary.
    Cancelled,
    /// The converged trial did not satisfy all independent physical gates.
    Balance {
        /// Nodal impulse residual norm [N s].
        impulse_n_s: f64,
        /// Midpoint displacement residual norm [m].
        kinematic_m: f64,
        /// Kinetic plus surface energy residual [J].
        energy_j: f64,
    },
}
impl core::fmt::Display for SurfaceDynamicsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "capillary surface dynamics: {self:?}")
    }
}
impl std::error::Error for SurfaceDynamicsError {}
impl From<CapillarySurfaceError> for SurfaceDynamicsError {
    fn from(value: CapillarySurfaceError) -> Self {
        Self::Surface(value)
    }
}

/// Borrowed fixed material connectivity and conserved reference nodal masses.
///
/// Coordinates and velocities remain in the existing `SecondOrderState`, in
/// Cartesian xyz order. No parallel trajectory or mesh representation is kept.
#[derive(Clone, Copy, Debug)]
pub struct CapillarySurfaceDynamics<'a> {
    masses_kg: &'a [f64],
    faces: &'a [[usize; 3]],
}

impl<'a> CapillarySurfaceDynamics<'a> {
    /// Bind existing triangle indices and strictly positive reference masses.
    /// Explicit topology caps bound every matrix-free callback's work.
    pub fn try_new(
        masses_kg: &'a [f64],
        faces: &'a [[usize; 3]],
        maximum_vertices: usize,
        maximum_faces: usize,
    ) -> Result<Self, SurfaceDynamicsError> {
        if masses_kg.len() < 3
            || faces.is_empty()
            || masses_kg.len() > maximum_vertices
            || faces.len() > maximum_faces
            || masses_kg.len().checked_mul(3).is_none()
            || masses_kg.iter().any(|m| !m.is_finite() || *m <= 0.0)
            || !masses_kg.iter().sum::<f64>().is_finite()
        {
            return Err(SurfaceDynamicsError::Invalid(
                "bounded nonempty topology and finite positive reference masses required",
            ));
        }
        for (face, indices) in faces.iter().enumerate() {
            for &vertex in indices {
                if vertex >= masses_kg.len() {
                    return Err(CapillarySurfaceError::InvalidVertexIndex { face, vertex }.into());
                }
            }
            if indices[0] == indices[1] || indices[1] == indices[2] || indices[2] == indices[0] {
                return Err(CapillarySurfaceError::DegenerateTriangle.into());
            }
        }
        Ok(Self { masses_kg, faces })
    }

    /// Initialize the original structural state with physical force/mass
    /// acceleration. Later accelerations retain native Newmark semantics;
    /// the accepted position and velocity satisfy the midpoint balance.
    pub fn initial_state(
        &self,
        time_s: f64,
        positions_m: &[[f64; 3]],
        velocities_m_s: &[[f64; 3]],
        tension: SurfaceTensionField<'_>,
        external_forces_n: &[[f64; 3]],
    ) -> Result<SecondOrderState, SurfaceDynamicsError> {
        if !time_s.is_finite()
            || !vectors(positions_m, self.masses_kg.len())
            || !vectors(velocities_m_s, self.masses_kg.len())
            || !vectors(external_forces_n, self.masses_kg.len())
        {
            return Err(SurfaceDynamicsError::Invalid(
                "finite per-node initial positions, velocities, forces and time required",
            ));
        }
        let surface =
            assemble_surface_tension_step(positions_m, positions_m, self.faces, tension, tension)?;
        let acceleration: Vec<_> = surface
            .forces_n
            .iter()
            .zip(external_forces_n)
            .zip(self.masses_kg)
            .flat_map(|((a, b), m)| (0..3).map(move |d| (a[d] + b[d]) / m))
            .collect();
        if acceleration.iter().any(|x| !x.is_finite()) {
            return Err(SurfaceDynamicsError::Invalid(
                "initial acceleration is unrepresentable",
            ));
        }
        let q: Vec<_> = positions_m.iter().flatten().copied().collect();
        let v: Vec<_> = velocities_m_s.iter().flatten().copied().collect();
        kinetic(self.masses_kg, &v)?;
        Ok(SecondOrderState::new(time_s, &q, &v, &acceleration))
    }

    /// Advance the SAME structural trajectory through a conservative capillary
    /// step. Tensions are prescribed endpoint fields; their change is external
    /// free-energy supply, never thermal deposition. External forces are held
    /// over the interval and perform force-times-displacement work.
    ///
    /// The analytic endpoint tangent is composed with the native midpoint
    /// kinematics. Mass scaling preconditions the existing matrix-free Krylov
    /// solve. Face kernels allocate nothing; no dense mesh Jacobian is formed.
    /// All physical gates are re-evaluated before publishing, including when
    /// the native residual meets its relative tolerance. Prior history is not
    /// copied during a trial. Cancellation is polled at most every 64 faces.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        state: &mut SecondOrderState,
        tension_before: SurfaceTensionField<'_>,
        tension_after: SurfaceTensionField<'_>,
        external_forces_n: &[[f64; 3]],
        dt_s: f64,
        controls: SurfaceDynamicsControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SurfaceDynamicsReport, SurfaceDynamicsError> {
        if cancelled() {
            return Err(SurfaceDynamicsError::Cancelled);
        }
        let dimension = 3 * self.masses_kg.len();
        if state.q.len() != dimension
            || state.v.len() != dimension
            || state.a.len() != dimension
            || state
                .q
                .iter()
                .chain(&state.v)
                .chain(&state.a)
                .any(|x| !x.is_finite())
            || !vectors(external_forces_n, self.masses_kg.len())
            || !dt_s.is_finite()
            || dt_s <= 0.0
            || !(dt_s * dt_s).is_finite()
            || !(4.0 / (dt_s * dt_s)).is_finite()
            || controls.solve.max_newton_iterations == 0
            || controls.solve.newton.linear_restart > controls.maximum_krylov_restart
            || [
                controls.impulse_tolerance_n_s,
                controls.kinematic_tolerance_m,
                controls.energy_tolerance_j,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(SurfaceDynamicsError::Invalid(
                "finite complete state, positive representable step and bounded physical solve controls required",
            ));
        }
        controls
            .solve
            .newton
            .validate()
            .map_err(|e| SurfaceDynamicsError::Time(TimeSolveError::NewtonSetup(e)))?;
        tension_before.validate(self.faces.len())?;
        tension_after.validate(self.faces.len())?;
        let before: Vec<_> = state
            .q
            .chunks_exact(3)
            .map(|v| [v[0], v[1], v[2]])
            .collect();
        let adapter = MidpointSurface {
            model: self,
            initial: &state.q,
            tension_before,
            tension_after,
            cancelled: RefCell::new(&mut cancelled),
            stopped: Cell::new(false),
        };
        // Validate physical start geometry before any native trial, retaining
        // its typed geometric error. The final evaluation below uses endpoints.
        let _ = adapter.surface(&before, &before)?;
        let mut candidate = SecondOrderState {
            t: state.t,
            q: state.q.clone(),
            v: state.v.clone(),
            a: state.a.clone(),
            steps: state.steps,
            history: Vec::new(),
        };
        let forcing: Vec<_> = external_forces_n.iter().flatten().copied().collect();
        let method = OperatorGeneralizedAlpha::new(dimension, dt_s, 1.0, controls.solve);
        let native =
            method.step_controlled(&mut candidate, &adapter, &forcing, &mut || adapter.cancel());
        if adapter.cancel() {
            return Err(SurfaceDynamicsError::Cancelled);
        }
        let native = native.map_err(SurfaceDynamicsError::Time)?;
        let after: Vec<_> = candidate
            .q
            .chunks_exact(3)
            .map(|v| [v[0], v[1], v[2]])
            .collect();
        let surface = adapter.surface(&before, &after)?;
        let mut impulse = 0.0_f64;
        let mut kinematic = 0.0_f64;
        let mut change = Sum::default();
        let mut external = Sum::default();
        for i in 0..dimension {
            if i % 192 == 0 && adapter.cancel() {
                return Err(SurfaceDynamicsError::Cancelled);
            }
            let dx = candidate.q[i] - state.q[i];
            let dv = candidate.v[i] - state.v[i];
            let velocity = 0.5 * state.v[i] + 0.5 * candidate.v[i];
            let momentum = self.masses_kg[i / 3] * dv;
            impulse =
                impulse.hypot(momentum - dt_s * (surface.forces_n[i / 3][i % 3] + forcing[i]));
            kinematic = kinematic.hypot(dx - dt_s * velocity);
            change.add(velocity * momentum);
            external.add(forcing[i] * dx);
        }
        let change = change.value();
        let external = external.value();
        let mut energy = Sum::default();
        for v in [
            change,
            surface.energy.energy_change_j,
            -external,
            -surface.energy.coefficient_supply_j,
        ] {
            energy.add(v);
        }
        let energy = energy.value();
        if !impulse.is_finite()
            || !kinematic.is_finite()
            || !energy.is_finite()
            || impulse > controls.impulse_tolerance_n_s
            || kinematic > controls.kinematic_tolerance_m
            || energy.abs() > controls.energy_tolerance_j
        {
            return Err(SurfaceDynamicsError::Balance {
                impulse_n_s: impulse,
                kinematic_m: kinematic,
                energy_j: energy,
            });
        }
        let report = SurfaceDynamicsReport {
            native,
            surface: surface.energy,
            kinetic_before_j: kinetic(self.masses_kg, &state.v)?,
            kinetic_after_j: kinetic(self.masses_kg, &candidate.v)?,
            kinetic_change_j: change,
            external_work_j: external,
            impulse_residual_n_s: impulse,
            kinematic_residual_m: kinematic,
            energy_residual_j: energy,
        };
        if adapter.cancel() {
            return Err(SurfaceDynamicsError::Cancelled);
        }
        state.q = candidate.q;
        state.v = candidate.v;
        state.a = candidate.a;
        state.t = candidate.t;
        state.steps = candidate.steps;
        state.history.extend(candidate.history);
        Ok(report)
    }
}

fn vectors(v: &[[f64; 3]], n: usize) -> bool {
    v.len() == n && v.iter().flatten().all(|x| x.is_finite())
}
fn kinetic(mass: &[f64], velocity: &[f64]) -> Result<f64, SurfaceDynamicsError> {
    let mut sum = Sum::default();
    for (v, m) in velocity.chunks_exact(3).zip(mass) {
        for &component in v {
            sum.add(0.5 * m * component * component);
        }
    }
    let value = sum.value();
    if !value.is_finite() {
        return Err(SurfaceDynamicsError::Invalid(
            "kinetic energy is unrepresentable",
        ));
    }
    Ok(value)
}

struct MidpointSurface<'a, 'm, C> {
    model: &'a CapillarySurfaceDynamics<'m>,
    initial: &'a [f64],
    tension_before: SurfaceTensionField<'a>,
    tension_after: SurfaceTensionField<'a>,
    cancelled: RefCell<&'a mut C>,
    stopped: Cell<bool>,
}
impl<C: FnMut() -> bool> MidpointSurface<'_, '_, C> {
    fn cancel(&self) -> bool {
        if !self.stopped.get() && (self.cancelled.borrow_mut())() {
            self.stopped.set(true);
        }
        self.stopped.get()
    }
    fn face(&self, indices: [usize; 3], midpoint: &[f64]) -> ([Vector; 3], [Vector; 3]) {
        let before = indices.map(|i| core::array::from_fn(|d| self.initial[3 * i + d]));
        let after = indices.map(|i| {
            core::array::from_fn(|d| {
                self.initial[3 * i + d] + 2.0 * (midpoint[3 * i + d] - self.initial[3 * i + d])
            })
        });
        (before, after)
    }
    fn surface(
        &self,
        before: &[Vector],
        after: &[Vector],
    ) -> Result<CapillaryPatchResponse, SurfaceDynamicsError> {
        let mut forces = vec![[0.0; 3]; before.len()];
        let mut energy = [Sum::default(); 7];
        for (face, &indices) in self.model.faces.iter().enumerate() {
            if face % 64 == 0 && self.cancel() {
                return Err(SurfaceDynamicsError::Cancelled);
            }
            let response = triangle_surface_tension_step(
                indices.map(|i| before[i]),
                indices.map(|i| after[i]),
                self.tension_before.at(face),
                self.tension_after.at(face),
            )?;
            for (node, force) in indices.into_iter().zip(response.forces_n) {
                for d in 0..3 {
                    forces[node][d] += force[d];
                }
            }
            for (sum, value) in energy.iter_mut().zip(response.energy.values()) {
                sum.add(value);
            }
        }
        let energy = energy.map(Sum::value);
        if forces
            .iter()
            .flatten()
            .chain(energy.iter())
            .any(|x| !x.is_finite())
        {
            return Err(CapillarySurfaceError::Unrepresentable.into());
        }
        Ok(CapillaryPatchResponse {
            forces_n: forces,
            energy: SurfaceEnergyChange::from_values(energy),
        })
    }
}
impl<C: FnMut() -> bool> SecondOrderProblem for MidpointSurface<'_, '_, C> {
    fn dimension(&self) -> usize {
        self.initial.len()
    }
    fn validate_state(&self, q: &[f64], v: &[f64]) -> Result<(), String> {
        if q.iter().chain(v).any(|x| !x.is_finite()) {
            Err("nonfinite surface position or velocity".into())
        } else {
            Ok(())
        }
    }
    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        for (i, (out, v)) in output.iter_mut().zip(input).enumerate() {
            *out = self.model.masses_kg[i / 3] * v;
        }
    }
    fn damping_apply(&self, _: &[f64], output: &mut [f64]) {
        output.fill(0.0);
    }
    fn internal_force(&self, q: &[f64], output: &mut [f64]) {
        output.fill(0.0);
        for (face, &indices) in self.model.faces.iter().enumerate() {
            if face % 64 == 0 && self.cancel() {
                output.fill(f64::NAN);
                return;
            }
            let (before, after) = self.face(indices, q);
            let Ok(response) = triangle_surface_tension_step(
                before,
                after,
                self.tension_before.at(face),
                self.tension_after.at(face),
            ) else {
                output.fill(f64::NAN);
                return;
            };
            for (node, force) in indices.into_iter().zip(response.forces_n) {
                for d in 0..3 {
                    output[3 * node + d] -= force[d];
                }
            }
        }
    }
    fn resistance(&self, q: &[f64], _: &[f64], output: &mut [f64]) {
        self.internal_force(q, output);
    }
    fn tangent_apply(&self, q: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(0.0);
        for (face, &indices) in self.model.faces.iter().enumerate() {
            if face % 64 == 0 && self.cancel() {
                output.fill(f64::NAN);
                return;
            }
            let (before, after) = self.face(indices, q);
            let direction = indices.map(|i| core::array::from_fn(|d| 2.0 * direction[3 * i + d]));
            let Ok(tangent) = triangle_surface_tension_step_tangent(
                before,
                after,
                self.tension_before.at(face),
                self.tension_after.at(face),
                direction,
            ) else {
                output.fill(f64::NAN);
                return;
            };
            for (node, force) in indices.into_iter().zip(tangent) {
                for d in 0..3 {
                    output[3 * node + d] -= force[d];
                }
            }
        }
    }
    fn preconditioner_apply(
        &self,
        _: &[f64],
        weights: SecondOrderOperatorWeights,
        _: usize,
        _: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        for (i, (out, r)) in output.iter_mut().zip(residual).enumerate() {
            *out = r / (weights.mass * self.model.masses_kg[i / 3]);
        }
    }
}

#[cfg(test)]
mod tests;
