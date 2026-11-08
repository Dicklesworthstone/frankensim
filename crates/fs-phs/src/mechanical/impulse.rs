//! Instantaneous original-coordinate loading without resetting material memory.

use super::LinearMechanicalSystem;
use crate::interconnection::accumulation::Sum;
use crate::{PhsError, storage_energy_change};

/// Admission of the actual stored event, independent of background momentum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MechanicalImpulseControls {
    /// Componentwise absolute allowance in original generalized impulse units.
    pub absolute_momentum_tolerance: f64,
    /// Relative allowance against max(|requested impulse|, |stored change|).
    /// Must lie in [0,1); the pre-existing momentum does not enlarge it.
    pub relative_momentum_tolerance: f64,
    /// Absolute allowance on constitutive storage change minus supplied work.
    /// This is joules when original coordinates/forces use conjugate SI units.
    pub work_tolerance: f64,
}

/// Accepted zero-duration event in the unchanged original [q,p,z] coordinates.
#[derive(Clone, Debug)]
pub struct MechanicalImpulse {
    /// Complete candidate; displacement and every Maxwell history retain bits.
    pub state_after: Vec<f64>,
    /// Full admitted inverse mass applied to the actual momentum change.
    pub velocity_change: Vec<f64>,
    /// Actual difference of the stored canonical momenta.
    pub momentum_change: Vec<f64>,
    /// Actual momentum change minus the supplied original-coordinate impulse.
    pub momentum_residual: Vec<f64>,
    /// Constitutive storage change, evaluated without subtracting total energies.
    /// Elastic and memory states are fixed, so this is the kinetic increment.
    pub kinetic_energy_change: f64,
    /// Signed impulse work from average physical velocity; braking is negative.
    pub supplied_work: f64,
    /// Work of each original port; empty for original-coordinate impulses.
    /// The total retains low words before rounding these individual observables.
    pub port_work: Vec<f64>,
    /// Exactly zero: finite viscosity and Maxwell relaxation act over time.
    pub dissipated: f64,
    /// Independently evaluated constitutive increment minus supplied work.
    pub work_residual: f64,
}

impl LinearMechanicalSystem {
    /// Apply simultaneous prescribed impulses in original generalized coordinates.
    ///
    /// At fixed q and material history, p+ = p- + J and v+ - v- = M^-1 J.
    /// Finite stiffness, viscosity and gyroscopic forces supply no impulse or
    /// heat at zero duration. Work is J^T (v-+v+)/2, including extraction.
    /// Neither physical time nor forcing phase is advanced. Continue an existing
    /// exact/prepared owner from the returned complete state at the event time.
    ///
    /// This applies already resolved loading; it does not infer contact,
    /// restitution, constraint reactions or event time. Do not divide an impulse
    /// by a sample period or replay it as a finite-duration force.
    ///
    /// The candidate is returned only after actual momentum and independent
    /// storage/work admission. Refusal never changes the model or caller state.
    /// Dense O(n^2) contractions use the prepared full inverse; there is no new
    /// mass factorization or integrator. Intermediate arithmetic must be finite.
    pub fn apply_impulse(
        &self,
        state: &[f64],
        impulse: &[f64],
        controls: MechanicalImpulseControls,
    ) -> Result<MechanicalImpulse, PhsError> {
        if impulse.len() != self.coordinates || impulse.iter().any(|v| !v.is_finite()) {
            return Err(invalid("complete finite original mechanical impulse"));
        }
        self.impulse_with(state, controls, |i| Sum::dot(&[impulse[i]], &[1.0]), None)
    }

    /// Apply integrated inputs through the original reciprocal force map P.
    ///
    /// All ports act simultaneously: J=P*J_port. Each signed port work is
    /// J_port[p] times its average conjugate velocity, including cross-port
    /// terms. Opposing large sources retain a weak net impulse/work through
    /// compensated contractions. They are never applied as sequential kicks.
    /// State, history, work and refusal semantics match [`Self::apply_impulse`].
    pub fn apply_port_impulse(
        &self,
        state: &[f64],
        impulse: &[f64],
        controls: MechanicalImpulseControls,
    ) -> Result<MechanicalImpulse, PhsError> {
        let ports = self.system.port_dim();
        if impulse.len() != ports || impulse.iter().any(|v| !v.is_finite()) {
            return Err(invalid("complete finite mechanical port impulse"));
        }
        let (_, _, map) = self.system.structure();
        self.impulse_with(
            state,
            controls,
            |i| {
                let row = (self.coordinates + i) * ports;
                Sum::dot(&map[row..row + ports], impulse)
            },
            Some(impulse),
        )
    }

    fn impulse_with(
        &self,
        state: &[f64],
        controls: MechanicalImpulseControls,
        source: impl Fn(usize) -> Option<Sum>,
        port_impulse: Option<&[f64]>,
    ) -> Result<MechanicalImpulse, PhsError> {
        self.admit_state(state)?;
        if !controls.absolute_momentum_tolerance.is_finite()
            || controls.absolute_momentum_tolerance < 0.0
            || !controls.relative_momentum_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.relative_momentum_tolerance)
            || !controls.work_tolerance.is_finite()
            || controls.work_tolerance < 0.0
        {
            return Err(invalid("finite mechanical impulse admission controls"));
        }
        let n = self.coordinates;
        let mut state_after = state.to_vec();
        let mut requested = Vec::with_capacity(n);
        let mut changes = Vec::with_capacity(n);
        let mut means = Vec::with_capacity(n);
        let mut momentum_change = Vec::with_capacity(n);
        let mut momentum_residual = Vec::with_capacity(n);
        for i in 0..n {
            let impulse = source(i).ok_or_else(arithmetic)?;
            let supplied = impulse.value().ok_or_else(arithmetic)?;
            // An exact zero event retains signed zero and every existing bit.
            if impulse.components().ok_or_else(arithmetic)? != (0.0, 0.0) {
                let mut next = impulse;
                next.product(state[n + i], 1.0).ok_or_else(arithmetic)?;
                state_after[n + i] = next.value().ok_or_else(arithmetic)?;
            }
            let mut change = Sum::default();
            change
                .product(state_after[n + i], 1.0)
                .ok_or_else(arithmetic)?;
            change.product(state[n + i], -1.0).ok_or_else(arithmetic)?;
            let actual = change.value().ok_or_else(arithmetic)?;
            let mut residual = change;
            residual.add_scaled(impulse, -1.0).ok_or_else(arithmetic)?;
            let residual = residual.value().ok_or_else(arithmetic)?;
            let allowance = controls.absolute_momentum_tolerance
                + controls.relative_momentum_tolerance * actual.abs().max(supplied.abs());
            if !allowance.is_finite() || residual.abs() > allowance {
                return Err(invalid(
                    "stored mechanical momentum cannot resolve the impulse",
                ));
            }
            let mut mean = Sum::default();
            if state[n + i] == state_after[n + i] {
                mean.product(state[n + i], 1.0).ok_or_else(arithmetic)?;
            } else {
                mean.half_product(state[n + i], 1.0)
                    .ok_or_else(arithmetic)?;
                mean.half_product(state_after[n + i], 1.0)
                    .ok_or_else(arithmetic)?;
            }
            requested.push(impulse);
            changes.push(change);
            means.push(mean);
            momentum_change.push(actual);
            momentum_residual.push(residual);
        }
        let mut velocity_change = Vec::with_capacity(n);
        let mut mean_velocity = Vec::with_capacity(n);
        let mut supplied = Sum::default();
        for i in 0..n {
            let mut change = Sum::default();
            let mut mean = Sum::default();
            for j in 0..n {
                let inverse = self.inverse_mass[i * n + j];
                change
                    .add_scaled(changes[j], inverse)
                    .ok_or_else(arithmetic)?;
                mean.add_scaled(means[j], inverse).ok_or_else(arithmetic)?;
            }
            if port_impulse.is_none() {
                supplied
                    .add_product(mean, requested[i])
                    .ok_or_else(arithmetic)?;
            }
            velocity_change.push(change.value().ok_or_else(arithmetic)?);
            mean_velocity.push(mean);
        }
        let mut port_work = Vec::new();
        if let Some(inputs) = port_impulse {
            let (_, _, map) = self.system.structure();
            port_work.reserve(inputs.len());
            for (port, &input) in inputs.iter().enumerate() {
                let mut output = Sum::default();
                for (i, &mean) in mean_velocity.iter().enumerate() {
                    output
                        .add_scaled(mean, map[(n + i) * inputs.len() + port])
                        .ok_or_else(arithmetic)?;
                }
                let mut work = Sum::default();
                work.add_scaled(output, input).ok_or_else(arithmetic)?;
                supplied.add_scaled(work, 1.0).ok_or_else(arithmetic)?;
                port_work.push(work.value().ok_or_else(arithmetic)?);
            }
        }
        let d = state.len();
        let kinetic_energy_change = storage_energy_change(
            self.system.storage(),
            state,
            &state_after,
            &mut vec![0.0; d],
            &mut vec![0.0; d],
        );
        let mut residual = Sum::default();
        residual
            .product(kinetic_energy_change, 1.0)
            .ok_or_else(arithmetic)?;
        residual.add_scaled(supplied, -1.0).ok_or_else(arithmetic)?;
        let work_residual = residual.value().ok_or_else(arithmetic)?;
        if work_residual.abs() > controls.work_tolerance {
            return Err(invalid(
                "mechanical impulse storage and original source work disagree",
            ));
        }
        Ok(MechanicalImpulse {
            state_after,
            velocity_change,
            momentum_change,
            momentum_residual,
            kinetic_energy_change,
            supplied_work: supplied.value().ok_or_else(arithmetic)?,
            port_work,
            dissipated: 0.0,
            work_residual,
        })
    }
}

fn invalid(what: &'static str) -> PhsError {
    PhsError::Unadmitted { what }
}
fn arithmetic() -> PhsError {
    invalid("mechanical impulse or work arithmetic is unrepresentable")
}
