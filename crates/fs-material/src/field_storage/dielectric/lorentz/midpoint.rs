//! Source-owned Lorentz condensation and recovery for coupled quasistatic fields.
//!
//! The complete original fixed-condition material is retained. Local poles
//! contribute positive endpoint permittivity and a history involving both P0
//! and J0. The field owner solves once, then recovers those same original states.
//! This does not add Ohmic conduction, thermal parameter work or a Maxwell solver.

use super::oscillator::{midpoint_field_energy, midpoint_field_increment};
use super::{
    FieldStorageError, LorentzError, LorentzState, PreparedLorentzMidpoint,
    ResolvedLorentzMaterial, finite, invalid,
};
use fs_math::dd::Dd;

/// Source-bank recovery for one midpoint field interval. All energies are
/// densities [J/m3]; the spatial owner retains the original geometry weights.
#[derive(Clone, Debug, PartialEq)]
pub struct LorentzMidpointStep<const N: usize> {
    /// Original P/J endpoint states, in the source material's pole order.
    pub states_after: Vec<LorentzState<N>>,
    pub displacement_before_c_m2: [f64; N],
    pub displacement_after_c_m2: [f64; N],
    pub stored_before_j_m3: f64,
    pub stored_after_j_m3: f64,
    /// Stable change of instantaneous field, restoring and kinetic storage.
    pub stored_change_j_m3: f64,
    /// Original midpoint electric work `Em dot (D1-D0)`.
    pub supplied_work_j_m3: f64,
    /// Independent nonnegative polarization friction; excludes Ohmic heat.
    pub heat_j_m3: f64,
    /// Stored change plus heat minus supplied work, without balancing a term.
    pub energy_residual_j_m3: f64,
}

/// Derived coefficient cache for one original material bank and fixed interval.
/// Preparation occurs once; field cells reuse it for histories and recovery.
/// It owns no cloned source record. The caller retains the complete material
/// and its identity alongside this immutable cache, as for prepared pole steps.
#[derive(Clone, Debug)]
pub struct PreparedLorentzMaterialMidpoint {
    epsilon_infinity: f64,
    epsilon_midpoint: f64,
    duration_s: f64,
    poles: Vec<PreparedLorentzMidpoint>,
}

impl PreparedLorentzMaterialMidpoint {
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    pub const fn midpoint_permittivity_f_m(&self) -> f64 {
        self.epsilon_midpoint
    }

    /// No coefficient preparation, allocation or source query in a warm cell.
    pub fn midpoint_history_c_m2<const N: usize>(
        &self,
        electric_before_v_m: [f64; N],
        states: &[LorentzState<N>],
    ) -> Result<[f64; N], FieldStorageError> {
        history(
            self.poles.iter().copied().map(Ok),
            electric_before_v_m,
            states,
        )
    }

    /// Recover with exactly the same cached coefficients as the field solve.
    /// Only the returned at-most-16-pole state vector allocates.
    pub fn step_midpoint<const N: usize>(
        &self,
        electric_before_v_m: [f64; N],
        electric_after_v_m: [f64; N],
        states: &[LorentzState<N>],
    ) -> Result<LorentzMidpointStep<N>, FieldStorageError> {
        recover(
            self.epsilon_infinity,
            self.poles.iter().copied().map(Ok),
            electric_before_v_m,
            electric_after_v_m,
            states,
        )
    }
}

fn interval(dt_s: f64) -> Result<(), FieldStorageError> {
    if dt_s.is_finite() && dt_s > 0.0 {
        Ok(())
    } else {
        Err(invalid(
            "positive finite Lorentz midpoint interval required",
        ))
    }
}

fn admit<const N: usize>(
    poles: usize,
    field: &[f64; N],
    states: &[LorentzState<N>],
) -> Result<(), FieldStorageError> {
    if !(1..=3).contains(&N)
        || states.len() != poles
        || field.iter().any(|v| !v.is_finite())
        || states.iter().any(|state| {
            state
                .polarization_c_m2
                .iter()
                .chain(&state.current_a_m2)
                .any(|v| !v.is_finite())
        })
    {
        return Err(invalid(
            "one finite 1D/2D/3D P/J state per original Lorentz pole required",
        ));
    }
    Ok(())
}

fn displacement<const N: usize>(
    epsilon: f64,
    electric: [f64; N],
    states: &[LorentzState<N>],
) -> Result<[f64; N], FieldStorageError> {
    let mut result = [0.0; N];
    for i in 0..N {
        let mut sum = Dd::from_f64(epsilon) * Dd::from_f64(electric[i]);
        for state in states {
            sum = sum + Dd::from_f64(state.polarization_c_m2[i]);
        }
        result[i] = finite(sum.to_f64())?;
    }
    Ok(result)
}

impl ResolvedLorentzMaterial {
    /// Total original displacement `epsilon_infinity*E + sum P` [C/m2].
    /// Current is retained in the state, but is not an additional displacement.
    pub fn displacement_c_m2<const N: usize>(
        &self,
        electric_v_m: [f64; N],
        states: &[LorentzState<N>],
    ) -> Result<[f64; N], FieldStorageError> {
        admit(self.poles.len(), &electric_v_m, states)?;
        displacement(self.epsilon_infinity, electric_v_m, states)
    }

    /// Full instantaneous field plus original polarization restoring/kinetic
    /// energy. It cannot be replaced by effective complex epsilon times E².
    pub fn stored_energy_density_j_m3<const N: usize>(
        &self,
        electric_v_m: [f64; N],
        states: &[LorentzState<N>],
    ) -> Result<f64, FieldStorageError> {
        admit(self.poles.len(), &electric_v_m, states)?;
        let mut sum = Dd::from_f64(midpoint_field_energy(self.epsilon_infinity, electric_v_m)?);
        for (pole, state) in self.poles.iter().zip(states) {
            sum = sum + Dd::from_f64(pole.stored_energy_j_m3(*state)?);
        }
        finite(sum.to_f64())
    }

    /// Prepare the bounded original bank once for repeated spatial cells and
    /// steps of this interval. No source snapshot clone, refresh or fitting.
    pub fn prepare_midpoint(
        &self,
        dt_s: f64,
    ) -> Result<PreparedLorentzMaterialMidpoint, FieldStorageError> {
        interval(dt_s)?;
        let mut poles = Vec::new();
        poles
            .try_reserve_exact(self.poles.len())
            .map_err(|_| invalid("Lorentz midpoint coefficient allocation failed"))?;
        let mut epsilon = Dd::from_f64(self.epsilon_infinity);
        for pole in &self.poles {
            let prepared = pole.prepare_midpoint(dt_s)?;
            epsilon = epsilon + Dd::from_f64(prepared.permittivity_increment_f_m()?);
            poles.push(prepared);
        }
        Ok(PreparedLorentzMaterialMidpoint {
            epsilon_infinity: self.epsilon_infinity,
            epsilon_midpoint: finite(epsilon.to_f64())?,
            duration_s: dt_s,
            poles,
        })
    }

    /// Positive endpoint coefficient after local midpoint condensation [F/m].
    /// Empty banks return epsilon_infinity exactly. This convenience query
    /// allocates nothing; repeated cell work should use `prepare_midpoint`.
    pub fn midpoint_permittivity_f_m(&self, dt_s: f64) -> Result<f64, FieldStorageError> {
        interval(dt_s)?;
        let mut sum = Dd::from_f64(self.epsilon_infinity);
        for pole in &self.poles {
            sum = sum + Dd::from_f64(pole.prepare_midpoint(dt_s)?.permittivity_increment_f_m()?);
        }
        finite(sum.to_f64())
    }

    /// Allocation-free history in `D1 = epsilon_midpoint*E1 + history`.
    /// Both P0 and J0 contribute. Recovery must use this same source and dt.
    /// Repeated cell work can reuse the source-owned prepared bank instead.
    pub fn midpoint_history_c_m2<const N: usize>(
        &self,
        dt_s: f64,
        electric_before_v_m: [f64; N],
        states: &[LorentzState<N>],
    ) -> Result<[f64; N], FieldStorageError> {
        interval(dt_s)?;
        history(
            self.poles.iter().map(|pole| pole.prepare_midpoint(dt_s)),
            electric_before_v_m,
            states,
        )
    }

    /// Recover the original bank after a coupled midpoint endpoint field solve.
    /// Only the at-most-16-pole output vector allocates. Failure leaves input
    /// states/source unchanged and publishes no partial bank. The spatial owner
    /// owns convergence, geometry integration, and energy-residual admission.
    /// Midpoint is passive and second order, not exact or L-stable; refinement
    /// must resolve resonance phase, damping and the physical field history.
    pub fn step_midpoint<const N: usize>(
        &self,
        dt_s: f64,
        electric_before_v_m: [f64; N],
        electric_after_v_m: [f64; N],
        states: &[LorentzState<N>],
    ) -> Result<LorentzMidpointStep<N>, FieldStorageError> {
        interval(dt_s)?;
        recover(
            self.epsilon_infinity,
            self.poles.iter().map(|pole| pole.prepare_midpoint(dt_s)),
            electric_before_v_m,
            electric_after_v_m,
            states,
        )
    }
}

fn history<const N: usize>(
    poles: impl ExactSizeIterator<Item = Result<PreparedLorentzMidpoint, LorentzError>>,
    electric_before: [f64; N],
    states: &[LorentzState<N>],
) -> Result<[f64; N], FieldStorageError> {
    admit(poles.len(), &electric_before, states)?;
    let mut history = [Dd::ZERO; N];
    for (pole, before) in poles.zip(states) {
        let contribution = pole?.history_c_m2(electric_before, *before)?;
        for i in 0..N {
            history[i] = history[i] + Dd::from_f64(contribution[i]);
        }
    }
    let mut result = [0.0; N];
    for i in 0..N {
        result[i] = finite(history[i].to_f64())?;
    }
    Ok(result)
}

fn recover<const N: usize>(
    epsilon: f64,
    poles: impl ExactSizeIterator<Item = Result<PreparedLorentzMidpoint, LorentzError>>,
    electric_before: [f64; N],
    electric_after: [f64; N],
    states: &[LorentzState<N>],
) -> Result<LorentzMidpointStep<N>, FieldStorageError> {
    admit(poles.len(), &electric_before, states)?;
    admit(poles.len(), &electric_after, states)?;
    let field_change = midpoint_field_increment(epsilon, electric_before, electric_after)?;
    let mut change = Dd::from_f64(field_change);
    let mut work = change;
    let mut heat = Dd::ZERO;
    let mut stored_before = Dd::from_f64(midpoint_field_energy(epsilon, electric_before)?);
    let mut stored_after = Dd::from_f64(midpoint_field_energy(epsilon, electric_after)?);
    let mut states_after = Vec::new();
    states_after
        .try_reserve_exact(poles.len())
        .map_err(|_| invalid("Lorentz midpoint state allocation failed"))?;
    for (pole, before) in poles.zip(states) {
        let step = pole?.advance(*before, electric_before, electric_after)?;
        change = change + Dd::from_f64(step.stored_change_j_m3);
        work = work + Dd::from_f64(step.supplied_work_j_m3);
        heat = heat + Dd::from_f64(step.heat_j_m3);
        stored_before = stored_before + Dd::from_f64(step.stored_before_j_m3);
        stored_after = stored_after + Dd::from_f64(step.stored_after_j_m3);
        states_after.push(step.state);
    }
    Ok(LorentzMidpointStep {
        displacement_before_c_m2: displacement(epsilon, electric_before, states)?,
        displacement_after_c_m2: displacement(epsilon, electric_after, &states_after)?,
        states_after,
        stored_before_j_m3: finite(stored_before.to_f64())?,
        stored_after_j_m3: finite(stored_after.to_f64())?,
        stored_change_j_m3: finite(change.to_f64())?,
        supplied_work_j_m3: finite(work.to_f64())?,
        heat_j_m3: finite(heat.to_f64())?,
        energy_residual_j_m3: finite((change + heat - work).to_f64())?,
    })
}
