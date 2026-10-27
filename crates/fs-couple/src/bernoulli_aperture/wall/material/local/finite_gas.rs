//! A finite thermochemical background and the continuing native DG wall owner.
//!
//! One frozen-initial-medium acoustic sample is followed by a conservative
//! material refresh. The existing scalar feedback solves the work transferred
//! from the finite background into the actual wave/diffusion storage. No gas
//! kinetics, mass-flow law, thermal solver or acoustic propagator is duplicated.

use super::{MaterialWallDg, RegionalMaterialWallDg};
use super::update::{GasRefresh, prepare_gas};
use crate::acoustic_realize::AcousticRealizeError;
use crate::cavity_phs::piston::thermal::feedback;
use fs_duct::dg::{DgDuct, DgGasMaterialUpdate, DgPropagator, DgState, DgStepOutput, PortInput};
use fs_material::gas_field::{AxialGasField, PrimitiveGasState};
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::composition::GasCompositionStep;
use fs_matdb::QueryPoint;

/// Select the original realized loss owner; this never introduces or drops memory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FiniteGasDgLoss {
    /// No gas diffusion memory; existing structural and boundary losses remain active.
    Lossless,
    /// Existing uniform circular diffusion realization, retaining its pole identities.
    UniformViscothermal { maximum_relative_error: f64 },
}

/// Explicit physical work, partition-step and nonlinear acceptance limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGasDgSolve {
    /// Signed work INTO the acoustic system [J]. Ordered and containing zero.
    pub work_bracket_j: [f64; 2],
    /// Allowed difference between gas work and actual acoustic parameter work [J].
    pub work_tolerance_j: f64,
    /// Absolute gas-owner AND combined endpoint energy tolerance [J].
    pub energy_tolerance_j: f64,
    pub maximum_iterations: usize,
    /// Limits freezing the medium during the acoustic sample; not an error estimate.
    pub maximum_temperature_increment_k: f64,
    pub loss: FiniteGasDgLoss,
}

/// One original gas/thermal/flow owner's complete candidate at a trial work.
// The callback owns its physical process and never mutates accepted external owners.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasDgBackground<S> {
    pub gas: ClosedGasParcel,
    /// External heat into this gas inventory [J], excluding transported enthalpy.
    pub supplied_heat_j: f64,
    /// Net inlet minus outlet enthalpy [J], from actual original mass-flow accounts.
    /// Zero for a closed inventory; never counted again as heat.
    pub transported_enthalpy_j: f64,
    /// Retained original source report and any uncommitted companion endpoints.
    pub source: S,
}

/// A failure publishes none of the gas, material binding, state or stepping owner.
#[derive(Debug)]
pub enum FiniteGasDgError<E> {
    Invalid(&'static str),
    Background(E),
    Acoustic(AcousticRealizeError),
    OutsideBracket { lower_residual_j: f64, upper_residual_j: f64 },
    NotConverged { iterations: usize, residual_j: f64 },
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for FiniteGasDgError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-gas DG coupling: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for FiniteGasDgError<E> {}

/// One accepted physical sample. All instantaneous observations use the final gas.
#[derive(Debug)]
pub struct FiniteGasDgStep<S> {
    /// Original exact held-port impulses/work, with final endpoint pressure and
    /// ledger updated after the zero-duration material transfer.
    pub acoustic: DgStepOutput,
    /// Actual native storage transfer, including every retained diffusion memory.
    pub refresh: DgGasMaterialUpdate,
    pub background: S,
    pub supplied_heat_j: f64,
    pub transported_enthalpy_j: f64,
    /// Actual extensive gas storage change, allowing conserved native mass transfer.
    pub gas_energy_change_j: f64,
    /// Work debited by the gas owner; its difference from refresh.parameter_work_j
    /// is explicitly bounded below rather than hidden as heat.
    pub work_by_gas_j: f64,
    pub work_residual_j: f64,
    pub gas_energy_residual_j: f64,
    /// Original unpartitioned acoustic/wall/boundary loss. It leaves this combined
    /// accounting boundary; it is NOT automatically deposited into the gas.
    pub acoustic_dissipation_j: f64,
    /// Measured gas + native storage + original loss - heat - enthalpy - port work.
    pub combined_energy_residual_j: f64,
    pub trials: usize,
}

struct Candidate<S> {
    duct: DgDuct,
    state: DgState,
    propagator: DgPropagator,
    gas: ClosedGasParcel,
    step: FiniteGasDgStep<S>,
}

fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn check<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), FiniteGasDgError<E>> {
    if cancelled() { Err(FiniteGasDgError::Cancelled) } else { Ok(()) }
}
fn acoustic<E>(error: AcousticRealizeError) -> FiniteGasDgError<E> {
    match error {
        AcousticRealizeError::Cancelled => FiniteGasDgError::Cancelled,
        error => FiniteGasDgError::Acoustic(error),
    }
}
fn native<E>(error: fs_duct::DuctError) -> FiniteGasDgError<E> {
    match error {
        fs_duct::DuctError::Cancelled { .. } => FiniteGasDgError::Cancelled,
        error => acoustic(AcousticRealizeError::Nonlinear(error.to_string())),
    }
}

/// Difference specific energies before multiplying, preserving formation references.
fn stored_change(old: &ClosedGasParcel, next: &ClosedGasParcel) -> f64 {
    let u = old.state().specific_internal_energy_j_kg();
    next.mass_kg() * (next.state().specific_internal_energy_j_kg() - u)
        + (next.mass_kg() - old.mass_kg()) * u
}

fn field<E>(original: &AxialGasField, gas: &ClosedGasParcel)
    -> Result<AxialGasField, FiniteGasDgError<E>>
{
    let mut context = QueryPoint::new();
    for (name, value) in gas.point().axes() {
        if name == "temperature" || name == "pressure" { continue; }
        context = match gas.point().axis_quantities().get(name) {
            Some(quantity) => context.with_quantity(name.clone(), quantity.clone(), *value),
            None => context.with(name.clone(), *value),
        }.map_err(|_| FiniteGasDgError::Invalid("gas source coordinates cannot be retained"))?;
    }
    let primitive = PrimitiveGasState {
        temperature_k: gas.state().state().temperature,
        composition: gas.model().thermodynamics().mole_composition().clone(),
    };
    let result = AxialGasField::new(
        original.stations_m().to_vec(), gas.state().state().pressure,
        vec![primitive; original.stations_m().len()], gas.model().clone(),
        context, original.phase(),
    ).map_err(|error| acoustic(AcousticRealizeError::Nonlinear(error.to_string())))?;
    if result.uniform_state() != Some(gas.state()) {
        return Err(FiniteGasDgError::Invalid(
            "field must resolve the identical complete endpoint gas and phase/source policy",
        ));
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn prepare<S, E, C: FnMut() -> bool>(
    duct: &DgDuct, state: &DgState, propagator: &DgPropagator,
    gas: &ClosedGasParcel, ports: &[PortInput], controls: FiniteGasDgSolve,
    mut evolve: impl FnMut(&ClosedGasParcel, f64, f64, &mut C)
        -> Result<FiniteGasDgBackground<S>, E>,
    cancelled: &mut C,
) -> Result<Candidate<S>, FiniteGasDgError<E>> {
    check(cancelled)?;
    let [lo, hi] = controls.work_bracket_j;
    if !lo.is_finite() || !hi.is_finite() || lo > hi || lo > 0.0 || hi < 0.0
        || !positive(controls.work_tolerance_j) || !positive(controls.energy_tolerance_j)
        || !(1..=128).contains(&controls.maximum_iterations)
        || !positive(controls.maximum_temperature_increment_k)
    {
        return Err(FiniteGasDgError::Invalid("finite ordered work bracket, tolerances and step limits required"));
    }
    if let FiniteGasDgLoss::UniformViscothermal { maximum_relative_error } = controls.loss
        && !positive(maximum_relative_error)
    {
        return Err(FiniteGasDgError::Invalid("positive diffusion error admission required"));
    }
    let original = duct.gas_field().ok_or(FiniteGasDgError::Invalid(
        "finite gas coupling requires an original complete source field",
    ))?;
    if original.uniform_state() != Some(gas.state()) {
        return Err(FiniteGasDgError::Invalid(
            "finite gas and the uniform acoustic source must have one identical complete initial state",
        ));
    }
    let radius = duct.uniform_circular_radius_m().ok_or(FiniteGasDgError::Invalid(
        "finite uniform background requires the original circular fixed geometry",
    ))?;
    let volume = core::f64::consts::PI * radius * radius * duct.total_length();
    if !positive(volume) || (gas.volume_m3() / volume - 1.0).abs() > 64.0 * f64::EPSILON {
        return Err(FiniteGasDgError::Invalid("finite background volume must equal the actual duct volume"));
    }
    let mut marched = state.clone();
    // A foreign state/propagator refuses before running the physical callback.
    propagator.sync_nodal(&mut marched).map_err(native)?;
    duct.port_pressures(&marched).map_err(native)?;
    let mut acoustic_step = propagator.step(&mut marched, ports).map_err(native)?;
    let dissipation = marched.cumulative_dissipated_j() - state.cumulative_dissipated_j();
    if !dissipation.is_finite() {
        return Err(FiniteGasDgError::Invalid("native acoustic loss increment is unrepresentable"));
    }
    check(cancelled)?;
    let maximum_work = lo.abs().max(hi.abs());
    let candidate = feedback::solve_scalar(
        controls.work_bracket_j, controls.work_tolerance_j, controls.maximum_iterations,
        |work, trials, cancelled| {
            let background = evolve(gas, propagator.dt(), work, cancelled)
                .map_err(FiniteGasDgError::Background)?;
            check(cancelled)?;
            if background.gas.volume_m3() != gas.volume_m3()
                || !background.supplied_heat_j.is_finite()
                || !background.transported_enthalpy_j.is_finite()
                || (background.gas.state().state().temperature - gas.state().state().temperature).abs()
                    > controls.maximum_temperature_increment_k
            {
                return Err(FiniteGasDgError::Invalid(
                    "background changed fixed volume, exceeded temperature increment or returned nonfinite energy",
                ));
            }
            let gas_change = stored_change(gas, &background.gas);
            let gas_residual = ((gas_change + work) - background.supplied_heat_j)
                - background.transported_enthalpy_j;
            if !gas_change.is_finite() || !gas_residual.is_finite()
                || gas_residual.abs() > controls.energy_tolerance_j
            {
                return Err(FiniteGasDgError::Invalid("background endpoint does not close heat, enthalpy and trial work"));
            }
            let replacement = field(original, &background.gas)?;
            let refresh_kind = match controls.loss {
                FiniteGasDgLoss::Lossless => GasRefresh::Lossless,
                FiniteGasDgLoss::UniformViscothermal { maximum_relative_error } =>
                    GasRefresh::UniformDiffusion(maximum_relative_error),
            };
            let (next_duct, mut next_state, next_propagator, refresh) = prepare_gas(
                duct, &marched, propagator, &replacement, refresh_kind, maximum_work, cancelled,
            ).map_err(acoustic)?;
            next_propagator.sync_nodal(&mut next_state).map_err(native)?;
            let residual = work - refresh.parameter_work_j;
            // Each trial starts from the same marched wave and original gas.
            let step = FiniteGasDgStep {
                acoustic: acoustic_step.clone(), refresh, background: background.source,
                supplied_heat_j: background.supplied_heat_j,
                transported_enthalpy_j: background.transported_enthalpy_j,
                gas_energy_change_j: gas_change, work_by_gas_j: work,
                work_residual_j: residual, gas_energy_residual_j: gas_residual,
                acoustic_dissipation_j: dissipation, combined_energy_residual_j: 0.0, trials,
            };
            Ok((Candidate {
                duct: next_duct, state: next_state, propagator: next_propagator,
                gas: background.gas, step,
            }, residual))
        }, cancelled,
    ).map_err(|error| match error {
        feedback::ScalarError::Evaluation(error) => error,
        feedback::ScalarError::Invalid(message) => FiniteGasDgError::Invalid(message),
        feedback::ScalarError::Cancelled => FiniteGasDgError::Cancelled,
        feedback::ScalarError::OutsideBracket { lower, upper } =>
            FiniteGasDgError::OutsideBracket { lower_residual_j: lower, upper_residual_j: upper },
        feedback::ScalarError::NotConverged { iterations, residual } =>
            FiniteGasDgError::NotConverged { iterations, residual_j: residual },
    })?;
    let mut candidate = candidate;
    let residual = ((candidate.step.gas_energy_change_j
        + (candidate.state.stored_energy_j() - state.stored_energy_j())
        + dissipation) - acoustic_step.ledger.port_work_step_j
        - candidate.step.supplied_heat_j) - candidate.step.transported_enthalpy_j;
    if !residual.is_finite() || residual.abs() > controls.energy_tolerance_j {
        return Err(FiniteGasDgError::Invalid("independent gas plus wave/wall endpoint energy balance failed"));
    }
    // Impulses and work retain their actual acoustic substep; endpoint pressures
    // and the global ledger belong to the refreshed accepted state at the same clock.
    let pressures = candidate.duct.port_pressures(&candidate.state).map_err(native)?;
    for output in &mut acoustic_step.port_outputs {
        output.end_conjugate = pressures.iter().find(|(id, _)| *id == output.id)
            .ok_or(FiniteGasDgError::Invalid("refreshed native port map lost an original port"))?.1;
    }
    let ledger = &mut acoustic_step.ledger;
    ledger.stored_energy_j = candidate.state.stored_energy_j();
    ledger.cumulative_parameter_work_j = candidate.state.cumulative_parameter_work_j();
    ledger.energy_balance_residual_j = ledger.stored_energy_j - candidate.state.initial_energy_j()
        - ledger.cumulative_port_work_j - ledger.cumulative_parameter_work_j
        + ledger.cumulative_dissipated_j;
    let scale = (candidate.state.initial_energy_j() + ledger.cumulative_port_work_j.abs()
        + ledger.cumulative_parameter_work_j.abs()).max(1e-300);
    ledger.relative_energy_drift = ledger.energy_balance_residual_j.abs() / scale;
    candidate.step.acoustic = acoustic_step;
    candidate.step.combined_energy_residual_j = residual;
    check(cancelled)?;
    Ok(candidate)
}

fn heated<C: FnMut() -> bool>(
    old: &ClosedGasParcel, heat_j: f64, work: f64, solve: GasEnergySolve, cancelled: &mut C,
) -> Result<FiniteGasDgBackground<GasCompositionStep>, GasEnergyError> {
    let mut gas = old.clone();
    let source = gas.replace_model_with_heat_work_volume(
        old.model().clone(), old.volume_m3(), heat_j, work, solve, cancelled,
    )?;
    Ok(FiniteGasDgBackground { gas, supplied_heat_j: heat_j, transported_enthalpy_j: 0.0, source })
}

impl MaterialWallDg {
    /// Advance the native sample and a finite source gas atomically.
    ///
    /// The repeatable callback receives (original gas, native dt, work OUT of gas,
    /// cancellation). It MUST propose from that original snapshot and return its
    /// actual gas endpoint, heat and separately transported enthalpy. It must not
    /// mutate accepted companions: return their candidates in `source`, then
    /// publish them only after this call succeeds. Callback/source, field phase,
    /// diffusion, propagator, energy and cancellation refusals preserve all four
    /// accepted owners. The heat convenience method supplies the original caloric owner.
    ///
    /// Uniform well-mixed gas in the actual fixed circular volume only. This is
    /// a first-order partition: no mean-flow momentum, pressure prestress, changing
    /// geometry, wall thermal evolution or spontaneous sound from uniform heating.
    /// All native wall and diffusion histories survive. Native dissipative losses
    /// remain explicit exported energy; their unknown gas/wall split is not guessed.
    /// Timestep refinement and the callback's own source/transport gates remain required.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_controlled<S, E, C: FnMut() -> bool>(
        &mut self, state: &mut DgState, propagator: &mut DgPropagator,
        gas: &mut ClosedGasParcel, ports: &[PortInput], controls: FiniteGasDgSolve,
        evolve: impl FnMut(&ClosedGasParcel, f64, f64, &mut C)
            -> Result<FiniteGasDgBackground<S>, E>,
        cancelled: &mut C,
    ) -> Result<FiniteGasDgStep<S>, FiniteGasDgError<E>> {
        let candidate = prepare(&self.duct, state, propagator, gas, ports, controls, evolve, cancelled)?;
        self.duct = candidate.duct;
        *state = candidate.state;
        *propagator = candidate.propagator;
        *gas = candidate.gas;
        Ok(candidate.step)
    }

    /// Closed-inventory heating with the existing complete NASA-9/EOS energy owner.
    /// Supplied heat and reciprocal parameter work remain distinct thermodynamic terms.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_heat_controlled<C: FnMut() -> bool>(
        &mut self, state: &mut DgState, propagator: &mut DgPropagator,
        gas: &mut ClosedGasParcel, ports: &[PortInput], heat_j: f64,
        gas_solve: GasEnergySolve, controls: FiniteGasDgSolve, cancelled: &mut C,
    ) -> Result<FiniteGasDgStep<GasCompositionStep>, FiniteGasDgError<GasEnergyError>> {
        self.step_with_finite_gas_controlled(state, propagator, gas, ports, controls,
            |old, _, work, cancel| heated(old, heat_j, work, gas_solve, cancel), cancelled)
    }
}

impl RegionalMaterialWallDg {
    /// Same finite-gas transaction, preserving every original regional wall source
    /// and hereditary branch. See `MaterialWallDg::step_with_finite_gas_controlled`.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_controlled<S, E, C: FnMut() -> bool>(
        &mut self, state: &mut DgState, propagator: &mut DgPropagator,
        gas: &mut ClosedGasParcel, ports: &[PortInput], controls: FiniteGasDgSolve,
        evolve: impl FnMut(&ClosedGasParcel, f64, f64, &mut C)
            -> Result<FiniteGasDgBackground<S>, E>,
        cancelled: &mut C,
    ) -> Result<FiniteGasDgStep<S>, FiniteGasDgError<E>> {
        let candidate = prepare(&self.duct, state, propagator, gas, ports, controls, evolve, cancelled)?;
        self.duct = candidate.duct;
        *state = candidate.state;
        *propagator = candidate.propagator;
        *gas = candidate.gas;
        Ok(candidate.step)
    }

    /// Original caloric heating and reciprocal work on the same regional wall state.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_heat_controlled<C: FnMut() -> bool>(
        &mut self, state: &mut DgState, propagator: &mut DgPropagator,
        gas: &mut ClosedGasParcel, ports: &[PortInput], heat_j: f64,
        gas_solve: GasEnergySolve, controls: FiniteGasDgSolve, cancelled: &mut C,
    ) -> Result<FiniteGasDgStep<GasCompositionStep>, FiniteGasDgError<GasEnergyError>> {
        self.step_with_finite_gas_controlled(state, propagator, gas, ports, controls,
            |old, _, work, cancel| heated(old, heat_j, work, gas_solve, cancel), cancelled)
    }
}

