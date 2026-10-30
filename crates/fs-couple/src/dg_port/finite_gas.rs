//! Finite source-gas feedback on the original complete DG port load.
//!
//! One original acoustic step is followed by a scalar energy feedback solve.
//! Every trial starts from the same marched bore, branch, channel and terminal
//! histories and the same original thermochemical inventory. The original gas
//! owner pays the actual native bore AND finite-channel parameter work. No
//! source inversion, momentum transfer, acoustic integrator or root solver is
//! duplicated here.

/// Original multi-parcel transport through complete continuing acoustic loads.
pub mod transport;

use super::{DgGasRefresh, DgPortGasUpdate, DgPortLoad, DgTerminalChannelKind};
use crate::acoustic_realize::AcousticRealizeError;
pub use crate::bernoulli_aperture::wall::material::{
    FiniteGasDgBackground, FiniteGasDgError, FiniteGasDgLoss, FiniteGasDgSolve,
};
use crate::bernoulli_aperture::wall::material::{
    admit_finite_gas_controls, finite_gas_field, finite_gas_stored_change, heated_finite_gas,
};
use crate::cavity_phs::piston::thermal::feedback;
use fs_duct::dg::{DgBoundary, DgStepOutput};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::composition::GasCompositionStep;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};

/// Actual gas and represented acoustic endpoints of one accepted sample.
/// Rational exterior exchange remains signed boundary work, not invented storage.
#[derive(Debug)]
pub struct FiniteGasDgPortStep<S> {
    /// Original interval work and pressure impulses, with current endpoint
    /// pressure and complete bore/branch/channel ledger after the gas refresh.
    pub acoustic: DgStepOutput,
    /// Original native work reports for the bore and each finite channel.
    pub refresh: DgPortGasUpdate,
    /// Original background report and any still-private companion candidates.
    pub background: S,
    pub supplied_heat_j: f64,
    pub transported_enthalpy_j: f64,
    pub gas_energy_change_j: f64,
    /// Signed work OUT of the finite gas inventory and INTO represented motion.
    pub work_by_gas_j: f64,
    pub work_residual_j: f64,
    pub gas_energy_residual_j: f64,
    /// Original loss exported from this accounting boundary. Channel viscosity
    /// and numerical loss retain their distinct original reports; neither is
    /// guessed to be heat deposited in this well-mixed inventory.
    pub acoustic_dissipation_j: f64,
    /// Signed original rational terminal and optional driven-mouth exchange [J].
    /// Already excluded from net acoustic port work; never subtract it twice.
    pub terminal_work_j: f64,
    /// Gas plus represented storage change plus original loss, less actual net
    /// acoustic port work, supplied heat and transported enthalpy [J].
    pub combined_energy_residual_j: f64,
    pub trials: usize,
}

pub(crate) struct FiniteGasDgPortCandidate<S> {
    pub load: DgPortLoad,
    pub gas: ClosedGasParcel,
    pub step: FiniteGasDgPortStep<S>,
}

pub(crate) fn check<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), FiniteGasDgError<E>> {
    if cancelled() {
        Err(FiniteGasDgError::Cancelled)
    } else {
        Ok(())
    }
}

pub(crate) fn acoustic<E>(error: AcousticRealizeError) -> FiniteGasDgError<E> {
    match error {
        AcousticRealizeError::Cancelled => FiniteGasDgError::Cancelled,
        AcousticRealizeError::Duct(fs_duct::DuctError::Cancelled { .. }) => {
            FiniteGasDgError::Cancelled
        }
        error => FiniteGasDgError::Acoustic(error),
    }
}

fn native<E>(error: fs_duct::DuctError) -> FiniteGasDgError<E> {
    acoustic(AcousticRealizeError::Duct(error))
}

/// Admission shared with moving sources, before any mechanics or gas callback.
pub(crate) fn admit_initial<E>(
    load: &DgPortLoad,
    gas: &ClosedGasParcel,
    controls: FiniteGasDgSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), FiniteGasDgError<E>> {
    check(cancelled)?;
    admit_finite_gas_controls(controls)?;
    let volume = admitted_volume(load, gas, cancelled)?;
    if (gas.volume_m3() / volume - 1.0).abs() > 64.0 * f64::EPSILON {
        return Err(FiniteGasDgError::Invalid(
            "finite inventory volume must equal the main bore plus its finite terminal passages",
        ));
    }
    Ok(())
}

/// Admit one original source and return its actual main-gas geometric volume.
/// A deliberately shared parcel may span multiple complete loads; the caller
/// checks their total volume once instead of assigning its mass to each load.
pub(crate) fn admitted_volume<E>(
    load: &DgPortLoad,
    gas: &ClosedGasParcel,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, FiniteGasDgError<E>> {
    check(cancelled)?;
    let original = load.duct().gas_field().ok_or(FiniteGasDgError::Invalid(
        "finite gas port coupling requires the original complete source field",
    ))?;
    if original.uniform_state() != Some(gas.state()) {
        return Err(FiniteGasDgError::Invalid(
            "finite gas and the uniform main bore must share one complete initial source state",
        ));
    }
    let radius = load
        .duct()
        .uniform_circular_radius_m()
        .ok_or(FiniteGasDgError::Invalid(
            "finite gas port coupling requires the original uniform circular fixed main bore",
        ))?;
    let mut volume = core::f64::consts::PI * radius * radius * load.duct().total_length();
    for channel in load.terminal_channels() {
        check(cancelled)?;
        // A finite branch and its terminal passage have an independent source.
        if load
            .branches()
            .iter()
            .any(|branch| branch.port() == channel.input().port)
        {
            continue;
        }
        if channel.gas() != gas.state().state() {
            return Err(FiniteGasDgError::Invalid(
                "main-bore terminal channels must share the finite inventory's complete gas",
            ));
        }
        volume += match channel.input().kind {
            DgTerminalChannelKind::Closed => 0.0,
            DgTerminalChannelKind::Slit {
                width_m,
                height_m,
                length_m,
                ..
            } => width_m * height_m * length_m,
            DgTerminalChannelKind::Circular {
                radius_m, length_m, ..
            } => core::f64::consts::PI * radius_m * radius_m * length_m,
        };
    }
    if !volume.is_finite() || volume <= 0.0 {
        return Err(FiniteGasDgError::Invalid(
            "finite inventory geometric volume must be positive and representable",
        ));
    }
    load.port_response().map_err(acoustic)?;
    check(cancelled)?;
    Ok(volume)
}

pub(crate) fn refresh_kind(loss: FiniteGasDgLoss) -> DgGasRefresh {
    match loss {
        FiniteGasDgLoss::Lossless => DgGasRefresh::Lossless,
        FiniteGasDgLoss::UniformViscothermal {
            maximum_relative_error,
        } => DgGasRefresh::UniformDiffusion {
            maximum_relative_error,
        },
    }
}

/// The original finite gas must close its actual caloric/transport/work account
/// before any acoustic region consumes the proposed medium.
pub(crate) fn admit_background<S, E>(
    gas: &ClosedGasParcel,
    background: &FiniteGasDgBackground<S>,
    work: f64,
    controls: FiniteGasDgSolve,
) -> Result<(f64, f64), FiniteGasDgError<E>> {
    if background.gas.volume_m3() != gas.volume_m3()
        || !background.supplied_heat_j.is_finite()
        || !background.transported_enthalpy_j.is_finite()
        || (background.gas.state().state().temperature - gas.state().state().temperature).abs()
            > controls.maximum_temperature_increment_k
    {
        return Err(FiniteGasDgError::Invalid(
            "finite background changed volume or exceeded thermal/source limits",
        ));
    }
    let gas_change = finite_gas_stored_change(gas, &background.gas);
    let gas_residual =
        ((gas_change + work) - background.supplied_heat_j) - background.transported_enthalpy_j;
    if !gas_change.is_finite()
        || !gas_residual.is_finite()
        || gas_residual.abs() > controls.energy_tolerance_j
    {
        return Err(FiniteGasDgError::Invalid(
            "finite gas endpoint does not close its actual heat, enthalpy and work",
        ));
    }
    Ok((gas_change, gas_residual))
}

pub(crate) fn feedback_error<E>(
    error: feedback::ScalarError<FiniteGasDgError<E>>,
) -> FiniteGasDgError<E> {
    match error {
        feedback::ScalarError::Evaluation(error) => error,
        feedback::ScalarError::Invalid(message) => FiniteGasDgError::Invalid(message),
        feedback::ScalarError::Cancelled => FiniteGasDgError::Cancelled,
        feedback::ScalarError::OutsideBracket { lower, upper } => {
            FiniteGasDgError::OutsideBracket {
                lower_residual_j: lower,
                upper_residual_j: upper,
            }
        }
        feedback::ScalarError::NotConverged {
            iterations,
            residual,
        } => FiniteGasDgError::NotConverged {
            iterations,
            residual_j: residual,
        },
    }
}

/// Complete a gas feedback proposal after the original source has advanced its
/// private acoustic candidate exactly once. The original mechanics need not be
/// cloneable and are never repeated inside the scalar work solve.
#[allow(clippy::too_many_arguments)]
pub(crate) fn prepare_after_step<S, E, C: FnMut() -> bool>(
    original: &DgPortLoad,
    marched: &DgPortLoad,
    mut output: DgStepOutput,
    gas: &ClosedGasParcel,
    controls: FiniteGasDgSolve,
    mut evolve: impl FnMut(&ClosedGasParcel, f64, f64, &mut C) -> Result<FiniteGasDgBackground<S>, E>,
    cancelled: &mut C,
) -> Result<FiniteGasDgPortCandidate<S>, FiniteGasDgError<E>> {
    admit_initial(original, gas, controls, cancelled)?;
    if original.state().step_count().checked_add(1) != Some(marched.state().step_count())
        || output.ledger.step != marched.state().step_count()
        || original.duct().operator_identity() != marched.duct().operator_identity()
        || original.propagator().dt().to_bits() != marched.propagator().dt().to_bits()
        || original.duct().gas_field() != marched.duct().gas_field()
        || original.branches().len() != marched.branches().len()
        || original.terminal_channels().len() != marched.terminal_channels().len()
    {
        return Err(FiniteGasDgError::Invalid(
            "finite gas feedback requires one original-source acoustic step on the same load",
        ));
    }
    for (old, next) in original.branches().iter().zip(marched.branches()) {
        check(cancelled)?;
        if old.port() != next.port()
            || old.duct().operator_identity() != next.duct().operator_identity()
            || old.state().step_count().checked_add(1) != Some(next.state().step_count())
        {
            return Err(FiniteGasDgError::Invalid(
                "finite branch lost its original one-step history",
            ));
        }
    }
    for (old, next) in original
        .terminal_channels()
        .iter()
        .zip(marched.terminal_channels())
    {
        check(cancelled)?;
        if old.input() != next.input() || old.gas() != next.gas() {
            return Err(FiniteGasDgError::Invalid(
                "finite channel changed before gas feedback",
            ));
        }
    }
    let dissipation = marched.cumulative_dissipated_j() - original.cumulative_dissipated_j();
    let terminal_work_j = marched
        .last_terminal_frames()
        .iter()
        .fold(0.0, |sum, frame| sum + frame.load_work_j)
        + marched
            .last_driven_mouth_frame()
            .map_or(0.0, |frame| frame.load_work_j);
    if !dissipation.is_finite()
        || !terminal_work_j.is_finite()
        || !output.ledger.stored_energy_j.is_finite()
        || (output.ledger.stored_energy_j - marched.stored_energy_j()).abs()
            > controls.energy_tolerance_j
        || !output.ledger.port_work_step_j.is_finite()
    {
        return Err(FiniteGasDgError::Invalid(
            "finite gas feedback requires the actual complete acoustic storage and exchanges",
        ));
    }
    let original_field = original
        .duct()
        .gas_field()
        .ok_or(FiniteGasDgError::Invalid(
            "finite gas feedback lost the admitted original source field",
        ))?;
    let [lo, hi] = controls.work_bracket_j;
    let maximum_work = lo.abs().max(hi.abs());
    let refresh_kind = refresh_kind(controls.loss);
    let mut candidate = feedback::solve_scalar(
        controls.work_bracket_j,
        controls.work_tolerance_j,
        controls.maximum_iterations,
        |work, trials, cancelled| {
            let background = evolve(gas, original.propagator().dt(), work, cancelled)
                .map_err(FiniteGasDgError::Background)?;
            check(cancelled)?;
            let (gas_change, gas_residual) = admit_background(gas, &background, work, controls)?;
            let replacement = finite_gas_field(original_field, &background.gas)?;
            let mut load = marched.clone();
            let refresh = load
                .rebind_gas_field_controlled(&replacement, refresh_kind, maximum_work, cancelled)
                .map_err(acoustic)?;
            let residual = work - refresh.parameter_work_j;
            Ok((
                FiniteGasDgPortCandidate {
                    load,
                    gas: background.gas,
                    step: FiniteGasDgPortStep {
                        acoustic: output.clone(),
                        refresh,
                        background: background.source,
                        supplied_heat_j: background.supplied_heat_j,
                        transported_enthalpy_j: background.transported_enthalpy_j,
                        gas_energy_change_j: gas_change,
                        work_by_gas_j: work,
                        work_residual_j: residual,
                        gas_energy_residual_j: gas_residual,
                        acoustic_dissipation_j: dissipation,
                        terminal_work_j,
                        combined_energy_residual_j: 0.0,
                        trials,
                    },
                },
                residual,
            ))
        },
        cancelled,
    )
    .map_err(feedback_error)?;

    // Check extensive endpoints independently of the scalar residual. The
    // acoustic row's net port work already excludes the rational exterior.
    let residual = ((candidate.step.gas_energy_change_j
        + (candidate.load.stored_energy_j() - original.stored_energy_j())
        + dissipation)
        - output.ledger.port_work_step_j
        - candidate.step.supplied_heat_j)
        - candidate.step.transported_enthalpy_j;
    if !residual.is_finite() || residual.abs() > controls.energy_tolerance_j {
        return Err(FiniteGasDgError::Invalid(
            "finite gas plus complete acoustic storage, loss and boundary work do not close",
        ));
    }
    refresh_output(&mut candidate.load, &mut output, cancelled)?;
    candidate.step.acoustic = output;
    candidate.step.combined_energy_residual_j = residual;
    check(cancelled)?;
    Ok(candidate)
}

/// Keep the original interval integrals while refreshing endpoint pressure and
/// the complete current ledger after a zero-duration native gas intervention.
pub(crate) fn refresh_output<E>(
    load: &mut DgPortLoad,
    output: &mut DgStepOutput,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), FiniteGasDgError<E>> {
    // Interval observations remain from the actual frozen-medium step. Only
    // endpoint pressure and current extensive storage follow the accepted gas.
    load.propagator
        .sync_nodal(&mut load.state)
        .map_err(native)?;
    let main_pressures = load.duct().port_pressures(load.state()).map_err(native)?;
    for port in &mut output.port_outputs {
        check(cancelled)?;
        port.end_conjugate = if let Some(branch) = load
            .branches()
            .iter()
            .find(|branch| branch.port() == port.id)
        {
            let DgBoundary::Port(outlet) = branch.duct().outlet() else {
                return Err(FiniteGasDgError::Invalid(
                    "finite branch lost its original outlet",
                ));
            };
            branch
                .duct()
                .port_pressures(branch.state())
                .map_err(native)?
                .into_iter()
                .find(|(id, _)| *id == outlet)
                .map(|(_, pressure)| pressure)
        } else {
            main_pressures
                .iter()
                .find(|(id, _)| *id == port.id)
                .map(|(_, pressure)| *pressure)
        }
        .ok_or(FiniteGasDgError::Invalid(
            "finite gas refresh lost an original external port",
        ))?;
    }
    let initial_energy = load
        .branches()
        .iter()
        .fold(load.state().initial_energy_j(), |sum, branch| {
            sum + branch.state().initial_energy_j()
        });
    let ledger = &mut output.ledger;
    ledger.stored_energy_j = load.stored_energy_j();
    ledger.cumulative_parameter_work_j = load.cumulative_parameter_work_j();
    ledger.energy_balance_residual_j = ledger.stored_energy_j
        - initial_energy
        - ledger.cumulative_port_work_j
        - ledger.cumulative_parameter_work_j
        + ledger.cumulative_dissipated_j;
    ledger.relative_energy_drift = ledger.energy_balance_residual_j.abs()
        / (initial_energy
            + ledger.cumulative_port_work_j.abs()
            + ledger.cumulative_parameter_work_j.abs())
        .max(1e-300);
    if !ledger.energy_balance_residual_j.is_finite() || !ledger.relative_energy_drift.is_finite() {
        return Err(FiniteGasDgError::Invalid(
            "refreshed complete acoustic ledger is not finite",
        ));
    }
    check(cancelled)
}

impl DgPortLoad {
    /// Advance one original held-flow sample and its finite thermochemical
    /// inventory atomically, including actual finite-terminal momentum work.
    ///
    /// The callback proposes from the ORIGINAL gas on every trial and receives
    /// (original gas, native dt, work OUT of gas, cancellation). Return actual
    /// external heat, net transported enthalpy and private companion endpoints;
    /// publish companions only after this method succeeds. The existing gas
    /// owner must close that heat/enthalpy/work account before acoustic refresh.
    ///
    /// The well-mixed fixed volume is the original uniform circular main bore
    /// PLUS every finite slit/circular terminal passage using its gas. Finite
    /// branches and their passages keep their independent gas and are excluded
    /// from this inventory. Existing equal-pressure junction/exterior admission
    /// still applies. The bracket magnitude also bounds the sum of absolute
    /// native parameter-work contributions, so opposing transfers cannot hide
    /// an excessive intervention.
    ///
    /// This is a first-order frozen-initial-medium partition. It does not infer
    /// mean mass transport from acoustic flow, thermalize numerical loss, solve
    /// static pressure prestress or move geometry. The callback supplies any
    /// actual mass/species transport. Original losses leave this accounting
    /// boundary with their original classifications. Refine timestep and source
    /// controls; numerical closure is not measured acoustic validation.
    pub fn step_with_finite_gas_controlled<S, E, C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        inlet_flow_m3_s: f64,
        controls: FiniteGasDgSolve,
        evolve: impl FnMut(&ClosedGasParcel, f64, f64, &mut C) -> Result<FiniteGasDgBackground<S>, E>,
        cancelled: &mut C,
    ) -> Result<FiniteGasDgPortStep<S>, FiniteGasDgError<E>> {
        admit_initial(self, gas, controls, cancelled)?;
        let mut marched = self.clone();
        let output = marched
            .preview_step(inlet_flow_m3_s)
            .map_err(acoustic)?
            .accept();
        let candidate =
            prepare_after_step(self, &marched, output, gas, controls, evolve, cancelled)?;
        *self = candidate.load;
        *gas = candidate.gas;
        Ok(candidate.step)
    }

    /// Original closed-inventory caloric heating and reciprocal acoustic work.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_heat_controlled<C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        inlet_flow_m3_s: f64,
        heat_j: f64,
        gas_solve: GasEnergySolve,
        controls: FiniteGasDgSolve,
        cancelled: &mut C,
    ) -> Result<FiniteGasDgPortStep<GasCompositionStep>, FiniteGasDgError<GasEnergyError>> {
        self.step_with_finite_gas_controlled(
            gas,
            inlet_flow_m3_s,
            controls,
            |old, _, work, cancel| heated_finite_gas(old, heat_j, work, gas_solve, cancel),
            cancelled,
        )
    }
}
