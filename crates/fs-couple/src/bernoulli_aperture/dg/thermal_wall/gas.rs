//! Reciprocal heat between the continuing DG gas and its finite material wall.
//!
//! The accepted mechanics and complete acoustic load march once with their
//! initial material. Every interface trial then starts from that same marched
//! state and the same original inventories. Gas work is resolved first and
//! wall work second; these zero-duration storage interventions telescope.
//! The shared conductance law uses actual internal-energy/enthalpy midpoints.

use super::{
    ApertureDg, DgApertureFrame, FiniteWallThermalSolve, accept_acoustic_candidate_observed, admit,
    admit_phase, body_candidate, interval_heat, invalid, phase_error, poll,
};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dynamic::ApertureDrive;
use crate::bernoulli_aperture::modal::pressure::MaterialPressureFrame;
use crate::bernoulli_aperture::tube::TubeDrive;
use crate::bernoulli_aperture::two_sided::LoadedAperture;
use crate::bernoulli_aperture::wall::material::heated_finite_gas;
use crate::cavity_phs::piston::thermal::layers::MaterialConductionLayers;
use crate::cavity_phs::piston::thermal::{ConductanceSample, HeatExchangeSolve, feedback};
use crate::dg_port::finite_gas::{
    FiniteGasDgError, FiniteGasDgPortStep, FiniteGasDgSolve, admit_initial as admit_gas,
    prepare_after_step, refresh_output,
};
use crate::dg_port::{DgPortLoad, DgPortWallUpdate, sum_terminal_work};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStepInput};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::composition::GasCompositionStep;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};

/// Separate native work, caloric, interface and total nested-work allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGasWallThermalSolve {
    pub gas: FiniteGasDgSolve,
    pub caloric: GasEnergySolve,
    pub wall: FiniteWallThermalSolve,
    /// Signed interface heat INTO gas and OUT of the finite body [J].
    pub exchange: HeatExchangeSolve,
    /// Across every exchange trial, cap the number of gas endpoint solves,
    /// wall endpoint solves, gas midpoint reconstructions and interface queries.
    /// Their internal inversion/source loops retain their own explicit budgets.
    /// Cold source/geometry admission and cheap body chart queries are excluded.
    pub maximum_material_trials: usize,
}

/// One accepted gas/body/mechanical/acoustic interval.
#[derive(Debug)]
pub struct FiniteGasWallApertureDgStep {
    /// Final gas AND wall endpoint; original interval work/impulses remain intact.
    pub coupled: DgApertureFrame,
    /// Gas-first intervention, before wall refresh. This retained intermediate
    /// report closes its own native gas/work account. Its acoustic observation
    /// and coupled.duct are two endpoints of ONE march, never additive storage.
    pub gas: FiniteGasDgPortStep<GasCompositionStep>,
    pub body_state: EquilibriumPhaseState,
    pub body_midpoint: EquilibriumPhaseState,
    pub gas_midpoint: ResolvedThermochemicalGasState,
    pub interface: ConductanceSample<ResolvedMaterialStatePoint>,
    /// Signed heat from body to gas [J].
    pub exchanged_heat_j: f64,
    pub wall: DgPortWallUpdate,
    /// Original main-wall Kelvin/Maxwell heat, deposited exactly once [J].
    pub wall_heat_j: f64,
    pub external_gas_heat_j: f64,
    pub external_body_heat_j: f64,
    pub body_storage_change_j: f64,
    /// Body storage + actual wall work + interface heat - its two inputs [J].
    pub body_energy_residual_j: f64,
    pub wall_work_residual_j: f64,
    /// Q - dt G(T_body_mid - T_gas_mid) [J].
    pub exchange_residual_j: f64,
    /// Original total mechanical/acoustic loss less deposited wall heat [J].
    pub exported_dissipation_j: f64,
    /// Both finite inventories and complete motion, with all external work [J].
    pub combined_energy_residual_j: f64,
    pub exchange_trials: usize,
    /// Gas endpoint proposals summed across ALL exchange trials.
    pub gas_trials: usize,
    /// Wall endpoint proposals summed across ALL exchange trials.
    pub body_trials: usize,
    pub material_trials: usize,
}

fn caloric_error(error: GasEnergyError) -> AcousticRealizeError {
    match error {
        GasEnergyError::Cancelled => AcousticRealizeError::Cancelled,
        error => AcousticRealizeError::Nonlinear(error.to_string()),
    }
}

fn gas_error(error: FiniteGasDgError<AcousticRealizeError>) -> AcousticRealizeError {
    match error {
        FiniteGasDgError::Cancelled => AcousticRealizeError::Cancelled,
        FiniteGasDgError::Acoustic(error) | FiniteGasDgError::Background(error) => error,
        error => AcousticRealizeError::Nonlinear(format!("finite gas/wall feedback: {error:?}")),
    }
}

fn exchange_error(error: feedback::ScalarError<AcousticRealizeError>) -> AcousticRealizeError {
    match error {
        feedback::ScalarError::Evaluation(error) => error,
        feedback::ScalarError::Cancelled => AcousticRealizeError::Cancelled,
        feedback::ScalarError::Invalid(what) => invalid(what),
        feedback::ScalarError::OutsideBracket { lower, upper } => AcousticRealizeError::Nonlinear(
            format!("gas/wall interface heat is outside its bracket: {lower} J, {upper} J"),
        ),
        feedback::ScalarError::NotConverged {
            iterations,
            residual,
        } => AcousticRealizeError::Nonlinear(format!(
            "gas/wall interface did not converge after {iterations} trials: {residual} J"
        )),
    }
}

#[derive(Default)]
struct Trials {
    gas: usize,
    body: usize,
    total: usize,
}

impl Trials {
    fn charge(&mut self, limit: usize) -> Result<(), AcousticRealizeError> {
        if self.total >= limit {
            return Err(invalid("joint gas/wall material-trial budget exhausted"));
        }
        self.total += 1;
        Ok(())
    }
}

fn admit_joint<A>(
    system: &ApertureDg<A>,
    gas: &ClosedGasParcel,
    body: UniformEnthalpyStepInput<'_>,
    layers: &MaterialConductionLayers<'_>,
    controls: FiniteGasWallThermalSolve,
    cold: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AcousticRealizeError> {
    poll(cancelled)?;
    admit(
        system,
        body,
        layers.material,
        layers.temperature_axis,
        controls.wall,
        cold,
        cancelled,
    )?;
    controls.caloric.validate().map_err(caloric_error)?;
    feedback::admit::<core::convert::Infallible>(controls.exchange).map_err(|_| {
        invalid("joint gas/wall heat requires a valid signed interface bracket and energy limits")
    })?;
    let temperature = gas.state().state().temperature;
    if controls.maximum_material_trials == 0
        || temperature < controls.caloric.minimum_temperature_k
        || temperature > controls.caloric.maximum_temperature_k
    {
        return Err(invalid(
            "joint gas/wall heat requires a material budget and initial caloric support",
        ));
    }
    if system.load.radiation_record().is_some() || !system.load.branches().is_empty() {
        return Err(invalid(
            "closed gas/wall heating cannot impose a mean-pressure jump on a computed exterior or independent branch; transport is required",
        ));
    }
    if system
        .load
        .duct()
        .gas_field()
        .is_none_or(|field| field.phase() != ThermochemicalPhase::FrozenGas)
    {
        return Err(invalid(
            "joint gas/wall heat requires the original explicit frozen-gas field",
        ));
    }
    admit_gas::<AcousticRealizeError>(&system.load, gas, controls.gas, cancelled)
        .map_err(gas_error)?;
    // Both source conductivities are required even for declared insulation.
    layers
        .resolve(gas.state(), body.initial)
        .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
    if cold {
        let mut load = system.load.clone();
        load.rebind_gas_field_controlled(
            system
                .load
                .duct()
                .gas_field()
                .expect("admitted original gas field"),
            controls.gas.loss.into(),
            controls.gas.work_bracket_j[0]
                .abs()
                .max(controls.gas.work_bracket_j[1].abs()),
            cancelled,
        )?;
    }
    poll(cancelled)
}

/// Cold admission without advancing or publishing any physical participant.
pub(crate) fn admit_initial(
    system: &ApertureDg,
    gas: &ClosedGasParcel,
    body: UniformEnthalpyStepInput<'_>,
    layers: &MaterialConductionLayers<'_>,
    controls: FiniteGasWallThermalSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AcousticRealizeError> {
    admit_joint(system, gas, body, layers, controls, true, cancelled)
}

struct Candidate {
    load: DgPortLoad,
    gas: ClosedGasParcel,
    step: FiniteGasWallApertureDgStep,
}

#[allow(clippy::too_many_arguments)]
fn finish(
    original: &DgPortLoad,
    marched: &DgPortLoad,
    initial_frame: DgApertureFrame,
    initial_energy_j: f64,
    gas: &ClosedGasParcel,
    body: UniformEnthalpyStepInput<'_>,
    layers: &MaterialConductionLayers<'_>,
    external_gas_heat_j: f64,
    wall_heat_j: f64,
    controls: FiniteGasWallThermalSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Candidate, AcousticRealizeError> {
    if !external_gas_heat_j.is_finite() {
        return Err(invalid("external gas heat must be finite"));
    }
    let mut counts = Trials::default();
    let mut accepted = feedback::solve_scalar(
        controls.exchange.heat_bracket_j,
        controls.exchange.source_tolerance_j,
        controls.exchange.maximum_iterations,
        |q, exchange_trials, cancelled| {
            let gas_heat = external_gas_heat_j + q;
            let body_heat = body.internal_heat_j - q;
            if !gas_heat.is_finite() || !body_heat.is_finite()
                || (q != 0.0 && (gas_heat == external_gas_heat_j || body_heat == body.internal_heat_j))
                || (external_gas_heat_j != 0.0 && gas_heat == q)
                || (body.internal_heat_j != 0.0 && body_heat == -q)
            {
                return Err(invalid("interface or external heat disappears in its finite inventory"));
            }
            // Both nested solves return to their original conserved inventories.
            // Only the zero-duration native interventions are ordered gas, wall.
            let gas_candidate = prepare_after_step(
                original, marched, initial_frame.duct.clone(), gas, controls.gas,
                |old, _, work, cancel| {
                    counts.charge(controls.maximum_material_trials)?;
                    counts.gas += 1;
                    heated_finite_gas(old, gas_heat, work, controls.caloric, cancel)
                        .map_err(caloric_error)
                },
                cancelled,
            ).map_err(gas_error)?;
            let mut wall_candidate = body_candidate(
                &gas_candidate.load,
                UniformEnthalpyStepInput { internal_heat_j: body_heat, ..body },
                layers.material, layers.temperature_axis, controls.wall, wall_heat_j,
                &mut || {
                    counts.charge(controls.maximum_material_trials)?;
                    counts.body += 1;
                    Ok(())
                },
                cancelled,
            )?;

            counts.charge(controls.maximum_material_trials)?;
            let midpoint_heat = 0.5 * gas_candidate.step.gas_energy_change_j;
            if gas_candidate.step.gas_energy_change_j != 0.0 && midpoint_heat == 0.0 {
                return Err(invalid("gas energy midpoint is unrepresentable"));
            }
            let midpoint = heated_finite_gas(
                gas, midpoint_heat, 0.0, controls.caloric, cancelled,
            ).map_err(caloric_error)?;
            let gas_midpoint = midpoint.gas.state().clone();
            let h0 = body.initial.specific_enthalpy_j_kg();
            let dh = wall_candidate.state.specific_enthalpy_j_kg() - h0;
            let hm = h0 + 0.5 * dh;
            if !hm.is_finite() || (dh != 0.0 && hm == h0) {
                return Err(invalid("body enthalpy midpoint is unrepresentable"));
            }
            let body_midpoint = body.curve.state_at_specific_enthalpy(hm).map_err(phase_error)?;
            admit_phase(body_midpoint, body.mass_kg / body.volume_m3)?;
            poll(cancelled)?;
            counts.charge(controls.maximum_material_trials)?;
            let interface = layers.resolve(&gas_midpoint, body_midpoint)
                .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
            let conductance = interface.conductance_w_k;
            let delta = body_midpoint.temperature_k() - gas_midpoint.state().temperature;
            let power = conductance * delta;
            let interface_heat = body.duration_s * power;
            let residual = q - interface_heat;
            if ![power, interface_heat, residual].iter().all(|x| x.is_finite())
                || (conductance != 0.0 && delta != 0.0 && power == 0.0)
                || (power != 0.0 && interface_heat == 0.0)
            {
                return Err(invalid("midpoint interface heat must be finite and representable"));
            }

            let mut frame = initial_frame.clone();
            refresh_output::<AcousticRealizeError>(
                &mut wall_candidate.load, &mut frame.duct, cancelled,
            ).map_err(gas_error)?;
            let gas_work = gas_candidate.step.refresh.parameter_work_j;
            let wall_work = wall_candidate.wall.parameter_work_j;
            frame.parameter_work_j = gas_work + wall_work;
            frame.stored_energy_j =
                frame.aperture.stored_energy_j + wall_candidate.load.stored_energy_j();
            frame.storage_change_j = frame.stored_energy_j - initial_energy_j;
            let intervention_residual = (wall_candidate.load.stored_energy_j()
                - marched.stored_energy_j()) - frame.parameter_work_j;
            let exported = frame.inferred_total_dissipation_j - wall_heat_j;
            let body_residual = ((wall_candidate.storage_change_j + wall_work + q)
                - body.internal_heat_j) - wall_heat_j;
            let combined = ((gas_candidate.step.gas_energy_change_j
                + wall_candidate.storage_change_j + frame.storage_change_j
                + exported + sum_terminal_work(&frame.terminals))
                - frame.upstream_work_j - frame.body_work_j
                - frame.aperture.mechanical_work_j - external_gas_heat_j)
                - body.internal_heat_j;
            let acoustic_residual = frame.balance_residual_j();
            let tolerance = controls.exchange.energy_tolerance_j;
            if ![frame.parameter_work_j, frame.stored_energy_j, frame.storage_change_j,
                intervention_residual, exported, body_residual, combined, acoustic_residual]
                .iter().all(|x| x.is_finite())
                || (wall_heat_j != 0.0 && exported == frame.inferred_total_dissipation_j)
                || (wall_heat_j > 0.0 && exported < -tolerance)
                || intervention_residual.abs() > tolerance
                || wall_candidate.body_energy_residual_j.abs()
                    > controls.wall.exchange.energy_tolerance_j
                || body_residual.abs() > tolerance
                || acoustic_residual.abs() > tolerance
                || combined.abs() > tolerance
            {
                return Err(invalid(
                    "finite gas, body, native interventions and coupled motion do not close their energy balances",
                ));
            }
            poll(cancelled)?;
            Ok((Candidate {
                load: wall_candidate.load,
                gas: gas_candidate.gas,
                step: FiniteGasWallApertureDgStep {
                    coupled: frame,
                    gas: gas_candidate.step,
                    body_state: wall_candidate.state,
                    body_midpoint,
                    gas_midpoint,
                    interface,
                    exchanged_heat_j: q,
                    wall: wall_candidate.wall,
                    wall_heat_j,
                    external_gas_heat_j,
                    external_body_heat_j: body.internal_heat_j,
                    body_storage_change_j: wall_candidate.storage_change_j,
                    body_energy_residual_j: body_residual,
                    wall_work_residual_j: wall_candidate.source_residual_j,
                    exchange_residual_j: residual,
                    exported_dissipation_j: exported,
                    combined_energy_residual_j: combined,
                    exchange_trials,
                    gas_trials: 0,
                    body_trials: 0,
                    material_trials: 0,
                },
            }, residual))
        },
        cancelled,
    ).map_err(exchange_error)?;
    accepted.step.gas_trials = counts.gas;
    accepted.step.body_trials = counts.body;
    accepted.step.material_trials = counts.total;
    poll(cancelled)?;
    Ok(accepted)
}

impl ApertureDg {
    /// Advance the original scalar source, complete DG load, closed finite gas
    /// and original finite main-wall body as one transaction.
    ///
    /// One interface heat is credited to gas and debited from body. Native
    /// gas and wall parameter work are paid by their respective caloric owners.
    /// Optional native solid heat is measured once from actual solved ports
    /// and removed once from the enclosing exported loss. Numerical, mechanical,
    /// gas and terminal loss are not relabeled as material heating.
    ///
    /// Fixed volume, mass/composition and a complete all-solid fixed-density
    /// body are required. This is the existing first-order frozen acoustic
    /// partition with midpoint interface conduction, not a resolved thermal
    /// profile, advection, convection, expansion, phase or mean-pressure flow
    /// model. Independent branches, upstream inventories and computed exterior
    /// radiation need their own transport closure and are explicitly refused.
    /// Refusal/cancellation leaves all inventories, material and acoustic
    /// histories, mechanics, terminal memory and clocks unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_wall_heat_controlled(
        &mut self,
        gas: &mut ClosedGasParcel,
        body_state: &mut EquilibriumPhaseState,
        body: UniformEnthalpyStepInput<'_>,
        layers: &MaterialConductionLayers<'_>,
        drive: TubeDrive,
        force_n: f64,
        external_gas_heat_j: f64,
        controls: FiniteGasWallThermalSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FiniteGasWallApertureDgStep, AcousticRealizeError> {
        if *body_state != body.initial {
            return Err(invalid(
                "joint wall input must begin at the accepted body state",
            ));
        }
        admit_joint(self, gas, body, layers, controls, false, cancelled)?;
        let initial_energy_j = self.stored_energy_j();
        let relation = self.load.port_response()?;
        let mechanical = self.aperture.preview_step_with_force(
            ApertureDrive {
                upstream_pressure_pa: drive.upstream_pressure_pa,
                incoming_pressure_pa: relation.history_pressure_pa / 2.0,
                body_flow_m3_s: drive.body_flow_m3_s,
            },
            force_n,
        )?;
        let mut marched = self.load.clone();
        let (frame, heat) = accept_acoustic_candidate_observed(
            &mut marched,
            relation,
            mechanical.frame,
            [MaterialPressureFrame::default(); 2],
            drive,
            |trial| interval_heat(trial, controls.wall, cancelled),
        )?;
        let candidate = finish(
            &self.load,
            &marched,
            frame,
            initial_energy_j,
            gas,
            body,
            layers,
            external_gas_heat_j,
            heat,
            controls,
            cancelled,
        )?;
        // Source rebind validates before its first mutation; publication below
        // has no fallible operation or cancellation poll.
        self.aperture.rebind_exact_gas_load(
            candidate.load.duct().inlet_gas(),
            candidate.load.impedance_pa_s_m3(),
        )?;
        self.aperture.accept_frame(mechanical);
        self.load = candidate.load;
        *gas = candidate.gas;
        *body_state = candidate.step.body_state;
        Ok(candidate.step)
    }
}

impl ApertureDg<LoadedAperture> {
    /// The same joint transaction for original single, paired or modal mechanics.
    /// Independent force/pressure work remains in its existing physical owner.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_wall_heat_controlled(
        &mut self,
        gas: &mut ClosedGasParcel,
        body_state: &mut EquilibriumPhaseState,
        body: UniformEnthalpyStepInput<'_>,
        layers: &MaterialConductionLayers<'_>,
        drive: TubeDrive,
        forces_n: [f64; 2],
        pressures_pa: [f64; 2],
        external_gas_heat_j: f64,
        controls: FiniteGasWallThermalSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FiniteGasWallApertureDgStep, AcousticRealizeError> {
        if *body_state != body.initial {
            return Err(invalid(
                "joint wall input must begin at the accepted body state",
            ));
        }
        admit_joint(self, gas, body, layers, controls, false, cancelled)?;
        let initial_energy_j = self.stored_energy_j();
        let relation = self.load.port_response()?;
        let mechanical = self.aperture.preview(
            ApertureDrive {
                upstream_pressure_pa: drive.upstream_pressure_pa,
                incoming_pressure_pa: relation.history_pressure_pa / 2.0,
                body_flow_m3_s: drive.body_flow_m3_s,
            },
            forces_n,
            pressures_pa,
        )?;
        let mut marched = self.load.clone();
        let (frame, heat) = accept_acoustic_candidate_observed(
            &mut marched,
            relation,
            mechanical.frame(),
            mechanical.spatial_pressure(),
            drive,
            |trial| interval_heat(trial, controls.wall, cancelled),
        )?;
        let candidate = finish(
            &self.load,
            &marched,
            frame,
            initial_energy_j,
            gas,
            body,
            layers,
            external_gas_heat_j,
            heat,
            controls,
            cancelled,
        )?;
        self.aperture.rebind_exact_gas_load(
            candidate.load.duct().inlet_gas(),
            candidate.load.impedance_pa_s_m3(),
        )?;
        self.aperture.accept(mechanical);
        self.load = candidate.load;
        *gas = candidate.gas;
        *body_state = candidate.step.body_state;
        Ok(candidate.step)
    }
}
