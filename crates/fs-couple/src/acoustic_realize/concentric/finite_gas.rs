//! Reciprocal finite-gas heat on the native layered-string transaction.
//!
//! The existing layered thermal network retains each annulus's caloric and
//! thermoelastic storage. The source gas parcel retains its composition,
//! internal-energy chart, EOS and transport. Aerodynamic work heats that parcel
//! once; a signed outer-surface transfer is debited from the solid and credited
//! to the gas. No new heat-capacity law, drag law or heat solver is introduced.

/// Conservative atmospheric species exchange on the continuing layered string.
pub mod flow;

use std::cell::{Cell, RefCell};

use super::{
    ConcentricStringError, ConcentricStringOhmicSource, ConcentricStringThermalFrame,
    ConcentricStringThermalSolve, ThermalConcentricStringRuntime, add_ohmic_heat,
    checkpoint, deposit, finite, invalid, positive,
};
use crate::acoustic_realize::{
    AcousticRealizeError, finite_gas::FiniteGasStringSolve,
    linear_accelerations, linear_string_member,
};
use crate::cavity_phs::piston::thermal::{
    ConductanceSample, HeatExchangeSolve, feedback,
    layers::{LayerConductionError, MaterialConductionLayers},
};
use fs_exec::Cx;
use fs_material::conductor::OhmicDrive;
use fs_material::phase::EquilibriumPhaseState;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::{
    ResolvedThermochemicalGasState,
    energy::{GasHeatStep, IsochoricGasParcel},
};
use fs_thermal::enthalpy::{
    DepositionError, DepositionSimulation, backward_euler::StepCoefficients,
};

/// Explicit budgets of the existing layer, source-gas and interface owners.
#[derive(Clone, Copy, Debug)]
pub struct ConcentricFiniteGasSolve {
    pub thermal: ConcentricStringThermalSolve,
    pub gas: FiniteGasStringSolve,
    pub exchange: HeatExchangeSolve,
}

/// One patch on the exposed lateral surface of the outermost annulus.
/// Internal annular contacts remain links in the original thermal network.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConcentricGasContact {
    pub area_m2: f64,
}

/// Actual accepted material/heat sources, preserving their concrete receipt types.
#[derive(Debug)]
pub struct ConcentricFiniteGasSource<S, R> {
    pub thermal: S,
    pub exchange: R,
    pub electrical: Option<ConcentricStringOhmicSource<()>>,
    pub gas: GasHeatStep,
    /// Accepted heat-exchange endpoint. In the `flow` composition this precedes
    /// species/enthalpy transport; it remains the source of the thermal receipts.
    pub next_gas: IsochoricGasParcel,
}

/// One accepted sample of the layered solid, motion, observer and source gas.
#[derive(Debug)]
pub struct ConcentricFiniteGasFrame<S, R> {
    /// Its solid subaccount includes -Q in external_heat_j, and treats air loss
    /// as an outgoing transfer. Both are internal in the combined account below.
    pub string: ConcentricStringThermalFrame<ConcentricFiniteGasSource<S, R>>,
    pub solid_to_gas_heat_j: f64,
    pub drag_heat_j: f64,
    pub external_gas_heat_j: f64,
    /// Independent thermal boundary/callback and optional electrical supply.
    /// Excludes aerodynamic/material damping and reversible/internal transfers.
    pub external_solid_heat_j: f64,
    pub conductance_w_k: f64,
    /// Q - dt G(T_outer,end - T_gas,end) [J].
    pub source_residual_j: f64,
    /// Work + independent supplies - changes of complete solid/mechanical/gas
    /// storage - authored loss - damping-partition roundoff residual [J].
    pub energy_balance_residual_j: f64,
    pub energy_tolerance_j: f64,
    /// Full coupled proposals; at most exchange.maximum_iterations + 2.
    pub trials: usize,
    /// Sum of actual thermal residual evaluations across ALL coupled proposals.
    pub thermal_trials: usize,
}

struct Candidate<'a, 't, S, R> {
    runtime: ThermalConcentricStringRuntime<'a>,
    thermal: DepositionSimulation<'t>,
    string: ConcentricStringThermalFrame<ConcentricFiniteGasSource<S, R>>,
    heat_j: f64,
    conductance_w_k: f64,
    external_solid_heat_j: f64,
    source_residual_j: f64,
    trials: usize,
}

// Shared by heat-only and conservative species-flow compositions. Keep the
// original outer material pin and lumped-domain check in one physical owner.
fn finite_gas_conductance<E: From<LayerConductionError>>(
    layers: &MaterialConductionLayers<'_>, binding: (fs_matdb::ClaimId, f64),
    gas: &ResolvedThermochemicalGasState, solid: EquilibriumPhaseState,
) -> Result<ConductanceSample<ResolvedMaterialStatePoint>, E> {
    let (selected, volume) = binding;
    let sample = layers.resolve(gas, solid).map_err(E::from)?;
    let property = sample.source.property_by_key(layers.conductivity)
        .ok_or_else(|| E::from(LayerConductionError::Invalid("resolved outer conductivity is absent")))?;
    if property.answer().receipt.selected != selected {
        return Err(E::from(LayerConductionError::Invalid("surface conductivity changed its original layer source pin")));
    }
    let biot = if sample.conductance_w_k == 0.0 { 0.0 }
        else { sample.conductance_w_k*(volume/layers.area_m2)
            /layers.area_m2/property.value_si() };
    if !biot.is_finite() || biot > fs_conduction::lumped::LUMPED_BIOT_CEILING {
        return Err(E::from(LayerConductionError::Invalid(
            "outer annulus exceeds the lumped Biot limit; use spatial thermal transport")));
    }
    Ok(sample)
}

impl<'a> ThermalConcentricStringRuntime<'a> {
    /// Continue unequal-temperature layers, optional DC drive and a finite
    /// source atmosphere in ONE transaction, preserving modal/observer history.
    ///
    /// Mechanics freezes accepted solid/gas coefficients for one exact sample.
    /// Its actual aerodynamic loss enters the source parcel once. For each
    /// proposed Q, the existing native thermal network solves all layer
    /// energies, including original-link conduction, per-layer material loss,
    /// reversible Helmholtz/entropy exchange and optional endpoint Joule heat.
    /// The native scalar feedback owner then enforces
    /// Q = dt G(T_outer,end,T_gas,end)(T_outer,end-T_gas,end).
    ///
    /// The same gas EOS/transport endpoint rebinds drag for the next sample and
    /// the current pressure observation. q/p and the previous accepted
    /// acceleration are retained; no pluck, extra epoch or elastic work is added.
    /// Total storage includes all network caloric cells, kinetic energy,
    /// axial/modal Helmholtz energy, EACH layer's T*s, and source gas energy.
    ///
    /// Both callbacks are repeatable proposals and must not mutate external
    /// physical history. resolve_thermal supplies only independent heat and
    /// original-link conductances; this owner supplies Joule heat, material loss,
    /// reversible heat and -Q. resolve_exchange receives the actual outer-layer
    /// caloric endpoint and source gas, retaining its original typed receipt.
    /// Its law owns source, geometry and lumped admission; the concrete
    /// conduction adapter below supplies source pin and pointwise Biot checks.
    ///
    /// The parcel must initially match the runtime's complete gas snapshot.
    /// Contact area fits both the true outer lateral surface and the exposed
    /// thermal cell area. An ambient-HTC boundary on that cell is refused to
    /// prevent double counting. Zero area requires zero conductance.
    ///
    /// Source bounds apply to every bracket/thermal proposal. Accepted solid
    /// and gas temperature increments have separate finite budgets. All
    /// failures preserve the three accepted owners. This is a first-order
    /// partitioned, fixed-geometry, well-mixed, fixed-volume model, without
    /// pressure work, chamber modes, spatial boundary layers or radiation load.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_isochoric_gas_exchange<'t, S, R, E>(
        &mut self, cx: &Cx<'_>, thermal: &mut DepositionSimulation<'t>,
        gas: &mut IsochoricGasParcel, controls: ConcentricFiniteGasSolve,
        energy_guess_j: &[f64], forces: &[f64], gas_heat_j: f64,
        drive: Option<OhmicDrive>, contact: ConcentricGasContact,
        resolve_thermal: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        resolve_exchange: impl Fn(&ResolvedThermochemicalGasState, EquilibriumPhaseState)
            -> Result<ConductanceSample<R>, E>,
    ) -> Result<ConcentricFiniteGasFrame<S, R>, ConcentricStringError<E>> {
        checkpoint(cx)?;
        controls.gas.gas.validate()?;
        feedback::admit::<core::convert::Infallible>(controls.exchange)
            .map_err(|_| invalid("invalid interface heat bracket, tolerance or work limit"))?;
        if !gas_heat_j.is_finite()
            || !positive(controls.gas.maximum_gas_temperature_increment_k)
            || self.gas != *gas.state().state()
        { return Err(invalid("finite gas heat, positive increment limit and identical gas snapshot required")); }
        self.admit_gas_contact(thermal, contact)?;
        let outer = *self.bindings.last().ok_or(invalid("missing outer caloric layer"))?;
        let cells: Vec<_> = self.bindings.iter().map(|binding| binding.cell).collect();
        let initial_gas_energy = gas.state().specific_internal_energy_j_kg();
        let initial_gas_temperature = gas.state().state().temperature;
        let thermal_trials = Cell::new(0usize);
        let solved = feedback::solve_scalar(
            controls.exchange.heat_bracket_j, controls.exchange.source_tolerance_j,
            controls.exchange.maximum_iterations,
            |q, trials, _cancelled| {
                let mut runtime = self.clone();
                let mut next_thermal = thermal.clone();
                let external_solid_heat = Cell::new(f64::NAN);
                let conductance = Cell::new(f64::NAN);
                // Air loss comes from the actual native mechanical advance,
                // not an independent damping estimator. At fixed Q it is
                // invariant across thermal/Newton probes: invert this SAME
                // gas source once, then retain the identical proposal.
                let gas_proposal = RefCell::new(None::<(IsochoricGasParcel, GasHeatStep)>);
                let string = runtime.step_with_resolved_thermal_network(
                    cx, &mut next_thermal, controls.thermal, energy_guess_j, forces,
                    |endpoint, source, dt, loss| {
                        checkpoint(cx)?;
                        self.admit_gas_contact(endpoint, contact)?;
                        let (next_gas, gas_step) = {
                            let mut cached = gas_proposal.borrow_mut();
                            if cached.is_none() {
                                let base = add_heat(gas_heat_j, loss.air_loss_j)?;
                                let supplied = add_heat(base, q)?;
                                let mut next = gas.clone();
                                let report = next.add_heat(supplied, controls.gas.gas,
                                    || cx.checkpoint().is_err())?;
                                *cached = Some((next, report));
                            }
                            cached.as_ref().expect("gas proposal was initialized").clone()
                        };
                        let base = resolve_thermal(endpoint).map_err(ConcentricStringError::Source)?;
                        let (mut coefficients, electrical) = if let Some(drive) = drive {
                            let driven = add_ohmic_heat(cx, endpoint, source,
                                self.string.length_m, dt, &cells, drive,
                                controls.thermal.total_energy_tolerance_j, base)?;
                            let StepCoefficients { sources_w, conductances_w_k, payload } = driven;
                            let ConcentricStringOhmicSource { drive, electrical, joule_heat_j,
                                layer_joule_heat_j, other_external_heat_j, transport } = payload;
                            (StepCoefficients { sources_w, conductances_w_k, payload: transport },
                                Some(ConcentricStringOhmicSource { drive, electrical, joule_heat_j,
                                    layer_joule_heat_j, other_external_heat_j, transport: () }))
                        } else { (base, None) };
                        let check = || cx.checkpoint().map_err(|_| DepositionError::Cancelled);
                        let external = finite(dt*endpoint.external_step_power_w(
                            &coefficients.sources_w, &coefficients.conductances_w_k, check)?)?;
                        deposit(&mut coefficients.sources_w, outer.cell, -q, dt,
                            controls.exchange.energy_tolerance_j)?;
                        let cell = endpoint.cells().get(&outer.cell)
                            .ok_or(invalid("thermal proposal omitted the outer layer"))?;
                        let phase = outer.caloric.curve().state_at_specific_enthalpy(cell.specific_h_j_kg)?;
                        let sample = resolve_exchange(next_gas.state(), phase)
                            .map_err(ConcentricStringError::Source)?;
                        if !sample.conductance_w_k.is_finite() || sample.conductance_w_k < 0.0
                            || (contact.area_m2 == 0.0 && sample.conductance_w_k != 0.0)
                        { return Err(invalid("finite nonnegative surface conductance and actual transfer area required")); }
                        external_solid_heat.set(external);
                        conductance.set(sample.conductance_w_k);
                        Ok(StepCoefficients {
                            sources_w: coefficients.sources_w,
                            conductances_w_k: coefficients.conductances_w_k,
                            payload: ConcentricFiniteGasSource { thermal: coefficients.payload,
                                exchange: sample.source, electrical, gas: gas_step, next_gas },
                        })
                    },
                )?;
                thermal_trials.set(thermal_trials.get().checked_add(string.thermal_trials)
                    .ok_or(invalid("coupled thermal work counter overflow"))?);
                let delta = finite(next_thermal.temperature_k(outer.cell)?
                    -string.source.next_gas.state().state().temperature)?;
                let power = finite(conductance.get()*delta)?;
                let duration = finite(next_thermal.time_s()-thermal.time_s())?;
                let source_heat = finite(duration*power)?;
                let residual = finite(q-source_heat)?;
                if (conductance.get() != 0.0 && delta != 0.0 && power == 0.0)
                    || (power != 0.0 && source_heat == 0.0)
                { return Err(invalid("surface heat law is unrepresentable at the accepted endpoint")); }
                Ok((Candidate { runtime, thermal: next_thermal, string, heat_j: q,
                    conductance_w_k: conductance.get(),
                    external_solid_heat_j: external_solid_heat.get(),
                    source_residual_j: residual, trials }, residual))
            },
            &mut || cx.checkpoint().is_err(),
        ).map_err(scalar_error)?;
        let Candidate { mut runtime, thermal: next_thermal, mut string, heat_j: q,
            conductance_w_k, external_solid_heat_j, source_residual_j, trials } = solved;
        let next_gas = &string.source.next_gas;
        let gas_temperature = next_gas.state().state().temperature;
        let delta = finite(next_thermal.temperature_k(outer.cell)?-gas_temperature)?;
        if (q > 0.0 && delta < 0.0) || (q < 0.0 && delta > 0.0) {
            return Err(invalid("surface heat crossed endpoint temperatures; tighten its source tolerance"));
        }
        if (gas_temperature-initial_gas_temperature).abs()
            > controls.gas.maximum_gas_temperature_increment_k
        { return Err(invalid("finite gas exceeded its per-sample temperature increment")); }

        self.refresh_finite_gas_medium(cx, forces, &mut runtime, &mut string.vibration,
            next_gas.state().state())?;

        let final_gas_energy = next_gas.state().specific_internal_energy_j_kg();
        let gas_storage = finite(gas.mass_kg()*(final_gas_energy-initial_gas_energy))?;
        let mechanical_storage = finite(string.mechanical_internal_energy_after_j
            -string.mechanical_internal_energy_before_j)?;
        let solid_storage = string.thermal.stored_energy_change_j;
        let loss = &string.vibration.dissipation;
        let drag_heat_j = loss.air_loss_j;
        let input = finite(string.vibration.acoustic.input_work_j+external_solid_heat_j+gas_heat_j)?;
        let residual = finite(input-(mechanical_storage+solid_storage+gas_storage
            +loss.authored_loss_j+loss.roundoff_residual_j))?;
        let scale = finite(input.abs()+external_solid_heat_j.abs()+gas_heat_j.abs()
            +q.abs()+drag_heat_j+mechanical_storage.abs()+solid_storage.abs()+gas_storage.abs()
            +gas.mass_kg()*(initial_gas_energy.abs()+final_gas_energy.abs()))?;
        let gas_tolerance = finite(gas.mass_kg()*controls.gas.gas.energy_tolerance_j_kg)?;
        // The deposited solid debit and native gas inversion must also close
        // their own requested interface allowance, separate from source-law Q.
        let exchange_residual = finite(string.external_heat_j-external_solid_heat_j+q
            +string.source.gas.residual_j)?;
        let exchange_roundoff = 64.0*f64::EPSILON*scale;
        if !positive(gas_tolerance) || exchange_residual.abs()
            > controls.exchange.energy_tolerance_j+exchange_roundoff
        { return Err(invalid("actual solid debit and source gas storage exceed the exchange energy allowance")); }
        let tolerance = finite(controls.thermal.total_energy_tolerance_j
            +string.roundoff_tolerance_j+gas_tolerance+controls.exchange.energy_tolerance_j
            +128.0*f64::EPSILON*scale)?;
        if residual.abs() > tolerance {
            return Err(invalid("independent layered mechanical, caloric and finite gas energy do not close"));
        }
        let accepted_gas = next_gas.clone();
        checkpoint(cx)?;
        // No fallible operation or destructor separates the three publications.
        let _old_runtime = core::mem::replace(self, runtime);
        let _old_thermal = core::mem::replace(thermal, next_thermal);
        let _old_gas = core::mem::replace(gas, accepted_gas);
        Ok(ConcentricFiniteGasFrame { string, solid_to_gas_heat_j: q, drag_heat_j,
            external_gas_heat_j: gas_heat_j, external_solid_heat_j, conductance_w_k,
            source_residual_j, energy_balance_residual_j: residual, energy_tolerance_j: tolerance,
            trials, thermal_trials: thermal_trials.get() })
    }

    /// Native series conduction bound to the original outer material source.
    ///
    /// Both thermal endpoints are actual accepted candidates. The original
    /// conductivity receipt and non-temperature conditions must match the
    /// annulus; gas conductivity comes from its same source parcel. The existing
    /// lumped convention Bi = G(V/A)/(A k_outer) must remain at or below 0.1.
    /// Effective gas/solid distances are explicit caller geometry, not a
    /// boundary-layer or radial-profile prediction. Zero conductance is admitted
    /// only with valid original material data.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_isochoric_gas_conduction<'t, S, E>(
        &mut self, cx: &Cx<'_>, thermal: &mut DepositionSimulation<'t>,
        gas: &mut IsochoricGasParcel, controls: ConcentricFiniteGasSolve,
        energy_guess_j: &[f64], forces: &[f64], gas_heat_j: f64,
        drive: Option<OhmicDrive>, layers: &MaterialConductionLayers<'_>,
        resolve_thermal: impl Fn(&DepositionSimulation<'t>)
            -> Result<StepCoefficients<S>, E>,
    ) -> Result<ConcentricFiniteGasFrame<S, ResolvedMaterialStatePoint>,
        ConcentricStringError<E>>
    where E: From<LayerConductionError>
    {
        checkpoint(cx)?;
        let binding = self.admit_finite_gas_conduction::<E>(layers)?;
        self.step_with_isochoric_gas_exchange(cx, thermal, gas, controls, energy_guess_j,
            forces, gas_heat_j, drive, ConcentricGasContact { area_m2: layers.area_m2 },
            resolve_thermal, |gas_end, solid_end|
                finite_gas_conductance(layers, binding, gas_end, solid_end))
    }

    fn admit_finite_gas_conduction<E>(&self, layers: &MaterialConductionLayers<'_>)
        -> Result<(fs_matdb::ClaimId, f64), ConcentricStringError<E>>
    {
        let section = self.source.section.section();
        let outer = section.layers().last().ok_or(invalid("missing outer material"))?;
        let same_conditions = layers.context.axes().iter().filter(|(axis, _)| axis.as_str() != "T")
            .map(|(axis, value)| (axis, value.to_bits()))
            .eq(outer.material.query_point().iter().filter(|(axis, _)| axis != "T")
                .map(|(axis, value)| (axis, value.to_bits())));
        if layers.material.content_hash() != outer.material.card_identity()
            || layers.temperature_axis != "T" || layers.liquid_fraction_axis.is_some()
            || layers.context.axis_quantities() != outer.material.axis_quantities()
            || !same_conditions
        { return Err(invalid("surface conduction must retain the actual outer card, source conditions and T axis")); }
        let selected = outer.material.property_by_key(layers.conductivity)
            .ok_or(invalid("original outer source lacks the declared conductivity"))?
            .answer().receipt.selected;
        let volume = section.layer_integrals().last().ok_or(invalid("missing outer annular geometry"))?
            .area_m2()*self.string.length_m;
        Ok((selected, volume))
    }

    // A candidate may contain a discarded heat-stage observation. Always use
    // the ORIGINAL accepted acceleration, so a physical sample advances the
    // dipole history once even when heat and species transport both rebind gas.
    fn refresh_finite_gas_medium<E>(
        &self, cx: &Cx<'_>, forces: &[f64], runtime: &mut Self,
        vibration: &mut super::MaterialStringFrame, gas: &fs_material::gas::GasState,
    ) -> Result<(), ConcentricStringError<E>> {
        checkpoint(cx)?;
        let mut member = linear_string_member(&runtime.string, None, None,
            gas, runtime.listener_m, runtime.sample_rate_hz, 1.0)?;
        if member.model.modes().len() != runtime.member.model.modes().len()
            || member.model.modes().iter().zip(runtime.member.model.modes()).any(|(a, b)|
                a.angular_frequency_rad_s.to_bits() != b.angular_frequency_rad_s.to_bits())
        { return Err(invalid("gas rebind changed the fixed mechanical energy basis")); }
        member.model.restore_states(runtime.states()).map_err(AcousticRealizeError::Modal)?;
        member.observer.previous_acceleration = self.member.observer.previous_acceleration.clone();
        let initial_acceleration = member.observer.previous_acceleration.is_none()
            .then(|| linear_accelerations(&self.member.model, forces));
        vibration.acoustic.observer_pressure_pa = member.observer.observe(
            linear_accelerations(&member.model, forces), initial_acceleration,
            member.model.sample_period_s())?;
        runtime.member = member;
        runtime.gas = *gas;
        Ok(())
    }

    fn admit_gas_contact<E>(&self, thermal: &DepositionSimulation<'_>,
        contact: ConcentricGasContact) -> Result<(), ConcentricStringError<E>>
    {
        let outer = self.bindings.last().ok_or(invalid("missing outer thermal layer"))?;
        let cell = thermal.cells().get(&outer.cell).ok_or(invalid("missing outer thermal cell"))?;
        let lateral_area = finite(2.0*core::f64::consts::PI
            *self.source.section.section().outer_radius_m()*self.string.length_m)?;
        if !contact.area_m2.is_finite() || contact.area_m2 < 0.0
            || contact.area_m2 > lateral_area
            || contact.area_m2 > thermal.exposed_area_m2(outer.cell)?
            || cell.boundary.htc_w_m2_k != 0.0
        { return Err(invalid("gas contact exceeds the exposed outer lateral area or overlaps ambient convection")); }
        Ok(())
    }
}

fn add_heat<E>(a: f64, b: f64) -> Result<f64, ConcentricStringError<E>> {
    let result = finite(a+b)?;
    if (a != 0.0 && result == b) || (b != 0.0 && result == a) {
        return Err(invalid("finite gas heat contribution disappears when combined"));
    }
    Ok(result)
}

fn scalar_error<E>(error: feedback::ScalarError<ConcentricStringError<E>>) -> ConcentricStringError<E> {
    match error {
        feedback::ScalarError::Evaluation(error) => error,
        feedback::ScalarError::Invalid(message) => invalid(message),
        feedback::ScalarError::Cancelled => ConcentricStringError::Acoustic(AcousticRealizeError::Cancelled),
        feedback::ScalarError::OutsideBracket { lower, upper } =>
            ConcentricStringError::ExchangeOutsideBracket { lower_residual_j: lower, upper_residual_j: upper },
        feedback::ScalarError::NotConverged { iterations, residual } =>
            ConcentricStringError::ExchangeNotConverged { iterations, residual_j: residual },
    }
}
