//! A thermally resistive deposit on a separately stored, nonvolatile substrate.
//!
//! The two original caloric charts retain distinct temperatures and enthalpy
//! datums. A caller-supplied conductance or original nonlinear passive heat law
//! represents the admitted support/deposit interface. No conductivity, thickness,
//! porosity, moving area or morphology is inferred from deposited mass.
//!
//! Both initial and accepted deposits must have positive mass. Bare nucleation
//! and complete extinction need an explicit interface-transition owner; this
//! operator refuses those boundaries instead of inventing a thermal seed.
//! The vapor reservoir retains the original thermostatted semantics. A finite
//! nonisothermal gas owner can use the heat-only kernel after its mass transfer.

/// Mass-owned coating dimensions and through-layer/contact thermal conductance.
pub mod geometry;

use super::super::{
    Endpoint, SolidSurfaceState, SolverConfig, SublimationError, SublimationLaw, SurfaceKinetics,
    SurfaceStep, VaporReservoir, finish_surface_step, finite, nonnegative, positive, properties,
    solve_endpoint_with_scale, transfer, validate_solver, validate_surface,
};
use core::cell::Cell;

/// Original constitutive heat flow between the separately stored bodies.
///
/// Positive power flows from support to deposit. The law must be continuous,
/// nondecreasing in support temperature and nonincreasing in deposit temperature
/// throughout its admitted domain. It must be passive: heat flows down the
/// temperature difference and is exactly zero at equal temperatures. Each
/// evaluated finite endpoint is checked for those sign/equality conditions;
/// finite sampling does not establish global monotonicity or source validity.
///
/// The law owns conductivity integrals, contact/radiative physics and any source
/// query context. Its own inner work and cancellation must be bounded if a query
/// is nontrivial. The outer solver charges each endpoint trial to its original
/// work allowance; it neither estimates a conductance nor invokes another
/// caloric solver. Temperatures are kelvin and the result is watts.
pub trait ResistiveInterfaceLaw {
    fn heat_rate_w(
        &self,
        substrate_temperature_k: f64,
        deposit_temperature_k: f64,
    ) -> Result<f64, SublimationError>;
}

impl<F: Fn(f64, f64) -> Result<f64, SublimationError> + ?Sized> ResistiveInterfaceLaw for F {
    fn heat_rate_w(
        &self,
        substrate_temperature_k: f64,
        deposit_temperature_k: f64,
    ) -> Result<f64, SublimationError> {
        self(substrate_temperature_k, deposit_temperature_k)
    }
}

#[derive(Clone, Copy)]
enum Interface<'a> {
    Conductance(f64),
    Nonlinear(&'a dyn ResistiveInterfaceLaw),
}

/// Independently stored substrate and existing deposit; geometry is caller-owned.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResistiveSurfaceState {
    pub substrate_mass_kg: f64,
    pub substrate_specific_enthalpy_j_kg: f64,
    pub deposit: SolidSurfaceState,
}

/// Fixed-mass heat exchange, including both separately represented caloric states.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResistiveHeatStep {
    pub state: ResistiveSurfaceState,
    pub substrate_temperature_k: f64,
    pub deposit_temperature_k: f64,
    /// Heat from substrate to deposit [J]; counted with opposite signs locally.
    pub interface_heat_j: f64,
    pub substrate_energy_change_j: f64,
    pub deposit_energy_change_j: f64,
    pub substrate_energy_residual_j: f64,
    pub deposit_energy_residual_j: f64,
    /// Actual total storage change minus the two external heat inputs [J].
    pub energy_residual_j: f64,
    /// All scalar balance trials, including refused chart/anchor trials.
    pub evaluations: usize,
}

/// One conservative phase/thermal endpoint with separate material temperatures.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResistiveSurfaceStep {
    pub state: ResistiveSurfaceState,
    /// Deposit-only receipt. Its supplied heat includes the interface heat;
    /// its residual includes deposit storage and the complete vapor enthalpy.
    pub phase: SurfaceStep,
    pub substrate_temperature_k: f64,
    pub interface_heat_j: f64,
    pub substrate_energy_change_j: f64,
    pub deposit_energy_change_j: f64,
    pub substrate_energy_residual_j: f64,
    /// Actual combined storage plus vapor energy minus external heat [J].
    pub energy_residual_j: f64,
    pub mass_residual_kg: f64,
    pub evaluations: usize,
}

/// Implicit two-body heat exchange at fixed positive substrate/deposit masses.
///
/// The original continuous nondecreasing `T_s(h_s)` and `T_d(h_d)` charts close
/// `Delta U_s = Q_s - q`, `Delta U_d = Q_d + q`, with
/// `q = duration_s * conductance_w_k * (T_s1 - T_d1)`. The deposit law admits
/// every trial deposit temperature; the substrate retains its independent
/// domain. Zero conductance gives independent caloric heat storage. At zero
/// duration external heat impulses remain permitted, but no contact heat flows.
///
/// Total energy eliminates the substrate enthalpy, leaving ONE bounded scalar
/// solve, without a frozen/effective heat capacity, nested nonlinear solve or
/// allocation. Cancellation and exhausted work publish no state.
#[allow(clippy::too_many_arguments)]
pub fn heat_resistive_surface<L: SublimationLaw + ?Sized>(
    law: &L,
    state: ResistiveSurfaceState,
    duration_s: f64,
    conductance_w_k: f64,
    substrate_heat_j: f64,
    deposit_heat_j: f64,
    config: SolverConfig,
    substrate_temperature: impl Fn(f64) -> Result<f64, SublimationError>,
    deposit_temperature: impl Fn(f64) -> Result<f64, SublimationError>,
    checkpoint: impl FnMut() -> Result<(), SublimationError>,
) -> Result<ResistiveHeatStep, SublimationError> {
    let (accepted, evaluations) = advance(
        law,
        state,
        None,
        duration_s,
        Interface::Conductance(conductance_w_k),
        substrate_heat_j,
        deposit_heat_j,
        config,
        substrate_temperature,
        deposit_temperature,
        checkpoint,
    )?;
    Ok(heat_step(accepted, evaluations))
}

/// The same original two-body caloric solve with an endpoint nonlinear heat law.
///
/// `q = duration_s * interface.heat_rate_w(T_s1, T_d1)` closes the same separate
/// energy balances. It does not freeze a secant conductance or differentiate the
/// supplied law. Total energy still eliminates the support enthalpy, leaving
/// one bounded scalar solve. At zero duration no interface query is made, and
/// authored heat impulses are stored independently. Both masses stay positive.
#[allow(clippy::too_many_arguments)]
pub fn heat_resistive_surface_with_interface<L: SublimationLaw + ?Sized>(
    law: &L,
    state: ResistiveSurfaceState,
    duration_s: f64,
    interface: &dyn ResistiveInterfaceLaw,
    substrate_heat_j: f64,
    deposit_heat_j: f64,
    config: SolverConfig,
    substrate_temperature: impl Fn(f64) -> Result<f64, SublimationError>,
    deposit_temperature: impl Fn(f64) -> Result<f64, SublimationError>,
    checkpoint: impl FnMut() -> Result<(), SublimationError>,
) -> Result<ResistiveHeatStep, SublimationError> {
    let (accepted, evaluations) = advance(
        law,
        state,
        None,
        duration_s,
        Interface::Nonlinear(interface),
        substrate_heat_j,
        deposit_heat_j,
        config,
        substrate_temperature,
        deposit_temperature,
        checkpoint,
    )?;
    Ok(heat_step(accepted, evaluations))
}

fn heat_step(accepted: Balance, evaluations: usize) -> ResistiveHeatStep {
    ResistiveHeatStep {
        state: accepted.state,
        substrate_temperature_k: accepted.substrate_temperature_k,
        deposit_temperature_k: accepted.endpoint.t,
        interface_heat_j: accepted.interface_heat_j,
        substrate_energy_change_j: accepted.substrate_energy_change_j,
        deposit_energy_change_j: accepted.deposit_energy_change_j,
        substrate_energy_residual_j: accepted.substrate_energy_residual_j,
        deposit_energy_residual_j: accepted.deposit_energy_residual_j,
        energy_residual_j: accepted.energy_residual_j,
        evaluations,
    }
}

/// Phase exchange at the deposit's temperature with finite substrate coupling.
///
/// For positive mass transfer `dm` into vapor, the deposit balance is
/// `m_d0*(h_d1-h_d0) + dm*L(T_d1) = Q_d + q`. The substrate balance and contact
/// relation are the same as [`heat_resistive_surface`]. The original molecular
/// and optional gas-film laws determine `dm`; a finite bath uses its original
/// frozen-temperature back-pressure relaxation. Every mass and energy account
/// derives from that one accepted transfer, and vapor enthalpy includes latent
/// heat once. This is a first-order two-node layer approximation.
#[allow(clippy::too_many_arguments)]
pub fn advance_resistive_surface_with_heat<L: SublimationLaw + ?Sized>(
    law: &L,
    kinetics: SurfaceKinetics,
    state: ResistiveSurfaceState,
    reservoir: VaporReservoir,
    duration_s: f64,
    conductance_w_k: f64,
    substrate_heat_j: f64,
    deposit_heat_j: f64,
    config: SolverConfig,
    substrate_temperature: impl Fn(f64) -> Result<f64, SublimationError>,
    deposit_temperature: impl Fn(f64) -> Result<f64, SublimationError>,
    checkpoint: impl FnMut() -> Result<(), SublimationError>,
) -> Result<ResistiveSurfaceStep, SublimationError> {
    advance_phase(
        law,
        kinetics,
        state,
        reservoir,
        duration_s,
        Interface::Conductance(conductance_w_k),
        substrate_heat_j,
        deposit_heat_j,
        config,
        substrate_temperature,
        deposit_temperature,
        checkpoint,
    )
}

/// Original phase exchange and nonlinear interfacial heat in the same endpoint.
///
/// Mass transfer, complete vapor enthalpy and latent heat retain the meanings of
/// [`advance_resistive_surface_with_heat`]. Only the contact relation changes to
/// the original [`ResistiveInterfaceLaw`] evaluated at both accepted material
/// temperatures. There is one scalar solve and one shared trial budget, with no
/// frozen/effective heat capacity or fitted conductance. This remains a first-
/// order two-node approximation; it supplies no spatial temperature profile.
#[allow(clippy::too_many_arguments)]
pub fn advance_resistive_surface_with_interface<L: SublimationLaw + ?Sized>(
    law: &L,
    kinetics: SurfaceKinetics,
    state: ResistiveSurfaceState,
    reservoir: VaporReservoir,
    duration_s: f64,
    interface: &dyn ResistiveInterfaceLaw,
    substrate_heat_j: f64,
    deposit_heat_j: f64,
    config: SolverConfig,
    substrate_temperature: impl Fn(f64) -> Result<f64, SublimationError>,
    deposit_temperature: impl Fn(f64) -> Result<f64, SublimationError>,
    checkpoint: impl FnMut() -> Result<(), SublimationError>,
) -> Result<ResistiveSurfaceStep, SublimationError> {
    advance_phase(
        law,
        kinetics,
        state,
        reservoir,
        duration_s,
        Interface::Nonlinear(interface),
        substrate_heat_j,
        deposit_heat_j,
        config,
        substrate_temperature,
        deposit_temperature,
        checkpoint,
    )
}

#[allow(clippy::too_many_arguments)]
fn advance_phase<L: SublimationLaw + ?Sized>(
    law: &L,
    kinetics: SurfaceKinetics,
    state: ResistiveSurfaceState,
    reservoir: VaporReservoir,
    duration_s: f64,
    interface: Interface<'_>,
    substrate_heat_j: f64,
    deposit_heat_j: f64,
    config: SolverConfig,
    substrate_temperature: impl Fn(f64) -> Result<f64, SublimationError>,
    deposit_temperature: impl Fn(f64) -> Result<f64, SublimationError>,
    mut checkpoint: impl FnMut() -> Result<(), SublimationError>,
) -> Result<ResistiveSurfaceStep, SublimationError> {
    let (accepted, evaluations) = advance(
        law,
        state,
        Some((kinetics, reservoir)),
        duration_s,
        interface,
        substrate_heat_j,
        deposit_heat_j,
        config,
        substrate_temperature,
        deposit_temperature,
        &mut checkpoint,
    )?;
    let mut phase = finish_surface_step(
        state.deposit,
        reservoir,
        accepted.endpoint,
        finite(deposit_heat_j + accepted.interface_heat_j)?,
    )?;
    phase.energy_residual_j = accepted.deposit_energy_residual_j;
    let mut mass_residual_kg =
        finite((phase.solid.mass_kg - state.deposit.mass_kg) + phase.transferred_mass_kg)?;
    if let (
        VaporReservoir::Finite {
            vapor_mass_kg: before,
            ..
        },
        VaporReservoir::Finite {
            vapor_mass_kg: after,
            ..
        },
    ) = (reservoir, phase.reservoir)
    {
        if phase.transferred_mass_kg != 0.0 && before == after {
            return Err(SublimationError::Invalid(
                "unrepresentable resistive vapor increment",
            ));
        }
        mass_residual_kg =
            finite(mass_residual_kg + ((after - before) - phase.transferred_mass_kg))?;
    }
    checkpoint()?;
    Ok(ResistiveSurfaceStep {
        state: accepted.state,
        phase,
        substrate_temperature_k: accepted.substrate_temperature_k,
        interface_heat_j: accepted.interface_heat_j,
        substrate_energy_change_j: accepted.substrate_energy_change_j,
        deposit_energy_change_j: accepted.deposit_energy_change_j,
        substrate_energy_residual_j: accepted.substrate_energy_residual_j,
        energy_residual_j: accepted.energy_residual_j,
        mass_residual_kg,
        evaluations,
    })
}

#[derive(Clone, Copy)]
struct Balance {
    endpoint: Endpoint,
    state: ResistiveSurfaceState,
    substrate_temperature_k: f64,
    interface_heat_j: f64,
    substrate_energy_change_j: f64,
    deposit_energy_change_j: f64,
    substrate_energy_residual_j: f64,
    deposit_energy_residual_j: f64,
    energy_residual_j: f64,
    energy_scale_j: f64,
}

#[allow(clippy::too_many_arguments)]
fn advance<L: SublimationLaw + ?Sized>(
    law: &L,
    state: ResistiveSurfaceState,
    phase: Option<(SurfaceKinetics, VaporReservoir)>,
    duration_s: f64,
    interface: Interface<'_>,
    substrate_heat_j: f64,
    deposit_heat_j: f64,
    config: SolverConfig,
    substrate_temperature: impl Fn(f64) -> Result<f64, SublimationError>,
    deposit_temperature: impl Fn(f64) -> Result<f64, SublimationError>,
    mut checkpoint: impl FnMut() -> Result<(), SublimationError>,
) -> Result<(Balance, usize), SublimationError> {
    if state.deposit.mass_kg == 0.0 {
        return Err(SublimationError::Invalid(
            "bare resistive coating needs an interface transition",
        ));
    }
    validate_surface(state.deposit, duration_s)?;
    validate_solver(config)?;
    if !positive(state.substrate_mass_kg)
        || !state.substrate_specific_enthalpy_j_kg.is_finite()
        || matches!(interface, Interface::Conductance(g) if !nonnegative(g))
        || !substrate_heat_j.is_finite()
        || !deposit_heat_j.is_finite()
    {
        return Err(SublimationError::Invalid(
            "invalid resistive substrate, conductance or heat",
        ));
    }
    checkpoint()?;
    if !positive(substrate_temperature(
        state.substrate_specific_enthalpy_j_kg,
    )?) {
        return Err(SublimationError::Invalid("invalid substrate temperature"));
    }
    let initial_temperature = deposit_temperature(state.deposit.specific_enthalpy_j_kg)?;
    let initial_properties = properties(law, initial_temperature)?;
    let initial_transfer = match phase {
        Some((kinetics, reservoir)) => transfer(
            initial_properties,
            initial_temperature,
            kinetics,
            state.deposit,
            reservoir,
            duration_s,
        )?,
        None => 0.0,
    };
    let initial_latent = finite(initial_transfer * initial_properties.latent_heat_j_kg)?;
    // Keep the original multiplication order and exact insulated branch for
    // authored constant conductance. A nonlinear law supplies POWER instead.
    let contact = match interface {
        Interface::Conductance(g) => Some(finite(duration_s * g)?),
        Interface::Nonlinear(_) => None,
    };
    let inactive_contact = duration_s == 0.0 || contact == Some(0.0);
    let total_heat = finite(substrate_heat_j + deposit_heat_j)?;
    let evaluations = Cell::new(0usize);
    let last = Cell::new(None::<Balance>);
    let invalid_interface = Cell::new(None::<SublimationError>);
    let mut evaluate = |h: f64| -> Result<Endpoint, SublimationError> {
        checkpoint()?;
        if evaluations.get() >= config.max_iterations {
            return Err(SublimationError::NonConverged);
        }
        evaluations.set(evaluations.get() + 1);
        finite(h)?;
        let t = deposit_temperature(h)?;
        let p = properties(law, t)?;
        let dm = match phase {
            Some((kinetics, reservoir)) => {
                transfer(p, t, kinetics, state.deposit, reservoir, duration_s)?
            }
            None => 0.0,
        };
        let mass = finite(state.deposit.mass_kg - dm)?;
        let sensible = finite(state.deposit.mass_kg * (h - state.deposit.specific_enthalpy_j_kg))?;
        let latent = finite(dm * p.latent_heat_j_kg)?;
        let needed = finite(sensible + latent - deposit_heat_j)?;
        let substrate_heat = if inactive_contact {
            substrate_heat_j
        } else {
            substrate_heat_j - needed
        };
        let substrate_h = if substrate_heat == 0.0 {
            state.substrate_specific_enthalpy_j_kg
        } else {
            finite(
                state.substrate_specific_enthalpy_j_kg
                    + finite(substrate_heat / state.substrate_mass_kg)?,
            )?
        };
        let substrate_t = substrate_temperature(substrate_h)?;
        if !positive(substrate_t) {
            return Err(SublimationError::Invalid("invalid substrate temperature"));
        }
        let interface = if let Some(contact) = contact {
            finite(contact * (substrate_t - t))?
        } else if duration_s == 0.0 {
            0.0
        } else {
            let Interface::Nonlinear(law) = interface else {
                unreachable!()
            };
            let rate = law.heat_rate_w(substrate_t, t)?;
            if !rate.is_finite()
                || (substrate_t == t && rate != 0.0)
                || (substrate_t > t && rate < 0.0)
                || (substrate_t < t && rate > 0.0)
            {
                let error = SublimationError::Invalid(
                    "nonlinear interface heat must be finite, passive and zero at equal temperatures",
                );
                // A broken heat law is not a caloric-domain bracket boundary.
                // Abort the shared solver immediately, then restore this exact
                // physical refusal rather than hiding it behind NonConverged.
                invalid_interface.set(Some(error));
                return Err(SublimationError::NonConverged);
            }
            finite(duration_s * rate)?
        };
        let substrate_storage = finite(
            state.substrate_mass_kg * (substrate_h - state.substrate_specific_enthalpy_j_kg),
        )?;
        let actual_dm = finite(mass - state.deposit.mass_kg)?;
        let deposit_storage = finite(sensible + actual_dm * h)?;
        let mass_roundoff_energy = finite((actual_dm + dm) * h)?;
        let substrate_residual = finite(substrate_storage - substrate_heat_j + interface)?;
        let deposit_residual = finite(needed - interface + mass_roundoff_energy)?;
        let energy_residual =
            finite(substrate_storage - substrate_heat_j + needed + mass_roundoff_energy)?;
        let endpoint = Endpoint {
            h,
            t,
            p,
            dm,
            residual: finite(needed - interface)?,
        };
        last.set(Some(Balance {
            endpoint,
            state: ResistiveSurfaceState {
                substrate_specific_enthalpy_j_kg: substrate_h,
                deposit: SolidSurfaceState {
                    mass_kg: mass,
                    specific_enthalpy_j_kg: h,
                    ..state.deposit
                },
                ..state
            },
            substrate_temperature_k: substrate_t,
            interface_heat_j: interface,
            substrate_energy_change_j: substrate_storage,
            deposit_energy_change_j: deposit_storage,
            substrate_energy_residual_j: substrate_residual,
            deposit_energy_residual_j: deposit_residual,
            energy_residual_j: energy_residual,
            energy_scale_j: substrate_heat_j
                .abs()
                .max(deposit_heat_j.abs())
                .max(substrate_storage.abs())
                .max(sensible.abs())
                .max(latent.abs())
                .max(interface.abs()),
        }));
        Ok(endpoint)
    };
    // Source-only allocation, unchanged deposit, and shared specific-enthalpy
    // allocation are bounded physical anchors. A refused predictor cannot spend
    // another budget or relax either original caloric domain.
    let h0 = state.deposit.specific_enthalpy_j_kg;
    let anchors = [
        h0 + (deposit_heat_j - initial_latent) / state.deposit.mass_kg,
        h0,
        h0 + (total_heat - initial_latent) / (state.substrate_mass_kg + state.deposit.mass_kg),
    ];
    let mut initial = None;
    let mut refusal = SublimationError::NonConverged;
    for (index, h) in anchors.into_iter().enumerate() {
        if !h.is_finite() || anchors[..index].contains(&h) {
            continue;
        }
        match evaluate(h) {
            Ok(value) => {
                initial = Some(value);
                break;
            }
            Err(SublimationError::Cancelled) => return Err(SublimationError::Cancelled),
            Err(SublimationError::NonConverged) => {
                return Err(invalid_interface
                    .get()
                    .unwrap_or(SublimationError::NonConverged));
            }
            Err(error) => refusal = error,
        }
    }
    let mut initial = initial.ok_or(refusal)?;
    let mut scale = substrate_heat_j
        .abs()
        .max(deposit_heat_j.abs())
        .max(initial_latent.abs());
    loop {
        let endpoint =
            solve_endpoint_with_scale(initial, state.deposit.mass_kg, scale, config, &mut evaluate)
                .map_err(|error| invalid_interface.get().unwrap_or(error))?;
        let accepted = last.get().ok_or(SublimationError::NonConverged)?;
        debug_assert_eq!(endpoint.h, accepted.endpoint.h);
        let tolerance = finite(
            config.absolute_energy_tolerance_j
                + config.relative_tolerance * accepted.energy_scale_j,
        )?;
        if accepted.energy_residual_j.abs() > tolerance {
            return Err(SublimationError::Invalid(
                "unrepresentable resistive storage increment",
            ));
        }
        if accepted.substrate_energy_residual_j.abs() <= tolerance
            && accepted.deposit_energy_residual_j.abs() <= tolerance
        {
            if !positive(accepted.state.deposit.mass_kg) {
                return Err(SublimationError::Invalid(
                    "resistive coating extinction needs an interface transition",
                ));
            }
            if endpoint.dm != 0.0 && accepted.state.deposit.mass_kg == state.deposit.mass_kg {
                return Err(SublimationError::Invalid(
                    "unrepresentable resistive deposit increment",
                ));
            }
            checkpoint()?;
            return Ok((accepted, evaluations.get()));
        }
        // A stiff contact's initial trial heat can greatly exceed the actual
        // finite stored heat. Polish against the accepted physical scale, using
        // the SAME work allowance, instead of admitting that predictor's scale.
        if endpoint.h == initial.h {
            return Err(SublimationError::NonConverged);
        }
        initial = endpoint;
        scale = accepted.energy_scale_j;
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
#[path = "resistance/nonlinear_tests.rs"]
mod nonlinear_tests;
