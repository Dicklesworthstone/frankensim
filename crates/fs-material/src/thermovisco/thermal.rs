//! Dissipation/heat feedback using the existing Maxwell update and thermal owner.
//!
//! One accepted mechanical history and one accepted caloric state are held fixed
//! through all within-step trials. Only physical dashpot dissipation is supplied
//! to transport. Backward-Euler numerical dissipation remains a separate loss.
//!
//! This rung requires temperature-INDEPENDENT spring energies and free strain,
//! a strictly increasing single-solid-phase caloric chart, and constant density.
//! Relaxation times and heat capacity may depend on temperature. Under those
//! restrictions, changes in the chart's h equal caloric internal-energy changes
//! at fixed pressure/density; its arbitrary reference cancels. Variable elastic
//! free energy, thermal eigenstrain, phase changes and volume work need a more
//! complete caloric law and are refused, not assigned invented heat sources.

use core::{convert::Infallible, fmt};

use super::{Energy, Point};
use crate::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState,
    PhaseStateError, UniformEnthalpyStep, UniformEnthalpyStepInput};

/// Hard work ceiling; a failed fixed point is not silently accepted.
pub const MAX_FEEDBACK_ITERATIONS: usize = 256;

/// Independent numerical gates, not a physical accuracy certificate.
#[derive(Clone, Copy, Debug)]
pub struct FeedbackConfig {
    /// Thermal proposals allowed in this step, in 1..=256.
    pub max_iterations: usize,
    /// Absolute endpoint-temperature consistency tolerance [K].
    pub temperature_tolerance_k: f64,
    /// Absolute allowance for heat-source and combined energy defects [J].
    pub energy_tolerance_j: f64,
    /// Fixed-point temperature relaxation in (0,1].
    pub relaxation: f64,
}

/// The original material, caloric or transport failure remains distinguishable.
#[derive(Debug)]
pub enum FeedbackError<E = Infallible> {
    /// Original Maxwell law refusal.
    Mechanical(super::Error),
    /// Original caloric-chart refusal.
    Caloric(PhaseStateError),
    /// Original transport-owner refusal, retaining its concrete type.
    Transport(E),
    /// Inconsistent state, unsupported physics or unrepresentable arithmetic.
    Invalid(&'static str),
    /// The old mechanical and caloric histories have not been changed.
    Cancelled,
    /// Work budget exhausted without a consistent physical endpoint.
    NotConverged {
        iterations: usize,
        temperature_defect_k: f64,
        heat_defect_j: f64,
        thermal_residual_j: f64,
        combined_residual_j: f64,
    },
}
impl<E: fmt::Display> fmt::Display for FeedbackError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Mechanical(e) => write!(f, "thermovisco heat feedback: {e}"),
            Self::Caloric(e) => write!(f, "thermovisco heat feedback: {e}"),
            Self::Transport(e) => write!(f, "thermovisco heat transport: {e}"),
            Self::Invalid(message) => write!(f, "thermovisco heat feedback: {message}"),
            Self::Cancelled => f.write_str("thermovisco heat feedback: cancelled"),
            Self::NotConverged { iterations, temperature_defect_k, heat_defect_j,
                thermal_residual_j, combined_residual_j } => write!(f,
                "thermovisco heat feedback did not converge in {iterations} proposals: dT={temperature_defect_k} K, dQ={heat_defect_j} J, thermal={thermal_residual_j} J, combined={combined_residual_j} J"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for FeedbackError<E> {}

/// One uniformly strained/isothermal specimen with fixed mass and volume.
///
/// The mechanical law and caloric chart are independently declared models;
/// binding them does not invent shared material-card provenance or qualify a
/// specimen. Transport owns the isothermal/Biot admission for its boundary law.
#[derive(Clone, Debug, PartialEq)]
pub struct CaloricPoint<'a> {
    point: Point<'a>,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    state: EquilibriumPhaseState,
    mass_kg: f64,
    volume_m3: f64,
    surface_area_m2: f64,
    time_s: f64,
}

/// A complete proposed endpoint and the actual accepted transport-owner result.
/// No external reservoir may be mutated until this proposal is accepted.
#[derive(Debug)]
pub struct FeedbackStep<'a, R> {
    pub next: CaloricPoint<'a>,
    /// Original Maxwell endpoint audit, in J/m3; not recomputed by this adapter.
    pub mechanical_energy: Energy,
    /// Final thermal proposal, including its concrete report and boundary heat.
    pub transport: UniformEnthalpyStep<R>,
    pub iterations: usize,
    /// Physical dissipation supplied to the final thermal proposal [J].
    pub applied_heat_j: f64,
    /// Mismatch against dissipation reevaluated at the final temperature [J].
    pub heat_defect_j: f64,
    pub temperature_defect_k: f64,
    /// m delta_h - final physical dissipation - external heat [J].
    pub thermal_residual_j: f64,
    /// delta_elastic + m delta_h + numerical_loss - mechanical_work - heat [J].
    pub combined_residual_j: f64,
}

impl<'a> CaloricPoint<'a> {
    /// Bind an existing mechanical history without resetting its branch strains.
    /// Its temperature must exactly match the supplied caloric coordinate.
    ///
    /// # Errors
    /// Refuses non-solid/variable-density caloric charts, temperature-dependent
    /// spring energy or free strain, inconsistent mass/volume, and invalid state.
    pub fn try_new(point: Point<'a>, curve: &'a EquilibriumEnthalpyPhaseCurve,
        specific_enthalpy_j_kg: f64, mass_kg: f64, volume_m3: f64,
        surface_area_m2: f64) -> Result<Self, FeedbackError> {
        for value in [mass_kg, volume_m3, surface_area_m2] {
            if !positive(value) { return Err(FeedbackError::Invalid("mass, volume and area must be finite and positive")); }
        }
        if curve.knots().len() > crate::phase::MAX_HEAT_CAPACITY_ENTHALPY_KNOTS {
            return Err(FeedbackError::Invalid("caloric chart exceeds the bounded knot count"));
        }
        let knots = &point.law.knots;
        let reference = knots[0];
        if knots.iter().any(|k| k.modulus_scale != reference.modulus_scale
            || k.free_strain != reference.free_strain) {
            return Err(FeedbackError::Invalid("temperature-dependent spring energy or eigenstrain needs a coupled caloric law"));
        }
        let density = curve.knots()[0].bulk_density_kg_m3;
        let represented_mass = product(volume_m3, density)?;
        if (represented_mass / mass_kg - 1.0).abs() > 32.0 * f64::EPSILON
            || curve.knots().iter().any(|k| k.liquid_mass_fraction != 0.0
                || k.bulk_density_kg_m3 != density)
            || curve.knots().windows(2).any(|pair| pair[0].temperature_k >= pair[1].temperature_k) {
            return Err(FeedbackError::Invalid("heat feedback needs a fixed-density, strictly increasing solid chart matching reference mass"));
        }
        let state = curve.state_at_specific_enthalpy(specific_enthalpy_j_kg)
            .map_err(FeedbackError::Caloric)?;
        if state.temperature_k() != point.temperature_k()
            || !point.stored_energy_j_m3().is_finite() {
            return Err(FeedbackError::Invalid("mechanical history and caloric temperature must match"));
        }
        Ok(Self { point, curve, state, mass_kg, volume_m3, surface_area_m2, time_s: 0.0 })
    }

    #[must_use]
    pub const fn mechanical(&self) -> &Point<'a> { &self.point }
    #[must_use]
    pub const fn caloric_state(&self) -> &EquilibriumPhaseState { &self.state }
    #[must_use]
    pub const fn caloric_curve(&self) -> &'a EquilibriumEnthalpyPhaseCurve { self.curve }
    #[must_use]
    pub const fn mass_kg(&self) -> f64 { self.mass_kg }
    #[must_use]
    pub const fn volume_m3(&self) -> f64 { self.volume_m3 }
    #[must_use]
    pub const fn time_s(&self) -> f64 { self.time_s }

    /// Propose one prescribed-strain step with implicit damping/heat feedback.
    ///
    /// Pass an existing thermal owner, e.g. `|input| environment.advance(cx,input)`.
    /// The callback must be side-effect-free with respect to physical history:
    /// EVERY trial uses the same old enthalpy, old Maxwell history and duration.
    /// Only its temperature guess changes. The final mechanical update uses the
    /// actual returned thermal temperature, never a stale fixed-point iterate.
    ///
    /// This is bounded fixed-point iteration, not unconditional convergence or
    /// stability. Temperature and physical-heat defects are independently gated.
    /// Numerical dissipation is never deposited as heat. Failure or cancellation
    /// returns no candidate; the original specimen and histories are immutable.
    pub fn advance<R, E>(&self, total_strain: f64, duration_s: f64,
        config: FeedbackConfig,
        mut transport: impl FnMut(UniformEnthalpyStepInput<'_>) -> Result<UniformEnthalpyStep<R>, E>,
        mut cancelled: impl FnMut() -> bool) -> Result<FeedbackStep<'a, R>, FeedbackError<E>> {
        poll(&mut cancelled)?;
        let final_time = self.time_s + duration_s;
        let material_tolerance = config.energy_tolerance_j / self.volume_m3;
        if config.max_iterations == 0 || config.max_iterations > MAX_FEEDBACK_ITERATIONS
            || !positive(config.temperature_tolerance_k) || !positive(config.energy_tolerance_j)
            || !positive(config.relaxation) || config.relaxation > 1.0
            || !positive(duration_s) || !positive(material_tolerance)
            || !total_strain.is_finite() || !final_time.is_finite() || final_time <= self.time_s {
            return Err(FeedbackError::Invalid("invalid feedback controls, strain or physical interval"));
        }
        let mut temperature = self.point.temperature_k();
        let mut defects = [f64::INFINITY; 4];
        for iteration in 1..=config.max_iterations {
            poll(&mut cancelled)?;
            let trial = self.point.advance(temperature, total_strain, duration_s, material_tolerance)
                .map_err(FeedbackError::Mechanical)?;
            let applied_heat_j = product(self.volume_m3, trial.energy.dissipated)?;
            let proposal = transport(UniformEnthalpyStepInput { curve: self.curve, initial: self.state,
                mass_kg: self.mass_kg, volume_m3: self.volume_m3, surface_area_m2: self.surface_area_m2,
                internal_heat_j: applied_heat_j, duration_s }).map_err(FeedbackError::Transport)?;
            poll(&mut cancelled)?;
            let checked = self.curve.state_at_specific_enthalpy(proposal.state.specific_enthalpy_j_kg())
                .map_err(FeedbackError::Caloric)?;
            if checked != proposal.state || !proposal.external_heat_j.is_finite()
                || !proposal.energy_residual_tolerance_j.is_finite()
                || proposal.energy_residual_tolerance_j < 0.0
                || proposal.energy_residual_tolerance_j > config.energy_tolerance_j {
                return Err(FeedbackError::Invalid("transport returned a foreign state or an inadmissible energy allowance"));
            }
            let thermal_change = product(self.mass_kg,
                proposal.state.specific_enthalpy_j_kg() - self.state.specific_enthalpy_j_kg())?;
            let transport_residual = finite(thermal_change - applied_heat_j - proposal.external_heat_j)?;
            if transport_residual.abs() > proposal.energy_residual_tolerance_j {
                return Err(FeedbackError::Invalid("transport proposal violates its own physical heat balance"));
            }
            let final_temperature = proposal.state.temperature_k();
            let mechanical = self.point.advance(final_temperature, total_strain, duration_s, material_tolerance)
                .map_err(FeedbackError::Mechanical)?;
            let actual_heat = product(self.volume_m3, mechanical.energy.dissipated)?;
            let thermal_residual = finite(thermal_change - actual_heat - proposal.external_heat_j)?;
            let combined_residual = finite(product(self.volume_m3, mechanical.energy.stored_change)?
                + thermal_change + product(self.volume_m3, mechanical.energy.numerical_dissipation)?
                - product(self.volume_m3, mechanical.energy.mechanical_work)? - proposal.external_heat_j)?;
            defects = [(final_temperature - temperature).abs(), (actual_heat - applied_heat_j).abs(),
                thermal_residual, combined_residual];
            if defects[0] <= config.temperature_tolerance_k && defects[1] <= config.energy_tolerance_j
                && defects[2].abs() <= config.energy_tolerance_j && defects[3].abs() <= config.energy_tolerance_j {
                poll(&mut cancelled)?;
                let next = Self { point: mechanical.point, curve: self.curve, state: proposal.state,
                    mass_kg: self.mass_kg, volume_m3: self.volume_m3,
                    surface_area_m2: self.surface_area_m2, time_s: final_time };
                return Ok(FeedbackStep { next, mechanical_energy: mechanical.energy, transport: proposal,
                    iterations: iteration, applied_heat_j, heat_defect_j: defects[1],
                    temperature_defect_k: defects[0], thermal_residual_j: defects[2], combined_residual_j: defects[3] });
            }
            temperature = finite((1.0 - config.relaxation) * temperature + config.relaxation * final_temperature)?;
        }
        poll(&mut cancelled)?;
        Err(FeedbackError::NotConverged { iterations: config.max_iterations,
            temperature_defect_k: defects[0], heat_defect_j: defects[1],
            thermal_residual_j: defects[2], combined_residual_j: defects[3] })
    }
}

fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn finite<E>(value: f64) -> Result<f64, FeedbackError<E>> {
    if value.is_finite() { Ok(value) } else { Err(FeedbackError::Invalid("unrepresentable feedback arithmetic")) }
}
fn product<E>(a: f64, b: f64) -> Result<f64, FeedbackError<E>> {
    let value = finite(a * b)?;
    if a != 0.0 && b != 0.0 && value == 0.0 {
        return Err(FeedbackError::Invalid("nonzero energy or mass underflowed"));
    }
    Ok(value)
}
fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), FeedbackError<E>> {
    if cancelled() { Err(FeedbackError::Cancelled) } else { Ok(()) }
}
