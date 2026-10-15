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

/// Work-conjugate uniaxial endpoint control, independent of the material law.
/// Stress is work per reference volume per total engineering strain, matching
/// the existing small-strain Maxwell owner; it is not a finite-strain load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MechanicalControl {
    /// Imposed total engineering strain. The material supplies the reaction.
    TotalStrain(f64),
    /// Imposed endpoint stress [Pa], positive in tension. Strain is solved from
    /// the existing affine constitutive tangent at EACH thermal trial. A changing
    /// load uses endpoint work, as in the existing backward-Euler material law.
    Stress {
        stress_pa: f64,
        /// Independent absolute endpoint stress residual gate [Pa].
        tolerance_pa: f64,
    },
}

impl MechanicalControl {
    fn validate<E>(self) -> Result<(), FeedbackError<E>> {
        match self {
            Self::TotalStrain(strain) if strain.is_finite() => Ok(()),
            Self::Stress { stress_pa, tolerance_pa }
                if stress_pa.is_finite() && positive(tolerance_pa) => Ok(()),
            _ => Err(FeedbackError::Invalid("control needs finite strain or stress and positive stress tolerance")),
        }
    }

    fn propose<'a, E>(self, old: &Point<'a>, temperature: f64, duration_s: f64,
        energy_tolerance: f64) -> Result<super::Step<'a>, FeedbackError<E>> {
        let strain = match self {
            Self::TotalStrain(strain) => strain,
            Self::Stress { stress_pa, .. } => {
                // This inversion changes the boundary condition, not the
                // constitutive update or its branch-history integration.
                let affine = old.linearize(temperature, duration_s)
                    .map_err(FeedbackError::Mechanical)?;
                finite(finite(stress_pa - affine.offset_pa)? / affine.tangent_pa)?
            }
        };
        let step = old.advance(temperature, strain, duration_s, energy_tolerance)
            .map_err(FeedbackError::Mechanical)?;
        if let Self::Stress { stress_pa, tolerance_pa } = self {
            let residual = finite(step.stress_pa - stress_pa)?;
            if residual.abs() > tolerance_pa {
                return Err(FeedbackError::StressResidual { target_pa: stress_pa,
                    actual_pa: step.stress_pa, tolerance_pa });
            }
        }
        Ok(step)
    }
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
    /// Actual constitutive stress missed the independently prescribed load gate.
    StressResidual { target_pa: f64, actual_pa: f64, tolerance_pa: f64 },
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
            Self::StressResidual { target_pa, actual_pa, tolerance_pa } => write!(f,
                "thermovisco load control: target={target_pa} Pa, actual={actual_pa} Pa, allowance={tolerance_pa} Pa"),
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
    /// The imposed boundary condition, not a fitted material parameter.
    pub control: MechanicalControl,
    /// Actual minus prescribed stress [Pa]; None for strain control.
    pub stress_residual_pa: Option<f64>,
    /// Endpoint boundary work [J]. For stress control this uses the prescribed
    /// load, not a possibly imperfect constitutive reaction. Under strain
    /// control it is the original law's mechanical work times reference volume.
    pub external_mechanical_work_j: f64,
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
    /// delta_elastic + m delta_h + numerical_loss - imposed work - heat [J].
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
        transport: impl FnMut(UniformEnthalpyStepInput<'_>) -> Result<UniformEnthalpyStep<R>, E>,
        cancelled: impl FnMut() -> bool) -> Result<FeedbackStep<'a, R>, FeedbackError<E>> {
        self.advance_controlled(MechanicalControl::TotalStrain(total_strain), duration_s,
            config, transport, cancelled)
    }

    /// Solve strain or stress control and heat feedback in the SAME iteration.
    /// Stress control enables creep, unloading and signed load reversal without
    /// resetting viscoelastic history. Every trial and the returned endpoint
    /// solve the affine stress law from the same old point. The original affine
    /// tangent is a fixed-temperature constitutive tangent, not a derivative of
    /// this coupled thermal solve. This method supplies no structural dynamics.
    ///
    /// Prescribed endpoint work enters the combined energy gate independently
    /// of the stress gate; an allowed force residual is not relabeled as heat.
    /// The original transport transaction, source-domain, work and cancellation
    /// limits apply to both control modes. A rejected proposal changes neither
    /// strain, branch memory, enthalpy nor the accepted clock.
    pub fn advance_controlled<R, E>(&self, control: MechanicalControl, duration_s: f64,
        config: FeedbackConfig,
        mut transport: impl FnMut(UniformEnthalpyStepInput<'_>) -> Result<UniformEnthalpyStep<R>, E>,
        mut cancelled: impl FnMut() -> bool) -> Result<FeedbackStep<'a, R>, FeedbackError<E>> {
        poll(&mut cancelled)?;
        control.validate()?;
        let final_time = self.time_s + duration_s;
        let material_tolerance = config.energy_tolerance_j / self.volume_m3;
        if config.max_iterations == 0 || config.max_iterations > MAX_FEEDBACK_ITERATIONS
            || !positive(config.temperature_tolerance_k) || !positive(config.energy_tolerance_j)
            || !positive(config.relaxation) || config.relaxation > 1.0
            || !positive(duration_s) || !positive(material_tolerance)
            || !final_time.is_finite() || final_time <= self.time_s {
            return Err(FeedbackError::Invalid("invalid feedback controls, strain or physical interval"));
        }
        let mut temperature = self.point.temperature_k();
        let mut defects = [f64::INFINITY; 4];
        for iteration in 1..=config.max_iterations {
            poll(&mut cancelled)?;
            let trial = control.propose(&self.point, temperature, duration_s, material_tolerance)?;
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
            let mechanical = control.propose(&self.point, final_temperature, duration_s, material_tolerance)?;
            let actual_heat = product(self.volume_m3, mechanical.energy.dissipated)?;
            let thermal_residual = finite(thermal_change - actual_heat - proposal.external_heat_j)?;
            let stress_residual_pa = match control {
                MechanicalControl::TotalStrain(_) => None,
                MechanicalControl::Stress { stress_pa, .. } => Some(finite(mechanical.stress_pa - stress_pa)?),
            };
            let work_density = match control {
                MechanicalControl::TotalStrain(_) => mechanical.energy.mechanical_work,
                MechanicalControl::Stress { stress_pa, .. } => product(stress_pa,
                    mechanical.point.total_strain() - self.point.total_strain())?,
            };
            let external_mechanical_work_j = product(self.volume_m3, work_density)?;
            let combined_residual = finite(product(self.volume_m3, mechanical.energy.stored_change)?
                + thermal_change + product(self.volume_m3, mechanical.energy.numerical_dissipation)?
                - external_mechanical_work_j - proposal.external_heat_j)?;
            defects = [(final_temperature - temperature).abs(), (actual_heat - applied_heat_j).abs(),
                thermal_residual, combined_residual];
            if defects[0] <= config.temperature_tolerance_k && defects[1] <= config.energy_tolerance_j
                && defects[2].abs() <= config.energy_tolerance_j && defects[3].abs() <= config.energy_tolerance_j {
                poll(&mut cancelled)?;
                let next = Self { point: mechanical.point, curve: self.curve, state: proposal.state,
                    mass_kg: self.mass_kg, volume_m3: self.volume_m3,
                    surface_area_m2: self.surface_area_m2, time_s: final_time };
                return Ok(FeedbackStep { next, control, stress_residual_pa, external_mechanical_work_j,
                    mechanical_energy: mechanical.energy, transport: proposal,
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
