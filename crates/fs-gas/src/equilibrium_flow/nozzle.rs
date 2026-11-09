//! Reacting equilibrium discharge from a stagnation reservoir.
//!
//! The original SP path supplies composition, density, direct Delta h and the
//! reacting sound speed at every pressure. Velocity follows v^2=-2 Delta h;
//! the shared converging-nozzle search selects back pressure or c_eq^2=v^2.
//! Formation energy and flow work already belong to h. Neither frozen gamma
//! nor a second reaction-heat term enters the flow law.
//!
//! This is the smooth, adiabatic, inviscid, instantaneous chemical-equilibrium
//! limit on a caller-admitted ordinary single-sonic branch. It does not infer
//! losses, a discharge coefficient, shocks, divergent geometry, finite-rate
//! chemistry or recipient mixing. Sampled branch checks are numerical evidence,
//! not a global choking/uniqueness certificate for arbitrary source laws.

use crate::nozzle::{PressureSearch, SearchError, Station};
use fs_math::{dd::Dd, det};
use fs_qty::{Amount, Area, MassFlowRate, Pressure, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::equilibrium::composition::{
    CompositionEquilibriumErrorV1, EquilibriumEntropySolveV1, EquilibriumExpansionErrorV1,
    EquilibriumThermodynamicStateV1, IdealGasEquilibriumV1,
    isentrope::{EquilibriumCaloricAnchorV1, EquilibriumIsentropeSampleV1, EquilibriumIsentropeV1},
};

/// Original material parcels admitted against this same reacting flow law.
#[cfg(feature = "material")]
pub mod material;

/// Finite back-pressure-aware depletion of the original reacting reservoir.
#[cfg(feature = "material")]
pub mod blowdown;

/// Bidirectional conservative exchange between two finite reacting volumes.
#[cfg(feature = "material")]
pub mod pair;

/// Whole-call pressure work ceiling, including preparation and re-publication.
pub const MAX_NOZZLE_PRESSURE_EVALUATIONS: usize = crate::nozzle::MAX_PRESSURE_EVALUATIONS;

/// Original source controls and finite ordinary-branch pressure search.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumNozzleControl {
    /// Original Gibbs/entropy tolerances and source-temperature/work bounds.
    /// Entropy tolerance is extensive for the supplied feed amount [J/K].
    pub isentropic: EquilibriumEntropySolveV1,
    /// Positive source-supported search floor [Pa], not an invented outlet.
    pub minimum_pressure_pa: f64,
    /// Descending scan ratio in [0.5,0.95].
    pub pressure_step_ratio: f64,
    /// Absolute allowance on c_eq^2-v^2 at a sonic throat [m2/s2].
    pub absolute_sonic_tolerance_m2_s2: f64,
    /// Additional sonic allowance relative to the larger of c_eq^2 and v^2.
    pub relative_sonic_tolerance: f64,
    /// Absolute width of the admitted sonic-pressure bracket [Pa].
    pub absolute_pressure_tolerance_pa: f64,
    /// Pressure width relative to the original stagnation pressure.
    pub relative_pressure_tolerance: f64,
    /// Absolute source Delta h + v^2/2 closure allowance [J/kg].
    pub absolute_energy_tolerance_j_kg: f64,
    /// Energy allowance relative to released enthalpy, not a formation datum.
    pub relative_energy_tolerance: f64,
    /// At most 256 pressure preparations/queries, including final publication.
    pub maximum_pressure_evaluations: usize,
}
impl EquilibriumNozzleControl {
    /// Retain explicit original source bounds while choosing numerical defaults.
    pub fn new(isentropic: EquilibriumEntropySolveV1, minimum_pressure: Pressure) -> Self {
        Self {
            isentropic,
            minimum_pressure_pa: minimum_pressure.value(),
            pressure_step_ratio: 0.75,
            absolute_sonic_tolerance_m2_s2: 1e-4,
            relative_sonic_tolerance: 1e-8,
            absolute_pressure_tolerance_pa: 1e-5,
            relative_pressure_tolerance: 1e-8,
            absolute_energy_tolerance_j_kg: 1e-5,
            relative_energy_tolerance: 1e-10,
            maximum_pressure_evaluations: 128,
        }
    }

    fn pressure_search(self) -> PressureSearch {
        PressureSearch {
            minimum_pressure_pa: self.minimum_pressure_pa,
            pressure_step_ratio: self.pressure_step_ratio,
            absolute_sonic_tolerance_m2_s2: self.absolute_sonic_tolerance_m2_s2,
            relative_sonic_tolerance: self.relative_sonic_tolerance,
            absolute_pressure_tolerance_pa: self.absolute_pressure_tolerance_pa,
            relative_pressure_tolerance: self.relative_pressure_tolerance,
        }
    }

    /// Native source controls are additionally validated by the original path.
    pub fn validate(self) -> Result<(), EquilibriumNozzleError> {
        self.pressure_search()
            .validate()
            .map_err(EquilibriumNozzleError::Invalid)?;
        if !positive(self.absolute_energy_tolerance_j_kg)
            || !self.relative_energy_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_energy_tolerance)
            || !(1..=MAX_NOZZLE_PRESSURE_EVALUATIONS).contains(&self.maximum_pressure_evaluations)
        {
            return Err(EquilibriumNozzleError::Invalid(
                "positive energy allowance and bounded nozzle source work required",
            ));
        }
        Ok(())
    }
}

/// No refused query publishes a partial outlet or consumes a physical reservoir.
#[derive(Clone, Debug, PartialEq)]
pub enum EquilibriumNozzleError {
    Invalid(&'static str),
    /// Original source, conservation, SP, derivative or numerical-work refusal.
    Chemistry {
        pressure_pa: f64,
        /// Attempted pressure preparations/queries, including this refusal.
        pressure_evaluations: usize,
        source: EquilibriumExpansionErrorV1,
    },
    WorkLimit {
        pressure_evaluations: usize,
    },
    NonMonotone {
        pressure_bracket_pa: [f64; 2],
        sonic_residuals_m2_s2: [f64; 2],
    },
    OutsidePressureBracket {
        minimum_pressure_pa: f64,
        back_pressure_pa: f64,
        sonic_residual_m2_s2: f64,
    },
    NotConverged {
        pressure_bracket_pa: [f64; 2],
        sonic_residual_m2_s2: f64,
    },
    /// A source-derived velocity, flow or actual energy balance is unresolved.
    Balance(&'static str),
    Cancelled,
}
impl core::fmt::Display for EquilibriumNozzleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reacting equilibrium converging nozzle: {self:?}")
    }
}
impl std::error::Error for EquilibriumNozzleError {}
type Error = EquilibriumNozzleError;

/// One complete chemical-equilibrium stream endpoint and its actual flow law.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumNozzleFlow {
    /// Original source-bearing equilibrium at the supplied stagnation T,p.
    pub inlet: EquilibriumThermodynamicStateV1,
    /// Original source-bearing back-pressure exit or sonic throat.
    pub outlet: EquilibriumThermodynamicStateV1,
    /// Final source-relative H/S increments and reacting derivatives. For a
    /// finite reservoir these retain the immutable tank-birth reference.
    pub sample: EquilibriumIsentropeSampleV1,
    /// Actual outlet-minus-stagnation specific enthalpy [J/kg]. This is relative
    /// to this flow's inlet, even when `sample` retains an earlier birth datum.
    pub specific_enthalpy_change_j_kg: f64,
    /// Exact original feed mass; a normalization inventory, not a flow rate.
    pub feed_mass_kg: f64,
    /// Authored effective flow area. No empirical discharge loss is inferred.
    pub effective_area: Area,
    /// Nonnegative pressure-driven mass flow through that area.
    pub mass_flow: MassFlowRate,
    pub mass_flux_kg_m2_s: f64,
    pub velocity_m_s: f64,
    /// Uses chemical-equilibrium acoustics, not frozen product Cp/Cv.
    pub mach_number: f64,
    pub choked: bool,
    pub back_pressure_pa: f64,
    /// Original Delta h plus actual kinetic energy [J/kg].
    pub stagnation_enthalpy_residual_j_kg: f64,
    /// Increasing positive pressures; present only for an admitted sonic root.
    pub sonic_pressure_bracket_pa: Option<[f64; 2]>,
    /// c_eq^2-v^2, absent for exact zero flow.
    pub sonic_residual_m2_s2: Option<f64>,
    /// Actual pressure preparations/queries, including final re-publication.
    pub pressure_evaluations: usize,
    /// Original composition iterations across preparation and every trial.
    pub composition_iterations: usize,
    /// Original temperature iterations across every SP trial, not brackets.
    pub temperature_iterations: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Trial {
    sample: EquilibriumIsentropeSampleV1,
    enthalpy_change_j_kg: f64,
    velocity_squared: f64,
    sonic: f64,
}
impl Station for Trial {
    fn pressure(&self) -> f64 {
        self.sample.pressure.value()
    }
    fn velocity_squared(&self) -> f64 {
        self.velocity_squared
    }
    fn sonic(&self) -> f64 {
        self.sonic
    }
}
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct Work {
    pressure: usize,
    composition: usize,
    temperature: usize,
}

/// Scalar accepted flow; allocating source receipts remains a separate final step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct ResolvedNozzle {
    accepted: Trial,
    bracket: Option<[f64; 2]>,
    velocity_m_s: f64,
    mass_flux_kg_m2_s: f64,
    mass_flow_kg_s: f64,
    mach_number: f64,
    energy_residual_j_kg: f64,
    work: Work,
}

/// Resolve a one-direction converging restriction from a stationary reservoir.
///
/// Feed supplies elemental ratios, not prescribed mass flow. The inlet is its
/// equilibrium at the supplied T,p. A positive back pressure is queried only if
/// it precedes the sonic limit. Vacuum still requires a supported positive sonic
/// station. Zero area or equal pressures returns the exact original inlet after
/// one preparation. Reverse pressure drive refuses; no opposite donor is invented.
/// All source work shares the original reusable SP workspace. Only the inlet and
/// final accepted outlet publish complete receipts; failed/cancelled trials do not.
pub fn equilibrium_nozzle_flow(
    chemistry: &IdealGasEquilibriumV1,
    stagnation_temperature: Temperature,
    stagnation_pressure: Pressure,
    feed: &[(SpeciesId, Amount)],
    back_pressure: Pressure,
    effective_area: Area,
    controls: EquilibriumNozzleControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<EquilibriumNozzleFlow, Error> {
    poll(&mut cancelled)?;
    let p0 = stagnation_pressure.value();
    let back = back_pressure.value();
    let area = effective_area.value();
    validate_flow_inputs(p0, back, area, controls)?;
    let mut path = EquilibriumIsentropeV1::new(
        chemistry, stagnation_temperature, stagnation_pressure, feed,
        controls.isentropic, &mut cancelled,
    ).map_err(|source| chemistry_error(p0, 1, source))?;
    let inlet = path.reference_state().clone();
    let mut stagnation = path.current_sample().expect("new isentrope has its reference");
    stagnation.composition_iterations = 0;
    let work = Work {
        pressure: 1,
        composition: inlet.equilibrium().composition_iterations(),
        temperature: 0,
    };
    let resolved = resolve_on_path(&mut path, stagnation, None, back, area, controls, work, &mut cancelled)?;
    publish_on_path(&mut path, inlet, stagnation, resolved, back, area, &mut cancelled)
}

fn validate_flow_inputs(p0: f64, back: f64, area: f64, controls: EquilibriumNozzleControl)
    -> Result<(), Error> {
    controls.validate()?;
    if !positive(p0) || !back.is_finite() || back < 0.0 || back > p0
        || !area.is_finite() || area < 0.0 {
        return Err(Error::Invalid(
            "positive stagnation pressure, nonnegative area and nonincreasing back pressure required",
        ));
    }
    let moving = area > 0.0 && back < p0;
    if moving && controls.minimum_pressure_pa >= p0 {
        return Err(Error::Invalid(
            "nozzle pressure floor must be below the stagnation pressure",
        ));
    }
    Ok(())
}

/// Reuse a genuine birth-isentrope sample supplied by the original reservoir.
/// Scalar pressure trials allocate no source receipts. The caller must keep the
/// same source/control binding that prepared this path; this seam is crate-local.
pub(crate) fn resolve_on_isentrope(
    path: &mut EquilibriumIsentropeV1<'_>,
    stagnation: EquilibriumIsentropeSampleV1,
    anchor: &EquilibriumCaloricAnchorV1,
    back: f64,
    area: f64,
    controls: EquilibriumNozzleControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ResolvedNozzle, Error> {
    resolve_on_path(path, stagnation, Some(anchor), back, area, controls, Work::default(), cancelled)
}

#[allow(clippy::too_many_arguments)]
fn resolve_on_path(
    path: &mut EquilibriumIsentropeV1<'_>,
    stagnation: EquilibriumIsentropeSampleV1,
    anchor: Option<&EquilibriumCaloricAnchorV1>,
    back: f64,
    area: f64,
    controls: EquilibriumNozzleControl,
    mut work: Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ResolvedNozzle, Error> {
    poll(cancelled)?;
    validate_flow_inputs(stagnation.pressure.value(), back, area, controls)?;
    if !stagnation.entropy_change_j_k.is_finite()
        || stagnation.entropy_change_j_k.abs() > controls.isentropic.entropy_tolerance_j_k {
        return Err(Error::Balance("stagnation state is outside the immutable source isentrope"));
    }
    if anchor.is_some_and(|anchor| anchor.sample() != stagnation)
        || (anchor.is_none() && stagnation.specific_enthalpy_change_j_kg != 0.0) {
        return Err(Error::Balance("nozzle caloric reference differs from its stagnation sample"));
    }
    let initial = trial(stagnation, 0.0, false)?;
    let moving = area > 0.0 && back < stagnation.pressure.value();
    if !moving {
        return finish_scalar(initial, area, None, controls, work, cancelled);
    }
    let mut current = initial;
    let (mut accepted, bracket) = crate::nozzle::resolve(
        initial, back, controls.pressure_search(),
        |pressure, cancelled| {
            current = query(path, pressure, anchor, controls, &mut work, cancelled)?;
            Ok(current)
        }, cancelled,
    ).map_err(search_error)?;
    // The search can select the earlier endpoint of its final bracket. Reuse
    // current work when possible; otherwise charge a genuine original SP query.
    if accepted != current {
        accepted = query(path, accepted.pressure(), anchor, controls, &mut work, cancelled)?;
    }
    if let Some(bracket) = bracket
        && accepted.sonic.abs() > controls.pressure_search().sonic_tolerance(&accepted) {
        return Err(Error::NotConverged {
            pressure_bracket_pa: bracket,
            sonic_residual_m2_s2: accepted.sonic,
        });
    }
    finish_scalar(accepted, area, bracket, controls, work, cancelled)
}

/// Publish a source-complete station only after its scalar pressure solve passed.
/// The inlet receipt is the original admitted reservoir, never a throat replay.
pub(crate) fn publish_on_path(
    path: &mut EquilibriumIsentropeV1<'_>,
    inlet: EquilibriumThermodynamicStateV1,
    stagnation: EquilibriumIsentropeSampleV1,
    resolved: ResolvedNozzle,
    back: f64,
    area: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<EquilibriumNozzleFlow, Error> {
    if inlet.equilibrium().temperature() != stagnation.temperature
        || inlet.equilibrium().pressure() != stagnation.pressure
        || inlet.response() != stagnation.response {
        return Err(Error::Balance("published inlet differs from the accepted stagnation state"));
    }
    let accepted = resolved.accepted;
    let outlet = if resolved.velocity_m_s == 0.0 {
        inlet.clone()
    } else {
        path.publish_current(&mut *cancelled)
            .map_err(|source| chemistry_error(accepted.pressure(), resolved.work.pressure, source))?
    };
    if outlet.equilibrium().temperature() != accepted.sample.temperature
        || outlet.equilibrium().pressure() != accepted.sample.pressure
        || outlet.response() != accepted.sample.response {
        return Err(Error::Balance(
            "published source state differs from the accepted nozzle trial",
        ));
    }
    poll(cancelled)?;
    Ok(EquilibriumNozzleFlow {
        inlet,
        outlet,
        sample: accepted.sample,
        specific_enthalpy_change_j_kg: accepted.enthalpy_change_j_kg,
        feed_mass_kg: path.feed_mass_kg(),
        effective_area: Area::new(area),
        mass_flow: MassFlowRate::new(resolved.mass_flow_kg_s),
        mass_flux_kg_m2_s: resolved.mass_flux_kg_m2_s,
        velocity_m_s: resolved.velocity_m_s,
        mach_number: resolved.mach_number,
        choked: resolved.bracket.is_some(),
        back_pressure_pa: back,
        stagnation_enthalpy_residual_j_kg: resolved.energy_residual_j_kg,
        sonic_pressure_bracket_pa: resolved.bracket,
        sonic_residual_m2_s2: (resolved.velocity_m_s > 0.0).then_some(accepted.sonic),
        pressure_evaluations: resolved.work.pressure,
        composition_iterations: resolved.work.composition,
        temperature_iterations: resolved.work.temperature,
    })
}

fn query(
    path: &mut EquilibriumIsentropeV1<'_>,
    pressure: f64,
    anchor: Option<&EquilibriumCaloricAnchorV1>,
    controls: EquilibriumNozzleControl,
    work: &mut Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Trial, Error> {
    poll(cancelled)?;
    if work.pressure >= controls.maximum_pressure_evaluations {
        return Err(Error::WorkLimit { pressure_evaluations: work.pressure });
    }
    work.pressure += 1;
    let sample = path.sample(Pressure::new(pressure), &mut *cancelled)
        .map_err(|source| chemistry_error(pressure, work.pressure, source))?;
    work.composition += sample.composition_iterations;
    work.temperature += sample.temperature_iterations;
    let dh = match anchor {
        Some(anchor) => path.specific_enthalpy_change_from(anchor, &mut *cancelled)
            .map_err(|source| chemistry_error(pressure, work.pressure, source))?,
        None => sample.specific_enthalpy_change_j_kg,
    };
    trial(sample, dh, true)
}

fn trial(sample: EquilibriumIsentropeSampleV1, dh: f64, moving: bool)
    -> Result<Trial, Error> {
    let velocity_squared = -2.0 * dh;
    let sound = sample.response.sound_speed_m_s;
    let sonic = sound * sound - velocity_squared;
    if !positive(sample.response.density_kg_m3) || !positive(sound)
        || !velocity_squared.is_finite() || velocity_squared < 0.0
        || !sonic.is_finite() || (moving && (dh >= 0.0 || velocity_squared == 0.0)) {
        return Err(Error::Balance(
            "expansion must resolve positive original enthalpy release and reacting acoustics",
        ));
    }
    Ok(Trial { sample, enthalpy_change_j_kg: dh, velocity_squared, sonic })
}

#[allow(clippy::too_many_arguments)]
fn finish_scalar(
    accepted: Trial,
    area: f64,
    bracket: Option<[f64; 2]>,
    controls: EquilibriumNozzleControl,
    work: Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ResolvedNozzle, Error> {
    let velocity = det::sqrt(accepted.velocity_squared);
    let flux = accepted.sample.response.density_kg_m3 * velocity;
    let rate = area * flux;
    let dh = accepted.enthalpy_change_j_kg;
    let residual = (Dd::from_f64(dh)
        + Dd::from_f64(0.5) * Dd::from_f64(velocity) * Dd::from_f64(velocity)).to_f64();
    let tolerance = controls.absolute_energy_tolerance_j_kg
        + controls.relative_energy_tolerance * dh.abs();
    let mach = velocity / accepted.sample.response.sound_speed_m_s;
    if !flux.is_finite() || !rate.is_finite() || !mach.is_finite()
        || !residual.is_finite() || !tolerance.is_finite() || residual.abs() > tolerance
        || (velocity > 0.0 && (flux == 0.0 || (area > 0.0 && rate == 0.0))) {
        return Err(Error::Balance(
            "unrepresentable nozzle flow or source stagnation-enthalpy balance",
        ));
    }
    poll(cancelled)?;
    Ok(ResolvedNozzle {
        accepted, bracket, mass_flow_kg_s: rate, mass_flux_kg_m2_s: flux,
        velocity_m_s: velocity, mach_number: mach, energy_residual_j_kg: residual, work,
    })
}

fn search_error(error: SearchError<Error>) -> Error {
    match error {
        SearchError::Invalid(reason) => Error::Invalid(reason),
        SearchError::Source(error) => error,
        SearchError::Cancelled => Error::Cancelled,
        SearchError::NonMonotone {
            pressure_bracket_pa,
            sonic_residuals_m2_s2,
        } => Error::NonMonotone {
            pressure_bracket_pa,
            sonic_residuals_m2_s2,
        },
        SearchError::OutsidePressureBracket {
            minimum_pressure_pa,
            back_pressure_pa,
            sonic_residual_m2_s2,
        } => Error::OutsidePressureBracket {
            minimum_pressure_pa,
            back_pressure_pa,
            sonic_residual_m2_s2,
        },
        SearchError::NotConverged {
            pressure_bracket_pa,
            sonic_residual_m2_s2,
        } => Error::NotConverged {
            pressure_bracket_pa,
            sonic_residual_m2_s2,
        },
    }
}
fn chemistry_error(
    pressure: f64,
    evaluations: usize,
    source: EquilibriumExpansionErrorV1,
) -> Error {
    if matches!(
        source,
        EquilibriumExpansionErrorV1::Equilibrium(CompositionEquilibriumErrorV1::Cancelled)
    ) {
        Error::Cancelled
    } else {
        Error::Chemistry {
            pressure_pa: pressure,
            pressure_evaluations: evaluations,
            source,
        }
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
