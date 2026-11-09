//! Original-enthalpy homogeneous-equilibrium converging-nozzle flow.
//!
//! The fixed overall composition follows its original PS path. Its source
//! enthalpy loss becomes kinetic energy, and the same equilibrating phase model
//! supplies density and sound speed. No gamma or frozen-phase Cp is substituted.
//! The admitted pressure lane must be a connected, ordinary single-sonic branch.
//! Phase-count crossings and acoustic stencil failures refuse: a bubble/dew kink
//! may choke without a smooth sonic root and needs a separate one-sided model.

use crate::nozzle::{PressureSearch, SearchError, Station};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::{
    ThermochemicalFlashError,
    equilibrium::{
        ResolvedThermochemicalPhaseState,
        acoustics::{
            PHASE_ACOUSTIC_PT_EVALUATIONS, PhaseEquilibriumAcousticControls,
            PhaseEquilibriumAcoustics,
        },
        entropy::{
            MAX_ENTROPY_FLASH_EVALUATIONS, PhaseEntropyControls, PhaseEntropySolution,
            prepared::PreparedPhaseEntropyPath,
        },
    },
};
use fs_math::det;
use fs_qty::{Area, MassFlowRate, Pressure, Temperature};

/// Bound on pressure queries, including the original stagnation state.
pub const MAX_NOZZLE_PRESSURE_EVALUATIONS: usize = crate::nozzle::MAX_PRESSURE_EVALUATIONS;
/// One preparation plus bounded PS and acoustic work for every pressure query.
pub const MAX_NOZZLE_PT_EVALUATIONS: usize = 1 + MAX_NOZZLE_PRESSURE_EVALUATIONS
    * (MAX_ENTROPY_FLASH_EVALUATIONS + PHASE_ACOUSTIC_PT_EVALUATIONS);

/// Original source policy and the caller's connected smooth expansion domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureNozzleControl {
    /// Original source bracket and entropy closure at every pressure.
    pub isentropic: PhaseEntropyControls,
    /// Native equilibrium derivative policy; phase policy must match isentropic.
    pub acoustics: PhaseEquilibriumAcousticControls,
    /// Lowest admitted pressure [Pa]. It is a search limit, not an invented exit.
    pub minimum_pressure_pa: f64,
    /// Each search pressure is this fraction of the preceding pressure.
    /// In [0.5,0.95]; sampled monotonicity is evidence, not a global root proof.
    pub pressure_step_ratio: f64,
    /// Absolute sonic residual allowance [m2/s2].
    pub absolute_sonic_tolerance_m2_s2: f64,
    /// Sonic allowance relative to the larger of velocity squared and sound speed squared.
    pub relative_sonic_tolerance: f64,
    /// Absolute width of the accepted sonic-pressure bracket [Pa].
    pub absolute_pressure_tolerance_pa: f64,
    /// Pressure-bracket width relative to the original stagnation pressure.
    pub relative_pressure_tolerance: f64,
    /// Absolute original stagnation-enthalpy balance allowance [J/kg].
    pub absolute_energy_tolerance_j_kg: f64,
    /// Energy allowance relative to actual kinetic energy, never absolute h.
    pub relative_energy_tolerance: f64,
    /// Whole-call pressure queries, including scans and root trials.
    pub maximum_pressure_evaluations: usize,
    /// Whole-call original PT work, including preparation and acoustic stencils.
    pub maximum_pt_evaluations: usize,
}
impl MixtureNozzleControl {
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

    #[must_use]
    pub fn new(temperature_bracket: [Temperature; 2], minimum_pressure: Pressure) -> Self {
        let acoustics = PhaseEquilibriumAcousticControls::default();
        let mut isentropic = PhaseEntropyControls::new(temperature_bracket);
        isentropic.flash_options = acoustics.phase.flash;
        isentropic.stability = acoustics.phase.stability;
        Self {
            isentropic,
            acoustics,
            minimum_pressure_pa: minimum_pressure.value(),
            pressure_step_ratio: 0.75,
            absolute_sonic_tolerance_m2_s2: 1e-4,
            relative_sonic_tolerance: 1e-8,
            absolute_pressure_tolerance_pa: 1e-5,
            relative_pressure_tolerance: 1e-8,
            absolute_energy_tolerance_j_kg: 1e-5,
            relative_energy_tolerance: 1e-10,
            maximum_pressure_evaluations: 128,
            maximum_pt_evaluations: MAX_NOZZLE_PT_EVALUATIONS,
        }
    }

    pub fn validate(self) -> Result<(), MixtureNozzleError> {
        self.isentropic
            .validate()
            .map_err(|e| source_error(MixtureNozzleStage::Isentropic, 0, e))?;
        self.acoustics
            .validate()
            .map_err(|e| source_error(MixtureNozzleStage::Acoustics, 0, e))?;
        if self.acoustics.phase != self.isentropic.phase_options() {
            return Err(Error::Invalid(
                "nozzle calorics and acoustics require identical phase policies",
            ));
        }
        if !positive(self.minimum_pressure_pa)
            || !(0.5..=0.95).contains(&self.pressure_step_ratio)
            || !positive(self.absolute_sonic_tolerance_m2_s2)
            || !positive(self.absolute_pressure_tolerance_pa)
            || !positive(self.absolute_energy_tolerance_j_kg)
            || [
                self.relative_sonic_tolerance,
                self.relative_pressure_tolerance,
                self.relative_energy_tolerance,
            ]
            .into_iter()
            .any(|x| !x.is_finite() || !(0.0..1.0).contains(&x))
            || !(1..=MAX_NOZZLE_PRESSURE_EVALUATIONS).contains(&self.maximum_pressure_evaluations)
            || !(1..=MAX_NOZZLE_PT_EVALUATIONS).contains(&self.maximum_pt_evaluations)
        {
            return Err(Error::Invalid(
                "finite positive nozzle bounds/tolerances and bounded work required",
            ));
        }
        Ok(())
    }
}

/// Original source stage that refused the expansion.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MixtureNozzleStage {
    Isentropic,
    Acoustics,
}

/// No partial nozzle state or inventory change is published on failure.
#[derive(Clone, Debug, PartialEq)]
pub enum MixtureNozzleError {
    Invalid(&'static str),
    Equilibrium {
        stage: MixtureNozzleStage,
        /// Successful original PT evaluations before the failing stage.
        completed_pt_evaluations: usize,
        /// The original refusal. PS work-limit errors retain attempted PT work.
        source: ThermochemicalFlashError,
    },
    WorkLimit {
        pressure_evaluations: usize,
        pt_evaluations: usize,
    },
    /// A smooth sonic root cannot be assumed across a phase-count kink.
    PhaseBoundary {
        pressure_pa: f64,
        initial_phases: usize,
        trial_phases: usize,
    },
    /// The sampled sonic residual violates the admitted ordinary branch.
    NonMonotone {
        pressure_bracket_pa: [f64; 2],
        sonic_residuals_m2_s2: [f64; 2],
    },
    /// Neither a back-pressure exit nor sonic root lies in the admitted interval.
    OutsidePressureBracket {
        minimum_pressure_pa: f64,
        back_pressure_pa: f64,
        sonic_residual_m2_s2: f64,
    },
    /// Floating-point pressure resolution was exhausted before both tolerances.
    NotConverged {
        pressure_bracket_pa: [f64; 2],
        sonic_residual_m2_s2: f64,
    },
    Cancelled,
}
impl core::fmt::Display for MixtureNozzleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "equilibrium mixture nozzle: {self:?}")
    }
}
impl std::error::Error for MixtureNozzleError {}
type Error = MixtureNozzleError;

/// One accepted steady nozzle state with the actual original energy/phase laws.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureNozzleFlow {
    pub inlet: ResolvedThermochemicalPhaseState,
    /// Isentropic exit/throat and datum-safe original enthalpy/entropy increments.
    pub outlet: PhaseEntropySolution,
    /// Native equilibrium acoustic evidence. Absent for exact zero flow.
    pub acoustics: Option<PhaseEquilibriumAcoustics>,
    /// Nonnegative flow through the supplied effective area.
    pub mass_flow: MassFlowRate,
    pub mass_flux_kg_m2_s: f64,
    pub velocity_m_s: f64,
    pub mach_number: f64,
    /// A smooth sonic limit was reached before the requested back pressure.
    pub choked: bool,
    pub back_pressure_pa: f64,
    /// Original Delta h plus the kinetic energy of the published velocity.
    pub stagnation_enthalpy_residual_j_kg: f64,
    /// Present only for a solved sonic limit, in increasing pressure order.
    pub sonic_pressure_bracket_pa: Option<[f64; 2]>,
    /// c_eq^2 - v^2; absent for exact zero flow. Acoustic uncertainty is separate.
    pub sonic_residual_m2_s2: Option<f64>,
    pub pressure_evaluations: usize,
    pub pt_evaluations: usize,
    /// Original EOS ceiling, not an observed EOS-call count.
    pub eos_evaluation_upper_bound: usize,
}

/// Resolve a one-direction steady adiabatic converging restriction.
///
/// The inlet is a stagnation state. Fixed overall composition, instantaneous
/// thermal/phase equilibrium and one common phase velocity define the HEM limit.
/// Upstream enthalpy already contains flow work: only its original integrated
/// decrease becomes kinetic energy. Effective area is explicit; no discharge
/// coefficient, friction, shaft work, heat transfer or recipient mixing is inferred.
///
/// The caller admits one connected, ordinary single-sonic pressure lane. A
/// descending scan brackets a smooth sonic root or reaches back pressure first.
/// Every sampled state must retain the inlet phase count, and each native
/// acoustic stencil must be admitted. Sampled monotonicity and finite phase
/// screening do not certify global choking/uniqueness for a nonclassical EOS.
/// Latent conversion within a two-phase lane is retained; phase-boundary kinks,
/// critical singularities, slip and finite-rate phase relaxation are unsupported.
///
/// Equal pressure or zero area returns the original inlet with zero flow after
/// one preparation and no acoustic work. Reverse pressure drive refuses; no
/// implicit check valve or opposite-direction mixture source is manufactured.
/// Vacuum back pressure is admitted without querying the EOS at zero pressure.
pub fn equilibrium_mixture_nozzle(
    model: &ThermochemicalPengRobinsonModel,
    stagnation_temperature: Temperature,
    stagnation_pressure: Pressure,
    back_pressure: Pressure,
    effective_area: Area,
    controls: MixtureNozzleControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixtureNozzleFlow, Error> {
    poll(&mut cancelled)?;
    controls.validate()?;
    let p0 = stagnation_pressure.value();
    let pb = back_pressure.value();
    let area = effective_area.value();
    let t0 = stagnation_temperature.value();
    let [tlo, thi] = controls.isentropic.temperature_bracket;
    if !positive(t0)
        || !positive(p0)
        || !pb.is_finite()
        || pb < 0.0
        || pb > p0
        || !area.is_finite()
        || area < 0.0
        || t0 < tlo.value()
        || t0 > thi.value()
    {
        return Err(Error::Invalid(
            "positive inlet pressure, nonnegative back pressure/area and admitted stagnation T required",
        ));
    }
    let moving = area > 0.0 && pb < p0;
    if moving && controls.minimum_pressure_pa >= p0 {
        return Err(Error::Invalid(
            "expansion pressure floor must lie below stagnation pressure",
        ));
    }
    let path = model
        .prepare_phase_entropy_path(
            stagnation_temperature,
            stagnation_pressure,
            controls.isentropic,
            &mut cancelled,
        )
        .map_err(|e| source_error(MixtureNozzleStage::Isentropic, 0, e))?;
    let mut work = Work {
        pressure: 0,
        pt: path.preparation_pt_evaluations(),
        eos: path.preparation_eos_evaluation_upper_bound(),
    };
    let mut context = Context {
        model,
        path: &path,
        controls,
        phases: None,
    };
    let initial = context.sample(p0, moving, &mut work, &mut cancelled)?;
    let inlet = initial.equilibrium.state().clone();
    context.phases = Some(phase_count(&inlet));
    if !moving {
        return finish(
            initial,
            inlet,
            pb,
            area,
            None,
            controls,
            work,
            &mut cancelled,
        );
    }
    let (accepted, bracket) = crate::nozzle::resolve(
        initial,
        pb,
        controls.pressure_search(),
        |pressure, cancelled| context.sample(pressure, true, &mut work, cancelled),
        &mut cancelled,
    )
    .map_err(search_error)?;
    finish(
        accepted,
        inlet,
        pb,
        area,
        bracket,
        controls,
        work,
        &mut cancelled,
    )
}

struct Work {
    pressure: usize,
    pt: usize,
    eos: usize,
}
impl Work {
    fn exhausted(&self) -> Error {
        Error::WorkLimit {
            pressure_evaluations: self.pressure,
            pt_evaluations: self.pt,
        }
    }
    fn charge(&mut self, pt: usize, eos: usize) -> Result<(), Error> {
        self.pt = self
            .pt
            .checked_add(pt)
            .ok_or(Error::Invalid("nozzle PT counter overflow"))?;
        self.eos = self
            .eos
            .checked_add(eos)
            .ok_or(Error::Invalid("nozzle EOS counter overflow"))?;
        Ok(())
    }
}
struct Sample {
    equilibrium: PhaseEntropySolution,
    acoustics: Option<PhaseEquilibriumAcoustics>,
    velocity_squared: f64,
    sonic: f64,
}
impl Station for Sample {
    fn pressure(&self) -> f64 {
        self.equilibrium.state().pressure_pa()
    }
    fn velocity_squared(&self) -> f64 {
        self.velocity_squared
    }
    fn sonic(&self) -> f64 {
        self.sonic
    }
}
struct Context<'m, 'p> {
    model: &'m ThermochemicalPengRobinsonModel,
    path: &'p PreparedPhaseEntropyPath<'m>,
    controls: MixtureNozzleControl,
    phases: Option<usize>,
}
impl Context<'_, '_> {
    fn sample(
        &self,
        pressure: f64,
        acoustic: bool,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Sample, Error> {
        poll(cancelled)?;
        if work.pressure == self.controls.maximum_pressure_evaluations {
            return Err(work.exhausted());
        }
        work.pressure += 1;
        let remaining = self.controls.maximum_pt_evaluations.saturating_sub(work.pt);
        let reserved = if acoustic {
            PHASE_ACOUSTIC_PT_EVALUATIONS
        } else {
            0
        };
        if remaining < reserved {
            return Err(work.exhausted());
        }
        let solved = self
            .path
            .solve_increment(
                Pressure::new(pressure),
                0.0,
                remaining - reserved,
                &mut *cancelled,
            )
            .map_err(|e| source_error(MixtureNozzleStage::Isentropic, work.pt, e))?;
        work.charge(solved.pt_evaluations(), solved.eos_evaluation_upper_bound())?;
        let phases = phase_count(solved.state());
        if let Some(initial) = self.phases
            && phases != initial
        {
            return Err(Error::PhaseBoundary {
                pressure_pa: pressure,
                initial_phases: initial,
                trial_phases: phases,
            });
        }
        let dh = solved.enthalpy_change_j_kg().ok_or(Error::Invalid(
            "nozzle requires original enthalpy increment",
        ))?;
        let velocity_squared = -2.0 * dh;
        if !velocity_squared.is_finite()
            || velocity_squared < 0.0
            || (self.phases.is_some() && dh >= 0.0)
        {
            return Err(Error::Invalid(
                "expansion must release positive finite original enthalpy",
            ));
        }
        let response = if acoustic {
            if self.controls.maximum_pt_evaluations.saturating_sub(work.pt)
                < PHASE_ACOUSTIC_PT_EVALUATIONS
            {
                return Err(work.exhausted());
            }
            let a = self
                .model
                .phase_equilibrium_acoustics(
                    Temperature::new(solved.state().temperature_k()),
                    Pressure::new(pressure),
                    self.controls.acoustics,
                    &mut *cancelled,
                )
                .map_err(|e| source_error(MixtureNozzleStage::Acoustics, work.pt, e))?;
            work.charge(a.pt_evaluations(), a.eos_evaluation_upper_bound())?;
            if a.source() != solved.state() {
                return Err(Error::Invalid(
                    "nozzle acoustics must retain the exact caloric source state",
                ));
            }
            Some(a)
        } else {
            None
        };
        let sonic = response.as_ref().map_or(0.0, |a| {
            let c = a.state().sound_speed_m_s();
            c * c - velocity_squared
        });
        if !sonic.is_finite() {
            return Err(Error::Invalid("unrepresentable nozzle sonic residual"));
        }
        Ok(Sample {
            equilibrium: solved,
            acoustics: response,
            velocity_squared,
            sonic,
        })
    }
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

#[allow(clippy::too_many_arguments)]
fn finish(
    sample: Sample,
    inlet: ResolvedThermochemicalPhaseState,
    back: f64,
    area: f64,
    bracket: Option<[f64; 2]>,
    controls: MixtureNozzleControl,
    work: Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<MixtureNozzleFlow, Error> {
    let velocity = det::sqrt(sample.velocity_squared);
    let flux = sample.equilibrium.state().density_kg_m3() * velocity;
    let mass_flow = area * flux;
    let dh = sample
        .equilibrium
        .enthalpy_change_j_kg()
        .ok_or(Error::Invalid(
            "nozzle requires original enthalpy increment",
        ))?;
    let residual = dh + 0.5 * velocity * velocity;
    let tolerance =
        controls.absolute_energy_tolerance_j_kg + controls.relative_energy_tolerance * dh.abs();
    let mach = sample
        .acoustics
        .as_ref()
        .map_or(0.0, |a| velocity / a.state().sound_speed_m_s());
    if !flux.is_finite()
        || !mass_flow.is_finite()
        || !residual.is_finite()
        || !mach.is_finite()
        || residual.abs() > tolerance
        || (velocity > 0.0 && area > 0.0 && mass_flow == 0.0)
    {
        return Err(Error::Invalid(
            "unrepresentable nozzle flow or original stagnation-enthalpy closure",
        ));
    }
    poll(cancelled)?;
    Ok(MixtureNozzleFlow {
        inlet,
        outlet: sample.equilibrium,
        sonic_residual_m2_s2: sample.acoustics.as_ref().map(|_| sample.sonic),
        acoustics: sample.acoustics,
        mass_flow: MassFlowRate::new(mass_flow),
        mass_flux_kg_m2_s: flux,
        velocity_m_s: velocity,
        mach_number: mach,
        choked: bracket.is_some(),
        back_pressure_pa: back,
        stagnation_enthalpy_residual_j_kg: residual,
        sonic_pressure_bracket_pa: bracket,
        pressure_evaluations: work.pressure,
        pt_evaluations: work.pt,
        eos_evaluation_upper_bound: work.eos,
    })
}
fn phase_count(state: &ResolvedThermochemicalPhaseState) -> usize {
    if state.phase().two_phase().is_some() {
        2
    } else {
        1
    }
}
fn source_error(
    stage: MixtureNozzleStage,
    completed: usize,
    source: ThermochemicalFlashError,
) -> Error {
    match source {
        ThermochemicalFlashError::Cancelled => Error::Cancelled,
        source => Error::Equilibrium {
            stage,
            completed_pt_evaluations: completed,
            source,
        },
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
