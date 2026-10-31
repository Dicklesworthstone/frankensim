//! Conserved-density/internal-energy closure with pure-fluid phase equilibrium.
//!
//! Each subcritical temperature uses the original PR saturation roots. The
//! convex volume closure admits a liquid, coexistence mixture or vapor before
//! evaluating calorics. A mechanically stable metastable density is therefore
//! not silently retained inside the dome. Solids and interfaces are excluded.
//! NASA-9 supplies the ideal reference once; PR supplies departures. No gas
//! transport or frozen-quality heat capacity is substituted for equilibrium.

/// Constant-mass, fixed-volume equilibrium heating with atomic publication.
pub mod parcel;
/// Explicit homogeneous-equilibrium pressure, entropy and acoustic response.
pub mod acoustics;

use fs_qty::{Density, Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::eos::IdealGasMolarVolumeQuantityV1 as MolarVolume;
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseRootV1;
use fs_thermochem::eos::peng_robinson::saturation::{
    PengRobinsonSaturatedMixtureV1, PengRobinsonSaturationErrorV1, PengRobinsonSaturationOptionsV1,
};
use fs_thermochem::mixture::FrozenIdealGasMixtureEvaluationV1;

use super::ThermochemicalPengRobinsonModel;
use crate::MaterialError;
use crate::thermochemical_gas::{failure, positive};

/// Explicit supported temperature/pressure envelope and whole-call work caps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PureFluidEquilibriumControls {
    /// Subcritical temperature bracket in the original NASA source domain.
    pub temperature_bracket: [Temperature; 2],
    /// Must bracket saturation throughout the temperature interval. The actual
    /// homogeneous pressure outside the dome need not lie in this bracket.
    pub pressure_bracket: [Pressure; 2],
    /// Absolute physical energy-balance tolerance [J/kg].
    pub energy_tolerance_j_kg: f64,
    /// Final Newton temperature correction or contracted bracket width [K].
    /// A short bracket never replaces energy-residual convergence.
    pub temperature_tolerance_k: f64,
    /// Whole-call caloric samples, including endpoints and increment reference.
    pub maximum_temperature_evaluations: usize,
    /// Whole-call EOS calls, including every inner coexistence iteration and
    /// each homogeneous density evaluation; no repeated per-search reset.
    pub maximum_eos_evaluations: usize,
    /// Original phase solver controls; its cap is further limited by remaining work.
    pub saturation_options: PengRobinsonSaturationOptionsV1,
}
impl PureFluidEquilibriumControls {
    /// Numerical defaults require explicit physical temperature/pressure bounds.
    #[must_use]
    pub fn new(temperature_bracket: [Temperature; 2], pressure_bracket: [Pressure; 2]) -> Self {
        Self {
            temperature_bracket,
            pressure_bracket,
            energy_tolerance_j_kg: 1e-6,
            temperature_tolerance_k: 1e-9,
            maximum_temperature_evaluations: 128,
            maximum_eos_evaluations: 8192,
            saturation_options: PengRobinsonSaturationOptionsV1 {
                log_fugacity_tolerance: 1e-12,
                relative_pressure_tolerance: 1e-12,
                ..PengRobinsonSaturationOptionsV1::default()
            },
        }
    }
    /// Validate numerical controls. Model critical temperature and original
    /// source coverage are admitted separately by the evaluating model.
    pub fn validate(self) -> Result<(), PureFluidEquilibriumError> {
        let [lo, hi] = self.temperature_bracket.map(|t| t.value());
        let [pl, pu] = self.pressure_bracket.map(|p| p.value());
        let phase = self.saturation_options;
        if !positive(lo)
            || !positive(hi)
            || lo >= hi
            || !positive(pl)
            || !positive(pu)
            || pl >= pu
            || !positive(self.energy_tolerance_j_kg)
            || !positive(self.temperature_tolerance_k)
            || !(2..=256).contains(&self.maximum_temperature_evaluations)
            || !(4..=65_536).contains(&self.maximum_eos_evaluations)
            || !(2..=512).contains(&phase.maximum_evaluations)
            || [
                phase.log_fugacity_tolerance,
                phase.relative_pressure_tolerance,
            ]
            .into_iter()
            .any(|x| !positive(x) || x >= 1.0)
        {
            return Err(E::Invalid(
                "finite ordered physical brackets, positive tolerances and bounded work required",
            ));
        }
        Ok(())
    }
}

/// A failed solve publishes no candidate or storage update.
#[derive(Clone, Debug, PartialEq)]
pub enum PureFluidEquilibriumError {
    /// Invalid controls, unsupported phase domain or unrepresentable calorics.
    Invalid(&'static str),
    /// Original material or reference-thermodynamic refusal.
    Material(MaterialError),
    /// Original EOS/coexistence refusal; unresolved phases never fall back to a
    /// homogeneous metastable answer.
    Phase(PengRobinsonSaturationErrorV1),
    /// Cancellation before a constitutive evaluation or result publication.
    Cancelled,
    /// Temperature endpoints do not bracket the requested energy.
    OutsideBracket {
        lower_residual_j_kg: f64,
        upper_residual_j_kg: f64,
    },
    /// Whole-call work exhausted without publishing a partial phase state.
    WorkLimit {
        temperature_evaluations: usize,
        eos_evaluations: usize,
    },
    /// Temperature stagnation or exhaustion without energy and thermal closure.
    NotConverged {
        temperature_evaluations: usize,
        eos_evaluations: usize,
        residual_j_kg: f64,
    },
}
impl core::fmt::Display for PureFluidEquilibriumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pure-fluid density/energy equilibrium: {self:?}")
    }
}
impl std::error::Error for PureFluidEquilibriumError {}
impl From<MaterialError> for PureFluidEquilibriumError {
    fn from(value: MaterialError) -> Self {
        Self::Material(value)
    }
}
use PureFluidEquilibriumError as E;

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), E> {
    if cancelled() {
        Err(E::Cancelled)
    } else {
        Ok(())
    }
}

/// Original fluid phase admitted by comparison with the resolved saturation
/// volumes. A homogeneous root inside those volumes is never used here.
#[derive(Clone, Debug, PartialEq)]
pub enum PureFluidEquilibriumPhase {
    /// Compressed liquid outside the liquid boundary of the dome.
    Liquid(PengRobinsonPhaseRootV1),
    /// Macroscopic liquid/vapor equilibrium, including exact saturated endpoints.
    Coexistence(PengRobinsonSaturatedMixtureV1),
    /// Vapor outside the vapor boundary of the dome.
    Vapor(PengRobinsonPhaseRootV1),
}
impl PureFluidEquilibriumPhase {
    fn residual_u(&self) -> f64 {
        match self {
            Self::Coexistence(mixture) => mixture.residual_internal_energy_j_mol(),
            Self::Liquid(root) | Self::Vapor(root) => root.residual_internal_energy_j_mol(),
        }
    }
    fn residual_h(&self) -> f64 {
        match self {
            Self::Coexistence(mixture) => mixture.residual_enthalpy_j_mol(),
            Self::Liquid(root) | Self::Vapor(root) => root.residual_enthalpy_j_mol(),
        }
    }
    fn residual_s(&self, temperature: f64) -> f64 {
        match self {
            Self::Coexistence(mixture) => mixture.residual_entropy_j_mol_k(),
            Self::Liquid(root) | Self::Vapor(root) => {
                root.residual_enthalpy_j_mol() / temperature
                    - R * root.log_fugacity_coefficients()[0]
            }
        }
    }
    fn complete_cv(&self, ideal_cv: f64) -> Result<f64, E> {
        let residual = match self {
            Self::Coexistence(mixture) => {
                for root in [mixture.saturation().liquid(), mixture.saturation().vapor()] {
                    if !positive(ideal_cv + root.residual_cv_j_mol_k()) {
                        return Err(E::Invalid(
                            "both equilibrium phases require positive original complete Cv",
                        ));
                    }
                }
                mixture
                    .equilibrium_isochoric_residual_cv_j_mol_k()
                    .map_err(E::Phase)?
            }
            Self::Liquid(root) | Self::Vapor(root) => root.residual_cv_j_mol_k(),
        };
        let complete = ideal_cv + residual;
        if !positive(complete) {
            return Err(E::Invalid(
                "positive original equilibrium caloric tangent required",
            ));
        }
        Ok(complete)
    }
}

/// Complete source-bound calorics at the conserved density. Acoustic consumers
/// explicitly select the homogeneous-equilibrium view; no gas transport or
/// acoustic relaxation policy is silently attached to this base state.
#[derive(Clone, Debug, PartialEq)]
pub struct PureFluidEquilibriumState {
    reference: FrozenIdealGasMixtureEvaluationV1,
    phase: PureFluidEquilibriumPhase,
    density_kg_m3: f64,
    specific_internal_energy_j_kg: f64,
    specific_enthalpy_j_kg: f64,
    specific_entropy_j_kg_k: f64,
    thermal_tangent_j_kg_k: f64,
}
impl PureFluidEquilibriumState {
    /// Original NASA-9 evaluation and source receipt at the accepted T,p.
    #[must_use]
    pub const fn reference(&self) -> &FrozenIdealGasMixtureEvaluationV1 {
        &self.reference
    }
    /// Original phase roots or conservative coexistence mixture.
    #[must_use]
    pub const fn phase(&self) -> &PureFluidEquilibriumPhase {
        &self.phase
    }
    /// Absolute equilibrium temperature [K].
    #[must_use]
    pub fn temperature_k(&self) -> f64 {
        self.reference.temperature().value()
    }
    /// Original EOS pressure [Pa].
    #[must_use]
    pub fn pressure_pa(&self) -> f64 {
        self.reference.pressure().value()
    }
    /// Exactly retained conserved mass density [kg/m³].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
    /// Vapor mass fraction (equal to mole fraction for this pure component).
    #[must_use]
    pub fn vapor_mass_fraction(&self) -> f64 {
        match &self.phase {
            PureFluidEquilibriumPhase::Liquid(_) => 0.0,
            PureFluidEquilibriumPhase::Vapor(_) => 1.0,
            PureFluidEquilibriumPhase::Coexistence(mixture) => mixture.vapor_fraction(),
        }
    }
    /// Full original-reference internal energy [J/kg].
    #[must_use]
    pub const fn specific_internal_energy_j_kg(&self) -> f64 {
        self.specific_internal_energy_j_kg
    }
    /// Full original-reference enthalpy [J/kg].
    #[must_use]
    pub const fn specific_enthalpy_j_kg(&self) -> f64 {
        self.specific_enthalpy_j_kg
    }
    /// Full original-reference entropy [J/(kg K)].
    #[must_use]
    pub const fn specific_entropy_j_kg_k(&self) -> f64 {
        self.specific_entropy_j_kg_k
    }
}

/// Accepted conservative state and actual source-derived energy closure/work.
#[derive(Clone, Debug, PartialEq)]
pub struct PureFluidEquilibriumEnergySolution {
    state: PureFluidEquilibriumState,
    residual_j_kg: f64,
    accepted_increment_j_kg: Option<f64>,
    temperature_evaluations: usize,
    eos_evaluations: usize,
}
impl PureFluidEquilibriumEnergySolution {
    /// Accepted complete equilibrium caloric state.
    #[must_use]
    pub const fn state(&self) -> &PureFluidEquilibriumState {
        &self.state
    }
    /// Consume the solution for atomic storage publication.
    #[must_use]
    pub fn into_state(self) -> PureFluidEquilibriumState {
        self.state
    }
    /// Evaluated absolute energy or increment minus its requested value [J/kg].
    #[must_use]
    pub const fn residual_j_kg(&self) -> f64 {
        self.residual_j_kg
    }
    /// For an increment solve, the independently evaluated NASA integral plus
    /// PR departure difference, before subtracting requested heat. None for an
    /// absolute-energy solve. A large formation datum is never subtracted.
    #[must_use]
    pub const fn accepted_increment_j_kg(&self) -> Option<f64> {
        self.accepted_increment_j_kg
    }
    /// Whole-call temperature samples including endpoints/reference.
    #[must_use]
    pub const fn temperature_evaluations(&self) -> usize {
        self.temperature_evaluations
    }
    /// Whole-call EOS evaluations, including all nested phase work.
    #[must_use]
    pub const fn eos_evaluations(&self) -> usize {
        self.eos_evaluations
    }
}

#[derive(Default)]
struct Work {
    temperatures: usize,
    eos: usize,
}
impl Work {
    fn exhausted(&self) -> E {
        E::WorkLimit {
            temperature_evaluations: self.temperatures,
            eos_evaluations: self.eos,
        }
    }
}
enum Target {
    Absolute(f64),
    Increment {
        temperature: Temperature,
        residual_u: f64,
        requested: f64,
    },
}
impl Target {
    fn requested(&self) -> f64 {
        match self {
            Self::Absolute(value) => *value,
            Self::Increment { requested, .. } => *requested,
        }
    }
}

impl ThermochemicalPengRobinsonModel {
    fn admit_pure_equilibrium(
        &self,
        density: f64,
        controls: PureFluidEquilibriumControls,
    ) -> Result<(), E> {
        controls.validate()?;
        if self.mechanics().components().len() != 1
            || self.mechanics().mole_fractions() != [1.0]
            || !positive(density)
            || !positive(self.mechanics().molar_mass().value() / density)
        {
            return Err(E::Invalid(
                "pure equilibrium requires one component and positive representable density/volume",
            ));
        }
        if controls.temperature_bracket[1].value()
            >= self.mechanics().components()[0]
                .critical_temperature
                .value()
        {
            return Err(E::Invalid(
                "pure-fluid equilibrium interval must lie below the original critical temperature",
            ));
        }
        Ok(())
    }

    fn equilibrium_sample(
        &self,
        density: f64,
        temperature: Temperature,
        controls: PureFluidEquilibriumControls,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PureFluidEquilibriumState, E> {
        poll(cancelled)?;
        let t = temperature.value();
        if !positive(t)
            || t < controls.temperature_bracket[0].value()
            || t > controls.temperature_bracket[1].value()
        {
            return Err(E::Invalid(
                "temperature outside declared subcritical equilibrium interval",
            ));
        }
        if work.temperatures >= controls.maximum_temperature_evaluations
            || controls.maximum_eos_evaluations.saturating_sub(work.eos) < 2
        {
            return Err(work.exhausted());
        }
        work.temperatures += 1;
        self.source()
            .thermodynamics()
            .evaluate_properties(temperature, controls.pressure_bracket[0])
            .map_err(|error| failure("equilibrium reference domain", error))?;
        let mut phase_options = controls.saturation_options;
        phase_options.maximum_evaluations = phase_options
            .maximum_evaluations
            .min(controls.maximum_eos_evaluations - work.eos);
        let saturation = match self.mechanics().pure_saturation_at_temperature_cancellable(
            temperature,
            controls.pressure_bracket,
            phase_options,
            &mut *cancelled,
        ) {
            Ok(state) => {
                work.eos += state.evaluations();
                state
            }
            Err(PengRobinsonSaturationErrorV1::Cancelled) => return Err(E::Cancelled),
            Err(error @ PengRobinsonSaturationErrorV1::NoConvergence { .. }) => {
                if let PengRobinsonSaturationErrorV1::NoConvergence { evaluations, .. } = &error {
                    work.eos += evaluations;
                }
                if work.eos >= controls.maximum_eos_evaluations {
                    return Err(work.exhausted());
                }
                return Err(E::Phase(error));
            }
            Err(error) => return Err(E::Phase(error)),
        };
        let molar_mass = self.mechanics().molar_mass().value();
        let volume = molar_mass / density;
        let vl = saturation.liquid().molar_volume_m3_per_mol();
        let vv = saturation.vapor().molar_volume_m3_per_mol();
        let (phase, pressure) = if volume >= vl && volume <= vv {
            let pressure = saturation.pressure();
            (
                PureFluidEquilibriumPhase::Coexistence(
                    saturation
                        .at_molar_volume(MolarVolume::new(volume))
                        .map_err(E::Phase)?,
                ),
                pressure,
            )
        } else {
            poll(cancelled)?;
            if work.eos >= controls.maximum_eos_evaluations {
                return Err(work.exhausted());
            }
            work.eos += 1;
            let homogeneous = self
                .mechanics()
                .evaluate_density_temperature(temperature, Density::new(density))
                .map_err(|error| E::Phase(error.into()))?;
            let root = homogeneous.root().clone();
            let phase = if volume < vl {
                PureFluidEquilibriumPhase::Liquid(root)
            } else {
                PureFluidEquilibriumPhase::Vapor(root)
            };
            (phase, homogeneous.pressure())
        };
        let reference = self
            .source()
            .thermodynamics()
            .evaluate(temperature, pressure)
            .map_err(|error| failure("equilibrium original reference", error))?;
        let properties = reference.properties();
        let u = (properties.u().value() + phase.residual_u()) / molar_mass;
        let h = (properties.h().value() + phase.residual_h()) / molar_mass;
        let s = (properties.s().value() + phase.residual_s(t)) / molar_mass;
        let tangent = phase.complete_cv(properties.cv().value())? / molar_mass;
        if ![u, h, s].into_iter().all(f64::is_finite) || !positive(tangent) {
            return Err(E::Invalid(
                "source-bound equilibrium calorics are unrepresentable",
            ));
        }
        Ok(PureFluidEquilibriumState {
            reference,
            phase,
            density_kg_m3: density,
            specific_internal_energy_j_kg: u,
            specific_enthalpy_j_kg: h,
            specific_entropy_j_kg_k: s,
            thermal_tangent_j_kg_k: tangent,
        })
    }

    /// Evaluate original fluid-phase equilibrium at retained density and T.
    /// Saturation must resolve even for a homogeneous endpoint: an unresolved
    /// phase comparison never becomes silent metastability or extrapolation.
    pub fn evaluate_pure_equilibrium_density_temperature(
        &self,
        density_kg_m3: f64,
        temperature: Temperature,
        controls: PureFluidEquilibriumControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PureFluidEquilibriumState, E> {
        self.admit_pure_equilibrium(density_kg_m3, controls)?;
        let state = self.equilibrium_sample(
            density_kg_m3,
            temperature,
            controls,
            &mut Work::default(),
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(state)
    }

    /// Recover T, p and phase amount from retained density and total specific
    /// internal energy. Kinetic/potential energy must be removed by its owner.
    /// The original formation/reference energy is included exactly once.
    /// Absolute-energy spacing too coarse for either requested tolerance
    /// refuses; use the source-anchored increment method for large datums.
    pub fn reconstruct_pure_density_energy(
        &self,
        density_kg_m3: f64,
        internal_energy_j_kg: f64,
        controls: PureFluidEquilibriumControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PureFluidEquilibriumEnergySolution, E> {
        self.admit_pure_equilibrium(density_kg_m3, controls)?;
        if !internal_energy_j_kg.is_finite() {
            return Err(E::Invalid("finite absolute internal energy required"));
        }
        self.solve_pure_energy(
            density_kg_m3,
            Target::Absolute(internal_energy_j_kg),
            controls,
            Work::default(),
            &mut cancelled,
        )
    }

    /// Recover equilibrium after signed fixed-volume heating/cooling. The
    /// original NASA energy increment is integrated before adding the change
    /// of PR departures; large formation energies are never subtracted.
    /// This traverses the dome boundary without imposing a homogeneous branch.
    pub fn reconstruct_pure_density_internal_energy_increment(
        &self,
        density_kg_m3: f64,
        initial_temperature: Temperature,
        change_j_kg: f64,
        controls: PureFluidEquilibriumControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PureFluidEquilibriumEnergySolution, E> {
        self.admit_pure_equilibrium(density_kg_m3, controls)?;
        if !change_j_kg.is_finite() {
            return Err(E::Invalid("finite internal-energy increment required"));
        }
        let mut work = Work::default();
        let initial = self.equilibrium_sample(
            density_kg_m3,
            initial_temperature,
            controls,
            &mut work,
            &mut cancelled,
        )?;
        if change_j_kg == 0.0 {
            poll(&mut cancelled)?;
            return Ok(PureFluidEquilibriumEnergySolution {
                state: initial,
                residual_j_kg: 0.0,
                accepted_increment_j_kg: Some(0.0),
                temperature_evaluations: work.temperatures,
                eos_evaluations: work.eos,
            });
        }
        let target = Target::Increment {
            temperature: initial_temperature,
            residual_u: initial.phase.residual_u(),
            requested: change_j_kg,
        };
        self.solve_pure_energy(density_kg_m3, target, controls, work, &mut cancelled)
    }

    fn evaluated_energy(
        &self,
        state: &PureFluidEquilibriumState,
        target: &Target,
        controls: PureFluidEquilibriumControls,
    ) -> Result<f64, E> {
        let value = match target {
            Target::Absolute(requested) => {
                let value = state.specific_internal_energy_j_kg();
                let spacing = |x: f64| (x.next_up() - x).max(x - x.next_down());
                let resolution = spacing(value).max(spacing(*requested));
                if resolution > controls.energy_tolerance_j_kg
                    || resolution / state.thermal_tangent_j_kg_k > controls.temperature_tolerance_k
                {
                    return Err(E::Invalid(
                        "absolute energy precision is unresolved; use a source-anchored energy increment",
                    ));
                }
                value
            }
            Target::Increment {
                temperature,
                residual_u,
                ..
            } => {
                let ideal = self
                    .source()
                    .thermodynamics()
                    .specific_internal_energy_change(*temperature, state.reference.temperature())
                    .map_err(|error| failure("equilibrium reference energy increment", error))?
                    .value();
                ideal
                    + (state.phase.residual_u() - residual_u)
                        / self.mechanics().molar_mass().value()
            }
        };
        if !value.is_finite() {
            return Err(E::Invalid("equilibrium energy balance is unrepresentable"));
        }
        Ok(value)
    }

    fn solve_pure_energy(
        &self,
        density: f64,
        target: Target,
        controls: PureFluidEquilibriumControls,
        mut work: Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PureFluidEquilibriumEnergySolution, E> {
        let requested = target.requested();
        let energy_tolerance = match target {
            Target::Absolute(_) => controls.energy_tolerance_j_kg,
            Target::Increment { .. } => controls.energy_tolerance_j_kg.min(0.5 * requested.abs()),
        };
        if !positive(energy_tolerance) {
            return Err(E::Invalid(
                "energy increment is below represented solve resolution",
            ));
        }
        let mut lower = self.equilibrium_sample(
            density,
            controls.temperature_bracket[0],
            controls,
            &mut work,
            cancelled,
        )?;
        let mut upper = self.equilibrium_sample(
            density,
            controls.temperature_bracket[1],
            controls,
            &mut work,
            cancelled,
        )?;
        let mut lower_value = self.evaluated_energy(&lower, &target, controls)?;
        let mut upper_value = self.evaluated_energy(&upper, &target, controls)?;
        if lower_value >= upper_value {
            return Err(E::Invalid(
                "represented equilibrium energy must increase across the admitted temperature interval",
            ));
        }
        let closed = |state: &PureFluidEquilibriumState, value: f64, width: f64| {
            let residual = value - requested;
            residual.abs() <= energy_tolerance
                && (residual.abs() / state.thermal_tangent_j_kg_k
                    <= controls.temperature_tolerance_k
                    || width <= controls.temperature_tolerance_k)
        };
        let publish = |state, value: f64, work: &Work| PureFluidEquilibriumEnergySolution {
            state,
            residual_j_kg: value - requested,
            accepted_increment_j_kg: if matches!(target, Target::Increment { .. }) {
                Some(value)
            } else {
                None
            },
            temperature_evaluations: work.temperatures,
            eos_evaluations: work.eos,
        };
        let width = upper.temperature_k() - lower.temperature_k();
        if closed(&lower, lower_value, width) {
            poll(cancelled)?;
            return Ok(publish(lower, lower_value, &work));
        }
        if closed(&upper, upper_value, width) {
            poll(cancelled)?;
            return Ok(publish(upper, upper_value, &work));
        }
        if lower_value > requested || upper_value < requested {
            return Err(E::OutsideBracket {
                lower_residual_j_kg: lower_value - requested,
                upper_residual_j_kg: upper_value - requested,
            });
        }
        let mut residual = lower_value - requested;
        let mut interior_iterations = 0;
        while work.temperatures < controls.maximum_temperature_evaluations {
            let lo = lower.temperature_k();
            let hi = upper.temperature_k();
            let (best, value) =
                if (lower_value - requested).abs() <= (upper_value - requested).abs() {
                    (&lower, lower_value)
                } else {
                    (&upper, upper_value)
                };
            let newton = best.temperature_k() - (value - requested) / best.thermal_tangent_j_kg_k;
            interior_iterations += 1;
            // Keep near-endpoint Newton corrections for smooth EOS convergence;
            // periodic bisection contracts the bracket even across a phase kink.
            let t =
                if interior_iterations % 4 != 0 && newton.is_finite() && newton > lo && newton < hi
                {
                    newton
                } else {
                    lo + 0.5 * (hi - lo)
                };
            if t <= lo || t >= hi {
                break;
            }
            let sample = self.equilibrium_sample(
                density,
                Temperature::new(t),
                controls,
                &mut work,
                cancelled,
            )?;
            let value = self.evaluated_energy(&sample, &target, controls)?;
            residual = value - requested;
            if closed(&sample, value, hi - lo) {
                poll(cancelled)?;
                return Ok(publish(sample, value, &work));
            }
            if residual < 0.0 {
                lower = sample;
                lower_value = value;
            } else {
                upper = sample;
                upper_value = value;
            }
        }
        Err(E::NotConverged {
            temperature_evaluations: work.temperatures,
            eos_evaluations: work.eos,
            residual_j_kg: residual,
        })
    }
}
