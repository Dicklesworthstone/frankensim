//! Conservative electron/heavy energy evolution at frozen populations and V.
//!
//! Independent external powers enter their explicitly selected pools. Elastic
//! electron-neutral exchange uses the SAME momentum-transfer cross sections as
//! material drag: `Qeh = G (Te-Th)`, where
//! `G = 3 kB V ne sum_n[n_n (me/m_n) km_n(Te)]`.
//! Banks (1966), equations 16–20, derives this stationary-heavy Maxwellian
//! closure: <https://ntrs.nasa.gov/api/citations/19660020172/downloads/19660020172.pdf>.
//! No electron-ion, inelastic, radiative or field power is invented. A user
//! may supply separately calculated electron/ion Joule power to these pools.
//!
//! Each substep lags G at its actual starting Te, then solves one monotone
//! implicit exchange with equal and opposite pool increments. Full/two-half
//! step comparison controls an ESTIMATED local error. All caloric inversions
//! use original-source energy DIFFERENCES from the initial state, so formation
//! energy is retained once and never subtracted as two large total energies.
//! Temperature/source limits, numerical work and the massive-target regime
//! are explicit. There is no clipping, equilibrium projection, species change,
//! two-temperature molecular reverse rate, or per-trial heap allocation.

use fs_qty::{Energy, Power, Temperature, Time};
use fs_thermochem::electron_collision::{
    AVOGADRO_PER_MOL, ELECTRON_MASS_KG, ElectronCollisionErrorV1, ElectronCollisionKindV1,
    ElectronCrossSectionV1, MaxwellianRateControlsV1, MaxwellianRateWeightingV1,
};
use fs_thermochem::equilibrium::composition::{
    FrozenIonizedGasStateErrorV1, FrozenIonizedGasStateV1, MAX_EQUILIBRIUM_SPECIES_V1,
};
use fs_thermochem::kinetics::{KineticsErrorV1, admit_caloric_interval};
use fs_thermochem::{SpeciesId, ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

mod solve;
pub(super) use solve::source_pool_calorics;
use solve::{Point, Problem};

/// Explicit physical domain and numerical error/work controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectronHeavyEnergySolve {
    /// Common minimum temperature admitted for every populated original source.
    pub minimum_temperature_k: f64,
    /// Common maximum temperature admitted for every populated original source.
    pub maximum_temperature_k: f64,
    /// Maximum charged-particle mole fraction, strictly between zero and one.
    /// This limits the chosen electron-neutral closure; it does not prove that
    /// omitted Coulomb collisions are negligible in a particular experiment.
    pub maximum_charged_mole_fraction: f64,
    /// Joint bound on `sqrt(me/m_target)` and
    /// `sqrt((Th/m_target)/(Te/me))`, strictly between zero and one. Both
    /// small recoil and slow target thermal motion are required by the
    /// stationary massive-target expansion; this is not an error certificate.
    pub maximum_target_thermal_speed_ratio: f64,
    /// Original-population charge admission at final state publication.
    pub relative_charge_tolerance: f64,
    /// Shared cross-section integration controls, including tail admission.
    pub rate_controls: MaxwellianRateControlsV1,
    /// Absolute per-pool source energy closure tolerance [J].
    pub energy_tolerance_j: f64,
    /// Absolute temperature scale in the estimated local error norm [K].
    pub temperature_tolerance_k: f64,
    /// Relative temperature and sensible-energy local error estimate tolerance.
    pub relative_tolerance: f64,
    /// First attempted step [s].
    pub initial_step_s: f64,
    /// Smallest retry step [s]; a shorter final remainder is permitted.
    pub minimum_step_s: f64,
    /// Largest attempted step [s].
    pub maximum_step_s: f64,
    /// Bound on accepted plus rejected steps.
    pub maximum_steps: usize,
    /// Bound on each monotone caloric/exchange solve.
    pub maximum_iterations: usize,
}
impl ElectronHeavyEnergySolve {
    /// Require explicit physical/source choices and supply editable numerical
    /// defaults. These defaults are integration controls, not measured accuracy.
    pub fn new(
        minimum: Temperature,
        maximum: Temperature,
        maximum_charged_mole_fraction: f64,
        maximum_target_thermal_speed_ratio: f64,
        rate_controls: MaxwellianRateControlsV1,
    ) -> Self {
        Self {
            minimum_temperature_k: minimum.value(),
            maximum_temperature_k: maximum.value(),
            maximum_charged_mole_fraction,
            maximum_target_thermal_speed_ratio,
            relative_charge_tolerance: 1e-10,
            rate_controls,
            energy_tolerance_j: 1e-9,
            temperature_tolerance_k: 1e-6,
            relative_tolerance: 1e-5,
            initial_step_s: 1e-6,
            minimum_step_s: 1e-18,
            maximum_step_s: 1e-3,
            maximum_steps: 100_000,
            maximum_iterations: 64,
        }
    }
    fn validate(self) -> Result<(), ElectronHeavyEnergyError> {
        if [
            self.minimum_temperature_k,
            self.maximum_temperature_k,
            self.energy_tolerance_j,
            self.temperature_tolerance_k,
            self.relative_tolerance,
            self.initial_step_s,
            self.minimum_step_s,
            self.maximum_step_s,
            self.maximum_target_thermal_speed_ratio,
            self.maximum_charged_mole_fraction,
            self.relative_charge_tolerance,
        ]
        .iter()
        .any(|v| !positive(*v))
            || self.minimum_temperature_k >= self.maximum_temperature_k
            || self.minimum_step_s > self.initial_step_s
            || self.initial_step_s > self.maximum_step_s
            || self.relative_tolerance >= 1.0
            || self.maximum_target_thermal_speed_ratio >= 1.0
            || self.maximum_charged_mole_fraction >= 1.0
            || self.relative_charge_tolerance >= 1.0
            || !(1..=1_000_000).contains(&self.maximum_steps)
            || !(4..=256).contains(&self.maximum_iterations)
        {
            return Err(ElectronHeavyEnergyError::Invalid(
                "energy integration controls",
            ));
        }
        Ok(())
    }
}

/// Atomic refusal; the caller's original immutable state is never changed.
#[derive(Clone, Debug, PartialEq)]
pub enum ElectronHeavyEnergyError {
    /// Invalid state, powers, source binding or physical/numerical controls.
    Invalid(&'static str),
    /// No exactly matching momentum source for a populated neutral target.
    MissingSource(SpeciesId),
    /// Original NASA-9 source refused the declared interval.
    Caloric(KineticsErrorV1),
    /// Original NASA-9 source refused a trial state/change.
    Species {
        species: SpeciesId,
        source: ThermochemErrorV1,
    },
    /// Original collision source, quadrature or tail could not support the rate.
    Collision {
        species: SpeciesId,
        source: ElectronCollisionErrorV1,
    },
    /// A temperature/energy or stationary-heavy trial exceeds the declared domain.
    OutsideDomain,
    /// A bounded monotone nonlinear solve failed to converge.
    Nonconvergence,
    /// The explicit time-step/work bound was exhausted.
    StepBudget,
    /// A needed arithmetic result is not finite/representable.
    Arithmetic(&'static str),
    /// A per-pool energy or exchange-sign audit exceeded its explicit tolerance.
    Balance(&'static str),
    /// Original frozen-state publication refused.
    State(FrozenIonizedGasStateErrorV1),
    /// Cancellation before final state publication.
    Cancelled,
}
impl core::fmt::Display for ElectronHeavyEnergyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "electron/heavy energy: {self:?}")
    }
}
impl std::error::Error for ElectronHeavyEnergyError {}

/// Actual per-pool balance and bounded-integration diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectronHeavyEnergyBalance {
    /// Independent external energy supplied to electrons [J].
    pub electron_deposited_energy_j: f64,
    /// Independent external energy supplied to heavy particles [J].
    pub heavy_deposited_energy_j: f64,
    /// Signed cumulative transfer from electrons to heavy particles [J].
    pub electron_to_heavy_energy_j: f64,
    /// Actual source electron energy change minus its target [J].
    pub electron_energy_residual_j: f64,
    /// Actual source heavy energy change minus its target [J].
    pub heavy_energy_residual_j: f64,
    /// Largest accepted full/two-half step error norm; an estimate, not proof.
    pub maximum_normalized_local_error: f64,
    /// Accepted two-half steps.
    pub accepted_steps: usize,
    /// Rejected full/two-half attempts.
    pub rejected_steps: usize,
}

/// Final source-bound state with original collision sources and energy account.
#[derive(Clone, Debug, PartialEq)]
pub struct ElectronHeavyEnergyResult {
    initial: FrozenIonizedGasStateV1,
    state: FrozenIonizedGasStateV1,
    sources: Vec<ElectronCrossSectionV1>,
    controls: ElectronHeavyEnergySolve,
    elapsed: Time,
    balance: ElectronHeavyEnergyBalance,
}
impl ElectronHeavyEnergyResult {
    /// Final actual populations, Th, Te, pressure and original species calorics.
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 {
        &self.state
    }
    /// Complete original frozen state, including the original energy reference.
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial
    }
    /// Original momentum cross sections, retaining source provenance/envelopes.
    pub fn momentum_sources(&self) -> &[ElectronCrossSectionV1] {
        &self.sources
    }
    /// Actual elapsed time [s].
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed
    }
    /// Original domain and numerical controls.
    pub const fn controls(&self) -> ElectronHeavyEnergySolve {
        self.controls
    }
    /// Per-pool source energy balance and step diagnostics.
    pub const fn balance(&self) -> ElectronHeavyEnergyBalance {
        self.balance
    }
    /// Signed cumulative elastic electron-to-heavy energy transfer [J].
    pub const fn exchanged_energy(&self) -> Energy {
        Energy::new(self.balance.electron_to_heavy_energy_j)
    }
}

/// Advance a fixed-composition, fixed-volume parcel with constant external
/// electron/heavy powers over one explicit interval. Positive power deposits
/// energy, negative power removes it. Every populated neutral target requires
/// one elastic momentum source when electrons exist; extra/duplicate sources
/// refuse. An absent electron population cannot receive electron power.
pub fn advance_frozen_energy(
    initial: &FrozenIonizedGasStateV1,
    momentum_sources: &[ElectronCrossSectionV1],
    electron_power: Power,
    heavy_power: Power,
    duration: Time,
    controls: ElectronHeavyEnergySolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ElectronHeavyEnergyResult, ElectronHeavyEnergyError> {
    controls.validate()?;
    poll(&mut cancelled)?;
    if !duration.value().is_finite()
        || duration.value() < 0.0
        || !electron_power.value().is_finite()
        || !heavy_power.value().is_finite()
    {
        return Err(ElectronHeavyEnergyError::Invalid(
            "finite powers and nonnegative duration required",
        ));
    }
    let problem = Problem::new(
        initial,
        momentum_sources,
        [heavy_power.value(), electron_power.value()],
        controls,
        &mut cancelled,
    )?;
    let mut point = Point {
        time: 0.0,
        exchange: 0.0,
        temperatures: [
            initial.heavy_temperature().value(),
            initial.electron_temperature().value(),
        ],
        capacities: [
            initial.heavy_heat_capacity_j_k(),
            initial.electron_heat_capacity_j_k(),
        ],
    };
    let mut accepted = 0;
    let mut rejected = 0;
    let mut maximum_error = 0.0_f64;
    let mut step = controls.initial_step_s;
    while point.time < duration.value() {
        poll(&mut cancelled)?;
        if accepted + rejected >= controls.maximum_steps {
            return Err(ElectronHeavyEnergyError::StepBudget);
        }
        let end = (point.time + step.min(controls.maximum_step_s)).min(duration.value());
        let middle = point.time + 0.5 * (end - point.time);
        if middle <= point.time || middle >= end {
            return Err(ElectronHeavyEnergyError::StepBudget);
        }
        let coupling = problem.coupling(&point, &mut cancelled)?;
        let trial = (|| {
            let full = problem.step(&point, end, coupling, &mut cancelled)?;
            let half = problem.step(&point, middle, coupling, &mut cancelled)?;
            let coupling_half = problem.coupling(&half, &mut cancelled)?;
            let two = problem.step(&half, end, coupling_half, &mut cancelled)?;
            Ok::<_, ElectronHeavyEnergyError>((full, two))
        })();
        let (full, two) = match trial {
            Ok(result) => result,
            Err(
                ElectronHeavyEnergyError::OutsideDomain | ElectronHeavyEnergyError::Nonconvergence,
            ) => {
                rejected += 1;
                step = 0.5 * (end - point.time);
                if step < controls.minimum_step_s {
                    return Err(ElectronHeavyEnergyError::StepBudget);
                }
                continue;
            }
            Err(error) => return Err(error),
        };
        let mut error = 0.0_f64;
        for pool in 0..2 {
            if two.capacities[pool] == 0.0 {
                continue;
            }
            let scale = controls.temperature_tolerance_k
                + controls.relative_tolerance
                    * two.temperatures[pool]
                        .abs()
                        .max(point.temperatures[pool].abs());
            error = error.max((two.temperatures[pool] - full.temperatures[pool]).abs() / scale);
        }
        let thermal_scale = (point.capacities[0] * point.temperatures[0])
            .min(point.capacities[1] * point.temperatures[1]);
        error = error.max(
            (two.exchange - full.exchange).abs()
                / (controls.energy_tolerance_j + controls.relative_tolerance * thermal_scale),
        );
        if !error.is_finite() {
            return Err(ElectronHeavyEnergyError::Arithmetic(
                "local energy error estimate",
            ));
        }
        if error > 1.0 {
            rejected += 1;
            step = (end - point.time) * (0.8 / fs_math::det::sqrt(error)).clamp(0.1, 0.8);
            if step < controls.minimum_step_s {
                return Err(ElectronHeavyEnergyError::StepBudget);
            }
            continue;
        }
        maximum_error = maximum_error.max(error);
        let accepted_step = end - point.time;
        point = two;
        accepted += 1;
        let factor = if error == 0.0 {
            2.0
        } else {
            (0.9 / fs_math::det::sqrt(error)).clamp(0.5, 2.0)
        };
        step = (accepted_step * factor)
            .min(controls.maximum_step_s)
            .max(controls.minimum_step_s);
    }
    problem.coupling(&point, &mut cancelled)?;
    let targets = problem.targets(point.time, point.exchange)?;
    let mut residuals = [0.0; 2];
    for pool in 0..2 {
        residuals[pool] = problem
            .calorics(pool, point.temperatures[pool], &mut cancelled)?
            .0
            - targets[pool];
        if !residuals[pool].is_finite() || residuals[pool].abs() > controls.energy_tolerance_j {
            return Err(ElectronHeavyEnergyError::Balance(
                "final per-pool source energy",
            ));
        }
    }
    let electron_temperature = if initial.electron_index().is_none() {
        point.temperatures[0]
    } else {
        point.temperatures[1]
    };
    let state = FrozenIonizedGasStateV1::new(
        initial.thermodynamics(),
        initial.amounts(),
        initial.volume(),
        Temperature::new(point.temperatures[0]),
        Temperature::new(electron_temperature),
        controls.relative_charge_tolerance,
        &mut cancelled,
    )
    .map_err(|error| match error {
        FrozenIonizedGasStateErrorV1::Cancelled => ElectronHeavyEnergyError::Cancelled,
        error => ElectronHeavyEnergyError::State(error),
    })?;
    poll(&mut cancelled)?;
    Ok(ElectronHeavyEnergyResult {
        initial: initial.clone(),
        state,
        sources: momentum_sources.to_vec(),
        controls,
        elapsed: duration,
        balance: ElectronHeavyEnergyBalance {
            electron_deposited_energy_j: electron_power.value() * duration.value(),
            heavy_deposited_energy_j: heavy_power.value() * duration.value(),
            electron_to_heavy_energy_j: point.exchange,
            electron_energy_residual_j: residuals[1],
            heavy_energy_residual_j: residuals[0],
            maximum_normalized_local_error: maximum_error,
            accepted_steps: accepted,
            rejected_steps: rejected,
        },
    })
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ElectronHeavyEnergyError> {
    if cancelled() {
        Err(ElectronHeavyEnergyError::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}

#[cfg(test)]
mod tests;
