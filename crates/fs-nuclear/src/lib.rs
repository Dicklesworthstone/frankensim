//! Conservative nuclear sources with explicit physical inventories and histories.
//!
//! The same accepted event count consumes finite fuel, releases prompt energy
//! and fills delayed-energy groups. Each group solves dE/dt = q F - lambda E
//! analytically for a constant fission rate, including exhaustion within a step.
//! Decay continues when the source is shut down or the fuel is exhausted.
//! Local heat, outward transportable radiation and nonrecoverable export remain
//! separate. A caller must supply and identify its own admitted nuclear data;
//! no isotope, energy yield, irradiation history or deposition fraction is guessed.
//!
//! The fission source is a reduced model, not neutron transport, criticality or
//! a radiation-transport solver. Its effective decay groups are energy reservoirs,
//! not delayed-neutron precursor populations. The separate [`decay`] owner evolves
//! source-identified radioactive parent/daughter populations and branch energy
//! on the shared `fs_la` propagator and `fs_thermochem` half-life laws. Thermal
//! material mass is a separate nonrelativistic storage coordinate; consumed fuel
//! nuclei must not be subtracted as lost bulk material mass.

#![deny(unsafe_code)]

use fs_math::decay::{attenuate, decayed_fraction};

/// Source-bound radioactive populations, daughter branches and integrated energy.
pub mod decay;
use std::sync::Arc;

/// Source coupling on the original material-backed thermal timeline.
#[cfg(feature = "thermal")]
pub mod thermal;

/// Multigroup neutron irradiation with competing fission and capture.
pub mod irradiation;

/// Identified, temperature-dependent microscopic neutron reaction data.
pub mod cross_sections;

/// Bounded processed ENDF/PENDF import and declared spectrum collapse.
pub mod endf;

/// Spatial multigroup neutron histories and composition-coupled depletion.
pub mod diffusion;

/// Fusion-neutron injection, kinetic-energy transport and disjoint heat/leakage.
pub mod fusion_neutrons;

/// Product-forming fission, capture and radioactive transmutation networks.
pub mod burnup;

/// Point-neutron kinetics with independent delayed-neutron precursor history.
pub mod kinetics;

/// Conservative deposition of emitted radiation using explicit transport responses.
pub mod radiation;

/// Disjoint energy destinations, in joules (or joules per fission in a model).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EnergyChannels {
    /// Energy actually thermalized in the source material.
    pub deposited_j: f64,
    /// Energy leaving this source for a separately owned radiation transport.
    pub radiation_j: f64,
    /// Nonrecoverable energy leaving the modeled system, such as neutrinos.
    pub escaped_j: f64,
}

impl EnergyChannels {
    /// Sum of the three disjoint accounts; validated model/state values are finite.
    #[must_use]
    pub fn total_j(self) -> f64 {
        sum(&[self.deposited_j, self.radiation_j, self.escaped_j])
    }

    fn validate(self) -> Result<(), FissionError> {
        if ![self.deposited_j, self.radiation_j, self.escaped_j]
            .into_iter().all(nonnegative) || !self.total_j().is_finite()
        {
            return Err(FissionError::InvalidInput("energy channels must be finite and nonnegative"));
        }
        Ok(())
    }

    fn scaled(self, factor: f64) -> Result<Self, FissionError> {
        let result = Self {
            deposited_j: finite(self.deposited_j * factor)?,
            radiation_j: finite(self.radiation_j * factor)?,
            escaped_j: finite(self.escaped_j * factor)?,
        };
        finite(result.total_j())?;
        Ok(result)
    }
}

/// One nonnegative exponential term in a source-supplied decay-energy model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayHeatGroup {
    /// Positive decay constant lambda [1/s].
    pub decay_constant_s_inv: f64,
    /// Integrated energy released by this group per fission [J/fission].
    /// For a supplied impulse-power fit alpha exp(-lambda t), this is alpha/lambda.
    pub energy_per_fission: EnergyChannels,
}

impl DecayHeatGroup {
    /// Admit a positive effective reservoir half-life [s] and its integrated
    /// energy per fission. This does not identify a physical nuclide or supply
    /// nuclear data. Under a continuing source, stored energy need not halve.
    pub fn from_half_life_s(
        half_life_s: f64,
        energy_per_fission: EnergyChannels,
    ) -> Result<Self, FissionError> {
        if !positive(half_life_s) {
            return Err(FissionError::InvalidInput("half-life must be positive and finite"));
        }
        energy_per_fission.validate()?;
        if !positive(energy_per_fission.total_j()) {
            return Err(FissionError::InvalidInput("decay groups require positive energy"));
        }
        Ok(Self {
            decay_constant_s_inv: finite(core::f64::consts::LN_2 / half_life_s)?,
            energy_per_fission,
        })
    }

    /// Shutdown half-life [s] of this effective energy reservoir, `ln(2)/lambda`.
    /// Public coefficients are checked here too; zero or invalid decay rates
    /// are not assigned a fictitious finite half-life.
    pub fn half_life_s(&self) -> Result<f64, FissionError> {
        if !positive(self.decay_constant_s_inv) {
            return Err(FissionError::InvalidInput("decay constant must be positive and finite"));
        }
        finite(core::f64::consts::LN_2 / self.decay_constant_s_inv)
    }
}

/// Immutable nuclear source data, kept distinct from chemical material laws.
#[derive(Clone, Debug, PartialEq)]
pub struct FissionModel {
    data_id: String,
    prompt_energy_per_fission: EnergyChannels,
    decay_groups: Vec<DecayHeatGroup>,
    total_energy_per_fission_j: f64,
}

impl FissionModel {
    /// Admit explicit, nonnegative source coefficients. The data identity names
    /// the caller's dataset, revision and applicability; it is not validation.
    pub fn new(
        data_id: String,
        prompt_energy_per_fission: EnergyChannels,
        decay_groups: Vec<DecayHeatGroup>,
    ) -> Result<Self, FissionError> {
        if data_id.trim().is_empty() {
            return Err(FissionError::InvalidInput("a nuclear source data identity is required"));
        }
        prompt_energy_per_fission.validate()?;
        let mut energy = Accumulator::new(prompt_energy_per_fission.total_j());
        for group in &decay_groups {
            group.energy_per_fission.validate()?;
            if !positive(group.decay_constant_s_inv) || !positive(group.energy_per_fission.total_j()) {
                return Err(FissionError::InvalidInput("decay groups require positive decay constants and energy"));
            }
            energy.add(group.energy_per_fission.total_j())?;
        }
        let total_energy_per_fission_j = energy.value();
        if !positive(total_energy_per_fission_j) {
            return Err(FissionError::InvalidInput("total fission energy must be positive and finite"));
        }
        Ok(Self { data_id, prompt_energy_per_fission, decay_groups, total_energy_per_fission_j })
    }

    /// Original source dataset identity.
    #[must_use]
    pub fn data_id(&self) -> &str { &self.data_id }

    /// Prompt energy destinations per accepted fission.
    #[must_use]
    pub const fn prompt_energy_per_fission(&self) -> EnergyChannels { self.prompt_energy_per_fission }

    /// Immutable, ordered delayed-energy groups.
    #[must_use]
    pub fn decay_groups(&self) -> &[DecayHeatGroup] { &self.decay_groups }

    /// Total accounted energy per fission, including outward and delayed energy.
    #[must_use]
    pub const fn total_energy_per_fission_j(&self) -> f64 { self.total_energy_per_fission_j }
}

/// Work and noncancelling roundoff limits for one source interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FissionBudget {
    /// Maximum admitted delayed groups; stepping is O(groups) in time and storage.
    pub max_decay_groups: usize,
    /// Absolute bound on aggregate local energy-account defects [J].
    pub energy_tolerance_j: f64,
    /// Absolute bound on represented fuel-debit defects [expected nuclei].
    pub fuel_tolerance_nuclei: f64,
}

impl FissionBudget {
    fn validate(self) -> Result<(), FissionError> {
        if !nonnegative(self.energy_tolerance_j) || !nonnegative(self.fuel_tolerance_nuclei) {
            return Err(FissionError::InvalidInput("fission tolerances must be finite and nonnegative"));
        }
        Ok(())
    }
}

/// Complete accepted source history. Clone is an independent physical checkpoint;
/// only immutable model data are shared. Counts are continuous expectations, not
/// stochastic integer nuclear events. No implicit equilibrium prehistory exists.
#[derive(Clone, Debug, PartialEq)]
pub struct FissionState {
    model: Arc<FissionModel>,
    remaining_fuel: Accumulator,
    fissions: Accumulator,
    delayed_energy_j: Vec<f64>,
    released: [Accumulator; 3],
}

impl FissionState {
    /// Bind finite fuel and explicit initial delayed energy to this exact model.
    /// Fresh fuel uses a zero reservoir for every group. Irradiated fuel needs
    /// an independently supplied history; temperature does not determine it.
    pub fn new(
        model: Arc<FissionModel>,
        remaining_fuel_nuclei: f64,
        delayed_energy_j: Vec<f64>,
    ) -> Result<Self, FissionError> {
        if !nonnegative(remaining_fuel_nuclei)
            || delayed_energy_j.len() != model.decay_groups.len()
            || !delayed_energy_j.iter().copied().all(nonnegative)
            || !sum(&delayed_energy_j).is_finite()
        {
            return Err(FissionError::InvalidInput("invalid fuel or source-bound delayed inventory"));
        }
        Ok(Self {
            model, remaining_fuel: Accumulator::new(remaining_fuel_nuclei),
            fissions: Accumulator::new(0.0), delayed_energy_j,
            released: [Accumulator::new(0.0); 3],
        })
    }

    /// Immutable source data retained by this history.
    #[must_use]
    pub fn model(&self) -> &FissionModel { &self.model }

    /// Remaining expected fuel nuclei; a compensated debit retains sub-ULP burns.
    #[must_use]
    pub fn remaining_fuel_nuclei(&self) -> f64 { self.remaining_fuel.value() }

    /// Fission events accepted since this state was constructed.
    #[must_use]
    pub fn cumulative_fissions(&self) -> f64 { self.fissions.value() }

    /// Remaining energy in each original model group [J].
    #[must_use]
    pub fn delayed_energy_j(&self) -> &[f64] { &self.delayed_energy_j }

    /// Total surviving delayed energy [J].
    #[must_use]
    pub fn stored_delayed_energy_j(&self) -> f64 { sum(&self.delayed_energy_j) }

    /// Energy actually released since construction, separate from initial storage.
    #[must_use]
    pub fn cumulative_released(&self) -> EnergyChannels {
        EnergyChannels {
            deposited_j: self.released[0].value(),
            radiation_j: self.released[1].value(),
            escaped_j: self.released[2].value(),
        }
    }

    /// Current decay-only power [W], independent of the prescribed prompt source.
    pub fn decay_power_w(&self) -> Result<EnergyChannels, FissionError> {
        let mut channels = [Accumulator::new(0.0); 3];
        for (group, &stored) in self.model.decay_groups.iter().zip(&self.delayed_energy_j) {
            let power = finite(stored * group.decay_constant_s_inv)?;
            add_channels(&mut channels, partition(power, group.energy_per_fission)?)?;
        }
        Ok(channel_values(channels))
    }

    /// Prepare an independent candidate for one constant-rate interval. Fuel
    /// exhaustion stops fission at N/F inside the interval; existing and newly
    /// produced reservoirs continue decaying for the rest of it. The report
    /// contains integrated joules, not an endpoint-power rectangle.
    ///
    /// Coefficients are a prescribed source reduction. Numerical energy closure
    /// does not bound the error of a fitted nuclear dataset or thermal splitting.
    pub fn trial(
        &self,
        duration_s: f64,
        fissions_per_s: f64,
        budget: FissionBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, FissionStep), FissionError> {
        checkpoint()?;
        budget.validate()?;
        if !nonnegative(duration_s) || !nonnegative(fissions_per_s) {
            return Err(FissionError::InvalidInput("duration and fission rate must be finite and nonnegative"));
        }
        if self.model.decay_groups.len() > budget.max_decay_groups {
            return Err(FissionError::Budget("delayed groups"));
        }
        let mut next = self.clone();
        let fuel = self.remaining_fuel.value();
        let (fissions, active_s, exhausted) = if duration_s == 0.0 || fissions_per_s == 0.0 || fuel == 0.0 {
            (0.0, 0.0, false)
        } else if duration_s >= fuel / fissions_per_s {
            let active = fuel / fissions_per_s;
            if !positive(active) { return Err(FissionError::Arithmetic); }
            (fuel, active, true)
        } else {
            (finite(fissions_per_s * duration_s)?, duration_s, false)
        };
        if fissions == 0.0 && active_s > 0.0 {
            return Err(FissionError::Arithmetic);
        }
        if exhausted {
            next.remaining_fuel = Accumulator::new(0.0);
        } else {
            next.remaining_fuel.add(-fissions)?;
        }
        if !nonnegative(next.remaining_fuel.value()) {
            return Err(FissionError::Arithmetic);
        }
        let fuel_residual_nuclei = finite(sum(&[
            next.remaining_fuel.hi, -self.remaining_fuel.hi,
            next.remaining_fuel.lo, -self.remaining_fuel.lo, fissions,
        ]))?;
        let fission_counter_error = next.fissions.add(fissions)?;
        if fuel_residual_nuclei.abs() + fission_counter_error.abs() > budget.fuel_tolerance_nuclei {
            return Err(FissionError::Budget("fuel accounting roundoff"));
        }
        let prompt = self.model.prompt_energy_per_fission.scaled(fissions)?;
        let mut released = [Accumulator::new(0.0); 3];
        add_channels(&mut released, prompt)?;
        let mut storage_change = Accumulator::new(0.0);
        let mut energy_error = Accumulator::new(0.0);
        let mut generated = Accumulator::new(prompt.total_j());
        let cooling_s = (duration_s - active_s).max(0.0);
        for (index, (group, &initial)) in self.model.decay_groups.iter()
            .zip(&self.delayed_energy_j).enumerate()
        {
            checkpoint()?;
            let q = group.energy_per_fission.total_j();
            let input = finite(fissions * q)?;
            generated.add(input)?;
            let exposure = decay_exposure(group.decay_constant_s_inv, duration_s)?;
            let cooling_exposure = decay_exposure(group.decay_constant_s_inv, cooling_s)?;
            let (retained_input, released_input) = source_energy(input, group.decay_constant_s_inv, active_s)?;
            let final_energy = finite(attenuate(initial, exposure)
                + attenuate(retained_input, cooling_exposure))?;
            // Stable integrals retain tiny heat releases when subtracting two
            // large reservoir totals would round them to zero.
            let output = finite(sum(&[
                initial * decayed_fraction(exposure), released_input,
                retained_input * decayed_fraction(cooling_exposure),
            ]))?;
            let destinations = partition(output, group.energy_per_fission)?;
            add_channels(&mut released, destinations)?;
            next.delayed_energy_j[index] = final_energy;
            storage_change.add(final_energy - initial)?;
            energy_error.add(sum(&[final_energy, -initial, destinations.total_j(), -input]).abs())?;
        }
        let released = channel_values(released);
        finite(released.total_j())?;
        for (accumulator, value) in next.released.iter_mut()
            .zip([released.deposited_j, released.radiation_j, released.escaped_j])
        {
            energy_error.add(accumulator.add(value)?.abs())?;
        }
        let stored_delayed_change_j = storage_change.value();
        // Independently anchor total generation to the same immutable Q and
        // accepted count, including roundoff in the component products/sum.
        let generated_energy_j = finite(fissions * self.model.total_energy_per_fission_j)?;
        energy_error.add((generated.value() - generated_energy_j).abs())?;
        let energy_residual_j = finite(sum(&[
            stored_delayed_change_j, released.total_j(), -generated_energy_j,
        ]))?;
        energy_error.add(energy_residual_j.abs())?;
        let energy_roundoff_j = energy_error.value();
        if energy_roundoff_j > budget.energy_tolerance_j {
            return Err(FissionError::Budget("energy accounting roundoff"));
        }
        checkpoint()?;
        Ok((next, FissionStep {
            duration_s, fissions, fission_duration_s: active_s, released,
            generated_energy_j, stored_delayed_change_j,
            energy_residual_j, energy_roundoff_j, fuel_residual_nuclei,
        }))
    }

    /// Atomically accept a trial. A refusal or cancellation leaves self intact.
    pub fn advance(
        &mut self, duration_s: f64, fissions_per_s: f64, budget: FissionBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<FissionStep, FissionError> {
        let (next, report) = self.trial(duration_s, fissions_per_s, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }
}

/// One accepted interval's extensive source, storage and export accounting.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FissionStep {
    /// Actual interval length [s].
    pub duration_s: f64,
    /// Accepted fissions, including any finite-inventory clipping.
    pub fissions: f64,
    /// Duration for which the prescribed fission rate actually acted [s].
    pub fission_duration_s: f64,
    /// Integrated energy released during this interval [J].
    pub released: EnergyChannels,
    /// Prompt plus newly produced delayed energy [J].
    pub generated_energy_j: f64,
    /// Signed change in surviving delayed storage [J].
    pub stored_delayed_change_j: f64,
    /// delta_delayed + released - generated [J].
    pub energy_residual_j: f64,
    /// Noncancelling sum of local and aggregate accounting defects [J].
    pub energy_roundoff_j: f64,
    /// delta_remaining_fuel + fissions [expected nuclei].
    pub fuel_residual_nuclei: f64,
}

/// Refusal without an accepted physical candidate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FissionError {
    /// Invalid source data, inventory or interval input.
    InvalidInput(&'static str),
    /// A declared work or accounting limit was exceeded.
    Budget(&'static str),
    /// A required physical quantity is not representable as finite f64.
    Arithmetic,
    /// Caller checkpoint requested cancellation.
    Cancelled,
}

impl core::fmt::Display for FissionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput(detail) => write!(f, "invalid fission input: {detail}"),
            Self::Budget(detail) => write!(f, "fission budget exceeded: {detail}"),
            Self::Arithmetic => f.write_str("fission arithmetic is not representable"),
            Self::Cancelled => f.write_str("fission advance cancelled"),
        }
    }
}
impl std::error::Error for FissionError {}

fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn nonnegative(value: f64) -> bool { value.is_finite() && value >= 0.0 }
fn finite(value: f64) -> Result<f64, FissionError> {
    if value.is_finite() { Ok(value) } else { Err(FissionError::Arithmetic) }
}

// Error-free two-sum, with a retained low part, prevents tiny repeated fuel
// debits and cumulative heat credits from disappearing against a large total.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Accumulator { hi: f64, lo: f64 }
impl Accumulator {
    const fn new(value: f64) -> Self { Self { hi: value, lo: 0.0 } }
    fn value(self) -> f64 { self.hi + self.lo }
    fn add(&mut self, value: f64) -> Result<f64, FissionError> {
        let old = *self;
        let s = finite(self.hi + value)?;
        let b = s - self.hi;
        let e = (self.hi - (s - b)) + (value - b);
        let correction = finite(self.lo + e)?;
        let hi = finite(s + correction)?;
        let lo = finite(correction - (hi - s))?;
        *self = Self { hi, lo };
        finite(sum(&[hi, -old.hi, lo, -old.lo, -value]))
    }
}

fn sum(values: &[f64]) -> f64 {
    let (mut total, mut correction) = (0.0_f64, 0.0_f64);
    for &value in values {
        let next = total + value;
        correction += if total.abs() >= value.abs() {
            (total - next) + value
        } else { (value - next) + total };
        total = next;
    }
    total + correction
}

fn add_channels(target: &mut [Accumulator; 3], values: EnergyChannels) -> Result<(), FissionError> {
    for (accumulator, value) in target.iter_mut()
        .zip([values.deposited_j, values.radiation_j, values.escaped_j])
    {
        accumulator.add(value)?;
    }
    Ok(())
}
fn channel_values(values: [Accumulator; 3]) -> EnergyChannels {
    EnergyChannels { deposited_j: values[0].value(), radiation_j: values[1].value(), escaped_j: values[2].value() }
}
fn partition(total: f64, per_fission: EnergyChannels) -> Result<EnergyChannels, FissionError> {
    let scale = per_fission.total_j();
    EnergyChannels {
        deposited_j: per_fission.deposited_j / scale,
        radiation_j: per_fission.radiation_j / scale,
        escaped_j: per_fission.escaped_j / scale,
    }.scaled(total)
}

// Overflow of an exposure is a resolved complete-decay limit, not overflow
// of a physical energy account. An underflowed positive exposure still refuses.
fn decay_exposure(lambda: f64, duration: f64) -> Result<f64, FissionError> {
    if duration == 0.0 { return Ok(0.0); }
    let x = lambda * duration;
    if x == 0.0 { return Err(FissionError::Arithmetic); }
    Ok(x)
}

// phi=(1-exp(-x))/x is the retained fraction of uniformly produced energy.
// Its complement needs a series to avoid cancellation for short intervals.
// Return energies rather than only fractions: input/x can remain representable
// even when x overflows. In that case both finite divisors exceed one, so the
// sequential divisions cannot overflow or prematurely erase a surviving result.
fn source_energy(input: f64, lambda: f64, duration: f64) -> Result<(f64, f64), FissionError> {
    if duration == 0.0 { return Ok((input, 0.0)); }
    let x = decay_exposure(lambda, duration)?;
    if x < 0.01 {
        let released = x * (0.5 + x * (-1.0 / 6.0 + x * (1.0 / 24.0
            + x * (-1.0 / 120.0 + x * (1.0 / 720.0
            + x * (-1.0 / 5040.0 + x / 40320.0))))));
        Ok((input * (1.0 - released), input * released))
    } else if x.is_infinite() {
        let retained = input / lambda / duration;
        Ok((retained, input - retained))
    } else {
        let retained = decayed_fraction(x) / x;
        Ok((input * retained, input * (1.0 - retained)))
    }
}

#[cfg(test)]
mod tests;

