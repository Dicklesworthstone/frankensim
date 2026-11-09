//! Pressure-driven exchange between two finite reacting, rigid gas volumes.
//!
//! Both tanks retain the same original elemental ratios and complete source bank.
//! Molecular composition, temperature and entropy may differ. The donor follows
//! its original depletion isentrope; its integrated U loss enters the receiver,
//! whose actual changing-density state is closed by the original UV owner.
//! Receiver inflow establishes a fresh accepted entropy reference. Jet energy is
//! thermalized once through total flow enthalpy; no separate kinetic heat enters.
//!
//! This is the ideal, well-mixed, instantaneous-equilibrium limit. It does not
//! mix different elemental feeds, conduct heat between equal-pressure tanks,
//! store port momentum, or model finite-rate reaction/valve losses.

mod exchange;

use super::{EquilibriumNozzleControl, EquilibriumNozzleError, EquilibriumNozzleFlow,
    publish_on_path, resolve_on_isentrope};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::equilibrium::blowdown::{
    EquilibriumBlowdownError, EquilibriumBlowdownSolve, EquilibriumBlowdownStep,
    EquilibriumGasReservoir,
    midpoint::{self, OutflowMidpointError, OutflowMidpointSolve, OutflowRate},
};
use fs_qty::{Amount, Area, Pressure, Temperature};
use fs_thermochem::{SpeciesId, equilibrium::composition::{
    CompositionEquilibriumErrorV1, EquilibriumExpansionErrorV1, IdealGasEquilibriumV1,
    isentrope::{EquilibriumIsentropeV1, density::EquilibriumDensityErrorV1},
    volume::{EquilibriumInternalEnergySolveV1, EquilibriumVolumeErrorV1},
}};
use std::sync::Arc;

/// A physical initial tank; mass follows from the original equilibrium EOS.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumTankInitialState {
    pub volume_m3: f64,
    pub temperature: Temperature,
    pub pressure: Pressure,
}

/// Immutable thermochemical and physical-time policies for the connected pair.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumGasPairControl {
    /// Original donor finite-energy and mass/EOS allowances.
    pub reservoir: EquilibriumBlowdownSolve,
    /// Original ordinary-branch reacting nozzle controls.
    pub nozzle: EquilibriumNozzleControl,
    /// Original receiving UV root, on the retained feed basis.
    pub receiver: EquilibriumInternalEnergySolveV1,
    /// Pair total entropy may fall by at most this numerical allowance [J/K].
    pub entropy_tolerance_j_k: f64,
    /// Actual total elemental-amount discrepancy, relative to each original row.
    pub relative_element_tolerance: f64,
    /// Fraction of current donor inventory per interval, in (0,0.5].
    pub maximum_outgoing_fraction: f64,
    /// At most 128 physical midpoint mass-rate evaluations.
    pub maximum_mass_iterations: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EquilibriumGasPairDirection {
    LeftToRight,
    RightToLeft,
}

/// Original physical reservoirs. Sources and elemental ratios cannot diverge.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumGasPair {
    chemistry: Arc<IdealGasEquilibriumV1>,
    feed: Vec<(SpeciesId, Amount)>,
    tanks: [EquilibriumGasReservoir; 2],
    reference: Arc<PairReference>,
    // Last independently admitted immutable-reference mass/U/S balances.
    reference_balance: [f64; 3],
    reference_element_residual: f64,
    controls: EquilibriumGasPairControl,
}

#[derive(Clone, Debug, PartialEq)]
struct PairReference {
    tanks: [EquilibriumGasReservoir; 2],
    right_minus_left_h: f64,
    right_minus_left_s: f64,
}

/// Complete conservative pair interval, including its physical midpoint flow.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumGasPairStep {
    pub duration_s: f64,
    pub direction: Option<EquilibriumGasPairDirection>,
    /// Original integrated donor energy/mass account; absent for exact rest.
    pub donor: Option<EquilibriumBlowdownStep>,
    /// Actual source-bearing midpoint, with finite receiver back pressure.
    pub midpoint_nozzle: Option<EquilibriumNozzleFlow>,
    /// Actual pair mass change from represented endpoint inventories [kg].
    pub mass_residual_kg: f64,
    /// Largest actual relative elemental-amount change from the INITIAL pair.
    pub maximum_relative_element_residual: f64,
    /// Receiver energy gain minus the original finite donor U loss [J].
    /// Source references cancel before summation; no rounded absolute subtraction.
    pub energy_residual_j: f64,
    /// Total actual pair entropy increment, including irreversible receiving [J/K].
    pub entropy_production_j_k: f64,
    /// Actual total mass change from the immutable INITIAL pair [kg].
    pub mass_change_since_reference_kg: f64,
    /// Direct original-source total U-U_initial [J], independently re-evaluated.
    /// The same energy allowance applies; it is not multiplied by step count.
    pub energy_change_since_reference_j: f64,
    /// Direct original-source total S-S_initial [J/K].
    pub entropy_change_since_reference_j_k: f64,
    /// Actual transferred mass minus dt times the published midpoint rate [kg].
    pub mass_quadrature_residual_kg: f64,
    pub mass_iterations: usize,
    pub density_pressure_queries: usize,
    pub receiver_pressure_queries: usize,
    pub receiver_composition_iterations: usize,
    pub nozzle_pressure_queries: usize,
}

/// Both tentative endpoints and their one shared transfer account.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumGasPairCandidate {
    endpoint: EquilibriumGasPair,
    step: EquilibriumGasPairStep,
}
impl EquilibriumGasPairCandidate {
    pub const fn endpoint(&self) -> &EquilibriumGasPair { &self.endpoint }
    pub const fn step(&self) -> &EquilibriumGasPairStep { &self.step }
    pub fn into_parts(self) -> (EquilibriumGasPair, EquilibriumGasPairStep) {
        (self.endpoint, self.step)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum EquilibriumGasPairError {
    Invalid(&'static str),
    Reservoir(EquilibriumBlowdownError),
    Receiver(EquilibriumVolumeErrorV1),
    Nozzle(EquilibriumNozzleError),
    InventoryLimit { maximum_duration_s: f64 },
    PressureCrossing { donor_pressure_pa: f64, receiver_pressure_pa: f64 },
    NoConvergence { iterations: usize, residual_kg: f64 },
    MassBalance { residual_kg: f64 },
    ElementBalance { maximum_relative_residual: f64 },
    EnergyBalance { residual_j: f64 },
    Entropy { production_j_k: f64 },
    Cancelled,
}
type Error = EquilibriumGasPairError;
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite reacting two-volume exchange: {self:?}")
    }
}
impl std::error::Error for Error {}
impl From<EquilibriumBlowdownError> for Error {
    fn from(error: EquilibriumBlowdownError) -> Self {
        match error {
            EquilibriumBlowdownError::Cancelled
            | EquilibriumBlowdownError::Thermochemistry(EquilibriumExpansionErrorV1::Equilibrium(
                CompositionEquilibriumErrorV1::Cancelled))
            | EquilibriumBlowdownError::Density(EquilibriumDensityErrorV1::Thermochemistry(
                EquilibriumExpansionErrorV1::Equilibrium(CompositionEquilibriumErrorV1::Cancelled)))
                => Self::Cancelled,
            other => Self::Reservoir(other),
        }
    }
}
impl From<EquilibriumExpansionErrorV1> for Error {
    fn from(error: EquilibriumExpansionErrorV1) -> Self {
        EquilibriumBlowdownError::Thermochemistry(error).into()
    }
}
impl From<EquilibriumVolumeErrorV1> for Error {
    fn from(error: EquilibriumVolumeErrorV1) -> Self {
        match error {
            EquilibriumVolumeErrorV1::Equilibrium(CompositionEquilibriumErrorV1::Cancelled)
            | EquilibriumVolumeErrorV1::Response(EquilibriumExpansionErrorV1::Equilibrium(
                CompositionEquilibriumErrorV1::Cancelled)) => Self::Cancelled,
            other => Self::Receiver(other),
        }
    }
}
impl From<EquilibriumNozzleError> for Error {
    fn from(error: EquilibriumNozzleError) -> Self {
        if matches!(error, EquilibriumNozzleError::Cancelled) { Self::Cancelled }
        else { Self::Nozzle(error) }
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}
fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }

impl EquilibriumGasPair {
    /// Assemble two original reservoirs with one exact source and elemental feed.
    /// Feed amounts specify the retained normalization, not either physical mass.
    /// No old outflow entropy is imposed on a tank after it receives gas.
    pub fn new(
        chemistry: Arc<IdealGasEquilibriumV1>,
        feed: Vec<(SpeciesId, Amount)>,
        left: EquilibriumTankInitialState,
        right: EquilibriumTankInitialState,
        controls: EquilibriumGasPairControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if !positive(controls.entropy_tolerance_j_k)
            || !positive(controls.relative_element_tolerance)
            || controls.relative_element_tolerance >= 0.01
            || !positive(controls.maximum_outgoing_fraction)
            || controls.maximum_outgoing_fraction > 0.5
            || !(1..=128).contains(&controls.maximum_mass_iterations)
            || controls.receiver.volume.equilibrium != controls.nozzle.isentropic.equilibrium
            || controls.receiver.minimum_temperature_k != controls.nozzle.isentropic.minimum_temperature_k
            || controls.receiver.maximum_temperature_k != controls.nozzle.isentropic.maximum_temperature_k {
            return Err(Error::Invalid("bounded pair work, positive balance allowances and one original source/Gibbs policy required"));
        }
        let make = |initial: EquilibriumTankInitialState, cancel: &mut _| {
            EquilibriumGasReservoir::new(
                chemistry.clone(), feed.clone(), initial.volume_m3,
                initial.temperature, initial.pressure, controls.nozzle.isentropic, cancel,
            )
        };
        let tanks = [make(left, &mut cancelled)?, make(right, &mut cancelled)?];
        // Original owners admit their complete policy even if this first state
        // has no pressure drive. Prepared work is numerical, never a new tank.
        for tank in &tanks {
            let _ = tank.prepare_depletion(controls.reservoir, &mut cancelled)?;
            super::validate_flow_inputs(
                tank.isentropic_sample().pressure.value(),
                tank.isentropic_sample().pressure.value(), 0.0, controls.nozzle,
            )?;
        }
        let basis = tanks[0].mass_kg() / tanks[0].physical_amount_scale();
        let _ = fs_thermochem::equilibrium::composition::volume::energy_path::EquilibriumEnergyVolumePathV1::new(
            &chemistry, left.temperature,
            fs_qty::Volume::new(left.volume_m3 * basis / tanks[0].mass_kg()),
            &feed, controls.receiver, &mut cancelled,
        )?;
        let right_path = EquilibriumIsentropeV1::new(
            &chemistry, right.temperature, right.pressure, &feed,
            controls.nozzle.isentropic, &mut cancelled,
        )?;
        let reference = Arc::new(PairReference {
            right_minus_left_h: right_path.specific_enthalpy_change_from(
                tanks[0].caloric_anchor(), &mut cancelled)?,
            right_minus_left_s: right_path.specific_entropy_change_from(
                tanks[0].caloric_anchor(), &mut cancelled)?,
            tanks: tanks.clone(),
        });
        poll(&mut cancelled)?;
        Ok(Self { chemistry, feed, tanks, reference, reference_balance: [0.0; 3], reference_element_residual: 0.0, controls })
    }
    pub const fn tanks(&self) -> &[EquilibriumGasReservoir; 2] { &self.tanks }
    pub const fn controls(&self) -> EquilibriumGasPairControl { self.controls }

    /// Advance one ideal restriction, choosing direction from actual pressures.
    ///
    /// The shared midpoint owner closes dm=dt*mdot at both evolving inventories.
    /// Every trial receiver absorbs the donor's exact finite U loss to that mass.
    /// Thus midpoint approximates only elapsed time, never the energy transfer.
    /// An endpoint pressure crossing, exhausted work, source failure or final
    /// cancellation publishes neither tank. Refine the interval; no clipping.
    pub fn advance(
        &self, effective_area: Area, duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EquilibriumGasPairCandidate, Error> {
        poll(&mut cancelled)?;
        if !duration_s.is_finite() || duration_s <= 0.0
            || !effective_area.value().is_finite() || effective_area.value() < 0.0 {
            return Err(Error::Invalid("positive finite time and finite nonnegative effective area required"));
        }
        let [left, right] = self.tanks.each_ref().map(|tank| tank.isentropic_sample().pressure.value());
        if effective_area.value() == 0.0 || left == right {
            poll(&mut cancelled)?;
            return Ok(EquilibriumGasPairCandidate { endpoint: self.clone(), step: EquilibriumGasPairStep {
                duration_s, direction: None, donor: None, midpoint_nozzle: None,
                mass_residual_kg: 0.0, maximum_relative_element_residual: self.reference_element_residual,
                energy_residual_j: 0.0, entropy_production_j_k: 0.0,
                mass_change_since_reference_kg: self.reference_balance[0],
                energy_change_since_reference_j: self.reference_balance[1],
                entropy_change_since_reference_j_k: self.reference_balance[2],
                mass_quadrature_residual_kg: 0.0, mass_iterations: 0,
                density_pressure_queries: 0, receiver_pressure_queries: 0,
                receiver_composition_iterations: 0, nozzle_pressure_queries: 0,
            } });
        }
        let donor = usize::from(right > left);
        let mut work = exchange::Workspace::new(self, donor, &mut cancelled)?;
        let old = &self.tanks[donor];
        let receiver = &self.tanks[1 - donor];
        let first = resolve_on_isentrope(
            work.depletion.isentrope(), old.isentropic_sample(), old.caloric_anchor(),
            receiver.isentropic_sample().pressure.value(), effective_area.value(),
            self.controls.nozzle, &mut cancelled,
        )?;
        work.nozzle_queries += first.work.pressure;
        let quadrature = midpoint::solve(
            old.mass_kg(), first.mass_flow_kg_s,
            OutflowMidpointSolve {
                duration_s,
                maximum_outgoing_fraction: self.controls.maximum_outgoing_fraction,
                mass_tolerance_kg: self.controls.reservoir.mass_tolerance_kg,
                maximum_iterations: self.controls.maximum_mass_iterations,
            },
            |mass, mut cancel| -> Result<OutflowRate, Error> {
                let trial = work.trial(mass, &mut cancel)?;
                let back = trial.receiver.sample.pressure.value();
                if trial.donor.pressure.value() <= back {
                    // A crossed trial only bounds the outward root. Publication
                    // still requires the complete endpoint pressure inequality.
                    return Ok(OutflowRate { mass_flow_kg_s: 0.0, derivative_per_s: None });
                }
                let anchor = work.depletion.isentrope().caloric_anchor(&mut cancel)?;
                let flow = resolve_on_isentrope(
                    work.depletion.isentrope(), trial.donor, &anchor, back,
                    effective_area.value(), self.controls.nozzle, &mut cancel,
                )?;
                work.nozzle_queries += flow.work.pressure;
                let derivative_per_s = work.rate_derivative(
                    &trial, flow.mass_flow_kg_s, flow.velocity_m_s,
                    flow.accepted.sample.response.density_kg_m3,
                    flow.accepted.sample.response.sound_speed_m_s, flow.bracket.is_some(),
                );
                Ok(OutflowRate { mass_flow_kg_s: flow.mass_flow_kg_s, derivative_per_s })
            },
            &mut cancelled,
        ).map_err(|error| match error {
            OutflowMidpointError::Invalid(reason) => Error::Invalid(reason),
            OutflowMidpointError::Source(source) => source,
            OutflowMidpointError::InventoryLimit { maximum_duration_s } => Error::InventoryLimit { maximum_duration_s },
            OutflowMidpointError::NoConvergence { iterations, mass_residual_kg } => Error::NoConvergence { iterations, residual_kg: mass_residual_kg },
            OutflowMidpointError::Cancelled => Error::Cancelled,
        })?;
        let midpoint = work.trial(quadrature.midpoint_mass_kg, &mut cancelled)?;
        let back = midpoint.receiver.sample.pressure.value();
        if midpoint.donor.pressure.value() <= back {
            return Err(Error::Invalid("positive pair transfer has no outward midpoint pressure drive"));
        }
        let anchor = work.depletion.isentrope().caloric_anchor(&mut cancelled)?;
        let inlet = work.depletion.isentrope().publish_current(&mut cancelled)?;
        let flow = resolve_on_isentrope(
            work.depletion.isentrope(), midpoint.donor, &anchor, back,
            effective_area.value(), self.controls.nozzle, &mut cancelled,
        )?;
        work.nozzle_queries += flow.work.pressure;
        let midpoint_nozzle = publish_on_path(
            work.depletion.isentrope(), inlet, midpoint.donor, flow, back,
            effective_area.value(), &mut cancelled,
        )?;
        let residual = quadrature.expelled_mass_kg - duration_s * midpoint_nozzle.mass_flow.value();
        if !residual.is_finite() || residual.abs() > self.controls.reservoir.mass_tolerance_kg {
            return Err(Error::NoConvergence { iterations: quadrature.iterations, residual_kg: residual });
        }
        let final_trial = work.trial(quadrature.remaining_mass_kg, &mut cancelled)?;
        work.finish(final_trial, midpoint_nozzle, duration_s, residual,
            quadrature.iterations, &mut cancelled)
    }
}
