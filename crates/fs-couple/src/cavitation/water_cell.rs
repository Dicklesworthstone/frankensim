//! Finite rigid-volume water storage with conservative mass and energy transport.
//!
//! The accepted coordinates are total water mass M and internal energy U.
//! Original `IsochoricWater` resolves v=V/M and u=U/M with IF97 Regions 1, 2
//! and 4. Its actual source state, finite phase masses/volumes and inversion
//! residuals remain visible; an approximate root never overwrites conserved U.
//! The original source supplies all sensible and latent calorics.
//!
//! A prescribed integrated stream carries mass, source enthalpy and entropy.
//! Pressure-flow work is already inside enthalpy. The boundary is rigid, so
//! no p*dV work is added. Withdrawal is explicitly well mixed and uses the
//! beginning-of-interval source enthalpy, not a claimed implicit endpoint.
//! This first-order held-donor partition requires interval refinement. The
//! caller's hydraulic owner supplies flow amounts/direction; this module is
//! not a valve, choking or momentum solver.
//!
//! Each complete update, including a paired donor/receiver transfer, publishes
//! only after source recovery, conservation checks and the final cancellation
//! poll. Total water, including vapor, is retained. The equilibrium approximation
//! excludes finite-rate nucleation, metastability, noncondensable gas, interface
//! energy and thermal gradients. The source domain is 273.16..=623.15 K; its
//! explicit temperature/pressure window may refuse an unsupported inventory.
//! Entropy receipts retain the raw IF97 values. Supplied heat has no authored
//! boundary temperature here, so an unaccounted entropy change is not certified
//! entropy production.

use fs_material::phase::water::flow::across::WaterFlowState;
use fs_material::phase::water::isochoric::{
    IsochoricWater, IsochoricWaterControls, IsochoricWaterError, IsochoricWaterState,
};
use fs_qty::Temperature;
use fs_thermochem::water::MassSpecificVolume;
use fs_thermochem::water::if97::If97State;

/// Frozen-phase acoustic observation of an accepted equilibrium cell.
#[path = "water_cell/acoustics.rs"]
pub mod acoustics;

/// Reciprocal finite-wall heat storage and native conducting-layer feedback.
#[path = "water_cell/thermal.rs"]
pub mod thermal;

/// Source-inversion work, gross transport and absolute physical budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterCellControls {
    pub equilibrium: IsochoricWaterControls,
    /// Sum of incoming and outgoing mass divided by initial cell mass.
    /// Strictly in (0,1); opposing streams cannot hide excessive movement.
    pub maximum_relative_mass_exchange: f64,
    /// Caller work bound, in 1..=256; every incoming parcel polls cancellation.
    pub maximum_inflow_parcels: usize,
    pub mass_tolerance_kg: f64,
    /// Budget for conserved first-law residual and each source recovery defect.
    pub energy_tolerance_j: f64,
    /// Budget for source and separately reconstructed phase volumes minus V.
    pub volume_tolerance_m3: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum WaterCellError {
    Invalid(&'static str),
    Source(IsochoricWaterError),
    ExcessiveMassExchange {
        relative: f64,
        maximum: f64,
    },
    Budget {
        quantity: &'static str,
        residual: f64,
        tolerance: f64,
    },
    Unrepresentable(&'static str),
    ClockMismatch,
    Cancelled,
}
impl From<IsochoricWaterError> for WaterCellError {
    fn from(error: IsochoricWaterError) -> Self {
        match error {
            IsochoricWaterError::Cancelled => Self::Cancelled,
            other => Self::Source(other),
        }
    }
}
impl core::fmt::Display for WaterCellError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite IF97 water cell: {self:?}")
    }
}
impl std::error::Error for WaterCellError {}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), WaterCellError> {
    if cancelled() {
        Err(WaterCellError::Cancelled)
    } else {
        Ok(())
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn finite(value: f64, what: &'static str) -> Result<f64, WaterCellError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(WaterCellError::Unrepresentable(what))
    }
}
fn budget(quantity: &'static str, residual: f64, tolerance: f64) -> Result<(), WaterCellError> {
    if !residual.is_finite() || residual.abs() > tolerance {
        Err(WaterCellError::Budget {
            quantity,
            residual,
            tolerance,
        })
    } else {
        Ok(())
    }
}

impl WaterCellControls {
    fn validate(self) -> Result<(), WaterCellError> {
        if !positive(self.maximum_relative_mass_exchange)
            || self.maximum_relative_mass_exchange >= 1.0
            || !(1..=256).contains(&self.maximum_inflow_parcels)
            || !positive(self.mass_tolerance_kg)
            || !positive(self.energy_tolerance_j)
            || !positive(self.volume_tolerance_m3)
            || !self.equilibrium.volume_tolerance_m3_kg.is_finite()
            || self.equilibrium.volume_tolerance_m3_kg < 0.0
            || !self.equilibrium.internal_energy_tolerance_j_kg.is_finite()
            || self.equilibrium.internal_energy_tolerance_j_kg < 0.0
        {
            return Err(WaterCellError::Invalid(
                "bounded gross exchange, parcel count and positive absolute budgets required",
            ));
        }
        Ok(())
    }
    fn source_controls(self, mass: f64) -> IsochoricWaterControls {
        let mut source = self.equilibrium;
        source.volume_tolerance_m3_kg = source
            .volume_tolerance_m3_kg
            .min(self.volume_tolerance_m3 / mass);
        source.internal_energy_tolerance_j_kg = source
            .internal_energy_tolerance_j_kg
            .min(self.energy_tolerance_j / mass);
        source
    }
}

/// Original water caloric owner retained by a transported parcel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WaterParcelSource {
    /// A directly evaluated original IF97 liquid or steam state.
    Single(If97State),
    /// Original IF97 pressure-device state, including dense/regional charts.
    Flow(WaterFlowState),
    /// A finite cell's original recovered equilibrium image.
    Equilibrium(IsochoricWaterState),
}
impl WaterParcelSource {
    pub fn specific_enthalpy_j_kg(self) -> f64 {
        match self {
            Self::Single(s) => s.specific_enthalpy().value(),
            Self::Flow(s) => s.specific_enthalpy_j_kg(),
            Self::Equilibrium(s) => s.specific_enthalpy_j_kg(),
        }
    }
    pub fn specific_entropy_j_kg_k(self) -> f64 {
        match self {
            Self::Single(s) => s.specific_entropy().value(),
            Self::Flow(s) => s.specific_entropy_j_kg_k(),
            Self::Equilibrium(s) => s.specific_entropy_j_kg_k(),
        }
    }
}

/// Positive integrated incoming mass with its actual source enthalpy/entropy.
/// External parcels represent a prescribed upstream owner. Use
/// `transfer_water` when that owner is another finite cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterParcel {
    mass_kg: f64,
    enthalpy_j: f64,
    entropy_j_k: f64,
    source: WaterParcelSource,
}
impl WaterParcel {
    pub fn new(mass_kg: f64, source: WaterParcelSource) -> Result<Self, WaterCellError> {
        if !positive(mass_kg) {
            return Err(WaterCellError::Invalid(
                "positive finite parcel mass required",
            ));
        }
        let h = source.specific_enthalpy_j_kg();
        let s = source.specific_entropy_j_kg_k();
        let enthalpy_j = finite(mass_kg * h, "parcel enthalpy")?;
        let entropy_j_k = finite(mass_kg * s, "parcel entropy")?;
        if (h != 0.0 && enthalpy_j == 0.0) || (s != 0.0 && entropy_j_k == 0.0) {
            return Err(WaterCellError::Unrepresentable("parcel calorics"));
        }
        Ok(Self {
            mass_kg,
            enthalpy_j,
            entropy_j_k,
            source,
        })
    }
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    pub const fn enthalpy_j(&self) -> f64 {
        self.enthalpy_j
    }
    pub const fn entropy_j_k(&self) -> f64 {
        self.entropy_j_k
    }
    pub const fn source(&self) -> WaterParcelSource {
        self.source
    }
}

/// Conserved inventory and its separately retained original-source image.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterCellState {
    time_s: f64,
    volume_m3: f64,
    mass_kg: f64,
    internal_energy_j: f64,
    source: IsochoricWaterState,
    source_internal_energy_j: f64,
    source_energy_residual_j: f64,
    source_volume_residual_m3: f64,
    phase_masses_kg: [f64; 2],
    phase_volumes_m3: [f64; 2],
}
impl WaterCellState {
    pub const fn source_state(&self) -> IsochoricWaterState {
        self.source
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    /// Conserved storage coordinate, not rounded back to an approximate root.
    pub const fn internal_energy_j(&self) -> f64 {
        self.internal_energy_j
    }
    pub const fn source_internal_energy_j(&self) -> f64 {
        self.source_internal_energy_j
    }
    pub const fn source_energy_residual_j(&self) -> f64 {
        self.source_energy_residual_j
    }
    pub const fn source_volume_residual_m3(&self) -> f64 {
        self.source_volume_residual_m3
    }
    /// [liquid, vapor], both finite masses drawn from the same water inventory.
    pub const fn phase_masses_kg(&self) -> [f64; 2] {
        self.phase_masses_kg
    }
    /// [liquid, vapor] volumes from the retained original phase endpoints.
    pub const fn phase_volumes_m3(&self) -> [f64; 2] {
        self.phase_volumes_m3
    }
    pub fn vapor_volume_fraction(&self) -> f64 {
        self.phase_volumes_m3[1] / (self.phase_volumes_m3[0] + self.phase_volumes_m3[1])
    }
}

/// One accepted complete-interval exchange. Flow work is inside enthalpy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterCellStep {
    pub initial: WaterCellState,
    pub state: WaterCellState,
    pub incoming_mass_kg: f64,
    pub incoming_enthalpy_j: f64,
    pub outgoing: Option<WaterParcel>,
    pub supplied_heat_j: f64,
    /// Exactly zero for this rigid control volume.
    pub boundary_work_j: f64,
    pub mass_change_kg: f64,
    pub internal_energy_change_j: f64,
    pub mass_residual_kg: f64,
    pub energy_residual_j: f64,
    /// Actual source-image delta U minus the same incoming/outgoing H and Q.
    /// This includes the difference of the two exposed inversion defects.
    pub source_first_law_residual_j: f64,
    pub source_entropy_change_j_k: f64,
    pub advected_entropy_j_k: f64,
    /// Source delta S minus advected S. Heat-boundary entropy is unspecified.
    pub unaccounted_entropy_change_j_k: f64,
    pub temperature_iterations: usize,
    pub pressure_iterations: usize,
}

/// Uniform equilibrium water in one fixed physical volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumWaterCell {
    controls: WaterCellControls,
    state: WaterCellState,
}

impl EquilibriumWaterCell {
    pub fn new(
        volume_m3: f64,
        mass_kg: f64,
        internal_energy_j: f64,
        time_s: f64,
        controls: WaterCellControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WaterCellError> {
        poll(&mut cancelled)?;
        admit_initial(volume_m3, mass_kg, internal_energy_j, time_s, controls)?;
        let (state, _, _) = recover(
            volume_m3,
            mass_kg,
            internal_energy_j,
            time_s,
            controls,
            None,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(Self { controls, state })
    }

    pub fn at_temperature(
        volume_m3: f64,
        mass_kg: f64,
        temperature: Temperature,
        time_s: f64,
        controls: WaterCellControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WaterCellError> {
        poll(&mut cancelled)?;
        admit_initial(volume_m3, mass_kg, 0.0, time_s, controls)?;
        let law = IsochoricWater::new(
            MassSpecificVolume::new(volume_m3 / mass_kg),
            controls.source_controls(mass_kg),
            &mut cancelled,
        )?;
        let source = law.state_at_temperature(temperature, &mut cancelled)?;
        let energy = finite(
            mass_kg * source.specific_internal_energy_j_kg(),
            "initial internal energy",
        )?;
        let state = image(volume_m3, mass_kg, energy, time_s, source, controls)?;
        poll(&mut cancelled)?;
        Ok(Self { controls, state })
    }

    pub const fn state(&self) -> WaterCellState {
        self.state
    }
    pub const fn controls(&self) -> WaterCellControls {
        self.controls
    }

    /// Apply prescribed integrated streams and heat at an increasing physical
    /// time. Withdrawal samples the initial well-mixed donor once. No caller
    /// state changes on an inadmissible inventory, solve failure or cancellation.
    pub fn advance(
        &mut self,
        inflows: &[WaterParcel],
        withdrawn_mass_kg: f64,
        supplied_heat_j: f64,
        end_time_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterCellStep, WaterCellError> {
        let step = self.propose(
            inflows,
            withdrawn_mass_kg,
            supplied_heat_j,
            end_time_s,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        self.state = step.state;
        Ok(step)
    }

    fn propose(
        &self,
        inflows: &[WaterParcel],
        withdrawn_mass_kg: f64,
        supplied_heat_j: f64,
        end_time_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WaterCellStep, WaterCellError> {
        poll(cancelled)?;
        if inflows.len() > self.controls.maximum_inflow_parcels
            || !withdrawn_mass_kg.is_finite()
            || withdrawn_mass_kg < 0.0
            || !supplied_heat_j.is_finite()
            || !end_time_s.is_finite()
            || end_time_s <= self.state.time_s
            || !positive(end_time_s - self.state.time_s)
        {
            return Err(WaterCellError::Invalid(
                "bounded inflows, nonnegative withdrawal, finite heat and increasing time required",
            ));
        }
        let (mut incoming_mass, mut incoming_h, mut incoming_s) = (0.0, 0.0, 0.0);
        for parcel in inflows {
            poll(cancelled)?;
            incoming_mass = finite(incoming_mass + parcel.mass_kg, "incoming mass sum")?;
            incoming_h = finite(incoming_h + parcel.enthalpy_j, "incoming enthalpy sum")?;
            incoming_s = finite(incoming_s + parcel.entropy_j_k, "incoming entropy sum")?;
        }
        let relative = (incoming_mass + withdrawn_mass_kg) / self.state.mass_kg;
        if !relative.is_finite() || relative > self.controls.maximum_relative_mass_exchange {
            return Err(WaterCellError::ExcessiveMassExchange {
                relative,
                maximum: self.controls.maximum_relative_mass_exchange,
            });
        }
        let outgoing = if withdrawn_mass_kg == 0.0 {
            None
        } else {
            Some(WaterParcel::new(
                withdrawn_mass_kg,
                WaterParcelSource::Equilibrium(self.state.source),
            )?)
        };
        let outgoing_h = outgoing.map_or(0.0, |p| p.enthalpy_j);
        let outgoing_s = outgoing.map_or(0.0, |p| p.entropy_j_k);
        let mass_change = finite(incoming_mass - withdrawn_mass_kg, "requested mass change")?;
        let energy_change = finite(
            incoming_h - outgoing_h + supplied_heat_j,
            "requested energy change",
        )?;
        let mass = finite(self.state.mass_kg + mass_change, "updated mass")?;
        let energy = finite(
            self.state.internal_energy_j + energy_change,
            "updated internal energy",
        )?;
        if mass <= 0.0
            || (mass_change != 0.0 && mass == self.state.mass_kg)
            || (energy_change != 0.0 && energy == self.state.internal_energy_j)
        {
            return Err(WaterCellError::Unrepresentable(
                "positive updated inventory or lost exchange",
            ));
        }
        let (state, temperature_iterations, pressure_iterations) =
            if mass == self.state.mass_kg && energy == self.state.internal_energy_j {
                let mut same = self.state;
                same.time_s = end_time_s;
                (same, 0, 0)
            } else {
                recover(
                    self.state.volume_m3,
                    mass,
                    energy,
                    end_time_s,
                    self.controls,
                    Some(self.state.source.law()),
                    cancelled,
                )?
            };
        let actual_mass_change = state.mass_kg - self.state.mass_kg;
        let actual_energy_change = state.internal_energy_j - self.state.internal_energy_j;
        let mass_residual = actual_mass_change - mass_change;
        let energy_residual = actual_energy_change - energy_change;
        budget(
            "conserved mass",
            mass_residual,
            self.controls.mass_tolerance_kg,
        )?;
        budget(
            "conserved energy",
            energy_residual,
            self.controls.energy_tolerance_j,
        )?;
        let source_residual = finite(
            state.source_internal_energy_j - self.state.source_internal_energy_j - energy_change,
            "source first-law residual",
        )?;
        let entropy_change = finite(
            state.mass_kg * state.source.specific_entropy_j_kg_k()
                - self.state.mass_kg * self.state.source.specific_entropy_j_kg_k(),
            "source entropy change",
        )?;
        let advected_entropy = finite(incoming_s - outgoing_s, "advected entropy")?;
        let unaccounted_entropy = finite(entropy_change - advected_entropy, "unaccounted entropy")?;
        poll(cancelled)?;
        Ok(WaterCellStep {
            initial: self.state,
            state,
            incoming_mass_kg: incoming_mass,
            incoming_enthalpy_j: incoming_h,
            outgoing,
            supplied_heat_j,
            boundary_work_j: 0.0,
            mass_change_kg: actual_mass_change,
            internal_energy_change_j: actual_energy_change,
            mass_residual_kg: mass_residual,
            energy_residual_j: energy_residual,
            source_first_law_residual_j: source_residual,
            source_entropy_change_j_k: entropy_change,
            advected_entropy_j_k: advected_entropy,
            unaccounted_entropy_change_j_k: unaccounted_entropy,
            temperature_iterations,
            pressure_iterations,
        })
    }
}

/// A completed two-cell transaction with one identical transported parcel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterCellTransfer {
    pub parcel: WaterParcel,
    pub donor: WaterCellStep,
    pub receiver: WaterCellStep,
    pub total_mass_residual_kg: f64,
    pub total_energy_residual_j: f64,
}

/// Transfer prescribed well-mixed water between two finite cells atomically.
/// `supplied_heat_j` is [donor, receiver]. Both accepted clocks must match.
/// This does not choose a flow amount or supply missing pump/kinetic work.
pub fn transfer_water(
    donor: &mut EquilibriumWaterCell,
    receiver: &mut EquilibriumWaterCell,
    mass_kg: f64,
    supplied_heat_j: [f64; 2],
    end_time_s: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<WaterCellTransfer, WaterCellError> {
    poll(&mut cancelled)?;
    if donor.state.time_s != receiver.state.time_s {
        return Err(WaterCellError::ClockMismatch);
    }
    if !positive(mass_kg) {
        return Err(WaterCellError::Invalid("positive transfer mass required"));
    }
    let donor_step = donor.propose(&[], mass_kg, supplied_heat_j[0], end_time_s, &mut cancelled)?;
    let parcel = donor_step
        .outgoing
        .ok_or(WaterCellError::Invalid("missing donor parcel"))?;
    let receiver_step = receiver.propose(
        &[parcel],
        0.0,
        supplied_heat_j[1],
        end_time_s,
        &mut cancelled,
    )?;
    let mass_residual = (donor_step.state.mass_kg + receiver_step.state.mass_kg)
        - (donor.state.mass_kg + receiver.state.mass_kg);
    let energy_residual = (donor_step.state.internal_energy_j
        + receiver_step.state.internal_energy_j)
        - (donor.state.internal_energy_j + receiver.state.internal_energy_j)
        - (supplied_heat_j[0] + supplied_heat_j[1]);
    budget(
        "two-cell mass",
        mass_residual,
        finite(
            donor.controls.mass_tolerance_kg + receiver.controls.mass_tolerance_kg,
            "paired mass budget",
        )?,
    )?;
    budget(
        "two-cell energy",
        energy_residual,
        finite(
            donor.controls.energy_tolerance_j + receiver.controls.energy_tolerance_j,
            "paired energy budget",
        )?,
    )?;
    poll(&mut cancelled)?;
    donor.state = donor_step.state;
    receiver.state = receiver_step.state;
    Ok(WaterCellTransfer {
        parcel,
        donor: donor_step,
        receiver: receiver_step,
        total_mass_residual_kg: mass_residual,
        total_energy_residual_j: energy_residual,
    })
}

fn admit_initial(
    volume: f64,
    mass: f64,
    energy: f64,
    time: f64,
    controls: WaterCellControls,
) -> Result<(), WaterCellError> {
    controls.validate()?;
    if !positive(volume) || !positive(mass) || !energy.is_finite() || !time.is_finite() {
        return Err(WaterCellError::Invalid(
            "positive finite volume/mass and finite energy/time required",
        ));
    }
    Ok(())
}

fn recover(
    volume: f64,
    mass: f64,
    energy: f64,
    time: f64,
    controls: WaterCellControls,
    previous: Option<&IsochoricWater>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(WaterCellState, usize, usize), WaterCellError> {
    poll(cancelled)?;
    let specific_volume = volume / mass;
    let specific_energy = energy / mass;
    if !positive(specific_volume)
        || !specific_energy.is_finite()
        || (energy != 0.0 && specific_energy == 0.0)
    {
        return Err(WaterCellError::Unrepresentable(
            "specific conserved coordinates",
        ));
    }
    let source_controls = controls.source_controls(mass);
    let law = if let Some(law) = previous.filter(|law| {
        law.fixed_specific_volume_m3_kg() == specific_volume && *law.controls() == source_controls
    }) {
        *law
    } else {
        IsochoricWater::new(
            MassSpecificVolume::new(specific_volume),
            source_controls,
            &mut *cancelled,
        )?
    };
    let solved = law.state_at_specific_internal_energy(specific_energy, &mut *cancelled)?;
    let state = image(volume, mass, energy, time, solved.into_state(), controls)?;
    poll(cancelled)?;
    Ok((
        state,
        solved.temperature_iterations(),
        solved.pressure_iterations(),
    ))
}

fn image(
    volume: f64,
    mass: f64,
    energy: f64,
    time: f64,
    source: IsochoricWaterState,
    controls: WaterCellControls,
) -> Result<WaterCellState, WaterCellError> {
    let source_energy = finite(
        mass * source.specific_internal_energy_j_kg(),
        "source internal energy",
    )?;
    let source_energy_residual = source_energy - energy;
    let source_volume_residual = mass * source.specific_volume_m3_kg() - volume;
    budget(
        "source energy recovery",
        source_energy_residual,
        controls.energy_tolerance_j,
    )?;
    budget(
        "source volume recovery",
        source_volume_residual,
        controls.volume_tolerance_m3,
    )?;
    let vapor_mass = finite(mass * source.vapor_mass_fraction(), "vapor mass")?;
    let liquid_mass = mass - vapor_mass;
    let phase_masses = [liquid_mass, vapor_mass];
    let volumes = if let Some(pair) = source.saturation_pair() {
        [
            liquid_mass * pair.liquid().specific_volume().value(),
            vapor_mass * pair.vapor().specific_volume().value(),
        ]
    } else if source.liquid_mass_fraction() == 1.0 {
        [mass * source.specific_volume_m3_kg(), 0.0]
    } else {
        [0.0, mass * source.specific_volume_m3_kg()]
    };
    for (index, fraction) in [source.liquid_mass_fraction(), source.vapor_mass_fraction()]
        .into_iter()
        .enumerate()
    {
        if !phase_masses[index].is_finite()
            || phase_masses[index] < 0.0
            || !volumes[index].is_finite()
            || volumes[index] < 0.0
            || (fraction > 0.0 && (phase_masses[index] == 0.0 || volumes[index] == 0.0))
        {
            return Err(WaterCellError::Unrepresentable("finite phase inventory"));
        }
    }
    budget(
        "phase mass reconstruction",
        liquid_mass + vapor_mass - mass,
        controls.mass_tolerance_kg,
    )?;
    budget(
        "phase volume reconstruction",
        volumes[0] + volumes[1] - volume,
        controls.volume_tolerance_m3,
    )?;
    Ok(WaterCellState {
        time_s: time,
        volume_m3: volume,
        mass_kg: mass,
        internal_energy_j: energy,
        source,
        source_internal_energy_j: source_energy,
        source_energy_residual_j: source_energy_residual,
        source_volume_residual_m3: source_volume_residual,
        phase_masses_kg: phase_masses,
        phase_volumes_m3: volumes,
    })
}

#[cfg(test)]
#[path = "water_cell/tests.rs"]
mod tests;
