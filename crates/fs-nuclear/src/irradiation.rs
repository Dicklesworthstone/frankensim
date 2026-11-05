//! Fission and competing capture under a supplied multigroup neutron flux.
//!
//! A homogeneous parent inventory obeys `dN/dt = -(k_f + k_c) N`, with
//! `k_r = sum_g sigma_r,g * phi_g`. The incident scalar flux is prescribed;
//! there is no transport, multiplication, self-shielding or criticality solve.
//! Microscopic coefficients are either explicitly fixed or interpolated from
//! identified temperature tables. The actual material temperature selects the
//! coefficients frozen for each interval; no Doppler broadening, unresolved
//! resonance reconstruction or spectrum-dependent self-shielding is inferred.
//!
//! Fission and capture debit the same parent inventory. Capture products and
//! fission fragments leave that tracked inventory but remain material: this is
//! not a daughter-production chain and does not subtract bulk thermal mass.
//! Prompt capture energy is supplied explicitly, in the same three destinations
//! as fission energy. Delayed fission energy uses the exact exponential source
//! history, including the coincident depletion/decay-rate limit.

use std::sync::Arc;

use fs_math::decay::{attenuate, decayed_fraction};

use crate::cross_sections::{CrossSectionTable, NeutronReaction};

use super::{
    Accumulator, EnergyChannels, FissionBudget, FissionError, FissionModel, FissionState,
    add_channels, channel_values, finite, nonnegative, partition, positive, sum,
};

/// Identified nuclide, without a claim that an evaluated dataset exists for it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NuclideId {
    atomic_number: u16,
    mass_number: u16,
    isomer: u16,
}

impl NuclideId {
    /// `isomer = 0` denotes the ground state; positive indices are caller-defined.
    pub fn new(atomic_number: u16, mass_number: u16, isomer: u16) -> Result<Self, FissionError> {
        if atomic_number == 0 || mass_number < atomic_number {
            return Err(FissionError::InvalidInput(
                "invalid nuclide atomic/mass number",
            ));
        }
        Ok(Self {
            atomic_number,
            mass_number,
            isomer,
        })
    }

    #[must_use]
    pub const fn atomic_number(self) -> u16 {
        self.atomic_number
    }
    #[must_use]
    pub const fn mass_number(self) -> u16 {
        self.mass_number
    }
    #[must_use]
    pub const fn isomer(self) -> u16 {
        self.isomer
    }
}

/// Immutable, explicitly identified ascending energy-group boundaries [J].
#[derive(Clone, Debug, PartialEq)]
pub struct NeutronGroups {
    data_id: String,
    edges_j: Arc<[f64]>,
}

impl NeutronGroups {
    /// Each flux/cross-section entry belongs to the corresponding adjacent edges.
    pub fn new(data_id: String, edges_j: Vec<f64>) -> Result<Self, FissionError> {
        if data_id.trim().is_empty()
            || edges_j.len() < 2
            || !edges_j.iter().copied().all(nonnegative)
            || !edges_j.windows(2).all(|pair| pair[0] < pair[1])
        {
            return Err(FissionError::InvalidInput(
                "identified ascending neutron energy groups are required",
            ));
        }
        Ok(Self {
            data_id,
            edges_j: edges_j.into(),
        })
    }

    #[must_use]
    pub fn data_id(&self) -> &str {
        &self.data_id
    }
    #[must_use]
    pub fn edges_j(&self) -> &[f64] {
        &self.edges_j
    }
    #[must_use]
    pub fn len(&self) -> usize {
        self.edges_j.len() - 1
    }
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        false
    }
}

/// Energy-integrated scalar neutron flux in each admitted group [1/(m² s)].
/// This is an absolute flux, not a normalized spectrum or a surface current.
#[derive(Clone, Debug, PartialEq)]
pub struct GroupFlux {
    groups: Arc<NeutronGroups>,
    values_m2_s: Arc<[f64]>,
}

impl GroupFlux {
    pub fn new(groups: Arc<NeutronGroups>, values_m2_s: Vec<f64>) -> Result<Self, FissionError> {
        if values_m2_s.len() != groups.len() || !values_m2_s.iter().copied().all(nonnegative) {
            return Err(FissionError::InvalidInput(
                "flux must provide one finite nonnegative value per group",
            ));
        }
        Ok(Self {
            groups,
            values_m2_s: values_m2_s.into(),
        })
    }

    #[must_use]
    pub fn groups(&self) -> &Arc<NeutronGroups> {
        &self.groups
    }
    #[must_use]
    pub fn values_m2_s(&self) -> &[f64] {
        &self.values_m2_s
    }
}

/// Microscopic cross sections per parent nucleus, in square metres.
/// Capture here means removal without fission; scattering is not absorption.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct GroupCrossSections {
    pub fission_m2: f64,
    pub capture_m2: f64,
}

/// Explicit data supplied for one fissionable parent and one energy-group grid.
/// Dataset identities and temperature bounds are declarations, not validation.
#[derive(Clone, Debug, PartialEq)]
pub struct IrradiationData {
    pub data_id: String,
    pub nuclide: NuclideId,
    pub groups: Arc<NeutronGroups>,
    pub cross_sections: Vec<GroupCrossSections>,
    pub temperature_min_k: f64,
    pub temperature_max_k: f64,
    pub fission_model: Arc<FissionModel>,
    /// Immediate energy per capture [J/capture]. No delayed capture heating is inferred.
    pub capture_energy_per_event: EnergyChannels,
}

/// Temperature-dependent fission and capture data on the same physical parent.
/// Both tables must cover the entire declared model temperature domain and
/// carry identical nuclide and energy-group identities. Event energies retain
/// their original owners; a table does not infer deposition or decay-heat data.
#[derive(Clone, Debug, PartialEq)]
pub struct TabulatedIrradiationData {
    pub data_id: String,
    pub fission: Arc<CrossSectionTable>,
    pub capture: Arc<CrossSectionTable>,
    pub temperature_min_k: f64,
    pub temperature_max_k: f64,
    pub fission_model: Arc<FissionModel>,
    pub capture_energy_per_event: EnergyChannels,
}

#[derive(Clone, Debug, PartialEq)]
struct ReactionTables {
    fission: Arc<CrossSectionTable>,
    capture: Arc<CrossSectionTable>,
}

/// Admitted immutable reaction coefficients. There is no built-in isotope data.
#[derive(Clone, Debug, PartialEq)]
pub struct IrradiationModel {
    data: IrradiationData,
    tables: Option<ReactionTables>,
}

impl IrradiationModel {
    pub fn new(data: IrradiationData) -> Result<Self, FissionError> {
        if data.data_id.trim().is_empty()
            || data.cross_sections.len() != data.groups.len()
            || !positive(data.temperature_min_k)
            || !positive(data.temperature_max_k)
            || data.temperature_max_k < data.temperature_min_k
            || !data
                .cross_sections
                .iter()
                .all(|xs| nonnegative(xs.fission_m2) && nonnegative(xs.capture_m2))
        {
            return Err(FissionError::InvalidInput(
                "invalid microscopic reaction data or temperature domain",
            ));
        }
        data.capture_energy_per_event.validate()?;
        Ok(Self { data, tables: None })
    }

    /// Bind processed temperature tables without changing the inventory or
    /// fission-energy solver. A finite lower-domain reference row is retained
    /// for the legacy `cross_sections` observation; actual rates ALWAYS evaluate
    /// the tables at the temperature passed to `rate_constants`.
    pub fn new_tabulated(data: TabulatedIrradiationData) -> Result<Self, FissionError> {
        if data.fission.reaction() != NeutronReaction::Fission
            || data.capture.reaction() != NeutronReaction::Capture
            || data.fission.nuclide() != data.capture.nuclide()
            || data.fission.groups() != data.capture.groups()
            || !positive(data.temperature_min_k)
            || !positive(data.temperature_max_k)
            || data.temperature_min_k > data.temperature_max_k
        {
            return Err(FissionError::InvalidInput(
                "inconsistent tabulated irradiation identities or domain",
            ));
        }
        for table in [&data.fission, &data.capture] {
            table.validate_temperature(data.temperature_min_k)?;
            table.validate_temperature(data.temperature_max_k)?;
        }
        let mut reference = Vec::with_capacity(data.fission.groups().len());
        for group in 0..data.fission.groups().len() {
            reference.push(GroupCrossSections {
                fission_m2: data
                    .fission
                    .cross_section_at_group(data.temperature_min_k, group)?,
                capture_m2: data
                    .capture
                    .cross_section_at_group(data.temperature_min_k, group)?,
            });
        }
        let mut model = Self::new(IrradiationData {
            data_id: data.data_id,
            nuclide: data.fission.nuclide(),
            groups: Arc::clone(data.fission.groups()),
            cross_sections: reference,
            temperature_min_k: data.temperature_min_k,
            temperature_max_k: data.temperature_max_k,
            fission_model: data.fission_model,
            capture_energy_per_event: data.capture_energy_per_event,
        })?;
        model.tables = Some(ReactionTables {
            fission: data.fission,
            capture: data.capture,
        });
        Ok(model)
    }

    #[must_use]
    pub fn data_id(&self) -> &str {
        &self.data.data_id
    }
    #[must_use]
    pub const fn nuclide(&self) -> NuclideId {
        self.data.nuclide
    }
    #[must_use]
    pub fn groups(&self) -> &Arc<NeutronGroups> {
        &self.data.groups
    }
    /// Fixed coefficients, or the lower-domain reference row for a tabulated
    /// model. The reference temperature is exposed separately; actual reaction
    /// evaluation must use `rate_constants`, which evaluates the current law.
    #[must_use]
    pub fn cross_sections(&self) -> &[GroupCrossSections] {
        &self.data.cross_sections
    }
    /// The retained row is temperature-independent for a fixed model. In a
    /// tabulated model it is the reference at this explicitly reported Kelvin
    /// temperature; use `rate_constants` for all physical reaction evaluations.
    #[must_use]
    pub fn cross_section_reference_temperature_k(&self) -> Option<f64> {
        self.tables.as_ref().map(|_| self.temperature_min_k())
    }
    /// Original immutable fission/capture tables, when this model uses them.
    #[must_use]
    pub fn temperature_tables(&self) -> Option<(&Arc<CrossSectionTable>, &Arc<CrossSectionTable>)> {
        self.tables
            .as_ref()
            .map(|tables| (&tables.fission, &tables.capture))
    }
    #[must_use]
    pub const fn temperature_min_k(&self) -> f64 {
        self.data.temperature_min_k
    }
    #[must_use]
    pub const fn temperature_max_k(&self) -> f64 {
        self.data.temperature_max_k
    }
    #[must_use]
    pub fn fission_model(&self) -> &FissionModel {
        &self.data.fission_model
    }
    #[must_use]
    pub const fn capture_energy_per_event(&self) -> EnergyChannels {
        self.data.capture_energy_per_event
    }

    /// Evaluate the microscopic law at the actual temperature, then collapse
    /// it against the supplied absolute group flux. These rates remain frozen
    /// within one inventory trial; changing thermal histories require refinement.
    /// Grid identity and every edge must match, even for an all-zero flux.
    pub fn rate_constants(
        &self,
        flux: &GroupFlux,
        temperature_k: f64,
        budget: IrradiationBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<ReactionRateConstants, FissionError> {
        checkpoint()?;
        budget.source.validate()?;
        if self.data.groups.len() > budget.max_neutron_groups {
            return Err(FissionError::Budget("neutron groups"));
        }
        if self.data.fission_model.decay_groups.len() > budget.source.max_decay_groups {
            return Err(FissionError::Budget("delayed groups"));
        }
        if !positive(temperature_k)
            || temperature_k < self.temperature_min_k()
            || temperature_k > self.temperature_max_k()
        {
            return Err(FissionError::InvalidInput(
                "temperature is outside the cross-section data domain",
            ));
        }
        if !Arc::ptr_eq(&self.data.groups, &flux.groups) && self.data.groups != flux.groups {
            return Err(FissionError::InvalidInput(
                "neutron flux and cross-section energy groups do not match",
            ));
        }
        if let Some(tables) = &self.tables {
            let result = ReactionRateConstants {
                fission_s_inv: tables.fission.rate_constant(
                    flux,
                    temperature_k,
                    budget.max_neutron_groups,
                    &mut checkpoint,
                )?,
                capture_s_inv: tables.capture.rate_constant(
                    flux,
                    temperature_k,
                    budget.max_neutron_groups,
                    &mut checkpoint,
                )?,
            };
            finite(result.total_s_inv())?;
            checkpoint()?;
            return Ok(result);
        }
        let mut fission = Accumulator::new(0.0);
        let mut capture = Accumulator::new(0.0);
        for (xs, &value) in self.data.cross_sections.iter().zip(flux.values_m2_s.iter()) {
            checkpoint()?;
            fission.add(rate_product(xs.fission_m2, value)?)?;
            capture.add(rate_product(xs.capture_m2, value)?)?;
        }
        let result = ReactionRateConstants {
            fission_s_inv: fission.value(),
            capture_s_inv: capture.value(),
        };
        finite(result.total_s_inv())?;
        checkpoint()?;
        Ok(result)
    }
}

/// Bounded neutron-group work plus the original source/accounting budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IrradiationBudget {
    pub max_neutron_groups: usize,
    pub source: FissionBudget,
}

/// Per-parent reaction rate constants for the supplied flux [1/s].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ReactionRateConstants {
    pub fission_s_inv: f64,
    pub capture_s_inv: f64,
}

impl ReactionRateConstants {
    #[must_use]
    pub fn total_s_inv(self) -> f64 {
        sum(&[self.fission_s_inv, self.capture_s_inv])
    }
}

/// Instantaneous expected events per second at the accepted parent inventory.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ReactionRates {
    pub fissions_per_s: f64,
    pub captures_per_s: f64,
}

/// Complete accepted irradiation history, with the original fission reservoirs.
/// Clones share only immutable model/grid data, so trials never alias inventory.
#[derive(Clone, Debug, PartialEq)]
pub struct IrradiationState {
    model: Arc<IrradiationModel>,
    source: FissionState,
    captures: Accumulator,
}

impl IrradiationState {
    pub fn new(
        model: Arc<IrradiationModel>,
        parent_nuclei: f64,
        delayed_energy_j: Vec<f64>,
    ) -> Result<Self, FissionError> {
        let source = FissionState::new(
            Arc::clone(&model.data.fission_model),
            parent_nuclei,
            delayed_energy_j,
        )?;
        Ok(Self {
            model,
            source,
            captures: Accumulator::new(0.0),
        })
    }

    #[must_use]
    pub fn model(&self) -> &IrradiationModel {
        &self.model
    }
    #[must_use]
    pub fn remaining_parent_nuclei(&self) -> f64 {
        self.source.remaining_fuel_nuclei()
    }
    #[must_use]
    pub fn cumulative_fissions(&self) -> f64 {
        self.source.cumulative_fissions()
    }
    #[must_use]
    pub fn cumulative_captures(&self) -> f64 {
        self.captures.value()
    }
    #[must_use]
    pub fn delayed_energy_j(&self) -> &[f64] {
        self.source.delayed_energy_j()
    }
    #[must_use]
    pub fn stored_delayed_energy_j(&self) -> f64 {
        self.source.stored_delayed_energy_j()
    }
    /// All released reaction energy since construction, including capture energy.
    #[must_use]
    pub fn cumulative_released(&self) -> EnergyChannels {
        self.source.cumulative_released()
    }
    pub fn decay_power_w(&self) -> Result<EnergyChannels, FissionError> {
        self.source.decay_power_w()
    }

    pub fn reaction_rates(
        &self,
        flux: &GroupFlux,
        temperature_k: f64,
        budget: IrradiationBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<ReactionRates, FissionError> {
        let constants = self
            .model
            .rate_constants(flux, temperature_k, budget, checkpoint)?;
        Ok(ReactionRates {
            fissions_per_s: rate_product(self.remaining_parent_nuclei(), constants.fission_s_inv)?,
            captures_per_s: rate_product(self.remaining_parent_nuclei(), constants.capture_s_inv)?,
        })
    }

    /// Instantaneous prompt fission/capture plus accepted decay power [W].
    pub fn power_w(
        &self,
        flux: &GroupFlux,
        temperature_k: f64,
        budget: IrradiationBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<EnergyChannels, FissionError> {
        let rates = self.reaction_rates(flux, temperature_k, budget, &mut checkpoint)?;
        let mut power = [Accumulator::new(0.0); 3];
        add_channels(
            &mut power,
            self.source
                .model
                .prompt_energy_per_fission
                .scaled(rates.fissions_per_s)?,
        )?;
        add_channels(
            &mut power,
            self.model
                .data
                .capture_energy_per_event
                .scaled(rates.captures_per_s)?,
        )?;
        for (group, &stored) in self
            .source
            .model
            .decay_groups
            .iter()
            .zip(&self.source.delayed_energy_j)
        {
            checkpoint()?;
            add_channels(
                &mut power,
                partition(
                    finite(stored * group.decay_constant_s_inv)?,
                    group.energy_per_fission,
                )?,
            )?;
        }
        let result = channel_values(power);
        finite(result.total_j())?;
        checkpoint()?;
        Ok(result)
    }

    /// Analytic frozen-flux trial. The same integrated reactions debit the
    /// parent and create energy. No endpoint-rate rectangle or decay substep is
    /// used. Numerical account closure does not bound nuclear-data uncertainty,
    /// spectral collapse error or error from holding an external flux fixed.
    pub fn trial(
        &self,
        duration_s: f64,
        flux: &GroupFlux,
        temperature_k: f64,
        budget: IrradiationBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, IrradiationStep), FissionError> {
        checkpoint()?;
        if !nonnegative(duration_s) {
            return Err(FissionError::InvalidInput(
                "irradiation duration must be finite and nonnegative",
            ));
        }
        let rates = self
            .model
            .rate_constants(flux, temperature_k, budget, &mut checkpoint)?;
        if duration_s == 0.0 {
            checkpoint()?;
            return Ok((self.clone(), IrradiationStep::default()));
        }
        let removal_rate = rates.total_s_inv();
        let exposure = finite(removal_rate * duration_s)?;
        if removal_rate > 0.0 && exposure == 0.0 {
            return Err(FissionError::Arithmetic);
        }
        let parent = self.remaining_parent_nuclei();
        let removed = finite(parent * decayed_fraction(exposure))?;
        let (fissions, captures) = if removal_rate == 0.0 {
            (0.0, 0.0)
        } else {
            (
                reaction_share(removed, rates.fission_s_inv, removal_rate)?,
                reaction_share(removed, rates.capture_s_inv, removal_rate)?,
            )
        };
        let mut next = self.clone();
        let mut parent_error = Accumulator::new(sum(&[fissions, captures, -removed]).abs());
        if exposure < 0.5 {
            // Small repeated removals must not disappear against a large parent.
            parent_error.add(next.source.remaining_fuel.add(-fissions)?.abs())?;
            parent_error.add(next.source.remaining_fuel.add(-captures)?.abs())?;
        } else {
            // For substantial burn, subtraction would erase the surviving tail.
            next.source.remaining_fuel = Accumulator::new(scaled_decay(parent, exposure)?);
        }
        if !nonnegative(next.remaining_parent_nuclei()) {
            return Err(FissionError::Arithmetic);
        }
        let parent_residual_nuclei = finite(sum(&[
            next.source.remaining_fuel.hi,
            -self.source.remaining_fuel.hi,
            next.source.remaining_fuel.lo,
            -self.source.remaining_fuel.lo,
            fissions,
            captures,
        ]))?;
        parent_error.add(parent_residual_nuclei.abs())?;
        parent_error.add(next.source.fissions.add(fissions)?.abs())?;
        parent_error.add(next.captures.add(captures)?.abs())?;
        let parent_roundoff_nuclei = parent_error.value();
        if parent_roundoff_nuclei > budget.source.fuel_tolerance_nuclei {
            return Err(FissionError::Budget("parent inventory accounting roundoff"));
        }
        let fission_prompt = self
            .source
            .model
            .prompt_energy_per_fission
            .scaled(fissions)?;
        let capture_prompt = self.model.data.capture_energy_per_event.scaled(captures)?;
        let mut released = [Accumulator::new(0.0); 3];
        add_channels(&mut released, fission_prompt)?;
        add_channels(&mut released, capture_prompt)?;
        let mut generated = Accumulator::new(fission_prompt.total_j());
        generated.add(capture_prompt.total_j())?;
        let mut storage_change = Accumulator::new(0.0);
        let mut energy_error = Accumulator::new(0.0);
        for (index, (group, &initial)) in self
            .source
            .model
            .decay_groups
            .iter()
            .zip(&self.source.delayed_energy_j)
            .enumerate()
        {
            checkpoint()?;
            let input = finite(fissions * group.energy_per_fission.total_j())?;
            generated.add(input)?;
            let decay = finite(group.decay_constant_s_inv * duration_s)?;
            if decay == 0.0 {
                return Err(FissionError::Arithmetic);
            }
            let (stored_input, released_input) = produced_energy(input, exposure, decay)?;
            let final_energy = finite(scaled_decay(initial, decay)? + stored_input)?;
            let output = finite(initial * decayed_fraction(decay) + released_input)?;
            let destinations = partition(output, group.energy_per_fission)?;
            add_channels(&mut released, destinations)?;
            next.source.delayed_energy_j[index] = final_energy;
            energy_error.add(storage_change.add(final_energy - initial)?.abs())?;
            energy_error
                .add(sum(&[final_energy, -initial, destinations.total_j(), -input]).abs())?;
        }
        let released = channel_values(released);
        finite(released.total_j())?;
        for (accumulator, value) in next.source.released.iter_mut().zip([
            released.deposited_j,
            released.radiation_j,
            released.escaped_j,
        ]) {
            energy_error.add(accumulator.add(value)?.abs())?;
        }
        let fission_energy_j = finite(fissions * self.source.model.total_energy_per_fission_j)?;
        let capture_energy_j =
            finite(captures * self.model.data.capture_energy_per_event.total_j())?;
        let generated_energy_j = finite(sum(&[fission_energy_j, capture_energy_j]))?;
        energy_error.add((generated.value() - generated_energy_j).abs())?;
        let stored_delayed_change_j = storage_change.value();
        let energy_residual_j = finite(sum(&[
            stored_delayed_change_j,
            released.total_j(),
            -generated_energy_j,
        ]))?;
        energy_error.add(energy_residual_j.abs())?;
        let energy_roundoff_j = energy_error.value();
        if energy_roundoff_j > budget.source.energy_tolerance_j {
            return Err(FissionError::Budget("energy accounting roundoff"));
        }
        checkpoint()?;
        Ok((
            next,
            IrradiationStep {
                duration_s,
                fissions,
                captures,
                released,
                fission_energy_j,
                capture_energy_j,
                generated_energy_j,
                stored_delayed_change_j,
                energy_residual_j,
                energy_roundoff_j,
                parent_residual_nuclei,
                parent_roundoff_nuclei,
            },
        ))
    }

    /// Publish only a fully accepted trial; cancellation/refusal preserves self.
    pub fn advance(
        &mut self,
        duration_s: f64,
        flux: &GroupFlux,
        temperature_k: f64,
        budget: IrradiationBudget,
        checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<IrradiationStep, FissionError> {
        let (next, report) = self.trial(duration_s, flux, temperature_k, budget, checkpoint)?;
        *self = next;
        Ok(report)
    }
}

/// Integrated expected reactions and their extensive energy accounts.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IrradiationStep {
    pub duration_s: f64,
    pub fissions: f64,
    pub captures: f64,
    pub released: EnergyChannels,
    pub fission_energy_j: f64,
    pub capture_energy_j: f64,
    pub generated_energy_j: f64,
    pub stored_delayed_change_j: f64,
    /// `delta_delayed + released - fission_energy - capture_energy` [J].
    pub energy_residual_j: f64,
    /// Noncancelling local and aggregate energy defects, not model uncertainty [J].
    pub energy_roundoff_j: f64,
    /// `delta_parent + fissions + captures` [expected nuclei].
    pub parent_residual_nuclei: f64,
    pub parent_roundoff_nuclei: f64,
}

fn rate_product(a: f64, b: f64) -> Result<f64, FissionError> {
    let product = finite(a * b)?;
    if product == 0.0 && a > 0.0 && b > 0.0 {
        Err(FissionError::Arithmetic)
    } else {
        Ok(product)
    }
}

fn reaction_share(removed: f64, rate: f64, total: f64) -> Result<f64, FissionError> {
    let fraction = rate / total;
    let result = if fraction < f64::MIN_POSITIVE && rate > 0.0 {
        // A subnormal fraction can lose many bits before multiplication by N.
        (removed / total) * rate
    } else {
        removed * fraction
    };
    finite(result)
}

// The shared deterministic scalar owner preserves representable late tails;
// this physical owner admits the exposure and refuses nonfinite results.
fn scaled_decay(value: f64, exposure: f64) -> Result<f64, FissionError> {
    finite(attenuate(value, exposure))
}

// phi(x) = (1-exp(-x))/x, with a separately stable complement.
fn dimensionless_factors(x: f64) -> (f64, f64) {
    if x < 0.01 {
        let complement = x
            * (0.5
                + x * (-1.0 / 6.0
                    + x * (1.0 / 24.0
                        + x * (-1.0 / 120.0
                            + x * (1.0 / 720.0 + x * (-1.0 / 5040.0 + x / 40320.0))))));
        (1.0 - complement, complement)
    } else {
        let phi = decayed_fraction(x) / x;
        (phi, 1.0 - phi)
    }
}

// For exponential parent depletion x=k*t and delayed decay y=lambda*t,
// retained/input = exp(-min(x,y))*phi(|x-y|)/phi(x). This has the correct
// repeated-rate limit. The released complement needs independent treatment
// when y is tiny; subtracting stored energy from its input would erase heat.
fn produced_energy(input: f64, x: f64, y: f64) -> Result<(f64, f64), FissionError> {
    if input == 0.0 {
        return Ok((0.0, 0.0));
    }
    let (phi_x, complement_x) = dimensionless_factors(x);
    if y < 0.01 {
        let released_fraction = if x < 0.02 {
            // Integral exp(-x*s)*(1-exp(-y*(1-s))) ds on s in [0,1]:
            // y * sum_{n>=0} (-1)^n h_n(x,y)/(n+2)!, where h_n is
            // x^n + x^(n-1)y + ... + y^n. Eight fixed terms have a
            // remainder below f64 resolution on this admitted small domain.
            let (mut homogeneous, mut y_power, mut factorial) = (1.0, 1.0, 2.0);
            let mut series = 0.5;
            for n in 1..=8 {
                y_power *= y;
                homogeneous = x * homogeneous + y_power;
                factorial *= f64::from(n + 2);
                let term = homogeneous / factorial;
                series += if n % 2 == 0 { term } else { -term };
            }
            y * series / phi_x
        } else {
            // y*(phi(y)-phi(x))/((x-y)*phi(x)); complements retain the
            // small numerator while x-y stays at least 0.01 in this branch.
            let (_, complement_y) = dimensionless_factors(y);
            // Reassociate x*phi(x)=1-exp(-x), so a very large x cannot
            // overflow the ratio or underflow y/x before it is rescaled.
            (y * (x / (x - y))) * (complement_x - complement_y) / decayed_fraction(x)
        };
        if !(0.0..=1.0).contains(&released_fraction) {
            return Err(FissionError::Arithmetic);
        }
        Ok((
            finite(input * (1.0 - released_fraction))?,
            finite(input * released_fraction)?,
        ))
    } else {
        let (phi_difference, _) = dimensionless_factors((x - y).abs());
        let log_retained = phi_difference.ln() - phi_x.ln() - x.min(y);
        if !log_retained.is_finite() || log_retained > 0.0 {
            return Err(FissionError::Arithmetic);
        }
        Ok((
            scaled_decay(input, -log_retained)?,
            finite(input * decayed_fraction(-log_retained))?,
        ))
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
#[path = "irradiation/temperature_tests.rs"]
mod temperature_tests;
