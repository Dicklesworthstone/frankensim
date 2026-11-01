//! Adaptive open molecular networks on the original simultaneous candidates.
//!
//! Only accepted fine steps enter the held-bath account. Bath enthalpy and
//! entropy use the unchanged original reservoir state over the whole window;
//! a per-step relative caloric chart is never integrated as boundary energy.

use super::*;
use crate::material::transport::parcel::adaptive::{
    self as time, ParcelMolecularAdvanceControls, ParcelMolecularAdvanceError as AdvanceError,
    ParcelMolecularIntegratedHeat, ParcelMolecularIntegratedTransfer,
};

/// Independent state and integrated-stream accuracy, with one bounded clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReservoirMolecularAdvanceControls {
    pub state: ParcelMolecularAdvanceControls,
    /// Absolute local error of each directional species stream [kg].
    pub species_transfer_tolerance_kg: f64,
    /// Absolute local error of each directional additional-heat stream [J].
    /// Internal links use the native donor-enthalpy correction; baths use
    /// reduced heat beyond their fixed partial-enthalpy streams.
    pub heat_transfer_tolerance_j: f64,
}

/// Accepted species streams on one fixed held-bath link, in canonical order.
#[derive(Clone, Debug, PartialEq)]
pub struct ReservoirMolecularIntegratedSpecies {
    pub id: u64,
    pub parcel: u64,
    pub species: SpeciesId,
    pub mass_into_gas_kg: f64,
    pub mass_out_of_gas_kg: f64,
    /// Independently compensated signed sum, retaining cancellation accuracy.
    pub net_mass_into_gas_kg: f64,
}

/// Accepted reduced heat beyond the fixed bath's partial-enthalpy streams.
#[derive(Clone, Debug, PartialEq)]
pub struct ReservoirMolecularIntegratedHeat {
    pub id: u64,
    pub parcel: u64,
    pub heat_into_gas_j: f64,
    pub heat_out_of_gas_j: f64,
    pub net_heat_into_gas_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ReservoirMolecularAdvanceReport {
    pub interval_s: [f64; 2],
    pub accepted_intervals: usize,
    pub accepted_substeps: usize,
    pub rejected_intervals: usize,
    pub evaluations: usize,
    pub minimum_accepted_substep_s: f64,
    pub maximum_accepted_substep_s: f64,
    /// Species inventory, temperature, species stream, additional-heat stream.
    pub maximum_accepted_error_ratios: [f64; 4],
    pub suggested_step_s: f64,
    pub transfers: Vec<ParcelMolecularIntegratedTransfer>,
    pub heat: Vec<ParcelMolecularIntegratedHeat>,
    pub reservoir_species: Vec<ReservoirMolecularIntegratedSpecies>,
    pub reservoir_heat: Vec<ReservoirMolecularIntegratedHeat>,
    pub accumulated_energy_residual_j: f64,
    /// Actual original-to-final storage minus the accepted held-bath streams.
    pub window: GasReservoirWindowAudit,
    pub last_substep: ReservoirMolecularStep,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ReservoirMolecularAdvanceCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    report: ReservoirMolecularAdvanceReport,
}
impl ReservoirMolecularAdvanceCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn report(&self) -> &ReservoirMolecularAdvanceReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, ClosedGasParcel>,
        ReservoirMolecularAdvanceReport,
    ) {
        (self.parcels, self.report)
    }
}

#[derive(Default)]
struct Directional {
    gross: [time::Sum; 2],
    net: time::Sum,
}
impl Directional {
    fn add(&mut self, value: f64) -> Result<(), AdvanceError> {
        self.net.add(value)?;
        self.gross[usize::from(value < 0.0)].add(value.abs())
    }
}

#[derive(Default)]
struct Accepted {
    native: time::Accepted,
    species: BTreeMap<(u64, u64, SpeciesId), Directional>,
    heat: BTreeMap<(u64, u64), Directional>,
}
impl Accepted {
    fn record(
        &mut self,
        step: &ReservoirMolecularStep,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AdvanceError> {
        self.native.record_native(&step.exchange, cancelled)?;
        for link in &step.links {
            time::poll(cancelled)?;
            for (species, &mass) in step.species.iter().zip(&link.mass_into_gas_kg) {
                time::poll(cancelled)?;
                if mass != 0.0 {
                    self.species
                        .entry((link.id, link.parcel, species.clone()))
                        .or_default()
                        .add(mass)?;
                }
            }
            self.heat
                .entry((link.id, link.parcel))
                .or_default()
                .add(link.reduced_heat_into_gas_j)?;
        }
        if self.native.entry_count() + self.species.len() + self.heat.len() > MAX_SNAPSHOT_TRANSFERS
        {
            return Err(AdvanceError::Invalid(
                "combined internal and reservoir window accounts exceed the native work bound",
            ));
        }
        Ok(())
    }
}

// Compare directional streams individually, so steady throughflow or reversal
// cannot disappear behind a small endpoint change or a cancellation of faces.
#[derive(Default)]
struct Streams {
    species: BTreeMap<(u8, u64, SpeciesId), [time::Sum; 2]>,
    heat: BTreeMap<(u8, u64), [time::Sum; 2]>,
}
impl Streams {
    fn record(
        &mut self,
        step: &ReservoirMolecularStep,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AdvanceError> {
        for (kind, id, masses, heat) in step
            .internal_links
            .iter()
            .map(|link| {
                (
                    0,
                    link.id,
                    &link.mass_left_to_right_kg,
                    link.heat_correction_left_to_right_j,
                )
            })
            .chain(step.links.iter().map(|link| {
                (
                    1,
                    link.id,
                    &link.mass_into_gas_kg,
                    link.reduced_heat_into_gas_j,
                )
            }))
        {
            time::poll(cancelled)?;
            for (species, &mass) in step.species.iter().zip(masses) {
                time::poll(cancelled)?;
                self.species.entry((kind, id, species.clone())).or_default()
                    [usize::from(mass < 0.0)]
                .add(mass.abs())?;
            }
            self.heat.entry((kind, id)).or_default()[usize::from(heat < 0.0)].add(heat.abs())?;
        }
        if self.species.len() + self.heat.len() > MAX_SNAPSHOT_TRANSFERS {
            return Err(AdvanceError::Invalid(
                "temporal stream comparison exceeds the native work bound",
            ));
        }
        Ok(())
    }
}

fn stream_errors(
    coarse: &ReservoirMolecularStep,
    first: &ReservoirMolecularStep,
    second: &ReservoirMolecularStep,
    controls: ReservoirMolecularAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 2], AdvanceError> {
    let (mut a, mut b) = (Streams::default(), Streams::default());
    a.record(coarse, cancelled)?;
    b.record(first, cancelled)?;
    b.record(second, cancelled)?;
    fn compare<K: Ord>(
        a: &BTreeMap<K, [time::Sum; 2]>,
        b: &BTreeMap<K, [time::Sum; 2]>,
        absolute: f64,
        relative: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, AdvanceError> {
        let mut maximum = 0.0_f64;
        for key in a.keys().chain(b.keys()) {
            time::poll(cancelled)?;
            for direction in 0..2 {
                let x = a.get(key).map_or(0.0, |s| s[direction].total());
                let y = b.get(key).map_or(0.0, |s| s[direction].total());
                let allowance = absolute + relative * x.abs().max(y.abs());
                let error = (x - y).abs() / allowance;
                if !positive(allowance) || !error.is_finite() {
                    return Err(AdvanceError::Invalid(
                        "unrepresentable temporal stream error scale",
                    ));
                }
                maximum = maximum.max(error);
            }
        }
        Ok(maximum)
    }
    Ok([
        compare(
            &a.species,
            &b.species,
            controls.species_transfer_tolerance_kg,
            controls.state.relative_tolerance,
            cancelled,
        )?,
        compare(
            &a.heat,
            &b.heat,
            controls.heat_transfer_tolerance_j,
            controls.state.relative_tolerance,
            cancelled,
        )?,
    ])
}

impl<L: TransportLaw> ParcelMolecularTransport<L> {
    /// Evolve internal continuum links and fixed held baths on one native clock.
    ///
    /// Full/fine comparison controls local state and directional stream errors.
    /// Only the fine path is accepted; it is never extrapolated. The original
    /// source domains and gross withdrawals remain native. Final mass, species,
    /// energy and bath-inclusive entropy use fixed whole-window budgets.
    /// This is a local temporal estimator, not a certified global error bound.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_network_candidate<C: FnMut() -> bool>(
        &mut self,
        original: &BTreeMap<u64, ClosedGasParcel>,
        internal_links: &[ParcelMolecularLink],
        reservoir_links: &[ReservoirMolecularLink<'_>],
        interval_s: [f64; 2],
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        physical: ParcelMolecularSolve,
        controls: ReservoirMolecularAdvanceControls,
        mut cancelled: C,
    ) -> Result<ReservoirMolecularAdvanceCandidate, AdvanceError> {
        time::poll(&mut cancelled)?;
        if !positive(controls.species_transfer_tolerance_kg)
            || !positive(controls.heat_transfer_tolerance_j)
        {
            return Err(AdvanceError::Invalid(
                "positive finite species and heat stream accuracy scales required",
            ));
        }
        let n = self.face.gas().species_count();
        if reservoir_links
            .len()
            .checked_mul(n + 2)
            .and_then(|count| {
                internal_links
                    .len()
                    .checked_mul(n + 1)
                    .and_then(|internal| count.checked_add(internal))
            })
            .is_none_or(|count| count > MAX_SNAPSHOT_TRANSFERS)
        {
            return Err(AdvanceError::Invalid(
                "combined continuum and reservoir graph exceeds the native work bound",
            ));
        }
        let mut held = BTreeMap::new();
        for link in reservoir_links {
            time::poll(&mut cancelled)?;
            if held.insert(link.id, link).is_some() {
                return Err(AdvanceError::Invalid("duplicate held-bath link ID"));
            }
        }
        let mut accepted = Accepted::default();
        let temporal = time::run(
            original,
            interval_s,
            endpoints,
            controls.state,
            &mut |state, dt, cancel: &mut C| {
                self.network_candidate(
                    state,
                    internal_links,
                    reservoir_links,
                    dt,
                    endpoints,
                    physical,
                    cancel,
                )
                .map(ReservoirMolecularCandidate::into_parts)
            },
            &mut |step, cancel: &mut C| accepted.record(step, cancel),
            &mut |coarse, first, second, cancel: &mut C| {
                stream_errors(coarse, first, second, controls, cancel)
            },
            &mut cancelled,
        )?;
        let mut by_bath: BTreeMap<u64, Vec<ReservoirSpeciesTransfer>> = BTreeMap::new();
        for ((id, _, species), sums) in &accepted.species {
            time::poll(&mut cancelled)?;
            by_bath
                .entry(*id)
                .or_default()
                .push(ReservoirSpeciesTransfer {
                    species: species.clone(),
                    mass_into_gas_kg: sums.net.total(),
                });
        }
        let exchanges: Vec<_> = held
            .values()
            .map(|link| GasReservoirExchange {
                parcel: link.parcel,
                catalog: link.reservoir.model(),
                state: link.reservoir.state(),
                species: by_bath.get(&link.id).map_or(&[], Vec::as_slice),
                reduced_heat_into_gas_j: accepted
                    .heat
                    .get(&(link.id, link.parcel))
                    .map_or(0.0, |s| s.net.total()),
            })
            .collect();
        let window = transport::audit_reservoir_window(
            original,
            &temporal.parcels,
            &exchanges,
            physical.physical,
            &mut cancelled,
        )
        .map_err(time::window_error)?;
        let (energy, transfers, heat) = accepted.native.into_parts();
        if energy.abs() > physical.physical.exchange.energy_tolerance_j {
            return Err(AdvanceError::Invalid(
                "accepted native residuals exceeded the fixed whole-window energy budget",
            ));
        }
        let reservoir_species = accepted
            .species
            .into_iter()
            .map(
                |((id, parcel, species), sums)| ReservoirMolecularIntegratedSpecies {
                    id,
                    parcel,
                    species,
                    mass_into_gas_kg: sums.gross[0].total(),
                    mass_out_of_gas_kg: sums.gross[1].total(),
                    net_mass_into_gas_kg: sums.net.total(),
                },
            )
            .collect();
        let reservoir_heat = accepted
            .heat
            .into_iter()
            .map(|((id, parcel), sums)| ReservoirMolecularIntegratedHeat {
                id,
                parcel,
                heat_into_gas_j: sums.gross[0].total(),
                heat_out_of_gas_j: sums.gross[1].total(),
                net_heat_into_gas_j: sums.net.total(),
            })
            .collect();
        time::poll(&mut cancelled)?;
        Ok(ReservoirMolecularAdvanceCandidate {
            parcels: temporal.parcels,
            report: ReservoirMolecularAdvanceReport {
                interval_s,
                accepted_intervals: temporal.accepted_intervals,
                accepted_substeps: 2 * temporal.accepted_intervals,
                rejected_intervals: temporal.rejected_intervals,
                evaluations: temporal.evaluations,
                minimum_accepted_substep_s: temporal.minimum_accepted_substep_s,
                maximum_accepted_substep_s: temporal.maximum_accepted_substep_s,
                maximum_accepted_error_ratios: temporal.maximum_accepted_error_ratios,
                suggested_step_s: temporal.suggested_step_s,
                transfers,
                heat,
                reservoir_species,
                reservoir_heat,
                accumulated_energy_residual_j: energy,
                window,
                last_substep: temporal.last_substep,
            },
        })
    }

    /// Publish the complete gas map and exact represented clock together.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_network_to<C: FnMut() -> bool>(
        &mut self,
        parcels: &mut BTreeMap<u64, ClosedGasParcel>,
        time_s: &mut f64,
        target_s: f64,
        internal_links: &[ParcelMolecularLink],
        reservoir_links: &[ReservoirMolecularLink<'_>],
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        physical: ParcelMolecularSolve,
        controls: ReservoirMolecularAdvanceControls,
        mut cancelled: C,
    ) -> Result<ReservoirMolecularAdvanceReport, AdvanceError> {
        let candidate = self.advance_network_candidate(
            parcels,
            internal_links,
            reservoir_links,
            [*time_s, target_s],
            endpoints,
            physical,
            controls,
            &mut cancelled,
        )?;
        time::poll(&mut cancelled)?;
        let (next, report) = candidate.into_parts();
        *parcels = next;
        *time_s = target_s;
        Ok(report)
    }
}
