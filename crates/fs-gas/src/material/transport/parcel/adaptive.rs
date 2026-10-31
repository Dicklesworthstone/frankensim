//! Whole-window temporal accuracy over the original finite gas snapshots.
//! Only the two-step fine trajectory enters accepted states and accounts.
//! The full/fine difference is a first-order local error estimate, not a
//! certified global error or an extrapolated thermodynamic state.

use super::{
    ClosedGasParcel, GasEnergySolve, GasTransportError, GasTransportPaths, MAX_SNAPSHOT_PARCELS,
    MAX_SNAPSHOT_TRANSFERS, ParcelMolecularCandidate, ParcelMolecularError, ParcelMolecularLink,
    ParcelMolecularSolve, ParcelMolecularStep, ParcelMolecularTransport, SpeciesId, TransportLaw,
    positive,
};
use fs_material::thermochemical_gas::energy::GasEnergyError;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::network::snapshot::SnapshotNozzleError;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::{
    diffusion::network::DiffusionNetworkError,
    network::{
        SpeciesNetworkError,
        transport::{
            GasTransportWindowAudit, audit_closed_window, membrane::MembraneNetworkError,
            porous::PorousNetworkError,
        },
    },
};
use std::collections::BTreeMap;

/// Reuse the native gross species/enthalpy account, including reverse directions.
pub use fs_material::thermochemical_gas::energy::parcel::transfer::species::diffusion::network::adaptive::DiffusionNetworkIntegratedTransfer
    as ParcelMolecularIntegratedTransfer;

/// Local accuracy and work limits; physical budgets remain independently fixed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParcelMolecularAdvanceControls {
    /// Requested comparison-step cap [s], rounded through endpoint clocks.
    /// Accepted states use two substeps with their actual represented durations.
    pub maximum_step_s: f64,
    /// Retry floor [s]; a shorter exact final remainder may still be attempted.
    pub minimum_step_s: f64,
    /// All native full/fine probes, including refusals; in 3..=100,000.
    pub maximum_evaluations: usize,
    /// Per-parcel, per-species absolute inventory error [kg].
    pub species_tolerance_kg: f64,
    /// Per-parcel absolute endpoint-temperature error [K].
    pub temperature_tolerance_k: f64,
    /// Multiplies the largest coarse/fine CHANGE from the accepted state.
    /// In [0,1); absolute temperatures/carrier inventories cannot hide changes.
    pub relative_tolerance: f64,
}

/// Accepted additional energy beyond donor partial enthalpy, grouped by pair.
/// These are the native combined heat corrections, not a second enthalpy flux.
#[derive(Clone, Debug, PartialEq)]
pub struct ParcelMolecularIntegratedHeat {
    pub parcels: [u64; 2],
    pub left_to_right_j: f64,
    pub right_to_left_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ParcelMolecularAdvanceReport {
    pub interval_s: [f64; 2],
    pub accepted_intervals: usize,
    pub accepted_substeps: usize,
    pub rejected_intervals: usize,
    pub evaluations: usize,
    pub minimum_accepted_substep_s: f64,
    pub maximum_accepted_substep_s: f64,
    /// Species and temperature order, over accepted comparisons only.
    pub maximum_accepted_error_ratios: [f64; 2],
    pub suggested_step_s: f64,
    /// All accepted native and continuum species streams, in physical-edge order.
    pub transfers: Vec<ParcelMolecularIntegratedTransfer>,
    pub heat: Vec<ParcelMolecularIntegratedHeat>,
    /// Compensated accepted native energy residual; checked against ONE budget.
    pub accumulated_energy_residual_j: f64,
    /// Independently checked original-to-final source-qualified physical storage.
    pub window: GasTransportWindowAudit,
    /// Detailed source reports for the last accepted half step, not all history.
    pub last_substep: ParcelMolecularStep,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ParcelMolecularAdvanceCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    report: ParcelMolecularAdvanceReport,
}
impl ParcelMolecularAdvanceCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn report(&self) -> &ParcelMolecularAdvanceReport {
        &self.report
    }
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, ParcelMolecularAdvanceReport) {
        (self.parcels, self.report)
    }
}

/// No refusal publishes a partially advanced physical map or caller clock.
#[derive(Clone, Debug, PartialEq)]
pub enum ParcelMolecularAdvanceError {
    Invalid(&'static str),
    Physical {
        time_s: f64,
        step_s: f64,
        error: ParcelMolecularError,
    },
    Accuracy {
        time_s: f64,
        step_s: f64,
        error_ratios: [f64; 2],
    },
    TimeResolution {
        time_s: f64,
        target_s: f64,
    },
    WorkLimit {
        time_s: f64,
        accepted_intervals: usize,
        evaluations: usize,
    },
    Window(GasTransportError),
    Cancelled,
}
impl core::fmt::Display for ParcelMolecularAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive finite molecular gas: {self:?}")
    }
}
impl std::error::Error for ParcelMolecularAdvanceError {}
use ParcelMolecularAdvanceError as Error;

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn physical_error(error: ParcelMolecularError, time_s: f64, step_s: f64) -> Error {
    if error == ParcelMolecularError::Cancelled {
        Error::Cancelled
    } else {
        Error::Physical {
            time_s,
            step_s,
            error,
        }
    }
}
fn retryable(error: &Error) -> bool {
    matches!(
        error,
        Error::Accuracy { .. }
            | Error::Physical {
                error: ParcelMolecularError::Transport(
                    GasTransportError::Exchange(SpeciesNetworkError::Withdrawal { .. })
                        | GasTransportError::Exchange(SpeciesNetworkError::Parcel {
                            error: GasEnergyError::OutsideBracket { .. }
                                | GasEnergyError::NotConverged { .. },
                            ..
                        })
                        | GasTransportError::Nozzle(SnapshotNozzleError::EquilibriumCrossed { .. })
                        | GasTransportError::Porous(PorousNetworkError::EquilibriumCrossed { .. })
                        | GasTransportError::Membrane(
                            MembraneNetworkError::EquilibriumCrossed { .. }
                        )
                        | GasTransportError::Diffusion(
                            DiffusionNetworkError::EntropyDecrease { .. }
                        )
                ),
                ..
            }
    )
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), Error> {
        let (value, lost) = fs_math::eft::two_sum(self.value, term);
        self.value = value;
        self.correction += lost;
        if [term, self.value, self.correction, self.total()]
            .iter()
            .all(|v| v.is_finite())
        {
            Ok(())
        } else {
            Err(Error::Invalid(
                "unrepresentable accepted whole-window account",
            ))
        }
    }
    fn total(&self) -> f64 {
        self.value + self.correction
    }
}
#[derive(Default)]
struct Accepted {
    energy: Sum,
    transfers: BTreeMap<(u64, u64, SpeciesId), [Sum; 2]>,
    heat: BTreeMap<[u64; 2], [Sum; 2]>,
}
impl Accepted {
    fn record(
        &mut self,
        step: &ParcelMolecularStep,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        let native = &step.transport.exchange;
        self.energy.add(native.total_energy_residual_j)?;
        for row in &native.fluxes {
            poll(cancelled)?;
            let key = (
                row.transfer.donor,
                row.transfer.receiver,
                row.transfer.species.clone(),
            );
            if !self.transfers.contains_key(&key)
                && self.transfers.len() + self.heat.len() == MAX_SNAPSHOT_TRANSFERS
            {
                return Err(Error::Invalid(
                    "whole-window transfer accounts exceed the native work bound",
                ));
            }
            let sum = self.transfers.entry(key).or_default();
            sum[0].add(row.transfer.mass_kg)?;
            sum[1].add(row.carried_enthalpy_j)?;
        }
        for row in &native.heat {
            poll(cancelled)?;
            if !self.heat.contains_key(&row.parcels)
                && self.transfers.len() + self.heat.len() == MAX_SNAPSHOT_TRANSFERS
            {
                return Err(Error::Invalid(
                    "whole-window heat accounts exceed the native work bound",
                ));
            }
            let sums = self.heat.entry(row.parcels).or_default();
            sums[usize::from(row.heat_left_to_right_j < 0.0)]
                .add(row.heat_left_to_right_j.abs())?;
        }
        Ok(())
    }
}

fn masses(parcel: &ClosedGasParcel) -> Result<BTreeMap<SpeciesId, f64>, ParcelMolecularError> {
    let model = parcel.model().thermodynamics();
    let molar: Vec<_> = model.components().iter().map(|s| s.molar_mass()).collect();
    let fractions = model
        .mole_composition()
        .to_mass_fractions(&molar)
        .map_err(ParcelMolecularError::Composition)?;
    model
        .components()
        .iter()
        .zip(fractions.fractions())
        .map(|(source, fraction)| {
            let mass = parcel.mass_kg() * fraction;
            if positive(mass) {
                Ok((source.species().clone(), mass))
            } else {
                Err(ParcelMolecularError::Invalid(
                    "positive temporal species inventory is unrepresentable",
                ))
            }
        })
        .collect()
}
fn ratio(
    old: f64,
    coarse: f64,
    fine: f64,
    absolute: f64,
    relative: f64,
) -> Result<f64, ParcelMolecularError> {
    let allowance = absolute + relative * (coarse - old).abs().max((fine - old).abs());
    let ratio = (coarse - fine).abs() / allowance;
    if positive(allowance) && ratio.is_finite() {
        Ok(ratio)
    } else {
        Err(ParcelMolecularError::Invalid(
            "unrepresentable temporal error scale",
        ))
    }
}
fn errors(
    old: &BTreeMap<u64, ClosedGasParcel>,
    coarse: &BTreeMap<u64, ClosedGasParcel>,
    fine: &BTreeMap<u64, ClosedGasParcel>,
    controls: ParcelMolecularAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 2], ParcelMolecularError> {
    let mut result = [0.0_f64; 2];
    for (&id, original) in old {
        super::poll(cancelled)?;
        let (a, b) = (&coarse[&id], &fine[&id]);
        let (m0, ma, mb) = (masses(original)?, masses(a)?, masses(b)?);
        for name in m0.keys().chain(ma.keys()).chain(mb.keys()) {
            result[0] = result[0].max(ratio(
                m0.get(name).copied().unwrap_or(0.0),
                ma.get(name).copied().unwrap_or(0.0),
                mb.get(name).copied().unwrap_or(0.0),
                controls.species_tolerance_kg,
                controls.relative_tolerance,
            )?);
        }
        result[1] = result[1].max(ratio(
            original.state().state().temperature,
            a.state().state().temperature,
            b.state().state().temperature,
            controls.temperature_tolerance_k,
            controls.relative_tolerance,
        )?);
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn probe<C: FnMut() -> bool>(
    state: &BTreeMap<u64, ClosedGasParcel>,
    time_s: f64,
    dt: f64,
    accepted: usize,
    evaluations: &mut usize,
    maximum: usize,
    step: &mut impl FnMut(
        &BTreeMap<u64, ClosedGasParcel>,
        f64,
        &mut C,
    ) -> Result<ParcelMolecularCandidate, ParcelMolecularError>,
    cancelled: &mut C,
) -> Result<ParcelMolecularCandidate, Error> {
    poll(cancelled)?;
    if *evaluations == maximum {
        return Err(Error::WorkLimit {
            time_s,
            accepted_intervals: accepted,
            evaluations: *evaluations,
        });
    }
    *evaluations += 1;
    step(state, dt, cancelled).map_err(|e| physical_error(e, time_s, dt))
}

impl<L: TransportLaw> ParcelMolecularTransport<L> {
    /// Complete a positive represented interval with original full/fine probes.
    /// Fixed faces, native paths and source context apply throughout this window.
    ///
    /// Only typed withdrawal, pressure-order, entropy or caloric bracket/root
    /// timestep refusals and local error excess may halve the comparison span.
    /// Actual material/source-domain refusals remain hard. Every probe retains the native source
    /// requalification and physical gates; only accepted fine steps enter reports.
    /// The original-to-final audit and accumulated native residual use the same
    /// absolute physical budgets without scaling by the number of substeps.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_candidate<C: FnMut() -> bool>(
        &mut self,
        original: &BTreeMap<u64, ClosedGasParcel>,
        links: &[ParcelMolecularLink],
        paths: GasTransportPaths<'_>,
        interval_s: [f64; 2],
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        physical: ParcelMolecularSolve,
        controls: ParcelMolecularAdvanceControls,
        mut cancelled: C,
    ) -> Result<ParcelMolecularAdvanceCandidate, Error> {
        poll(&mut cancelled)?;
        let [start, target] = interval_s;
        if !interval_s.iter().all(|t| t.is_finite())
            || !positive(target - start)
            || ![
                controls.minimum_step_s,
                controls.maximum_step_s,
                controls.species_tolerance_kg,
                controls.temperature_tolerance_k,
            ]
            .iter()
            .all(|v| positive(*v))
            || controls.minimum_step_s > controls.maximum_step_s
            || !(3..=100_000).contains(&controls.maximum_evaluations)
            || !controls.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.relative_tolerance)
            || original.is_empty()
            || original.len() > MAX_SNAPSHOT_PARCELS
            || !original.keys().eq(endpoints.keys())
        {
            return Err(Error::Invalid(
                "forward interval, bounded graph/work, exact endpoint map and positive accuracy scales required",
            ));
        }
        for (&id, solve) in endpoints {
            poll(&mut cancelled)?;
            solve.validate().map_err(|error| {
                physical_error(
                    ParcelMolecularError::Transport(GasTransportError::Exchange(
                        SpeciesNetworkError::Parcel { id, error },
                    )),
                    start,
                    0.0,
                )
            })?;
            if !(solve.minimum_temperature_k..=solve.maximum_temperature_k)
                .contains(&original[&id].state().state().temperature)
            {
                return Err(Error::Invalid(
                    "original gas temperature outside endpoint caloric bracket",
                ));
            }
        }
        let mut step = |state: &BTreeMap<u64, ClosedGasParcel>, dt, cancelled: &mut C| {
            self.step_candidate(state, links, paths, dt, endpoints, physical, cancelled)
        };
        let mut state = original.clone();
        let (mut time, mut ceiling) = (start, controls.maximum_step_s);
        let (mut accepted, mut rejected, mut evaluations) = (0, 0, 0);
        let (mut minimum, mut maximum) = (f64::INFINITY, 0.0_f64);
        let mut maxima = [0.0_f64; 2];
        let mut ledger = Accepted::default();
        let mut last = None;
        while time < target {
            poll(&mut cancelled)?;
            let mut end = if ceiling >= target - time {
                target
            } else {
                time + ceiling
            };
            // Coalesce only an adjacent representable final tail BEFORE either
            // physical trajectory is sampled; never snap an already solved clock.
            if end < target && end.next_up() == target {
                end = target;
            }
            let span = end - time;
            let middle = time + 0.5 * span;
            let (h1, h2) = (middle - time, end - middle);
            if !(time < middle && middle < end && end <= target) || h1 + h2 != span {
                return Err(Error::TimeResolution {
                    time_s: time,
                    target_s: end,
                });
            }
            let attempted = (|| {
                let coarse = probe(
                    &state,
                    time,
                    span,
                    accepted,
                    &mut evaluations,
                    controls.maximum_evaluations,
                    &mut step,
                    &mut cancelled,
                )?;
                let first = probe(
                    &state,
                    time,
                    h1,
                    accepted,
                    &mut evaluations,
                    controls.maximum_evaluations,
                    &mut step,
                    &mut cancelled,
                )?;
                let second = probe(
                    first.parcels(),
                    middle,
                    h2,
                    accepted,
                    &mut evaluations,
                    controls.maximum_evaluations,
                    &mut step,
                    &mut cancelled,
                )?;
                let mut ratios = errors(
                    &state,
                    coarse.parcels(),
                    second.parcels(),
                    controls,
                    &mut cancelled,
                )
                .map_err(|e| physical_error(e, time, span))?;
                // Existing pore-controller convention for unequal represented halves.
                let (u, v) = (h1 / span, h2 / span);
                let factor = (u * u + v * v) / (2.0 * u * v);
                for ratio in &mut ratios {
                    *ratio *= factor;
                }
                if ratios.iter().any(|r| !r.is_finite() || *r > 1.0) {
                    return Err(Error::Accuracy {
                        time_s: time,
                        step_s: span,
                        error_ratios: ratios,
                    });
                }
                Ok((first, second, ratios))
            })();
            match attempted {
                Ok((first, second, ratios)) => {
                    ledger.record(first.step(), &mut cancelled)?;
                    ledger.record(second.step(), &mut cancelled)?;
                    let (next, report) = second.into_parts();
                    state = next;
                    last = Some(report);
                    accepted += 1;
                    minimum = minimum.min(h1.min(h2));
                    maximum = maximum.max(h1.max(h2));
                    for (max, ratio) in maxima.iter_mut().zip(ratios) {
                        *max = max.max(ratio);
                    }
                    time = end;
                    ceiling = if ratios.iter().all(|r| *r <= 0.25) {
                        span + span.min(controls.maximum_step_s - span)
                    } else {
                        span
                    };
                    ceiling = ceiling.clamp(controls.minimum_step_s, controls.maximum_step_s);
                }
                Err(error) => {
                    let smaller = (0.5 * span).max(controls.minimum_step_s);
                    if retryable(&error) && smaller < span {
                        rejected += 1;
                        ceiling = smaller;
                    } else {
                        return Err(error);
                    }
                }
            }
        }
        let window = audit_closed_window(original, &state, physical.physical, &mut cancelled)
            .map_err(|e| {
                if e == GasTransportError::Cancelled {
                    Error::Cancelled
                } else {
                    Error::Window(e)
                }
            })?;
        let energy = ledger.energy.total();
        if energy.abs() > physical.physical.exchange.energy_tolerance_j {
            return Err(Error::Invalid(
                "accepted native residuals exceeded the fixed whole-window energy budget",
            ));
        }
        poll(&mut cancelled)?;
        Ok(ParcelMolecularAdvanceCandidate {
            parcels: state,
            report: ParcelMolecularAdvanceReport {
                interval_s,
                accepted_intervals: accepted,
                accepted_substeps: 2 * accepted,
                rejected_intervals: rejected,
                evaluations,
                minimum_accepted_substep_s: minimum,
                maximum_accepted_substep_s: maximum,
                maximum_accepted_error_ratios: maxima,
                suggested_step_s: ceiling,
                accumulated_energy_residual_j: energy,
                window,
                transfers: ledger
                    .transfers
                    .into_iter()
                    .map(
                        |((donor, receiver, species), sums)| ParcelMolecularIntegratedTransfer {
                            donor,
                            receiver,
                            species,
                            mass_kg: sums[0].total(),
                            carried_enthalpy_j: sums[1].total(),
                        },
                    )
                    .collect(),
                heat: ledger
                    .heat
                    .into_iter()
                    .map(|(parcels, sums)| ParcelMolecularIntegratedHeat {
                        parcels,
                        left_to_right_j: sums[0].total(),
                        right_to_left_j: sums[1].total(),
                    })
                    .collect(),
                last_substep: last
                    .ok_or(Error::Invalid("whole gas window accepted no interval"))?,
            },
        })
    }

    /// Atomically publish the complete proposed gas map and its exact target clock.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to<C: FnMut() -> bool>(
        &mut self,
        parcels: &mut BTreeMap<u64, ClosedGasParcel>,
        time_s: &mut f64,
        target_s: f64,
        links: &[ParcelMolecularLink],
        paths: GasTransportPaths<'_>,
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        physical: ParcelMolecularSolve,
        controls: ParcelMolecularAdvanceControls,
        mut cancelled: C,
    ) -> Result<ParcelMolecularAdvanceReport, Error> {
        let candidate = self.advance_candidate(
            parcels,
            links,
            paths,
            [*time_s, target_s],
            endpoints,
            physical,
            controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        let (next, report) = candidate.into_parts();
        *parcels = next;
        *time_s = target_s;
        Ok(report)
    }
}
