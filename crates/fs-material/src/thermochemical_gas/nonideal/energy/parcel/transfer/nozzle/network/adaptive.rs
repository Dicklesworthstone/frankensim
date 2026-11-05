//! Adaptive evolution of the original pressure-driven real-gas network.
//!
//! Compare one original interval with two original half intervals, and accept
//! only the conservative fine trajectory. No extrapolated parcel, clipped
//! stream, frozen-pressure history or second caloric solver is introduced.
//! Energy error uses U - sum_i m_i u_i^ideal(T_start), evaluated with finite
//! NASA increments and the retained PR departure. Each species' arbitrary
//! formation datum cancels BEFORE arithmetic, including when mixtures change.
//! Species inventories and temperature have independent local error controls.
//! Estimates concern temporal discretization, not model error or global bounds.

use super::{
    GasEnergyError, GasEnergySolve, IsochoricRealGasParcel, RealGasNozzleNetworkCandidate,
    RealGasNozzleNetworkError, RealGasNozzleNetworkLink, RealGasNozzleNetworkStep,
    SpeciesNetworkError, SpeciesNetworkSolve, Sum, advance, inventory, positive,
};
use fs_qty::Temperature;
use std::collections::{BTreeMap, BTreeSet};

/// Bounded physical work and local accuracy, separate from endpoint root gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasNozzleNetworkTimeSolve {
    /// Full comparison spans; the accepted physical half steps are shorter.
    pub initial_dt_s: f64,
    pub minimum_dt_s: f64,
    pub maximum_dt_s: f64,
    /// Includes inaccurate and infeasible trials, at most 65,536.
    pub maximum_trials: usize,
    /// Two slots per accepted comparison, at most 131,072.
    pub maximum_steps: usize,
    /// Absolute tolerance for each node's reference-cancelled storage [J].
    pub energy_absolute_tolerance_j: f64,
    /// Absolute tolerance for total AND individual species inventories [kg].
    pub mass_absolute_tolerance_kg: f64,
    pub temperature_absolute_tolerance_k: f64,
    /// Scales CHANGES since the accepted state, never absolute energy datums,
    /// original inventory magnitudes or absolute Kelvin. Must lie in [0,1).
    pub relative_tolerance: f64,
    /// Whole-window sum of absolute accepted local caloric residuals [J].
    pub cumulative_energy_tolerance_j: f64,
    /// Whole-window total and individual species conservation tolerance [kg].
    pub cumulative_mass_tolerance_kg: f64,
}

/// Coarse/fine discrepancies before the first-order estimator correction.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RealGasNozzleNetworkTimeError {
    pub maximum_energy_difference_j: f64,
    pub maximum_mass_difference_kg: f64,
    pub maximum_temperature_difference_k: f64,
    /// Corrected for unequal representable half spans. At most one is accepted.
    pub normalized: f64,
}

/// Accepted flow integrated along the two canonical directions of one path.
/// Enthalpy retains the original source convention and may be negative.
/// Separate gross accounts preserve physical recirculation and flow reversals.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasNozzlePathIntegral {
    pub id: u64,
    pub parcels: [u64; 2],
    pub mass_left_to_right_kg: f64,
    pub mass_right_to_left_kg: f64,
    pub enthalpy_left_to_right_j: f64,
    pub enthalpy_right_to_left_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RealGasNozzleNetworkEvolution {
    pub duration_s: f64,
    pub comparisons: usize,
    /// Includes failed or rejected original interval calls.
    pub evaluated_steps: usize,
    /// Only accepted fine half steps contribute flow or caloric accounts.
    pub accepted_steps: usize,
    pub accuracy_rejections: usize,
    pub feasibility_rejections: usize,
    pub maximum_accepted_normalized_error: f64,
    pub last_error: Option<RealGasNozzleNetworkTimeError>,
    pub suggested_dt_s: f64,
    pub paths: Vec<RealGasNozzlePathIntegral>,
    /// Sum of original local caloric residuals on the accepted trajectory.
    pub energy_residual_j: f64,
    pub accumulated_absolute_energy_residual_j: f64,
    pub total_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RealGasNozzleNetworkEvolutionCandidate {
    parcels: BTreeMap<u64, IsochoricRealGasParcel>,
    report: RealGasNozzleNetworkEvolution,
}
impl RealGasNozzleNetworkEvolutionCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, IsochoricRealGasParcel> {
        &self.parcels
    }
    pub const fn report(&self) -> &RealGasNozzleNetworkEvolution {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, IsochoricRealGasParcel>,
        RealGasNozzleNetworkEvolution,
    ) {
        (self.parcels, self.report)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum RealGasNozzleNetworkTimeFailure {
    Invalid(&'static str),
    Network(RealGasNozzleNetworkError),
    Budget(&'static str),
    Accuracy {
        time_s: f64,
        span_s: f64,
        error: RealGasNozzleNetworkTimeError,
    },
}
impl From<RealGasNozzleNetworkError> for RealGasNozzleNetworkTimeFailure {
    fn from(error: RealGasNozzleNetworkError) -> Self {
        Self::Network(error)
    }
}
impl core::fmt::Display for RealGasNozzleNetworkTimeFailure {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "real-gas nozzle evolution: {self:?}")
    }
}
impl std::error::Error for RealGasNozzleNetworkTimeFailure {}
use RealGasNozzleNetworkTimeFailure as Error;

fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    super::check(cancelled).map_err(Into::into)
}

fn validate(duration: f64, controls: RealGasNozzleNetworkTimeSolve) -> Result<(), Error> {
    if !positive(duration)
        || ![
            controls.initial_dt_s,
            controls.minimum_dt_s,
            controls.maximum_dt_s,
            controls.energy_absolute_tolerance_j,
            controls.mass_absolute_tolerance_kg,
            controls.temperature_absolute_tolerance_k,
            controls.cumulative_energy_tolerance_j,
            controls.cumulative_mass_tolerance_kg,
        ]
        .into_iter()
        .all(positive)
        || controls.minimum_dt_s > controls.initial_dt_s
        || controls.initial_dt_s > controls.maximum_dt_s
        || !controls.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&controls.relative_tolerance)
        || !(1..=65_536).contains(&controls.maximum_trials)
        || !(2..=131_072).contains(&controls.maximum_steps)
    {
        return Err(Error::Invalid(
            "positive bounded time, ordered spans, absolute errors and finite relative tolerance in [0,1) required",
        ));
    }
    Ok(())
}

fn retryable(error: &RealGasNozzleNetworkError) -> bool {
    matches!(
        error,
        RealGasNozzleNetworkError::Withdrawal { .. }
            | RealGasNozzleNetworkError::EquilibriumCrossing { .. }
            | RealGasNozzleNetworkError::Exchange(SpeciesNetworkError::Withdrawal { .. })
    )
}

// Retain the original EOS departure already present in the accepted material;
// comparison needs no new cubic, partial-caloric query or material resolution.
fn relative_energy(parcel: &IsochoricRealGasParcel, anchor_t: f64) -> Result<f64, GasEnergyError> {
    let ideal = parcel
        .model()
        .source()
        .thermodynamics()
        .specific_internal_energy_change(
            Temperature::new(anchor_t),
            Temperature::new(parcel.state().temperature_k()),
        )
        .map_err(|error| crate::thermochemical_gas::failure("network time energy", error))?
        .value();
    let residual = parcel.state().root().residual_internal_energy_j_mol()
        / parcel.model().mechanics().molar_mass().value();
    let value = parcel.mass_kg() * (ideal + residual);
    if !value.is_finite() {
        return Err(GasEnergyError::Invalid(
            "network comparison energy is unrepresentable",
        ));
    }
    Ok(value)
}

fn ratio(old: f64, coarse: f64, fine: f64, absolute: f64, relative: f64) -> Result<f64, Error> {
    let scale = absolute + relative * (coarse - old).abs().max((fine - old).abs());
    let value = (fine - coarse).abs() / scale;
    if !positive(scale) || !value.is_finite() {
        return Err(Error::Invalid(
            "network time-error scale is unrepresentable",
        ));
    }
    Ok(value)
}

#[allow(clippy::too_many_arguments)]
fn compare(
    old: &BTreeMap<u64, IsochoricRealGasParcel>,
    coarse: &BTreeMap<u64, IsochoricRealGasParcel>,
    fine: &BTreeMap<u64, IsochoricRealGasParcel>,
    first_s: f64,
    second_s: f64,
    controls: RealGasNozzleNetworkTimeSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RealGasNozzleNetworkTimeError, Error> {
    let a = first_s / (first_s + second_s);
    let b = second_s / (first_s + second_s);
    // Euler local defects scale as dt^2. Fine/(coarse-fine) =
    // (a^2+b^2)/(2*a*b); equal halves recover the familiar factor one.
    let correction = (a * a + b * b) / (2.0 * a * b);
    if !positive(correction) {
        return Err(Error::Invalid(
            "unrepresentable half-interval error estimate",
        ));
    }
    let mut error = RealGasNozzleNetworkTimeError::default();
    for (&id, before) in old {
        check(cancelled)?;
        let c = &coarse[&id];
        let f = &fine[&id];
        let t0 = before.state().temperature_k();
        let tc = c.state().temperature_k();
        let tf = f.state().temperature_k();
        let material = |e| {
            Error::Network(RealGasNozzleNetworkError::Exchange(
                SpeciesNetworkError::Parcel { id, error: e },
            ))
        };
        let u0 = relative_energy(before, t0).map_err(material)?;
        let uc = relative_energy(c, t0).map_err(material)?;
        let uf = relative_energy(f, t0).map_err(material)?;
        error.maximum_energy_difference_j = error.maximum_energy_difference_j.max((uf - uc).abs());
        error.maximum_temperature_difference_k =
            error.maximum_temperature_difference_k.max((tf - tc).abs());
        error.normalized = error
            .normalized
            .max(ratio(
                u0,
                uc,
                uf,
                controls.energy_absolute_tolerance_j,
                controls.relative_tolerance,
            )?)
            .max(ratio(
                t0,
                tc,
                tf,
                controls.temperature_absolute_tolerance_k,
                controls.relative_tolerance,
            )?);
        let initial = inventory(before.model().source(), before.mass_kg()).map_err(material)?;
        let coarse = inventory(c.model().source(), c.mass_kg()).map_err(material)?;
        let fine = inventory(f.model().source(), f.mass_kg()).map_err(material)?;
        error.maximum_mass_difference_kg = error
            .maximum_mass_difference_kg
            .max((f.mass_kg() - c.mass_kg()).abs());
        error.normalized = error.normalized.max(ratio(
            before.mass_kg(),
            c.mass_kg(),
            f.mass_kg(),
            controls.mass_absolute_tolerance_kg,
            controls.relative_tolerance,
        )?);
        for species in initial
            .keys()
            .chain(coarse.keys())
            .chain(fine.keys())
            .collect::<BTreeSet<_>>()
        {
            let m0 = initial.get(species).copied().unwrap_or(0.0);
            let mc = coarse.get(species).copied().unwrap_or(0.0);
            let mf = fine.get(species).copied().unwrap_or(0.0);
            error.maximum_mass_difference_kg =
                error.maximum_mass_difference_kg.max((mf - mc).abs());
            error.normalized = error.normalized.max(ratio(
                m0,
                mc,
                mf,
                controls.mass_absolute_tolerance_kg,
                controls.relative_tolerance,
            )?);
        }
    }
    error.normalized *= correction;
    if !error.normalized.is_finite() {
        return Err(Error::Invalid("network time error is unrepresentable"));
    }
    Ok(error)
}

#[derive(Default)]
struct PathAccount {
    parcels: [u64; 2],
    mass: [Sum; 2],
    enthalpy: [Sum; 2],
}

fn account(
    step: &RealGasNozzleNetworkStep,
    paths: &mut BTreeMap<u64, PathAccount>,
    signed_energy: &mut Sum,
    absolute_energy: &mut Sum,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    for row in &step.links {
        check(cancelled)?;
        let endpoints = [row.donor.min(row.receiver), row.donor.max(row.receiver)];
        let direction = usize::from(row.donor != endpoints[0]);
        let path = paths.entry(row.id).or_insert_with(|| PathAccount {
            parcels: endpoints,
            ..PathAccount::default()
        });
        path.mass[direction].add(row.transferred_mass_kg)?;
        path.enthalpy[direction].add(row.transferred_mass_kg * row.flow.upstream_enthalpy_j_kg)?;
    }
    for parcel in &step.exchange.parcels {
        check(cancelled)?;
        signed_energy.add(parcel.energy_residual_j)?;
        absolute_energy.add(parcel.energy_residual_j.abs())?;
    }
    Ok(())
}

fn conservation(
    original: &BTreeMap<u64, IsochoricRealGasParcel>,
    final_state: &BTreeMap<u64, IsochoricRealGasParcel>,
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64), Error> {
    let mut mass = Sum::default();
    let mut species = BTreeMap::<_, Sum>::new();
    for (&id, initial) in original {
        check(cancelled)?;
        let final_node = &final_state[&id];
        mass.add(final_node.mass_kg() - initial.mass_kg())?;
        let material = |error| {
            Error::Network(RealGasNozzleNetworkError::Exchange(
                SpeciesNetworkError::Parcel { id, error },
            ))
        };
        let a = inventory(initial.model().source(), initial.mass_kg()).map_err(material)?;
        let b = inventory(final_node.model().source(), final_node.mass_kg()).map_err(material)?;
        for id in a.keys().chain(b.keys()).collect::<BTreeSet<_>>() {
            species
                .entry(id.clone())
                .or_default()
                .add(b.get(id).copied().unwrap_or(0.0) - a.get(id).copied().unwrap_or(0.0))?;
        }
    }
    let maximum = species.values().map(|s| s.value.abs()).fold(0.0, f64::max);
    if mass.value.abs() > tolerance || maximum > tolerance {
        return Err(Error::Budget(
            "whole-window conserved species or total mass",
        ));
    }
    Ok((mass.value, maximum))
}

/// Evolve an unforced, closed pressure network for one physical duration.
///
/// All rates and donor enthalpies are resampled on every ORIGINAL interval.
/// Only original inventory-limit and pressure-crossing refusals trigger
/// subdivision; invalid inputs, unsupported sources, caloric/root failures and
/// cancellation propagate. Failed comparisons contribute work but no streams.
/// The full window is transactional, including its final conservation audit.
/// A final span below the floor may be tried once; it still must pass accuracy.
///
/// Work is at most three bounded network calls per trial. Only three candidate
/// maps and one aggregate per physical path are retained, independent of the
/// number of accepted steps. The caller receives the next suggested span for
/// continuation. No latent-phase flash, momentum, externally forced heat/species
/// history or global error enclosure is supplied by this closed-network owner.
#[allow(clippy::too_many_arguments)]
pub fn advance_adaptive(
    parcels: &BTreeMap<u64, IsochoricRealGasParcel>,
    links: &[RealGasNozzleNetworkLink],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    network: SpeciesNetworkSolve,
    controls: RealGasNozzleNetworkTimeSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RealGasNozzleNetworkEvolutionCandidate, Error> {
    check(&mut cancelled)?;
    validate(duration_s, controls)?;
    // The original owner performs full input/source admission before publishing
    // any trial. Retain its parcel/row bounds before making the working clone.
    if parcels.is_empty()
        || parcels.len() > super::MAX_SNAPSHOT_PARCELS
        || links.len() > super::MAX_SNAPSHOT_TRANSFERS
    {
        return Err(Error::Invalid(
            "adaptive nozzle network exceeds original work bounds",
        ));
    }
    let mut next = parcels.clone();
    let mut time = 0.0;
    let mut span = controls.initial_dt_s;
    let mut paths = BTreeMap::new();
    let mut signed_energy = Sum::default();
    let mut absolute_energy = Sum::default();
    let mut report = RealGasNozzleNetworkEvolution {
        duration_s,
        comparisons: 0,
        evaluated_steps: 0,
        accepted_steps: 0,
        accuracy_rejections: 0,
        feasibility_rejections: 0,
        maximum_accepted_normalized_error: 0.0,
        last_error: None,
        suggested_dt_s: span,
        paths: Vec::new(),
        energy_residual_j: 0.0,
        accumulated_absolute_energy_residual_j: 0.0,
        total_mass_residual_kg: 0.0,
        maximum_species_mass_residual_kg: 0.0,
    };
    while time < duration_s {
        check(&mut cancelled)?;
        if report.comparisons == controls.maximum_trials {
            return Err(Error::Budget("adaptive nozzle comparisons"));
        }
        if controls.maximum_steps - report.accepted_steps < 2 {
            return Err(Error::Budget("accepted nozzle half steps"));
        }
        let end = if span >= duration_s - time {
            duration_s
        } else {
            time + span
        };
        let middle = time + 0.5 * (end - time);
        if !(time < middle && middle < end) {
            return Err(Error::Invalid(
                "adaptive nozzle half times are unrepresentable",
            ));
        }
        let whole = end - time;
        let first_dt = middle - time;
        let second_dt = end - middle;
        report.comparisons += 1;
        let mut run = |state: &BTreeMap<u64, IsochoricRealGasParcel>, dt| {
            report.evaluated_steps += 1;
            advance(state, links, dt, endpoints, network, &mut cancelled)
        };
        let trial = (|| -> Result<_, RealGasNozzleNetworkError> {
            let coarse = run(&next, whole)?;
            let first = run(&next, first_dt)?;
            let second = run(first.parcels(), second_dt)?;
            Ok((coarse, first, second))
        })();
        let (coarse, first, second): (
            RealGasNozzleNetworkCandidate,
            RealGasNozzleNetworkCandidate,
            RealGasNozzleNetworkCandidate,
        ) = match trial {
            Ok(value) => value,
            Err(error) if retryable(&error) && whole > controls.minimum_dt_s => {
                report.feasibility_rejections += 1;
                span = (0.5 * whole).max(controls.minimum_dt_s);
                continue;
            }
            Err(error) => return Err(error.into()),
        };
        let estimate = compare(
            &next,
            coarse.parcels(),
            second.parcels(),
            first_dt,
            second_dt,
            controls,
            &mut cancelled,
        )?;
        report.last_error = Some(estimate);
        if estimate.normalized > 1.0 {
            report.accuracy_rejections += 1;
            if whole <= controls.minimum_dt_s {
                return Err(Error::Accuracy {
                    time_s: time,
                    span_s: whole,
                    error: estimate,
                });
            }
            span = (whole * (0.9 / estimate.normalized.sqrt()).clamp(0.1, 0.5))
                .max(controls.minimum_dt_s);
            continue;
        }
        account(
            first.step(),
            &mut paths,
            &mut signed_energy,
            &mut absolute_energy,
            &mut cancelled,
        )?;
        account(
            second.step(),
            &mut paths,
            &mut signed_energy,
            &mut absolute_energy,
            &mut cancelled,
        )?;
        if absolute_energy.value > controls.cumulative_energy_tolerance_j {
            return Err(Error::Budget("whole-window accumulated caloric residual"));
        }
        report.accepted_steps += 2;
        report.maximum_accepted_normalized_error = report
            .maximum_accepted_normalized_error
            .max(estimate.normalized);
        next = second.into_parts().0;
        time = end;
        let growth = if estimate.normalized == 0.0 {
            2.0
        } else {
            (0.9 / estimate.normalized.sqrt()).clamp(0.5, 2.0)
        };
        span = (whole * growth).clamp(controls.minimum_dt_s, controls.maximum_dt_s);
    }
    let (mass, species) = conservation(
        parcels,
        &next,
        controls.cumulative_mass_tolerance_kg,
        &mut cancelled,
    )?;
    report.total_mass_residual_kg = mass;
    report.maximum_species_mass_residual_kg = species;
    report.energy_residual_j = signed_energy.value;
    report.accumulated_absolute_energy_residual_j = absolute_energy.value;
    report.suggested_dt_s = span;
    report.paths = paths
        .into_iter()
        .map(|(id, path)| RealGasNozzlePathIntegral {
            id,
            parcels: path.parcels,
            mass_left_to_right_kg: path.mass[0].value,
            mass_right_to_left_kg: path.mass[1].value,
            enthalpy_left_to_right_j: path.enthalpy[0].value,
            enthalpy_right_to_left_j: path.enthalpy[1].value,
        })
        .collect();
    check(&mut cancelled)?;
    Ok(RealGasNozzleNetworkEvolutionCandidate {
        parcels: next,
        report,
    })
}

#[cfg(test)]
mod tests;
