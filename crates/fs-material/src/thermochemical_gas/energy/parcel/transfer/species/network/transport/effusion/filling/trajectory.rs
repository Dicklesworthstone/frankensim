//! Adaptive transition from exact vacuum to ordinary reciprocal effusion.
//!
//! The first interval compares one full capture against a half capture followed
//! by a reciprocal half step. Receiver backflow therefore participates in the
//! error estimate immediately. Only the refined pair is accepted. The existing
//! network trajectory then owns every remaining interval; no duplicate rate,
//! thermochemistry, pressure floor or continuing vacuum-outlet model is used.
//!
//! This retains the parent owner's thermalized, fixed rigid-volume assumption.
//! It does not resolve the initially non-Maxwellian molecular transit. Error
//! estimates are local; source-supported endpoint pressures can limit refinement.

use super::super::{
    EffusionNetworkLink, face_error,
    reservoir::{self, ReservoirEffusionCandidate, ReservoirEffusionStep},
    trajectory::{
        self as native, EffusionAdvanceControls, EffusionAdvanceError, EffusionAdvanceStep,
    },
};
use super::{
    EffusionFillBalance, EffusionFillCandidate, EffusionFillError, EffusionFillLink,
    EffusionFillStep, EffusionReceiver, audit_pair, step_candidate,
};
use crate::thermochemical_gas::energy::parcel::transfer::species::effusion::positive_product_ratio;
use crate::thermochemical_gas::energy::{
    GasEnergySolve,
    parcel::{
        ClosedGasParcel,
        evacuated::EvacuatedGasVolume,
        transfer::species::network::transport::{
            GasTransportError, GasTransportPaths, GasTransportSolve,
        },
    },
};
use fs_math::dd::Dd;
use fs_qty::Temperature;
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use std::collections::BTreeMap;

#[derive(Clone, Debug, PartialEq)]
pub enum EffusionFillAdvanceError {
    Invalid(&'static str),
    /// The original incident beam and the receiver's temperature bracket imply
    /// a smaller initial step than the caller permits. No source was extrapolated.
    InitialPressureLimit {
        maximum_step_s: f64,
        minimum_step_s: f64,
    },
    Fill {
        duration_s: f64,
        error: EffusionFillError,
    },
    /// Native errors use elapsed time and work counts for the COMPLETE window.
    Advance(EffusionAdvanceError),
    Cancelled,
}
impl core::fmt::Display for EffusionFillAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive finite molecular filling: {self:?}")
    }
}
impl std::error::Error for EffusionFillAdvanceError {}
use EffusionFillAdvanceError as Error;

#[derive(Clone, Debug, PartialEq)]
pub struct EffusionFillAdvanceStep {
    pub duration_s: f64,
    /// Initial full/half probes plus all evaluations in the native continuation.
    pub attempted_steps: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub maximum_accepted_error: f64,
    /// Accepted first half of the initial interval. For an exactly inactive
    /// aperture this is the one unchanged full-window transaction instead.
    pub first_fill: EffusionFillStep,
    /// Accepted second half, with BOTH real finite gases and reciprocal beams.
    pub first_reciprocal: Option<ReservoirEffusionStep>,
    pub continuation: Option<EffusionAdvanceStep>,
    /// Independent original donor + vacuum to final pair mass/energy/entropy.
    pub balance: EffusionFillBalance,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffusionFillAdvanceCandidate {
    donor: ClosedGasParcel,
    receiver: EffusionReceiver,
    step: EffusionFillAdvanceStep,
}
impl EffusionFillAdvanceCandidate {
    pub const fn donor(&self) -> &ClosedGasParcel {
        &self.donor
    }
    pub const fn receiver(&self) -> &EffusionReceiver {
        &self.receiver
    }
    pub const fn step(&self) -> &EffusionFillAdvanceStep {
        &self.step
    }
    pub fn into_parts(self) -> (ClosedGasParcel, EffusionReceiver, EffusionFillAdvanceStep) {
        (self.donor, self.receiver, self.step)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn advance_error(error: EffusionAdvanceError) -> Error {
    if error == EffusionAdvanceError::Cancelled {
        Error::Cancelled
    } else {
        Error::Advance(error)
    }
}
fn fill_error(error: EffusionFillError, duration_s: f64) -> Error {
    if error == EffusionFillError::Cancelled {
        Error::Cancelled
    } else {
        Error::Fill { duration_s, error }
    }
}
fn native_error(error: GasTransportError, elapsed_s: f64, duration_s: f64) -> Error {
    if error == GasTransportError::Cancelled {
        Error::Cancelled
    } else {
        advance_error(EffusionAdvanceError::Step {
            elapsed_s,
            duration_s,
            error,
        })
    }
}
fn budget(elapsed_s: f64, attempted_steps: usize, accepted_intervals: usize) -> Error {
    advance_error(EffusionAdvanceError::WorkBudget {
        elapsed_s,
        attempted_steps,
        accepted_intervals,
    })
}
fn retryable(error: &Error) -> bool {
    match error {
        Error::Fill {
            error: EffusionFillError::EntropyDecrease { .. },
            ..
        } => true,
        Error::Fill {
            error: EffusionFillError::Transport(error),
            ..
        }
        | Error::Advance(EffusionAdvanceError::Step { error, .. }) => {
            native::retryable_native(error)
        }
        // Initial receiver temperature/composition is set by the incident
        // beam. Reducing time cannot repair an unsupported caloric source.
        _ => false,
    }
}
fn count_attempt(attempted: &mut usize, controls: EffusionAdvanceControls) -> Result<(), Error> {
    if *attempted == controls.maximum_attempts {
        return Err(budget(0.0, *attempted, 0));
    }
    *attempted += 1;
    Ok(())
}

fn pair(
    candidate: &EffusionFillCandidate,
    link: EffusionFillLink<'_>,
) -> Result<BTreeMap<u64, ClosedGasParcel>, Error> {
    let receiver = candidate.receiver().gas().ok_or(Error::Invalid(
        "a positive first-fill probe lost its receiver inventory",
    ))?;
    Ok(BTreeMap::from([
        (link.donor, candidate.donor().clone()),
        (link.receiver, receiver.clone()),
    ]))
}

fn add(total: &mut Dd, value: Dd) -> Result<(), Error> {
    *total = *total + value;
    if !value.hi.is_finite()
        || !value.lo.is_finite()
        || !total.hi.is_finite()
        || !total.lo.is_finite()
    {
        return Err(Error::Invalid(
            "initial effusion error coordinate is unrepresentable",
        ));
    }
    Ok(())
}

// Relative physical storage is also the net transported molecular energy of
// each rigid chamber, because no external heat/work is present. Use exact
// original source increments, never rounded absolute enthalpy differences.
fn energy_coordinate(
    original: &ClosedGasParcel,
    gas: &ClosedGasParcel,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let inventory = gas
        .species_inventory_kg()
        .map_err(|_| Error::Invalid("unrepresentable filled species inventory"))?;
    let sources = original.model().thermodynamics().components();
    let mut result = Dd::ZERO;
    for (id, mass) in inventory {
        poll(cancelled)?;
        let index = sources
            .binary_search_by(|s| s.species().cmp(&id))
            .map_err(|_| Error::Invalid("filled endpoint lost its original species source"))?;
        let source = &sources[index];
        let change = source
            .molar_internal_energy_change(
                Temperature::new(original.state().state().temperature),
                Temperature::new(gas.state().state().temperature),
            )
            .map_err(|_| Error::Invalid("initial error estimate left original caloric support"))?
            .value();
        add(
            &mut result,
            Dd::from_f64(mass) * Dd::from_f64(change) / Dd::from_f64(source.molar_mass().value()),
        )?;
    }
    Ok(result.to_f64())
}

fn initial_error(
    original: &ClosedGasParcel,
    coarse: &EffusionFillCandidate,
    first: &EffusionFillCandidate,
    fine: &ReservoirEffusionCandidate,
    link: EffusionFillLink<'_>,
    controls: EffusionAdvanceControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let coarse_pair = pair(coarse, link)?;
    let mut error = native::state_error(
        &BTreeMap::from([(link.donor, original.clone())]),
        &coarse_pair,
        fine.parcels(),
        controls,
        cancelled,
    )
    .map_err(advance_error)?;
    for (&id, a) in &coarse_pair {
        let x = energy_coordinate(original, a, cancelled)?;
        let y = energy_coordinate(original, &fine.parcels()[&id], cancelled)?;
        error = error.max(
            native::ratio(
                x,
                y,
                x.abs().max(y.abs()),
                controls.stream_energy_tolerance_j,
                controls.relative_tolerance,
            )
            .map_err(advance_error)?,
        );
    }
    // Independent gross directions retain receiver backflow, including when
    // opposite streams nearly cancel in the final species inventory.
    let mut streams: BTreeMap<(u64, u64, SpeciesId), [Dd; 2]> = BTreeMap::new();
    for (probe, step) in [(0, coarse.step()), (1, first.step())] {
        for (species, &mass) in &step.captured_species_kg {
            poll(cancelled)?;
            add(
                &mut streams
                    .entry((link.donor, link.receiver, species.clone()))
                    .or_default()[probe],
                Dd::from_f64(mass),
            )?;
        }
    }
    for row in &fine.step().exchange.fluxes {
        poll(cancelled)?;
        let row = &row.transfer;
        add(
            &mut streams
                .entry((row.donor, row.receiver, row.species.clone()))
                .or_default()[1],
            Dd::from_f64(row.mass_kg),
        )?;
    }
    for values in streams.values() {
        poll(cancelled)?;
        let (a, b) = (values[0].to_f64(), values[1].to_f64());
        error = error.max(
            native::ratio(
                a,
                b,
                a.max(b),
                controls.species_mass_tolerance_kg,
                controls.relative_tolerance,
            )
            .map_err(advance_error)?,
        );
    }
    Ok(error)
}

fn first_rows(step: &EffusionFillStep) -> usize {
    let native = &step.transport;
    1 + step.captured_species_kg.len()
        + native.exchange.boundaries.len()
        + native.exchange.parcels.len()
        + native
            .outlets
            .iter()
            .map(|o| 1 + o.flux.species.len())
            .sum::<usize>()
}

// Only a numerical startup ceiling: p <= n_dot*R*T_max*dt/V. The original
// donor pressure is already source-admitted, while a coarse vacuum-only fill
// could exceed it by orders of magnitude in a small receiver. Start below that
// pressure without evaluating an unsupported gas, then let actual reciprocal
// flow and ordinary error control choose every subsequent interval. This is
// never a pressure floor or a restriction on later thermal transpiration.
fn initial_pressure_limit(
    donor: &ClosedGasParcel,
    receiver: &EvacuatedGasVolume,
    link: EffusionFillLink<'_>,
    receiver_solve: GasEnergySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<f64>, Error> {
    receiver_solve
        .validate()
        .map_err(|e| fill_error(EffusionFillError::Receiver(e), 0.0))?;
    let beam = donor
        .effusive_flux_to_vacuum(link.interface, link.area_m2, &mut *cancelled)
        .map_err(|e| native_error(face_error(link.id, e).into(), 0.0, 0.0))?;
    if beam.molar_rate_mol_s == 0.0 {
        return Ok(None);
    }
    let pressure = donor
        .state()
        .state()
        .pressure
        .min(link.interface.support().pressure_pa.1);
    let ceiling = positive_product_ratio(
        &[pressure, receiver.volume_m3()],
        &[
            beam.molar_rate_mol_s,
            R,
            receiver_solve.maximum_temperature_k,
        ],
    )
    .map_err(|e| native_error(face_error(link.id, e).into(), 0.0, 0.0))?;
    // Leave representational room for the original composition normalization
    // and EOS reconstruction. Actual source admission still decides acceptance.
    let ceiling = (ceiling * (1.0 - 512.0 * f64::EPSILON)).next_down();
    if !ceiling.is_finite() || ceiling <= 0.0 {
        return Err(Error::Invalid(
            "initial receiver pressure ceiling is unrepresentable",
        ));
    }
    Ok(Some(ceiling))
}

fn continuation_error(error: EffusionAdvanceError, elapsed: f64, attempted: usize) -> Error {
    advance_error(match error {
        EffusionAdvanceError::Step {
            elapsed_s,
            duration_s,
            error,
        } => EffusionAdvanceError::Step {
            elapsed_s: elapsed + elapsed_s,
            duration_s,
            error,
        },
        EffusionAdvanceError::MinimumStep {
            elapsed_s,
            duration_s,
            normalized_error,
        } => EffusionAdvanceError::MinimumStep {
            elapsed_s: elapsed + elapsed_s,
            duration_s,
            normalized_error,
        },
        EffusionAdvanceError::WorkBudget {
            elapsed_s,
            attempted_steps,
            accepted_intervals,
        } => EffusionAdvanceError::WorkBudget {
            elapsed_s: elapsed + elapsed_s,
            attempted_steps: attempted + attempted_steps,
            accepted_intervals: 1 + accepted_intervals,
        },
        other => other,
    })
}

/// Advance one original finite donor and one exactly evacuated finite receiver.
///
/// One accuracy/work policy governs initialization AND subsequent reciprocal
/// flow. Accepted first-fill history, actual source identities, positive trace
/// species and auxiliary query axes are retained. Final closed mass, energy and
/// entropy are checked against the original donor plus zero receiver inventory
/// using one fixed budget. Any refusal publishes neither endpoint.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate(
    donor: &ClosedGasParcel,
    receiver: &EvacuatedGasVolume,
    link: EffusionFillLink<'_>,
    duration_s: f64,
    donor_solve: GasEnergySolve,
    receiver_solve: GasEnergySolve,
    physical: GasTransportSolve,
    controls: EffusionAdvanceControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<EffusionFillAdvanceCandidate, Error> {
    poll(&mut cancelled)?;
    controls.validate().map_err(advance_error)?;
    if !duration_s.is_finite() || duration_s <= 0.0 || link.donor == link.receiver {
        return Err(Error::Invalid("positive finite filling duration required"));
    }
    let links = [EffusionNetworkLink {
        id: link.id,
        parcels: [link.donor, link.receiver],
        interface: link.interface,
        area_m2: link.area_m2,
    }];
    let paths = GasTransportPaths {
        effusion: &links,
        ..GasTransportPaths::default()
    };
    let endpoints = BTreeMap::from([(link.donor, donor_solve), (link.receiver, receiver_solve)]);
    let mut span = duration_s.min(controls.initial_step_s);
    if let Some(ceiling) =
        initial_pressure_limit(donor, receiver, link, receiver_solve, &mut cancelled)?
    {
        if ceiling < span {
            if ceiling < controls.minimum_step_s {
                return Err(Error::InitialPressureLimit {
                    maximum_step_s: ceiling,
                    minimum_step_s: controls.minimum_step_s,
                });
            }
            span = ceiling;
        }
    }
    let mut attempted = 0;
    let mut rejected = 0;

    let (first, second, maximum_error) = loop {
        poll(&mut cancelled)?;
        let half = 0.5 * span;
        if half <= 0.0 || half + half != span {
            return Err(Error::Invalid(
                "first-fill half interval is unrepresentable",
            ));
        }
        let trial = (|| {
            count_attempt(&mut attempted, controls)?;
            let coarse = step_candidate(
                donor,
                receiver,
                link,
                span,
                donor_solve,
                receiver_solve,
                physical,
                &mut cancelled,
            )
            .map_err(|e| fill_error(e, span))?;
            if coarse.receiver().gas().is_none() {
                return Ok((coarse, None, 0.0));
            }
            count_attempt(&mut attempted, controls)?;
            let first = step_candidate(
                donor,
                receiver,
                link,
                half,
                donor_solve,
                receiver_solve,
                physical,
                &mut cancelled,
            )
            .map_err(|e| fill_error(e, half))?;
            let start = pair(&first, link)?;
            count_attempt(&mut attempted, controls)?;
            let second = reservoir::step_candidate(
                &start,
                paths,
                &[],
                &[],
                half,
                &endpoints,
                physical,
                &mut cancelled,
            )
            .map_err(|e| native_error(e, half, half))?;
            let error = initial_error(
                donor,
                &coarse,
                &first,
                &second,
                link,
                controls,
                &mut cancelled,
            )?;
            Ok::<_, Error>((first, Some(second), error))
        })();
        match trial {
            Ok(result) if result.2 <= 1.0 => break result,
            Ok((_, _, error)) => {
                rejected += 1;
                let smaller = (0.5 * span).max(controls.minimum_step_s);
                if smaller >= span {
                    return Err(advance_error(EffusionAdvanceError::MinimumStep {
                        elapsed_s: 0.0,
                        duration_s: span,
                        normalized_error: error,
                    }));
                }
                span = smaller;
            }
            Err(error) if retryable(&error) => {
                rejected += 1;
                let smaller = (0.5 * span).max(controls.minimum_step_s);
                if smaller >= span {
                    return Err(error);
                }
                span = smaller;
            }
            Err(error) => return Err(error),
        }
    };

    if second.is_none() {
        // An exactly zero beam at fixed geometry/source state remains inactive.
        // Retain an actual full-duration receipt, including source admission.
        let accepted = if span == duration_s {
            first
        } else {
            count_attempt(&mut attempted, controls)?;
            step_candidate(
                donor,
                receiver,
                link,
                duration_s,
                donor_solve,
                receiver_solve,
                physical,
                &mut cancelled,
            )
            .map_err(|e| fill_error(e, duration_s))?
        };
        if accepted.receiver().gas().is_some() {
            return Err(Error::Invalid(
                "an inactive molecular beam changed under duration alone",
            ));
        }
        let (donor, receiver, first_fill) = accepted.into_parts();
        let balance = first_fill.balance;
        poll(&mut cancelled)?;
        return Ok(EffusionFillAdvanceCandidate {
            donor,
            receiver,
            step: EffusionFillAdvanceStep {
                duration_s,
                attempted_steps: attempted,
                accepted_intervals: 1,
                rejected_intervals: rejected,
                maximum_accepted_error: 0.0,
                first_fill,
                first_reciprocal: None,
                continuation: None,
                balance,
            },
        });
    }

    let (_, _, first_fill) = first.into_parts();
    let (mut current, first_reciprocal) = second.unwrap().into_parts();
    let mut retained = first_rows(&first_fill) + native::retained_rows(&first_reciprocal);
    let mut accepted_intervals = 1;
    let mut maximum_accepted_error = maximum_error;
    let remaining = duration_s - span;
    let continuation = if remaining > 0.0 {
        let available = controls.maximum_attempts - attempted;
        let intervals = controls.maximum_accepted_intervals - 1;
        if available < 3 || intervals == 0 {
            return Err(budget(span, attempted, 1));
        }
        let next_controls = EffusionAdvanceControls {
            initial_step_s: span.clamp(controls.minimum_step_s, controls.maximum_step_s),
            maximum_attempts: available,
            maximum_accepted_intervals: intervals,
            ..controls
        };
        let rows_available = native::MAX_EFFUSION_RETAINED_ROWS
            .checked_sub(retained)
            .filter(|&n| n > 0)
            .ok_or_else(|| budget(span, attempted, 1))?;
        let next = native::advance_candidate_with_row_limit(
            &current,
            paths,
            &[],
            &[],
            &[],
            remaining,
            &endpoints,
            physical,
            next_controls,
            rows_available,
            &mut cancelled,
        )
        .map_err(|e| continuation_error(e, span, attempted))?;
        let (next, report) = next.into_parts();
        current = next;
        attempted += report.attempted_steps;
        accepted_intervals += report.accepted_intervals;
        rejected += report.rejected_intervals;
        maximum_accepted_error = maximum_accepted_error.max(report.maximum_accepted_error);
        for step in &report.accepted_steps {
            poll(&mut cancelled)?;
            retained = retained
                .checked_add(native::retained_rows(step))
                .ok_or_else(|| budget(duration_s, attempted, accepted_intervals))?;
        }
        Some(report)
    } else {
        None
    };
    if retained > native::MAX_EFFUSION_RETAINED_ROWS {
        return Err(budget(duration_s, attempted, accepted_intervals));
    }
    let balance = audit_pair(
        donor,
        (link.donor, &current[&link.donor]),
        Some((link.receiver, &current[&link.receiver])),
        physical,
        &mut cancelled,
    )
    .map_err(|e| fill_error(e, duration_s))?;
    let donor = current
        .remove(&link.donor)
        .ok_or(Error::Invalid("filled donor identity missing"))?;
    let receiver = EffusionReceiver::Filled(
        current
            .remove(&link.receiver)
            .ok_or(Error::Invalid("filled receiver identity missing"))?,
    );
    poll(&mut cancelled)?;
    Ok(EffusionFillAdvanceCandidate {
        donor,
        receiver,
        step: EffusionFillAdvanceStep {
            duration_s,
            attempted_steps: attempted,
            accepted_intervals,
            rejected_intervals: rejected,
            maximum_accepted_error,
            first_fill,
            first_reciprocal: Some(first_reciprocal),
            continuation,
            balance,
        },
    })
}

#[cfg(test)]
mod tests;
