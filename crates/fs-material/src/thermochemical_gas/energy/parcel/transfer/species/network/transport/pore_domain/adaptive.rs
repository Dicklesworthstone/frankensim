//! Local temporal accuracy for the original finite pore diffusion owner.
//!
//! Ordinary and adsorbing pores share one step-doubling controller. Every probe
//! is one original physical interval; only the two fine intervals contribute
//! accepted state and transfers. Fixed geometry, source domains and physical
//! balances remain owned by the native candidate. No extrapolation is applied.

use super::*;

/// Absolute component scales and a relative allowance on actual state changes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PoreTimeControls {
    pub gas_species_tolerance_kg: f64,
    pub gas_temperature_tolerance_k: f64,
    /// Multiplies the largest coarse/fine change from the accepted starting state.
    /// Large carrier inventories cannot hide a small transported species change.
    pub relative_tolerance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PoreAdvanceReport {
    pub interval_s: [f64; 2],
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    /// Includes failed and rejected probes as well as accepted fine intervals.
    pub evaluations: usize,
    /// Gas species inventory and gas temperature, respectively.
    pub maximum_accepted_error_ratios: [f64; 2],
    pub suggested_step_s: f64,
    /// Whole-window balances and gross transfers from accepted fine steps only.
    pub physical: PoreDomainStep,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PoreAdvanceCandidate {
    domain: PoreDomain,
    report: PoreAdvanceReport,
}
impl PoreAdvanceCandidate {
    pub const fn domain(&self) -> &PoreDomain {
        &self.domain
    }
    pub const fn report(&self) -> &PoreAdvanceReport {
        &self.report
    }
    pub fn into_parts(self) -> (PoreDomain, PoreAdvanceReport) {
        (self.domain, self.report)
    }
}

/// Original physical refusal or bounded temporal-control failure. The ordinary
/// gas owner uses two error components; the shared adsorbing driver uses four.
#[derive(Clone, Debug, PartialEq)]
pub enum PoreAdvanceError<const N: usize = 2> {
    Invalid(&'static str),
    Physical {
        time_s: f64,
        step_s: f64,
        error: PoreDomainError,
    },
    Accuracy {
        time_s: f64,
        step_s: f64,
        error_ratios: [f64; N],
    },
    TimeResolution {
        time_s: f64,
        target_s: f64,
    },
    WorkLimit {
        time_s: f64,
        accepted_substeps: usize,
        evaluations: usize,
    },
    Cancelled,
}
impl<const N: usize> core::fmt::Display for PoreAdvanceError<N> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive pore domain: {self:?}")
    }
}
impl<const N: usize> std::error::Error for PoreAdvanceError<N> {}
impl<const N: usize> PoreAdvanceError<N> {
    fn physical(error: PoreDomainError, time_s: f64, step_s: f64) -> Self {
        if error == PoreDomainError::Cancelled {
            Self::Cancelled
        } else {
            Self::Physical {
                time_s,
                step_s,
                error,
            }
        }
    }
    fn retryable(&self) -> bool {
        match self {
            Self::Physical {
                error: PoreDomainError::Transport { error, .. },
                ..
            } => super::retryable(error),
            Self::Physical {
                error: PoreDomainError::Adsorption { error, .. },
                ..
            } => adsorption::retryable(error),
            Self::Accuracy { .. } => true,
            _ => false,
        }
    }
}

fn poll<const N: usize>(cancelled: &mut impl FnMut() -> bool) -> Result<(), PoreAdvanceError<N>> {
    if cancelled() {
        Err(PoreAdvanceError::Cancelled)
    } else {
        Ok(())
    }
}

pub(super) struct TimeWindow<D, S, const N: usize> {
    pub state: D,
    pub last: S,
    pub accepted: usize,
    pub rejected: usize,
    pub evaluations: usize,
    pub minimum: f64,
    pub maximum: f64,
    pub maxima: [f64; N],
    pub suggested: f64,
}

#[allow(clippy::too_many_arguments)]
fn probe<D, S, C: FnMut() -> bool, const N: usize>(
    state: &D,
    time_s: f64,
    dt: f64,
    accepted: usize,
    evaluations: &mut usize,
    maximum: usize,
    advance: &mut impl FnMut(&D, f64, &mut C) -> Result<(D, S), PoreDomainError>,
    cancelled: &mut C,
) -> Result<(D, S), PoreAdvanceError<N>> {
    poll(cancelled)?;
    if *evaluations == maximum {
        return Err(PoreAdvanceError::WorkLimit {
            time_s,
            accepted_substeps: 2 * accepted,
            evaluations: *evaluations,
        });
    }
    *evaluations += 1;
    advance(state, dt, cancelled).map_err(|e| PoreAdvanceError::physical(e, time_s, dt))
}

/// The two pore owners differ only in their native state, error components and
/// accepted account reducer. Their physical laws are never implemented here.
#[allow(clippy::too_many_arguments)]
pub(super) fn run<D: Clone, S, C: FnMut() -> bool, const N: usize>(
    original: &D,
    target_s: f64,
    controls: PoreDomainSolve,
    time: impl Fn(&D) -> f64,
    mut advance: impl FnMut(&D, f64, &mut C) -> Result<(D, S), PoreDomainError>,
    mut errors: impl FnMut(&D, &D, &D, &mut C) -> Result<[f64; N], PoreAdvanceError<N>>,
    mut record: impl FnMut(&S, &mut C) -> Result<(), PoreAdvanceError<N>>,
    cancelled: &mut C,
) -> Result<TimeWindow<D, S, N>, PoreAdvanceError<N>> {
    poll(cancelled)?;
    let start = time(original);
    if !target_s.is_finite()
        || target_s <= start
        || !positive(target_s - start)
        || !positive(controls.minimum_step_s)
        || !positive(controls.maximum_step_s)
        || controls.minimum_step_s > controls.maximum_step_s
        || !(2..=MAX_PORE_SUBSTEPS).contains(&controls.maximum_substeps)
        || !(3..=MAX_PORE_ATTEMPTS).contains(&controls.maximum_attempts)
    {
        return Err(PoreAdvanceError::Invalid(
            "forward target and bounded positive comparison/evaluation controls required",
        ));
    }
    let mut state = original.clone();
    let mut ceiling = controls.maximum_step_s;
    let (mut accepted, mut rejected, mut evaluations) = (0, 0, 0);
    let (mut minimum, mut maximum) = (f64::INFINITY, 0.0_f64);
    let mut maxima = [0.0_f64; N];
    let mut last = None;
    while time(&state) < target_s {
        poll(cancelled)?;
        let now = time(&state);
        if 2 * accepted + 2 > controls.maximum_substeps {
            return Err(PoreAdvanceError::WorkLimit {
                time_s: now,
                accepted_substeps: 2 * accepted,
                evaluations,
            });
        }
        let end = if ceiling >= target_s - now {
            target_s
        } else {
            now + ceiling
        };
        let span = end - now;
        let middle = now + 0.5 * span;
        if !(now < middle && middle < end) {
            return Err(PoreAdvanceError::TimeResolution {
                time_s: now,
                target_s: end,
            });
        }
        let (h1, h2) = (middle - now, end - middle);
        let attempted = (|| {
            let coarse = probe(
                &state,
                now,
                span,
                accepted,
                &mut evaluations,
                controls.maximum_attempts,
                &mut advance,
                cancelled,
            )?;
            let first = probe(
                &state,
                now,
                h1,
                accepted,
                &mut evaluations,
                controls.maximum_attempts,
                &mut advance,
                cancelled,
            )?;
            if time(&first.0).to_bits() != middle.to_bits() {
                return Err(PoreAdvanceError::TimeResolution {
                    time_s: now,
                    target_s: middle,
                });
            }
            let second = probe(
                &first.0,
                middle,
                h2,
                accepted,
                &mut evaluations,
                controls.maximum_attempts,
                &mut advance,
                cancelled,
            )?;
            if time(&coarse.0).to_bits() != end.to_bits()
                || time(&second.0).to_bits() != end.to_bits()
            {
                return Err(PoreAdvanceError::TimeResolution {
                    time_s: now,
                    target_s: end,
                });
            }
            let mut ratios = errors(&state, &coarse.0, &second.0, cancelled)?;
            // For first-order local error C*h^2, the fine error is C*(h1^2+h2^2)
            // and coarse-minus-fine is 2*C*h1*h2. Use actual represented halves.
            let (u, v) = (h1 / span, h2 / span);
            let factor = (u * u + v * v) / (2.0 * u * v);
            for ratio in &mut ratios {
                *ratio *= factor;
                if !ratio.is_finite() || *ratio < 0.0 {
                    return Err(PoreAdvanceError::Invalid(
                        "unrepresentable pore temporal error ratio",
                    ));
                }
            }
            if ratios.iter().any(|r| *r > 1.0) {
                return Err(PoreAdvanceError::Accuracy {
                    time_s: now,
                    step_s: span,
                    error_ratios: ratios,
                });
            }
            Ok((first, second, ratios))
        })();
        match attempted {
            Ok((first, (next, step), ratios)) => {
                record(&first.1, cancelled)?;
                record(&step, cancelled)?;
                state = next;
                last = Some(step);
                accepted += 1;
                minimum = minimum.min(h1.min(h2));
                maximum = maximum.max(h1.max(h2));
                for (max, r) in maxima.iter_mut().zip(ratios) {
                    *max = max.max(r);
                }
                ceiling = if ratios.iter().all(|r| *r <= 0.25) {
                    span + span.min(controls.maximum_step_s - span)
                } else {
                    span
                };
                ceiling = ceiling.clamp(controls.minimum_step_s, controls.maximum_step_s);
            }
            Err(error) => {
                let smaller = (0.5 * span).max(controls.minimum_step_s);
                if error.retryable() && smaller < span {
                    rejected += 1;
                    ceiling = smaller;
                } else {
                    return Err(error);
                }
            }
        }
    }
    poll(cancelled)?;
    Ok(TimeWindow {
        state,
        last: last.expect("positive duration accepted fine steps"),
        accepted,
        rejected,
        evaluations,
        minimum,
        maximum,
        maxima,
        suggested: ceiling,
    })
}

fn errors(
    old: &PoreDomain,
    coarse: &PoreDomain,
    fine: &PoreDomain,
    controls: PoreTimeControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 2], PoreAdvanceError> {
    let mut result = [0.0_f64; 2];
    let ratio = |before: f64, a: f64, b: f64, absolute: f64| {
        let allowance =
            absolute + controls.relative_tolerance * (a - before).abs().max((b - before).abs());
        let value = (a - b).abs() / allowance;
        if !positive(allowance) || !value.is_finite() {
            Err(PoreAdvanceError::Invalid(
                "unrepresentable pore component allowance",
            ))
        } else {
            Ok(value)
        }
    };
    for (&id, before) in &old.parcels {
        poll(cancelled)?;
        let (a, b) = (&coarse.parcels[&id], &fine.parcels[&id]);
        let get = |gas: &ClosedGasParcel| {
            inventory(gas.model(), gas.mass_kg())
                .map_err(|e| PoreAdvanceError::physical(gas_error(e), old.time_s, 0.0))
        };
        let (m0, ma, mb) = (get(before)?, get(a)?, get(b)?);
        // Include newly arrived species as well as original traces. Repeated
        // keys cost no allocation and cannot change the componentwise maximum.
        for species in m0.keys().chain(ma.keys()).chain(mb.keys()) {
            result[0] = result[0].max(ratio(
                m0.get(species).copied().unwrap_or(0.0),
                ma.get(species).copied().unwrap_or(0.0),
                mb.get(species).copied().unwrap_or(0.0),
                controls.gas_species_tolerance_kg,
            )?);
        }
        result[1] = result[1].max(ratio(
            before.state().state().temperature,
            a.state().state().temperature,
            b.state().state().temperature,
            controls.gas_temperature_tolerance_k,
        )?);
    }
    Ok(result)
}

impl PoreDomain {
    /// Propose a complete window with local species/temperature accuracy.
    /// The comparison span is bounded by maximum_step_s; its minimum is a retry
    /// floor, so a shorter final remainder still undergoes an accuracy check.
    /// Maximum substeps counts accepted fine steps and maximum attempts counts
    /// every native probe. Source/domain failures are fatal, and no partial
    /// trajectory is published. This is first-order local error control, not a
    /// global error bound or spatial-refinement certificate.
    pub fn step_candidate_adaptive_to<C: FnMut() -> bool>(
        &self,
        target_s: f64,
        faces: &[PoreFace<'_>],
        endpoints: &[GasEnergySolve],
        controls: PoreDomainSolve,
        accuracy: PoreTimeControls,
        mut cancelled: C,
    ) -> Result<PoreAdvanceCandidate, PoreAdvanceError> {
        poll(&mut cancelled)?;
        if !positive(accuracy.gas_species_tolerance_kg)
            || !positive(accuracy.gas_temperature_tolerance_k)
            || !accuracy.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&accuracy.relative_tolerance)
        {
            return Err(PoreAdvanceError::Invalid(
                "positive pore temporal scales and relative tolerance in [0,1) required",
            ));
        }
        let mut ledger: BTreeMap<(usize, SpeciesId), [Sum; 3]> = BTreeMap::new();
        let run = run(
            self,
            target_s,
            controls,
            |state| state.time_s,
            |state, dt, cancelled| {
                // Prevent an inner retry from hiding a failed coarse/fine probe.
                state
                    .step_candidate(
                        dt,
                        faces,
                        endpoints,
                        PoreDomainSolve {
                            maximum_step_s: dt,
                            minimum_step_s: dt,
                            maximum_substeps: 1,
                            maximum_attempts: 1,
                            ..controls
                        },
                        cancelled,
                    )
                    .map(PoreDomainCandidate::into_parts)
            },
            |old, coarse, fine, cancelled| errors(old, coarse, fine, accuracy, cancelled),
            |step, cancelled| {
                for row in &step.faces {
                    poll(cancelled)?;
                    let sums = ledger.entry((row.face, row.species.clone())).or_default();
                    for (sum, value) in sums.iter_mut().zip([
                        row.mass_left_to_right_kg,
                        row.mass_right_to_left_kg,
                        row.enthalpy_left_to_right_j,
                    ]) {
                        sum.add(value).map_err(|e| {
                            PoreAdvanceError::physical(gas_error(e), self.time_s, 0.0)
                        })?;
                    }
                }
                Ok(())
            },
            &mut cancelled,
        )?;
        let (mass, species, energy, entropy, pore_mass) = window_balances(
            &self.parcels,
            &run.state.parcels,
            self.cells.len(),
            controls.transport,
            &mut cancelled,
        )
        .map_err(|e| PoreAdvanceError::physical(e, target_s, 0.0))?;
        let mut step = run.last;
        step.duration_s = target_s - self.time_s;
        step.accepted_substeps = 2 * run.accepted;
        step.attempts = run.evaluations;
        step.minimum_accepted_step_s = run.minimum;
        step.maximum_accepted_step_s = run.maximum;
        step.total_mass_residual_kg = mass;
        step.maximum_species_mass_residual_kg = species;
        step.total_energy_residual_j = energy;
        step.entropy_change_j_k = entropy;
        step.pore_mass_change_kg = pore_mass;
        step.faces = ledger
            .into_iter()
            .map(|((face, species), sums)| PoreFaceTransfer {
                face,
                species,
                mass_left_to_right_kg: sums[0].value,
                mass_right_to_left_kg: sums[1].value,
                enthalpy_left_to_right_j: sums[2].value,
            })
            .collect();
        poll(&mut cancelled)?;
        Ok(PoreAdvanceCandidate {
            domain: run.state,
            report: PoreAdvanceReport {
                interval_s: [self.time_s, target_s],
                accepted_intervals: run.accepted,
                rejected_intervals: run.rejected,
                evaluations: run.evaluations,
                maximum_accepted_error_ratios: run.maxima,
                suggested_step_s: run.suggested,
                physical: step,
            },
        })
    }

    /// Publish all pore and reservoir inventories only after the complete
    /// accuracy-controlled trajectory, whole-window audit and cancellation gate.
    pub fn advance_adaptive_to(
        &mut self,
        target_s: f64,
        faces: &[PoreFace<'_>],
        endpoints: &[GasEnergySolve],
        controls: PoreDomainSolve,
        accuracy: PoreTimeControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PoreAdvanceReport, PoreAdvanceError> {
        let candidate = self.step_candidate_adaptive_to(
            target_s,
            faces,
            endpoints,
            controls,
            accuracy,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        let (next, report) = candidate.into_parts();
        *self = next;
        Ok(report)
    }
}
