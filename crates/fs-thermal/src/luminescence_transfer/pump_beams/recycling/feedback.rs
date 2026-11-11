//! Safeguarded source closure on original finite-population optical stages.

use super::*;

#[derive(Clone, Copy)]
struct Source {
    cell: usize,
    state: usize,
    upper_rate: f64,
}

struct Trial<const B: usize> {
    history: History<B>,
    residual: Vec<f64>,
    norm: f64,
    consistent: bool,
}

impl<const B: usize> SaturablePhotonRecycling<B> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn frozen_history(
        &self,
        frozen: &mut SaturablePumpBeams<B>,
        initial: &[[f64; 3]],
        incident: &[[f64; B]],
        local_pumps: &[LocalPhotonPump],
        duration: f64,
        limits: SaturableRecyclingLimits,
        evaluations: &mut usize,
        coupling: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<History<B>>, SingletTripletError> {
        let nonlinear = duration > 0.0
            && frozen.paths.iter().zip(incident).any(|(path, rates)| {
                path.segments.iter().any(|segment| {
                    (0..B).any(|lane| {
                        rates[lane] > 0.0
                            && segment.optical_depth_per_ground[lane]
                                .iter()
                                .any(|&x| x > 0.0)
                    })
                })
            });
        let mut attribution = SourceAttribution::new(initial.len(), self.external_paths);
        let Some(beams) = frozen.uniform_attributed(
            initial,
            incident,
            local_pumps,
            duration,
            1,
            nonlinear,
            limits.pump,
            evaluations,
            coupling,
            Some(&mut attribution),
            cancelled,
        )?
        else {
            return Ok(None);
        };
        Ok(Some(History {
            beams,
            external_absorbed: attribution
                .external_absorbed
                .into_iter()
                .map(values)
                .collect::<Result<_, _>>()?,
            recycled_absorbed: attribution
                .recycled_absorbed
                .into_iter()
                .map(values)
                .collect::<Result<_, _>>()?,
            external_events: attribution
                .external_events
                .into_iter()
                .map(values)
                .collect::<Result<_, _>>()?,
            recycled_events: attribution
                .recycled_events
                .into_iter()
                .map(values)
                .collect::<Result<_, _>>()?,
            external_relaxation: attribution
                .external_relaxation
                .into_iter()
                .map(values)
                .collect::<Result<_, _>>()?,
            recycled_relaxation: attribution
                .recycled_relaxation
                .into_iter()
                .map(values)
                .collect::<Result<_, _>>()?,
            external_passive: attribution
                .external_passive
                .into_iter()
                .map(|v| v.value())
                .collect::<Result<_, _>>()?,
            recycled_passive: attribution
                .recycled_passive
                .into_iter()
                .map(|v| v.value())
                .collect::<Result<_, _>>()?,
            source_allowance_j: 0.0,
        }))
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn feedback(
        &self,
        frozen: &mut SaturablePumpBeams<B>,
        initial: &[[f64; 3]],
        external: &[[f64; B]],
        local_pumps: &[LocalPhotonPump],
        duration: f64,
        limits: SaturableRecyclingLimits,
        evaluations: &mut usize,
        coupling: &mut usize,
        feedback: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<History<B>>, SingletTripletError> {
        let source_capacity =
            initial
                .len()
                .checked_mul(2)
                .ok_or(SingletTripletError::Unrepresentable(
                    "finite photon source count",
                ))?;
        admit_storage::<Source>(source_capacity, "finite photon source storage")?;
        admit_storage::<f64>(source_capacity, "finite photon source vectors")?;
        admit_storage::<[Option<usize>; 2]>(initial.len(), "finite photon source indices")?;
        let mut sources = Vec::new();
        let mut source_index = vec![[None; 2]; initial.len()];
        let mut x = Vec::new();
        for (cell, model) in frozen.models.iter().enumerate() {
            poll(cancelled)?;
            let count = inventory(initial[cell])?;
            for state in 0..2 {
                let upper_rate = positive_product(
                    count,
                    model.radiative_rates_s()[state],
                    "finite source rate bound",
                )?;
                if upper_rate > 0.0 {
                    source_index[cell][state] = Some(sources.len());
                    sources.push(Source {
                        cell,
                        state,
                        upper_rate,
                    });
                    x.push(initial[cell][state + 1] / count);
                }
            }
        }
        let n = sources.len();
        let mut inverse = secant::SecantInverse::try_new(
            n,
            limits
                .maximum_secant_history
                .min(limits.maximum_feedback_iterations),
        )?;
        let mut stage_evaluations = 0usize;
        // Inner optical/root error must leave room for emitted-source closure
        // and the final external-energy audit. This tightens the nested solve;
        // it never changes the caller's physical or temporal allowances.
        let mut inner_limits = limits;
        inner_limits.pump.column.kinetic.relative_tolerance *= 0.0625;
        let mut evaluate = |position: &[f64],
                            cancelled: &mut _|
         -> Result<Option<Trial<B>>, SingletTripletError> {
            if stage_evaluations == limits.maximum_feedback_iterations {
                return Ok(None);
            }
            stage_evaluations += 1;
            *feedback = feedback
                .checked_add(1)
                .ok_or(SingletTripletError::Unrepresentable(
                    "photon feedback work counter",
                ))?;
            poll(cancelled)?;
            let mut incident = external.to_vec();
            for emitted in &self.emission {
                let rate =
                    if let Some(index) = source_index[emitted.source_cell][emitted.source_state] {
                        positive_product(
                            position[index],
                            sources[index].upper_rate,
                            "trial radiative rate",
                        )?
                    } else {
                        0.0
                    };
                let mut lanes = [0.0; B];
                for lane in 0..B {
                    lanes[lane] = positive_product(
                        rate,
                        emitted.photon_fractions[lane],
                        "trial spectral source",
                    )?;
                }
                incident.push(lanes);
            }
            let Some(mut history) = self.frozen_history(
                frozen,
                initial,
                &incident,
                local_pumps,
                duration,
                inner_limits,
                evaluations,
                coupling,
                cancelled,
            )?
            else {
                return Ok(None);
            };
            let mut residual = Vec::with_capacity(n);
            let mut norm = 0.0_f64;
            let mut closure_energy = Sum::default();
            // Magnitudes of the terms that cancel in this stage's closure c_k
            // and in its composed-source mismatch d_k; they bound the rounding
            // between the two (admit_composed_source() states the chain).
            let mut cancelling = Sum::default();
            for (index, source) in sources.iter().enumerate() {
                let scale = positive_product(source.upper_rate, duration, "radiative count bound")?;
                let actual = history.beams.cells[source.cell].radiative_events[source.state];
                let supplied =
                    positive_product(position[index], scale, "supplied radiative count")?;
                let difference = finite(supplied - actual, "radiative source closure")?;
                closure_energy.add(positive_product(
                    difference.abs(),
                    self.models()[source.cell].excitation_energies_j()[source.state],
                    "radiative closure energy",
                )?)?;
                for count in [supplied, actual] {
                    cancelling.add(positive_product(
                        count,
                        self.models()[source.cell].excitation_energies_j()[source.state],
                        "closure rounding term",
                    )?)?;
                }
                let f = difference / scale;
                if !f.is_finite() {
                    return Err(SingletTripletError::Unrepresentable(
                        "normalized radiative closure",
                    ));
                }
                norm = norm.max(f.abs());
                residual.push(f);
            }
            for (path, beam) in self
                .emission
                .iter()
                .zip(&history.beams.beams[self.external_paths..])
            {
                let events =
                    history.beams.cells[path.source_cell].radiative_events[path.source_state];
                for lane in 0..B {
                    let count = positive_product(
                        events,
                        path.photon_fractions[lane],
                        "emitted closure photons",
                    )?;
                    for photons in [count, beam.incident_photons[lane]] {
                        cancelling.add(positive_product(
                            photons,
                            path.path.photon_energies_j[lane],
                            "closure rounding term",
                        )?)?;
                    }
                }
            }
            // Gross radiative traffic NEVER enlarges this source tolerance.
            // Only actual initial storage, external supply and bath exchange
            // can pay for an unresolved internal photon mismatch.
            let scale = self.external_scale(initial, &history.beams)?;
            let allowance = 0.0625 * limits.pump.column.kinetic.relative_tolerance * scale;
            let consistent = closure_energy.value()? <= allowance;
            // A returned trial carries a_k + r_k + f c_k: the allowance it
            // closed to, r_k = 4 eps x the cancelling magnitudes above, and
            // fraction_roundoff x its closure (admit_composed_source() derives
            // both). Neither extra term enters the admission just above. A
            // midpoint history sums all three over its corrector stages.
            history.source_allowance_j = allowance
                + 4.0 * f64::EPSILON * cancelling.value()?
                + self.fraction_roundoff * closure_energy.value()?;
            Ok(Some(Trial {
                history,
                residual,
                norm,
                consistent,
            }))
        };
        let Some(mut trial) = evaluate(&x, cancelled)? else {
            return Ok(None);
        };
        loop {
            poll(cancelled)?;
            if trial.consistent {
                return Ok(Some(trial.history));
            }
            // Inverse Broyden learns the slow return modes with bounded
            // O(n * history) storage and application. At identity this remains
            // ordinary physical source iteration; no extra native evaluations
            // are used to assemble a dense Jacobian.
            let mut direction = vec![0.0; n];
            inverse.apply(&trial.residual, &mut direction, cancelled)?;
            for value in &mut direction {
                *value = -*value;
            }
            let mut alpha = feasible_step(&x, &direction);
            let mut accepted = None;
            for _ in 0..8 {
                if alpha <= 0.0 {
                    break;
                }
                let candidate: Vec<_> = x
                    .iter()
                    .zip(&direction)
                    .map(|(&a, &d)| a + alpha * d)
                    .collect();
                if candidate
                    .iter()
                    .all(|&a| a.is_finite() && (0.0..=1.0).contains(&a))
                {
                    let Some(next) = evaluate(&candidate, cancelled)? else {
                        return Ok(None);
                    };
                    if next.consistent || next.norm < trial.norm * (1.0 - 1e-4 * alpha) {
                        accepted = Some((candidate, next));
                        break;
                    }
                }
                alpha *= 0.5;
            }
            let (next_x, next) = if let Some(next) = accepted {
                next
            } else {
                // Reset an unhelpful secant approximation. The fallback is
                // another admissible native source evaluation, never clipping
                // a physical population or an accepted photon account.
                inverse.reset();
                let candidate: Vec<_> = x
                    .iter()
                    .zip(&trial.residual)
                    .map(|(&a, &f)| a - f)
                    .collect();
                if !candidate
                    .iter()
                    .all(|&a| a.is_finite() && (0.0..=1.0).contains(&a))
                {
                    return Ok(None);
                }
                let Some(next) = evaluate(&candidate, cancelled)? else {
                    return Ok(None);
                };
                (candidate, next)
            };
            let displacement: Vec<_> = next_x.iter().zip(&x).map(|(&a, &b)| a - b).collect();
            let change: Vec<_> = next
                .residual
                .iter()
                .zip(&trial.residual)
                .map(|(&a, &b)| a - b)
                .collect();
            inverse.update(&displacement, &change, cancelled)?;
            x = next_x;
            trial = next;
        }
    }
}

fn feasible_step(x: &[f64], direction: &[f64]) -> f64 {
    let mut alpha = 1.0_f64;
    for (&value, &delta) in x.iter().zip(direction) {
        if !delta.is_finite() {
            return 0.0;
        }
        let bound = if delta > 0.0 {
            (1.0 - value) / delta
        } else if delta < 0.0 {
            -value / delta
        } else {
            f64::INFINITY
        };
        if bound < alpha {
            alpha = 0.95 * bound;
        }
    }
    alpha
}
