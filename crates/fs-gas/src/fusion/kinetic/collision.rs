//! Positive energy diffusion with explicit heat rewards for both Maxwellian baths.

use super::grid::log_ratio;
use super::{AlphaEnergyGrid, MAX_ALPHA_BINS};
use crate::fusion::solve::Sum;
use crate::fusion::{AVOGADRO, DtFusionError, poll, positive};
use fs_math::det::{exp, expm1, ln};
use fs_qty::{Energy, Temperature, Time};
use fs_thermochem::coulomb::CoulombErrorV1;
use fs_thermochem::coulomb::test_particle::PreparedTestParticleCoulombV1;
use fs_thermochem::fusion::BOLTZMANN_J_PER_K;

const MAX_EDGES: usize = MAX_ALPHA_BINS - 1;
const MAX_POISSON_TERMS: usize = 96;
const MAX_EXPOSURE: f64 = 8.0;
const POISSON_TAIL_TOLERANCE: f64 = 1e-18;
const MAX_APPLICATIONS: usize = 1_000_000;

/// Frozen birth/death generator for the actual source-bound Maxwellian baths.
///
/// Each bath separately satisfies discrete detailed balance against
/// `width(E) sqrt(E) exp(-E/kT_b)`. The Scharfetter-Gummel flux uses the full
/// energy-diffusion coefficient, D=B/2, at the interior arithmetic face. It is
/// not a drag-only characteristic or a terminal-cutoff thermalization model.
pub(crate) struct AlphaCollisionOperator {
    grid: AlphaEnergyGrid,
    lambda: f64,
    upward: [f64; MAX_ALPHA_BINS],
    downward: [f64; MAX_ALPHA_BINS],
    staying: [f64; MAX_ALPHA_BINS],
    reward_j_per_mol: [[f64; MAX_ALPHA_BINS]; 2],
    maximum_reward_j_per_mol: [f64; 2],
}

/// One atomic accepted frozen collision interval. Amounts are in moles.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct AlphaCollisionStep {
    pub(crate) amounts: [f64; MAX_ALPHA_BINS],
    /// Signed energy deposited into [heavy, electron] baths. Energy diffusion
    /// can heat a cold kinetic population, making a bath's deposit negative.
    pub(crate) deposited_j: [f64; 2],
    pub(crate) applications: usize,
    /// Sum of analytic omitted-Poisson-probability bounds. Tail mass is retained
    /// in the last Markov state, never discarded or used to rescale populations.
    /// This bounds polynomial truncation in total variation; it excludes grid
    /// error, source/model error and floating-point roundoff.
    pub(crate) numerical_tail_bound: f64,
    /// Bounds on heat-reward error from the same Poisson-tail approximation.
    pub(crate) heat_tail_bound_j: [f64; 2],
    pub(crate) number_residual_mol: f64,
    pub(crate) energy_residual_j: f64,
}

impl AlphaCollisionOperator {
    pub(crate) fn new(
        grid: &AlphaEnergyGrid,
        prepared: &PreparedTestParticleCoulombV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DtFusionError> {
        Self::from_diffusion(
            grid,
            prepared.temperatures(),
            |energy| {
                prepared
                    .evaluate(energy)
                    .map(|value| value.energy_diffusion_j2_s)
                    .map_err(|error| match error {
                        CoulombErrorV1::Cancelled => DtFusionError::Cancelled,
                        error => DtFusionError::Coulomb(error),
                    })
            },
            cancelled,
        )
    }

    fn from_diffusion(
        grid: &AlphaEnergyGrid,
        temperatures: [Temperature; 2],
        mut diffusion: impl FnMut(Energy) -> Result<[f64; 2], DtFusionError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DtFusionError> {
        poll(&mut cancelled)?;
        let thermal_energy = temperatures.map(|t| BOLTZMANN_J_PER_K * t.value());
        if thermal_energy.iter().any(|&energy| !positive(energy)) {
            return Err(DtFusionError::Invalid(
                "positive finite alpha collision bath temperatures required",
            ));
        }
        let mut up = [[0.0; MAX_EDGES]; 2];
        let mut down = [[0.0; MAX_EDGES]; 2];
        for edge in 0..grid.len() - 1 {
            poll(&mut cancelled)?;
            let low = grid.energies()[edge].value();
            let high = grid.energies()[edge + 1].value();
            let gap = high - low;
            let coefficient = diffusion(Energy::new(low + 0.5 * gap))?;
            let log_energy_ratio = log_ratio(high, low);
            for bath in 0..2 {
                if !coefficient[bath].is_finite() || coefficient[bath] < 0.0 {
                    return Err(DtFusionError::Arithmetic(
                        "finite nonnegative alpha energy diffusion",
                    ));
                }
                let w = gap / thermal_energy[bath] - 0.5 * log_energy_ratio;
                if !w.is_finite() {
                    return Err(DtFusionError::Arithmetic("alpha collision edge exponent"));
                }
                up[bath][edge] = edge_rate(coefficient[bath], gap, grid.widths_j()[edge], w)?;
                down[bath][edge] =
                    edge_rate(coefficient[bath], gap, grid.widths_j()[edge + 1], -w)?;
            }
        }
        Self::from_rates(grid, &up, &down, cancelled)
    }

    fn from_rates(
        grid: &AlphaEnergyGrid,
        up: &[[f64; MAX_EDGES]; 2],
        down: &[[f64; MAX_EDGES]; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DtFusionError> {
        poll(&mut cancelled)?;
        let mut upward = [0.0; MAX_ALPHA_BINS];
        let mut downward = [0.0; MAX_ALPHA_BINS];
        let mut staying = [0.0; MAX_ALPHA_BINS];
        let mut lambda = 0.0_f64;
        for edge in 0..grid.len() - 1 {
            poll(&mut cancelled)?;
            for bath in 0..2 {
                if [up[bath][edge], down[bath][edge]]
                    .iter()
                    .any(|rate| !rate.is_finite() || *rate < 0.0)
                {
                    return Err(DtFusionError::Arithmetic(
                        "finite nonnegative alpha collision rates",
                    ));
                }
                upward[edge] += up[bath][edge];
                downward[edge + 1] += down[bath][edge];
            }
        }
        for i in 0..grid.len() {
            let outgoing = upward[i] + downward[i];
            if !outgoing.is_finite() {
                return Err(DtFusionError::Arithmetic("alpha collision outgoing rate"));
            }
            lambda = lambda.max(outgoing);
        }
        let mut reward_j_per_mol = [[0.0; MAX_ALPHA_BINS]; 2];
        let mut maximum_reward_j_per_mol = [0.0_f64; 2];
        if lambda > 0.0 {
            for i in 0..grid.len() {
                // Form the diagonal from the same summed outgoing rate. This
                // leaves a nonnegative probability even at the maximum column.
                staying[i] = 1.0 - (upward[i] + downward[i]) / lambda;
                upward[i] /= lambda;
                downward[i] /= lambda;
            }
            for edge in 0..grid.len() - 1 {
                poll(&mut cancelled)?;
                let gap = grid.energies()[edge + 1].value() - grid.energies()[edge].value();
                for bath in 0..2 {
                    // N is mol; the generator's energy coordinate is J/alpha.
                    reward_j_per_mol[bath][edge] -= (gap * (up[bath][edge] / lambda)) * AVOGADRO;
                    reward_j_per_mol[bath][edge + 1] +=
                        (gap * (down[bath][edge] / lambda)) * AVOGADRO;
                }
            }
            for bath in 0..2 {
                for value in &reward_j_per_mol[bath][..grid.len()] {
                    if !value.is_finite() {
                        return Err(DtFusionError::Arithmetic("alpha collision heat reward"));
                    }
                    maximum_reward_j_per_mol[bath] =
                        maximum_reward_j_per_mol[bath].max(value.abs());
                }
            }
        } else {
            staying[..grid.len()].fill(1.0);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            grid: grid.clone(),
            lambda,
            upward,
            downward,
            staying,
            reward_j_per_mol,
            maximum_reward_j_per_mol,
        })
    }

    /// Positive uniformization of exp(Q dt), with the same Poisson weights for
    /// population states and their integrated two-bath Markov heat rewards.
    /// Success-path workspace is bounded stack storage. The input is immutable;
    /// cancellation, domain errors or exhausted work publish no partial result.
    pub(crate) fn advance(
        &self,
        amounts: &[f64; MAX_ALPHA_BINS],
        duration: Time,
        max_applications: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AlphaCollisionStep, DtFusionError> {
        poll(&mut cancelled)?;
        self.grid.admit_amounts(amounts)?;
        if !duration.value().is_finite() || duration.value() < 0.0 {
            return Err(DtFusionError::Invalid(
                "finite nonnegative alpha collision duration required",
            ));
        }
        if max_applications > MAX_APPLICATIONS {
            return Err(DtFusionError::Invalid(
                "alpha collision budget must not exceed one million applications",
            ));
        }
        let initial_number = number(amounts, self.grid.len())?;
        let initial_energy = energy(amounts, &self.grid)?;
        let mut result = AlphaCollisionStep {
            amounts: *amounts,
            deposited_j: [0.0; 2],
            applications: 0,
            numerical_tail_bound: 0.0,
            heat_tail_bound_j: [0.0; 2],
            number_residual_mol: 0.0,
            energy_residual_j: 0.0,
        };
        if duration.value() == 0.0 || self.lambda == 0.0 || initial_number == 0.0 {
            poll(&mut cancelled)?;
            return Ok(result);
        }
        let exposure = self.lambda * duration.value();
        if !exposure.is_finite() {
            return Err(DtFusionError::StepBudget);
        }
        if exposure == 0.0 {
            poll(&mut cancelled)?;
            return Ok(result);
        }
        let pieces = (exposure / MAX_EXPOSURE).ceil().max(1.0);
        if pieces > max_applications as f64 {
            return Err(DtFusionError::StepBudget);
        }
        let pieces = pieces as usize;
        let mut heat = [Sum::default(); 2];
        let mut tail = Sum::default();
        let mut heat_tail = [Sum::default(); 2];
        for piece in 0..pieces {
            poll(&mut cancelled)?;
            // Index-derived endpoints retain the exact total duration without
            // accumulating a separate floating-point time drift.
            let begin = duration.value() * (piece as f64 / pieces as f64);
            let end = duration.value() * ((piece + 1) as f64 / pieces as f64);
            let z = self.lambda * (end - begin);
            if z == 0.0 {
                continue;
            }
            let poisson = Poisson::new(z)?;
            if poisson.terms > max_applications - result.applications {
                return Err(DtFusionError::StepBudget);
            }
            let mut current = result.amounts;
            let mut accumulated = [Sum::default(); MAX_ALPHA_BINS];
            for i in 0..self.grid.len() {
                accumulated[i].add(poisson.weights[0] * current[i]);
            }
            let mut path_heat = [Sum::default(); 2];
            let mut mean_heat = [Sum::default(); 2];
            for k in 1..=poisson.terms {
                poll(&mut cancelled)?;
                let (next, reward) = self.apply(&current, &mut cancelled)?;
                result.applications += 1;
                for bath in 0..2 {
                    path_heat[bath].add(reward[bath]);
                    mean_heat[bath].add(poisson.weights[k] * path_heat[bath].value);
                }
                for i in 0..self.grid.len() {
                    accumulated[i].add(poisson.weights[k] * next[i]);
                }
                current = next;
            }
            result.amounts = accumulated.map(|sum| sum.value);
            for bath in 0..2 {
                heat[bath].add(mean_heat[bath].value);
                // Upper tail mass is lumped on P^K; its excess over the exact
                // tail can also move weight from P^0. Include both contributions
                // to the linear-growth reward bound, not just excess event count.
                // Population error from this interval also affects subsequent
                // rewards: its L1 bound is 2*N*tail, propagated through a
                // stochastic semigroup for the remaining expected event count.
                let remaining_exposure = self.lambda * (duration.value() - end);
                heat_tail[bath].add(
                    initial_number
                        * (self.maximum_reward_j_per_mol[bath]
                            * (poisson.terms as f64 * poisson.tail_bound
                                + poisson.excess_events_bound
                                + 2.0 * remaining_exposure * poisson.tail_bound)),
                );
            }
            tail.add(poisson.tail_bound);
        }
        self.grid.admit_amounts(&result.amounts)?;
        result.deposited_j = heat.map(|sum| sum.value);
        result.numerical_tail_bound = tail.value;
        result.heat_tail_bound_j = heat_tail.map(|sum| sum.value);
        result.number_residual_mol = number(&result.amounts, self.grid.len())? - initial_number;
        let final_energy = energy(&result.amounts, &self.grid)?;
        result.energy_residual_j =
            (final_energy - initial_energy) + result.deposited_j[0] + result.deposited_j[1];
        let rounding_allowance = 64.0 * f64::EPSILON * (result.applications as f64 + 1.0);
        let energy_scale = initial_energy
            .abs()
            .max(final_energy.abs())
            .max(result.deposited_j[0].abs())
            .max(result.deposited_j[1].abs());
        if result
            .deposited_j
            .iter()
            .chain(result.heat_tail_bound_j.iter())
            .chain([
                &result.number_residual_mol,
                &result.energy_residual_j,
                &result.numerical_tail_bound,
            ])
            .any(|value| !value.is_finite())
            || result.number_residual_mol.abs() > rounding_allowance * initial_number
            || result.energy_residual_j.abs()
                > (4.0 * rounding_allowance) * energy_scale
                    + result.heat_tail_bound_j[0]
                    + result.heat_tail_bound_j[1]
        {
            return Err(DtFusionError::Balance("alpha collision number or energy"));
        }
        poll(&mut cancelled)?;
        Ok(result)
    }

    fn apply(
        &self,
        current: &[f64; MAX_ALPHA_BINS],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<([f64; MAX_ALPHA_BINS], [f64; 2]), DtFusionError> {
        let mut next = [0.0; MAX_ALPHA_BINS];
        let mut heat = [Sum::default(); 2];
        for i in 0..self.grid.len() {
            if i % 32 == 0 {
                poll(cancelled)?;
            }
            let mut population = Sum::default();
            population.add(self.staying[i] * current[i]);
            if i > 0 {
                population.add(self.upward[i - 1] * current[i - 1]);
            }
            if i + 1 < self.grid.len() {
                population.add(self.downward[i + 1] * current[i + 1]);
            }
            if !population.value.is_finite() || population.value < 0.0 {
                return Err(DtFusionError::Arithmetic(
                    "positive alpha Markov transition",
                ));
            }
            next[i] = population.value;
            for bath in 0..2 {
                heat[bath].add(self.reward_j_per_mol[bath][i] * current[i]);
            }
        }
        Ok((next, heat.map(|sum| sum.value)))
    }
}

struct Poisson {
    weights: [f64; MAX_POISSON_TERMS],
    terms: usize,
    tail_bound: f64,
    excess_events_bound: f64,
}

impl Poisson {
    fn new(z: f64) -> Result<Self, DtFusionError> {
        // Subintervals target z<=8. Index-derived floating-point endpoints can
        // exceed that target by clock roundoff; retain a strict finite bound
        // without clipping their represented physical durations.
        if !positive(z) || z > 2.0 * MAX_EXPOSURE {
            return Err(DtFusionError::Arithmetic("bounded alpha Poisson exposure"));
        }
        let mut weights = [0.0; MAX_POISSON_TERMS];
        let mut p = exp(-z);
        let mut nonzero_weights = Sum::default();
        for k in 1..MAX_POISSON_TERMS {
            p *= z / k as f64;
            weights[k] = p;
            nonzero_weights.add(p);
            let next = p * (z / (k + 1) as f64);
            let ratio = z / (k + 2) as f64;
            if ratio < 1.0 {
                let tail_bound = next / (1.0 - ratio);
                if tail_bound <= POISSON_TAIL_TOLERANCE {
                    // Use a positive upper bound on the remaining tail rather
                    // than 1-sum(p_k), which can round negative. Assign its mass
                    // to the final Markov state and form p0 as the complement.
                    // The bounded exposure keeps this complement well resolved.
                    weights[k] += tail_bound;
                    nonzero_weights.add(tail_bound);
                    weights[0] = 1.0 - nonzero_weights.value;
                    if !positive(weights[0]) {
                        return Err(DtFusionError::Arithmetic("positive alpha Poisson weights"));
                    }
                    return Ok(Self {
                        weights,
                        terms: k,
                        tail_bound,
                        excess_events_bound: next / ((1.0 - ratio) * (1.0 - ratio)),
                    });
                }
            }
        }
        Err(DtFusionError::StepBudget)
    }
}

fn edge_rate(diffusion: f64, gap: f64, width: f64, w: f64) -> Result<f64, DtFusionError> {
    if diffusion == 0.0 {
        return Ok(0.0);
    }
    let factor = bernoulli(w);
    let direct = (diffusion / gap) * (factor / width);
    if positive(direct) {
        return Ok(direct);
    }
    // A tiny Bernoulli factor can compensate a large prefactor. Combine their
    // logarithms before exponentiating when separated arithmetic loses it.
    let log_factor = if w > 50.0 {
        ln(w) - w
    } else if w < -50.0 {
        ln(-w)
    } else {
        ln(factor)
    };
    let rate = exp(ln(diffusion) - ln(gap) - ln(width) + log_factor);
    if rate.is_finite() && rate >= 0.0 {
        Ok(rate)
    } else {
        Err(DtFusionError::Arithmetic(
            "representable alpha collision rate",
        ))
    }
}

fn bernoulli(w: f64) -> f64 {
    if w.abs() <= 1e-4 {
        let squared = w * w;
        1.0 - 0.5 * w + squared * (1.0 / 12.0 + squared * (-1.0 / 720.0 + squared / 30240.0))
    } else if w > 50.0 {
        let exponential = exp(-w);
        w * exponential / (1.0 - exponential)
    } else if w < -50.0 {
        -w / (1.0 - exp(w))
    } else {
        w / expm1(w)
    }
}

fn number(amounts: &[f64; MAX_ALPHA_BINS], len: usize) -> Result<f64, DtFusionError> {
    let mut total = Sum::default();
    for amount in &amounts[..len] {
        total.add(*amount);
    }
    if total.value.is_finite() {
        Ok(total.value)
    } else {
        Err(DtFusionError::Arithmetic("finite alpha population sum"))
    }
}

fn energy(amounts: &[f64; MAX_ALPHA_BINS], grid: &AlphaEnergyGrid) -> Result<f64, DtFusionError> {
    let mut total = Sum::default();
    for (amount, energy) in amounts.iter().zip(grid.energies()) {
        total.add((amount * energy.value()) * AVOGADRO);
    }
    if total.value.is_finite() {
        Ok(total.value)
    } else {
        Err(DtFusionError::Arithmetic("finite alpha kinetic energy sum"))
    }
}

#[cfg(test)]
mod tests;
