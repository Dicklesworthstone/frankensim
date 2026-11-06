//! Transient photon recycling between equal-energy luminescent stores.
//!
//! Every admitted photon has the same energy as every receiving excitation.
//! The caller must enforce this resonant condition and supply all receiving
//! material cells, including cells with no external pump. In joules,
//! `U_j' = P_j - U_j/tau_j + sum_i(phi_i R_ij U_i/tau_i)`.
//! The existing positive reaction/source owner integrates both the stores and
//! repeated emissions throughout the interval, including self-return cycles.
//! No steady-state inverse, leakage floor or finite bounce truncation is used.
//!
//! Rates and transport fractions stay fixed over one prepared interval. Photon
//! flight time is neglected on this kinetic clock; this model does not include
//! wavelength conversion, saturation, delayed transport or coherent emission.
//! Geometry, spectral admission and changing-temperature rates belong to the
//! caller. Numerical resources retain the reaction owner's 64-cell bound.

use super::ExcitationDecay;
use core::fmt;
use fs_la::decay::reaction::{ReactionEdge, source::ReactionSourcePlan};
use fs_la::decay::{DecayError, MAX_DECAY_EDGES, MAX_DECAY_NODES};
use fs_math::eft::two_sum;

/// Existing numerical work, size and accuracy limits, without another dependency
/// in clients of the thermal owner.
pub use fs_la::decay::DecayLimits as RecyclingLimits;

/// A directed absorption fraction per emitted photon. Self-return is permitted;
/// duplicate endpoints add without renormalization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RecyclingEdge {
    pub source_cell: usize,
    pub receiver_cell: usize,
    pub absorbed_fraction: f64,
}

/// Physical admission or a refusal from the existing positive evolution owner.
#[derive(Clone, Debug, PartialEq)]
pub enum RecyclingError {
    Invalid(&'static str),
    Unrepresentable(&'static str),
    Evolution(DecayError),
    /// Internal reemissions do not contribute to this external energy allowance.
    EnergyBalance {
        residual_j: f64,
        allowance_j: f64,
    },
}
impl fmt::Display for RecyclingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "resonant luminescence recycling: {self:?}")
    }
}
impl std::error::Error for RecyclingError {}
impl From<DecayError> for RecyclingError {
    fn from(value: DecayError) -> Self {
        Self::Evolution(value)
    }
}

/// Reusable equal-energy transport and decay interval, independent of inventory.
#[derive(Clone, Debug)]
pub struct ResonantRecyclingPlan {
    duration_s: f64,
    relative_tolerance: f64,
    quantum_yields: Vec<f64>,
    links: Vec<RecyclingEdge>,
    escape_fractions: Vec<f64>,
    unresolved_fractions: Vec<f64>,
    evolution: ReactionSourcePlan,
}

/// Integrated channels for a candidate interval. Every extensive quantity is J.
#[derive(Clone, Debug, PartialEq)]
pub struct RecyclingStep {
    pub stored_energy_j: Vec<f64>,
    pub pumped_energy_j: Vec<f64>,
    /// Emission by source, including photons that are subsequently recycled.
    /// This can exceed the external energy supplied to a closed recycling system.
    pub emitted_energy_j: Vec<f64>,
    pub nonradiative_heat_j: Vec<f64>,
    /// Integrated incoming photon energy by receiver, including repeated visits.
    pub reabsorbed_energy_j: Vec<f64>,
    pub escaped_energy_j: f64,
    /// Declared unresolved outgoing transport; it is never assigned to heat.
    pub unresolved_energy_j: f64,
    /// Initial store + pump - final store - dark heat - escape - unresolved.
    pub residual_j: f64,
}

impl ResonantRecyclingPlan {
    /// Prepare an interval whose absorption, escape and unresolved fractions
    /// close each emitter's outgoing partition within `fraction_roundoff`.
    /// That explicit nonnegative allowance must be smaller than one; fractions
    /// are never rescaled. Positive sinks erased by the feedback representation
    /// refuse. Numerical admission of each candidate uses `limits` separately.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        decays: &[ExcitationDecay],
        absorption_edges: &[RecyclingEdge],
        escape_fractions: &[f64],
        unresolved_fractions: &[f64],
        duration_s: f64,
        fraction_roundoff: f64,
        limits: RecyclingLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RecyclingError> {
        poll(&mut cancelled)?;
        let n = decays.len();
        cap(
            "recycling cells",
            n,
            limits.maximum_nodes.min(MAX_DECAY_NODES),
        )?;
        cap(
            "recycling links",
            absorption_edges.len(),
            limits.maximum_edges.min(MAX_DECAY_EDGES),
        )?;
        if n == 0
            || escape_fractions.len() != n
            || unresolved_fractions.len() != n
            || !nonnegative(duration_s)
            || !nonnegative(fraction_roundoff)
            || fraction_roundoff >= 1.0
            || !escape_fractions
                .iter()
                .chain(unresolved_fractions)
                .copied()
                .all(fraction)
        {
            return Err(RecyclingError::Invalid(
                "cells, duration or transport fractions",
            ));
        }
        let mut rates = Vec::with_capacity(n);
        let mut quantum_yields = Vec::with_capacity(n);
        let mut partition_errors = vec![Sum::default(); n];
        // Aggregate duplicate links in the same receiver/source order used by
        // the numerical owner before assessing its represented net sink.
        let mut feedback = vec![Sum::default(); n * n];
        for (i, decay) in decays.iter().enumerate() {
            poll(&mut cancelled)?;
            rates.push(finite(1.0 / decay.lifetime_s(), "inverse lifetime")?);
            quantum_yields.push(decay.quantum_yield());
            partition_errors[i].add(-1.0)?;
            partition_errors[i].add(escape_fractions[i])?;
            partition_errors[i].add(unresolved_fractions[i])?;
        }
        for link in absorption_edges {
            poll(&mut cancelled)?;
            if link.source_cell >= n || link.receiver_cell >= n || !fraction(link.absorbed_fraction)
            {
                return Err(RecyclingError::Invalid("absorption endpoint or fraction"));
            }
            partition_errors[link.source_cell].add(link.absorbed_fraction)?;
            let yield_per_removal = product(
                quantum_yields[link.source_cell],
                link.absorbed_fraction,
                "recycling yield",
            )?;
            feedback[link.receiver_cell * n + link.source_cell].add(yield_per_removal)?;
        }
        let mut edges = Vec::with_capacity(absorption_edges.len());
        for i in 0..n {
            poll(&mut cancelled)?;
            if partition_errors[i].total()?.abs() > fraction_roundoff {
                return Err(RecyclingError::Invalid(
                    "outgoing photon fractions do not sum to one",
                ));
            }
            let mut net_yield = Sum::default();
            net_yield.add(-1.0)?;
            for j in 0..n {
                let yield_per_removal = feedback[j * n + i].total()?;
                net_yield.add(yield_per_removal)?;
                if yield_per_removal > 0.0 {
                    edges.push(ReactionEdge {
                        parent: i,
                        daughter: j,
                        yield_per_removal,
                    });
                }
            }
            let phi = quantum_yields[i];
            // External channels are applied to resolved emitted joules later.
            // Do not underflow an unnecessary phi*escape intermediate here:
            // a large store can still make that outgoing energy representable.
            let positive_loss = phi < 1.0
                || (phi > 0.0 && (escape_fractions[i] > 0.0 || unresolved_fractions[i] > 0.0));
            if positive_loss && net_yield.total()? >= 0.0 {
                return Err(RecyclingError::Unrepresentable(
                    "positive loss erased by recycling yields",
                ));
            }
        }
        let evolution =
            ReactionSourcePlan::prepare(&rates, &edges, duration_s, limits, &mut cancelled)?;
        Ok(Self {
            duration_s,
            relative_tolerance: limits.relative_tolerance,
            quantum_yields,
            links: absorption_edges.to_vec(),
            escape_fractions: escape_fractions.to_vec(),
            unresolved_fractions: unresolved_fractions.to_vec(),
            evolution,
        })
    }

    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    pub fn cell_count(&self) -> usize {
        self.quantum_yields.len()
    }

    /// Evolve all material stores under constant admitted excitation power [W].
    pub fn advance(
        &self,
        stored_energy_j: &[f64],
        pump_power_w: &[f64],
    ) -> Result<RecyclingStep, RecyclingError> {
        self.advance_with_cancel(stored_energy_j, pump_power_w, || false)
    }

    /// Transactional candidate; cancelled or inadmissible states are not returned.
    pub fn advance_with_cancel(
        &self,
        stored_energy_j: &[f64],
        pump_power_w: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RecyclingStep, RecyclingError> {
        poll(&mut cancelled)?;
        let n = self.cell_count();
        if stored_energy_j.len() != n
            || pump_power_w.len() != n
            || !stored_energy_j
                .iter()
                .chain(pump_power_w)
                .copied()
                .all(nonnegative)
        {
            return Err(RecyclingError::Invalid(
                "one finite nonnegative store and pump per cell required",
            ));
        }
        let mut pumped = Vec::with_capacity(n);
        for &power in pump_power_w {
            poll(&mut cancelled)?;
            pumped.push(product(power, self.duration_s, "pumped energy")?);
        }
        let evolved = self
            .evolution
            .apply_with_cancel(stored_energy_j, &pumped, &mut cancelled)?;
        let mut step = RecyclingStep {
            stored_energy_j: evolved.populations,
            pumped_energy_j: pumped,
            emitted_energy_j: vec![0.0; n],
            nonradiative_heat_j: vec![0.0; n],
            reabsorbed_energy_j: vec![0.0; n],
            escaped_energy_j: 0.0,
            unresolved_energy_j: 0.0,
            residual_j: 0.0,
        };
        let mut escape = Sum::default();
        let mut unresolved = Sum::default();
        for i in 0..n {
            poll(&mut cancelled)?;
            let phi = self.quantum_yields[i];
            step.emitted_energy_j[i] = product(evolved.removals[i], phi, "emitted energy")?;
            step.nonradiative_heat_j[i] =
                product(evolved.removals[i], 1.0 - phi, "nonradiative heat")?;
            escape.add(product(
                step.emitted_energy_j[i],
                self.escape_fractions[i],
                "escaped energy",
            )?)?;
            unresolved.add(product(
                step.emitted_energy_j[i],
                self.unresolved_fractions[i],
                "unresolved energy",
            )?)?;
        }
        let mut incoming = vec![Sum::default(); n];
        for link in &self.links {
            poll(&mut cancelled)?;
            incoming[link.receiver_cell].add(product(
                step.emitted_energy_j[link.source_cell],
                link.absorbed_fraction,
                "reabsorbed energy",
            )?)?;
        }
        for (i, sum) in incoming.into_iter().enumerate() {
            step.reabsorbed_energy_j[i] = sum.total()?;
        }
        step.escaped_energy_j = escape.total()?;
        step.unresolved_energy_j = unresolved.total()?;
        step.residual_j = admit_energy_balance(stored_energy_j, &step, self.relative_tolerance)?;
        poll(&mut cancelled)?;
        Ok(step)
    }
}

// Internal photon traffic is excluded: a nearly closed, rapidly cycling graph
// must not use its much larger emission counts to hide errors in net energy.
fn admit_energy_balance(
    initial: &[f64],
    step: &RecyclingStep,
    tolerance: f64,
) -> Result<f64, RecyclingError> {
    let mut residual = Sum::default();
    let mut scale = Sum::default();
    for (i, &old) in initial.iter().enumerate() {
        for term in [
            old,
            -step.stored_energy_j[i],
            step.pumped_energy_j[i],
            -step.nonradiative_heat_j[i],
        ] {
            residual.add(term)?;
            scale.add(term.abs())?;
        }
    }
    for loss in [step.escaped_energy_j, step.unresolved_energy_j] {
        residual.add(-loss)?;
        scale.add(loss)?;
    }
    let residual_j = residual.total()?;
    let allowance_j = finite(tolerance * scale.total()?, "energy balance allowance")?;
    if residual_j.abs() > allowance_j {
        Err(RecyclingError::EnergyBalance {
            residual_j,
            allowance_j,
        })
    } else {
        Ok(residual_j)
    }
}

fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn fraction(value: f64) -> bool {
    nonnegative(value) && value <= 1.0
}
fn finite(value: f64, quantity: &'static str) -> Result<f64, RecyclingError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RecyclingError::Unrepresentable(quantity))
    }
}
fn product(a: f64, b: f64, quantity: &'static str) -> Result<f64, RecyclingError> {
    let result = finite(a * b, quantity)?;
    if a > 0.0 && b > 0.0 && result == 0.0 {
        Err(RecyclingError::Unrepresentable(quantity))
    } else {
        Ok(result)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RecyclingError> {
    if cancelled() {
        Err(DecayError::Cancelled.into())
    } else {
        Ok(())
    }
}
fn cap(quantity: &'static str, required: usize, maximum: usize) -> Result<(), RecyclingError> {
    if required > maximum {
        Err(DecayError::Budget {
            quantity,
            required: required as u64,
            maximum: maximum as u64,
        }
        .into())
    } else {
        Ok(())
    }
}
#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), RecyclingError> {
        let (value, error) = two_sum(self.value, term);
        self.value = finite(value, "partial compensated sum")?;
        self.correction = finite(self.correction + error, "compensated sum correction")?;
        Ok(())
    }
    fn total(self) -> Result<f64, RecyclingError> {
        finite(self.value + self.correction, "compensated sum")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn decay(tau: f64, phi: f64) -> ExcitationDecay {
        ExcitationDecay::try_new(tau, phi).unwrap()
    }
    fn edge(source_cell: usize, receiver_cell: usize, absorbed_fraction: f64) -> RecyclingEdge {
        RecyclingEdge {
            source_cell,
            receiver_cell,
            absorbed_fraction,
        }
    }
    fn prepare(
        decays: &[ExcitationDecay],
        links: &[RecyclingEdge],
        escape: &[f64],
        unresolved: &[f64],
        h: f64,
    ) -> ResonantRecyclingPlan {
        ResonantRecyclingPlan::prepare(
            decays,
            links,
            escape,
            unresolved,
            h,
            32.0 * f64::EPSILON,
            RecyclingLimits::default(),
            || false,
        )
        .unwrap()
    }
    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 3e-11 * expected.abs().max(1e-300),
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn g0_closed_self_return_retains_store_and_integrates_all_emissions() {
        let p = prepare(&[decay(2.0, 1.0)], &[edge(0, 0, 1.0)], &[0.0], &[0.0], 3.0);
        let stopped = p.advance(&[4.0], &[0.0]).unwrap();
        close(stopped.stored_energy_j[0], 4.0);
        close(stopped.emitted_energy_j[0], 6.0);
        let pumped = p.advance(&[4.0], &[2.0]).unwrap();
        close(pumped.stored_energy_j[0], 10.0);
        close(pumped.pumped_energy_j[0], 6.0);
        close(pumped.emitted_energy_j[0], 10.5);
        close(pumped.reabsorbed_energy_j[0], 10.5);
        assert_eq!(pumped.nonradiative_heat_j, [0.0]);
        assert_eq!(pumped.escaped_energy_j, 0.0);
        assert_eq!(pumped.unresolved_energy_j, 0.0);
        assert!(pumped.residual_j.abs() < 1e-12);
    }

    #[test]
    fn g0_two_cell_exchange_matches_independent_closed_form() {
        let h: f64 = 0.7;
        let p = prepare(
            &[decay(1.0, 1.0); 2],
            &[edge(0, 1, 1.0), edge(1, 0, 1.0)],
            &[0.0; 2],
            &[0.0; 2],
            h,
        );
        let s = p.advance(&[4.0, 0.0], &[0.0; 2]).unwrap();
        let e = (-2.0 * h).exp();
        close(s.stored_energy_j[0], 2.0 + 2.0 * e);
        close(s.stored_energy_j[1], 2.0 - 2.0 * e);
        close(s.emitted_energy_j[0], 2.0 * h + 1.0 - e);
        close(s.emitted_energy_j[1], 2.0 * h - 1.0 + e);
        close(s.reabsorbed_energy_j[0], s.emitted_energy_j[1]);
        close(s.reabsorbed_energy_j[1], s.emitted_energy_j[0]);
        close(s.stored_energy_j.iter().sum(), 4.0);
        assert_eq!(s.escaped_energy_j, 0.0);
    }

    #[test]
    fn g0_self_recycling_with_pump_keeps_escape_dark_and_unresolved_separate() {
        let h: f64 = 1.5;
        let p = prepare(&[decay(2.0, 0.8)], &[edge(0, 0, 0.6)], &[0.3], &[0.1], h);
        let s = p.advance(&[5.0], &[0.4]).unwrap();
        let net_rate = 0.5 * (1.0 - 0.8 * 0.6);
        let equilibrium = 0.4 / net_rate;
        let e = (-net_rate * h).exp();
        let removed = 0.5 * (equilibrium * h + (5.0 - equilibrium) * (1.0 - e) / net_rate);
        close(s.stored_energy_j[0], equilibrium + (5.0 - equilibrium) * e);
        close(s.emitted_energy_j[0], 0.8 * removed);
        close(s.nonradiative_heat_j[0], 0.2 * removed);
        close(s.reabsorbed_energy_j[0], 0.6 * 0.8 * removed);
        close(s.escaped_energy_j, 0.3 * 0.8 * removed);
        close(s.unresolved_energy_j, 0.1 * 0.8 * removed);
        assert!(s.residual_j.abs() < 1e-11);
    }

    #[test]
    fn g3_split_steps_match_every_integrated_recycling_channel() {
        let d = [decay(0.3, 0.8), decay(2.0, 0.6)];
        let links = [
            edge(0, 0, 0.2),
            edge(0, 1, 0.4),
            edge(1, 0, 0.5),
            edge(1, 1, 0.1),
        ];
        let all = prepare(&d, &links, &[0.3, 0.2], &[0.1, 0.2], 1.0)
            .advance(&[1.0, 3.0], &[2.0, 0.1])
            .unwrap();
        let half = prepare(&d, &links, &[0.3, 0.2], &[0.1, 0.2], 0.5);
        let a = half.advance(&[1.0, 3.0], &[2.0, 0.1]).unwrap();
        let b = half.advance(&a.stored_energy_j, &[2.0, 0.1]).unwrap();
        for i in 0..2 {
            close(all.stored_energy_j[i], b.stored_energy_j[i]);
            close(
                all.emitted_energy_j[i],
                a.emitted_energy_j[i] + b.emitted_energy_j[i],
            );
            close(
                all.nonradiative_heat_j[i],
                a.nonradiative_heat_j[i] + b.nonradiative_heat_j[i],
            );
            close(
                all.reabsorbed_energy_j[i],
                a.reabsorbed_energy_j[i] + b.reabsorbed_energy_j[i],
            );
        }
        close(
            all.escaped_energy_j,
            a.escaped_energy_j + b.escaped_energy_j,
        );
        close(
            all.unresolved_energy_j,
            a.unresolved_energy_j + b.unresolved_energy_j,
        );
    }

    #[test]
    fn g0_positive_loss_cannot_disappear_behind_recycled_energy() {
        let result = ResonantRecyclingPlan::prepare(
            &[decay(1.0, 1.0)],
            &[edge(0, 0, 1.0)],
            &[1e-20],
            &[0.0],
            1.0,
            1e-14,
            RecyclingLimits::default(),
            || false,
        );
        assert!(matches!(
            result,
            Err(RecyclingError::Unrepresentable(
                "positive loss erased by recycling yields"
            ))
        ));
        let mut s = prepare(&[decay(1.0, 1.0)], &[edge(0, 0, 1.0)], &[0.0], &[0.0], 1.0)
            .advance(&[1.0], &[0.0])
            .unwrap();
        s.emitted_energy_j[0] = 1e20;
        s.reabsorbed_energy_j[0] = 1e20;
        s.escaped_energy_j = 1.0;
        assert!(matches!(
            admit_energy_balance(&[1.0], &s, 2e-12),
            Err(RecyclingError::EnergyBalance { .. })
        ));
    }

    #[test]
    fn g0_positive_feedback_and_emission_underflow_refuse() {
        let lost_weight = ResonantRecyclingPlan::prepare(
            &[decay(1.0, 1e-300)],
            &[edge(0, 0, 1e-300)],
            &[1.0],
            &[0.0],
            1.0,
            1e-14,
            RecyclingLimits::default(),
            || false,
        );
        assert!(matches!(
            lost_weight,
            Err(RecyclingError::Unrepresentable("recycling yield"))
        ));
        let tiny = prepare(&[decay(1.0, 1e-100)], &[], &[1.0], &[0.0], 1.0);
        assert!(matches!(
            tiny.advance(&[1e-300], &[0.0]),
            Err(RecyclingError::Unrepresentable("emitted energy"))
        ));
        // An unused phi*escape coefficient may underflow even though scaling
        // the actual emitting store first resolves both outgoing channels.
        let scaled = prepare(&[decay(1.0, f64::from_bits(1))], &[], &[0.5], &[0.5], 1.0)
            .advance(&[1e300], &[0.0])
            .unwrap();
        assert!(scaled.emitted_energy_j[0] > 0.0);
        close(scaled.escaped_energy_j, 0.5 * scaled.emitted_energy_j[0]);
        close(scaled.unresolved_energy_j, 0.5 * scaled.emitted_energy_j[0]);
    }

    #[test]
    fn g0_zero_interval_shape_fraction_budget_and_cancellation_admission() {
        let d = [decay(1.0, 0.5)];
        let zero = prepare(&d, &[], &[1.0], &[0.0], 0.0)
            .advance(&[3.0], &[7.0])
            .unwrap();
        assert_eq!(zero.stored_energy_j, [3.0]);
        assert_eq!(zero.pumped_energy_j, [0.0]);
        assert_eq!(zero.emitted_energy_j, [0.0]);
        let p = prepare(&d, &[], &[1.0], &[0.0], 1.0);
        assert!(p.advance(&[], &[0.0]).is_err());
        assert!(p.advance(&[-1.0], &[0.0]).is_err());
        assert!(
            ResonantRecyclingPlan::prepare(
                &d,
                &[],
                &[0.5],
                &[0.0],
                1.0,
                1e-14,
                RecyclingLimits::default(),
                || false
            )
            .is_err()
        );
        assert!(
            ResonantRecyclingPlan::prepare(
                &d,
                &[edge(0, 1, 1.0)],
                &[0.0],
                &[0.0],
                1.0,
                1e-14,
                RecyclingLimits::default(),
                || false
            )
            .is_err()
        );
        let too_many = [decay(1.0, 1.0); MAX_DECAY_NODES + 1];
        assert!(matches!(
            ResonantRecyclingPlan::prepare(
                &too_many,
                &[],
                &[],
                &[],
                1.0,
                1e-14,
                RecyclingLimits::default(),
                || false
            ),
            Err(RecyclingError::Evolution(DecayError::Budget { .. }))
        ));
        assert!(matches!(
            ResonantRecyclingPlan::prepare(
                &d,
                &[],
                &[1.0],
                &[0.0],
                1.0,
                1e-14,
                RecyclingLimits::default(),
                || true
            ),
            Err(RecyclingError::Evolution(DecayError::Cancelled))
        ));
        assert!(matches!(
            p.advance_with_cancel(&[1.0], &[0.0], || true),
            Err(RecyclingError::Evolution(DecayError::Cancelled))
        ));
    }
}
