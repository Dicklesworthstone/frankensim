//! Coupled singlet/triplet populations and their integrated energy channels.
//!
//! Index 0 denotes the singlet and index 1 the triplet. With authored rates,
//! `Ns' = Ps - (ksr + ksnr + kst) Ns + kts Nt` and
//! `Nt' = Pt + kst Ns - (ktr + ktnr + kts) Nt`.
//! Forward intersystem crossing supplies `(Es-Et)` as heat per event; reverse
//! crossing absorbs that energy from the thermal bath. A triplet can therefore
//! feed delayed singlet emission without inventing an extra photon or its energy.
//!
//! The positive matrix exponential and integrated event counts belong to
//! `fs_la::decay::reaction`; equal rates, constant pumping and stiff intervals
//! use that existing numerical owner. Preparation is reusable across cells.
//! Rates remain constant over a prepared interval. The source owner must choose
//! thermally admissible reverse rates and resolve changing bath temperatures.
//! This dilute, first-order model excludes finite ground-state depletion,
//! annihilation, stimulated emission and coherent spin dynamics.

use core::fmt;
use fs_la::decay::reaction::{ReactionEdge, source::ReactionSourcePlan};
use fs_la::decay::{DecayError, DecayLimits};
use fs_math::eft::two_sum;

/// Spectral photon-number recycling between coupled excitation stores.
#[path = "luminescence_transfer/recycling.rs"]
pub mod recycling;

/// Exact dilute triplet-pair annihilation with explicit excitation and heat.
#[path = "luminescence_transfer/annihilation.rs"]
pub mod annihilation;

/// Physical admission or a refusal from the existing positive evolution owner.
#[derive(Clone, Debug, PartialEq)]
pub enum SingletTripletError {
    /// Nonfinite, negative or physically unordered authored data.
    Invalid(&'static str),
    /// A required count, channel fraction or energy exceeds numerical range.
    Unrepresentable(&'static str),
    /// Authored ground-state channels do not balance the represented evolution.
    PopulationBalance { residual: f64, allowance: f64 },
    /// External energy does not close without counting internal photon traffic.
    EnergyBalance { residual_j: f64, allowance_j: f64 },
    /// Numerical admission, cancellation or resource limits were not satisfied.
    Evolution(DecayError),
}
impl fmt::Display for SingletTripletError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "singlet/triplet luminescence: {self:?}")
    }
}
impl std::error::Error for SingletTripletError {}
impl From<DecayError> for SingletTripletError {
    fn from(value: DecayError) -> Self {
        Self::Evolution(value)
    }
}

/// Authored excitation energies and competing rates, ordered `[singlet, triplet]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SingletTripletDecay {
    energies_j: [f64; 2],
    radiative_rates_s: [f64; 2],
    nonradiative_rates_s: [f64; 2],
    crossing_rates_s: [f64; 2],
    total_rates_s: [f64; 2],
    fractions: [[f64; 3]; 2],
}

/// Prepared fixed-rate interval; no initial populations or pump rates are cached.
#[derive(Clone, Debug)]
pub struct SingletTripletPlan {
    model: SingletTripletDecay,
    duration_s: f64,
    relative_tolerance: f64,
    evolution: ReactionSourcePlan,
}

/// Candidate populations, integrated events and joules, ordered `[singlet, triplet]`.
#[derive(Clone, Debug, PartialEq)]
pub struct SingletTripletStep {
    /// Endpoint numbers of excitations, not concentrations or energy densities.
    pub populations: [f64; 2],
    /// Excitations admitted uniformly during the interval.
    pub injected_excitations: [f64; 2],
    /// Integrated transitions to the ground state through the radiative channel.
    pub radiative_events: [f64; 2],
    /// Integrated transitions to the ground state through the nonradiative channel.
    pub nonradiative_events: [f64; 2],
    /// Integrated `[singlet->triplet, triplet->singlet]` events, including revisits.
    pub crossing_events: [f64; 2],
    /// Endpoint excitation energy in each state, joules.
    pub stored_energy_j: [f64; 2],
    /// Energy admitted to the relaxed stores, after any pump Stokes partition.
    pub pumped_energy_j: [f64; 2],
    /// State energy released radiatively BEFORE photon-spectrum/Stokes partition.
    /// Photon escape and reabsorption are the transport owner's responsibility.
    pub radiative_decay_energy_j: [f64; 2],
    /// Local heat released through each state's nonradiative ground-state channel.
    pub nonradiative_heat_j: [f64; 2],
    /// Heat supplied to the bath by downhill singlet-to-triplet transitions.
    pub crossing_heat_j: f64,
    /// Heat withdrawn from the bath by reverse, uphill transitions.
    /// Kept separate so large reversible exchanges cannot cancel invisibly.
    pub reverse_crossing_heat_absorbed_j: f64,
    /// `old_store + pumped + bath_absorbed - new_store - radiative - heat`.
    /// The signed arithmetic residual is reported without correcting a channel.
    pub residual_j: f64,
}

impl SingletTripletDecay {
    /// Admit positive finite excitation energies with `Es >= Et`, and finite
    /// nonnegative rates [1/s]. Crossing rates are `[kST, kTS]`. A state with
    /// every outgoing rate zero is an explicitly stable store. No lifetime,
    /// branching yield, reverse-rate balance or material constant is inferred.
    pub fn try_new(
        excitation_energies_j: [f64; 2],
        radiative_rates_s: [f64; 2],
        nonradiative_rates_s: [f64; 2],
        crossing_rates_s: [f64; 2],
    ) -> Result<Self, SingletTripletError> {
        if !excitation_energies_j
            .iter()
            .all(|&x| x.is_finite() && x > 0.0)
            || excitation_energies_j[0] < excitation_energies_j[1]
            || !radiative_rates_s
                .iter()
                .chain(&nonradiative_rates_s)
                .chain(&crossing_rates_s)
                .all(|&x| nonnegative(x))
        {
            return Err(SingletTripletError::Invalid("energies or competing rates"));
        }
        let mut total_rates_s = [0.0; 2];
        let mut fractions = [[0.0; 3]; 2];
        for i in 0..2 {
            let channels = [
                radiative_rates_s[i],
                nonradiative_rates_s[i],
                crossing_rates_s[i],
            ];
            total_rates_s[i] = finite(channels.iter().sum(), "total removal rate")?;
            if total_rates_s[i] != 0.0 {
                for (j, rate) in channels.into_iter().enumerate() {
                    fractions[i][j] = rate / total_rates_s[i];
                    if rate > 0.0 && fractions[i][j] == 0.0 {
                        return Err(SingletTripletError::Unrepresentable(
                            "positive channel fraction",
                        ));
                    }
                }
                // The reaction owner represents its sink as lambda*(1-yield).
                // A rounded unit crossing yield must not erase a positive
                // radiative/dark sink while those events are still reported.
                // For nonunit yields, the interval's population balance below
                // assesses the actual error at the caller's requested accuracy;
                // a rate-relative cutoff would reject harmless short intervals.
                if channels[0] + channels[1] > 0.0 && fractions[i][2] == 1.0 {
                    return Err(SingletTripletError::Unrepresentable(
                        "ground-state sink erased by crossing fraction",
                    ));
                }
            }
        }
        Ok(Self {
            energies_j: excitation_energies_j,
            radiative_rates_s,
            nonradiative_rates_s,
            crossing_rates_s,
            total_rates_s,
            fractions,
        })
    }

    /// Authored per-excitation energies, joules, ordered `[singlet, triplet]`.
    pub const fn excitation_energies_j(&self) -> [f64; 2] {
        self.energies_j
    }
    /// Authored radiative rates, inverse seconds.
    pub const fn radiative_rates_s(&self) -> [f64; 2] {
        self.radiative_rates_s
    }
    /// Authored nonradiative rates, inverse seconds.
    pub const fn nonradiative_rates_s(&self) -> [f64; 2] {
        self.nonradiative_rates_s
    }
    /// Authored `[singlet->triplet, triplet->singlet]` rates, inverse seconds.
    pub const fn crossing_rates_s(&self) -> [f64; 2] {
        self.crossing_rates_s
    }

    /// Prepare the same interval once for many independent cells. Cancellation,
    /// work limits and numerical refusals retain the reaction owner's meaning.
    pub fn prepare(
        &self,
        duration_s: f64,
        limits: DecayLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SingletTripletPlan, SingletTripletError> {
        let mut edges = Vec::with_capacity(2);
        for i in 0..2 {
            if self.crossing_rates_s[i] > 0.0 {
                edges.push(ReactionEdge {
                    parent: i,
                    daughter: 1 - i,
                    yield_per_removal: self.fractions[i][2],
                });
            }
        }
        let evolution = ReactionSourcePlan::prepare(
            &self.total_rates_s,
            &edges,
            duration_s,
            limits,
            cancelled,
        )?;
        Ok(SingletTripletPlan {
            model: *self,
            duration_s,
            relative_tolerance: limits.relative_tolerance,
            evolution,
        })
    }
}

impl SingletTripletPlan {
    /// Interval length, seconds.
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Evolve fixed, dilute excitation counts under uniform production [1/s].
    /// Pump admission must account for the original photon's energy separately.
    pub fn advance(
        &self,
        populations: [f64; 2],
        pump_rates_s: [f64; 2],
    ) -> Result<SingletTripletStep, SingletTripletError> {
        self.advance_with_cancel(populations, pump_rates_s, || false)
    }

    /// Compute a candidate without mutation, with bounded reaction checkpoints.
    pub fn advance_with_cancel(
        &self,
        populations: [f64; 2],
        pump_rates_s: [f64; 2],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SingletTripletStep, SingletTripletError> {
        if !populations
            .iter()
            .chain(&pump_rates_s)
            .all(|&x| nonnegative(x))
        {
            return Err(SingletTripletError::Invalid(
                "populations or excitation source rates",
            ));
        }
        let mut injected = [0.0; 2];
        for i in 0..2 {
            injected[i] = finite(pump_rates_s[i] * self.duration_s, "injected count")?;
            if pump_rates_s[i] > 0.0 && self.duration_s > 0.0 && injected[i] == 0.0 {
                return Err(SingletTripletError::Unrepresentable(
                    "positive injected count",
                ));
            }
        }
        let result = self
            .evolution
            .apply_with_cancel(&populations, &injected, cancelled)?;
        let mut step = SingletTripletStep {
            populations: [result.populations[0], result.populations[1]],
            injected_excitations: injected,
            radiative_events: [0.0; 2],
            nonradiative_events: [0.0; 2],
            crossing_events: [0.0; 2],
            stored_energy_j: [0.0; 2],
            pumped_energy_j: [0.0; 2],
            radiative_decay_energy_j: [0.0; 2],
            nonradiative_heat_j: [0.0; 2],
            crossing_heat_j: 0.0,
            reverse_crossing_heat_absorbed_j: 0.0,
            residual_j: 0.0,
        };
        let mut balance = CompensatedSum::default();
        for i in 0..2 {
            let energy = self.model.energies_j[i];
            let fractions = self.model.fractions[i];
            step.radiative_events[i] =
                positive_product(result.removals[i], fractions[0], "radiative events")?;
            step.nonradiative_events[i] =
                positive_product(result.removals[i], fractions[1], "nonradiative events")?;
            step.crossing_events[i] =
                positive_product(result.removals[i], fractions[2], "crossing events")?;
            step.stored_energy_j[i] =
                positive_product(step.populations[i], energy, "stored energy")?;
            step.pumped_energy_j[i] = positive_product(injected[i], energy, "pumped energy")?;
            step.radiative_decay_energy_j[i] =
                positive_product(step.radiative_events[i], energy, "radiative energy")?;
            step.nonradiative_heat_j[i] =
                positive_product(step.nonradiative_events[i], energy, "nonradiative energy")?;
            for term in [
                positive_product(populations[i], energy, "initial energy")?,
                -step.stored_energy_j[i],
                step.pumped_energy_j[i],
                -step.radiative_decay_energy_j[i],
                -step.nonradiative_heat_j[i],
            ] {
                balance.add(term)?;
            }
        }
        admit_population_balance(
            &populations,
            &injected,
            &step.populations,
            &step.radiative_events,
            &step.nonradiative_events,
            self.relative_tolerance,
        )?;
        let gap = self.model.energies_j[0] - self.model.energies_j[1];
        step.crossing_heat_j = positive_product(step.crossing_events[0], gap, "crossing heat")?;
        step.reverse_crossing_heat_absorbed_j =
            positive_product(step.crossing_events[1], gap, "reverse crossing heat")?;
        balance.add(-step.crossing_heat_j)?;
        balance.add(step.reverse_crossing_heat_absorbed_j)?;
        step.residual_j = finite(balance.value + balance.correction, "energy residual")?;
        Ok(step)
    }
}

// Crossing events cancel exactly in this physical identity. Neither reversible
// traffic nor its enormous removal totals enters the admission scale: otherwise
// a nearly closed cycle could hide artificial ground-state photons in roundoff.
// This uses the same residual <= relative_tolerance * sum(abs(terms)) convention
// as the reaction owner, applied to the independently authored physical sinks.
fn admit_population_balance(
    initial: &[f64; 2],
    injected: &[f64; 2],
    final_populations: &[f64; 2],
    radiative: &[f64; 2],
    nonradiative: &[f64; 2],
    tolerance: f64,
) -> Result<(), SingletTripletError> {
    let mut residual = CompensatedSum::default();
    let mut scale = CompensatedSum::default();
    for i in 0..2 {
        for term in [
            final_populations[i],
            -initial[i],
            -injected[i],
            radiative[i],
            nonradiative[i],
        ] {
            residual.add(term)?;
            scale.add(term.abs())?;
        }
    }
    let residual = finite(residual.value + residual.correction, "population residual")?;
    let allowance = finite(
        tolerance * (scale.value + scale.correction),
        "population allowance",
    )?;
    if residual.abs() > allowance {
        Err(SingletTripletError::PopulationBalance {
            residual,
            allowance,
        })
    } else {
        Ok(())
    }
}

fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn finite(value: f64, quantity: &'static str) -> Result<f64, SingletTripletError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SingletTripletError::Unrepresentable(quantity))
    }
}
// Public event counts and energies must survive their own multiplication; a
// later large energy factor cannot justify silently returning a zero count.
fn positive_product(a: f64, b: f64, quantity: &'static str) -> Result<f64, SingletTripletError> {
    let product = finite(a * b, quantity)?;
    if a > 0.0 && b > 0.0 && product == 0.0 {
        Err(SingletTripletError::Unrepresentable(quantity))
    } else {
        Ok(product)
    }
}
#[derive(Default)]
struct CompensatedSum {
    value: f64,
    correction: f64,
}
impl CompensatedSum {
    fn add(&mut self, term: f64) -> Result<(), SingletTripletError> {
        let (value, error) = two_sum(self.value, term);
        self.value = finite(value, "partial compensated sum")?;
        self.correction = finite(self.correction + error, "compensated sum correction")?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn model(
        radiative: [f64; 2],
        nonradiative: [f64; 2],
        crossing: [f64; 2],
    ) -> SingletTripletDecay {
        SingletTripletDecay::try_new([4e-19, 3e-19], radiative, nonradiative, crossing).unwrap()
    }
    fn plan(model: SingletTripletDecay, h: f64) -> SingletTripletPlan {
        model.prepare(h, DecayLimits::default(), || false).unwrap()
    }
    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
            "actual={actual:e}, expected={expected:e}"
        );
    }
    fn balanced(step: &SingletTripletStep, scale: f64) {
        assert!(step.residual_j.abs() <= 3e-11 * scale, "{step:?}");
    }

    #[test]
    fn g0_crossing_and_phosphorescence_follow_independent_bateman_solution() {
        let h: f64 = 0.7;
        let step = plan(model([2.0, 0.1], [1.0, 0.1], [3.0, 0.0]), h)
            .advance([1.0, 0.0], [0.0; 2])
            .unwrap();
        let ns = (-6.0 * h).exp();
        let nt = 3.0 * ((-0.2 * h).exp() - ns) / 5.8;
        let crossings = 0.5 * (1.0 - ns);
        close(step.populations[0], ns, 3e-12);
        close(step.populations[1], nt, 3e-12);
        close(step.crossing_events[0], crossings, 3e-12);
        close(step.radiative_events[0], (1.0 - ns) / 3.0, 3e-12);
        close(step.radiative_events[1], 0.5 * (crossings - nt), 3e-12);
        close(step.crossing_heat_j, crossings * 1e-19, 3e-12);
        assert_eq!(step.reverse_crossing_heat_absorbed_j, 0.0);
        balanced(&step, 4e-19);
    }

    #[test]
    fn g1_equal_rates_and_tiny_intervals_preserve_positive_triplet_channels() {
        let m = model([1.0, 2.0], [0.0; 2], [1.0, 0.0]);
        for h in [0.4_f64, 1e-12] {
            let step = plan(m, h).advance([1.0, 0.0], [0.0; 2]).unwrap();
            close(step.populations[1], h * (-2.0 * h).exp(), 3e-12);
            assert!(step.radiative_events[1] > 0.0);
            if h < 1e-10 {
                close(step.radiative_events[1], h * h, 3e-12);
            }
            balanced(&step, 4e-19);
        }
    }

    #[test]
    fn g0_uniform_pump_has_continuous_production_ages() {
        let h: f64 = 2.0;
        let step = plan(model([2.0, 0.1], [1.0, 0.1], [3.0, 0.0]), h)
            .advance([0.0; 2], [5.0, 0.0])
            .unwrap();
        let es = (-6.0 * h).exp();
        let et = (-0.2 * h).exp();
        let ns = (5.0 / 6.0) * (1.0 - es);
        let nt = 2.5 * ((1.0 - et) / 0.2 - (et - es) / 5.8);
        let crossings = 2.5 * (h - (1.0 - es) / 6.0);
        close(step.populations[0], ns, 3e-12);
        close(step.populations[1], nt, 3e-12);
        close(step.crossing_events[0], crossings, 3e-12);
        close(step.radiative_events[1], 0.5 * (crossings - nt), 3e-12);
        close(step.pumped_energy_j[0], 5.0 * h * 4e-19, 3e-12);
        balanced(&step, 4e-18);
    }

    #[test]
    fn g0_reverse_crossing_pays_for_delayed_singlet_emission() {
        let step = plan(model([1.0, 0.0], [0.0; 2], [0.0, 2.0]), 1.0)
            .advance([0.0, 1.0], [0.0; 2])
            .unwrap();
        let nt = (-2.0_f64).exp();
        let ns = 2.0 * ((-1.0_f64).exp() - nt);
        close(step.populations[0], ns, 3e-12);
        close(step.crossing_events[1], 1.0 - nt, 3e-12);
        close(step.radiative_events[0], 1.0 - ns - nt, 3e-12);
        close(
            step.reverse_crossing_heat_absorbed_j,
            (1.0 - nt) * 1e-19,
            3e-12,
        );
        assert_eq!(step.crossing_heat_j, 0.0);
        assert_eq!(step.radiative_events[1], 0.0);
        balanced(&step, 3e-19);
    }

    #[test]
    fn g0_reversible_crossing_retains_both_thermal_exchanges() {
        let step = plan(model([0.0; 2], [0.0; 2], [1.0, 3.0]), 2.0)
            .advance([1.0, 0.0], [0.0; 2])
            .unwrap();
        close(step.populations[0], 0.75 + 0.25 * (-8.0_f64).exp(), 3e-12);
        close(step.populations.iter().sum(), 1.0, 3e-12);
        assert!(step.crossing_heat_j > 0.0 && step.reverse_crossing_heat_absorbed_j > 0.0);
        assert_eq!(step.radiative_events, [0.0; 2]);
        balanced(&step, 4e-19);
    }

    #[test]
    fn g3_constant_pump_split_steps_match_state_and_every_integrated_channel() {
        let m = model([2.0, 0.03], [0.5, 0.01], [3.0, 0.02]);
        let whole = plan(m, 3.0).advance([8.0, 2.0], [4.0, 0.2]).unwrap();
        let half = plan(m, 1.5);
        let a = half.advance([8.0, 2.0], [4.0, 0.2]).unwrap();
        let b = half.advance(a.populations, [4.0, 0.2]).unwrap();
        for i in 0..2 {
            close(whole.populations[i], b.populations[i], 1e-11);
            for (all, left, right) in [
                (
                    whole.radiative_events[i],
                    a.radiative_events[i],
                    b.radiative_events[i],
                ),
                (
                    whole.nonradiative_events[i],
                    a.nonradiative_events[i],
                    b.nonradiative_events[i],
                ),
                (
                    whole.crossing_events[i],
                    a.crossing_events[i],
                    b.crossing_events[i],
                ),
            ] {
                close(all, left + right, 1e-11);
            }
        }
        balanced(&whole, 1e-17);
    }

    #[test]
    fn g0_triplet_afterglow_remains_after_fast_singlet_exhaustion() {
        let m = model([5e5, 0.1], [0.0; 2], [5e5, 0.0]);
        let step = plan(m, 10.0).advance([1.0, 0.0], [0.0; 2]).unwrap();
        let expected_triplet = (5e5 / (1e6 - 0.1)) * (-1.0_f64).exp();
        close(step.populations[1], expected_triplet, 3e-11);
        close(step.radiative_events[0], 0.5, 3e-11);
        close(step.radiative_events[1], 0.5 - expected_triplet, 3e-11);
        balanced(&step, 4e-19);
    }

    #[test]
    fn g0_stable_store_zero_interval_and_cancelled_candidate() {
        let m = model([0.0; 2], [0.0; 2], [0.0; 2]);
        let stable = plan(m, 2.0).advance([1.0, 2.0], [3.0, 4.0]).unwrap();
        assert_eq!(stable.populations, [7.0, 10.0]);
        assert_eq!(stable.radiative_events, [0.0; 2]);
        balanced(&stable, 6e-18);
        let zero = plan(m, 0.0).advance([1.0, 2.0], [3.0, 4.0]).unwrap();
        assert_eq!(zero.populations, [1.0, 2.0]);
        assert_eq!(zero.injected_excitations, [0.0; 2]);
        assert!(matches!(
            m.prepare(1.0, DecayLimits::default(), || true),
            Err(SingletTripletError::Evolution(DecayError::Cancelled))
        ));
        assert!(matches!(
            plan(m, 1.0).advance_with_cancel([1.0; 2], [0.0; 2], || true),
            Err(SingletTripletError::Evolution(DecayError::Cancelled))
        ));
    }

    #[test]
    fn g0_fast_crossing_cannot_erase_or_hide_ground_state_channels() {
        for (radiative, dark) in [([1.0; 2], [0.0; 2]), ([0.0; 2], [1.0; 2])] {
            assert!(matches!(
                SingletTripletDecay::try_new([4.0, 3.0], radiative, dark, [1e20; 2]),
                Err(SingletTripletError::Unrepresentable(
                    "ground-state sink erased by crossing fraction"
                ))
            ));
        }
        // A nonzero sink with rate-relative roundoff can still resolve a short
        // physical interval. Admission is based on the actual candidate counts.
        let fast = model([1.0; 2], [0.0; 2], [1e9; 2]);
        let short = plan(fast, 1e-12).advance([1.0, 0.0], [0.0; 2]).unwrap();
        assert!(short.radiative_events.iter().sum::<f64>() > 0.0);
        // A candidate retaining all initial population while radiating one
        // excitation is invalid regardless of its unreported internal traffic.
        assert!(matches!(
            admit_population_balance(
                &[1.0, 0.0],
                &[0.0; 2],
                &[1.0, 0.0],
                &[0.25, 0.75],
                &[0.0; 2],
                2e-12,
            ),
            Err(SingletTripletError::PopulationBalance { .. })
        ));
    }

    #[test]
    fn g0_positive_event_and_energy_underflow_refuse() {
        // The radiative event count is ~1e-400, but multiplying by its 1e300 J
        // state energy would give resolved ~1e-100 J. Returning zero events and
        // zero light would therefore hide a physically representable channel.
        let tiny_branch =
            SingletTripletDecay::try_new([1e300; 2], [1e-300, 0.0], [1.0, 0.0], [0.0; 2]).unwrap();
        assert!(matches!(
            plan(tiny_branch, 1e-100).advance([1.0, 0.0], [0.0; 2]),
            Err(SingletTripletError::Unrepresentable("radiative events"))
        ));
        let tiny_energy =
            SingletTripletDecay::try_new([1e-300; 2], [0.0; 2], [0.0; 2], [0.0; 2]).unwrap();
        assert!(matches!(
            plan(tiny_energy, 1.0).advance([1e-100, 0.0], [0.0; 2]),
            Err(SingletTripletError::Unrepresentable("stored energy"))
        ));
    }

    #[test]
    fn g0_invalid_rates_energies_and_unrepresentable_results_refuse() {
        for bad in [-1.0, f64::INFINITY, f64::NAN] {
            assert!(
                SingletTripletDecay::try_new([4.0, 3.0], [bad, 1.0], [0.0; 2], [0.0; 2]).is_err()
            );
        }
        assert!(SingletTripletDecay::try_new([2.0, 3.0], [1.0; 2], [0.0; 2], [0.0; 2]).is_err());
        assert!(
            SingletTripletDecay::try_new([4.0, 3.0], [f64::MAX; 2], [f64::MAX; 2], [0.0; 2])
                .is_err()
        );
        assert!(
            SingletTripletDecay::try_new(
                [4.0, 3.0],
                [f64::from_bits(1), 0.0],
                [1e300, 0.0],
                [0.0; 2]
            )
            .is_err()
        );
        let p = plan(model([0.0; 2], [0.0; 2], [0.0; 2]), 2.0);
        assert!(p.advance([0.0; 2], [f64::MAX, 0.0]).is_err());
        assert!(p.advance([f64::NAN, 0.0], [0.0; 2]).is_err());
        let expensive =
            SingletTripletDecay::try_new([f64::MAX; 2], [0.0; 2], [0.0; 2], [0.0; 2]).unwrap();
        assert!(plan(expensive, 1.0).advance([2.0, 0.0], [0.0; 2]).is_err());
    }
}
