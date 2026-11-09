//! Conserved ground/singlet/triplet populations under a prescribed pump field.
//!
//! Each pump hazard acts on the CURRENT ground-state population, so
//! `G' = -(Ws + Wt) G + (ksr + ksnr) S + (ktr + ktnr) T`,
//! `S' = Ws G - (ksr + ksnr + kst) S + kts T`, and
//! `T' = Wt G + kst S - (ktr + ktnr + kts) T`.
//! Their sum is constant. Pumping cannot create extra emitters; shelving in a
//! long-lived triplet reduces the ground population available for absorption.
//! The existing positive reaction exponential integrates both populations and
//! every channel's actual events, including repeated excitation/decay cycles.
//!
//! Pump rates are per-ground-emitter hazards [1/s], NOT already absorbed watts
//! or constant excitation injection. Optical callers can supply cross section
//! times local photon flux; other excitation mechanisms must author their own
//! hazards and energy per event. Only actual pump events debit that energy.
//! The field and coefficients stay fixed over this interval: this is not a
//! solution for attenuation, electrical carrier transport or pump depletion.
//! Do not feed a fixed Beer absorption result into this model as if the ground
//! population could not change. Returned radiative joules still precede spectral
//! Stokes partition, escape and reabsorption.
//!
//! This incoherent rate model omits stimulated emission, excited-state pump
//! absorption, bleaching, trap distributions and pair annihilation. An optical
//! use requires the regime in which rapid pump relaxation justifies those
//! omissions; it is not a coherently driven resonant two-level atom. The
//! ground/singlet/triplet saturation mechanism is used, for example, by Fujiwara
//! et al., J. Cell Biol. 222, e202110160 (2023), Methods, "Triplet bottleneck
//! saturation", <https://doi.org/10.1083/jcb.202110160>. No coefficients are
//! inferred from that example or from a material name.

use super::{
    CompensatedSum, SingletTripletDecay, SingletTripletError, admit_population_balance, finite,
    nonnegative, positive_product,
};
use fs_la::decay::DecayLimits;
use fs_la::decay::reaction::{ReactionEdge, ReactionPlan};

/// Reusable fixed-field evolution of a conserved `[ground, singlet, triplet]`
/// inventory. Populations are counts, not number densities.
#[derive(Clone, Debug)]
pub struct SaturableSingletTripletPlan {
    model: SingletTripletDecay,
    pump_rates_s: [f64; 2],
    pump_energies_j: [f64; 2],
    pump_fractions: [f64; 2],
    duration_s: f64,
    relative_tolerance: f64,
    evolution: ReactionPlan,
}

/// Candidate populations and integrated physical channels over one interval.
/// No input is mutated, no population is clipped and no residual becomes heat.
#[derive(Clone, Debug, PartialEq)]
pub struct SaturableSingletTripletStep {
    /// Endpoint counts ordered `[ground, singlet, triplet]`.
    pub populations: [f64; 3],
    /// Actual ground-to-`[singlet, triplet]` transitions, including repumping.
    pub pump_events: [f64; 2],
    /// Original external pump energy actually consumed by those transitions.
    pub absorbed_pump_energy_j: [f64; 2],
    /// Pump energy above the admitted state's energy, deposited as local heat.
    pub pump_relaxation_heat_j: [f64; 2],
    /// Integrated radiative ground-state transitions from `[singlet, triplet]`.
    pub radiative_events: [f64; 2],
    /// Integrated nonradiative ground-state transitions from each excited state.
    pub nonradiative_events: [f64; 2],
    /// Integrated `[singlet -> triplet, triplet -> singlet]` transitions.
    pub crossing_events: [f64; 2],
    /// Endpoint excitation energies `[singlet, triplet]` relative to the ground.
    pub stored_energy_j: [f64; 2],
    /// State energy released radiatively before photon-spectrum conversion.
    pub radiative_decay_energy_j: [f64; 2],
    /// Local heat from each state's nonradiative ground-state transitions.
    pub nonradiative_heat_j: [f64; 2],
    /// Downhill singlet-to-triplet energy deposited in the thermal bath.
    pub crossing_heat_j: f64,
    /// Uphill triplet-to-singlet energy withdrawn from the thermal bath.
    pub reverse_crossing_heat_absorbed_j: f64,
    /// Signed final-minus-initial emitter count, without correction.
    pub population_residual: f64,
    /// `old excitation + absorbed pump - new excitation - radiation - net heat`.
    /// Ground inventory and reversible crossing traffic do not inflate its scale.
    pub residual_j: f64,
}

impl SaturableSingletTripletPlan {
    /// Prepare per-ground pump hazards `[G -> S, G -> T]` and the energy each
    /// actual pump event consumes [J]. Every active pump must supply at least
    /// its target state's energy. An inactive pump may use zero energy.
    /// Preparation/application retain the reaction owner's work and accuracy
    /// limits and bounded cancellation. Coefficients and pump field are fixed.
    pub fn prepare(
        model: SingletTripletDecay,
        pump_rates_s: [f64; 2],
        pump_energies_j: [f64; 2],
        duration_s: f64,
        limits: DecayLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SingletTripletError> {
        Self::prepare_impl(
            model,
            pump_rates_s,
            pump_energies_j,
            duration_s,
            limits,
            false,
            cancelled,
        )
    }

    /// Also retain the native population-time integral for each real state.
    /// This includes a stable ground state and does not divide removal counts
    /// by an artificial or vanishing rate. The extra preparation work remains
    /// subject to the original kinetic limits.
    pub fn prepare_with_population_integrals(
        model: SingletTripletDecay,
        pump_rates_s: [f64; 2],
        pump_energies_j: [f64; 2],
        duration_s: f64,
        limits: DecayLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SingletTripletError> {
        Self::prepare_impl(
            model,
            pump_rates_s,
            pump_energies_j,
            duration_s,
            limits,
            true,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_impl(
        model: SingletTripletDecay,
        pump_rates_s: [f64; 2],
        pump_energies_j: [f64; 2],
        duration_s: f64,
        limits: DecayLimits,
        with_population_integrals: bool,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SingletTripletError> {
        if !pump_rates_s.into_iter().all(nonnegative)
            || !pump_energies_j.into_iter().all(nonnegative)
            || (0..2).any(|i| pump_rates_s[i] > 0.0 && pump_energies_j[i] < model.energies_j[i])
        {
            return Err(SingletTripletError::Invalid(
                "ground-state pump hazards or event energies",
            ));
        }
        let ground_rate = finite(pump_rates_s[0] + pump_rates_s[1], "total ground pump rate")?;
        let mut pump_fractions = [0.0; 2];
        let mut edges = Vec::with_capacity(8);
        for i in 0..2 {
            if pump_rates_s[i] > 0.0 {
                pump_fractions[i] = pump_rates_s[i] / ground_rate;
                if pump_fractions[i] == 0.0 {
                    return Err(SingletTripletError::Unrepresentable(
                        "positive ground-state pump fraction",
                    ));
                }
                edges.push(ReactionEdge {
                    parent: 0,
                    daughter: i + 1,
                    yield_per_removal: pump_fractions[i],
                });
            }
            // Returning the emitter to ground is distinct from disposing of
            // its excitation energy. Preserve both authored physical channels.
            for fraction in &model.fractions[i][..2] {
                if *fraction > 0.0 {
                    edges.push(ReactionEdge {
                        parent: i + 1,
                        daughter: 0,
                        yield_per_removal: *fraction,
                    });
                }
            }
            if model.fractions[i][2] > 0.0 {
                edges.push(ReactionEdge {
                    parent: i + 1,
                    daughter: 2 - i,
                    yield_per_removal: model.fractions[i][2],
                });
            }
        }
        let rates = [ground_rate, model.total_rates_s[0], model.total_rates_s[1]];
        let evolution = if with_population_integrals {
            ReactionPlan::prepare_with_population_integrals(
                &rates, &edges, duration_s, limits, cancelled,
            )?
        } else {
            ReactionPlan::prepare(&rates, &edges, duration_s, limits, cancelled)?
        };
        Ok(Self {
            model,
            pump_rates_s,
            pump_energies_j,
            pump_fractions,
            duration_s,
            relative_tolerance: limits.relative_tolerance,
            evolution,
        })
    }

    /// Original excited-state energies and competing decay/crossing rates.
    pub const fn model(&self) -> &SingletTripletDecay {
        &self.model
    }
    /// Per-ground-emitter `[G -> S, G -> T]` pump hazards, inverse seconds.
    pub const fn pump_rates_s(&self) -> [f64; 2] {
        self.pump_rates_s
    }
    /// Energy consumed by each actual `[G -> S, G -> T]` pump event, joules.
    pub const fn pump_energies_j(&self) -> [f64; 2] {
        self.pump_energies_j
    }
    /// Prepared interval length, seconds.
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Population-time integrals `[ground, singlet, triplet]` in count-seconds.
    /// Requires `prepare_with_population_integrals`; inputs are unchanged.
    pub fn integrate_populations_with_cancel(
        &self,
        populations: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<[f64; 3], SingletTripletError> {
        let result = self
            .evolution
            .integrate_populations_with_cancel(&populations, cancelled)?;
        Ok([result[0], result[1], result[2]])
    }

    /// Advance finite nonnegative `[ground, singlet, triplet]` counts.
    pub fn advance(
        &self,
        populations: [f64; 3],
    ) -> Result<SaturableSingletTripletStep, SingletTripletError> {
        self.advance_with_cancel(populations, || false)
    }

    /// Compute a complete independent candidate with bounded row checkpoints.
    /// The caller owns thermal-bath admissibility and acceptance of this step.
    pub fn advance_with_cancel(
        &self,
        populations: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SaturableSingletTripletStep, SingletTripletError> {
        let result = self.evolution.apply_with_cancel(&populations, cancelled)?;
        let mut step = SaturableSingletTripletStep {
            populations: [
                result.populations[0],
                result.populations[1],
                result.populations[2],
            ],
            pump_events: [0.0; 2],
            absorbed_pump_energy_j: [0.0; 2],
            pump_relaxation_heat_j: [0.0; 2],
            radiative_events: [0.0; 2],
            nonradiative_events: [0.0; 2],
            crossing_events: [0.0; 2],
            stored_energy_j: [0.0; 2],
            radiative_decay_energy_j: [0.0; 2],
            nonradiative_heat_j: [0.0; 2],
            crossing_heat_j: 0.0,
            reverse_crossing_heat_absorbed_j: 0.0,
            population_residual: 0.0,
            residual_j: 0.0,
        };
        for i in 0..2 {
            let energy = self.model.energies_j[i];
            let fractions = self.model.fractions[i];
            step.pump_events[i] =
                positive_product(result.removals[0], self.pump_fractions[i], "pump events")?;
            step.absorbed_pump_energy_j[i] = positive_product(
                step.pump_events[i],
                self.pump_energies_j[i],
                "absorbed pump energy",
            )?;
            let excess = if self.pump_rates_s[i] > 0.0 {
                self.pump_energies_j[i] - energy
            } else {
                0.0
            };
            step.pump_relaxation_heat_j[i] =
                positive_product(step.pump_events[i], excess, "pump relaxation heat")?;
            step.radiative_events[i] =
                positive_product(result.removals[i + 1], fractions[0], "radiative events")?;
            step.nonradiative_events[i] =
                positive_product(result.removals[i + 1], fractions[1], "nonradiative events")?;
            step.crossing_events[i] =
                positive_product(result.removals[i + 1], fractions[2], "crossing events")?;
            step.stored_energy_j[i] =
                positive_product(step.populations[i + 1], energy, "stored excitation energy")?;
            step.radiative_decay_energy_j[i] =
                positive_product(step.radiative_events[i], energy, "radiative energy")?;
            step.nonradiative_heat_j[i] =
                positive_product(step.nonradiative_events[i], energy, "nonradiative heat")?;
        }
        admit_population_balance(
            &[populations[1], populations[2]],
            &step.pump_events,
            &[step.populations[1], step.populations[2]],
            &step.radiative_events,
            &step.nonradiative_events,
            self.relative_tolerance,
        )?;
        let mut residual = CompensatedSum::default();
        let mut population_scale = CompensatedSum::default();
        for i in 0..3 {
            residual.add(step.populations[i])?;
            residual.add(-populations[i])?;
            population_scale.add(step.populations[i])?;
            population_scale.add(populations[i])?;
        }
        step.population_residual = total(&residual, "finite emitter population residual")?;
        let allowance = self.relative_tolerance * total(&population_scale, "emitter inventory")?;
        if step.population_residual.abs() > allowance {
            return Err(SingletTripletError::PopulationBalance {
                residual: step.population_residual,
                allowance,
            });
        }
        let gap = self.model.energies_j[0] - self.model.energies_j[1];
        step.crossing_heat_j = positive_product(step.crossing_events[0], gap, "crossing heat")?;
        step.reverse_crossing_heat_absorbed_j =
            positive_product(step.crossing_events[1], gap, "reverse crossing heat")?;
        step.residual_j = self.energy_residual(populations, &step)?;
        Ok(step)
    }

    fn energy_residual(
        &self,
        initial: [f64; 3],
        step: &SaturableSingletTripletStep,
    ) -> Result<f64, SingletTripletError> {
        let mut balance = CompensatedSum::default();
        let mut scale = CompensatedSum::default();
        let mut heat = CompensatedSum::default();
        for i in 0..2 {
            for term in [
                positive_product(
                    initial[i + 1],
                    self.model.energies_j[i],
                    "initial excitation",
                )?,
                step.absorbed_pump_energy_j[i],
                -step.stored_energy_j[i],
                -step.radiative_decay_energy_j[i],
            ] {
                balance.add(term)?;
                scale.add(term.abs())?;
            }
            heat.add(step.pump_relaxation_heat_j[i])?;
            heat.add(step.nonradiative_heat_j[i])?;
        }
        heat.add(step.crossing_heat_j)?;
        heat.add(-step.reverse_crossing_heat_absorbed_j)?;
        // Only the NET bath exchange enters the allowance; nearly balanced
        // crossing cycles must not hide a physical energy-accounting error.
        let net_heat = total(&heat, "net finite emitter thermal exchange")?;
        balance.add(-net_heat)?;
        scale.add(net_heat.abs())?;
        let residual_j = total(&balance, "finite emitter energy residual")?;
        let allowance_j = self.relative_tolerance * total(&scale, "finite emitter energy scale")?;
        if residual_j.abs() > allowance_j {
            return Err(SingletTripletError::EnergyBalance {
                residual_j,
                allowance_j,
            });
        }
        Ok(residual_j)
    }
}

fn total(sum: &CompensatedSum, name: &'static str) -> Result<f64, SingletTripletError> {
    finite(sum.value + sum.correction, name)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 2e-10 * expected.abs().max(1e-300),
            "actual={actual:e}, expected={expected:e}"
        );
    }

    fn prepare(
        model: SingletTripletDecay,
        pump: [f64; 2],
        duration: f64,
    ) -> SaturableSingletTripletPlan {
        SaturableSingletTripletPlan::prepare(
            model,
            pump,
            [5e-19, 4e-19],
            duration,
            DecayLimits::default(),
            || false,
        )
        .unwrap()
    }

    #[test]
    fn g1_ground_depletion_matches_two_state_saturation_and_actual_pump_power() {
        let model =
            SingletTripletDecay::try_new([4e-19, 3e-19], [6.0, 0.0], [2.0, 0.0], [0.0; 2]).unwrap();
        let n = 2e12;
        let h: f64 = 0.7;
        for w in [0.001, 8.0, 800.0] {
            let step = prepare(model, [w, 0.0], h).advance([n, 0.0, 0.0]).unwrap();
            let steady = n * w / (w + 8.0);
            let singlet = steady * (1.0 - (-(w + 8.0) * h).exp());
            let residence = steady * (h - (1.0 - (-(w + 8.0) * h).exp()) / (w + 8.0));
            close(step.populations[0], n - singlet);
            close(step.populations[1], singlet);
            assert_eq!(step.populations[2], 0.0);
            close(step.radiative_events[0], 6.0 * residence);
            close(step.nonradiative_events[0], 2.0 * residence);
            close(step.pump_events[0], singlet + 8.0 * residence);
            close(step.absorbed_pump_energy_j[0], step.pump_events[0] * 5e-19);
            close(step.pump_relaxation_heat_j[0], step.pump_events[0] * 1e-19);
            assert!(step.pump_events[0] < n * w * h);
            assert!(step.populations[1] <= n);
            assert!(step.residual_j.abs() < 1e-17);
        }
    }

    #[test]
    fn g1_triplet_shelving_depletes_absorbers_without_inventing_emitters() {
        let model =
            SingletTripletDecay::try_new([4e-19, 3e-19], [0.0; 2], [0.0; 2], [3.0, 0.0]).unwrap();
        let n = 100.0;
        let h: f64 = 1.2;
        let survived = (-3.0 * h).exp();
        let step = prepare(model, [3.0, 0.0], h)
            .advance([n, 0.0, 0.0])
            .unwrap();
        close(step.populations[0], n * survived);
        close(step.populations[1], n * 3.0 * h * survived);
        close(step.populations[2], n * (1.0 - (1.0 + 3.0 * h) * survived));
        close(step.pump_events[0], n * (1.0 - survived));
        close(step.crossing_events[0], step.populations[2]);
        close(step.crossing_heat_j, step.populations[2] * 1e-19);
        let shelved = prepare(model, [1e5, 0.0], 1.0)
            .advance([0.0, 0.0, n])
            .unwrap();
        assert_eq!(shelved.populations, [0.0, 0.0, n]);
        assert_eq!(shelved.pump_events, [0.0; 2]);
        assert_eq!(shelved.absorbed_pump_energy_j, [0.0; 2]);
    }

    #[test]
    fn g1_saturated_store_takes_no_more_than_one_pump_quantum_per_emitter() {
        let model =
            SingletTripletDecay::try_new([4e-19, 3e-19], [0.0; 2], [0.0; 2], [0.0; 2]).unwrap();
        let step = prepare(model, [100.0, 0.0], 1.0)
            .advance([1e9, 0.0, 0.0])
            .unwrap();
        close(step.populations[1], 1e9);
        close(step.pump_events[0], 1e9);
        close(step.absorbed_pump_energy_j[0], 5e-10);
        assert_eq!(step.radiative_events, [0.0; 2]);
        let empty = prepare(model, [100.0; 2], 1.0).advance([0.0; 3]).unwrap();
        assert_eq!(empty.populations, [0.0; 3]);
        assert_eq!(empty.absorbed_pump_energy_j, [0.0; 2]);
    }

    #[test]
    fn g3_both_pumps_crossings_and_every_integrated_channel_compose() {
        let model =
            SingletTripletDecay::try_new([4e-19, 3e-19], [6.0, 0.2], [2.0, 0.3], [5.0, 0.8])
                .unwrap();
        let initial = [50.0, 10.0, 40.0];
        let whole = prepare(model, [7.0, 2.0], 1.0).advance(initial).unwrap();
        let first = prepare(model, [7.0, 2.0], 0.3).advance(initial).unwrap();
        let second = prepare(model, [7.0, 2.0], 0.7)
            .advance(first.populations)
            .unwrap();
        for i in 0..3 {
            close(whole.populations[i], second.populations[i]);
        }
        for (total, a, b) in [
            (whole.pump_events, first.pump_events, second.pump_events),
            (
                whole.absorbed_pump_energy_j,
                first.absorbed_pump_energy_j,
                second.absorbed_pump_energy_j,
            ),
            (
                whole.pump_relaxation_heat_j,
                first.pump_relaxation_heat_j,
                second.pump_relaxation_heat_j,
            ),
            (
                whole.radiative_events,
                first.radiative_events,
                second.radiative_events,
            ),
            (
                whole.nonradiative_events,
                first.nonradiative_events,
                second.nonradiative_events,
            ),
            (
                whole.crossing_events,
                first.crossing_events,
                second.crossing_events,
            ),
            (
                whole.radiative_decay_energy_j,
                first.radiative_decay_energy_j,
                second.radiative_decay_energy_j,
            ),
            (
                whole.nonradiative_heat_j,
                first.nonradiative_heat_j,
                second.nonradiative_heat_j,
            ),
        ] {
            for i in 0..2 {
                close(total[i], a[i] + b[i]);
            }
        }
        close(
            whole.crossing_heat_j,
            first.crossing_heat_j + second.crossing_heat_j,
        );
        close(
            whole.reverse_crossing_heat_absorbed_j,
            first.reverse_crossing_heat_absorbed_j + second.reverse_crossing_heat_absorbed_j,
        );
        close(whole.populations.iter().sum(), 100.0);
        assert!(whole.population_residual.abs() < 1e-9);
        assert!(whole.residual_j.abs() < 1e-25);
    }

    #[test]
    fn g3_pump_off_rejoins_original_afterglow_and_restores_ground() {
        let model =
            SingletTripletDecay::try_new([4e-19, 3e-19], [5.0, 0.1], [1.0, 0.4], [2.0, 0.3])
                .unwrap();
        let initial = [20.0, 30.0, 50.0];
        let dark = prepare(model, [0.0; 2], 1.0).advance(initial).unwrap();
        let original = model
            .prepare(1.0, DecayLimits::default(), || false)
            .unwrap()
            .advance([30.0, 50.0], [0.0; 2])
            .unwrap();
        for i in 0..2 {
            close(dark.populations[i + 1], original.populations[i]);
            close(dark.radiative_events[i], original.radiative_events[i]);
            close(dark.nonradiative_events[i], original.nonradiative_events[i]);
        }
        close(
            dark.populations[0],
            100.0 - original.populations.iter().sum::<f64>(),
        );
        assert_eq!(dark.absorbed_pump_energy_j, [0.0; 2]);
    }

    #[test]
    fn g3_nanosecond_rates_and_macroscopic_counts_preserve_si_scaling() {
        // Visible-photon state energies, a ~1 ns singlet, a ~20 us triplet,
        // and a microsecond exposure containing many singlet lifetimes.
        let model =
            SingletTripletDecay::try_new([4e-19, 3e-19], [6e8, 2e4], [2e8, 3e4], [5e7, 8e3])
                .unwrap();
        let plan = prepare(model, [7e7, 2e6], 1e-6);
        let unit = plan.advance([0.9, 0.05, 0.05]).unwrap();
        let bulk = plan.advance([9e19, 5e18, 5e18]).unwrap();
        for i in 0..3 {
            close(bulk.populations[i] / 1e20, unit.populations[i]);
        }
        for (large, small) in [
            (bulk.pump_events, unit.pump_events),
            (bulk.stored_energy_j, unit.stored_energy_j),
            (bulk.absorbed_pump_energy_j, unit.absorbed_pump_energy_j),
            (bulk.radiative_decay_energy_j, unit.radiative_decay_energy_j),
            (bulk.nonradiative_heat_j, unit.nonradiative_heat_j),
        ] {
            for i in 0..2 {
                close(large[i] / 1e20, small[i]);
            }
        }
        close(bulk.populations.iter().sum(), 1e20);
        assert!(bulk.populations.iter().all(|&count| count >= 0.0));
        assert!(bulk.population_residual.abs() < 1e-10 * 1e20);
        assert!(bulk.residual_j.abs() < 1e-9);
    }

    #[test]
    fn g0_energy_admission_zero_time_cancellation_and_budgets_refuse_atomically() {
        let model =
            SingletTripletDecay::try_new([4e-19, 3e-19], [1.0; 2], [0.0; 2], [0.0; 2]).unwrap();
        let plan = prepare(model, [5.0, 0.0], 0.0);
        let state = [10.0, 2.0, 3.0];
        let step = plan.advance(state).unwrap();
        assert_eq!(step.populations, state);
        assert_eq!(step.pump_events, [0.0; 2]);
        assert_eq!(step.radiative_events, [0.0; 2]);
        assert!(plan.advance([-1.0, 0.0, 0.0]).is_err());
        assert!(plan.advance([f64::NAN, 0.0, 0.0]).is_err());
        assert!(
            SaturableSingletTripletPlan::prepare(
                model,
                [1.0, 0.0],
                [3e-19, 0.0],
                1.0,
                DecayLimits::default(),
                || false
            )
            .is_err()
        );
        assert!(
            SaturableSingletTripletPlan::prepare(
                model,
                [1.0, 0.0],
                [5e-19, 0.0],
                1.0,
                DecayLimits::default(),
                || true
            )
            .is_err()
        );
        assert!(
            SaturableSingletTripletPlan::prepare(
                model,
                [1.0, 0.0],
                [5e-19, 0.0],
                1.0,
                DecayLimits {
                    maximum_nodes: 2,
                    ..DecayLimits::default()
                },
                || false
            )
            .is_err()
        );
        assert!(plan.advance_with_cancel(state, || true).is_err());
        assert_eq!(state, [10.0, 2.0, 3.0]);
    }
}
