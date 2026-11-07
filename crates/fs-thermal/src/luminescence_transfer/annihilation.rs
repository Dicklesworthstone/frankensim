//! Exact, spatially uniform triplet-triplet annihilation over a fixed interval.
//!
//! The authored pair-event coefficient `k` [m^3/s] means `R = k N_t^2 / V`
//! events/s, so `dN_t/dt = -2 R`. A coefficient defined by
//! `dn_t/dt = -k_loss n_t^2` must therefore be supplied as `k_loss / 2`.
//! Each event consumes TWO triplets and produces one singlet with probability
//! `phi`, otherwise only heat. The singlet subsequently evolves through the
//! ordinary linear owner; this step emits no photon. A positive singlet branch
//! requires `E_s <= 2 E_t` per event, irrespective of its average yield.
//!
//! This is a continuum mean-population model with a fixed active volume and
//! authored, constant encounter coefficient. It does not infer a diffusion
//! rate, spin statistical yield, ground-state saturation or spatial mixing.
//! The rational solution is exact for this isolated channel; combining it with
//! pumping, crossing, diffusion or decay requires a controlled coupling step.
//! Work is bounded and allocation-free. Arithmetic is deterministic IEEE;
//! binary scaling preserves channels across extreme intermediate exposures.

use super::{CompensatedSum, SingletTripletError, finite, nonnegative, positive_product};

/// Authored pair-event coefficient, active volume and singlet yield per pair.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TripletAnnihilation {
    pair_rate_m3_s: f64,
    active_volume_m3: f64,
    singlet_yield: f64,
}

/// Candidate populations and integrated annihilation channels on one interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TripletAnnihilationStep {
    /// Endpoint excitation counts, ordered `[singlet, triplet]`.
    pub populations: [f64; 2],
    /// Integrated pair events; each event removes two triplets.
    pub pair_events: f64,
    /// Singlets produced, equal to the authored yield times pair events.
    pub singlet_production: f64,
    /// Local heat from unsuccessful pairs and the productive branch's energy gap.
    pub nonradiative_heat_j: f64,
    /// Represented `old_store - new_store - heat`; never folded into a channel.
    pub residual_j: f64,
}

impl TripletAnnihilation {
    /// Admit finite `k >= 0`, `V > 0`, and singlet yield in `[0, 1]`.
    /// Counts passed to `advance` refer to this entire active volume, not density.
    pub fn try_new(
        pair_rate_m3_s: f64,
        active_volume_m3: f64,
        singlet_yield: f64,
    ) -> Result<Self, SingletTripletError> {
        if !nonnegative(pair_rate_m3_s)
            || !active_volume_m3.is_finite()
            || active_volume_m3 <= 0.0
            || !singlet_yield.is_finite()
            || !(0.0..=1.0).contains(&singlet_yield)
        {
            return Err(SingletTripletError::Invalid(
                "annihilation coefficient, volume or yield",
            ));
        }
        Ok(Self {
            pair_rate_m3_s,
            active_volume_m3,
            singlet_yield,
        })
    }

    /// Authored pair-event coefficient, cubic metres per second.
    pub const fn pair_rate_m3_s(&self) -> f64 {
        self.pair_rate_m3_s
    }
    /// Active volume containing the supplied excitation counts, cubic metres.
    pub const fn active_volume_m3(&self) -> f64 {
        self.active_volume_m3
    }
    /// Probability of producing one singlet per annihilated triplet pair.
    pub const fn singlet_yield(&self) -> f64 {
        self.singlet_yield
    }

    /// Advance `N_t(h) = N_t(0)/(1 + 2 k N_t(0) h/V)` without mutation.
    ///
    /// Energies are positive finite joules per excitation with `E_s >= E_t`,
    /// as in the linear singlet/triplet owner. Counts and seconds are finite
    /// and nonnegative. Positive channels that underflow are refused, including
    /// a positive pair count even if a later energy factor could rescale it.
    /// The signed residual exposes rounding of endpoint stores to the caller.
    pub fn advance(
        &self,
        energies: [f64; 2],
        populations: [f64; 2],
        duration_s: f64,
    ) -> Result<TripletAnnihilationStep, SingletTripletError> {
        let [es, et] = energies;
        if !energies.iter().all(|&x| x.is_finite() && x > 0.0)
            || es < et
            || !populations.into_iter().all(nonnegative)
            || !nonnegative(duration_s)
        {
            return Err(SingletTripletError::Invalid(
                "annihilation energies, populations or duration",
            ));
        }
        // Subtracting first avoids overflowing 2*Et, and admits the branch
        // eventwise: a tiny yield cannot make an uphill individual event valid.
        if self.singlet_yield > 0.0 && es - et > et {
            return Err(SingletTripletError::Invalid(
                "annihilation singlet exceeds two triplet energies",
            ));
        }
        let initial = [
            positive_product(populations[0], es, "initial singlet energy")?,
            positive_product(populations[1], et, "initial triplet energy")?,
        ];
        if self.pair_rate_m3_s == 0.0 || populations[1] == 0.0 || duration_s == 0.0 {
            return Ok(TripletAnnihilationStep {
                populations,
                pair_events: 0.0,
                singlet_production: 0.0,
                nonradiative_heat_j: 0.0,
                residual_j: 0.0,
            });
        }
        let (km, ke) = parts(self.pair_rate_m3_s);
        let (nm, ne) = parts(populations[1]);
        let (hm, he) = parts(duration_s);
        let (vm, ve) = parts(self.active_volume_m3);
        // z = zm * 2^ze = 2*k*Nt*h/V. Its mantissa stays in [1/2, 8),
        // even when z itself overflows or underflows binary64.
        let zm = km * nm * hm / vm;
        let ze = ke + ne + he - ve + 1;
        let z = scaled(zm, ze);
        let (triplets, pair_events) = if z <= 1.0 {
            let denominator = 1.0 + z;
            (
                populations[1] / denominator,
                scaled(nm * zm / (2.0 * denominator), ne + ze),
            )
        } else {
            let denominator = 1.0 + scaled(1.0 / zm, -ze);
            (
                scaled(nm / zm / denominator, ne - ze),
                populations[1] / (2.0 * denominator),
            )
        };
        if triplets == 0.0 || pair_events == 0.0 {
            return Err(SingletTripletError::Unrepresentable(
                "positive triplet endpoint or pair count",
            ));
        }
        let singlet_production = positive_product(
            pair_events,
            self.singlet_yield,
            "annihilation singlet production",
        )?;
        let new_populations = [
            finite(
                populations[0] + singlet_production,
                "annihilation singlet endpoint",
            )?,
            triplets,
        ];
        let dark_pairs = positive_product(
            pair_events,
            1.0 - self.singlet_yield,
            "nonproductive annihilation pairs",
        )?;
        let dark_energy = positive_product(dark_pairs, et, "nonproductive annihilation energy")?;
        let productive_gap = if self.singlet_yield > 0.0 {
            et - (es - et)
        } else {
            0.0
        };
        let mut heat = CompensatedSum::default();
        heat.add(dark_energy)?;
        heat.add(dark_energy)?;
        heat.add(positive_product(
            singlet_production,
            productive_gap,
            "productive annihilation heat",
        )?)?;
        let nonradiative_heat_j = finite(heat.value + heat.correction, "annihilation heat")?;
        let mut balance = CompensatedSum::default();
        for i in 0..2 {
            balance.add(initial[i])?;
            balance.add(-positive_product(
                new_populations[i],
                energies[i],
                "annihilation stored energy",
            )?)?;
        }
        balance.add(-nonradiative_heat_j)?;
        Ok(TripletAnnihilationStep {
            populations: new_populations,
            pair_events,
            singlet_production,
            nonradiative_heat_j,
            residual_j: finite(
                balance.value + balance.correction,
                "annihilation energy residual",
            )?,
        })
    }
}

// Exact binary decomposition of a finite positive input, including subnormals.
fn parts(value: f64) -> (f64, i32) {
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    if exponent == 0 {
        let (mantissa, exponent) = parts(value * 4_503_599_627_370_496.0);
        return (mantissa, exponent - 52);
    }
    (
        f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000),
        exponent - 1023,
    )
}

// Only bounded, positive mantissas reach this helper. Powers of two are exact;
// one final multiplication rounds a subnormal result without flushing it.
fn scaled(mut mantissa: f64, mut exponent: i32) -> f64 {
    while mantissa >= 2.0 {
        mantissa *= 0.5;
        exponent += 1;
    }
    while mantissa < 1.0 {
        mantissa *= 2.0;
        exponent -= 1;
    }
    if exponent > 1023 {
        f64::INFINITY
    } else if exponent < -1075 {
        0.0
    } else if exponent >= -1022 {
        mantissa * f64::from_bits(((exponent + 1023) as u64) << 52)
    } else {
        (mantissa * f64::from_bits(((exponent + 1076) as u64) << 52))
            * (1.0 / 9_007_199_254_740_992.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 8e-15 * expected.abs().max(f64::MIN_POSITIVE),
            "actual={actual:e}, expected={expected:e}"
        );
    }

    #[test]
    fn pure_triplets_follow_reciprocal_time_and_the_authored_pair_convention() {
        let model = TripletAnnihilation::try_new(0.25, 2.0, 0.5).unwrap();
        for duration in [0.125, 0.5, 2.0] {
            let step = model.advance([3.0, 2.0], [0.0, 8.0], duration).unwrap();
            close(
                1.0 / step.populations[1] - 1.0 / 8.0,
                2.0 * 0.25 * duration / 2.0,
            );
            close(step.pair_events, 0.5 * (8.0 - step.populations[1]));
            close(step.singlet_production, 0.5 * step.pair_events);
            close(step.populations[0], step.singlet_production);
            close(step.nonradiative_heat_j, 2.5 * step.pair_events);
            assert!(step.residual_j.abs() <= 8e-14);
        }
        let first = model.advance([3.0, 2.0], [0.0, 8.0], 0.5).unwrap();
        let second = model.advance([3.0, 2.0], first.populations, 0.5).unwrap();
        let whole = model.advance([3.0, 2.0], [0.0, 8.0], 1.0).unwrap();
        for i in 0..2 {
            close(second.populations[i], whole.populations[i]);
        }
        close(first.pair_events + second.pair_events, whole.pair_events);
        close(
            first.nonradiative_heat_j + second.nonradiative_heat_j,
            whole.nonradiative_heat_j,
        );
    }

    #[test]
    fn zero_rate_duration_or_triplets_are_exact_identity_steps() {
        for (rate, populations, duration) in [
            (0.0, [3.0, 8.0], 1.0),
            (0.25, [3.0, 8.0], 0.0),
            (0.25, [3.0, 0.0], 1.0),
        ] {
            let step = TripletAnnihilation::try_new(rate, 2.0, 0.5)
                .unwrap()
                .advance([3.0, 2.0], populations, duration)
                .unwrap();
            assert_eq!(step.populations, populations);
            assert_eq!(
                [
                    step.pair_events,
                    step.singlet_production,
                    step.nonradiative_heat_j,
                    step.residual_j
                ],
                [0.0; 4]
            );
        }
    }

    #[test]
    fn volume_controls_density_and_extensive_counts_scale_together() {
        let run = |volume, count| {
            TripletAnnihilation::try_new(0.25, volume, 0.5)
                .unwrap()
                .advance([3.0, 2.0], [0.0, count], 0.5)
                .unwrap()
        };
        let base = run(2.0, 8.0);
        let dilute = run(4.0, 8.0);
        assert!(dilute.populations[1] > base.populations[1]);
        let extensive = run(20.0, 80.0);
        for i in 0..2 {
            close(extensive.populations[i], 10.0 * base.populations[i]);
        }
        close(extensive.pair_events, 10.0 * base.pair_events);
        close(
            extensive.nonradiative_heat_j,
            10.0 * base.nonradiative_heat_j,
        );
    }

    #[test]
    fn singlet_branch_must_be_possible_per_event_even_at_tiny_yield() {
        let tiny_yield = TripletAnnihilation::try_new(1.0, 1.0, 1e-20).unwrap();
        assert!(matches!(
            tiny_yield.advance([4.1, 2.0], [0.0, 2.0], 1.0),
            Err(SingletTripletError::Invalid(_))
        ));
        let dark = TripletAnnihilation::try_new(1.0, 1.0, 0.0)
            .unwrap()
            .advance([4.1, 2.0], [0.0, 2.0], 1.0)
            .unwrap();
        close(dark.nonradiative_heat_j, 4.0 * dark.pair_events);
        assert_eq!(dark.singlet_production, 0.0);
        let resonant = TripletAnnihilation::try_new(1.0, 1.0, 1.0)
            .unwrap()
            .advance([4.0, 2.0], [0.0, 2.0], 1.0)
            .unwrap();
        assert_eq!(resonant.nonradiative_heat_j, 0.0);
        close(resonant.singlet_production, resonant.pair_events);
        let large_energy = TripletAnnihilation::try_new(1.0, 1.0, 1.0)
            .unwrap()
            .advance([f64::MAX, f64::MAX], [0.0, 0.5], 1.0)
            .unwrap();
        close(large_energy.nonradiative_heat_j, 0.125 * f64::MAX);
    }

    #[test]
    fn extreme_exposure_preserves_representable_tail_and_tiny_event_count() {
        let stiff = TripletAnnihilation::try_new(1e100, 1.0, 0.0)
            .unwrap()
            .advance([1.0, 1.0], [0.0, 1e300], 1.0)
            .unwrap();
        close(stiff.populations[1], 5e-101);
        close(stiff.pair_events, 5e299);
        close(stiff.nonradiative_heat_j, 1e300);
        let tiny = TripletAnnihilation::try_new(1e-300, 1e300, 0.0)
            .unwrap()
            .advance([1.0, 1.0], [0.0, 1e300], 1e-300)
            .unwrap();
        close(tiny.pair_events, 1e-300);
        close(tiny.nonradiative_heat_j, 2e-300);
        assert_eq!(tiny.populations[1], 1e300);
        // The public residual exposes energy that endpoint rounding cannot store.
        close(tiny.residual_j, -2e-300);
        let subnormal = TripletAnnihilation::try_new(f64::from_bits(2046_u64 << 52), 1.0, 0.0)
            .unwrap()
            .advance([1.0, 1.0], [0.0, f64::from_bits(3)], 281_474_976_710_656.0)
            .unwrap();
        assert_eq!(subnormal.pair_events, f64::from_bits(1));
        assert_eq!(subnormal.populations[1], f64::from_bits(2));
    }

    #[test]
    fn unrepresentable_events_and_endpoints_refuse_without_a_candidate() {
        let least = f64::from_bits(1);
        let slow = TripletAnnihilation::try_new(least, 1.0, 0.0).unwrap();
        assert!(matches!(
            slow.advance([1.0, 1.0], [0.0, 1.0], least),
            Err(SingletTripletError::Unrepresentable(_))
        ));
        let vanishing_branch = TripletAnnihilation::try_new(1.0, 1.0, least).unwrap();
        assert!(matches!(
            vanishing_branch.advance([1.0, 1.0], [0.0, 1.0], 0.5),
            Err(SingletTripletError::Unrepresentable(_))
        ));
        let model = TripletAnnihilation::try_new(1.0, 1.0, 1.0).unwrap();
        assert!(matches!(
            model.advance([2.0, 1.0], [f64::MAX, 1.0], 1.0),
            Err(SingletTripletError::Unrepresentable(_))
        ));
        assert!(matches!(
            model.advance([0.5, 0.5], [f64::MAX, f64::MAX], 1.0),
            Err(SingletTripletError::Unrepresentable(_))
        ));
    }

    #[test]
    fn invalid_coefficients_states_and_times_are_refused() {
        for (rate, volume, yield_) in [
            (-1.0, 1.0, 0.5),
            (f64::NAN, 1.0, 0.5),
            (1.0, 0.0, 0.5),
            (1.0, f64::INFINITY, 0.5),
            (1.0, 1.0, 1.1),
        ] {
            assert!(TripletAnnihilation::try_new(rate, volume, yield_).is_err());
        }
        let model = TripletAnnihilation::try_new(1.0, 1.0, 0.5).unwrap();
        for (energies, populations, duration) in [
            ([1.0, 2.0], [0.0, 1.0], 1.0),
            ([2.0, 1.0], [-1.0, 1.0], 1.0),
            ([2.0, 1.0], [0.0, 1.0], f64::NAN),
        ] {
            assert!(model.advance(energies, populations, duration).is_err());
        }
    }
}
