//! Conservative co-located stream reduction in the original liquid datum.
//!
//! A mean velocity alone loses kinetic energy. Retain the mass-weighted
//! central second moment `q = sum(m_i |v_i-v_mean|²)/(2 sum(m_i))` alongside
//! mean enthalpy and velocity. Use positive mass magnitudes for either stream
//! direction. Complete inflow mixing turns `dm*q` into heat; powered outflow
//! requires `r*q` actuator work for removed mass `r`. Both retain the same
//! term in their signed carried kinetic-energy account.
//! Pairwise central-moment combination avoids subtracting large kinetic
//! energies and is covariant under a common velocity or caloric-datum shift.

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct IncomingMixture {
    pub(crate) mass_kg: f64,
    pub(crate) velocity_m_s: [f64; 3],
    pub(crate) specific_enthalpy_j_kg: f64,
    pub(crate) heat_j: f64,
    pub(crate) specific_kinetic_variance_j_kg: f64,
}

impl IncomingMixture {
    pub(crate) fn new(
        mass_kg: f64,
        velocity_m_s: [f64; 3],
        specific_enthalpy_j_kg: f64,
        heat_j: f64,
    ) -> Option<Self> {
        (mass_kg.is_finite()
            && mass_kg > 0.0
            && velocity_m_s.into_iter().all(f64::is_finite)
            && specific_enthalpy_j_kg.is_finite()
            && heat_j.is_finite())
        .then_some(Self {
            mass_kg,
            velocity_m_s,
            specific_enthalpy_j_kg,
            heat_j,
            specific_kinetic_variance_j_kg: 0.0,
        })
    }

    /// One allocation-free merge. Neither source is changed on refusal.
    pub(crate) fn combine(&mut self, incoming: Self) -> Option<()> {
        let mass = self.mass_kg + incoming.mass_kg;
        if !mass.is_finite() || mass <= 0.0 {
            return None;
        }
        let fraction = incoming.mass_kg / mass;
        let retained = self.mass_kg / mass;
        let delta: [f64; 3] =
            core::array::from_fn(|d| incoming.velocity_m_s[d] - self.velocity_m_s[d]);
        let velocity = core::array::from_fn(|d| self.velocity_m_s[d] + fraction * delta[d]);
        let enthalpy = self.specific_enthalpy_j_kg
            + fraction * (incoming.specific_enthalpy_j_kg - self.specific_enthalpy_j_kg);
        let variance = retained * self.specific_kinetic_variance_j_kg
            + fraction * incoming.specific_kinetic_variance_j_kg
            + 0.5 * fraction * retained * delta.iter().map(|u| u * u).sum::<f64>();
        let heat = self.heat_j + incoming.heat_j;
        if velocity.iter().any(|u| !u.is_finite())
            || !enthalpy.is_finite()
            || !variance.is_finite()
            || !heat.is_finite()
        {
            return None;
        }
        *self = Self {
            mass_kg: mass,
            velocity_m_s: velocity,
            specific_enthalpy_j_kg: enthalpy,
            heat_j: heat,
            specific_kinetic_variance_j_kg: variance,
        };
        Some(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64, scale: f64) {
        assert!(
            (actual - expected).abs() <= 128.0 * f64::EPSILON * scale,
            "actual={actual:e}, expected={expected:e}, scale={scale:e}"
        );
    }

    fn streams() -> [IncomingMixture; 3] {
        [
            IncomingMixture::new(0.2, [0.8, -0.4, 0.2], 120_000.0, 2.0).unwrap(),
            IncomingMixture::new(0.3, [-0.6, 0.3, -0.1], 240_000.0, -1.0).unwrap(),
            IncomingMixture::new(0.1, [0.2, 0.5, -0.4], 180_000.0, 0.5).unwrap(),
        ]
    }

    fn reduced(sources: [IncomingMixture; 3], order: [usize; 3]) -> IncomingMixture {
        let mut result = sources[order[0]];
        for i in order.into_iter().skip(1) {
            result.combine(sources[i]).unwrap();
        }
        result
    }

    #[test]
    fn distinct_stream_momenta_and_both_energies_survive_one_native_mass_fraction() {
        let sources = streams();
        let mixed = reduced(sources, [0, 1, 2]);
        let total_mass: f64 = sources.iter().map(|s| s.mass_kg).sum();
        close(mixed.mass_kg, total_mass, total_mass);
        for represented_fraction in [0.37, 1.0, 1.0 + 2.0 * f64::EPSILON] {
            let dm = mixed.mass_kg * represented_fraction;
            for d in 0..3 {
                let momentum: f64 = sources
                    .iter()
                    .map(|s| represented_fraction * s.mass_kg * s.velocity_m_s[d])
                    .sum();
                close(dm * mixed.velocity_m_s[d], momentum, dm);
            }
            let enthalpy: f64 = sources
                .iter()
                .map(|s| represented_fraction * s.mass_kg * s.specific_enthalpy_j_kg)
                .sum();
            close(dm * mixed.specific_enthalpy_j_kg, enthalpy, enthalpy.abs());
            let kinetic: f64 = sources
                .iter()
                .map(|s| {
                    0.5 * represented_fraction
                        * s.mass_kg
                        * s.velocity_m_s.iter().map(|u| u * u).sum::<f64>()
                })
                .sum();
            let resolved = 0.5 * dm * mixed.velocity_m_s.iter().map(|u| u * u).sum::<f64>();
            close(
                resolved + dm * mixed.specific_kinetic_variance_j_kg,
                kinetic,
                kinetic,
            );
            assert!(kinetic - resolved > 0.05 * represented_fraction);
        }
        // Authored heat is not mass transport and is never fraction-scaled.
        assert_eq!(mixed.heat_j, 1.5);
    }

    #[test]
    fn opposite_jets_keep_the_heat_that_a_zero_mean_velocity_would_erase() {
        let left = IncomingMixture::new(0.3, [2.0, -1.0, 0.5], 100.0, 0.0).unwrap();
        let right = IncomingMixture::new(0.3, [-2.0, 1.0, -0.5], 300.0, 0.0).unwrap();
        let mut mixed = left;
        mixed.combine(right).unwrap();
        assert_eq!(mixed.velocity_m_s, [0.0; 3]);
        assert_eq!(mixed.specific_enthalpy_j_kg, 200.0);
        assert_eq!(mixed.specific_kinetic_variance_j_kg, 2.625);
        close(
            mixed.mass_kg * mixed.specific_kinetic_variance_j_kg,
            1.575,
            1.575,
        );
    }

    #[test]
    fn permutation_density_boost_and_caloric_datum_covariance() {
        let baseline = reduced(streams(), [0, 1, 2]);
        let orders = [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ];
        for scale in [1e-10, 0.4, 1.0, 1e10] {
            for boost in [[0.0; 3], [7.0, -3.0, 2.0]] {
                for datum in [0.0, -1e8, 1e8] {
                    let transformed = streams().map(|s| {
                        IncomingMixture::new(
                            scale * s.mass_kg,
                            core::array::from_fn(|d| s.velocity_m_s[d] + boost[d]),
                            s.specific_enthalpy_j_kg + datum,
                            scale * s.heat_j,
                        )
                        .unwrap()
                    });
                    for order in orders {
                        let actual = reduced(transformed, order);
                        close(actual.mass_kg, scale * baseline.mass_kg, scale);
                        for d in 0..3 {
                            close(
                                actual.velocity_m_s[d],
                                baseline.velocity_m_s[d] + boost[d],
                                8.0,
                            );
                        }
                        close(
                            actual.specific_kinetic_variance_j_kg,
                            baseline.specific_kinetic_variance_j_kg,
                            1.0,
                        );
                        close(
                            actual.specific_enthalpy_j_kg,
                            baseline.specific_enthalpy_j_kg + datum,
                            1e8,
                        );
                        close(actual.heat_j, scale * baseline.heat_j, 2.0 * scale);
                    }
                }
            }
        }
    }

    #[test]
    fn singleton_is_unchanged_and_overflow_refuses_before_modifying_the_group() {
        let source = IncomingMixture::new(0.2, [-0.0, 0.03, -0.02], -1e9, -0.0).unwrap();
        assert_eq!(source.mass_kg.to_bits(), 0.2_f64.to_bits());
        assert_eq!(source.velocity_m_s[0].to_bits(), (-0.0_f64).to_bits());
        assert_eq!(source.heat_j.to_bits(), (-0.0_f64).to_bits());
        assert_eq!(source.specific_kinetic_variance_j_kg, 0.0);
        assert!(IncomingMixture::new(0.0, [0.0; 3], 0.0, 0.0).is_none());
        assert!(IncomingMixture::new(1.0, [f64::NAN, 0.0, 0.0], 0.0, 0.0).is_none());
        for pair in [
            [IncomingMixture::new(f64::MAX, [0.0; 3], 0.0, 0.0).unwrap(); 2],
            [IncomingMixture::new(1.0, [0.0; 3], 0.0, f64::MAX).unwrap(); 2],
            [
                IncomingMixture::new(1.0, [f64::MAX, 0.0, 0.0], 0.0, 0.0).unwrap(),
                IncomingMixture::new(1.0, [-f64::MAX, 0.0, 0.0], 0.0, 0.0).unwrap(),
            ],
        ] {
            let mut mixed = pair[0];
            assert!(mixed.combine(pair[1]).is_none());
            assert_eq!(mixed, pair[0]);
        }
    }
}
