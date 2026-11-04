//! Incremental kinematics of the retained mean-drift gas source.
//! No thermal or phase closure is implemented here.

use fs_math::dd::Dd;

#[derive(Clone, Copy, Debug)]
pub(super) struct GasMotion {
    pub momentum: [f64; 3],
    pub velocity: [f64; 3],
    pub kinetic_change_j: f64,
    pub mixing_heat_j: f64,
}

/// Signed mass into gas, at its ORIGINAL donor's retained mean velocity.
#[derive(Clone, Copy, Debug)]
pub(super) struct GasStream {
    pub mass_kg: f64,
    pub velocity_m_s: [f64; 3],
}

/// One frozen-donor phase batch. Withdrawals use the original gas drift and
/// original inventory; incoming vapor cannot finance them. The retained old
/// gas and all incoming streams mix simultaneously. Weighted relative variance
/// deposits their physical loss once, with no subtraction of absolute kinetic
/// energies and no intermediate per-surface gas state.
pub(super) fn advance_many(
    old_mass: f64,
    new_mass: f64,
    momentum: [f64; 3],
    streams: &[GasStream],
) -> Option<GasMotion> {
    let old_velocity = velocity(old_mass, momentum)?;
    let mut condensed = Dd::ZERO;
    let mut impulse = [Dd::ZERO; 3];
    for stream in streams {
        if !stream.mass_kg.is_finite()
            || stream.velocity_m_s.iter().any(|u| !u.is_finite())
            || (stream.mass_kg < 0.0 && stream.velocity_m_s != old_velocity)
        {
            return None;
        }
        if stream.mass_kg < 0.0 {
            condensed = condensed - Dd::from_f64(stream.mass_kg);
        }
        for d in 0..3 {
            impulse[d] =
                impulse[d] + Dd::from_f64(stream.mass_kg) * Dd::from_f64(stream.velocity_m_s[d]);
        }
    }
    let mut retained = (Dd::from_f64(old_mass) - condensed).to_f64();
    if !retained.is_finite() || retained < 0.0 {
        return None;
    }
    let next_momentum = core::array::from_fn(|d| (Dd::from_f64(momentum[d]) + impulse[d]).to_f64());
    let next_velocity = velocity(new_mass, next_momentum)?;
    if (0..3).any(|d| impulse[d] != Dd::ZERO && next_momentum[d] == momentum[d]) {
        return None;
    }
    let mut relative_mean = [0.0; 3];
    let mut loss = Dd::ZERO;
    for stream in streams.iter().filter(|s| s.mass_kg > 0.0) {
        let relative =
            core::array::from_fn::<_, 3, _>(|d| stream.velocity_m_s[d] - old_velocity[d]);
        let combined = retained + stream.mass_kg;
        if !combined.is_finite() || combined <= 0.0 {
            return None;
        }
        if retained == 0.0 {
            relative_mean = relative;
        } else {
            let difference = core::array::from_fn(|d| relative[d] - relative_mean[d]);
            let reduced = retained.min(stream.mass_kg) * (retained.max(stream.mass_kg) / combined);
            loss = loss + Dd::from_f64(kinetic_energy(reduced, difference));
            for d in 0..3 {
                relative_mean[d] += (stream.mass_kg / combined) * difference[d];
            }
        }
        retained = combined;
    }
    // Exact changing-mass identity, using actual represented momentum/mass.
    // Delta K = (u0+u1)/2 . Delta P - (u0.u1)/2 Delta m.
    let dm = Dd::from_f64(new_mass) - Dd::from_f64(old_mass);
    let mut kinetic = Dd::ZERO;
    for d in 0..3 {
        let before = Dd::from_f64(momentum[d]) / Dd::from_f64(old_mass);
        let after = Dd::from_f64(next_momentum[d]) / Dd::from_f64(new_mass);
        let dp = Dd::from_f64(next_momentum[d]) - Dd::from_f64(momentum[d]);
        kinetic = kinetic + (before * Dd::from_f64(0.5) + after * Dd::from_f64(0.5)) * dp
            - dm * Dd::from_f64(0.5) * before * after;
    }
    let kinetic_change_j = kinetic.to_f64();
    let mixing_heat_j = loss.to_f64();
    [
        kinetic_change_j,
        mixing_heat_j,
        kinetic_energy(old_mass, old_velocity),
        kinetic_energy(new_mass, next_velocity),
    ]
    .into_iter()
    .all(f64::is_finite)
    .then_some(GasMotion {
        momentum: next_momentum,
        velocity: next_velocity,
        kinetic_change_j,
        mixing_heat_j,
    })
}

pub(super) fn velocity(mass: f64, momentum: [f64; 3]) -> Option<[f64; 3]> {
    if !mass.is_finite() || mass <= 0.0 {
        return None;
    }
    let result = momentum.map(|p| p / mass);
    (momentum.iter().chain(&result).all(|v| v.is_finite())
        && !(0..3).any(|d| momentum[d] != 0.0 && result[d] == 0.0))
    .then_some(result)
}

pub(super) fn kinetic_energy(mass: f64, velocity: [f64; 3]) -> f64 {
    // Scale before squaring; neither v^2 nor 0.5*m needs to be representable.
    let scale = mass.sqrt() * core::f64::consts::FRAC_1_SQRT_2;
    velocity.iter().map(|u| (scale * u).powi(2)).sum()
}

/// Signed accepted transfer is positive into gas. Actual represented endpoint
/// mass remains separate, so inventory rounding cannot silently disappear.
pub(super) fn advance(
    old_mass: f64,
    new_mass: f64,
    momentum: [f64; 3],
    transferred_mass: f64,
    donor_velocity: [f64; 3],
) -> Option<GasMotion> {
    let old_velocity = velocity(old_mass, momentum)?;
    if !transferred_mass.is_finite()
        || transferred_mass <= -old_mass
        || donor_velocity.iter().any(|u| !u.is_finite())
        || (transferred_mass < 0.0 && donor_velocity != old_velocity)
    {
        return None;
    }
    // Preserve the full product through the accumulation. Rounding the
    // impulse separately would throw away subnormal but representable changes
    // before the compensated sum can observe them.
    let impulse = core::array::from_fn::<_, 3, _>(|d| {
        Dd::from_f64(transferred_mass) * Dd::from_f64(donor_velocity[d])
    });
    let next_momentum = core::array::from_fn(|d| {
        (Dd::from_f64(momentum[d]) + impulse[d]).to_f64()
    });
    if (0..3).any(|d| impulse[d] != Dd::ZERO && next_momentum[d] == momentum[d]) {
        return None;
    }
    let next_velocity = velocity(new_mass, next_momentum)?;
    let mixing_heat_j = if transferred_mass > 0.0 {
        // The receiver and incoming masses may have a representable
        // reduced mass even when their raw sum overflows. Use normalized
        // weights, retaining the existing physical mixing law.
        let larger = old_mass.max(transferred_mass);
        let smaller = old_mass.min(transferred_mass);
        let ratio = smaller / larger;
        let reduced_mass = smaller / (1.0 + ratio);
        kinetic_energy(
            reduced_mass,
            core::array::from_fn(|d| donor_velocity[d] - old_velocity[d]),
        )
    } else {
        0.0
    };
    // Changing-mass identity on the represented endpoint momenta. The old
    // carried-plus-drift expression lost small velocity changes when absolute
    // velocities were large, and computed a rounded mass difference.
    let dm = Dd::from_f64(new_mass) - Dd::from_f64(old_mass);
    let mut kinetic = Dd::ZERO;
    for d in 0..3 {
        let before = Dd::from_f64(momentum[d]) / Dd::from_f64(old_mass);
        let after = Dd::from_f64(next_momentum[d]) / Dd::from_f64(new_mass);
        let dp = Dd::from_f64(next_momentum[d]) - Dd::from_f64(momentum[d]);
        kinetic = kinetic + (before * Dd::from_f64(0.5) + after * Dd::from_f64(0.5)) * dp
            - dm * Dd::from_f64(0.5) * before * after;
    }
    let kinetic_change_j = kinetic.to_f64();
    [
        kinetic_energy(old_mass, old_velocity),
        kinetic_energy(new_mass, next_velocity),
        kinetic_change_j,
        mixing_heat_j,
    ]
    .iter()
    .all(|v| v.is_finite())
    .then_some(GasMotion {
        momentum: next_momentum,
        velocity: next_velocity,
        kinetic_change_j,
        mixing_heat_j,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: f64, b: f64, tolerance: f64) {
        assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
    }

    #[test]
    fn evaporation_mixes_incoming_momentum_and_only_relative_energy() {
        let mass = 2.0;
        let incoming = 0.5;
        let velocity = [3.0, -1.0, 2.0];
        let donor = [-2.0, 4.0, 1.0];
        let initial = velocity.map(|u| mass * u);
        let next = advance(mass, mass + incoming, initial, incoming, donor).unwrap();
        for d in 0..3 {
            near(next.momentum[d], initial[d] + incoming * donor[d], 1e-15);
        }
        near(next.mixing_heat_j, 10.2, 1e-14);
        near(
            next.kinetic_change_j + next.mixing_heat_j,
            kinetic_energy(incoming, donor),
            1e-14,
        );
    }

    #[test]
    fn condensation_withdraws_the_actual_gas_mean_without_a_jet_recoil() {
        let initial = [6.0, -2.0, 4.0];
        let donor = [3.0, -1.0, 2.0];
        let next = advance(2.0, 1.5, initial, -0.5, donor).unwrap();
        assert_eq!(next.velocity, donor);
        assert_eq!(next.mixing_heat_j, 0.0);
        near(next.kinetic_change_j, -kinetic_energy(0.5, donor), 1e-14);
        assert!(advance(2.0, 1.5, initial, -0.5, [0.0; 3]).is_none());
    }

    #[test]
    fn common_boost_preserves_receiver_loss_and_transforms_stream_energy() {
        for boost in [[0.0; 3], [100.0, -200.0, 50.0], [-1e4, 1e4, 2e4]] {
            let velocity = core::array::from_fn(|d| [3.0, -1.0, 2.0][d] + boost[d]);
            let donor = core::array::from_fn(|d| [-2.0, 4.0, 1.0][d] + boost[d]);
            let next = advance(2.0, 2.5, velocity.map(|u| 2.0 * u), 0.5, donor).unwrap();
            near(next.mixing_heat_j, 10.2, 1e-12);
            near(
                next.kinetic_change_j + next.mixing_heat_j,
                kinetic_energy(0.5, donor),
                1e-7,
            );
        }
    }

    #[test]
    fn zero_flux_and_equal_drift_create_no_physical_mixing_heat() {
        let p = [6.0, -2.0, 4.0];
        let zero = advance(2.0, 2.0, p, 0.0, [0.0; 3]).unwrap();
        assert_eq!(zero.momentum, p);
        assert_eq!(zero.kinetic_change_j, 0.0);
        assert_eq!(zero.mixing_heat_j, 0.0);
        let equal = advance(2.0, 2.5, p, 0.5, [3.0, -1.0, 2.0]).unwrap();
        assert_eq!(equal.mixing_heat_j, 0.0);
        assert_eq!(equal.velocity, [3.0, -1.0, 2.0]);
    }

    #[test]
    fn represented_mass_error_remains_in_the_mechanical_balance() {
        let exact = advance(2.0, 2.5, [6.0, 0.0, 0.0], 0.5, [-2.0, 0.0, 0.0]).unwrap();
        let rounded = advance(2.0, 2.500001, [6.0, 0.0, 0.0], 0.5, [-2.0, 0.0, 0.0]).unwrap();
        assert_eq!(rounded.momentum, exact.momentum);
        assert_eq!(rounded.mixing_heat_j, exact.mixing_heat_j);
        assert_ne!(rounded.kinetic_change_j, exact.kinetic_change_j);
    }

    #[test]
    fn weighted_kinetics_survive_extreme_scales_and_invalid_states_refuse() {
        near(kinetic_energy(1e-300, [1e150, 0.0, 0.0]), 0.5, 1e-15);
        assert!(velocity(0.0, [0.0; 3]).is_none());
        assert!(velocity(1.0, [f64::INFINITY, 0.0, 0.0]).is_none());
        assert!(advance(1.0, 0.0, [1.0; 3], -1.0, [1.0; 3]).is_none());
        assert!(advance(1.0, 2.0, [0.0; 3], f64::NAN, [0.0; 3]).is_none());
        // Small but physically nonzero impulse cannot vanish into a large
        // existing momentum without a visible refusal.
        assert!(advance(1e20, 1e20, [1e30, 0.0, 0.0], 1e-12, [1.0, 0.0, 0.0]).is_none());
    }

    fn stream(mass_kg: f64, velocity_m_s: [f64; 3]) -> GasStream {
        GasStream {
            mass_kg,
            velocity_m_s,
        }
    }

    #[test]
    fn single_phase_matches_distributed_under_boost_and_inventory_rounding() {
        // Both public paths must use the same represented-momentum account,
        // including a nontrivial incoming relative velocity and a rounded
        // endpoint inventory. The thermal loss is frame-independent.
        for boost in [[0.0; 3], [1024.0, -4096.0, 256.0]] {
            let old: [f64; 3] =
                core::array::from_fn(|d| [3.0, -1.0, 2.0][d] + boost[d]);
            let donor: [f64; 3] =
                core::array::from_fn(|d| [-2.0, 4.0, 1.0][d] + boost[d]);
            let p = old.map(|u| 2.0 * u);
            for new_mass in [2.5, 2.500001] {
                let single = advance(2.0, new_mass, p, 0.5, donor).unwrap();
                let batch = advance_many(2.0, new_mass, p, &[stream(0.5, donor)]).unwrap();
                assert_eq!(single.momentum, batch.momentum);
                assert_eq!(single.velocity, batch.velocity);
                near(single.mixing_heat_j, batch.mixing_heat_j, 1e-12);
                near(single.kinetic_change_j, batch.kinetic_change_j, 1e-9);
            }
        }
    }

    #[test]
    fn mixing_reduced_mass_remains_representable_for_large_inventory() {
        // A mixture's reduced mass can fit even when adding both masses
        // directly would overflow. No artificial clamping of velocity.
        let old = 1.0e308;
        let incoming = 1.0e308;
        let expected = 0.5e308;
        let larger = old.max(incoming);
        let smaller = old.min(incoming);
        let reduced = smaller / (1.0 + smaller / larger);
        assert_eq!(reduced, expected);
        // Stable weighted energy remains finite without forming m*u*u
        // in ordinary precision.
        near(kinetic_energy(reduced, [1.0e-154, 0.0, 0.0]), 0.25, 1e-15);
    }

    #[test]
    fn single_phase_rejects_a_lost_nonzero_impulse() {
        let mass = 1e20;
        let momentum = [1e30, 0.0, 0.0];
        let transfer = 1e-12;
        assert!(advance(mass, mass, momentum, transfer, [1.0, 0.0, 0.0]).is_none());
        assert!(advance_many(mass, mass, momentum, &[stream(transfer, [1.0, 0.0, 0.0])]).is_none());
    }

    #[test]
    fn distributed_opposing_phase_closes_independent_rational_inventory() {
        let streams = [
            stream(-0.5, [3.0, -1.0, 2.0]),
            stream(0.25, [-2.0, 4.0, 1.0]),
            stream(0.75, [1.0, 2.0, -3.0]),
        ];
        let next = advance_many(2.0, 2.5, [6.0, -2.0, 4.0], &streams).unwrap();
        assert_eq!(next.momentum, [4.75, 1.0, 1.0]);
        assert_eq!(next.velocity, [1.9, 0.4, 0.4]);
        near(next.kinetic_change_j, -727.0 / 80.0, 2e-15);
        near(next.mixing_heat_j, 1077.0 / 80.0, 4e-15);
        near(
            next.kinetic_change_j + next.mixing_heat_j,
            35.0 / 8.0,
            4e-15,
        );
    }

    #[test]
    fn equal_gross_exchange_changes_momentum_and_heat_at_zero_net_mass() {
        let next = advance_many(
            2.0,
            2.0,
            [6.0, 0.0, 0.0],
            &[stream(-0.5, [3.0, 0.0, 0.0]), stream(0.5, [-1.0, 0.0, 0.0])],
        )
        .unwrap();
        assert_eq!(next.momentum, [4.0, 0.0, 0.0]);
        assert_eq!(next.kinetic_change_j, -5.0);
        near(next.mixing_heat_j, 3.0, 1e-15);
        let dispersion = advance_many(
            2.0,
            4.0,
            [0.0; 3],
            &[stream(1.0, [3.0, 0.0, 0.0]), stream(1.0, [-3.0, 0.0, 0.0])],
        )
        .unwrap();
        assert_eq!(dispersion.momentum, [0.0; 3]);
        assert_eq!(dispersion.kinetic_change_j, 0.0);
        near(dispersion.mixing_heat_j, 9.0, 4e-15);
    }

    #[test]
    fn distributed_mixing_keeps_relative_loss_under_boost_and_permutation() {
        for boost in [[0.0; 3], [1048576.0, -2097152.0, 524288.0]] {
            let old = core::array::from_fn(|d| [3.0, -1.0, 2.0][d] + boost[d]);
            let a = core::array::from_fn(|d| [-2.0, 4.0, 1.0][d] + boost[d]);
            let b = core::array::from_fn(|d| [1.0, 2.0, -3.0][d] + boost[d]);
            let streams = [stream(-0.5, old), stream(0.25, a), stream(0.75, b)];
            let p = old.map(|u| 2.0 * u);
            let next = advance_many(2.0, 2.5, p, &streams).unwrap();
            let reverse = advance_many(2.0, 2.5, p, &[streams[2], streams[1], streams[0]]).unwrap();
            assert_eq!(next.momentum, reverse.momentum);
            near(next.mixing_heat_j, 1077.0 / 80.0, 4e-15);
            near(next.mixing_heat_j, reverse.mixing_heat_j, 4e-15);
            // Exact transform Delta K' = Delta K + boost.Delta P + |boost|² Delta m/2.
            let transformed = -727.0 / 80.0
                + (0..3)
                    .map(|d| boost[d] * [-1.25, 3.0, -3.0][d])
                    .sum::<f64>()
                + 0.25 * boost.iter().map(|u| u * u).sum::<f64>();
            near(next.kinetic_change_j, transformed, 1e-12);
        }
    }

    #[test]
    fn common_gas_donor_is_not_financed_by_incoming_vapor() {
        assert!(
            advance_many(
                1.0,
                1.5,
                [2.0, 0.0, 0.0],
                &[stream(-1.5, [2.0, 0.0, 0.0]), stream(2.0, [0.0; 3]),]
            )
            .is_none()
        );
        assert!(
            advance_many(1.0, 0.5, [2.0, 0.0, 0.0], &[stream(-0.5, [1.0, 0.0, 0.0]),]).is_none()
        );
        let replaced = advance_many(
            1.0,
            1.0,
            [2.0, 0.0, 0.0],
            &[
                stream(-1.0, [2.0, 0.0, 0.0]),
                stream(0.5, [3.0, 0.0, 0.0]),
                stream(0.5, [-3.0, 0.0, 0.0]),
            ],
        )
        .unwrap();
        near(replaced.mixing_heat_j, 4.5, 2e-15);
        assert_eq!(replaced.kinetic_change_j, -2.0);
    }

    #[test]
    fn singleton_and_thermal_only_batches_retain_the_existing_physical_law() {
        let p = [6.0, -2.0, 4.0];
        for (dm, donor) in [(0.5, [-2.0, 4.0, 1.0]), (-0.5, [3.0, -1.0, 2.0])] {
            let single = advance(2.0, 2.0 + dm, p, dm, donor).unwrap();
            let batch = advance_many(2.0, 2.0 + dm, p, &[stream(dm, donor)]).unwrap();
            assert_eq!(single.momentum, batch.momentum);
            near(single.kinetic_change_j, batch.kinetic_change_j, 3e-15);
            near(single.mixing_heat_j, batch.mixing_heat_j, 3e-15);
        }
        let zero = advance_many(2.0, 2.0, p, &[stream(0.0, [9.0; 3])]).unwrap();
        assert_eq!(zero.momentum, p);
        assert_eq!(zero.kinetic_change_j, 0.0);
        assert_eq!(zero.mixing_heat_j, 0.0);
    }

    #[test]
    fn represented_mass_error_and_disappearing_impulse_cannot_hide() {
        let streams = [stream(0.5, [-2.0, 0.0, 0.0])];
        let exact = advance_many(2.0, 2.5, [6.0, 0.0, 0.0], &streams).unwrap();
        let rounded = advance_many(2.0, 2.500001, [6.0, 0.0, 0.0], &streams).unwrap();
        assert_eq!(exact.momentum, rounded.momentum);
        assert_eq!(exact.mixing_heat_j, rounded.mixing_heat_j);
        assert_ne!(exact.kinetic_change_j, rounded.kinetic_change_j);
        assert!(
            advance_many(
                1e20,
                1e20,
                [1e30, 0.0, 0.0],
                &[stream(1e-12, [1.0, 0.0, 0.0]),]
            )
            .is_none()
        );
        assert!(advance_many(1.0, 1.5, [0.0; 3], &[stream(0.5, [f64::NAN; 3])]).is_none());
    }

    #[test]
    fn distributed_weighted_kinetics_remain_finite_at_extreme_mass_scale() {
        let scale = 1e-200;
        let speed = 1e90;
        let next = advance_many(
            2.0 * scale,
            2.5 * scale,
            [6.0 * scale * speed, 0.0, 0.0],
            &[
                stream(-0.5 * scale, [3.0 * speed, 0.0, 0.0]),
                stream(scale, [-speed, 0.0, 0.0]),
            ],
        )
        .unwrap();
        near(next.mixing_heat_j / (scale * speed * speed), 4.8, 3e-15);
        near(
            next.kinetic_change_j / (scale * speed * speed),
            -6.55,
            3e-15,
        );
        // The mean-speed SUM overflows, although the mass-weighted kinetic
        // increment is representable. Halve each speed before adding it.
        let mass = 2e-320;
        let incoming = 1e-320;
        let p = [mass * 1.2e308, 0.0, 0.0];
        let u = velocity(mass, p).unwrap();
        let extreme = advance_many(mass, mass + incoming, p, &[stream(incoming, u)]).unwrap();
        assert_eq!(extreme.mixing_heat_j, 0.0);
        near(
            extreme.kinetic_change_j / kinetic_energy(incoming, u),
            1.0,
            1e-13,
        );
    }
}
