//! Incremental kinematics of the retained mean-drift gas source.
//! No thermal or phase closure is implemented here.

#[derive(Clone, Copy, Debug)]
pub(super) struct GasMotion {
    pub momentum: [f64; 3],
    pub velocity: [f64; 3],
    pub kinetic_change_j: f64,
    pub mixing_heat_j: f64,
}

pub(super) fn velocity(mass: f64, momentum: [f64; 3]) -> Option<[f64; 3]> {
    let result = momentum.map(|p| p / mass);
    (mass.is_finite()
        && mass > 0.0
        && momentum.iter().chain(&result).all(|v| v.is_finite())
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
    let next_momentum =
        core::array::from_fn(|d| momentum[d] + transferred_mass * donor_velocity[d]);
    let next_velocity = velocity(new_mass, next_momentum)?;
    let mixing_heat_j = if transferred_mass > 0.0 {
        let mixed_mass = old_mass + transferred_mass;
        if !mixed_mass.is_finite() || mixed_mass <= 0.0 {
            return None;
        }
        let reduced_mass =
            old_mass.min(transferred_mass) * (old_mass.max(transferred_mass) / mixed_mass);
        kinetic_energy(
            reduced_mass,
            core::array::from_fn(|d| donor_velocity[d] - old_velocity[d]),
        )
    } else {
        0.0
    };
    let actual_dm = new_mass - old_mass;
    let carried = actual_dm.signum() * kinetic_energy(actual_dm.abs(), old_velocity);
    let drift = (0..3)
        .map(|d| {
            let du = next_velocity[d] - old_velocity[d];
            (new_mass * du) * (old_velocity[d] + 0.5 * du)
        })
        .sum::<f64>();
    let kinetic_change_j = carried + drift;
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
    }
}
