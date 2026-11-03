use super::*;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn law(darcy: bool, wall: f64, rate: f64) -> WallDrag {
    if darcy {
        WallDrag::MovingDarcy {
            friction_factor: rate,
            hydraulic_diameter_m: 0.5,
            wall_velocity_m_s: wall,
        }
    } else {
        WallDrag::MovingLinear {
            rate_per_s: rate,
            wall_velocity_m_s: wall,
        }
    }
}

#[test]
fn moving_wall_exact_entrainment_braking_and_slip_heating() {
    for darcy in [false, true] {
        for (velocity, wall) in [(0.0, 0.8), (1.1, 0.3), (0.7, -0.4), (-0.6, -0.2)] {
            for rate in [0.0, 1e-8, 2.0, 1e6] {
                let (density, volume, dt) = (1.7, 0.3, 0.2);
                let momentum = density * velocity;
                let energy = 8.0 + 0.5 * density * velocity * velocity;
                let update = advance(
                    law(darcy, wall, rate),
                    momentum,
                    energy,
                    density,
                    velocity,
                    volume,
                    dt,
                )
                .unwrap();
                let initial_slip: f64 = velocity - wall;
                let slip = if darcy {
                    initial_slip / (1.0 + rate * initial_slip.abs() * dt)
                } else {
                    initial_slip * (-rate * dt).exp()
                };
                close(update.momentum / density, wall + slip, 4e-15);
                let impulse = (update.momentum - momentum) * volume;
                close(update.work_j, wall * impulse, 2e-15);
                close((update.energy - energy) * volume, update.work_j, 5e-15);
                let delta_kinetic =
                    0.5 * volume * (update.momentum * update.momentum - momentum * momentum)
                        / density;
                close(delta_kinetic + update.heating_j, update.work_j, 3e-15);
                close(
                    update.heating_j,
                    0.5 * density * volume * (initial_slip * initial_slip - slip * slip),
                    3e-15,
                );
                assert!(update.heating_j >= 0.0);
                assert!(update.momentum / density >= velocity.min(wall) - 1e-15);
                assert!(update.momentum / density <= velocity.max(wall) + 1e-15);
            }
        }
    }
}

#[test]
fn moving_wall_source_is_galilean_covariant_with_invariant_dissipation() {
    let (density, volume, dt, velocity, wall) = (1.25, 0.7, 0.2, -0.3, 0.8);
    for darcy in [false, true] {
        let baseline = advance(
            law(darcy, wall, 1.8),
            density * velocity,
            9.0 + 0.5 * density * velocity * velocity,
            density,
            velocity,
            volume,
            dt,
        )
        .unwrap();
        for boost in [-100.0, -2.0, 0.0, 5.0, 100.0] {
            let boosted_velocity = velocity + boost;
            let changed = advance(
                law(darcy, wall + boost, 1.8),
                density * boosted_velocity,
                9.0 + 0.5 * density * boosted_velocity * boosted_velocity,
                density,
                boosted_velocity,
                volume,
                dt,
            )
            .unwrap();
            close(changed.momentum, baseline.momentum + density * boost, 5e-14);
            close(
                changed.energy,
                baseline.energy + boost * baseline.momentum + 0.5 * density * boost * boost,
                3e-12,
            );
            close(changed.heating_j, baseline.heating_j, 3e-14);
            close(
                changed.work_j,
                baseline.work_j + boost * (baseline.momentum - density * velocity) * volume,
                3e-12,
            );
        }
    }
}

#[test]
fn moving_wall_no_slip_zero_coupling_and_stationary_limit_preserve_state() {
    for darcy in [false, true] {
        for (velocity, wall, rate, dt) in [
            (0.7, 0.7, 4.0, 0.2),
            (0.7, -0.3, 0.0, 0.2),
            (0.7, -0.3, 4.0, 0.0),
        ] {
            let p = 1.3 * velocity;
            let update = advance(law(darcy, wall, rate), p, 9.0, 1.3, velocity, 0.4, dt).unwrap();
            assert_eq!(update.momentum, p);
            assert_eq!(update.energy, 9.0);
            assert_eq!(update.work_j, 0.0);
            assert_eq!(update.heating_j, 0.0);
        }
        let fixed = if darcy {
            WallDrag::Darcy {
                friction_factor: 2.0,
                hydraulic_diameter_m: 0.5,
            }
        } else {
            WallDrag::Linear { rate_per_s: 2.0 }
        };
        let static_update = advance(fixed, 0.91, 9.0, 1.3, 0.7, 0.4, 0.2).unwrap();
        let moving_zero = advance(law(darcy, 0.0, 2.0), 0.91, 9.0, 1.3, 0.7, 0.4, 0.2).unwrap();
        assert_eq!(
            moving_zero.momentum.to_bits(),
            static_update.momentum.to_bits()
        );
        assert_eq!(moving_zero.energy.to_bits(), static_update.energy.to_bits());
        assert_eq!(
            moving_zero.heating_j.to_bits(),
            static_update.heating_j.to_bits()
        );
    }
}

#[test]
fn moving_wall_weak_drive_retains_small_representable_momentum() {
    for darcy in [false, true] {
        let update = advance(law(darcy, 1.0, 1e-18), 0.0, 1.0, 1.0, 0.0, 1.0, 1.0).unwrap();
        close(update.momentum, 1e-18, 1e-33);
        close(update.work_j, 1e-18, 1e-33);
        close(update.heating_j, 1e-18, 1e-33);
        assert!(update.momentum > 0.0);
    }
}
