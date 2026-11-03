use super::*;

fn near(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}
fn vector_near(a: Vec3, b: Vec3, relative: f64) {
    near(a.x, b.x, relative);
    near(a.y, b.y, relative);
    near(a.z, b.z, relative);
}
fn state() -> RocketState {
    RocketState {
        position_m: Vec3::new(2.0, -3.0, 5.0),
        velocity_m_s: Vec3::new(10.0, 40.0, -20.0),
        dry_mass_kg: 12.0,
        propellant_mass_kg: 8.0,
    }
}
fn forcing() -> ExhaustForcing {
    let w = Vec3::new(400.0, -30.0, 100.0);
    ExhaustForcing {
        mass_flow_kg_s: 2.0,
        relative_momentum_rate_n: w.scale(2.0),
        relative_energy_rate_w: 2.0 * (200_000.0 + 0.5 * w.norm_squared()),
    }
}

#[test]
fn physical_momentum_and_energy_match_independent_inertial_flux_quadrature() {
    let (s, f, dt) = (state(), forcing(), 2.0);
    let step = s.advance_with_exhaust(f, dt, || false).unwrap();
    let m0 = s.dry_mass_kg + s.propellant_mass_kg;
    let dm = f.mass_flow_kg_s * dt;
    let m1 = m0 - dm;
    let w = f.relative_momentum_rate_n.scale(f.mass_flow_kg_s.recip());
    let expected_dv = w.scale(-(m0 / m1).ln());
    vector_near(
        step.mechanics.state_after.velocity_m_s,
        s.velocity_m_s.add(expected_dv),
        2e-13,
    );

    // G1: independent Simpson integration of the physical particle moments.
    // Uses std log and 1024 time panels, not the production antiderivatives.
    let panels = 1024;
    let h = dt / f64::from(panels);
    let mut impulse = Vec3::ZERO;
    let (mut energy, mut mean_squared, mut mean_velocity) = (0.0, 0.0, Vec3::ZERO);
    for i in 0..=panels {
        let t = f64::from(i) * h;
        let mass = m0 - f.mass_flow_kg_s * t;
        let v = s.velocity_m_s.sub(w.scale((m0 / mass).ln()));
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        impulse = impulse.add(
            v.scale(f.mass_flow_kg_s)
                .add(f.relative_momentum_rate_n)
                .scale(weight),
        );
        energy += weight
            * (f.relative_energy_rate_w
                + v.dot(f.relative_momentum_rate_n)
                + 0.5 * f.mass_flow_kg_s * v.norm_squared());
        mean_velocity = mean_velocity.add(v.scale(weight));
        mean_squared += weight * v.norm_squared();
    }
    vector_near(
        step.exhaust_momentum_world_kg_m_s,
        impulse.scale(h / 3.0),
        2e-12,
    );
    near(step.exhaust_energy_world_j, energy * (h / 3.0), 2e-12);
    vector_near(
        step.mean_velocity_world_m_s,
        mean_velocity.scale(h / (3.0 * dt)),
        2e-12,
    );
    near(
        step.mean_speed_squared_m2_s2,
        mean_squared * (h / (3.0 * dt)),
        2e-12,
    );
    vector_near(
        step.mechanics
            .state_after
            .velocity_m_s
            .scale(m1)
            .add(step.exhaust_momentum_world_kg_m_s),
        s.velocity_m_s.scale(m0),
        2e-12,
    );
    let dk = 0.5 * m1 * step.mechanics.state_after.velocity_m_s.norm_squared()
        - 0.5 * m0 * s.velocity_m_s.norm_squared();
    near(
        dk + step.exhaust_energy_world_j,
        f.relative_energy_rate_w * dt,
        2e-13,
    );
    assert!(step.momentum_residual_kg_m_s.norm_squared() < 1e-22);
    assert!(step.kinetic_energy_residual_j.abs() < 1e-8);
}

#[test]
fn boosts_transform_the_exported_particle_moments_and_preserve_internal_energy() {
    let (s, f, dt) = (state(), forcing(), 2.0);
    let base = s.advance_with_exhaust(f, dt, || false).unwrap();
    let boost = Vec3::new(1000.0, -500.0, 800.0);
    let shifted = RocketState {
        velocity_m_s: s.velocity_m_s.add(boost),
        ..s
    }
    .advance_with_exhaust(f, dt, || false)
    .unwrap();
    let dm = base.mechanics.burned_propellant_mass_kg;
    assert_eq!(base.relative_energy_j, shifted.relative_energy_j);
    vector_near(
        shifted.mechanics.state_after.velocity_m_s,
        base.mechanics.state_after.velocity_m_s.add(boost),
        2e-13,
    );
    vector_near(
        shifted.mechanics.state_after.position_m,
        base.mechanics.state_after.position_m.add(boost.scale(dt)),
        2e-13,
    );
    vector_near(
        shifted.exhaust_momentum_world_kg_m_s,
        base.exhaust_momentum_world_kg_m_s.add(boost.scale(dm)),
        2e-13,
    );
    near(
        shifted.exhaust_energy_world_j,
        base.exhaust_energy_world_j
            + boost.dot(base.exhaust_momentum_world_kg_m_s)
            + 0.5 * dm * boost.norm_squared(),
        2e-13,
    );

    // The kinetic chart survives a chemical reference far larger than recoil.
    let signed = s
        .advance_with_exhaust(
            ExhaustForcing {
                relative_energy_rate_w: -1e20,
                ..f
            },
            dt,
            || false,
        )
        .unwrap();
    assert_eq!(
        signed.relative_kinetic_transport_energy_j,
        base.relative_kinetic_transport_energy_j
    );
    assert_eq!(
        signed.kinetic_energy_residual_j,
        base.kinetic_energy_residual_j
    );
    assert!(signed.exhaust_energy_world_j < 0.0);
}

#[test]
fn opposed_streams_lose_mass_and_energy_without_accelerating() {
    let s = state();
    let dt = 1.0;
    let f = ExhaustForcing {
        mass_flow_kg_s: 2.0,
        relative_momentum_rate_n: Vec3::ZERO,
        relative_energy_rate_w: 600_000.0,
    };
    let step = s.advance_with_exhaust(f, dt, || false).unwrap();
    assert_eq!(step.mechanics.state_after.velocity_m_s, s.velocity_m_s);
    assert_eq!(step.mechanics.state_after.propellant_mass_kg, 6.0);
    assert_eq!(step.mechanics.thrust_impulse_n_s, Vec3::ZERO);
    assert_eq!(
        step.exhaust_momentum_world_kg_m_s,
        s.velocity_m_s.scale(2.0)
    );
    near(
        step.exhaust_energy_world_j,
        f.relative_energy_rate_w + s.velocity_m_s.norm_squared(),
        1e-15,
    );
    assert_eq!(step.mean_speed_squared_m2_s2, s.velocity_m_s.norm_squared());
    assert_eq!(step.momentum_residual_kg_m_s, Vec3::ZERO);
    assert_eq!(step.kinetic_energy_residual_j, 0.0);
}

#[test]
fn split_steps_reproduce_exported_moments_through_burnout_and_coast() {
    let (s, f) = (state(), forcing());
    let whole = s.advance_with_exhaust(f, 6.0, || false).unwrap();
    assert_eq!(whole.mechanics.powered_duration_s, 4.0);
    assert_eq!(whole.mechanics.coast_duration_s, 2.0);
    for split in [0.5, 2.0, 4.0, 5.0] {
        let a = s.advance_with_exhaust(f, split, || false).unwrap();
        let b = a
            .mechanics
            .state_after
            .advance_with_exhaust(f, 6.0 - split, || false)
            .unwrap();
        assert_eq!(b.mechanics.state_after.propellant_mass_kg, 0.0);
        vector_near(
            b.mechanics.state_after.velocity_m_s,
            whole.mechanics.state_after.velocity_m_s,
            2e-13,
        );
        vector_near(
            b.mechanics.state_after.position_m,
            whole.mechanics.state_after.position_m,
            2e-13,
        );
        vector_near(
            a.exhaust_momentum_world_kg_m_s
                .add(b.exhaust_momentum_world_kg_m_s),
            whole.exhaust_momentum_world_kg_m_s,
            2e-13,
        );
        near(
            a.exhaust_energy_world_j + b.exhaust_energy_world_j,
            whole.exhaust_energy_world_j,
            2e-13,
        );
    }
}

#[test]
fn tiny_depletion_is_retained_and_unrepresentable_boost_recoil_is_exposed() {
    let s = RocketState {
        dry_mass_kg: 1e16,
        propellant_mass_kg: 1.0,
        position_m: Vec3::ZERO,
        velocity_m_s: Vec3::ZERO,
    };
    let f = ExhaustForcing {
        mass_flow_kg_s: 1e-4,
        relative_momentum_rate_n: Vec3::new(0.1, 0.0, 0.0),
        relative_energy_rate_w: 100.0,
    };
    let step = s.advance_with_exhaust(f, 0.01, || false).unwrap();
    assert!(step.mechanics.state_after.velocity_m_s.x < 0.0);
    assert!(step.mechanics.state_after.propellant_mass_kg < 1.0);
    assert!(step.relative_kinetic_transport_energy_j < 0.0);
    near(step.exhaust_momentum_world_kg_m_s.x, 0.001, 2e-16);
    assert_eq!(step.momentum_residual_kg_m_s, Vec3::ZERO);

    let boosted = RocketState {
        velocity_m_s: Vec3::new(1e20, 0.0, 0.0),
        ..s
    };
    let erased = boosted.advance_with_exhaust(f, 0.01, || false).unwrap();
    assert_eq!(
        erased.mechanics.state_after.velocity_m_s,
        boosted.velocity_m_s
    );
    assert!(erased.momentum_residual_kg_m_s.x.abs() > 0.0009);
    assert!(erased.kinetic_energy_residual_j.abs() > 9e16);

    // The relative speed squared overflows, but its tiny-time velocity
    // variance, actual recoil and transported energy are all representable.
    let extreme = RocketState {
        dry_mass_kg: 1.0,
        propellant_mass_kg: 1e-90,
        position_m: Vec3::ZERO,
        velocity_m_s: Vec3::ZERO,
    }
    .advance_with_exhaust(
        ExhaustForcing {
            mass_flow_kg_s: 1e-100,
            relative_momentum_rate_n: Vec3::new(1e100, 0.0, 0.0),
            relative_energy_rate_w: 1e300,
        },
        1.0,
        || false,
    )
    .unwrap();
    near(extreme.mean_speed_squared_m2_s2, 1e200 / 3.0, 2e-13);
    assert!(extreme.relative_kinetic_transport_energy_j < 0.0);
}

#[test]
fn rounded_inventory_loss_retains_its_advected_momentum_and_energy_defect() {
    let s = RocketState {
        position_m: Vec3::ZERO,
        velocity_m_s: Vec3::new(1e10, 0.0, 0.0),
        dry_mass_kg: 1.0,
        propellant_mass_kg: 1.0,
    };
    let step = s
        .advance_with_exhaust(
            ExhaustForcing {
                mass_flow_kg_s: 0.1,
                relative_momentum_rate_n: Vec3::ZERO,
                relative_energy_rate_w: 1.0,
            },
            1.0,
            || false,
        )
        .unwrap();
    let defect = step.mechanics.burned_propellant_mass_kg
        - (s.propellant_mass_kg - step.mechanics.state_after.propellant_mass_kg);
    assert_ne!(defect, 0.0);
    assert_eq!(step.mass_residual_kg, defect);
    assert_eq!(step.mechanics.state_after.velocity_m_s, s.velocity_m_s);
    assert_eq!(step.momentum_residual_kg_m_s.x, defect * 1e10);
    assert_eq!(step.kinetic_energy_residual_j, 0.5 * defect * 1e20);
    assert!(step.kinetic_energy_residual_j.abs() > 100.0);
}

#[test]
fn invalid_particle_flux_and_cancellation_return_no_state() {
    let (s, f) = (state(), forcing());
    for bad in [
        ExhaustForcing {
            mass_flow_kg_s: -1.0,
            ..f
        },
        ExhaustForcing {
            mass_flow_kg_s: 0.0,
            ..f
        },
        ExhaustForcing {
            relative_momentum_rate_n: Vec3::new(f64::NAN, 0.0, 0.0),
            ..f
        },
        ExhaustForcing {
            relative_energy_rate_w: f64::INFINITY,
            ..f
        },
    ] {
        assert!(s.advance_with_exhaust(bad, 1.0, || false).is_err());
    }
    assert!(s.advance_with_exhaust(f, -1.0, || false).is_err());
    assert_eq!(
        s.advance_with_exhaust(f, 1.0, || true),
        Err(RocketError::Cancelled)
    );
    let unchanged = s.advance_with_exhaust(f, 0.0, || false).unwrap();
    assert_eq!(unchanged.mechanics.state_after, s);
    assert_eq!(unchanged.exhaust_energy_world_j, 0.0);
    let coast = s
        .advance_with_exhaust(
            ExhaustForcing {
                mass_flow_kg_s: 0.0,
                relative_momentum_rate_n: Vec3::ZERO,
                relative_energy_rate_w: 0.0,
            },
            2.0,
            || false,
        )
        .unwrap();
    assert_eq!(coast.mechanics.state_after.velocity_m_s, s.velocity_m_s);
    assert_eq!(coast.exhaust_momentum_world_kg_m_s, Vec3::ZERO);
    assert_eq!(s, state());
}
