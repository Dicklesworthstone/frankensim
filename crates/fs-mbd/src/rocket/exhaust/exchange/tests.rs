//! Independent physical flux integrals and finite bidirectional limits.
use super::*;
use crate::rocket::exhaust::ExhaustForcing;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * actual.abs().max(expected.abs()).max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}
fn vector_near(actual: Vec3, expected: Vec3, tolerance: f64) {
    near(actual.x, expected.x, tolerance);
    near(actual.y, expected.y, tolerance);
    near(actual.z, expected.z, tolerance);
}
fn state() -> RocketState {
    RocketState {
        position_m: Vec3::new(2.0, -3.0, 5.0),
        velocity_m_s: Vec3::new(7.0, 12.0, -3.0),
        dry_mass_kg: 3.0,
        propellant_mass_kg: 5.0,
    }
}
fn stream(rate: f64, velocity: Vec3) -> ParticleStreamMoments {
    ParticleStreamMoments {
        mass_flow_kg_s: rate,
        relative_momentum_rate_n: velocity.scale(rate),
        relative_energy_rate_w: rate * (2000.0 + 0.5 * velocity.norm_squared()),
    }
}
fn forcing(incoming: f64, outgoing: f64) -> ParticleExchangeForcing {
    ParticleExchangeForcing {
        incoming: stream(incoming, Vec3::new(30.0, -12.0, 8.0)),
        outgoing: stream(outgoing, Vec3::new(-70.0, 8.0, 10.0)),
    }
}
fn environment() -> ExhaustEnvironment {
    ExhaustEnvironment {
        gravity_m_s2: Vec3::new(0.3, -2.0, 0.1),
        external_force_n: Vec3::new(3.0, -5.0, 0.8),
    }
}
fn advance(
    s: RocketState,
    f: ParticleExchangeForcing,
    e: ExhaustEnvironment,
    dt: f64,
) -> ParticleExchangeStep {
    s.advance_with_particle_exchange(f, e, dt, || false)
        .unwrap()
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let adjusted = value - self.correction;
        let next = self.value + adjusted;
        self.correction = (next - self.value) - adjusted;
        self.value = next;
    }
}

// Independently integrate q*v+P, v.P+q*v²/2, F+M*g, and their work.
// The velocity uses std log1p, not native mass/variance antiderivatives.
fn quadrature(
    s: RocketState,
    f: ParticleExchangeForcing,
    e: ExhaustEnvironment,
    dt: f64,
) -> [f64; 16] {
    let m0 = s.dry_mass_kg + s.propellant_mass_kg;
    let rate = f.incoming.mass_flow_kg_s - f.outgoing.mass_flow_kg_s;
    let force = e
        .external_force_n
        .add(f.incoming.relative_momentum_rate_n)
        .sub(f.outgoing.relative_momentum_rate_n);
    let mut integrals = [Sum::default(); 16];
    const PANELS: usize = 8192;
    let h = dt / PANELS as f64;
    for i in 0..=PANELS {
        let t = i as f64 * h;
        let coefficient = if rate == 0.0 {
            t / m0
        } else {
            (rate * t / m0).ln_1p() / rate
        };
        let v = s
            .velocity_m_s
            .add(force.scale(coefficient))
            .add(e.gravity_m_s2.scale(t));
        let pi = v
            .scale(f.incoming.mass_flow_kg_s)
            .add(f.incoming.relative_momentum_rate_n);
        let po = v
            .scale(f.outgoing.mass_flow_kg_s)
            .add(f.outgoing.relative_momentum_rate_n);
        let ei = v.dot(f.incoming.relative_momentum_rate_n)
            + 0.5 * f.incoming.mass_flow_kg_s * v.norm_squared();
        let eo = v.dot(f.outgoing.relative_momentum_rate_n)
            + 0.5 * f.outgoing.mass_flow_kg_s * v.norm_squared();
        let external = e.external_force_n.add(e.gravity_m_s2.scale(m0 + rate * t));
        let weight = if i == 0 || i == PANELS {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let values = [
            pi.x,
            pi.y,
            pi.z,
            po.x,
            po.y,
            po.z,
            ei,
            eo,
            external.x,
            external.y,
            external.z,
            external.dot(v),
            v.x,
            v.y,
            v.z,
            v.norm_squared(),
        ];
        for (sum, value) in integrals.iter_mut().zip(values) {
            sum.add(weight * value);
        }
    }
    integrals.map(|sum| sum.value * (h / 3.0))
}

#[test]
fn gain_loss_and_throughflow_match_independent_physical_quadrature() {
    let (s, e, dt) = (state(), environment(), 1.3);
    for (qi, qo) in [(1.4, 0.3), (0.2, 1.1), (0.8, 0.8)] {
        let f = forcing(qi, qo);
        let step = advance(s, f, e, dt);
        let q = quadrature(s, f, e, dt);
        vector_near(
            step.incoming.momentum_world_kg_m_s,
            Vec3::new(q[0], q[1], q[2]),
            3e-12,
        );
        vector_near(
            step.outgoing.momentum_world_kg_m_s,
            Vec3::new(q[3], q[4], q[5]),
            3e-12,
        );
        near(step.incoming.kinetic_transport_energy_j, q[6], 3e-12);
        near(step.outgoing.kinetic_transport_energy_j, q[7], 3e-12);
        vector_near(
            step.external_impulse_n_s,
            Vec3::new(q[8], q[9], q[10]),
            3e-12,
        );
        near(step.external_work_j, q[11], 3e-12);
        let mean = Vec3::new(q[12], q[13], q[14]).scale(dt.recip());
        vector_near(step.mean_velocity_world_m_s, mean, 3e-12);
        vector_near(
            step.state_after.position_m,
            s.position_m.add(mean.scale(dt)),
            3e-12,
        );
        near(step.mean_speed_squared_m2_s2, q[15] / dt, 3e-12);
        near(
            step.incoming.energy_world_j,
            f.incoming.relative_energy_rate_w * dt + q[6],
            3e-12,
        );
        near(
            step.outgoing.energy_world_j,
            f.outgoing.relative_energy_rate_w * dt + q[7],
            3e-12,
        );
        near(step.mass_residual_kg, 0.0, 3e-15);
        vector_near(step.momentum_residual_kg_m_s, Vec3::ZERO, 2e-11);
        near(step.kinetic_energy_residual_j, 0.0, 2e-9);
        near(
            step.kinetic_energy_change_j
                + step.retained_relative_energy_j
                + step.outgoing.energy_world_j
                - step.incoming.energy_world_j
                - step.external_work_j,
            0.0,
            2e-9,
        );
    }
    assert_eq!(s, state());
}

#[test]
fn equal_gross_flow_keeps_force_and_incoming_relative_kinetic_thermalization() {
    let s = state();
    let f = ParticleExchangeForcing {
        incoming: stream(2.0, Vec3::new(10.0, 0.0, 0.0)),
        outgoing: stream(2.0, Vec3::new(-5.0, 0.0, 0.0)),
    };
    let dt = 0.4;
    let step = advance(s, f, ExhaustEnvironment::ZERO, dt);
    assert_eq!(step.state_after.propellant_mass_kg, s.propellant_mass_kg);
    near(step.incoming.mass_kg, 0.8, 1e-15);
    near(step.outgoing.mass_kg, 0.8, 1e-15);
    vector_near(
        step.state_after.velocity_m_s,
        s.velocity_m_s.add(Vec3::new(1.5, 0.0, 0.0)),
        2e-14,
    );
    vector_near(
        step.state_after.position_m,
        s.position_m
            .add(s.velocity_m_s.scale(dt))
            .add(Vec3::new(0.3, 0.0, 0.0)),
        2e-14,
    );
    near(
        step.retained_relative_energy_j,
        0.5 * 0.8 * (100.0 - 25.0),
        2e-14,
    );
    let intake = advance(
        s,
        ParticleExchangeForcing {
            outgoing: ParticleStreamMoments::ZERO,
            ..f
        },
        ExhaustEnvironment::ZERO,
        dt,
    );
    near(
        intake.retained_relative_energy_j,
        0.8 * 2000.0 + 0.5 * 0.8 * 100.0,
        2e-14,
    );
    vector_near(
        intake.state_after.velocity_m_s,
        s.velocity_m_s.add(Vec3::new(10.0 * 1.1_f64.ln(), 0.0, 0.0)),
        2e-14,
    );
}

#[test]
fn one_way_exchange_matches_original_exhaust_and_held_subdivision_composes() {
    let (s, e, dt) = (state(), environment(), 1.3);
    let f = forcing(0.0, 0.7);
    let step = advance(s, f, e, dt);
    let original = s
        .advance_with_exhaust_forced(
            ExhaustForcing {
                mass_flow_kg_s: f.outgoing.mass_flow_kg_s,
                relative_momentum_rate_n: f.outgoing.relative_momentum_rate_n,
                relative_energy_rate_w: f.outgoing.relative_energy_rate_w,
            },
            e,
            dt,
            || false,
        )
        .unwrap();
    vector_near(
        step.state_after.velocity_m_s,
        original.mechanics.state_after.velocity_m_s,
        2e-13,
    );
    vector_near(
        step.state_after.position_m,
        original.mechanics.state_after.position_m,
        2e-13,
    );
    vector_near(
        step.outgoing.momentum_world_kg_m_s,
        original.exhaust_momentum_world_kg_m_s,
        2e-13,
    );
    near(
        step.outgoing.energy_world_j,
        original.exhaust_energy_world_j,
        2e-13,
    );
    near(step.external_work_j, original.external_work_j, 2e-13);
    for f in [forcing(1.4, 0.3), forcing(0.2, 1.1), forcing(0.8, 0.8)] {
        let whole = advance(s, f, e, dt);
        let first = advance(s, f, e, 0.5);
        let second = advance(first.state_after, f, e, 0.8);
        vector_near(
            second.state_after.velocity_m_s,
            whole.state_after.velocity_m_s,
            3e-13,
        );
        vector_near(
            second.state_after.position_m,
            whole.state_after.position_m,
            3e-13,
        );
        for (a, b, c) in [
            (first.incoming, second.incoming, whole.incoming),
            (first.outgoing, second.outgoing, whole.outgoing),
        ] {
            near(a.mass_kg + b.mass_kg, c.mass_kg, 3e-14);
            vector_near(
                a.momentum_world_kg_m_s.add(b.momentum_world_kg_m_s),
                c.momentum_world_kg_m_s,
                3e-13,
            );
            near(a.energy_world_j + b.energy_world_j, c.energy_world_j, 3e-13);
        }
        vector_near(
            first.external_impulse_n_s.add(second.external_impulse_n_s),
            whole.external_impulse_n_s,
            3e-13,
        );
        near(
            first.external_work_j + second.external_work_j,
            whole.external_work_j,
            3e-13,
        );
    }
}

#[test]
fn closed_coast_zero_duration_and_late_refusals_publish_no_partial_state() {
    let (s, e, dt) = (state(), environment(), 0.7);
    let closed = ParticleExchangeForcing {
        incoming: ParticleStreamMoments::ZERO,
        outgoing: ParticleStreamMoments::ZERO,
    };
    let step = advance(s, closed, e, dt);
    let acceleration = e.external_force_n.scale(1.0 / 8.0).add(e.gravity_m_s2);
    vector_near(
        step.state_after.velocity_m_s,
        s.velocity_m_s.add(acceleration.scale(dt)),
        2e-14,
    );
    vector_near(
        step.state_after.position_m,
        s.position_m
            .add(s.velocity_m_s.scale(dt))
            .add(acceleration.scale(0.5 * dt * dt)),
        2e-14,
    );
    assert_eq!(step.mean_velocity_world_m_s, Vec3::ZERO);
    assert_eq!(step.mean_speed_squared_m2_s2, 0.0);
    near(
        step.external_work_j,
        step.external_impulse_n_s
            .dot(s.velocity_m_s.add(acceleration.scale(0.5 * dt))),
        2e-14,
    );
    let huge = RocketState {
        velocity_m_s: Vec3::new(1e200, 0.0, 0.0),
        ..s
    };
    let zero = advance(huge, forcing(1.0, 0.5), e, 0.0);
    assert_eq!(zero.state_after, huge);
    assert_eq!(zero.kinetic_energy_change_j, 0.0);
    assert!(matches!(
        s.advance_with_particle_exchange(forcing(0.0, 6.0), e, 1.0, || false),
        Err(RocketError::InvalidInput(_))
    ));
    let invalid = ParticleExchangeForcing {
        incoming: ParticleStreamMoments {
            relative_momentum_rate_n: Vec3::new(1.0, 0.0, 0.0),
            ..ParticleStreamMoments::ZERO
        },
        ..closed
    };
    assert!(matches!(
        s.advance_with_particle_exchange(invalid, e, 1.0, || false),
        Err(RocketError::InvalidInput(_))
    ));
    assert!(matches!(
        s.advance_with_particle_exchange(forcing(-1.0, 0.0), e, 1.0, || false),
        Err(RocketError::InvalidInput(_))
    ));
    assert!(matches!(
        s.advance_with_particle_exchange(closed, e, f64::NAN, || false),
        Err(RocketError::NonFiniteInput(_))
    ));
    let mut calls = 0;
    assert_eq!(
        s.advance_with_particle_exchange(forcing(1.0, 0.5), e, 1.0, || {
            calls += 1;
            calls >= 6
        }),
        Err(RocketError::Cancelled)
    );
    assert!(calls >= 6);
    assert_eq!(s, state());
}

#[test]
fn galilean_stream_work_and_caloric_datum_transforms_remain_separate() {
    let (s, f, e, dt) = (state(), forcing(1.1, 0.3), environment(), 0.9);
    let original = advance(s, f, e, dt);
    let boost = Vec3::new(120.0, -35.0, 50.0);
    let translation = Vec3::new(1e4, -2e4, 3e4);
    let shifted = advance(
        RocketState {
            position_m: s.position_m.add(translation),
            velocity_m_s: s.velocity_m_s.add(boost),
            ..s
        },
        f,
        e,
        dt,
    );
    vector_near(
        shifted.state_after.velocity_m_s,
        original.state_after.velocity_m_s.add(boost),
        3e-13,
    );
    vector_near(
        shifted.state_after.position_m,
        original
            .state_after
            .position_m
            .add(translation)
            .add(boost.scale(dt)),
        3e-13,
    );
    for (a, b) in [
        (original.incoming, shifted.incoming),
        (original.outgoing, shifted.outgoing),
    ] {
        vector_near(
            b.momentum_world_kg_m_s,
            a.momentum_world_kg_m_s.add(boost.scale(a.mass_kg)),
            3e-13,
        );
        near(
            b.energy_world_j,
            a.energy_world_j
                + boost.dot(a.momentum_world_kg_m_s)
                + 0.5 * a.mass_kg * boost.norm_squared(),
            3e-13,
        );
    }
    near(
        shifted.external_work_j,
        original.external_work_j + boost.dot(original.external_impulse_n_s),
        3e-13,
    );
    near(
        shifted.retained_relative_energy_j,
        original.retained_relative_energy_j,
        0.0,
    );
    let datum = 1e12;
    let mut changed = f;
    changed.incoming.relative_energy_rate_w += f.incoming.mass_flow_kg_s * datum;
    changed.outgoing.relative_energy_rate_w += f.outgoing.mass_flow_kg_s * datum;
    let reference = advance(s, changed, e, dt);
    assert_eq!(reference.state_after, original.state_after);
    assert_eq!(
        reference.incoming.kinetic_transport_energy_j,
        original.incoming.kinetic_transport_energy_j
    );
    assert_eq!(
        reference.outgoing.kinetic_transport_energy_j,
        original.outgoing.kinetic_transport_energy_j
    );
    assert_eq!(
        reference.kinetic_energy_residual_j,
        original.kinetic_energy_residual_j
    );
    near(
        reference.retained_relative_energy_j,
        original.retained_relative_energy_j
            + datum * (original.incoming.mass_kg - original.outgoing.mass_kg),
        3e-15,
    );
}

#[test]
fn tiny_mass_changes_preserve_force_cancellation_and_expose_endpoint_rounding() {
    for gain in [false, true] {
        let s = RocketState {
            position_m: Vec3::ZERO,
            velocity_m_s: Vec3::ZERO,
            dry_mass_kg: 1.0,
            propellant_mass_kg: 1.0,
        };
        let rate = 1e-10;
        let flux = ParticleStreamMoments {
            mass_flow_kg_s: rate,
            ..ParticleStreamMoments::ZERO
        };
        let f = if gain {
            ParticleExchangeForcing {
                incoming: flux,
                outgoing: ParticleStreamMoments::ZERO,
            }
        } else {
            ParticleExchangeForcing {
                incoming: ParticleStreamMoments::ZERO,
                outgoing: flux,
            }
        };
        let high = if gain { 2.0 + rate } else { 2.0 };
        let x = rate / high;
        let step = advance(
            s,
            f,
            ExhaustEnvironment {
                external_force_n: Vec3::new(high, 0.0, 0.0),
                gravity_m_s2: Vec3::new(-1.0, 0.0, 0.0),
            },
            1.0,
        );
        assert!((step.state_after.velocity_m_s.x / (0.5 * x) - 1.0).abs() < 2e-6);
        let expected_mean = x / if gain { 3.0 } else { 6.0 };
        let expected_squared = x * x * if gain { 2.0 / 15.0 } else { 1.0 / 20.0 };
        assert!((step.mean_velocity_world_m_s.x / expected_mean - 1.0).abs() < 2e-6);
        assert!((step.mean_speed_squared_m2_s2 / expected_squared - 1.0).abs() < 5e-6);
    }
    let tiny = RocketState {
        position_m: Vec3::ZERO,
        velocity_m_s: Vec3::ZERO,
        dry_mass_kg: 1.0,
        propellant_mass_kg: 1e-299,
    };
    let tiny_flux = ParticleExchangeForcing {
        incoming: ParticleStreamMoments {
            mass_flow_kg_s: 2e-300,
            ..ParticleStreamMoments::ZERO
        },
        outgoing: ParticleStreamMoments {
            mass_flow_kg_s: 1e-300,
            ..ParticleStreamMoments::ZERO
        },
    };
    // F/(qi-qo) would overflow, but its finite acceleration does not need it.
    let tiny_step = advance(
        tiny,
        tiny_flux,
        ExhaustEnvironment {
            external_force_n: Vec3::new(1e20, 0.0, 0.0),
            gravity_m_s2: Vec3::ZERO,
        },
        1.0,
    );
    near(tiny_step.state_after.velocity_m_s.x, 1e20, 2e-14);
    // time/mass underflows here; the physical force*time/mass does not.
    let massive = RocketState {
        dry_mass_kg: 1e300,
        propellant_mass_kg: 0.0,
        ..tiny
    };
    let closed = ParticleExchangeForcing {
        incoming: ParticleStreamMoments::ZERO,
        outgoing: ParticleStreamMoments::ZERO,
    };
    let brief = advance(
        massive,
        closed,
        ExhaustEnvironment {
            external_force_n: Vec3::new(1e300, 0.0, 0.0),
            gravity_m_s2: Vec3::ZERO,
        },
        1e-100,
    );
    assert!((brief.state_after.velocity_m_s.x / 1e-100 - 1.0).abs() < 2e-14);
    near(brief.external_work_j, 5e99, 3e-14);
    let boosted = RocketState {
        position_m: Vec3::ZERO,
        velocity_m_s: Vec3::new(1e16, 0.0, 0.0),
        dry_mass_kg: 1.0,
        propellant_mass_kg: 1.0,
    };
    let through = ParticleExchangeForcing {
        incoming: ParticleStreamMoments {
            mass_flow_kg_s: 1.0,
            relative_momentum_rate_n: Vec3::new(1.0, 0.0, 0.0),
            relative_energy_rate_w: 0.5,
        },
        outgoing: ParticleStreamMoments {
            mass_flow_kg_s: 1.0,
            ..ParticleStreamMoments::ZERO
        },
    };
    let lost = advance(boosted, through, ExhaustEnvironment::ZERO, 1.0);
    assert_eq!(lost.state_after.velocity_m_s, boosted.velocity_m_s);
    near(lost.momentum_residual_kg_m_s.x, -1.0, 0.0);
    let rounded = advance(state(), forcing(0.1, 0.0), ExhaustEnvironment::ZERO, 0.1);
    assert_ne!(rounded.mass_residual_kg, 0.0);
    near(
        rounded.mass_residual_kg,
        rounded.state_after.propellant_mass_kg - state().propellant_mass_kg
            + rounded.outgoing.mass_kg
            - rounded.incoming.mass_kg,
        0.0,
    );
}
