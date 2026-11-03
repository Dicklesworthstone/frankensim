//! Independent forced-trajectory integrals and finite-mass mechanical limits.
use crate::Vec3;
use crate::rocket::{
    RocketError, RocketState,
    exhaust::{ExhaustEnvironment, ExhaustForcing, ExhaustStep},
};

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
fn forcing() -> ExhaustForcing {
    ExhaustForcing {
        mass_flow_kg_s: 1.0,
        relative_momentum_rate_n: Vec3::new(90.0, -40.0, 15.0),
        relative_energy_rate_w: 30_000.0,
    }
}
fn environment() -> ExhaustEnvironment {
    ExhaustEnvironment {
        gravity_m_s2: Vec3::new(0.2, -1.3, 0.05),
        external_force_n: Vec3::new(2.0, -4.0, 1.0),
    }
}
fn advance(s: &RocketState, f: ExhaustForcing, e: ExhaustEnvironment, dt: f64) -> ExhaustStep {
    s.advance_with_exhaust_forced(f, e, dt, || false).unwrap()
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
struct Quadrature {
    momentum: Vec3,
    kinetic_export: f64,
    external_impulse: Vec3,
    work: f64,
    mean_velocity: Vec3,
    mean_speed_squared: f64,
    displacement: Vec3,
    final_velocity: Vec3,
}

// Simpson integrates physical q(v+w), q(v.w+|v|²/2), (F+M g), and
// (F+M g).v directly. The independently evaluated velocity uses std log1p;
// no production mass integrals, variance series or exhaust accounts are reused.
// Split at burnout so an endpoint mass-flux discontinuity cannot pollute quadrature.
fn quadrature(s: RocketState, f: ExhaustForcing, e: ExhaustEnvironment, dt: f64) -> Quadrature {
    let q = f.mass_flow_kg_s;
    let b = if q > 0.0 {
        dt.min(s.propellant_mass_kg / q)
    } else {
        0.0
    };
    let c = dt - b;
    let m0 = s.dry_mass_kg + s.propellant_mass_kg;
    let m1 = m0 - q * b;
    let w = if b > 0.0 {
        f.relative_momentum_rate_n.scale(q.recip())
    } else {
        Vec3::ZERO
    };
    let a = if b > 0.0 {
        e.external_force_n.scale(q.recip()).sub(w)
    } else {
        Vec3::ZERO
    };
    let burn_velocity = |t: f64| {
        s.velocity_m_s
            .add(a.scale(-(-q * t / m0).ln_1p()))
            .add(e.gravity_m_s2.scale(t))
    };
    let vb = burn_velocity(b);
    let coast_acceleration = e.external_force_n.scale(m1.recip()).add(e.gravity_m_s2);
    let mut integrals = [Sum::default(); 15];
    const PANELS: usize = 16384;
    for (duration, emitting) in [(b, true), (c, false)] {
        if duration == 0.0 {
            continue;
        }
        let h = duration / PANELS as f64;
        for i in 0..=PANELS {
            let t = i as f64 * h;
            let v = if emitting {
                burn_velocity(t)
            } else {
                vb.add(coast_acceleration.scale(t))
            };
            let mass = if emitting { m0 - q * t } else { m1 };
            let external = e.external_force_n.add(e.gravity_m_s2.scale(mass));
            let p = if emitting {
                v.add(w).scale(q)
            } else {
                Vec3::ZERO
            };
            let kinetic = if emitting {
                q * (v.dot(w) + 0.5 * v.norm_squared())
            } else {
                0.0
            };
            let mean = if emitting { v } else { Vec3::ZERO };
            let mean_squared = if emitting { v.norm_squared() } else { 0.0 };
            let weight = h / 3.0
                * if i == 0 || i == PANELS {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
            for (sum, value) in integrals.iter_mut().zip([
                p.x,
                p.y,
                p.z,
                kinetic,
                external.x,
                external.y,
                external.z,
                external.dot(v),
                mean.x,
                mean.y,
                mean.z,
                mean_squared,
                v.x,
                v.y,
                v.z,
            ]) {
                sum.add(weight * value);
            }
        }
    }
    let value = |i: usize| integrals[i].value;
    Quadrature {
        momentum: Vec3::new(value(0), value(1), value(2)),
        kinetic_export: value(3),
        external_impulse: Vec3::new(value(4), value(5), value(6)),
        work: value(7),
        mean_velocity: if b > 0.0 {
            Vec3::new(value(8), value(9), value(10)).scale(b.recip())
        } else {
            Vec3::ZERO
        },
        mean_speed_squared: if b > 0.0 { value(11) / b } else { 0.0 },
        displacement: Vec3::new(value(12), value(13), value(14)),
        final_velocity: vb.add(coast_acceleration.scale(c)),
    }
}

#[test]
fn forced_burn_and_coast_match_independent_particle_flux_and_work_quadrature() {
    let (s, f, e) = (state(), forcing(), environment());
    for dt in [2.0, 6.5] {
        let expected = quadrature(s, f, e, dt);
        let step = advance(&s, f, e, dt);
        vector_near(
            step.mechanics.state_after.velocity_m_s,
            expected.final_velocity,
            3e-12,
        );
        vector_near(
            step.mechanics.state_after.position_m,
            s.position_m.add(expected.displacement),
            3e-12,
        );
        vector_near(step.exhaust_momentum_world_kg_m_s, expected.momentum, 3e-12);
        near(
            step.relative_kinetic_transport_energy_j,
            expected.kinetic_export,
            3e-12,
        );
        vector_near(step.external_impulse_n_s, expected.external_impulse, 3e-12);
        near(step.external_work_j, expected.work, 3e-12);
        vector_near(step.mean_velocity_world_m_s, expected.mean_velocity, 3e-12);
        near(
            step.mean_speed_squared_m2_s2,
            expected.mean_speed_squared,
            3e-12,
        );
        near(
            step.relative_energy_j,
            f.relative_energy_rate_w * dt.min(5.0),
            1e-15,
        );
        let m1 = s.dry_mass_kg + step.mechanics.state_after.propellant_mass_kg;
        let delta_k = 0.5 * m1 * expected.final_velocity.norm_squared()
            - 0.5 * 8.0 * s.velocity_m_s.norm_squared();
        near(step.kinetic_energy_change_j, delta_k, 3e-12);
        vector_near(
            step.mechanics
                .state_after
                .velocity_m_s
                .scale(m1)
                .sub(s.velocity_m_s.scale(8.0))
                .add(step.exhaust_momentum_world_kg_m_s)
                .sub(step.external_impulse_n_s),
            Vec3::ZERO,
            2e-10,
        );
        assert!((delta_k + expected.kinetic_export - expected.work).abs() < 1e-8);
        assert!(step.momentum_residual_kg_m_s.norm_squared() < 1e-18);
        assert!(step.kinetic_energy_residual_j.abs() < 1e-8);
    }
}

#[test]
fn force_balancing_relative_thrust_keeps_velocity_while_mass_and_energy_escape() {
    let s = RocketState {
        velocity_m_s: Vec3::new(3.0, -2.0, 1.0),
        ..state()
    };
    let w = Vec3::new(12.0, -6.0, 2.0);
    let f = ExhaustForcing {
        mass_flow_kg_s: 0.5,
        relative_momentum_rate_n: w.scale(0.5),
        relative_energy_rate_w: 100.0,
    };
    let e = ExhaustEnvironment {
        gravity_m_s2: Vec3::ZERO,
        external_force_n: f.relative_momentum_rate_n,
    };
    let step = advance(&s, f, e, 2.0);
    vector_near(
        step.mechanics.state_after.velocity_m_s,
        s.velocity_m_s,
        2e-14,
    );
    vector_near(
        step.mechanics.state_after.position_m,
        s.position_m.add(s.velocity_m_s.scale(2.0)),
        2e-14,
    );
    near(step.mechanics.burned_propellant_mass_kg, 1.0, 1e-15);
    vector_near(step.mean_velocity_world_m_s, s.velocity_m_s, 2e-14);
    near(step.mean_speed_squared_m2_s2, 14.0, 2e-14);
    vector_near(
        step.exhaust_momentum_world_kg_m_s,
        s.velocity_m_s.add(w),
        2e-14,
    );
    vector_near(step.external_impulse_n_s, w, 2e-14);
    near(step.external_work_j, 50.0, 2e-14);
    near(step.relative_kinetic_transport_energy_j, 57.0, 2e-14);
    near(step.kinetic_energy_change_j, -7.0, 2e-14);
}

#[test]
fn gravity_uses_changing_mass_and_its_correlation_with_velocity() {
    let s = RocketState {
        velocity_m_s: Vec3::new(1.0, 4.0, -2.0),
        ..state()
    };
    let g = Vec3::new(0.0, -2.0, 0.5);
    let f = ExhaustForcing {
        relative_momentum_rate_n: Vec3::ZERO,
        ..forcing()
    };
    let e = ExhaustEnvironment {
        gravity_m_s2: g,
        external_force_n: Vec3::ZERO,
    };
    let dt = 2.0;
    let step = advance(&s, f, e, dt);
    let displacement = s.velocity_m_s.scale(dt).add(g.scale(0.5 * dt * dt));
    vector_near(
        step.mechanics.state_after.velocity_m_s,
        s.velocity_m_s.add(g.scale(dt)),
        2e-14,
    );
    vector_near(
        step.external_impulse_n_s,
        g.scale(8.0 * dt - 0.5 * dt * dt),
        2e-14,
    );
    let expected_work = s.velocity_m_s.dot(g) * (8.0 * dt - 0.5 * dt * dt)
        + g.norm_squared() * (4.0 * dt * dt - dt * dt * dt / 3.0);
    near(step.external_work_j, expected_work, 3e-14);
    // Average mass times displacement discards the M(t)-v(t) correlation.
    let average_mass_work = 7.0 * g.dot(displacement);
    assert!((step.external_work_j - average_mass_work).abs() > 1.0);
    vector_near(
        step.mean_velocity_world_m_s,
        s.velocity_m_s.add(g.scale(dt / 2.0)),
        2e-14,
    );
    near(
        step.mean_speed_squared_m2_s2,
        s.velocity_m_s.norm_squared()
            + dt * s.velocity_m_s.dot(g)
            + dt * dt * g.norm_squared() / 3.0,
        2e-14,
    );
    assert!(step.kinetic_energy_residual_j.abs() < 1e-11);
}

#[test]
fn closed_coast_and_split_burnout_preserve_external_impulse_and_work() {
    let (s, f, e) = (state(), forcing(), environment());
    let coast = advance(
        &s,
        ExhaustForcing {
            mass_flow_kg_s: 0.0,
            relative_momentum_rate_n: Vec3::ZERO,
            relative_energy_rate_w: 0.0,
        },
        e,
        2.0,
    );
    let acceleration = e.external_force_n.scale(1.0 / 8.0).add(e.gravity_m_s2);
    let displacement = s.velocity_m_s.scale(2.0).add(acceleration.scale(2.0));
    vector_near(
        coast.mechanics.state_after.velocity_m_s,
        s.velocity_m_s.add(acceleration.scale(2.0)),
        2e-14,
    );
    vector_near(coast.external_impulse_n_s, acceleration.scale(16.0), 2e-14);
    near(
        coast.external_work_j,
        acceleration.scale(8.0).dot(displacement),
        2e-14,
    );
    assert_eq!(coast.exhaust_momentum_world_kg_m_s, Vec3::ZERO);
    assert_eq!(coast.exhaust_energy_world_j, 0.0);
    assert_eq!(coast.mean_velocity_world_m_s, Vec3::ZERO);
    assert_eq!(coast.mean_speed_squared_m2_s2, 0.0);
    assert!(coast.kinetic_energy_residual_j.abs() < 1e-11);
    let whole = advance(&s, f, e, 6.5);
    assert_eq!(whole.mechanics.powered_duration_s, 5.0);
    assert_eq!(whole.mechanics.coast_duration_s, 1.5);
    for split in [0.125, 2.5, 5.0, 5.75] {
        let a = advance(&s, f, e, split);
        let b = advance(&a.mechanics.state_after, f, e, 6.5 - split);
        vector_near(
            b.mechanics.state_after.position_m,
            whole.mechanics.state_after.position_m,
            3e-12,
        );
        vector_near(
            b.mechanics.state_after.velocity_m_s,
            whole.mechanics.state_after.velocity_m_s,
            3e-12,
        );
        assert_eq!(b.mechanics.state_after.propellant_mass_kg, 0.0);
        vector_near(
            a.exhaust_momentum_world_kg_m_s
                .add(b.exhaust_momentum_world_kg_m_s),
            whole.exhaust_momentum_world_kg_m_s,
            3e-12,
        );
        near(
            a.exhaust_energy_world_j + b.exhaust_energy_world_j,
            whole.exhaust_energy_world_j,
            3e-12,
        );
        vector_near(
            a.external_impulse_n_s.add(b.external_impulse_n_s),
            whole.external_impulse_n_s,
            3e-12,
        );
        near(
            a.external_work_j + b.external_work_j,
            whole.external_work_j,
            3e-12,
        );
    }
}

#[test]
fn forced_galilean_crosswalk_transforms_work_and_keeps_the_source_energy_datum() {
    let (s, f, e, dt) = (state(), forcing(), environment(), 6.5);
    let base = advance(&s, f, e, dt);
    let boost = Vec3::new(20.0, -10.0, 5.0);
    let moved = advance(
        &RocketState {
            velocity_m_s: s.velocity_m_s.add(boost),
            ..s
        },
        f,
        e,
        dt,
    );
    let dm = base.mechanics.burned_propellant_mass_kg;
    vector_near(
        moved.mechanics.state_after.position_m,
        base.mechanics.state_after.position_m.add(boost.scale(dt)),
        3e-12,
    );
    vector_near(
        moved.exhaust_momentum_world_kg_m_s,
        base.exhaust_momentum_world_kg_m_s.add(boost.scale(dm)),
        3e-12,
    );
    near(
        moved.relative_kinetic_transport_energy_j,
        base.relative_kinetic_transport_energy_j
            + boost.dot(base.exhaust_momentum_world_kg_m_s)
            + 0.5 * dm * boost.norm_squared(),
        3e-12,
    );
    assert_eq!(moved.external_impulse_n_s, base.external_impulse_n_s);
    near(
        moved.external_work_j,
        base.external_work_j + boost.dot(base.external_impulse_n_s),
        3e-12,
    );
    vector_near(
        moved.mean_velocity_world_m_s,
        base.mean_velocity_world_m_s.add(boost),
        3e-12,
    );
    near(
        moved.mean_speed_squared_m2_s2,
        base.mean_speed_squared_m2_s2
            + 2.0 * boost.dot(base.mean_velocity_world_m_s)
            + boost.norm_squared(),
        3e-12,
    );
    let shifted = advance(
        &s,
        ExhaustForcing {
            relative_energy_rate_w: -1e100,
            ..f
        },
        e,
        dt,
    );
    assert_eq!(
        shifted.relative_kinetic_transport_energy_j,
        base.relative_kinetic_transport_energy_j
    );
    assert_eq!(shifted.external_work_j, base.external_work_j);
    assert_eq!(
        shifted.kinetic_energy_residual_j,
        base.kinetic_energy_residual_j
    );
    assert!(shifted.exhaust_energy_world_j < 0.0);
    let free = s.advance_with_exhaust(f, 2.0, || false).unwrap();
    let zero = advance(&s, f, ExhaustEnvironment::ZERO, 2.0);
    assert_eq!(zero.mechanics, free.mechanics);
    vector_near(
        zero.exhaust_momentum_world_kg_m_s,
        free.exhaust_momentum_world_kg_m_s,
        3e-12,
    );
    near(
        zero.exhaust_energy_world_j,
        free.exhaust_energy_world_j,
        3e-12,
    );
}

#[test]
fn tiny_mass_fraction_retains_variance_and_does_not_form_force_divided_by_rate() {
    let tiny = RocketState {
        dry_mass_kg: 1.0,
        propellant_mass_kg: 1e-290,
        position_m: Vec3::ZERO,
        velocity_m_s: Vec3::ZERO,
    };
    let f = ExhaustForcing {
        mass_flow_kg_s: 1e-300,
        relative_momentum_rate_n: Vec3::new(1e-295, 0.0, 0.0),
        relative_energy_rate_w: 1e-280,
    };
    let e = ExhaustEnvironment {
        gravity_m_s2: Vec3::ZERO,
        external_force_n: Vec3::new(1e10, 0.0, 0.0),
    };
    assert!((e.external_force_n.x / f.mass_flow_kg_s).is_infinite());
    let step = advance(&tiny, f, e, 1.0);
    near(step.mechanics.state_after.velocity_m_s.x, 1e10, 3e-14);
    near(step.mean_velocity_world_m_s.x, 5e9, 3e-14);
    near(step.mean_speed_squared_m2_s2, 1e20 / 3.0, 3e-14);
    assert!(step.exhaust_energy_world_j > 0.0);
    assert!(step.mechanics.state_after.propellant_mass_kg < tiny.propellant_mass_kg);

    // Initially balanced force/gravity leaves a small quadratic drift as mass
    // drains. Its leading <v>=q/12 and <v²>=q²/80 are independent small-q limits.
    let s = RocketState {
        dry_mass_kg: 1.0,
        propellant_mass_kg: 1.0,
        position_m: Vec3::ZERO,
        velocity_m_s: Vec3::ZERO,
    };
    let q = 1e-7;
    let small = advance(
        &s,
        ExhaustForcing {
            mass_flow_kg_s: q,
            relative_momentum_rate_n: Vec3::ZERO,
            relative_energy_rate_w: q,
        },
        ExhaustEnvironment {
            gravity_m_s2: Vec3::new(-1.0, 0.0, 0.0),
            external_force_n: Vec3::new(2.0, 0.0, 0.0),
        },
        1.0,
    );
    assert!((small.mean_velocity_world_m_s.x / (q / 12.0) - 1.0).abs() < 2e-7);
    assert!((small.mean_speed_squared_m2_s2 / (q * q / 80.0) - 1.0).abs() < 2e-7);
}

#[test]
fn forced_admission_and_cancellation_return_no_partial_mechanical_state() {
    let (s, f, e) = (state(), forcing(), environment());
    let original = s;
    for environment in [
        ExhaustEnvironment {
            gravity_m_s2: Vec3::new(f64::NAN, 0.0, 0.0),
            ..e
        },
        ExhaustEnvironment {
            external_force_n: Vec3::new(0.0, f64::INFINITY, 0.0),
            ..e
        },
    ] {
        assert!(
            s.advance_with_exhaust_forced(f, environment, 1.0, || false)
                .is_err()
        );
    }
    assert_eq!(
        s.advance_with_exhaust_forced(f, e, 1.0, || true),
        Err(RocketError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        s.advance_with_exhaust_forced(f, e, 1.0, || {
            polls += 1;
            polls >= 3
        }),
        Err(RocketError::Cancelled)
    );
    assert_eq!(s, original);
    let zero = advance(&s, f, e, 0.0);
    assert_eq!(zero.mechanics.state_after, s);
    assert_eq!(zero.external_impulse_n_s, Vec3::ZERO);
    assert_eq!(zero.external_work_j, 0.0);
}
